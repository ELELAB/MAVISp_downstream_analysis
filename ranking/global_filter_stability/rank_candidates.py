#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ---------- Helpers ----------
def normalize_colname(s: str) -> str:
    """Lowercase, replace NBSP/ZWSP, collapse whitespace/newlines."""
    if s is None:
        return ""
    s = str(s).replace("\u00A0", " ").replace("\u200B", "")
    s = re.sub(r"\s+", " ", s, flags=re.UNICODE)
    return s.strip().lower()


def resolve_column(df: pd.DataFrame, target_name: str) -> str:
    """Resolve a column name robustly against MAVISp headers."""
    if target_name in df.columns:
        return target_name
    norm_map = {normalize_colname(c): c for c in df.columns}
    key = normalize_colname(target_name)
    if key in norm_map:
        return norm_map[key]
    raise KeyError(f"Column not found: '{target_name}'. Available:\n  - " + "\n  - ".join(df.columns))


# ---------- Config parsing ----------
class FeatureSpec:
    def __init__(self, kind: str, col: str, weight: float):
        self.kind = kind  # "NUM" or "CAT"
        self.col = col
        self.weight = weight


class NumSpec(FeatureSpec):
    def __init__(self, col: str, weight: float, direction: str = "high", normalize: str = "minmax"):
        super().__init__("NUM", col, weight)
        self.direction = direction  # "high" or "low"
        self.normalize = normalize  # "minmax" or "none"


class CatSpec(FeatureSpec):
    def __init__(self, col: str, weight: float, order: List[str], rest: str = "last"):
        super().__init__("CAT", col, weight)
        self.order = order  # best -> worst
        self.rest = rest    # "last" or "zero"


def parse_config_line(line: str) -> Optional[FeatureSpec]:
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    # Split head (kind: col) and tail (semicolon-separated key=val)
    head, *rest = [x.strip() for x in line.split(";", maxsplit=10)]
    # head like: NUM: Column name
    m = re.match(r"^(NUM|CAT)\s*:\s*(.+)$", head, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Invalid config line (head): {line}")
    kind = m.group(1).upper()
    col = m.group(2).strip()

    params = {}
    for piece in rest:
        if not piece:
            continue
        if "=" not in piece:
            raise ValueError(f"Invalid parameter (missing '='): {piece} in line: {line}")
        k, v = piece.split("=", 1)
        params[k.strip().lower()] = v.strip()

    weight = float(params.get("weight", "1.0"))
    if kind == "NUM":
        direction = params.get("direction", "high").lower()
        if direction not in {"high", "low"}:
            raise ValueError(f"direction must be high|low in: {line}")
        normalize = params.get("normalize", "minmax").lower()
        if normalize not in {"minmax", "none"}:
            raise ValueError(f"normalize must be minmax|none in: {line}")
        return NumSpec(col=col, weight=weight, direction=direction, normalize=normalize)

    else:  # CAT
        order_str = params.get("order", "")
        if not order_str:
            raise ValueError(f"CAT requires order=... in: {line}")
        order = [x.strip() for x in order_str.split(">") if x.strip()]
        rest = params.get("rest", "last").lower()
        if rest not in {"last", "zero"}:
            raise ValueError(f"rest must be last|zero in: {line}")
        return CatSpec(col=col, weight=weight, order=order, rest=rest)


def parse_config_file(path: Path) -> List[FeatureSpec]:
    specs: List[FeatureSpec] = []
    with path.open("r", encoding="utf-8") as fh:
        for i, raw in enumerate(fh, 1):
            try:
                spec = parse_config_line(raw)
                if spec:
                    specs.append(spec)
            except Exception as e:
                raise ValueError(f"Config parse error at line {i}: {e}") from e
    if not specs:
        raise ValueError("Config file produced no features.")
    return specs


# ---------- Scoring ----------
def minmax_normalize(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    xmin = np.nanmin(x.values)
    xmax = np.nanmax(x.values)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax == xmin:
        # Constant or non-finite -> neutral 0.5
        return pd.Series([0.5] * len(x), index=x.index)
    return (x - xmin) / (xmax - xmin)


def apply_num_feature(series: pd.Series, spec: NumSpec) -> pd.Series:
    s = series.copy()
    if spec.normalize == "minmax":
        s = minmax_normalize(s)
    else:
        # none: try to coerce to float and clip to [0,1] if it looks like probs;
        # otherwise pass-through then clip to [0,1] after direction.
        s = pd.to_numeric(s, errors="coerce")
        # If many NaNs or values outside [0,1], just fill NaN with 0.5 later.
    # Direction handling
    if spec.direction == "low":
        s = 1.0 - s
    # Fill NaNs with neutral
    s = s.fillna(0.5)
    s = s.clip(lower=0.0, upper=1.0)
    return s * float(spec.weight)


def apply_cat_feature(series: pd.Series, spec: CatSpec) -> pd.Series:
    # Map listed categories: best->1.0, worst->0.0
    n = max(len(spec.order) - 1, 1)
    scores: Dict[str, float] = {}
    for i, val in enumerate(spec.order):
        score = 1.0 - (i / n) if n > 0 else 1.0
        scores[val.strip().casefold()] = score

    def map_val(v: str) -> float:
        if v is None:
            return 0.0
        key = str(v).strip().casefold()
        if key in scores:
            return scores[key]
        return 0.0 if spec.rest == "zero" else 0.0  # both rest behaviors map to 0.0; keep simple/strict

    mapped = series.apply(map_val)
    mapped = mapped.fillna(0.0).clip(0.0, 1.0)
    return mapped * float(spec.weight)


# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Rank MAVISp variants for experimental validation.\n"
            "Pre-filter: AlphaMissense pathogenic AND |DeMaSk score| >= threshold (both columns configurable).\n"
            "Ranking: user-defined features from a config TXT (numeric or categorical), weighted and normalized."
        )
    )
    p.add_argument("-i", "--input", required=True, help="Input MAVISp CSV")
    p.add_argument("-o", "--output", required=True, help="Output ranked CSV")
    p.add_argument("--config", required=True, help="TXT config with feature scoring rules")

    # Pre-filter columns/params
    p.add_argument("--am-class-col", default="AlphaMissense classification",
                   help="Column with AlphaMissense classification (default: 'AlphaMissense classification').")
    p.add_argument("--am-class-value", default="pathogenic",
                   help="Value to require in AlphaMissense classification (default: 'pathogenic').")
    p.add_argument("--demask-col", default="DeMaSk delta fitness",
                   help="Column for DeMaSk score (e.g., 'DeMaSk delta fitness').")
    p.add_argument("--demask-threshold", type=float, default=0.25,
                   help="Absolute DeMaSk threshold; keep rows with |score| >= threshold (default: 0.25).")

    # Optional: include contributions/tiebreakers in output
    p.add_argument("--mutation-col", default="Mutation", help="Mutation column name for display")
    p.add_argument("--quote-all", action="store_true", help="Quote all fields in the output CSV")
    p.add_argument("--top", type=int, default=0, help="Print top N to screen (0 = skip).")
    p.add_argument("--debug-cols", action="store_true", help="Show columns (original -> normalized) and exit.")
    return p.parse_args()


def main():
    args = parse_args()

    in_path = Path(args.input)
    cfg_path = Path(args.config)
    out_path = Path(args.output)

    # Load data
    try:
        df = pd.read_csv(in_path, dtype=str, keep_default_na=False, na_values=[],
                         encoding="utf-8-sig", engine="python")
    except Exception as e:
        print(f"[ERROR] Failed to read CSV: {e}", file=sys.stderr)
        sys.exit(1)

    if args.debug_cols:
        print("[DEBUG] Columns (original -> normalized):")
        for c in df.columns:
            print(f"  - {c}  ->  {normalize_colname(c)}")
        return

    # Resolve key columns
    try:
        col_am_class = resolve_column(df, args.am_class_col)
    except KeyError as e:
        print(f"[ERROR] {e}", file=sys.stderr); sys.exit(1)
    try:
        col_demask = resolve_column(df, args.demask_col)
    except KeyError as e:
        print(f"[ERROR] {e}", file=sys.stderr); sys.exit(1)
    try:
        col_mut = resolve_column(df, args.mutation_col)
    except KeyError:
        col_mut = None  # optional

    # Pre-filter
    amc = df[col_am_class].astype(str).str.strip().str.casefold()
    mask_am = (amc == str(args.am_class_value).strip().casefold())

    def as_float(x):
        try:
            return float(str(x).strip())
        except Exception:
            return np.nan

    demask_vals = df[col_demask].apply(as_float)
    mask_dm = demask_vals.abs() >= float(args.demask_threshold)

    mask = mask_am & mask_dm
    pre = df[mask].copy()
    print(f"[INFO] Pre-filter kept {len(pre)} / {len(df)} rows "
          f"(AlphaMissense == '{args.am_class_value}', |{args.demask_col}| >= {args.demask_threshold}).")

    if len(pre) == 0:
        quoting = csv.QUOTE_ALL if args.quote_all else csv.QUOTE_MINIMAL
        pre.to_csv(out_path, index=False, encoding="utf-8-sig", quoting=quoting)
        print(f"[INFO] No candidates after pre-filter. Saved empty file: {out_path}")
        return

    # Parse config
    try:
        specs = parse_config_file(cfg_path)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    # Resolve feature columns
    resolved_specs: List[Tuple[FeatureSpec, str]] = []
    for spec in specs:
        try:
            col = resolve_column(pre, spec.col)
            resolved_specs.append((spec, col))
        except KeyError as e:
            print(f"[ERROR] {e}", file=sys.stderr); sys.exit(1)

    # Compute per-feature contributions
    contrib_cols: Dict[str, pd.Series] = {}
    total_score = pd.Series([0.0] * len(pre), index=pre.index, dtype=float)

    for spec, col in resolved_specs:
        safe_name = f"feat_{col}"
        if spec.kind == "NUM":
            s = apply_num_feature(pre[col], spec)  # 0..1 * weight
        else:
            s = apply_cat_feature(pre[col], spec)  # 0..1 * weight
        contrib_cols[safe_name] = s
        total_score = total_score + s.fillna(0.0)

    # Build output
    out = pre.copy()
    for name, series in contrib_cols.items():
        out[name] = series
    out["RankingScore"] = total_score

    # Rank descending by score
    out = out.sort_values(by=["RankingScore"], ascending=False)
    out["Rank"] = range(1, len(out) + 1)

    # Save
    quoting = csv.QUOTE_ALL if args.quote_all else csv.QUOTE_MINIMAL
    try:
        out.to_csv(out_path, index=False, encoding="utf-8-sig", quoting=quoting)
    except Exception as e:
        print(f"[ERROR] Failed to write output: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Wrote ranked candidates to: {out_path}")
    if args.top > 0:
        cols_to_show = [c for c in [col_mut, "RankingScore"] if c] + list(contrib_cols.keys())
        head = out[cols_to_show].head(args.top)
        print("[TOP] Preview:")
        with pd.option_context("display.max_columns", None, "display.width", 200):
            print(head.to_string(index=False))


if __name__ == "__main__":
    main()

