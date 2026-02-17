#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------- Utilities ----------
def normalize_colname(s: str) -> str:
    """Lowercase, replace NBSP/ZWSP, collapse whitespace/newlines for robust column matching."""
    if s is None:
        return ""
    s = str(s).replace("\u00A0", " ").replace("\u200B", "")
    s = re.sub(r"\s+", " ", s)
    return s.strip().lower()


def resolve_column(df: pd.DataFrame, target_name: str) -> str:
    """Resolve a column name robustly against MAVISp-style headers."""
    if target_name in df.columns:
        return target_name
    norm_map = {normalize_colname(c): c for c in df.columns}
    key = normalize_colname(target_name)
    if key in norm_map:
        return norm_map[key]
    raise KeyError(f"Column not found: '{target_name}'. Available:\n  - " + "\n  - ".join(df.columns))


def parse_config_for_order(cfg_path: Optional[Path]) -> List[str]:
    """
    Read a ranking config (the same .txt you used) and return the feature column
    names in the order they appear. This is purely for ordering/selection.
    """
    if cfg_path is None:
        return []
    cols: List[str] = []
    for i, raw in enumerate(cfg_path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        # Accept "NUM: <col> ; ..." or "CAT: <col> ; ..."
        m = re.match(r"^(NUM|CAT)\s*:\s*(.+?)(;|$)", line, flags=re.IGNORECASE)
        if not m:
            # ignore unrelated lines
            continue
        cols.append(m.group(2).strip())
    return cols


def split_list(s: Optional[str]) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()] if s else []


def to_float(series: pd.Series) -> pd.Series:
    def cast(x):
        try:
            return float(str(x).strip())
        except Exception:
            return np.nan
    return series.apply(cast)


def minmax_0_1(arr: np.ndarray) -> np.ndarray:
    mask = np.isfinite(arr)
    if not mask.any():
        return np.zeros_like(arr)
    vmin = np.nanmin(arr[mask])
    vmax = np.nanmax(arr[mask])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax == vmin:
        out = np.zeros_like(arr)
        out[mask] = 0.5  # neutral if constant
        return out
    out = (arr - vmin) / (vmax - vmin)
    out[~mask] = np.nan
    return out


# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Plot a heatmap of the top-N ranked variants from the ranked CSV.\n"
            "Features on X-axis, variants on Y-axis, color-coded with a bar."
        )
    )
    p.add_argument("-i", "--input", required=True, help="Ranked CSV produced by rank_candidates_for_validation.py")
    p.add_argument("-o", "--output", required=True, help="Output image (e.g., .png, .pdf, .svg)")
    p.add_argument("--top", type=int, default=30, help="Number of top variants to plot (default: 30)")
    p.add_argument("--mutation-col", default="Mutation", help="Column used to label variants (default: 'Mutation')")
    p.add_argument("--rank-col", default="Rank", help="Rank column (default: 'Rank'); set empty '' to hide")
    p.add_argument("--score-col", default="RankingScore", help="Total score column (default: 'RankingScore'); set empty '' to hide")

    # Feature selection/order
    g = p.add_mutually_exclusive_group(required=False)
    g.add_argument("--features", default=None,
                   help="Comma-separated list of feature columns to show (either 'feat_*' columns or raw).")
    g.add_argument("--features-from-config", default=None,
                   help="Path to the ranking config .txt to use feature order from (will look for matching 'feat_<col>' first).")

    # Which data to plot
    p.add_argument("--use-contribs", action="store_true",
                   help="Prefer 'feat_<col>' contribution columns if present. (Recommended: matches ranking inputs.)")
    p.add_argument("--normalize-raw", choices=["minmax", "none"], default="minmax",
                   help="Normalization for raw columns when not using contributions (default: minmax).")

    # Figure tuning
    p.add_argument("--figsize", default=None,
                   help='Figure size "W,H" in inches (e.g., 10,8). Defaults to auto based on N and K.')
    p.add_argument("--dpi", type=int, default=200, help="Output DPI (default: 200)")
    p.add_argument("--cmap", default="viridis", help="Matplotlib colormap (default: viridis)")
    p.add_argument("--annotate", action="store_true", help="Overlay numeric values (rounded) on each cell")
    p.add_argument("--na-color", default="#f0f0f0", help="Color for missing values (default: #f0f0f0)")
    p.add_argument("--title", default="Top ranked variants – feature heatmap", help="Plot title")

    # Debug
    p.add_argument("--debug-cols", action="store_true", help="Print columns (original -> normalized) and exit.")
    return p.parse_args()


# ---------- Main ----------
def main():
    args = parse_args()

    df = pd.read_csv(args.input, dtype=str, keep_default_na=False, na_values=[], encoding="utf-8-sig", engine="python")

    if args.debug_cols:
        print("[DEBUG] Columns (original -> normalized):")
        for c in df.columns:
            print(f"  - {c}  ->  {normalize_colname(c)}")
        return

    # Resolve label columns
    try:
        col_mut = resolve_column(df, args.mutation_col)
    except KeyError as e:
        print(f"[ERROR] {e}", file=sys.stderr); sys.exit(1)

    col_rank = None
    if args.rank_col:
        try:
            col_rank = resolve_column(df, args.rank_col)
        except KeyError:
            col_rank = None

    col_score = None
    if args.score_col:
        try:
            col_score = resolve_column(df, args.score_col)
        except KeyError:
            col_score = None

    # Make sure RankingScore is numeric for sorting
    if col_score:
        df["_score_float"] = to_float(df[col_score])
    else:
        df["_score_float"] = 0.0

    # Pick top N by RankingScore (desc); if no score col, keep first N
    if col_score:
        df_sorted = df.sort_values(by="_score_float", ascending=False)
    else:
        df_sorted = df.copy()

    topN = max(1, int(args.top))
    dft = df_sorted.head(topN).copy()

    # Decide feature columns to plot
    requested_features: List[str] = []
    if args.features:
        requested_features = split_list(args.features)
    elif args.features_from_config:
        requested_features = parse_config_for_order(Path(args.features_from_config))
    else:
        # fallback: infer all "feat_" contribution columns in file, preserving input order
        requested_features = [c for c in dft.columns if c.startswith("feat_")]

    if not requested_features:
        print("[ERROR] No features selected or found. "
              "Pass --features, or --features-from-config, or ensure 'feat_*' columns exist.", file=sys.stderr)
        sys.exit(1)

    # Build list of actual columns to plot, preferring contributions if requested
    feat_cols: List[Tuple[str, str]] = []  # (display_name, column_name_in_df)
    for name in requested_features:
        # If user gave raw config names but we're using contributions, try 'feat_<resolved>'
        display = name
        try:
            col_raw = resolve_column(dft, name)
        except KeyError:
            col_raw = None

        col_feat = None
        if args.use_contribs:
            # contributions are named feat_<original_header>
            if col_raw:
                candidate = f"feat_{col_raw}"
            else:
                candidate = f"feat_{name}"  # best effort
            if candidate in dft.columns:
                col_feat = candidate

        if col_feat:
            feat_cols.append((display, col_feat))
        elif col_raw:
            feat_cols.append((display, col_raw))
        else:
            print(f"[WARN] Feature column not found (raw or feat_): {name}")

    if not feat_cols:
        print("[ERROR] None of the requested features are present.", file=sys.stderr)
        sys.exit(1)

    # Build matrix
    rows = []
    ylabels = []
    for _, row in dft.iterrows():
        # Label: Rank. Mutation (Score)
        parts = []
        if col_rank:
            parts.append(str(row[col_rank]).strip())
        if col_mut:
            parts.append(str(row[col_mut]).strip())
        label = ". ".join(parts) if parts else ""
        if col_score:
            try:
                s = float(row[col_score])
                label = f"{label} ({s:.3f})" if label else f"{s:.3f}"
            except Exception:
                pass
        ylabels.append(label)

        vals = []
        for disp, col in feat_cols:
            v = row.get(col, "")
            try:
                vals.append(float(v))
            except Exception:
                vals.append(np.nan)
        rows.append(vals)

    M = np.array(rows, dtype=float)  # shape: (N, K)

    # If plotting raw columns (not feat_), normalize columns to 0..1 for a fair heatmap
    plotting_contribs = all(col.startswith("feat_") for _, col in feat_cols)
    if not plotting_contribs and args.normalize_raw == "minmax":
        for j in range(M.shape[1]):
            M[:, j] = minmax_0_1(M[:, j])

    # Figure sizing
    K = M.shape[1]
    if args.figsize:
        try:
            W, H = [float(x.strip()) for x in args.figsize.split(",")]
            figsize = (W, H)
        except Exception:
            print("[WARN] Invalid --figsize; using auto.", file=sys.stderr)
            figsize = None
    else:
        # heuristic: width scales with K, height with N
        figsize = (max(8, 0.6 * K + 3), max(4, 0.4 * topN + 2))

    # Plot
    plt.rcParams["figure.dpi"] = args.dpi
    fig, ax = plt.subplots(figsize=figsize)
    cmap = plt.get_cmap(args.cmap)
    cmap = cmap.copy()
    try:
        cmap.set_bad(color=args.na_color)
    except Exception:
        pass

    im = ax.imshow(M, aspect="auto", interpolation="nearest", cmap=cmap)

    # Axes ticks/labels
    ax.set_yticks(np.arange(topN))
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_xticks(np.arange(K))
    ax.set_xticklabels([disp for disp, _ in feat_cols], rotation=45, ha="right", fontsize=9)

    ax.set_title(args.title, pad=12, fontsize=12)
    ax.set_xlabel("Features")
    ax.set_ylabel("Top variants")

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel("Scaled value" if not plotting_contribs else "Weighted contribution", rotation=270, labelpad=12)

    # Optional annotations
    if args.annotate:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7)

    fig.tight_layout()

    # Save
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved heatmap to: {out_path}")


if __name__ == "__main__":
    main()

