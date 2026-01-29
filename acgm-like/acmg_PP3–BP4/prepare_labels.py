#!/usr/bin/env python3
"""
prepare_labels_pole_pold1.py

Prepare a labeled dataset (P/LP vs B/LB) from MAVISp CSV files for POLE/POLD1,
using ClinVar and optionally POLED_DB functional data.

Examples
--------

ClinVar only, all review levels:

    python prepare_labels_pole_pold1.py \
        -i POLE-ensemble_mode.csv POLD1-ensemble_mode.csv \
        -o pole_pold1_labels_clinvar_all.csv \
        --use_clinvar

ClinVar only, require at least 2-star confidence:

    python prepare_labels_pole_pold1.py \
        -i POLE-ensemble_mode.csv POLD1-ensemble_mode.csv \
        -o pole_pold1_labels_clinvar_2plus.csv \
        --use_clinvar \
        --min_clinvar_stars 2

ClinVar (≥1 star) + POLED_DB as fallback for variants without ClinVar label:

    python prepare_labels_pole_pold1.py \
        -i POLE-ensemble_mode.csv POLD1-ensemble_mode.csv \
        -o pole_pold1_labels_clinvar_poledb.csv \
        --use_clinvar \
        --min_clinvar_stars 1 \
        --use_poledb
"""

import argparse
import os
import re
import sys
from typing import Optional

import pandas as pd


# ---------------------- Helpers ------------------------------------------------


def infer_gene_from_filename(path: str) -> str:
    """
    Infer gene name from filename, e.g.:
    'POLE-ensemble_mode.csv' -> 'POLE'
    'POLD1-simple_mode.csv'  -> 'POLD1'
    """
    base = os.path.basename(path)
    m = re.match(r"([A-Za-z0-9]+)", base)
    if m:
        return m.group(1)
    return base


def normalize(s: Optional[str]) -> str:
    """Lowercase & strip, handle NaN/None safely."""
    if pd.isna(s):
        return ""
    return str(s).strip().lower()


# ---- ClinVar mapping ---------------------------------------------------------

PATHOGENIC_TERMS = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic, likely oncogenic",   # POLE / cancer-specific
    "pathogenic/likely pathogenic",   # POLD1 form
}

BENIGN_TERMS = {
    "benign",
    "likely benign",
    "benign/likely benign",
    "likely benign, benign",          # POLD1 combination form
}


def parse_clinvar_stars(review_status_value) -> int:
    """
    Parse MAVISp 'ClinVar Review Status' which can be:
      '0', '1', '2', '1, 1', '1, 1, 2', '2, 2', etc.

    We interpret this as a comma-separated list of star ratings per submitter
    and take the maximum as the effective star level.
    """
    if pd.isna(review_status_value):
        return 0

    s = str(review_status_value).strip()
    if not s:
        return 0

    parts = [p.strip() for p in s.split(",")]
    stars = []
    for p in parts:
        if not p:
            continue
        try:
            stars.append(int(p))
        except ValueError:
            continue

    if not stars:
        return 0

    return max(stars)


def classify_clinvar_interpretation(interpretation: str) -> Optional[str]:
    """
    Map ClinVar Interpretation to coarse labels:

    - 'P/LP'  (pathogenic / likely pathogenic)
    - 'B/LB'  (benign / likely benign)
    - None    (VUS, conflicting, drug response, etc.)

    Uses exact matches based on dictionaries from POLE + POLD1 MAVISp files.
    """
    s = normalize(interpretation)

    if s in PATHOGENIC_TERMS:
        return "P/LP"
    if s in BENIGN_TERMS:
        return "B/LB"

    # Everything else (uncertain, conflicting, drug response, etc.) -> no label
    return None


# ---- POLED_DB mapping --------------------------------------------------------


def classify_poledb(exp_class: str) -> Optional[str]:
    """
    Map POLED_DB experimental classification to P/LP or B/LB.

    In your MAVISp POLE and POLD1 files, the column
      'Experimental data classification (miscellaneous, Functional (POLED_DB))'
    currently has only:

        'damaging' -> P/LP
        'neutral'  -> B/LB

    We map those directly and keep a small set of fallbacks
    in case future files add more textual descriptions.
    """
    s = normalize(exp_class)
    if not s:
        return None

    # Standard POLED_DB categories
    if s == "damaging":
        return "P/LP"
    if s == "neutral":
        return "B/LB"

    # Conservative fallback (for future-proofing)
    if any(kw in s for kw in [
        "loss of function",
        "lof",
        "deleterious",
        "reduced function",
        "reduced activity",
        "strongly decreased",
        "non-functional",
        "impair",
        "inactivate",
    ]):
        return "P/LP"

    if any(kw in s for kw in [
        "benign",
        "wt-like",
        "wildtype-like",
        "neutral effect",
        "no effect",
        "no functional effect",
        "similar to wildtype",
        "retained activity",
    ]):
        return "B/LB"

    return None


# ---------------------- Main script -------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare labeled P/LP vs B/LB dataset from MAVISp POLE/POLD1 CSVs "
            "using ClinVar and optionally POLED_DB."
        )
    )
    parser.add_argument(
        "-i", "--input", nargs="+", required=True,
        help="Input MAVISp CSV files (e.g. POLE-ensemble_mode.csv POLD1-ensemble_mode.csv)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output CSV with labels (includes all MAVISp columns plus label columns)."
    )
    parser.add_argument(
        "--use_clinvar", action="store_true",
        help="Use ClinVar Interpretation + Review Status to assign labels."
    )
    parser.add_argument(
        "--min_clinvar_stars", type=int, default=0,
        help=(
            "Minimum ClinVar review status stars to accept (0–4). "
            "Uses the maximum star across submitters in 'ClinVar Review Status'. "
            "Default: 0 (no filtering)."
        )
    )
    parser.add_argument(
        "--use_poledb", action="store_true",
        help=(
            "Use POLED_DB experimental classification as an additional label source "
            "for variants without a ClinVar label (or if ClinVar is not used)."
        )
    )
    parser.add_argument(
        "--keep_unlabeled", action="store_true",
        help=(
            "If set, keep unlabeled variants in the output. "
            "By default, only labeled variants (P/LP or B/LB) are written."
        )
    )

    args = parser.parse_args()

    if not args.use_clinvar and not args.use_poledb:
        print(
            "⚠️  You did not select any label source (--use_clinvar or --use_poledb). "
            "Nothing to do.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Load and pool input files
    frames = []
    for path in args.input:
        if not os.path.exists(path):
            print(f"❌ Input file not found: {path}", file=sys.stderr)
            sys.exit(1)
        try:
            df = pd.read_csv(path, low_memory=False)
        except Exception as e:
            print(f"❌ Failed to read {path}: {e}", file=sys.stderr)
            sys.exit(1)

        gene = infer_gene_from_filename(path)
        df["Gene"] = gene
        frames.append(df)

    if not frames:
        print("❌ No input data loaded.", file=sys.stderr)
        sys.exit(1)

    data = pd.concat(frames, ignore_index=True)

    # Check required columns
    required_cols = []
    if args.use_clinvar:
        required_cols += ["ClinVar Interpretation", "ClinVar Review Status"]
    if args.use_poledb:
        required_cols += ["Experimental data classification (miscellaneous, Functional (POLED_DB))"]

    missing = [c for c in required_cols if c not in data.columns]
    if missing:
        print("❌ Missing required columns in input data:", ", ".join(missing), file=sys.stderr)
        sys.exit(1)

    # Initialize label columns
    data["MAVISp_label"] = pd.NA        # "P/LP" or "B/LB"
    data["MAVISp_label_source"] = pd.NA # "ClinVar" or "POLED_DB"

    # ---- ClinVar-based labels ----
    if args.use_clinvar:
        data["ClinVar_stars"] = data["ClinVar Review Status"].apply(parse_clinvar_stars)
        clinvar_mask = data["ClinVar_stars"] >= args.min_clinvar_stars

        for idx, row in data[clinvar_mask].iterrows():
            label = classify_clinvar_interpretation(row["ClinVar Interpretation"])
            if label is not None:
                data.at[idx, "MAVISp_label"] = label
                data.at[idx, "MAVISp_label_source"] = "ClinVar"

    # ---- POLED_DB-based labels (fallback / additional) ----
    if args.use_poledb:
        col_poledb = "Experimental data classification (miscellaneous, Functional (POLED_DB))"

        for idx, row in data.iterrows():
            # If ClinVar already labeled it, keep ClinVar as primary source
            if pd.notna(row["MAVISp_label"]):
                continue
            label = classify_poledb(row[col_poledb])
            if label is not None:
                data.at[idx, "MAVISp_label"] = label
                data.at[idx, "MAVISp_label_source"] = "POLED_DB"

    # --- Summary BEFORE dropping unlabeled ---
    total = len(data)
    labeled_mask = data["MAVISp_label"].notna()
    labeled = labeled_mask.sum()
    num_p = (data["MAVISp_label"] == "P/LP").sum()
    num_b = (data["MAVISp_label"] == "B/LB").sum()
    num_clinvar = (data["MAVISp_label_source"] == "ClinVar").sum()
    num_poledb = (data["MAVISp_label_source"] == "POLED_DB").sum()

    print(f"📊 Total variants read: {total}")
    print(f"📊 Labeled variants: {labeled} (P/LP: {num_p}, B/LB: {num_b})")
    if args.use_clinvar:
        print(f"   - From ClinVar : {num_clinvar}")
        if "ClinVar_stars" in data.columns:
            cv_labeled = data[data["MAVISp_label_source"] == "ClinVar"]["ClinVar_stars"]
            if not cv_labeled.empty:
                print("   - ClinVar star distribution (labeled):")
                print(cv_labeled.value_counts().sort_index().to_string())
    if args.use_poledb:
        print(f"   - From POLED_DB: {num_poledb}")

    # Optionally drop unlabeled variants
    if not args.keep_unlabeled:
        data = data[labeled_mask].copy()

    # Add numeric binary label: 1 = P/LP, 0 = B/LB
    def to_binary(lbl):
        if lbl == "P/LP":
            return 1
        if lbl == "B/LB":
            return 0
        return pd.NA

    data["MAVISp_label_binary"] = data["MAVISp_label"].apply(to_binary)

    # Write output
    try:
        data.to_csv(args.output, index=False)
    except Exception as e:
        print(f"❌ Failed to write output CSV {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"✅ Wrote labeled dataset to: {args.output}")
    print("   (Includes all original MAVISp columns plus: Gene, ClinVar_stars (if used), "
          "MAVISp_label, MAVISp_label_source, MAVISp_label_binary)")


if __name__ == "__main__":
    main()

