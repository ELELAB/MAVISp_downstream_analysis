#!/usr/bin/env python3
"""
prepare_features_from_yaml.py

Use a YAML config (feature_categories.yaml) to define ACMG/ClinGen
computational categories and build a feature table from the labeled
MAVISp dataset (POLE + POLD1).

This script:
  - loads a labeled CSV (with MAVISp_label / MAVISp_label_binary)
  - loads feature_categories.yaml (your ontology)
  - selects the columns for each category
  - engineers summary interaction DDG features:
      * max_local_int_ddg_protein  (max |ΔΔG| across protein–protein interfaces)
      * max_local_int_ddg_dna      (max |ΔΔG| across DNA interfaces)
  - writes an output CSV with:
      identifiers + all selected columns + engineered features + labels

Example:

  python prepare_features_from_yaml.py \
      -i pole_pold1_labels_clinvar.csv \
      -c config_features.yaml \
      -o pole_pold1_features_for_calibration.csv
"""

import argparse
import sys
import pandas as pd

try:
    import yaml
except ImportError:
    yaml = None


def load_yaml_config(path):
    if yaml is None:
        raise ImportError(
            "PyYAML is not installed. Install it via:\n\n"
            "   pip install pyyaml\n"
        )
    with open(path, "r") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare feature table from labeled MAVISp dataset using feature_categories.yaml."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input labeled CSV (output of prepare_labels_pole_pold1.py)."
    )
    parser.add_argument(
        "-c", "--config", required=True,
        help="YAML config file defining categories and column mapping (feature_categories.yaml)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output CSV with features + labels."
    )
    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # Load input data
    # -------------------------------------------------------------------------
    try:
        df = pd.read_csv(args.input, low_memory=False)
    except Exception as e:
        print(f"❌ Could not read input file {args.input}: {e}", file=sys.stderr)
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Load YAML config
    # -------------------------------------------------------------------------
    try:
        config = load_yaml_config(args.config)
    except Exception as e:
        print(f"❌ Could not read YAML config {args.config}: {e}", file=sys.stderr)
        sys.exit(1)

    categories = config.get("categories", {})
    options = config.get("options", {})
    ignore_missing = bool(options.get("ignore_missing_columns", True))
    verbose = bool(options.get("verbose", False))

    # -------------------------------------------------------------------------
    # 1) Figure out which columns we want from each category
    # -------------------------------------------------------------------------
    selected_cols_by_category = {}
    all_columns_from_categories = set()

    for cat_name, cat_info in categories.items():
        cols = cat_info.get("columns", []) or []
        existing = []
        missing = []

        for c in cols:
            if c in df.columns:
                existing.append(c)
                all_columns_from_categories.add(c)
            else:
                missing.append(c)

        selected_cols_by_category[cat_name] = existing

        if missing and not ignore_missing:
            print(
                f"❌ Missing columns for category '{cat_name}': {', '.join(missing)}",
                file=sys.stderr,
            )
            sys.exit(1)
        elif missing and verbose and missing:
            print(
                f"⚠️  Category '{cat_name}': ignoring missing columns: {', '.join(missing)}",
                file=sys.stderr,
            )

        if verbose:
            print(f"✅ Category '{cat_name}' will use columns:")
            for c in existing:
                print(f"   - {c}")

    # -------------------------------------------------------------------------
    # 2) Engineer summary interaction DDG features
    #
    # We treat interaction_functional columns as:
    #   - protein–protein interfaces:
    #       "Local Int. (Binding with" AND NOT "With DNA"
    #   - DNA interfaces:
    #       "Local Int. With DNA"
    #
    # For each row we compute:
    #   max_local_int_ddg_protein = max(|ΔΔG| over all protein–protein columns)
    #   max_local_int_ddg_dna     = max(|ΔΔG| over all DNA columns)
    # -------------------------------------------------------------------------
    interaction_cols = selected_cols_by_category.get("interaction_functional", [])

    # Protein–protein ΔΔG columns (any partner, FoldX or Rosetta), not DNA
    pp_ddg_cols = [
        c for c in interaction_cols
        if c.startswith("Local Int. (Binding with") and "With DNA" not in c
    ]

    # DNA ΔΔG columns
    dna_ddg_cols = [
        c for c in interaction_cols
        if c.startswith("Local Int. With DNA")
    ]

    if verbose:
        print("Interaction summary features:")
        print(f"  Protein–protein ΔΔG columns: {len(pp_ddg_cols)}")
        for c in pp_ddg_cols:
            print(f"    - {c}")
        print(f"  DNA ΔΔG columns: {len(dna_ddg_cols)}")
        for c in dna_ddg_cols:
            print(f"    - {c}")

    # Compute max |ΔΔG| over each group
    engineered_cols = []

    if pp_ddg_cols:
        # Convert subset to numeric safely
        pp_df = df[pp_ddg_cols].apply(pd.to_numeric, errors="coerce")
        df["max_local_int_ddg_protein"] = pp_df.abs().max(axis=1)
        engineered_cols.append("max_local_int_ddg_protein")
    else:
        # no pp cols selected in config or present in df
        df["max_local_int_ddg_protein"] = pd.NA
        engineered_cols.append("max_local_int_ddg_protein")

    if dna_ddg_cols:
        dna_df = df[dna_ddg_cols].apply(pd.to_numeric, errors="coerce")
        df["max_local_int_ddg_dna"] = dna_df.abs().max(axis=1)
        engineered_cols.append("max_local_int_ddg_dna")
    else:
        df["max_local_int_ddg_dna"] = pd.NA
        engineered_cols.append("max_local_int_ddg_dna")

    # -------------------------------------------------------------------------
    # 3) Build final feature table
    # -------------------------------------------------------------------------
    # Identifier + label columns that we want to keep if present
    id_cols = ["Gene", "Mutation", "HGVSp"]
    label_cols = ["MAVISp_label", "MAVISp_label_binary"]

    final_cols = []
    for col in id_cols + list(all_columns_from_categories) + engineered_cols + label_cols:
        if col in df.columns and col not in final_cols:
            final_cols.append(col)

    feature_df = df[final_cols].copy()

    # -------------------------------------------------------------------------
    # 4) Save output
    # -------------------------------------------------------------------------
    try:
        feature_df.to_csv(args.output, index=False)
    except Exception as e:
        print(f"❌ Failed to write output file {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

    # -------------------------------------------------------------------------
    # 5) Report summary
    # -------------------------------------------------------------------------
    n_rows, n_cols = feature_df.shape
    print("==== FEATURE TABLE SUMMARY ====")
    print(f"Input labeled variants: {len(df)}")
    print(f"Output feature rows:    {n_rows}")
    print(f"Output feature columns: {n_cols}")
    print(f"Categories used:        {', '.join(categories.keys())}")
    print(f"Engineered features:    {', '.join(engineered_cols)}")
    print(f"Output written to:      {args.output}")
    print("================================")


if __name__ == "__main__":
    main()

