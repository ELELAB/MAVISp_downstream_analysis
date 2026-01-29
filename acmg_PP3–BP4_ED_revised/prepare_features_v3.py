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
  - engineers combined stability feature:
      * stability_combined_ddg = mean of FoldX5 & RaSP ΔΔG (md_25) per variant
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
    #   Protein–protein interfaces (any partner, FoldX or Rosetta), not DNA:
    #       "Local Int. (Binding with" AND NOT "With DNA"
    #
    #   DNA interfaces:
    #       "Local Int. With DNA"
    # -------------------------------------------------------------------------
    interaction_cols = selected_cols_by_category.get("interaction_functional", [])

    pp_ddg_cols = [
        c for c in interaction_cols
        if c.startswith("Local Int. (Binding with") and "With DNA" not in c
    ]

    dna_ddg_cols = [
        c for c in interaction_cols
        if c.startswith("Local Int. With DNA")
    ]

    engineered_cols = []

    if verbose:
        print("Interaction summary features:")
        print(f"  Protein–protein ΔΔG columns: {len(pp_ddg_cols)}")
        for c in pp_ddg_cols:
            print(f"    - {c}")
        print(f"  DNA ΔΔG columns: {len(dna_ddg_cols)}")
        for c in dna_ddg_cols:
            print(f"    - {c}")

    if pp_ddg_cols:
        pp_df = df[pp_ddg_cols].apply(pd.to_numeric, errors="coerce")
        df["max_local_int_ddg_protein"] = pp_df.abs().max(axis=1)
        engineered_cols.append("max_local_int_ddg_protein")
    else:
        df["max_local_int_ddg_protein"] = pd.NA
        engineered_cols.append("max_local_int_ddg_protein")

    if dna_ddg_cols:
        dna_df = df[dna_ddg_cols].apply(pd.to_numeric, errors="coerce")

        # Legacy absolute summary across all DNA interfaces (kept for backward compatibility):
        #   max_local_int_ddg_dna = max |ΔΔG| across all DNA-bound structures.
        df["max_local_int_ddg_dna"] = dna_df.abs().max(axis=1)
        engineered_cols.append("max_local_int_ddg_dna")

        # New: sign-aware summaries for DNA interactions
        #
        #   • max_local_int_ddg_dna_pos — strongest destabilizing effect
        #       (max positive ΔΔG across all DNA interfaces)
        #
        #   • max_local_int_ddg_dna_neg — strongest stabilizing effect
        #       (max |negative ΔΔG| across all DNA interfaces, reported as a positive magnitude)
        #
        dna_pos = dna_df.where(dna_df > 0)
        dna_neg = dna_df.where(dna_df < 0)

        # Strongest destabilization (if any positives; otherwise NaN)
        df["max_local_int_ddg_dna_pos"] = dna_pos.max(axis=1)

        # Strongest stabilization, reported as |ΔΔG| of the most negative value
        df["max_local_int_ddg_dna_neg"] = (-dna_neg).max(axis=1)

        engineered_cols.append("max_local_int_ddg_dna_pos")
        engineered_cols.append("max_local_int_ddg_dna_neg")
    else:
        df["max_local_int_ddg_dna"] = pd.NA
        df["max_local_int_ddg_dna_pos"] = pd.NA
        df["max_local_int_ddg_dna_neg"] = pd.NA
        engineered_cols.append("max_local_int_ddg_dna")
        engineered_cols.append("max_local_int_ddg_dna_pos")
        engineered_cols.append("max_local_int_ddg_dna_neg")

    # -------------------------------------------------------------------------
    # 3) Engineer combined stability feature (FoldX + RaSP)
    #
    # stability_combined_ddg = mean of available:
    #   - Stability (FoldX5, kcal/mol) [md_25]
    #   - Stability (RaSP, kcal/mol) [md_25]
    # per variant (row-wise), ignoring NaNs.
    # -------------------------------------------------------------------------
    stability_cols = []
    foldx_col = "Stability (FoldX5, kcal/mol) [md_25]"
    rasp_col = "Stability (RaSP, kcal/mol) [md_25]"

    for c in (foldx_col, rasp_col):
        if c in df.columns:
            stability_cols.append(c)

    if stability_cols:
        stab_df = df[stability_cols].apply(pd.to_numeric, errors="coerce")
        df["stability_combined_ddg"] = stab_df.mean(axis=1)
        engineered_cols.append("stability_combined_ddg")
    else:
        df["stability_combined_ddg"] = pd.NA
        engineered_cols.append("stability_combined_ddg")

    # -------------------------------------------------------------------------
    # 4) Build final feature table
    # -------------------------------------------------------------------------
    id_cols = ["Gene", "Mutation", "HGVSp"]
    label_cols = ["MAVISp_label", "MAVISp_label_binary"]

    final_cols = []
    for col in id_cols + list(all_columns_from_categories) + engineered_cols + label_cols:
        if col in df.columns and col not in final_cols:
            final_cols.append(col)

    feature_df = df[final_cols].copy()

    # -------------------------------------------------------------------------
    # 5) Save output
    # -------------------------------------------------------------------------
    try:
        feature_df.to_csv(args.output, index=False)
    except Exception as e:
        print(f"❌ Failed to write output file {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

    # -------------------------------------------------------------------------
    # 6) Report summary
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

