#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import csv

def main():
    parser = argparse.ArgumentParser(
        description="Summarize PP3/BP4/PM1_like and output per-category CSVs."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input MAVISp CSV file"
    )
    parser.add_argument(
        "-o", "--out-prefix",
        help="Output prefix (default: input filename without extension)"
    )

    args = parser.parse_args()

    input_path = args.input
    if args.out_prefix:
        out_prefix = args.out_prefix
    else:
        base = os.path.basename(input_path)
        out_prefix = os.path.splitext(base)[0]

    # Read CSV
    df = pd.read_csv(input_path, dtype=str)

    # Check required columns
    required_cols = ["Mutation", "PP3_BP4", "PM1_like"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in input: {', '.join(missing)}")

    # Normalize strings for safe comparison
    pp3_bp4_series = df["PP3_BP4"].astype(str).str.strip().str.upper()
    pm1_series = df["PM1_like"].astype(str).str.strip().str.upper()

    # Masks
    mask_pp3 = pp3_bp4_series == "PP3"
    mask_bp4 = pp3_bp4_series == "BP4"
    mask_pm1_true = pm1_series == "TRUE"
    mask_pm1_false = pm1_series == "FALSE"

    # Counts
    n_pp3 = mask_pp3.sum()
    n_bp4 = mask_bp4.sum()
    n_pm1_true = mask_pm1_true.sum()
    n_pm1_false = mask_pm1_false.sum()

    print("===== SUMMARY =====")
    print(f"PP3 variants (PP3_BP4 == 'PP3'):      {n_pp3}")
    print(f"BP4 variants (PP3_BP4 == 'BP4'):      {n_bp4}")
    print(f"PM1_like TRUE (PM1_like == True):     {n_pm1_true}")
    print(f"PM1_like FALSE (PM1_like == False):   {n_pm1_false}")

    # Prepare minimal columns
    cols_to_keep = ["Mutation", "PP3_BP4", "PM1_like"]

    df_pp3 = df[mask_pp3][cols_to_keep]
    df_bp4 = df[mask_bp4][cols_to_keep]
    df_pm1_true = df[mask_pm1_true][cols_to_keep]

    # Write CSVs, all values quoted
    df_pp3.to_csv(f"{out_prefix}_PP3.csv", index=False,
                  quoting=csv.QUOTE_ALL)
    df_bp4.to_csv(f"{out_prefix}_BP4.csv", index=False,
                  quoting=csv.QUOTE_ALL)
    df_pm1_true.to_csv(f"{out_prefix}_PM1_like_true.csv", index=False,
                       quoting=csv.QUOTE_ALL)

    # Optional: also write a small summary table
    summary_df = pd.DataFrame({
        "Category": ["PP3", "BP4", "PM1_like_TRUE", "PM1_like_FALSE"],
        "Count": [n_pp3, n_bp4, n_pm1_true, n_pm1_false]
    })
    summary_df.to_csv(f"{out_prefix}_ACMG_summary_counts.csv",
                      index=False, quoting=csv.QUOTE_ALL)

    print(f"\nSaved:")
    print(f"  {out_prefix}_PP3.csv")
    print(f"  {out_prefix}_BP4.csv")
    print(f"  {out_prefix}_PM1_like_true.csv")
    print(f"  {out_prefix}_ACMG_summary_counts.csv")

if __name__ == "__main__":
    main()

