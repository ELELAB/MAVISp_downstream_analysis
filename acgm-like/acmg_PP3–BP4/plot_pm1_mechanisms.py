#!/usr/bin/env python3
import argparse
import os
import csv

import pandas as pd
import matplotlib.pyplot as plt


def get_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter variants where PP3_BP4 == 'PP3' and PM1_like == True, "
            "annotate mechanisms of action, and visualize counts + heatmap."
        )
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input MAVISp CSV file (e.g. POLE_vus_acmg.csv)"
    )
    parser.add_argument(
        "-o", "--out-prefix",
        help="Output prefix (default: input name without extension)"
    )
    return parser.parse_args()


def normalize_series(series):
    return series.astype(str).str.strip().str.lower()


def extract_partner_from_col(col_name):
    """
    From 'Local Int. classification (PCNA_9B8T_AC) [md]'
    -> 'PCNA'
    """
    try:
        start = col_name.index("(") + 1
        end = col_name.index(")", start)
        inside = col_name[start:end]  # e.g. 'PCNA_9B8T_AC'
        return inside.split("_")[0]
    except ValueError:
        return "UNKNOWN"


def main():
    args = get_args()
    input_path = args.input
    if args.out_prefix:
        out_prefix = args.out_prefix
    else:
        base = os.path.basename(input_path)
        out_prefix = os.path.splitext(base)[0]

    # ---- Read CSV ----
    df = pd.read_csv(input_path, dtype=str)

    # ---- Required cols for filtering ----
    for col in ["PP3_BP4", "PM1_like", "Mutation"]:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' not found in input file.")

    pp3 = normalize_series(df["PP3_BP4"]) == "pp3"
    pm1_true = normalize_series(df["PM1_like"]) == "true"
    mask = pp3 & pm1_true

    df_filt = df[mask].copy()
    if df_filt.empty:
        print("No variants with PP3_BP4 == 'PP3' and PM1_like == True.")
        return

    print(f"Filtered variants: {df_filt.shape[0]}")

    # ---- Collect column groups for mechanisms ----

    # Stability classification _md_25
    stab_col = None
    for c in df.columns:
        if "Stability classification" in c and "RaSP" in c and "[md_25]" in c:
            stab_col = c
            break
    if stab_col is None:
        raise ValueError("Could not find 'Stability classification, (RaSP, FoldX) [md_25]' column.")

    # Local Int. With DNA classification
    dna_class_cols = [
        c for c in df.columns
        if "Local Int. With DNA classification" in c
    ]

    # Local Int. (protein–protein) classification (no DNA)
    local_class_cols = [
        c for c in df.columns
        if "Local Int. classification (" in c and "With DNA" not in c
    ]

    # PTM / Long-range / functional site columns
    ptm_stab_col = "PTM effect in stability [md]"
    ptm_func_col = "PTM effect in function [md]"
    allo_col = "AlloSigma2-PSN classification - pockets and interfaces [md]"
    func_cofactor_col = "Functional sites (cofactor)"
    func_active_col = "Functional sites (active site)"

    for col in [ptm_stab_col, ptm_func_col, allo_col,
                func_cofactor_col, func_active_col]:
        if col not in df.columns:
            print(f"Warning: expected column '{col}' not found. It will be treated as empty.")

    # ---- Mechanism annotation ----

    broad_mechs = ["Stability",
                   "Local Int.",
                   "Local Int. DNA",
                   "PTM stability",
                   "PTM function",
                   "Long Range",
                   "Functional sites"]

    mech_flags = {m: [] for m in broad_mechs}
    local_partners_list = []  # e.g. 'PCNA;MCM2'

    for _, row in df_filt.iterrows():
        row_mechs = set()
        local_partners = set()

        # Stability
        val_stab = str(row.get(stab_col, "")).strip().lower()
        if val_stab == "destabilizing":
            row_mechs.add("Stability")

        # Local Int. with DNA
        has_local_dna = False
        for c in dna_class_cols:
            v = str(row.get(c, "")).strip().lower()
            if v in ("destabilizing", "stabilizing"):
                has_local_dna = True
                break
        if has_local_dna:
            row_mechs.add("Local Int. DNA")

        # Local Int. protein–protein
        has_local_int = False
        for c in local_class_cols:
            v = str(row.get(c, "")).strip().lower()
            if v == "destabilizing":
                has_local_int = True
                partner = extract_partner_from_col(c)
                local_partners.add(partner)
        if has_local_int:
            row_mechs.add("Local Int.")

        # PTM stability
        val_ptm_stab = str(row.get(ptm_stab_col, "")).strip().lower()
        if val_ptm_stab == "damaging":
            row_mechs.add("PTM stability")

        # PTM function
        val_ptm_func = str(row.get(ptm_func_col, "")).strip().lower()
        if val_ptm_func == "damaging":
            row_mechs.add("PTM function")

        # Long Range
        val_allo = str(row.get(allo_col, "")).strip().lower()
        if val_allo == "damaging":
            row_mechs.add("Long Range")

        # Functional sites
        val_cof = str(row.get(func_cofactor_col, "")).strip().lower()
        val_act = str(row.get(func_active_col, "")).strip().lower()
        if val_cof == "damaging" or val_act == "damaging":
            row_mechs.add("Functional sites")

        # Store flags
        for m in broad_mechs:
            mech_flags[m].append(1 if m in row_mechs else 0)

        if local_partners:
            local_partners_list.append(";".join(sorted(local_partners)))
        else:
            local_partners_list.append("")

    # Add mechanism columns to df_filt
    for m in broad_mechs:
        colname = f"Mechanism_{m.replace(' ', '_')}"
        df_filt[colname] = mech_flags[m]

    df_filt["Local_Int_partners"] = local_partners_list

    # Also add a compact list of mechanisms per variant
    mech_list_compact = []
    for i in range(df_filt.shape[0]):
        this_list = [m for m in broad_mechs if mech_flags[m][i] == 1]
        mech_list_compact.append(";".join(this_list))
    df_filt["Mechanisms_list"] = mech_list_compact

    # ---- Save filtered table ----
    filtered_csv = f"{out_prefix}_PP3_PM1like_mechanisms.csv"
    df_filt.to_csv(filtered_csv, index=False, quoting=csv.QUOTE_ALL)
    print(f"Saved filtered variants with mechanisms to: {filtered_csv}")

    # ---- Summary counts per mechanism ----
    counts = {m: sum(mech_flags[m]) for m in broad_mechs}
    summary_df = pd.DataFrame({
        "Mechanism": list(counts.keys()),
        "Count": list(counts.values())
    }).sort_values("Count", ascending=False)

    summary_csv = f"{out_prefix}_PP3_PM1like_mechanism_counts.csv"
    summary_df.to_csv(summary_csv, index=False, quoting=csv.QUOTE_ALL)
    print(f"Saved mechanism counts to: {summary_csv}")

    # ---- Bar plot ----
    plt.figure(figsize=(6, 4))
    plt.bar(summary_df["Mechanism"], summary_df["Count"])
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Number of variants")
    plt.title("Mechanisms for variants with PP3 & PM1_like")
    plt.tight_layout()

    png_bar = f"{out_prefix}_PP3_PM1like_mechanism_counts.png"
    pdf_bar = f"{out_prefix}_PP3_PM1like_mechanism_counts.pdf"
    plt.savefig(png_bar, dpi=300)
    plt.savefig(pdf_bar)
    plt.close()
    print(f"Saved bar plots:\n  {png_bar}\n  {pdf_bar}")

    # ---- Heatmap: variants x mechanisms ----
    # Build matrix: rows = variants, cols = mechanisms
    mech_matrix = pd.DataFrame(
        {m: mech_flags[m] for m in broad_mechs},
        index=df_filt["Mutation"]
    )

    # Optionally sort rows by total number of mechanisms (descending)
    mech_matrix["_sum"] = mech_matrix.sum(axis=1)
    mech_matrix = mech_matrix.sort_values("_sum", ascending=False)
    mech_matrix = mech_matrix.drop(columns=["_sum"])

    plt.figure(figsize=(8, max(4, 0.25 * mech_matrix.shape[0])))
    plt.imshow(mech_matrix.values, aspect="auto")
    plt.colorbar(label="Mechanism present (1) / absent (0)")
    plt.xticks(range(len(broad_mechs)), broad_mechs, rotation=45, ha="right")

    # Show mutation labels on y-axis only if not crazy many
    if mech_matrix.shape[0] <= 60:
        plt.yticks(range(mech_matrix.shape[0]), mech_matrix.index)
    else:
        plt.yticks([])

    plt.xlabel("Mechanism")
    plt.ylabel("Mutation")
    plt.title("Mechanism heatmap for variants with PP3 & PM1_like")
    plt.tight_layout()

    png_heat = f"{out_prefix}_PP3_PM1like_mechanism_heatmap.png"
    pdf_heat = f"{out_prefix}_PP3_PM1like_mechanism_heatmap.pdf"
    plt.savefig(png_heat, dpi=300)
    plt.savefig(pdf_heat)
    plt.close()
    print(f"Saved heatmap:\n  {png_heat}\n  {pdf_heat}")


if __name__ == "__main__":
    main()

