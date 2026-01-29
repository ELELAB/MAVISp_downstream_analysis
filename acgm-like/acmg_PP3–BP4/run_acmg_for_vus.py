#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import yaml

PP3_FEATURES = [
    "REVEL score",
    "AlphaMissense pathogenicity score",
    "EVE score",
    "GEMME Score",
    "DeMaSk delta fitness"
]

LOWER_IS_WORSE = {"GEMME Score", "DeMaSk delta fitness"}

def load_thresholds(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)["thresholds"]

def is_vus(clinvar_interp):
    if isinstance(clinvar_interp, str):
        c = clinvar_interp.lower()
        return (
            "uncertain" in c or
            c == "vus" or
            "variant of uncertain significance" in c
        )
    return False

def assign_pp3_bp4(row, thr):
    pp3_votes = 0
    bp4_votes = 0

    for feat in PP3_FEATURES:
        if feat not in row or feat not in thr:
            continue

        x = row[feat]
        if pd.isna(x):
            continue

        p25 = thr[feat]["percentile_pathogenic_25th"]
        b75 = thr[feat]["percentile_benign_75th"]
        invert = feat in LOWER_IS_WORSE

        if invert:
            # lower is worse
            if x <= p25:
                pp3_votes += 1
            if x >= b75:
                bp4_votes += 1
        else:
            # higher is worse
            if x >= p25:
                pp3_votes += 1
            if x <= b75:
                bp4_votes += 1

    if pp3_votes >= 2 and bp4_votes == 0:
        return "PP3_supporting"
    elif bp4_votes >= 2 and pp3_votes == 0:
        return "BP4_supporting"
    else:
        return "None"

def assign_pm1_like(row):
    # PM1: strong critical functional sites
    if str(row.get("Functional sites (active site)", "")).lower() == "yes":
        return "PM1"
    if str(row.get("Functional sites (cofactor)", "")).lower() == "yes":
        return "PM1"

    allo = str(row.get("AlloSigma2-PSN classification - active sites [md]", "")).lower()
    if allo in {"strong", "critical", "direct"}:
        return "PM1"

    # PM1-supporting
    hits = 0

    if "max_local_int_ddg_protein" in row and pd.notna(row["max_local_int_ddg_protein"]):
        if abs(row["max_local_int_ddg_protein"]) > 1.0:
            hits += 1

    if "max_local_int_ddg_dna" in row and pd.notna(row["max_local_int_ddg_dna"]):
        if abs(row["max_local_int_ddg_dna"]) > 1.0:
            hits += 1

    allo2 = str(row.get("AlloSigma2-PSN classification - pockets and interfaces [md]", "")).lower()
    if allo2 in {"strong", "impact", "moderate"}:
        hits += 1

    if str(row.get("is site part of phospho-SLiM [md]", "")).lower() == "yes":
        hits += 1

    if str(row.get("PTM effect in function [md]", "")).lower() in {"impact", "strong"}:
        hits += 1

    if hits >= 1:
        return "PM1_supporting"

    return "None"

def main():
    ap = argparse.ArgumentParser(description="ACMG PP3/BP4 and PM1-like annotation for ClinVar VUS.")
    ap.add_argument("-i", "--input", required=True, help="MAVISp CSV file (POLE+POLD1)")
    ap.add_argument("-t", "--thresholds", required=True, help="calibrated_thresholds.yaml")
    ap.add_argument("-o", "--output", required=True, help="Output CSV")

    args = ap.parse_args()

    df = pd.read_csv(args.input)
    thr = load_thresholds(args.thresholds)

    # Select only VUS
    vus_df = df[df["ClinVar Interpretation"].apply(is_vus)].copy()

    print(f"Total variants in file: {len(df)}")
    print(f"Total ClinVar VUS: {len(vus_df)}")

    vus_df["PP3_BP4"] = vus_df.apply(lambda r: assign_pp3_bp4(r, thr), axis=1)
    vus_df["PM1_like"] = vus_df.apply(assign_pm1_like, axis=1)

    # Optional interpretive flag
    def combine(r):
        if r["PP3_BP4"] != "None" and r["PM1_like"] != "None":
            return "PP3 + PM1_supporting"
        if r["PP3_BP4"] != "None":
            return r["PP3_BP4"]
        if r["PM1_like"] != "None":
            return r["PM1_like"]
        return "No_evidence"
    vus_df["Potential_evidence_gain"] = vus_df.apply(combine, axis=1)

    vus_df.to_csv(args.output, index=False)
    print(f"Results written to {args.output}")

if __name__ == "__main__":
    main()

