#!/usr/bin/env python3
import argparse
import pandas as pd
import yaml
import numpy as np

###############################################################################
# Helper functions
###############################################################################

def load_thresholds(yaml_path):
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)
    return data.get("thresholds", {})

def is_vus(clinvar_interp):
    """
    Define VUS using ClinVar interpretation text.
    """
    if pd.isna(clinvar_interp):
        return False
    text = str(clinvar_interp).lower()
    return any([
        "uncertain" in text,
        "vus" in text,
        "conflicting" in text,
        "association" in text,
        "risk factor" in text,
        "not provided" in text
    ])

def pp3_bp4_call(value, feature_name, thr_info):
    """
    PP3/BP4 assignment for a single predictor.
    """
    if pd.isna(value):
        return None

    pos25 = thr_info.get("percentile_pathogenic_25th", None)
    neg75 = thr_info.get("percentile_benign_75th", None)

    if pos25 is None or neg75 is None:
        return None

    # Predictors where lower = more pathogenic
    negative_direction = {"GEMME Score", "DeMaSk delta fitness"}

    if feature_name in negative_direction:
        if value <= pos25:
            return "PP3"
        elif value >= neg75:
            return "BP4"
    else:
        if value >= pos25:
            return "PP3"
        elif value <= neg75:
            return "BP4"

    return None


###############################################################################
# PM1-like logic
###############################################################################

def pm1_like(row, pm1_rules):
    """
    PM1-like TRUE if ANY structural/functional damaging evidence applies.
    """
    # (A) Functional sites
    for col in pm1_rules["functional_damaging"]:
        val = row.get(col, None)
        if isinstance(val, str) and val.lower() == "damaging":
            return True

    # (B) AlloSigma2 damaging
    for col in pm1_rules["allosigma_cols"]:
        val = row.get(col, None)
        if isinstance(val, str) and val.lower() == "damaging":
            return True

    # (C) Local interactions (protein)
    for col in pm1_rules["localint_cols"]:
        val = row.get(col, None)
        if isinstance(val, str) and val.lower() == "destabilizing":
            return True

    # (D) Local interactions (DNA)
    for col in pm1_rules["localint_dna_cols"]:
        val = row.get(col, None)
        if isinstance(val, str) and val.lower() in { 
            "destabilizing", "stabilizing"
         }:
            return True

    # (E) PTM functional effect
    ptm_func = row.get("PTM effect in function [md]", None)
    if isinstance(ptm_func, str) and ptm_func.lower() == "potentially_damaging":
        return True

    # (F) PTM stability effect
    ptm_stab = row.get("PTM effect in stability [md]", None)
    if isinstance(ptm_stab, str) and ptm_stab.lower() in {
        "damaging", "potentially_stabilizing"
    }:
        return True

    # (G) Stability consensus
    stab_cons = row.get("Stability classification, (RaSP, FoldX) [md_25]", None)
    if isinstance(stab_cons, str) and stab_cons.lower() == "destabilizing":
        return True

    return False


###############################################################################
# Main
###############################################################################

def main():
    parser = argparse.ArgumentParser(description="ACMG PP3/BP4 and PM1-like annotation for ClinVar VUS.")
    parser.add_argument(
        "-i", "--input", nargs="+", required=True,
        help="MAVISp CSV file(s) (POLE + POLD1)."
    )
    parser.add_argument(
        "-t", "--thresholds", required=True,
        help="calibrated_thresholds.yaml"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output CSV"
    )
    args = parser.parse_args()

    ###########################################################################
    # 1. Load MAVISp data
    ###########################################################################
    print("Loading MAVISp CSV files...")
    df = pd.concat(
        [pd.read_csv(path, low_memory=False) for path in args.input],
        ignore_index=True
    )
    print(f"Combined variants: {df.shape[0]}")

    ###########################################################################
    # 2. Load thresholds
    ###########################################################################
    thresholds = load_thresholds(args.thresholds)
    print(f"Loaded thresholds for {len(thresholds)} predictors.")

    ###########################################################################
    # 3. Extract VUS
    ###########################################################################
    df["is_VUS"] = df["ClinVar Interpretation"].apply(is_vus)
    vus_df = df[df["is_VUS"]].copy()
    print(f"ClinVar VUS detected: {vus_df.shape[0]}")

    if vus_df.empty:
        print("No VUS found. Exiting.")
        return

    ###########################################################################
    # 4. PP3/BP4 predictors
    ###########################################################################
    pp3bp4_predictors = [
        "REVEL score",
        "AlphaMissense pathogenicity score",
        "EVE score",
        "GEMME Score",
        "DeMaSk delta fitness",
    ]

    # Ensure predictors are numeric
    for pred in pp3bp4_predictors:
        if pred in vus_df.columns:
            vus_df[pred] = pd.to_numeric(vus_df[pred], errors="coerce")

    # Apply PP3/BP4 per predictor
    for pred in pp3bp4_predictors:
        if pred not in vus_df.columns:
            print(f"WARNING: predictor missing: {pred}")
            continue
        thr = thresholds.get(pred, None)
        if thr is None:
            print(f"WARNING: no threshold for {pred}")
            continue

        vus_df[f"{pred}_ACMG"] = vus_df[pred].apply(
            lambda x: pp3_bp4_call(x, pred, thr)
        )

    # Combine overall PP3/BP4
    def combine_pp3bp4(row):
        calls = [row.get(f"{p}_ACMG") for p in pp3bp4_predictors]
        if "PP3" in calls:
            return "PP3"
        if "BP4" in calls:
            return "BP4"
        return None

    vus_df["PP3_BP4"] = vus_df.apply(combine_pp3bp4, axis=1)

    ###########################################################################
    # 5. PM1-like rules
    ###########################################################################

    pm1_rules = {
        "functional_damaging": [
            "Functional sites (active site)",
            "Functional sites (cofactor)",
        ],
        "allosigma_cols": [
            "AlloSigma2-PSN classification - pockets and interfaces [md]",
            "AlloSigma2-PSN classification - active sites [md]",
            "AlloSigma2-PSN classification - cofactor sites [md]",
        ],
        "localint_cols": [
            "Local Int. classification (MCM5_7PLO) [md]",
            "Local Int. classification (PCNA_9B8T_AD) [md]",
            "Local Int. classification (PCNA_9B8T_AB) [md]",
            "Local Int. classification (CDC45_7PLO) [md]",
            "Local Int. classification (POLE2_7PLO) [md]",
            "Local Int. classification (MCM2_7PLO) [md]",
            "Local Int. classification (PCNA_9B8T_AC) [md]",
            "Local Int. classification (LASP1_AFmulti_Huri) [md]",
            "Local Int. classification (PCNA_AFmulti) [md]",
            "Local Int. classification (POLD4_AFmulti) [md]",
            "Local Int. classification (POLD2_6TNY) [md]",
        ],
        "localint_dna_cols": [
            "Local Int. With DNA classification (9F6I_27-1175_retromut_arrest) [md]",
            "Local Int. With DNA classification (9B8T_24-1198_retromut_closed) [md]",
            "Local Int. With DNA classification (9B8S_27-1198_retromut_open) [md]",
            "Local Int. With DNA classification (9F6K_27-1175_frayed) [md]",
            "Local Int. With DNA classification (9F6L_27-1175_mismatch_excision) [md]",
            "Local Int. With DNA classification (6TNY_APT) [md]",
            "Local Int. With DNA classification (6S1M_APT) [md]",
            "Local Int. With DNA classification (6S1O_APT) [md]",
            "Local Int. With DNA classification (6S1N_APT) [md]",
            "Local Int. With DNA classification (6TNZ_APT) [md]",
        ],
    }

    vus_df["PM1_like"] = vus_df.apply(lambda row: pm1_like(row, pm1_rules), axis=1)

    ###########################################################################
    # 6. Output
    ###########################################################################
    vus_df.to_csv(args.output, index=False)
    print(f"Saved ACMG evidence for VUS → {args.output}")


if __name__ == "__main__":
    main()

