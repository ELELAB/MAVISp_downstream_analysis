#!/usr/bin/env python3
import argparse
import csv
import re
import sys
from pathlib import Path
import pandas as pd

AA3_TO_1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
    "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
    "Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Sec":"U","Pyl":"O"
}

MUT_COL_CANDIDATES = [
    "Mutation", "HGVSp", "vep_HGVSp", "AF_mut", "Variant", "variant", "HGVSP", "HGVSp_short"
]

PLDDT_COL_CANDIDATES = [
    # exacts seen across your files
    "AlphaFold2 model pLDDT score",
    "pLDDT", "pLDDT (AF)", "pLDDT_alphafold", "pLDDT alphafold",
    "pLDDT, alphafold", "alphafold_pLDDT", "AF_pLDDT"
]

def guess_column(df: pd.DataFrame, candidates, fuzzy_substr=None):
    if fuzzy_substr is None:
        fuzzy_substr = []
    cols_lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in cols_lower:
            return cols_lower[c.lower()]
    for c in df.columns:
        cl = c.lower()
        if all(s in cl for s in fuzzy_substr):
            return c
    return None

def to_short_mutation(s: str) -> str:
    if not isinstance(s, str):
        return s
    s = s.strip()
    if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]", s):
        return s
    s_clean = re.sub(r"^\w+:\s*p\.", "", s)   # 'NP_...:p.Trp62Cys' -> 'Trp62Cys'
    s_clean = re.sub(r"^p\.", "", s_clean)    # 'p.Trp62Cys' -> 'Trp62Cys'
    s_clean = s_clean.replace("=", "")
    m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", s_clean)
    if m and m.group(1) in AA3_TO_1 and m.group(3) in AA3_TO_1:
        return f"{AA3_TO_1[m.group(1)]}{m.group(2)}{AA3_TO_1[m.group(3)]}"
    m = re.search(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", s_clean)
    if m and m.group(1) in AA3_TO_1 and m.group(3) in AA3_TO_1:
        return f"{AA3_TO_1[m.group(1)]}{m.group(2)}{AA3_TO_1[m.group(3)]}"
    m = re.fullmatch(r"[pP]\.([ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY])", s)
    if m:
        return m.group(1)
    return s

def main():
    parser = argparse.ArgumentParser(
        description="Filter mechanistic-indicator variants by AlphaFold pLDDT from a MAVISp CSV."
    )
    parser.add_argument("-i", "--mavisp", required=True, help="Path to MAVISp CSV (ensemble_mode.csv).")
    parser.add_argument("-m", "--mechanistic", required=True, help="Path to mechanistic indicators CSV.")
    parser.add_argument("-t", "--threshold", required=True, type=float,
                        help="pLDDT threshold: variants with pLDDT < threshold are flagged/removed.")
    parser.add_argument("-o", "--outprefix", default="mech_plddt", help="Output prefix.")
    parser.add_argument("--sep", default=",", help="CSV separator for inputs (default ',').")
    parser.add_argument("--encoding", default="utf-8", help="CSV encoding (default utf-8).")
    parser.add_argument("--plddt-col", default=None,
                        help="Override pLDDT column name (e.g., 'AlphaFold2 model pLDDT score').")
    parser.add_argument("--drop-na", action="store_true",
                        help="Also drop rows with missing pLDDT from the filtered mechanistic file.")
    args = parser.parse_args()

    try:
        df_mav = pd.read_csv(args.mavisp, sep=args.sep, encoding=args.encoding)
    except Exception as e:
        print(f"[ERROR] Could not read MAVISp CSV: {e}", file=sys.stderr); sys.exit(1)
    try:
        df_mech = pd.read_csv(args.mechanistic, sep=args.sep, encoding=args.encoding)
    except Exception as e:
        print(f"[ERROR] Could not read mechanistic CSV: {e}", file=sys.stderr); sys.exit(1)

    mut_col_mav = guess_column(df_mav, MUT_COL_CANDIDATES)
    if mut_col_mav is None:
        print(f"[ERROR] Mutation column not found in MAVISp CSV. Tried: {MUT_COL_CANDIDATES}", file=sys.stderr)
        sys.exit(1)

    if args.plddt_col:
        plddt_col = args.plddt_col
        if plddt_col not in df_mav.columns:
            print(f"[ERROR] Provided --plddt-col '{plddt_col}' not in MAVISp CSV.", file=sys.stderr)
            sys.exit(1)
    else:
        # prefer exact matches incl. your 'AlphaFold2 model pLDDT score'
        plddt_col = guess_column(df_mav, PLDDT_COL_CANDIDATES) \
                    or guess_column(df_mav, [], fuzzy_substr=["plddt"])
        if plddt_col is None:
            print("[ERROR] pLDDT column not found (looked for typical names and '*plddt*').", file=sys.stderr)
            sys.exit(1)

    if "Mutation" not in df_mech.columns:
        mech_mut = guess_column(df_mech, ["Mutation","variant","Variant","HGVSp","vep_HGVSp"])
        if mech_mut is None:
            print("[ERROR] Mechanistic CSV needs a 'Mutation' column (or equivalent).", file=sys.stderr)
            sys.exit(1)
        df_mech = df_mech.rename(columns={mech_mut: "Mutation"})

    # Normalize keys
    df_mech["_mut_key"] = df_mech["Mutation"].map(to_short_mutation)
    df_mav["_mut_key"]   = df_mav[mut_col_mav].map(to_short_mutation)

    # Extract pLDDT and coerce to numeric
    df_mav_sub = df_mav[["_mut_key", plddt_col]].drop_duplicates("_mut_key").rename(columns={plddt_col: "pLDDT"})
    df_mav_sub["pLDDT"] = pd.to_numeric(df_mav_sub["pLDDT"], errors="coerce")

    merged = df_mech.merge(df_mav_sub, on="_mut_key", how="left")

    # classify
    thr = float(args.threshold)
    has_plddt = merged["pLDDT"].notna()
    low_mask = has_plddt & (merged["pLDDT"] < thr)

    low_out = merged[low_mask].copy()
    # filtered: drop low; optionally drop NA as well
    keep_mask = ~low_mask if not args.drop_na else (~low_mask) & merged["pLDDT"].notna()
    filtered_mech = merged[keep_mask].drop(columns=["_mut_key"])

    # Arrange columns for low list
    low_cols = ["Mutation", "pLDDT"] + [c for c in merged.columns if c not in ("Mutation","pLDDT","_mut_key")]
    low_out = low_out[low_cols]

    # Write with quotes preserved
    out_low = f"{args.outprefix}_low_plddt_variants.csv"
    out_mech = f"{args.outprefix}_mechanistic_filtered.csv"
    low_out.to_csv(out_low, index=False, quoting=csv.QUOTE_ALL)
    filtered_mech.to_csv(out_mech, index=False, quoting=csv.QUOTE_ALL)

    # Summary
    total_mech = len(df_mech)
    matched = merged["pLDDT"].notna().sum()
    n_low = low_mask.sum()
    print(f"[INFO] Mechanistic variants total: {total_mech}", file=sys.stderr)
    print(f"[INFO] Matched pLDDT from MAVISp:  {matched}", file=sys.stderr)
    print(f"[INFO] Below threshold (< {thr}):  {n_low}", file=sys.stderr)
    print(f"[OK]  Wrote: {out_low}", file=sys.stderr)
    print(f"[OK]  Wrote: {out_mech}", file=sys.stderr)

if __name__ == "__main__":
    pd.options.mode.copy_on_write = True
    main()

