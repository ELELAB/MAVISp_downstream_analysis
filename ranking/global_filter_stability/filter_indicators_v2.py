#!/usr/bin/env python3
import pandas as pd
import argparse
import csv
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Filter MAVISp CSV by mutations from mechanistic indicators file, optionally dropping selected columns."
    )
    parser.add_argument(
        "-m", "--mech", required=True,
        help="Mechanistic indicators CSV (with 'Mutation' column)"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Full MAVISp CSV file"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output filtered MAVISp CSV"
    )
    parser.add_argument(
        "-l", "--outlist", required=True,
        help="Output mutation list file (TSV format)"
    )
    parser.add_argument(
        "-d", "--drop-cols", nargs="+", default=[],
        help=(
            "Columns to remove from the final output CSV. "
            "Pass each column name as its own quoted argument. "
            "Example: -d \"Col A\" \"Col B with, comma\""
        )
    )

    args = parser.parse_args()

    # use drop-cols as provided (no splitting on commas!)
    drop_cols = [c.strip() for c in args.drop_cols]

    # 1. Read mechanistic indicators file
    mech_df = pd.read_csv(args.mech, quotechar='"')
    if "Mutation" not in mech_df.columns:
        sys.exit("❌ Error: Column 'Mutation' not found in mechanistic indicators file.")

    mut_list = (
        mech_df["Mutation"]
        .dropna()
        .astype(str)
        .str.strip()
        .unique()
        .tolist()
    )

    # 2. Read full MAVISp CSV
    full_df = pd.read_csv(args.input, quotechar='"')
    if "Mutation" not in full_df.columns:
        sys.exit("❌ Error: Column 'Mutation' not found in MAVISp CSV file.")

    # 3. Filter rows by mutation list
    filtered_df = full_df[full_df["Mutation"].isin(mut_list)].copy()

    # 4. Drop user-specified columns, if they exist
    if drop_cols:
        existing_cols = [c for c in drop_cols if c in filtered_df.columns]
        missing_cols = [c for c in drop_cols if c not in filtered_df.columns]

        if missing_cols:
            print("⚠️  Warning: these columns were not found and were skipped:")
            for c in missing_cols:
                print(f"   - {c}")

        if existing_cols:
            filtered_df.drop(columns=existing_cols, inplace=True)
            print(f"🗑️  Dropped {len(existing_cols)} column(s): {', '.join(existing_cols)}")

    # 5. Save filtered MAVISp CSV (quoted everywhere)
    filtered_df.to_csv(args.output, index=False, quoting=csv.QUOTE_ALL)

    # 6. Write mutation list + found/not found info
    found = set(filtered_df["Mutation"].unique())
    with open(args.outlist, "w") as fh:
        fh.write("Mutation\tfound_in_MAVISp\n")
        for mut in mut_list:
            fh.write(f"{mut}\t{'yes' if mut in found else 'no'}\n")

    # 7. Summary to terminal
    not_found = [m for m in mut_list if m not in found]
    print("\n📊 SUMMARY")
    print("────────────────────────────────────────")
    print(f"Total mutations in mechanistic file: {len(mut_list)}")
    print(f"Matched in MAVISp file: {len(found)}")
    print(f"Output written to: {args.output}")
    print(f"Mutation list saved as: {args.outlist}")
    if not_found:
        print("\n⚠️  Mutations not found in MAVISp file:")
        for m in not_found:
            print(f"   - {m}")
    else:
        print("\n✅ All mutations were found!")

if __name__ == "__main__":
    main()

