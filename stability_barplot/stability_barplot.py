#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Stability plotting script with grouped methods and per-chunk plots.
Generates CSVs with only the columns used for each plot and saves figures
with the legend outside the plotting area and x-labels at 45°.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import re
import matplotlib.cm as cm


# ============================================================
# REGEX TO PARSE STABILITY COLUMNS
# ============================================================

ddg_ensemble_regex = re.compile(r"^Stability \((.+?), kcal/mol\) \[(.+?)\]$")
std_ensemble_regex = re.compile(r"^Stability \((.+?), kcal/mol[, ]*st\.?\s*dev\.?\) \[(.+?)\]$")
classification_regex = re.compile(r"^Stability classification, *\((?P<methods>[^)]+)\) *\[(?P<tag>[^\]]+)\]$")

# === simple_mode regex ===
ddg_simple_regex = re.compile(r"^Stability\s*\((?P<method>[^,]+),\s*(?P<ensemble>[^,]+),\s*kcal/mol\)$")
std_simple_regex = re.compile(r"^Stability\s*\((?P<method>[^,]+),\s*(?P<ensemble>[^,]+),\s*kcal/mol,\s*st\. dev\.\)$")
classification_simple_regex = re.compile(r"^Stability classification,\s*(?P<ensemble>[^,]+),\s*\((?P<methods>[^)]+)\)$")



def extract_columns(df,fmt):
    parsed = {}
    classification = {}

    for col in df.columns:
        clean = col.strip()

        # ===== ENSEMBLE MODE =====
        if fmt == "ensemble":
            m_ddg = ddg_ensemble_regex.match(clean)
            if m_ddg:
                method = m_ddg.group(1).strip()
                ensemble = m_ddg.group(2).strip()
                parsed.setdefault(method, {}).setdefault(ensemble, {})['ddg'] = col
                continue

            m_std = std_ensemble_regex.match(clean)
            if m_std:
                method = m_std.group(1).strip()
                ensemble = m_std.group(2).strip()
                parsed.setdefault(method, {}).setdefault(ensemble, {})['std'] = col
                continue

            m_class = classification_regex.match(clean)
            if m_class:
                ensemble = m_class.group("tag").strip()
                classification.setdefault(ensemble, []).append(col)
                continue

        # ===== SIMPLE MODE =====
        elif fmt == "simple":
            m_ddg = ddg_simple_regex.match(clean)
            if m_ddg:
                method = m_ddg.group("method").strip()
                ensemble = m_ddg.group("ensemble").strip()
                parsed.setdefault(method, {}).setdefault(ensemble, {})['ddg'] = col
                continue

            m_std = std_simple_regex.match(clean)
            if m_std:
                method = m_std.group("method").strip()
                ensemble = m_std.group("ensemble").strip()
                parsed.setdefault(method, {}).setdefault(ensemble, {})['std'] = col
                continue

            m_class = classification_simple_regex.match(clean)
            if m_class:
                ensemble = m_class.group("ensemble").strip()
                classification.setdefault(ensemble, []).append(col)
                continue

    return parsed, classification


def detect_format(df):
    cols = df.columns
    for col in cols:
        if ddg_ensemble_regex.match(col) or std_ensemble_regex.match(col) or classification_regex.match(col):
            return "ensemble"
        if ddg_simple_regex.match(col) or std_simple_regex.match(col) or classification_simple_regex.match(col):
            return "simple"
    return "unknown"


def extract_method_name(col_name):
    # ensemble_mode
    m_ddg = ddg_ensemble_regex.match(col_name)
    if m_ddg:
        return m_ddg.group(1).strip()

    m_std = std_ensemble_regex.match(col_name)
    if m_std:
        return m_std.group(1).strip()

    m_class = classification_regex.match(col_name)
    if m_class:
        return [x.strip() for x in m_class.group("methods").split(",")]

    # simple_mode
    m_ddg_s = ddg_simple_regex.match(col_name)
    if m_ddg_s:
        return m_ddg_s.group("method").strip()
    
    m_std_s = std_simple_regex.match(col_name)
    if m_std_s:
        return m_std_s.group(1).strip()

    m_class_s = classification_simple_regex.match(col_name)
    if m_class_s:
        return [x.strip() for x in m_class_s.group("methods").split(",")]

    return col_name


def dynamic_errorbar_style(n_bars, base_capsize=5, base_linewidth=2, base_capthick=2):
    """
    Reduces the size/line width of the error bars based on the total number of bars.
    """
    scale = max(0.3, min(1.0, 20 / n_bars))  
    # 20 = ideal number of bars before starting reducing

    return {
        "capsize": base_capsize * scale,
        "error_kw": {
            "elinewidth": base_linewidth * scale,
            "capthick": base_capthick * scale,
            "ecolor": "black"
        }
    }

def filter_destabilizing(df, ensemble, classification_cols, target_methods):
    """
    Keeps only the rows where at least one classification column for the ensemble 
    contains 'Destabilizing' for the target methods.
    """
    if ensemble not in classification_cols:
        return df

    cols = classification_cols[ensemble]
    mask = pd.Series(False, index=df.index)

    for col in cols:
        for method in target_methods:
            # check if the column is in method
            if method in col:
                # check for destabilizing 
                mask |= df[col].fillna("").str.contains("Destabilizing")

    return df[mask]


# ============================================================
# PLOTTING FUNCTIONS
# ============================================================
def get_palette_color(cmap_name, index):
    cmap = cm.get_cmap(cmap_name)
    return cmap(index)
    

def plot_chunk(chunk, method, ensemble, ddg_col, std_col, args, idx):
    """Plot a single method + ensemble per chunk, using palette-based colors only."""
    
    chunk = chunk.copy()
    ddg = pd.to_numeric(chunk[ddg_col], errors="coerce")
    std = pd.to_numeric(chunk[std_col], errors="coerce") if std_col else None
    labels = chunk["Mutation"]

    # =============================
    #  Auto-coloring from palette
    # =============================
    method_lower = method.lower()

    if "foldx" in method_lower:
        bar_color = get_palette_color(args.palette, 0)
    elif "rosetta" in method_lower:
        bar_color = get_palette_color(args.palette, 1)
    elif "rasp" in method_lower:
        bar_color = get_palette_color(args.palette, 2)
    else:
        bar_color = get_palette_color(args.palette, 3)

    # =============================
    #  AF label colors (unchanged)
    # =============================
    if args.color_AF and "AlphaFold2 model pLDDT score" in chunk.columns:
        colors = chunk["AlphaFold2 model pLDDT score"].apply(
            lambda x: 'orange' if x < 50 else
                      ('yellow' if 50 <= x < 70 else
                       ('cyan' if 70 <= x < 90 else 'blue'))
        )
    else:
        colors = ['black'] * len(chunk)

    fig, ax = plt.subplots(figsize=(args.figure_width/2.54, args.figure_height/2.54))
    n_bars = len(labels)
    style = dynamic_errorbar_style(n_bars)

    # =============================
    #  VERTICAL MODE
    # =============================
    if args.barh:
        ax.barh(labels, ddg, xerr=std, color=bar_color, **style)
        ax.invert_yaxis()
        ax.set_xlabel("ΔΔG (kcal/mol)", fontsize=args.y_title)
        ax.set_ylabel("Mutation", fontsize=args.x_title)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(
            labels,
            fontsize=args.x_labels,
            fontweight=args.x_labels_weight
        )

        for tick, c in zip(ax.get_yticklabels(), colors):
            tick.set_color(c)

        ax.set_title(f"{method} — {ensemble} — chunk {idx}")
        return fig

    # =============================
    #  HORIZONTAL MODE
    # =============================
    ax.bar(labels, ddg, yerr=std, color=bar_color, **style)
    ax.set_ylabel("ΔΔG (kcal/mol)", fontsize=args.y_title)
    ax.set_xlabel("Mutation", fontsize=args.x_title)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(
        labels,
        rotation=45,
        ha='right',
        fontsize=args.x_labels,
        fontweight=args.x_labels_weight
    )

    for tick, c in zip(ax.get_xticklabels(), colors):
        tick.set_color(c)

    ax.set_title(f"{method} — {ensemble} — chunk {idx}")
    ax.legend([method], bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.subplots_adjust(right=0.75)

    return fig

def plot_grouped_methods(chunk, parsed, ensemble, args, idx, specific_cols=None):
    """
    Grouped-methods plot with palette-based auto-coloring.
    """
    mutations = chunk["Mutation"].tolist()
    data, errors, legend_labels, colors_bars = [], [], [], []

    # -------------------------------------------
    # Extract methods and assign sequential colors
    # -------------------------------------------
    method_list = []

    if specific_cols:
        for col in specific_cols:
            method_list.append(extract_method_name(col))
    else:
        for method, ensembles in parsed.items():
            if ensemble in ensembles:
                method_list.append(method)

    # Assign dynamically colors
    color_map = {
        method_list[i]: get_palette_color(args.palette, i)
        for i in range(len(method_list))
    }

    # -------------------------------------------
    # Collect values
    # -------------------------------------------
    if specific_cols:
        for col in specific_cols:
            ddg_vals = pd.to_numeric(chunk[col], errors="coerce").tolist()
            method_name = extract_method_name(col)

            # find std col
            std_col = None
            if method_name in parsed:
                for ens, v in parsed[method_name].items():
                    if v.get("ddg") == col:
                        std_col = v.get("std")
                        break

            std_vals = pd.to_numeric(chunk[std_col], errors="coerce").tolist() if std_col else [0]*len(ddg_vals)

            data.append(ddg_vals)
            errors.append(std_vals)
            legend_labels.append(method_name)
            colors_bars.append(color_map[method_name])

    else:
        for method, ensembles in parsed.items():
            if ensemble not in ensembles:
                continue

            ddg_col = ensembles[ensemble].get("ddg")
            std_col = ensembles[ensemble].get("std", None)

            ddg_vals = pd.to_numeric(chunk[ddg_col], errors="coerce").tolist()
            std_vals = pd.to_numeric(chunk[std_col], errors="coerce").tolist() if std_col else [0]*len(ddg_vals)

            data.append(ddg_vals)
            errors.append(std_vals)
            legend_labels.append(method)
            colors_bars.append(color_map[method])

    # -------------------------------------------
    # Plotting
    # -------------------------------------------
    n_methods = len(legend_labels)
    x = np.arange(len(mutations))
    bar_width = min(0.8 / n_methods, args.bars_width)
    offsets = (np.arange(n_methods) - (n_methods - 1) / 2) * bar_width

    fig, ax = plt.subplots(figsize=(args.figure_width/2.54, args.figure_height/2.54))
    n_bars = len(mutations) * n_methods
    style = dynamic_errorbar_style(n_bars)
    # =============================
    #  AF label colors (unchanged)
    # =============================
    if args.color_AF and "AlphaFold2 model pLDDT score" in chunk.columns:
        colors = chunk["AlphaFold2 model pLDDT score"].apply(
            lambda x: 'orange' if x < 50 else
                      ('yellow' if 50 <= x < 70 else
                       ('cyan' if 70 <= x < 90 else 'blue'))
        )
    else:
        colors = ['black'] * len(chunk)

    for i in range(n_methods):
        if args.barh:
            ax.barh(x + offsets[i], data[i], xerr=errors[i],
                    color=colors_bars[i], height=bar_width,
                    label=legend_labels[i], **style)
        else:
            ax.bar(x + offsets[i], data[i], yerr=errors[i],
                   color=colors_bars[i], width=bar_width,
                   label=legend_labels[i], **style)

    # Labels
    if args.barh:
        ax.set_yticks(x)
        ax.set_yticklabels(
            mutations,
            fontsize=args.x_labels,
            fontweight=args.x_labels_weight
        )

        ax.set_xlabel("ΔΔG (kcal/mol)", fontsize=args.y_title)
        ax.set_ylabel("Mutation", fontsize=args.x_title)
        ax.invert_yaxis()
        for tick, c in zip(ax.get_yticklabels(), colors):
            tick.set_color(c)
    else:
        ax.set_xticks(x)
        ax.set_xticklabels(
            mutations,
            rotation=45,
            ha='right',
            fontsize=args.x_labels,
            fontweight=args.x_labels_weight
        )

        ax.set_ylabel("ΔΔG (kcal/mol)", fontsize=args.y_title)
        ax.set_xlabel("Mutation", fontsize=args.x_title)
        for tick, c in zip(ax.get_xticklabels(), colors):
            tick.set_color(c)

    ax.set_title(f"Grouped methods — ensemble {ensemble} — chunk {idx}")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.subplots_adjust(right=0.75)

    return fig, legend_labels, data, errors

# ============================================================
# FILE LOADING
# ============================================================

def load_and_combine(input_folder, grouped_columns=None):
    filename_regex = re.compile(r"^(?P<protein>[^-]+)-(?P<mode>simple|ensemble)_.*\.(csv|tsv)$")

    dfs = {"simple": [], "ensemble": []}
    filenames = {"simple": [], "ensemble": []}

    for fname in os.listdir(input_folder):
        if not (fname.endswith(".csv") or fname.endswith(".tsv")):
            continue

        match = filename_regex.match(fname)
        if not match:
            raise ValueError(
                f"Invalid file name format: '{fname}'. Expected '{{PROTEIN}}-{{simple|ensemble}}_*.csv'"
            )

        protein = match.group("protein")
        mode = match.group("mode")

        fpath = os.path.join(input_folder, fname)
        print(f"Loading {fname}")

        df = pd.read_csv(fpath) if fname.endswith(".csv") else pd.read_csv(fpath, sep="\t")
        df["protein"] = protein

        dfs[mode].append(df)
        filenames[mode].append(fname)

    combined = {}

    for mode, list_of_dfs in dfs.items():
        if not list_of_dfs:
            continue

        if grouped_columns:
            valid_dfs = []
            excluded = []

            for df, fname in zip(list_of_dfs, filenames[mode]):
                if all(col in df.columns for col in grouped_columns):
                    valid_dfs.append(df)
                else:
                    excluded.append(fname)

            if len(valid_dfs) == 0:
                raise ValueError(
                    f"ERROR: None of the '{mode}' files contain all columns required by --grouped_columns: "
                    f"{grouped_columns}"
                )

            if excluded:
                print("\n WARNING: The following files were excluded because "
                      "they do NOT contain all grouped_columns:")
                for f in excluded:
                    print(f"  - {f}")
                print("Proceeding with only the valid files.\n")

            combined[mode] = pd.concat(valid_dfs, ignore_index=True)

        else:
            combined[mode] = pd.concat(list_of_dfs, ignore_index=True)

    if not combined:
        raise ValueError("No valid {PROTEIN}-{simple|ensemble}_*.csv files found.")

    return combined



# ============================================================
# FILESYSTEM UTILS
# ============================================================
def ensure(path):
    if not os.path.exists(path):
        os.makedirs(path)


# ============================================================
# MAIN SCRIPT
# ============================================================
def main():
    parser = argparse.ArgumentParser()

    # Input / output
    parser.add_argument("-i", "--input_folder",
                        dest="input_folder",
                        required=True,
                        help="Path to the folder containing the input CSV files.")

    parser.add_argument("-o", "--output_folder",
                        dest = "output_folder",
                        required = False,
                        default = "output",
                        help = "Path to the folder where output figures and CSVs will be saved.")
    
    parser.add_argument("-m","--mutation_list",
                        dest="mutation_list",
                        required=False,
                        default=None,
                        help="Optional CSV file with columns 'protein,Mutation' used to filter the dataset before plotting.")

    # Chunking & plotting
    parser.add_argument("-ch", "--chunk_size",
                        dest="chunk_size",
                        type=int,
                        default=50,
                        required=False,
                        help="Number of mutations per chunk when plotting barplots. Default: 10.")

    parser.add_argument("-g", "--grouped_columns",
                        dest="grouped_columns",
                        nargs="+",
                        default=None,
                        required=False,
                        help="Specific columns to include in grouped methods plots. Default: None (all methods).")

    parser.add_argument("-w", "--figure_width",
                        dest="figure_width",
                        type=float,
                        default=18,
                        required=False,
                        help="Figure width in cm. Default: 18.")

    parser.add_argument("-he", "--figure_height",
                        dest="figure_height",
                        type=float,
                        default=21.5,
                        required=False,
                        help="Figure height in cm. Default: 21.5.")

    parser.add_argument("-bw", "--bars_width",
                        dest="bars_width",
                        type=float,
                        default=0.45,
                        required=False,
                        help="Maximum width of bars in the barplot. Default: 0.45.")
    
    parser.add_argument("-xf","--x_font_labels",
                       dest="x_labels",
                       type=int,
                       default=7,
                       required=False,
                       help="Font size of x-axis labels (mutations). Default: 7.")

    parser.add_argument("-xw","--weight_of_xlabels",
                        dest="x_labels_weight",
                        type=str,
                        default="regular",
                        required=False,
                        help="Font weight of x-axis labels. Options: "
                             "ultralight, light, normal, regular, book, medium, "
                             "semibold, demibold, bold, heavy, extra bold, black. "
                             "Default: regular.")

    parser.add_argument("-xt","--x_title_font",
                        dest="x_title",
                        type=int,
                        default=15,
                        required=False,
                        help="Font size of x-axis title (Mutation). Default: 15.")

    parser.add_argument("-yt","--y_title_font",
                        dest="y_title",
                        type=int,
                        default=15,
                        required=False,
                        help="Font size of y-axis title (ΔΔG). Default: 15.")


    # Orientation
    parser.add_argument("-v", "--vertical_plot",
                        dest="barh",
                        default=False,
                        required=False,
                        action="store_true",
                        help="Plot bars vertically (horizontal plot is default).")

    # Destabilizing mutations flags
    parser.add_argument("-dr", "--only_rosetta_destabilizing",
                        dest="dest_rosetta",
                        default=False,
                        required=False,
                        action="store_true",
                        help="Plot only the destabilizing mutations according to FoldX and Rosetta.")

    parser.add_argument("-dp", "--only_rasp_destabilizing",
                        dest="dest_rasp",
                        default=False,
                        required=False,
                        action="store_true",
                        help="Plot only the destabilizing mutations according to FoldX and RaSP.")

    # Color options
    parser.add_argument("-c", "--color_accordingly_AF",
                        dest="color_AF",
                        default=False,
                        required=False,
                        action="store_true",
                        help="Color x-labels according to the AlphaFold2 pLDDT score. Default: False.")

    parser.add_argument(
                        "-p", "--palette",
                        dest="palette",
                        type=str,
                        default="tab10",
                        required=False,
                        help=(
                            "Color palette to automatically assign colors to methods (FoldX, Rosetta, RaSP). "
                            "The selected palette will generate N distinct colors based on the number of methods.\n\n"
                            "Available palettes (matplotlib):\n"
                            "  - tab10\n"
                            "  - tab20\n"
                            "  - Set1, Set2, Set3\n"
                            "  - Pastel1, Pastel2\n"
                            "  - Dark2\n"
                            "  - Accent\n"
                            "  - Paired\n"
                            "  - viridis, plasma, inferno, magma, cividis\n"
                            "  - coolwarm, bwr, seismic\n"
                            "  - turbo\n\n"
                            "Default: tab10"
                        )
                    )

    args = parser.parse_args()
    combined = load_and_combine(args.input_folder,args.grouped_columns)
    # -----------------------------------------
    # -----------------------------------------
    # OPTIONAL MUTATION+PROTEIN FILTERING
    # -----------------------------------------
    for mode,df in combined.items():
        output_folder = os.path.join(args.output_folder,mode)
        if args.mutation_list:
            if not os.path.isfile(args.mutation_list):
                raise FileNotFoundError(f"Mutation list file not found: {args.mutation_list}")

            print(f"Applying mutation+protein filter from: {args.mutation_list}")
            mdf = pd.read_csv(args.mutation_list)

            # Validate columns
            if "Mutation" not in mdf.columns or "protein" not in mdf.columns:
                raise ValueError("The mutation list CSV must contain 'protein' and 'Mutation' columns.")

            # Build set of allowed pairs
            allowed_pairs = set(
                tuple(x) for x in mdf[["protein", "Mutation"]].astype(str).to_numpy()
            )

            print(f"Allowed (protein, Mutation) pairs: {len(allowed_pairs)}")

            # Ensure df has needed columns
            if "protein" not in df.columns:
                raise ValueError("Input dataframe does not contain a 'protein' column required for filtering.")
            if "Mutation" not in df.columns:
                raise ValueError("Input dataframe does not contain a 'Mutation' column required for filtering.")

            # Apply filtering
            df = df[
                df.apply(
                    lambda row: (str(row["protein"]), str(row["Mutation"])) in allowed_pairs,
                    axis=1
                )
            ]

            print(f"Filtered dataframe now has {len(df)} rows")

            if len(df) == 0:
                raise ValueError("After filtering, no rows remain. Check your mutation list.")

        parsed,classification = extract_columns(df,mode)
        
        ensure(output_folder)

        for ensemble in classification:
            if args.dest_rosetta:
                df = filter_destabilizing(df, ensemble, classification, ["Rosetta", "FoldX"])
            if args.dest_rasp:
                df = filter_destabilizing(df, ensemble, classification, ["RaSP", "FoldX"])

        N_chunks = int(np.ceil(len(df) / args.chunk_size))
        for idx in range(N_chunks):
            chunk = df.iloc[idx*args.chunk_size:(idx+1)*args.chunk_size]

            for method, ensembles in parsed.items():
                for ensemble in ensembles:
                    ddg_col = ensembles[ensemble].get("ddg")
                    if not ddg_col:
                        continue
                    std_col = ensembles[ensemble].get("std", None)
                    # Standard plot
                    out_dir = os.path.join(output_folder, ensemble, method)
                    ensure(out_dir)
                    fig = plot_chunk(chunk, method, ensemble, ddg_col, std_col, args, idx)
                    outfile = os.path.join(out_dir, f"{method}_{ensemble}_chunk_{idx}.png")
                    fig.savefig(outfile, dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    print("Saved:", outfile)

            # Grouped methods plot
            for ensemble in set([e for m in parsed for e in parsed[m]]):
                grouped_dir = os.path.join(output_folder, ensemble, "grouped_methods")
                ensure(grouped_dir)
                fig2, legend_labels, data_vals, errors_vals = plot_grouped_methods(
                    chunk, parsed, ensemble, args, idx, specific_cols=args.grouped_columns
                )
                outfile2 = os.path.join(grouped_dir, f"grouped_{ensemble}_chunk_{idx}.png")
                fig2.savefig(outfile2, dpi=300, bbox_inches='tight')
                plt.close(fig2)
                print("Saved grouped:", outfile2)

                # Save CSV
                csv_cols = ["Mutation"]
                if args.grouped_columns:
                    for col in args.grouped_columns:
                        csv_cols.append(col)
                        method_name = extract_method_name(col)
                        # robust std detection
                        std_col = None
                        if method_name in parsed:
                            for ens, v in parsed[method_name].items():
                                if "std" in v and v.get("ddg") == col:
                                    std_col = v["std"]
                                    break
                        if std_col:
                            csv_cols.append(std_col)
                else:
                    for method_name, ensembles_dict in parsed.items():
                        if ensemble in ensembles_dict:
                            csv_cols.append(ensembles_dict[ensemble]["ddg"])
                            if "std" in ensembles_dict[ensemble]:
                                csv_cols.append(ensembles_dict[ensemble]["std"])
                csv_file = os.path.join(grouped_dir, f"grouped_{ensemble}_chunk_{idx}.csv")
                chunk[csv_cols].to_csv(csv_file, index=False)


if __name__ == "__main__":
    main()