#!/usr/bin/env python3
"""
calibrate_threshold.py

Given:
  - a feature table (from prepare_features_from_yaml.py)
  - MAVISp_label_binary (1 = pathogenic, 0 = benign)

This script:
  - identifies continuous numeric features
  - computes ROC curves and AUC for each feature
  - finds optimal threshold via Youden’s J
  - computes percentile-based thresholds for pathogenic/benign enrichment
  - generates distribution plots for P/LP vs B/LB
  - outputs calibrated thresholds to calibrated_thresholds.yaml

Usage:

  python calibrate_threshold.py \
      -i pole_pold1_features_for_calibration.csv \
      -o calibrated_thresholds.yaml \
      -p calibration_plots

"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from scipy.stats import mannwhitneyu

try:
    import yaml
except ImportError:
    yaml = None


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def save_yaml(data, path):
    if yaml is None:
        raise ImportError(
            "PyYAML is not installed. Install it via:\n\n"
            "   pip install pyyaml\n"
        )
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def compute_youden_threshold(y, scores):
    """
    Computes the threshold that maximizes Youden's J = TPR - FPR.
    Returns:
        threshold_opt, fpr, tpr, auc_value
    """
    fpr, tpr, thresholds = roc_curve(y, scores)
    J = tpr - fpr
    ix = np.argmax(J)
    return thresholds[ix], fpr, tpr, auc(fpr, tpr)


def percentile_thresholds(pos_values, neg_values):
    """
    Returns percentile-based thresholds:
      pathogenic_threshold = 25th percentile of P/LP
      benign_threshold     = 75th percentile of B/LB

    Works only if we have enough numeric values on both sides.
    """
    pos = pd.to_numeric(pos_values, errors="coerce").dropna()
    neg = pd.to_numeric(neg_values, errors="coerce").dropna()

    if len(pos) < 3 or len(neg) < 3:
        return None, None

    pathogenic_thr = float(np.percentile(pos, 25))
    benign_thr = float(np.percentile(neg, 75))
    return pathogenic_thr, benign_thr


def plot_distributions(feature, pos_values, neg_values, outdir):
    """
    Simple histogram-based distribution plot.
    """
    plt.figure(figsize=(6, 4))
    plt.hist(
        pos_values,
        bins=30,
        density=True,
        alpha=0.6,
        label="P/LP",
    )
    plt.hist(
        neg_values,
        bins=30,
        density=True,
        alpha=0.6,
        label="B/LB",
    )
    plt.title(f"Distribution: {feature}")
    plt.xlabel(feature)
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    fname = feature.replace(" ", "_").replace("/", "_")
    plt.savefig(os.path.join(outdir, f"{fname}_dist.png"), dpi=120)
    plt.close()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate feature thresholds using labeled MAVISp variants."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Feature table CSV with MAVISp_label_binary.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output YAML file with calibrated thresholds.",
    )
    parser.add_argument(
        "-p",
        "--plotdir",
        required=True,
        help="Directory to save distribution + ROC plots.",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Load input
    # ------------------------------------------------------------------
    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"❌ Cannot read input {args.input}: {e}", file=sys.stderr)
        sys.exit(1)

    if "MAVISp_label_binary" not in df.columns:
        print("❌ Column MAVISp_label_binary missing in input table.", file=sys.stderr)
        sys.exit(1)

    y = df["MAVISp_label_binary"]
    # ensure numeric labels
    y = pd.to_numeric(y, errors="coerce")

    if not os.path.exists(args.plotdir):
        os.makedirs(args.plotdir)

    thresholds_out = {}
    n_features_used = 0

    # Columns to skip explicitly
    skip_cols = {
        "MAVISp_label",
        "MAVISp_label_binary",
        "Mutation",
        "HGVSp",
        "Gene",
    }

    # ------------------------------------------------------------------
    # Iterate over candidate features
    # ------------------------------------------------------------------
    for feature in df.columns:
        if feature in skip_cols:
            continue

        series = df[feature]

        # Try to coerce to numeric
        Xnum = pd.to_numeric(series, errors="coerce")

        # Drop NaN
        valid_mask = Xnum.notna() & y.notna()
        X = Xnum[valid_mask]
        y_valid = y[valid_mask]

        # Need both classes present
        if len(X) < 6:
            if feature in ("max_local_int_ddg_protein", "max_local_int_ddg_dna"):
               print(f"Skipping {feature}: only {len(X)} non-NaN values")
            continue
        if not ((y_valid == 0).any() and (y_valid == 1).any()):
            continue

        # Skip if binary/boolean/constant (<=2 unique values)
        unique_vals = pd.unique(X)
        if len(unique_vals) <= 2:
            # This is likely binary or constant; not suitable for continuous calibration
            continue

        # Split into P/LP vs B/LB arrays
        pos_values = X[y_valid == 1]
        neg_values = X[y_valid == 0]

        if len(pos_values) < 3 or len(neg_values) < 3:
           if feature in ("max_local_int_ddg_protein", "max_local_int_ddg_dna"):
             print(f"Skipping {feature}: only {len(pos_values)} pos / {len(neg_values)} neg") 
           continue

        # ------------------------------------------------------------------
        # ROC + Youden threshold
        # ------------------------------------------------------------------
        try:
            youden_thr, fpr, tpr, auc_value = compute_youden_threshold(y_valid, X)
        except Exception as e:
            print(f"⚠️ Skipping {feature}: cannot compute ROC ({e})", file=sys.stderr)
            continue

        # ------------------------------------------------------------------
        # Percentile thresholds
        # ------------------------------------------------------------------
        p_thr, b_thr = percentile_thresholds(pos_values, neg_values)
        if p_thr is None:
            # not enough data coverage
            continue

        # ------------------------------------------------------------------
        # Distribution plot
        # ------------------------------------------------------------------
        try:
            plot_distributions(feature, pos_values, neg_values, args.plotdir)
        except Exception as e:
            print(f"⚠️ Could not plot distributions for {feature}: {e}", file=sys.stderr)

        # ------------------------------------------------------------------
        # ROC curve plot
        # ------------------------------------------------------------------
        try:
            plt.figure(figsize=(6, 4))
            plt.plot(fpr, tpr, label=f"AUC={auc_value:.3f}")
            plt.plot([0, 1], [0, 1], "k--")
            plt.title(f"ROC: {feature}")
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.legend()
            plt.tight_layout()
            fname = feature.replace(" ", "_").replace("/", "_")
            plt.savefig(os.path.join(args.plotdir, f"{fname}_ROC.png"), dpi=120)
            plt.close()
        except Exception as e:
            print(f"⚠️ Could not plot ROC for {feature}: {e}", file=sys.stderr)

        # ------------------------------------------------------------------
        # Mann–Whitney U test
        # ------------------------------------------------------------------
        try:
            mw_p = float(
                mannwhitneyu(pos_values, neg_values, alternative="two-sided").pvalue
            )
        except Exception:
            mw_p = float("nan")

        # ------------------------------------------------------------------
        # Save thresholds + metrics for this feature
        # ------------------------------------------------------------------
        thresholds_out[feature] = {
            "auc": float(auc_value),
            "youden_threshold": float(youden_thr),
            "percentile_pathogenic_25th": float(p_thr),
            "percentile_benign_75th": float(b_thr),
            "mannwhitney_u_pvalue": mw_p,
            "n_positive": int(len(pos_values)),
            "n_negative": int(len(neg_values)),
        }

        n_features_used += 1

    # ----------------------------------------------------------------------
    # Write YAML output
    # ----------------------------------------------------------------------
    print(f"Calibration complete. Features processed: {n_features_used}")
    print(f"Writing thresholds to {args.output}")

    save_yaml({"thresholds": thresholds_out}, args.output)

    print("Done.")


if __name__ == "__main__":
    main()

