import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


DOT_PLOT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = DOT_PLOT_DIR / "tests" / "data"

DATASETS = {
    "simple_mode": {
        "input": DATA_DIR / "simple_mode" / "SMPD1-simple_mode.csv",
        "dictionary": DATA_DIR / "simple_mode" / "dictionary.csv",
        "expected": DATA_DIR / "simple_mode" / "expected" / "mechanistic_indicators_out.csv",
    },
    "ensemble_mode": {
        "input": DATA_DIR / "ensemble_mode" / "TP53-ensemble_mode.csv",
        "dictionary": DATA_DIR / "ensemble_mode" / "dictionary.csv",
        "expected": DATA_DIR / "ensemble_mode" / "expected" / "mechanistic_indicators_out.csv",
    },
}


def run_cli(tmp_path, dataset_key, extra_args=None):
    dataset = DATASETS[dataset_key]
    output_prefix = tmp_path / f"{dataset_key}_dot_plot"
    cmd = [
        sys.executable,
        str(DOT_PLOT_DIR / "dot_plot.py"),
        "-i",
        str(dataset["input"]),
        "-v",
        str(dataset["dictionary"]),
        "-o",
        str(output_prefix),
    ]
    if extra_args:
        cmd.extend(extra_args)

    subprocess.run(cmd, cwd=tmp_path, check=True)
    return {
        "pdf": output_prefix.with_suffix(".pdf"),
        "log": tmp_path / "log.txt",
        "mechanistic": tmp_path / "mechanistic_indicators_out.csv",
    }


def assert_csvs_equal(actual, expected):
    actual_df = pd.read_csv(actual)
    expected_df = pd.read_csv(expected)
    pd.testing.assert_frame_equal(
        actual_df.sort_index(axis=1), expected_df.sort_index(axis=1), check_like=True
    )


@pytest.mark.integration
def test_simple_mode_cli(tmp_path):
    paths = run_cli(tmp_path, "simple_mode")

    assert paths["pdf"].exists() and paths["pdf"].stat().st_size > 0
    assert paths["log"].exists()
    assert paths["mechanistic"].exists()

    assert_csvs_equal(paths["mechanistic"], DATASETS["simple_mode"]["expected"])

    log_content = paths["log"].read_text()
    assert "AlphaMissense" in log_content


@pytest.mark.integration
def test_ensemble_mode_cli(tmp_path):
    paths = run_cli(tmp_path, "ensemble_mode")

    assert paths["pdf"].exists()
    assert paths["mechanistic"].exists()
    assert_csvs_equal(paths["mechanistic"], DATASETS["ensemble_mode"]["expected"])


@pytest.mark.integration
def test_vep_filter_eve_produces_empty_result(tmp_path):
    paths = run_cli(tmp_path, "simple_mode", extra_args=["-vep", "eve"])
    actual_df = pd.read_csv(paths["mechanistic"])
    assert actual_df.empty
