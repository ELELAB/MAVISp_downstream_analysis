import numpy as np
import pandas as pd
import pytest

from dot_plot import dot_plot as dp


def test_convert_to_float_handles_multiple_values():
    assert dp.convert_to_float("0.5, 0.7,0.9") == pytest.approx(0.7)


def test_convert_to_float_passthrough_for_non_string():
    value = np.nan
    assert np.isnan(dp.convert_to_float(value))


def test_map_clinvar_categories_only_updates_clinvar_rows():
    df = pd.DataFrame(
        {
            "ClinVar Interpretation": ["pathogenic", "benign"],
            "Mutation sources": ["ClinVar;cosmic", "literature"],
        },
        index=["mut1", "mut2"],
    )

    clinvar_dict = {"pathogenic": "pathogenic"}
    mapped = dp.map_clinvar_categories(df.copy(), clinvar_dict)

    assert mapped.loc["mut1", "ClinVar Category"] == "pathogenic"
    assert pd.isna(mapped.loc["mut2", "ClinVar Category"])


def _build_summary_and_df():
    summary = pd.DataFrame(
        {
            "Stability classification": [1, 0],
            "Mutation": ["mut1", "mut2"],
        },
        index=["mut1", "mut2"],
    )

    classification = pd.DataFrame(
        {
            "AlphaMissense classification": [1, 0],
            "GEMME predicted consequence": [2, 1],
            "REVEL": [1, 2],
            "EVE classification (25% Uncertain)": [1, 2],
            "DeMaSk predicted consequence": [7, 2],
        },
        index=["mut1", "mut2"],
    )
    return summary, classification


def test_filter_vep_summary_honours_vep_choice():
    summary, classification = _build_summary_and_df()
    filtered = dp.filter_vep_summary(summary, classification, "alphamissense", False)

    assert list(filtered.index) == ["mut1"]
    assert "MAVISp Effects" in filtered.columns


def test_filter_vep_summary_with_demask_filter():
    summary, classification = _build_summary_and_df()
    filtered = dp.filter_vep_summary(summary, classification, "none", True)
    assert list(filtered.index) == ["mut1"]


def test_effect_summary_aggregates_categories():
    df = pd.DataFrame(
        {
            "Stability classification, test": [1, 0],
            "Local Int. classification, test": [0, 1],
            "PTM effect in regulation": [0, 1],
            "AlloSigMA 2 predicted consequence - pockets and interfaces": [1, 0],
            "Functional sites (active site)": [0, 1],
            "Loss of disulfide bridge": [1, 0],
            "Predicted de-novo disulfide bridge": [0, 1],
        },
        index=["mut1", "mut2"],
    )

    enriched = dp.effect_summary(df.copy())
    assert enriched.loc["mut1", "MAVISp Effects"] == "Stability, Long Range, Disulfide bridges"
    assert "Local Int." in enriched.loc["mut2", "MAVISp Effects"]
    assert "Disulfide bridges" in enriched.loc["mut1", "MAVISp Effects"]
