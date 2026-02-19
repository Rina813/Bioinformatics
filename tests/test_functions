import os
import sys

# Add project root (MyProjectFolder) to path so imports work in pytest
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import pandas as pd
import pytest

from my_module.functions import (
    clean_sample_name,
    add_base_sample_column,
    summarize_duplicates,
    filter_by_target,
)


def test_clean_sample_name_removes_numeric_suffix():
    assert clean_sample_name("KC_sample1_1") == "KC_sample1"
    assert clean_sample_name("KC_sample2_2") == "KC_sample2"
    assert clean_sample_name("Control") == "Control"


def test_add_base_sample_column_creates_expected_column():
    df = pd.DataFrame({
        "Sample": ["KC_sample1_1", "KC_sample1_2", "Control"],
        "Cq": [20.0, 21.0, 18.5],
    })

    df_out = add_base_sample_column(df)

    assert "base_sample" in df_out.columns
    assert df_out.loc[0, "base_sample"] == "KC_sample1"
    assert df_out.loc[1, "base_sample"] == "KC_sample1"
    assert df_out.loc[2, "base_sample"] == "Control"


def test_summarize_duplicates_computes_mean_and_count():
    df = pd.DataFrame({
        "Sample": ["KC_sample1_1", "KC_sample1_2", "Control_1"],
        "Cq": [20.0, 22.0, 18.0],
    })

    df = add_base_sample_column(df)
    summary = summarize_duplicates(df)

    row = summary[summary["base_sample"] == "KC_sample1"].iloc[0]
    assert pytest.approx(row["mean_cq"]) == 21.0
    assert row["n_reps"] == 2

    row_control = summary[summary["base_sample"] == "Control"].iloc[0]
    assert pytest.approx(row_control["mean_cq"]) == 18.0
    assert row_control["n_reps"] == 1


def test_filter_by_target_returns_only_requested_gene():
    df = pd.DataFrame({
        "Target": ["Actb", "Actb", "Stat3"],
        "Sample": ["s1", "s2", "s1"],
        "Cq": [20.0, 21.0, 24.0],
    })

    df_actb = filter_by_target(df, "Actb")
    df_stat3 = filter_by_target(df, "Stat3")

    assert set(df_actb["Target"].unique()) == {"Actb"}
    assert set(df_stat3["Target"].unique()) == {"Stat3"}
