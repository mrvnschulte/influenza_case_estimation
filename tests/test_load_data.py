import numpy as np
import pandas as pd
import pytest

from influenza_case_estimation.load_data import (
    load_strain_fractions,
    load_symptomatic_illnesses,
    load_test_positivity,
    load_weighted_ili_rate,
)


def write_cdc_csv(tmp_path, filename, data):
    """Write a CDC-CSV with the expected first row by the loader."""
    path = tmp_path / filename
    pd.DataFrame(data).to_csv(path, index=False)

    content = path.read_text()
    path.write_text("CDC data export\n" + content)

    return path


class TestLoadWeightedIliRate:
    def test_loads_national_records_as_proportions(self, tmp_path):
        path = write_cdc_csv(
            tmp_path,
            "ili.csv",
            {
                "REGION TYPE": ["National", "HHS Region"],
                "YEAR": [2021, 2021],
                "WEEK": [2, 1],
                "% WEIGHTED ILI": [1.5, 9.0],
            },
        )

        result = load_weighted_ili_rate(path)

        assert result.columns.tolist() == [
            "year",
            "week",
            "weighted_ili_rate",
        ]
        assert result[["year", "week"]].values.tolist() == [[2021, 2]]
        assert result.loc[0, "weighted_ili_rate"] == pytest.approx(0.015)

    def test_sorts_by_year_and_week(self, tmp_path):
        path = write_cdc_csv(
            tmp_path,
            "ili.csv",
            {
                "REGION TYPE": ["National", "National", "National"],
                "YEAR": [2021, 2020, 2021],
                "WEEK": [2, 52, 1],
                "% WEIGHTED ILI": [2.0, 3.0, 1.0],
            },
        )

        result = load_weighted_ili_rate(path)

        assert result[["year", "week"]].values.tolist() == [
            [2020, 52],
            [2021, 1],
            [2021, 2],
        ]

    def test_raises_when_no_national_records_exist(self, tmp_path):
        path = write_cdc_csv(
            tmp_path,
            "ili.csv",
            {
                "REGION TYPE": ["HHS Region"],
                "YEAR": [2021],
                "WEEK": [1],
                "% WEIGHTED ILI": [1.5],
            },
        )

        with pytest.raises(ValueError, match="No national records"):
            load_weighted_ili_rate(path)


class TestLoadTestPositivity:
    def test_loads_national_records_as_proportions(self, tmp_path):
        path = write_cdc_csv(
            tmp_path,
            "positivity.csv",
            {
                "REGION TYPE": ["National", "HHS Region"],
                "YEAR": [2021, 2021],
                "WEEK": [2, 1],
                "PERCENT POSITIVE": [18.0, 50.0],
            },
        )

        result = load_test_positivity(path)

        assert result.columns.tolist() == [
            "year",
            "week",
            "test_positivity",
        ]
        assert result.loc[0, "test_positivity"] == pytest.approx(0.18)

    @pytest.mark.parametrize("percent_positive", [-1.0, 101.0])
    def test_raises_for_positivity_outside_valid_range(
        self, tmp_path, percent_positive
    ):
        path = write_cdc_csv(
            tmp_path,
            "positivity.csv",
            {
                "REGION TYPE": ["National"],
                "YEAR": [2021],
                "WEEK": [1],
                "PERCENT POSITIVE": [percent_positive],
            },
        )

        with pytest.raises(ValueError, match="between 0 and 1"):
            load_test_positivity(path)

    def test_raises_when_no_national_records_exist(self, tmp_path):
        path = write_cdc_csv(
            tmp_path,
            "positivity.csv",
            {
                "REGION TYPE": ["HHS Region"],
                "YEAR": [2021],
                "WEEK": [1],
                "PERCENT POSITIVE": [10.0],
            },
        )

        with pytest.raises(ValueError, match="No national records"):
            load_test_positivity(path)


class TestLoadSymptomaticIllnesses:
    def test_loads_season_estimates(self, tmp_path):
        path = tmp_path / "burden.csv"
        pd.DataFrame(
            {
                "Flu Season": ["2020-2021", "2021-2022"],
                "Symptomatic Illnesses": ["1000000", "2000000"],
                "Other column": [1, 2],
            }
        ).to_csv(path, index=False)

        result = load_symptomatic_illnesses(path)

        assert result.columns.tolist() == [
            "flu_season",
            "symptomatic_illnesses",
        ]
        assert result["flu_season"].tolist() == ["2020-2021", "2021-2022"]
        np.testing.assert_array_equal(
            result["symptomatic_illnesses"].to_numpy(),
            [1_000_000, 2_000_000],
        )

    def test_raises_for_non_numeric_estimates(self, tmp_path):
        path = tmp_path / "burden.csv"
        pd.DataFrame(
            {
                "Flu Season": ["2020-2021"],
                "Symptomatic Illnesses": ["not a number"],
            }
        ).to_csv(path, index=False)

        with pytest.raises((ValueError, TypeError)):
            load_symptomatic_illnesses(path)


class TestLoadStrainFractions:
    @staticmethod
    def _write_strain_file(tmp_path, data):
        return write_cdc_csv(tmp_path, "strains.csv", data)

    def test_returns_expected_columns_and_fractions(self, tmp_path):
        path = self._write_strain_file(
            tmp_path,
            {
                "REGION TYPE": ["National"],
                "YEAR": [2021],
                "WEEK": [40],
                "A (2009 H1N1)": [10],
                "A (H3)": [30],
                "A (Subtyping not Performed)": [0],
                "B": [0],
                "BVic": [20],
                "BYam": [40],
            },
        )

        result = load_strain_fractions(path)

        assert result.columns.tolist() == [
            "year",
            "week",
            "a_h1n1",
            "a_h3",
            "b_victoria",
            "b_yamagata",
        ]
        np.testing.assert_allclose(
            result.loc[
                0, ["a_h1n1", "a_h3", "b_victoria", "b_yamagata"]
            ].to_numpy(dtype=float),
            [0.1, 0.3, 0.2, 0.4],
        )

    def test_distributes_untyped_cases_using_weekly_typed_shares(
        self, tmp_path
    ):
        path = self._write_strain_file(
            tmp_path,
            {
                "REGION TYPE": ["National"],
                "YEAR": [2021],
                "WEEK": [40],
                "A (2009 H1N1)": [1],
                "A (H3)": [3],
                "A (Subtyping not Performed)": [4],
                "B": [0],
                "BVic": [0],
                "BYam": [0],
            },
        )

        result = load_strain_fractions(path)

        # Die vier nicht typisierten A-Fälle werden im Verhältnis 1:3 verteilt.
        np.testing.assert_allclose(
            result.loc[0, ["a_h1n1", "a_h3"]].to_numpy(dtype=float),
            [0.25, 0.75],
        )

    def test_ignores_h5_and_h3n2v_cases(self, tmp_path):
        path = self._write_strain_file(
            tmp_path,
            {
                "REGION TYPE": ["National"],
                "YEAR": [2021],
                "WEEK": [40],
                "A (2009 H1N1)": [1],
                "A (H3)": [1],
                "A (Subtyping not Performed)": [0],
                "B": [0],
                "BVic": [0],
                "BYam": [0],
                "A (H5)": [100],
                "A (H3N2v)": [100],
            },
        )

        result = load_strain_fractions(path)

        np.testing.assert_allclose(
            result.loc[
                0, ["a_h1n1", "a_h3", "b_victoria", "b_yamagata"]
            ].to_numpy(dtype=float),
            [0.5, 0.5, 0.0, 0.0],
        )

    def test_rows_without_assignable_cases_have_missing_fractions(
        self, tmp_path
    ):
        path = self._write_strain_file(
            tmp_path,
            {
                "REGION TYPE": ["National"],
                "YEAR": [2021],
                "WEEK": [40],
                "A (2009 H1N1)": [0],
                "A (H3)": [0],
                "A (Subtyping not Performed)": [0],
                "B": [0],
                "BVic": [0],
                "BYam": [0],
            },
        )

        result = load_strain_fractions(path)

        assert (
            result[["a_h1n1", "a_h3", "b_victoria", "b_yamagata"]]
            .isna()
            .all(axis=None)
        )

    def test_raises_when_no_national_records_exist(self, tmp_path):
        path = self._write_strain_file(
            tmp_path,
            {
                "REGION TYPE": ["HHS Region"],
                "YEAR": [2021],
                "WEEK": [40],
                "A (2009 H1N1)": [1],
                "A (H3)": [1],
                "A (Subtyping not Performed)": [0],
                "B": [0],
                "BVic": [0],
                "BYam": [0],
            },
        )

        with pytest.raises(ValueError, match="No national records"):
            load_strain_fractions(path)
