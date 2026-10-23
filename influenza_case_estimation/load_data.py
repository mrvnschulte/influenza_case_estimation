"""
Functions for loading and preparing national influenza surveillance and burden
data from CDC CSV files.
"""

from pathlib import Path

import numpy as np
import pandas as pd


def load_weighted_ili_rate(file_path: str | Path) -> pd.DataFrame:
    """
    Load national weekly weighted ILI rates.

    The source file reports weighted ILI as a percentage. This function
    converts it to a decimal proportion, for example 1.90509 becomes
    0.0190509.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the ``ILINet.csv`` file.

    Returns
    -------
    pandas.DataFrame
        Data frame with the columns ``year``, ``week``, and
        ``weighted_ili_rate``. The rate is represented as a decimal
        proportion rather than a percentage.

    Raises
    ------
    ValueError
        If no national records are found in the input file.
    """
    ili_data = pd.read_csv(file_path, skiprows=1)

    national_data = ili_data.loc[
        ili_data["REGION TYPE"].eq("National"),
        ["YEAR", "WEEK", "% WEIGHTED ILI"],
    ].copy()

    if national_data.empty:
        raise ValueError("No national records found in the ILINet file.")

    national_data.columns = ["year", "week", "weighted_ili_rate"]
    national_data["weighted_ili_rate"] = (
        pd.to_numeric(national_data["weighted_ili_rate"], errors="raise") / 100
    )

    return national_data.sort_values(["year", "week"]).reset_index(drop=True)


def load_test_positivity(file_path: str | Path) -> pd.DataFrame:
    """
    Load national weekly influenza test positivity rates from clinical labs.

    The source file reports the overall positivity rate in percent. This
    function converts it to a decimal proportion, for example 18.2344 becomes
    0.182344.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the ``ICL_NREVSS_Clinical_Labs.csv`` file.

    Returns
    -------
    pandas.DataFrame
        Data frame with the columns ``year``, ``week``, and
        ``test_positivity``. The positivity rate is represented as a decimal
        proportion between 0 and 1.

    Raises
    ------
    ValueError
        If no national records are found, the positivity values are not
        numeric, or converted positivity rates are outside the inclusive
        range from 0 to 1.

    Examples
    --------
    >>> positivity = load_test_positivity("ICL_NREVSS_Clinical_Labs.csv")
    >>> positivity.loc[0, "test_positivity"]
    0.00674072
    """
    clinical_labs = pd.read_csv(file_path, skiprows=1)

    national_data = clinical_labs.loc[
        clinical_labs["REGION TYPE"].eq("National"),
        ["YEAR", "WEEK", "PERCENT POSITIVE"],
    ].copy()

    if national_data.empty:
        raise ValueError(
            "No national records found in the clinical labs file."
        )

    national_data.columns = ["year", "week", "test_positivity"]
    national_data["test_positivity"] = (
        pd.to_numeric(national_data["test_positivity"], errors="raise") / 100
    )

    if not national_data["test_positivity"].between(0, 1).all():
        raise ValueError(
            "Converted test positivity values must be between 0 and 1."
        )

    return national_data.sort_values(["year", "week"]).reset_index(drop=True)


def load_symptomatic_illnesses(file_path: str | Path) -> pd.DataFrame:
    """
    Load estimated total symptomatic influenza illnesses by flu season.

    The estimates correspond to the point estimates published by the CDC.
    Uncertainty intervals included in the source file are not returned.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the
        ``flu-disease-burden-past-season-rounded-estimates-all-ages.csv` file.

    Returns
    -------
    pandas.DataFrame
        Data frame with the columns ``flu_season`` and
        ``symptomatic_illnesses``. Each row contains the estimated total
        number of symptomatic influenza illnesses for one flu season.

    Raises
    ------
    ValueError
        If the input file contains missing or non-numeric illness estimates.
    """
    burden_data = pd.read_csv(file_path)

    result = burden_data[["Flu Season", "Symptomatic Illnesses"]].copy()
    result.columns = ["flu_season", "symptomatic_illnesses"]

    result["symptomatic_illnesses"] = pd.to_numeric(
        result["symptomatic_illnesses"],
        errors="raise",
    )

    return result


# def load_strain_fractions(file_path: str | Path) -> pd.DataFrame:
#     """
#     Load national weekly influenza subtype and lineage fractions.

#     Fractions are calculated from the counts of positive laboratory tests
#     reported for each influenza subtype or lineage. Each row represents one
#     epidemiological week, and each strain column contains its fraction of all
#     reported positive tests in that week.

#     The returned strain categories are:

#     - ``a_h1n1pdm09``
#     - ``a_h3``
#     - ``a_unsubtyped``
#     - ``b``
#     - ``b_victoria``
#     - ``b_yamagata``
#     - ``h3n2v``
#     - ``a_h5``

#     Parameters
#     ----------
#     file_path : str or pathlib.Path
#         Path to the ``ICL_NREVSS_Public_Health_Labs.csv`` file.

#     Returns
#     -------
#     pandas.DataFrame
#         Data frame with columns ``year`` and ``week``, followed by one column
#         per influenza subtype or lineage. Strain columns contain decimal
#         fractions between 0 and 1. Rows without positive tests are assigned
#         missing values for all strain fractions.

#     Raises
#     ------
#     ValueError
#         If no national records are found in the input file.
#     """
#     laboratory_data = pd.read_csv(file_path, skiprows=1)

#     strain_columns = {
#         "A (2009 H1N1)": "a_h1n1pdm09",
#         "A (H3)": "a_h3",
#         "A (Subtyping not Performed)": "a_unsubtyped",
#         "B": "b",
#         "BVic": "b_victoria",
#         "BYam": "b_yamagata",
#         "H3N2v": "h3n2v",
#         "A (H5)": "a_h5",
#     }

#     national_data = laboratory_data.loc[
#         laboratory_data["REGION TYPE"].eq("National"),
#         ["YEAR", "WEEK", *strain_columns],
#     ].copy()

#     if national_data.empty:
#         raise ValueError(
#             "No national records found in the public health laboratories file."
#         )

#     national_data = national_data.rename(
#         columns={
#             "YEAR": "year",
#             "WEEK": "week",
#             **strain_columns,
#         }
#     )

#     fraction_columns = list(strain_columns.values())

#     for column in fraction_columns:
#         national_data[column] = pd.to_numeric(
#             national_data[column],
#             errors="raise",
#         )

#     total_positive_tests = national_data[fraction_columns].sum(axis=1)

#     national_data[fraction_columns] = national_data[fraction_columns].div(
#         total_positive_tests.replace(0, pd.NA),
#         axis=0,
#     )

#     return national_data.sort_values(["year", "week"]).reset_index(drop=True)


def load_strain_fractions(file_path: str | Path) -> pd.DataFrame:
    """
    Load national weekly fractions for the four influenza model strains.

    The source file contains laboratory-confirmed influenza-positive tests by
    subtype and lineage. This function returns fractions for A(H1N1), A(H3),
    B/Victoria, and B/Yamagata only.

    Unsubtyped influenza-A cases are distributed between A(H1N1) and A(H3)
    according to their proportions among typed influenza-A cases in the same
    epidemiological week. Similarly, influenza-B cases without lineage
    information are distributed between B/Victoria and B/Yamagata according
    to their proportions among lineaged influenza-B cases in that week.

    If a week contains unassigned A or B cases but no corresponding typed
    subtype or lineage cases, all unassigned cases are allocated to the
    subtype or lineage with the larger total number of typed cases in the
    complete influenza season. Ties are resolved in favour of A(H1N1) and
    B/Victoria, respectively.

    H3N2v and A(H5) cases are ignored and do not contribute to the returned
    fractions or their denominators.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the ``ICL_NREVSS_Public_Health_Labs.csv`` file.

    Returns
    -------
    pandas.DataFrame
        Data frame with the columns ``year``, ``week``, ``a_h1n1``, ``a_h3``,
        ``b_victoria``, and ``b_yamagata``. The strain columns contain decimal
        fractions between 0 and 1 and sum to 1 for weeks with at least one
        assigned case. Rows without assignable cases contain missing values
        for all strain fractions.

    Raises
    ------
    ValueError
        If no national records are found in the input file or if required
        count columns contain non-numeric values.

    Examples
    --------
    >>> fractions = load_strain_fractions(
    ...     "ICL_NREVSS_Public_Health_Labs.csv"
    ... )
    >>> fractions.columns.tolist()
    ['year', 'week', 'a_h1n1', 'a_h3', 'b_victoria', 'b_yamagata']
    """
    laboratory_data = pd.read_csv(file_path, skiprows=1)

    source_columns = {
        "A (2009 H1N1)": "a_h1n1",
        "A (H3)": "a_h3",
        "A (Subtyping not Performed)": "a_unsubtyped",
        "B": "b_unlineaged",
        "BVic": "b_victoria",
        "BYam": "b_yamagata",
    }

    national_data = laboratory_data.loc[
        laboratory_data["REGION TYPE"].eq("National"),
        ["YEAR", "WEEK", *source_columns],
    ].copy()

    if national_data.empty:
        raise ValueError(
            "No national records found in the public health laboratories file."
        )

    national_data = national_data.rename(
        columns={
            "YEAR": "year",
            "WEEK": "week",
            **source_columns,
        }
    )

    count_columns = list(source_columns.values())
    national_data[count_columns] = national_data[count_columns].apply(
        pd.to_numeric,
        errors="raise",
    )

    season_start_year = np.where(
        national_data["week"] >= 40,
        national_data["year"],
        national_data["year"] - 1,
    )

    national_data["flu_season"] = (
        season_start_year.astype(str)
        + "-"
        + (season_start_year + 1).astype(str)
    )

    seasonal_totals = national_data.groupby("flu_season")[
        ["a_h1n1", "a_h3", "b_victoria", "b_yamagata"]
    ].transform("sum")

    # Bei Gleichstand wird A(H1N1) bzw. B/Victoria bevorzugt.
    national_data["season_prefers_h1n1"] = (
        seasonal_totals["a_h1n1"] >= seasonal_totals["a_h3"]
    )
    national_data["season_prefers_victoria"] = (
        seasonal_totals["b_victoria"] >= seasonal_totals["b_yamagata"]
    )

    typed_a_total = national_data["a_h1n1"] + national_data["a_h3"]
    typed_b_total = national_data["b_victoria"] + national_data["b_yamagata"]

    h1n1_weekly_share = np.divide(
        national_data["a_h1n1"],
        typed_a_total,
        out=np.zeros(len(national_data), dtype=float),
        where=typed_a_total.to_numpy() > 0,
    )

    victoria_weekly_share = np.divide(
        national_data["b_victoria"],
        typed_b_total,
        out=np.zeros(len(national_data), dtype=float),
        where=typed_b_total.to_numpy() > 0,
    )

    h1n1_share = np.where(
        typed_a_total > 0,
        h1n1_weekly_share,
        national_data["season_prefers_h1n1"].astype(float),
    )

    victoria_share = np.where(
        typed_b_total > 0,
        victoria_weekly_share,
        national_data["season_prefers_victoria"].astype(float),
    )

    result = national_data[["year", "week"]].copy()

    result["a_h1n1"] = (
        national_data["a_h1n1"] + national_data["a_unsubtyped"] * h1n1_share
    )
    result["a_h3"] = national_data["a_h3"] + national_data["a_unsubtyped"] * (
        1 - h1n1_share
    )
    result["b_victoria"] = (
        national_data["b_victoria"]
        + national_data["b_unlineaged"] * victoria_share
    )
    result["b_yamagata"] = national_data["b_yamagata"] + national_data[
        "b_unlineaged"
    ] * (1 - victoria_share)

    model_strains = [
        "a_h1n1",
        "a_h3",
        "b_victoria",
        "b_yamagata",
    ]

    model_total = result[model_strains].sum(axis=1)

    result[model_strains] = result[model_strains].div(
        model_total.replace(0, np.nan),
        axis=0,
    )

    return result.sort_values(["year", "week"]).reset_index(drop=True)
