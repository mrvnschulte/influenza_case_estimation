"""
Functions for loading and preparing national influenza surveillance and burden
data from CDC CSV files.
"""

from pathlib import Path

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
        raise ValueError("No national records found in the clinical labs file.")

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

    result = burden_data[
        ["Flu Season", "Symptomatic Illnesses"]
    ].copy()
    result.columns = ["flu_season", "symptomatic_illnesses"]

    result["symptomatic_illnesses"] = pd.to_numeric(
        result["symptomatic_illnesses"],
        errors="raise",
    )

    return result


def load_strain_fractions(file_path: str | Path) -> pd.DataFrame:
    """
    Load national weekly influenza subtype and lineage fractions.

    Fractions are calculated from the counts of positive laboratory tests
    reported for each influenza subtype or lineage. Each row represents one
    epidemiological week, and each strain column contains its fraction of all
    reported positive tests in that week.

    The returned strain categories are:

    - ``a_h1n1pdm09``
    - ``a_h3``
    - ``a_unsubtyped``
    - ``b``
    - ``b_victoria``
    - ``b_yamagata``
    - ``h3n2v``
    - ``a_h5``

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the ``ICL_NREVSS_Public_Health_Labs.csv`` file.

    Returns
    -------
    pandas.DataFrame
        Data frame with columns ``year`` and ``week``, followed by one column
        per influenza subtype or lineage. Strain columns contain decimal
        fractions between 0 and 1. Rows without positive tests are assigned
        missing values for all strain fractions.

    Raises
    ------
    ValueError
        If no national records are found in the input file.
    """
    laboratory_data = pd.read_csv(file_path, skiprows=1)

    strain_columns = {
        "A (2009 H1N1)": "a_h1n1pdm09",
        "A (H3)": "a_h3",
        "A (Subtyping not Performed)": "a_unsubtyped",
        "B": "b",
        "BVic": "b_victoria",
        "BYam": "b_yamagata",
        "H3N2v": "h3n2v",
        "A (H5)": "a_h5",
    }

    national_data = laboratory_data.loc[
        laboratory_data["REGION TYPE"].eq("National"),
        ["YEAR", "WEEK", *strain_columns],
    ].copy()

    if national_data.empty:
        raise ValueError(
            "No national records found in the public health laboratories file."
        )

    national_data = national_data.rename(
        columns={
            "YEAR": "year",
            "WEEK": "week",
            **strain_columns,
        }
    )

    fraction_columns = list(strain_columns.values())

    for column in fraction_columns:
        national_data[column] = pd.to_numeric(
            national_data[column],
            errors="raise",
        )

    total_positive_tests = national_data[fraction_columns].sum(axis=1)

    national_data[fraction_columns] = national_data[fraction_columns].div(
        total_positive_tests.replace(0, pd.NA),
        axis=0,
    )

    return national_data.sort_values(["year", "week"]).reset_index(drop=True)