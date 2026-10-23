"""
Estimates weekly symptomatic influenza cases per strain for the 2017-2018
influenza season and exports the results to a CSV file.

Weekly total cases are estimated from weighted ILI rates and influenza test
positivity. These cases are distributed among A(H1N1), A(H3), B/Victoria, and
B/Yamagata using weekly strain fractions from public health laboratory data.
"""

import numpy as np
import pandas as pd

from influenza_case_estimation import (
    load_strain_fractions,
    load_symptomatic_illnesses,
    load_test_positivity,
    load_weighted_ili_rate,
    weekly_cases_per_strain,
    weekly_influenza_activity,
    weekly_influenza_cases,
)


def is_in_flu_season(
    data: pd.DataFrame,
    flu_season: str,
) -> pd.Series:
    """
    Return a mask for epidemiological weeks in one influenza season.

    A season such as ``2017-2018`` starts in calendar week 40 of the first
    year and ends in calendar week 39 of the following year.

    Parameters
    ----------
    data : pandas.DataFrame
        Data frame containing the columns ``year`` and ``week``.
    flu_season : str
        Influenza season in the format ``"YYYY-YYYY"``.

    Returns
    -------
    pandas.Series
        Boolean mask identifying rows that belong to the requested season.

    Raises
    ------
    ValueError
        If ``flu_season`` does not describe two consecutive years.
    """
    try:
        start_year, end_year = map(int, flu_season.split("-"))
    except ValueError as error:
        raise ValueError(
            "flu_season must use the format 'YYYY-YYYY', "
            "for example '2017-2018'."
        ) from error

    if end_year != start_year + 1:
        raise ValueError(
            "A flu season must cover consecutive years, "
            "for example '2017-2018'."
        )

    return ((data["year"] == start_year) & (data["week"] >= 40)) | (
        (data["year"] == end_year) & (data["week"] <= 39)
    )


# File locations and influenza season
data_folder = "../data"

ili_file = f"{data_folder}/ILINet.csv"
positivity_file = f"{data_folder}/ICL_NREVSS_Clinical_Labs.csv"
public_health_labs_file = f"{data_folder}/ICL_NREVSS_Public_Health_Labs.csv"
burden_file = (
    f"{data_folder}/"
    "flu-disease-burden-past-season-rounded-estimates-all-ages.csv"
)

flu_season = "2017-2018"
output_file = f"{data_folder}/weekly_cases_per_strain_{flu_season}.csv"

strain_columns = [
    "a_h1n1",
    "a_h3",
    "b_victoria",
    "b_yamagata",
]

output_strain_columns = [
    "A(H1N1)",
    "A(H3)",
    "B/Victoria",
    "B/Yamagata",
]


# Load surveillance and burden data
ili_rates = load_weighted_ili_rate(ili_file)
test_positivity = load_test_positivity(positivity_file)
strain_fractions = load_strain_fractions(public_health_labs_file)
symptomatic_illnesses = load_symptomatic_illnesses(burden_file)

season_total = symptomatic_illnesses.loc[
    symptomatic_illnesses["flu_season"].eq(flu_season),
    "symptomatic_illnesses",
]

if season_total.empty:
    raise ValueError(
        f"No symptomatic illness estimate found for {flu_season}."
    )

total_cases = int(round(season_total.iloc[0]))


# Combine weekly data for the selected season
weekly_data = ili_rates.merge(
    test_positivity,
    on=["year", "week"],
    how="inner",
    validate="one_to_one",
)

weekly_data = weekly_data.loc[is_in_flu_season(weekly_data, flu_season)].copy()

if weekly_data.empty:
    raise ValueError(f"No surveillance data found for season {flu_season}.")

weekly_data = weekly_data.merge(
    strain_fractions,
    on=["year", "week"],
    how="left",
    validate="one_to_one",
)

weekly_data = weekly_data.sort_values(["year", "week"]).reset_index(drop=True)

missing_strain_columns = set(strain_columns) - set(weekly_data.columns)

if missing_strain_columns:
    raise ValueError(
        "Strain columns are missing from the public health laboratory data: "
        f"{sorted(missing_strain_columns)}"
    )

if weekly_data[strain_columns].isna().any().any():
    missing_weeks = weekly_data.loc[
        weekly_data[strain_columns].isna().any(axis=1),
        ["year", "week"],
    ]

    raise ValueError(
        "No assignable strain fractions for the following weeks: "
        f"{missing_weeks.to_dict(orient='records')}"
    )

fraction_sums = weekly_data[strain_columns].sum(axis=1)

if not np.allclose(
    fraction_sums.to_numpy(),
    1.0,
    rtol=0.0,
    atol=1e-10,
):
    raise ValueError(
        "The four strain fractions do not sum to 1 in every week."
    )


# Estimate weekly cases and distribute them among strains
weekly_data["activity"] = weekly_influenza_activity(
    weekly_data["weighted_ili_rate"].to_numpy(),
    weekly_data["test_positivity"].to_numpy(),
)

weekly_data["estimated_cases"] = weekly_influenza_cases(
    total_cases,
    weekly_data["activity"].to_numpy(),
)

cases_per_strain = weekly_cases_per_strain(
    weekly_data["estimated_cases"].to_numpy(),
    weekly_data[strain_columns].to_numpy(dtype=float),
)


# Create and save CSV output
result = weekly_data[["year", "week", "estimated_cases"]].copy()

result["calendar_week"] = (
    result["year"].astype(str) + "-W" + result["week"].astype(str).str.zfill(2)
)

result[output_strain_columns] = cases_per_strain

if not np.allclose(
    result[output_strain_columns].sum(axis=1).to_numpy(),
    result["estimated_cases"].to_numpy(),
    rtol=0.0,
    atol=1e-6,
):
    raise ValueError(
        "The four strain-specific case counts do not sum to total cases."
    )

result = result[
    [
        "year",
        "week",
        "calendar_week",
        "estimated_cases",
        *output_strain_columns,
    ]
]

result.to_csv(
    output_file,
    index=False,
    float_format="%.6f",
)

print(f"CSV file saved: {output_file}")
print(f"Total estimated cases: {result['estimated_cases'].sum():,}")

for strain in output_strain_columns:
    print(f"{strain}: {result[strain].sum():,.0f}")
