"""
Estimate weekly symptomatic influenza cases per variant and create stacked
plots for multiple influenza seasons.

The total number of symptomatic illnesses is distributed across weeks using
the product of weighted ILI rates and influenza test positivity. Weekly totals
are then distributed among A(H1N1), A(H3), B/Victoria, and B/Yamagata using
the strain fractions from public-health laboratory surveillance data.
"""

from pathlib import Path

import matplotlib.pyplot as plt
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
        Boolean mask identifying rows belonging to the requested season.

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


def plot_flu_season(
    flu_season: str,
    ili_rates: pd.DataFrame,
    test_positivity: pd.DataFrame,
    strain_fractions: pd.DataFrame,
    symptomatic_illnesses: pd.DataFrame,
) -> None:
    """
    Estimate weekly cases per variant and create a stacked area plot.

    Parameters
    ----------
    flu_season : str
        Influenza season in the format ``"YYYY-YYYY"``.
    ili_rates : pandas.DataFrame
        National weekly weighted ILI rates.
    test_positivity : pandas.DataFrame
        National weekly influenza test positivity rates.
    strain_fractions : pandas.DataFrame
        National weekly fractions for the four model strains.
    symptomatic_illnesses : pandas.DataFrame
        Estimated total symptomatic influenza illnesses by season.

    Raises
    ------
    ValueError
        If required burden, surveillance, or strain-fraction data is missing.
    """
    season_total = symptomatic_illnesses.loc[
        symptomatic_illnesses["flu_season"].eq(flu_season),
        "symptomatic_illnesses",
    ]

    if season_total.empty:
        raise ValueError(
            f"No symptomatic illness estimate found for {flu_season}."
        )

    total_cases = int(round(season_total.iloc[0]))

    weekly_data = ili_rates.merge(
        test_positivity,
        on=["year", "week"],
        how="inner",
        validate="one_to_one",
    )

    weekly_data = weekly_data.loc[
        is_in_flu_season(weekly_data, flu_season)
    ].copy()

    if weekly_data.empty:
        raise ValueError(
            f"No surveillance data found for season {flu_season}."
        )

    weekly_data = weekly_data.merge(
        strain_fractions,
        on=["year", "week"],
        how="left",
        validate="one_to_one",
    )

    weekly_data = weekly_data.sort_values(["year", "week"]).reset_index(
        drop=True
    )

    missing_strain_columns = set(strain_columns) - set(weekly_data.columns)

    if missing_strain_columns:
        raise ValueError(
            "Strain columns are missing from the public health laboratory "
            f"data: {sorted(missing_strain_columns)}"
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

    weekly_data["activity"] = weekly_influenza_activity(
        weekly_data["weighted_ili_rate"].to_numpy(),
        weekly_data["test_positivity"].to_numpy(),
    )

    weekly_data["estimated_cases"] = weekly_influenza_cases(
        total_cases,
        weekly_data["activity"].to_numpy(),
    )

    weekly_cases_by_strain = weekly_cases_per_strain(
        weekly_data["estimated_cases"].to_numpy(),
        weekly_data[strain_columns].to_numpy(dtype=float),
    )

    if not np.allclose(
        weekly_cases_by_strain.sum(axis=1),
        weekly_data["estimated_cases"].to_numpy(),
        rtol=0.0,
        atol=1e-6,
    ):
        raise ValueError(
            "The sum of variant-specific cases does not equal total cases."
        )

    weekly_data["calendar_week"] = (
        weekly_data["year"].astype(str)
        + "-W"
        + weekly_data["week"].astype(str).str.zfill(2)
    )

    x = np.arange(len(weekly_data))

    fig = plt.figure()

    plt.stackplot(
        x,
        *weekly_cases_by_strain.T,
        labels=strain_labels,
        colors=strain_colors,
        alpha=0.85,
    )

    plt.xlabel("Calendar week", fontsize=fontsize_text)
    plt.ylabel(
        "Estimated cases",
        fontsize=fontsize_text,
    )

    plt.xlim(x[0], x[-1])
    plt.ylim(bottom=0)
    plt.xticks(
        x[::8],
        weekly_data.loc[x[::8], "calendar_week"],
        rotation=90,
        fontsize=fontsize_ticks,
    )
    plt.yticks(fontsize=fontsize_ticks)
    plt.ticklabel_format(
        axis="y",
        style="sci",
        scilimits=(y_axis_exponent, y_axis_exponent),
        useMathText=True,
    )
    plt.gca().yaxis.get_offset_text().set_fontsize(fontsize_ticks)
    plt.grid(axis="y", alpha=alpha_grid)

    fig.legend(
        title="Strain",
        title_fontproperties={
            "weight": "bold",
            "size": fontsize_text,
        },
        fontsize=fontsize_ticks,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.05),
        ncol=2,
    )

    plt.tight_layout()

    output_file = (
        f"{folder}/weekly_influenza_cases_by_strain_{flu_season}."
        f"{fileformat}"
    )

    plt.savefig(
        output_file,
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.show()
    plt.close()

    print(f"Plot saved to: {output_file}")
    print(
        "Allocated total cases: " f"{weekly_data['estimated_cases'].sum():,}"
    )


# Plot settings
dpi = 1200
fileformat = "pdf"
folder = "../plots"
fontsize_text = 20
fontsize_ticks = 16
alpha_grid = 0.7
y_axis_exponent = 6

data_folder = "../data"

flu_seasons = [
    "2015-2016",
    "2016-2017",
    "2017-2018",
    "2018-2019",
    "2019-2020",
    "2020-2021",
    "2021-2022",
    "2022-2023",
    "2023-2024",
    "2024-2025",
]

strain_columns = [
    "a_h1n1",
    "a_h3",
    "b_victoria",
    "b_yamagata",
]

strain_labels = [
    "A(H1N1)",
    "A(H3)",
    "B/Victoria",
    "B/Yamagata",
]

strain_colors = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
]

ili_file = f"{data_folder}/ILINet.csv"
positivity_file = f"{data_folder}/ICL_NREVSS_Clinical_Labs.csv"
public_health_labs_file = f"{data_folder}/ICL_NREVSS_Public_Health_Labs.csv"
burden_file = (
    f"{data_folder}/"
    "flu-disease-burden-past-season-rounded-estimates-all-ages.csv"
)


# Load data and create plots
Path(folder).mkdir(parents=True, exist_ok=True)

ili_rates = load_weighted_ili_rate(ili_file)
test_positivity = load_test_positivity(positivity_file)
strain_fractions = load_strain_fractions(public_health_labs_file)
symptomatic_illnesses = load_symptomatic_illnesses(burden_file)

for flu_season in flu_seasons:
    print()
    print(f"Processing season: {flu_season}")

    try:
        plot_flu_season(
            flu_season,
            ili_rates,
            test_positivity,
            strain_fractions,
            symptomatic_illnesses,
        )
    except ValueError as error:
        print(f"Skipped {flu_season}: {error}")
