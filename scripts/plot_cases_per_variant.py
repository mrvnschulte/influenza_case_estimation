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


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIRECTORY = PROJECT_ROOT / "data"
PLOTS_DIRECTORY = PROJECT_ROOT / "plots"

ILI_FILE = DATA_DIRECTORY / "ILINet.csv"
POSITIVITY_FILE = DATA_DIRECTORY / "ICL_NREVSS_Clinical_Labs.csv"
PUBLIC_HEALTH_LABS_FILE = (
    DATA_DIRECTORY / "ICL_NREVSS_Public_Health_Labs.csv"
)
BURDEN_FILE = (
    DATA_DIRECTORY
    / "flu-disease-burden-past-season-rounded-estimates-all-ages.csv"
)

FLU_SEASONS = [
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

STRAIN_LABELS = {
    "a_h1n1pdm09": "A(H1N1)pdm09",
    "a_h3": "A(H3)",
    "a_unsubtyped": "A (not subtyped)",
    "b": "B",
    "b_victoria": "B/Victoria",
    "b_yamagata": "B/Yamagata",
    "h3n2v": "H3N2v",
    "a_h5": "A(H5)",
}

STRAIN_COLORS = (
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
)


def is_in_flu_season(data: pd.DataFrame, flu_season: str) -> pd.Series:
    """
    Return a mask for epidemiological weeks in a flu season.

    A season such as ``2024-2025`` starts in calendar week 40 of the first
    year and ends in calendar week 39 of the following year.
    """
    try:
        start_year, end_year = map(int, flu_season.split("-"))
    except ValueError as error:
        raise ValueError(
            "flu_season must use the format 'YYYY-YYYY', "
            "for example '2024-2025'."
        ) from error

    if end_year != start_year + 1:
        raise ValueError(
            "A flu season must cover consecutive years, "
            "for example '2024-2025'."
        )

    return (
        ((data["year"] == start_year) & (data["week"] >= 40))
        | ((data["year"] == end_year) & (data["week"] <= 39))
    )


def plot_flu_season(
    flu_season: str,
    ili_rates: pd.DataFrame,
    test_positivity: pd.DataFrame,
    strain_fractions: pd.DataFrame,
    symptomatic_illnesses: pd.DataFrame,
) -> None:
    """Estimate weekly cases per strain and save a plot for one flu season."""
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

    # Fehlende Variantendaten sollen keine Wochen unbemerkt entfernen.
    weekly_data = weekly_data.merge(
        strain_fractions,
        on=["year", "week"],
        how="left",
        validate="one_to_one",
    )

    weekly_data = weekly_data.sort_values(
        ["year", "week"]
    ).reset_index(drop=True)

    strain_columns = list(STRAIN_LABELS)

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
            "Strain data from public health labs is missing for at least one "
            f"week: {missing_weeks.to_dict(orient='records')}"
        )

    weekly_data["activity"] = weekly_influenza_activity(
        weekly_data["weighted_ili_rate"],
        weekly_data["test_positivity"],
    )

    weekly_data["estimated_cases"] = weekly_influenza_cases(
        total_cases,
        weekly_data["activity"].to_numpy(),
    )

    weekly_cases_by_strain = pd.DataFrame(
        weekly_cases_per_strain(
            weekly_data["estimated_cases"].to_numpy(),
            weekly_data[strain_columns].to_numpy(dtype=float),
        ),
        columns=strain_columns,
        index=weekly_data.index,
    )

    if not np.allclose(
        weekly_cases_by_strain.sum(axis=1).to_numpy(),
        weekly_data["estimated_cases"].to_numpy(),
        rtol=0.0,
        atol=1e-6,
    ):
        raise ValueError(
            "The sum of cases per strain does not equal the total cases."
        )

    weekly_data["calendar_week"] = (
        weekly_data["year"].astype(str)
        + "-W"
        + weekly_data["week"].astype(str).str.zfill(2)
    )

    PLOTS_DIRECTORY.mkdir(parents=True, exist_ok=True)

    x = weekly_data.index.to_numpy()

    figure, axis = plt.subplots(figsize=(12, 6))

    axis.stackplot(
        x,
        *weekly_cases_by_strain.to_numpy().T,
        labels=list(STRAIN_LABELS.values()),
        colors=STRAIN_COLORS,
        alpha=0.85,
        baseline="zero",
    )

    # Oberer Rand der Stapelflächen entspricht der geschätzten Gesamtfallzahl.
    axis.plot(
        x,
        weekly_data["estimated_cases"],
        color="black",
        linewidth=1.1,
        label="Total estimated cases",
    )

    axis.set_title(
        "Estimated Weekly Symptomatic Influenza Cases per Strain, "
        f"{flu_season} Season"
    )
    axis.set_xlabel("Calendar week")
    axis.set_ylabel("Estimated symptomatic influenza cases")
    axis.set_ylim(bottom=0)
    axis.grid(axis="y", alpha=0.3)
    axis.ticklabel_format(style="plain", axis="y")
    axis.legend(title="Strain", loc="upper left", ncol=2)

    tick_positions = x[::4]
    axis.set_xticks(tick_positions)
    axis.set_xticklabels(
        weekly_data.iloc[tick_positions]["calendar_week"],
        rotation=45,
        ha="right",
    )

    output_file = (
        PLOTS_DIRECTORY
        / f"weekly_influenza_cases_by_strain_{flu_season}.png"
    )

    figure.tight_layout()
    figure.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"Plot saved to: {output_file}")
    print(
        "Allocated total cases: "
        f"{weekly_data['estimated_cases'].sum():,}"
    )


def main() -> None:
    """Load all data once and create plots for all configured flu seasons."""
    ili_rates = load_weighted_ili_rate(ILI_FILE)
    test_positivity = load_test_positivity(POSITIVITY_FILE)
    strain_fractions = load_strain_fractions(PUBLIC_HEALTH_LABS_FILE)
    symptomatic_illnesses = load_symptomatic_illnesses(BURDEN_FILE)

    for flu_season in FLU_SEASONS:
        print(f"\nProcessing season: {flu_season}")

        try:
            plot_flu_season(
                flu_season=flu_season,
                ili_rates=ili_rates,
                test_positivity=test_positivity,
                strain_fractions=strain_fractions,
                symptomatic_illnesses=symptomatic_illnesses,
            )
        except ValueError as error:
            print(f"Skipped {flu_season}: {error}")


if __name__ == "__main__":
    main()