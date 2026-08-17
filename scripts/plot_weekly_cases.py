"""
Create a plot of estimated weekly symptomatic influenza cases for the
2024-2025 influenza season.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from influenza_case_estimation import (
    load_symptomatic_illnesses,
    load_test_positivity,
    load_weighted_ili_rate,
    weekly_influenza_activity,
    weekly_influenza_cases,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIRECTORY = PROJECT_ROOT / "data"
PLOTS_DIRECTORY = PROJECT_ROOT / "plots"

ILI_FILE = DATA_DIRECTORY / "ILINet.csv"
POSITIVITY_FILE = DATA_DIRECTORY / "ICL_NREVSS_Clinical_Labs.csv"
BURDEN_FILE = (
    DATA_DIRECTORY
    / "flu-disease-burden-past-season-rounded-estimates-all-ages.csv"
)

FLU_SEASON = "2024-2025"


def is_in_2024_2025_season(data: pd.DataFrame) -> pd.Series:
    """
    Return a mask for epidemiological weeks in the 2024-2025 flu season.

    The season starts with calendar week 40 of 2024 and ends with calendar
    week 39 of 2025.
    """
    return (
        ((data["year"] == 2024) & (data["week"] >= 40))
        | ((data["year"] == 2025) & (data["week"] <= 39))
    )


def main() -> None:
    """Load surveillance data, estimate weekly cases, and save a line plot."""
    ili_rates = load_weighted_ili_rate(ILI_FILE)
    test_positivity = load_test_positivity(POSITIVITY_FILE)
    symptomatic_illnesses = load_symptomatic_illnesses(BURDEN_FILE)

    season_total = symptomatic_illnesses.loc[
        symptomatic_illnesses["flu_season"].eq(FLU_SEASON),
        "symptomatic_illnesses",
    ]

    if season_total.empty:
        raise ValueError(f"No symptomatic illness estimate found for {FLU_SEASON}.")

    total_cases = int(round(season_total.iloc[0]))

    weekly_data = ili_rates.merge(
        test_positivity,
        on=["year", "week"],
        how="inner",
        validate="one_to_one",
    )

    weekly_data = weekly_data.loc[
        is_in_2024_2025_season(weekly_data)
    ].copy()

    if weekly_data.empty:
        raise ValueError(f"No surveillance data found for season {FLU_SEASON}.")

    weekly_data = weekly_data.sort_values(["year", "week"]).reset_index(drop=True)

    weekly_data["activity"] = weekly_influenza_activity(
        weekly_data["weighted_ili_rate"],
        weekly_data["test_positivity"],
    )

    weekly_data["estimated_cases"] = weekly_influenza_cases(
        total_cases,
        weekly_data["activity"].to_numpy(),
    )

    weekly_data["calendar_week"] = (
        weekly_data["year"].astype(str)
        + "-W"
        + weekly_data["week"].astype(str).str.zfill(2)
    )

    PLOTS_DIRECTORY.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(12, 6))

    axis.plot(
        weekly_data["calendar_week"],
        weekly_data["estimated_cases"],
        color="tab:blue",
        linewidth=2,
    )
    axis.fill_between(
        weekly_data["calendar_week"],
        weekly_data["estimated_cases"],
        color="tab:blue",
        alpha=0.2,
    )

    axis.set_title(
        "Estimated Weekly Symptomatic Influenza Cases, 2024–2025 Season"
    )
    axis.set_xlabel("Calendar week")
    axis.set_ylabel("Estimated symptomatic influenza cases")
    axis.grid(axis="y", alpha=0.3)
    axis.ticklabel_format(style="plain", axis="y")

    tick_positions = weekly_data.index[::4]
    axis.set_xticks(tick_positions)
    axis.set_xticklabels(
        weekly_data.loc[tick_positions, "calendar_week"],
        rotation=45,
        ha="right",
    )

    output_file = PLOTS_DIRECTORY / "weekly_influenza_cases_2024_2025.png"

    figure.tight_layout()
    figure.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"Plot saved to: {output_file}")
    print(f"Allocated total cases: {weekly_data['estimated_cases'].sum():,}")


if __name__ == "__main__":
    main()





