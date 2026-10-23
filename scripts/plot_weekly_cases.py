"""
Estimate weekly symptomatic influenza cases and create a line plot for one
influenza season.

The seasonal total number of symptomatic influenza illnesses is distributed
across epidemiological weeks according to the product of weighted ILI rates
and influenza test positivity.
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


# Plot settings
dpi = 1200
fileformat = "pdf"
folder = "../plots"
fontsize_text = 20
fontsize_ticks = 16
alpha_grid = 0.7
y_axis_exponent = 6

data_folder = "../data"

flu_season = "2017-2018"

ili_file = f"{data_folder}/ILINet.csv"
positivity_file = f"{data_folder}/ICL_NREVSS_Clinical_Labs.csv"
burden_file = (
    f"{data_folder}/"
    "flu-disease-burden-past-season-rounded-estimates-all-ages.csv"
)


# Load and prepare data
Path(folder).mkdir(parents=True, exist_ok=True)

ili_rates = load_weighted_ili_rate(ili_file)
test_positivity = load_test_positivity(positivity_file)
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

weekly_data = ili_rates.merge(
    test_positivity,
    on=["year", "week"],
    how="inner",
    validate="one_to_one",
)

weekly_data = weekly_data.loc[is_in_flu_season(weekly_data, flu_season)].copy()

if weekly_data.empty:
    raise ValueError(f"No surveillance data found for season {flu_season}.")

weekly_data = weekly_data.sort_values(["year", "week"]).reset_index(drop=True)

weekly_data["activity"] = weekly_influenza_activity(
    weekly_data["weighted_ili_rate"].to_numpy(),
    weekly_data["test_positivity"].to_numpy(),
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


# Create plot
x = weekly_data.index.to_numpy()

plt.figure()

plt.plot(
    x,
    weekly_data["estimated_cases"],
    color="tab:blue",
)

plt.fill_between(
    x,
    weekly_data["estimated_cases"],
    color="tab:blue",
    alpha=0.2,
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

plt.tight_layout()

output_file = f"{folder}/weekly_influenza_cases_{flu_season}.{fileformat}"

plt.savefig(
    output_file,
    dpi=dpi,
    bbox_inches="tight",
)
plt.show()
plt.close()

print(f"Plot saved to: {output_file}")
print("Allocated total cases: " f"{weekly_data['estimated_cases'].sum():,}")
