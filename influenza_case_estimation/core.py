"""
Core functions for estimating weekly influenza activity, case counts, and
strain-specific case distributions from surveillance data.
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray


def weekly_influenza_activity(
    ili_rate: ArrayLike,
    test_positivity: ArrayLike,
) -> NDArray[np.floating]:
    """
    Calculate weekly influenza activity from ILI rates and test positivity.

    Weekly activity is calculated as the element-wise product of the
    influenza-like illness (ILI) rate and the influenza test positivity rate.

    Parameters
    ----------
    ili_rate : ArrayLike
        One-dimensional array of weekly ILI rates. Values must be finite and
        non-negative.
    test_positivity : ArrayLike
        One-dimensional array of weekly influenza test positivity rates.
        Values must be finite and non-negative. Must have the same length as
        `ili_rate`.

    Returns
    -------
    numpy.ndarray
        One-dimensional floating-point array containing weekly influenza
        activity values.

    Raises
    ------
    ValueError
        If either input is not one-dimensional, inputs have different shapes,
        an input contains non-finite values, or an input contains negative
        values.

    Examples
    --------
    >>> weekly_influenza_activity([0.1, 0.2], [0.1, 0.25])
    array([0.01, 0.05])
    """
    ili_rate = np.asarray(ili_rate, dtype=float)
    test_positivity = np.asarray(test_positivity, dtype=float)

    if ili_rate.ndim != 1 or test_positivity.ndim != 1:
        raise ValueError("Inputs must be one-dimensional arrays.")
    if ili_rate.shape != test_positivity.shape:
        raise ValueError("ili_rate and test_positivity must have the same length.")
    if not np.all(np.isfinite(ili_rate)) or not np.all(np.isfinite(test_positivity)):
        raise ValueError("Inputs must contain only finite values.")
    if np.any(ili_rate < 0) or np.any(test_positivity < 0):
        raise ValueError("Inputs must not contain negative values.")

    return ili_rate * test_positivity


def weekly_influenza_cases(
    total_cases: int,
    weekly_activity: ArrayLike,
) -> NDArray[np.integer]:
    """
    Allocate a total number of influenza cases proportionally across weeks.

    The allocation is based on each week's share of total influenza activity.
    Cases are converted to integers using the largest remainder method:
    values are first rounded down, and remaining cases are assigned to weeks
    with the largest fractional remainders. This guarantees that the sum of
    returned weekly cases equals `total_cases`.

    Parameters
    ----------
    total_cases : int
        Total number of influenza cases to allocate. Must be non-negative.
    weekly_activity : ArrayLike
        One-dimensional array of weekly influenza activity values. Values must
        be finite and non-negative, and their sum must be greater than zero.

    Returns
    -------
    numpy.ndarray
        One-dimensional integer array containing the estimated number of
        influenza cases allocated to each week. Its sum equals `total_cases`.

    Raises
    ------
    ValueError
        If `total_cases` is negative, `weekly_activity` is not
        one-dimensional, contains non-finite or negative values, or has a
        sum of zero.

    Examples
    --------
    >>> weekly_influenza_cases(100, [0.1, 0.2, 0.3])
    array([17, 33, 50])
    """
    if total_cases < 0:
        raise ValueError("total_cases must not be negative.")

    weekly_activity = np.asarray(weekly_activity, dtype=float)

    if weekly_activity.ndim != 1:
        raise ValueError("weekly_activity must be one-dimensional.")
    if not np.all(np.isfinite(weekly_activity)) or np.any(weekly_activity < 0):
        raise ValueError("weekly_activity must contain finite, non-negative values.")

    activity_sum = weekly_activity.sum()
    if activity_sum == 0:
        raise ValueError("weekly_activity must have a sum greater than zero.")

    expected_cases = total_cases * weekly_activity / activity_sum
    cases = np.floor(expected_cases).astype(int)

    remainder = total_cases - cases.sum()
    
    if remainder > 0:
        largest_fractions = np.argsort(expected_cases - cases)[-remainder:]
        cases[largest_fractions] += 1

    return cases


def weekly_cases_per_strain(
    weekly_cases: ArrayLike,
    strain_fraction: ArrayLike,
) -> NDArray[np.floating]:
    """
    Calculate weekly influenza case counts for multiple viral strains.

    Each row of `strain_fraction` represents one week, and each column
    represents a viral strain. The function multiplies the total case count
    for each week by the corresponding fraction for every strain.

    Parameters
    ----------
    weekly_cases : ArrayLike
        One-dimensional array of estimated total influenza case counts per
        week. Values must be finite and non-negative.
    strain_fraction : ArrayLike
        Two-dimensional array of strain fractions. Rows correspond to weeks
        and columns correspond to viral strains. Values must be finite and
        lie between 0 and 1 inclusive. The number of rows must equal the
        number of elements in `weekly_cases`.

    Returns
    -------
    numpy.ndarray
        Two-dimensional floating-point array of estimated case counts per
        week and strain. Rows correspond to weeks and columns correspond to
        strains.

    Raises
    ------
    ValueError
        If `weekly_cases` is not one-dimensional, `strain_fraction` is not
        two-dimensional, the number of weeks does not match, either input
        contains non-finite values, `weekly_cases` contains negative values,
        or `strain_fraction` contains values outside the range from 0 to 1.

    Examples
    --------
    >>> weekly_cases_per_strain(
    ...     [10, 20, 30],
    ...     [[0.4, 0.6], [0.5, 0.5], [0.6, 0.4]],
    ... )
    array([[ 4.,  6.],
           [10., 10.],
           [18., 12.]])
    """
    weekly_cases = np.asarray(weekly_cases, dtype=float)
    strain_fraction = np.asarray(strain_fraction, dtype=float)

    if weekly_cases.ndim != 1:
        raise ValueError("weekly_cases must be one-dimensional.")
    if strain_fraction.ndim != 2:
        raise ValueError("strain_fraction must be two-dimensional.")
    if weekly_cases.shape[0] != strain_fraction.shape[0]:
        raise ValueError(
            "weekly_cases must have the same number of elements as "
            "strain_fraction has rows."
        )
    if not np.all(np.isfinite(weekly_cases)):
        raise ValueError("weekly_cases must contain only finite values.")
    if not np.all(np.isfinite(strain_fraction)):
        raise ValueError("strain_fraction must contain only finite values.")
    if np.any(weekly_cases < 0):
        raise ValueError("weekly_cases must not contain negative values.")
    if np.any((strain_fraction < 0) | (strain_fraction > 1)):
        raise ValueError("strain_fraction values must be between 0 and 1.")

    return weekly_cases[:, np.newaxis] * strain_fraction