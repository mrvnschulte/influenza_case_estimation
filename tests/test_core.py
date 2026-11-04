import numpy as np
import pytest

from influenza_case_estimation.core import (
    weekly_cases_per_strain,
    weekly_influenza_activity,
    weekly_influenza_cases,
)


class TestWeeklyInfluenzaActivity:
    def test_multiplies_inputs_elementwise(self):
        result = weekly_influenza_activity([0.1, 0.2], [0.1, 0.25])
        np.testing.assert_allclose(result, [0.01, 0.05])

    def test_returns_ndarray(self):
        result = weekly_influenza_activity([0.1], [0.2])
        assert isinstance(result, np.ndarray)

    def test_zero_input_gives_zero_activity(self):
        result = weekly_influenza_activity([0.0, 0.2], [0.5, 0.0])
        np.testing.assert_array_equal(result, [0.0, 0.0])

    @pytest.mark.parametrize(
        "ili_rate, test_positivity",
        [
            ([[0.1, 0.2]], [0.1, 0.2]),
            ([0.1, 0.2], [[0.1, 0.2]]),
        ],
    )
    def test_raises_for_non_1d_input(self, ili_rate, test_positivity):
        with pytest.raises(ValueError, match="one-dimensional"):
            weekly_influenza_activity(ili_rate, test_positivity)

    def test_raises_for_different_input_lengths(self):
        with pytest.raises(ValueError, match="same length"):
            weekly_influenza_activity([0.1, 0.2], [0.1])

    @pytest.mark.parametrize("invalid_value", [-0.1, np.nan, np.inf])
    def test_raises_for_invalid_ili_rate(self, invalid_value):
        with pytest.raises(ValueError):
            weekly_influenza_activity([invalid_value], [0.1])

    @pytest.mark.parametrize("invalid_value", [-0.1, np.nan, np.inf])
    def test_raises_for_invalid_test_positivity(self, invalid_value):
        with pytest.raises(ValueError):
            weekly_influenza_activity([0.1], [invalid_value])


class TestWeeklyInfluenzaCases:
    def test_allocates_cases_proportionally(self):
        result = weekly_influenza_cases(100, [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(result, [17, 33, 50])

    def test_returns_integer_ndarray(self):
        result = weekly_influenza_cases(10, [1.0, 1.0])
        assert isinstance(result, np.ndarray)
        assert np.issubdtype(result.dtype, np.integer)

    def test_allocated_cases_sum_to_total(self):
        result = weekly_influenza_cases(101, [1.0, 2.0, 3.0])
        assert result.sum() == 101

    def test_zero_total_returns_zeros(self):
        result = weekly_influenza_cases(0, [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(result, [0, 0, 0])

    def test_zero_activity_weeks_receive_no_cases(self):
        result = weekly_influenza_cases(10, [0.0, 1.0, 0.0])
        np.testing.assert_array_equal(result, [0, 10, 0])

    def test_largest_remainder_allocation(self):
        # Expected counts [1.5, 1.5]; the tie is resolved deterministically.
        result = weekly_influenza_cases(3, [1.0, 1.0])
        np.testing.assert_array_equal(result, [1, 2])

    def test_raises_for_negative_total(self):
        with pytest.raises(ValueError, match="total_cases"):
            weekly_influenza_cases(-1, [1.0])

    def test_raises_for_non_1d_activity(self):
        with pytest.raises(ValueError, match="one-dimensional"):
            weekly_influenza_cases(10, [[1.0, 2.0]])

    @pytest.mark.parametrize("activity", [[-1.0, 2.0], [np.nan], [np.inf]])
    def test_raises_for_invalid_activity(self, activity):
        with pytest.raises(ValueError):
            weekly_influenza_cases(10, activity)

    def test_raises_for_zero_total_activity(self):
        with pytest.raises(ValueError, match="sum greater than zero"):
            weekly_influenza_cases(10, [0.0, 0.0])


class TestWeeklyCasesPerStrain:
    def test_multiplies_cases_by_strain_fractions(self):
        result = weekly_cases_per_strain(
            [10, 20, 30],
            [[0.4, 0.6], [0.5, 0.5], [0.6, 0.4]],
        )
        np.testing.assert_allclose(
            result,
            [[4.0, 6.0], [10.0, 10.0], [18.0, 12.0]],
        )

    def test_returns_ndarray(self):
        result = weekly_cases_per_strain([10], [[0.25, 0.75]])
        assert isinstance(result, np.ndarray)

    def test_output_shape(self):
        result = weekly_cases_per_strain(
            [10, 20],
            [[0.25, 0.25, 0.25, 0.25], [0.1, 0.2, 0.3, 0.4]],
        )
        assert result.shape == (2, 4)

    def test_strain_cases_sum_to_weekly_cases_when_fractions_sum_to_one(self):
        weekly_cases = np.array([100.0, 200.0])
        fractions = np.array([[0.2, 0.8], [0.3, 0.7]])

        result = weekly_cases_per_strain(weekly_cases, fractions)

        np.testing.assert_allclose(result.sum(axis=1), weekly_cases)

    def test_zero_cases_give_zero_strain_cases(self):
        result = weekly_cases_per_strain([0, 0], [[0.5, 0.5], [0.2, 0.8]])
        np.testing.assert_array_equal(result, np.zeros((2, 2)))

    def test_raises_for_non_1d_weekly_cases(self):
        with pytest.raises(ValueError, match="one-dimensional"):
            weekly_cases_per_strain([[10, 20]], [[0.5, 0.5]])

    def test_raises_for_non_2d_fractions(self):
        with pytest.raises(ValueError, match="two-dimensional"):
            weekly_cases_per_strain([10, 20], [0.5, 0.5])

    def test_raises_when_week_counts_do_not_match(self):
        with pytest.raises(ValueError, match="same number of elements"):
            weekly_cases_per_strain([10, 20], [[0.5, 0.5]])

    @pytest.mark.parametrize("cases", [[-1], [np.nan], [np.inf]])
    def test_raises_for_invalid_weekly_cases(self, cases):
        with pytest.raises(ValueError):
            weekly_cases_per_strain(cases, [[0.5, 0.5]])

    @pytest.mark.parametrize(
        "fractions",
        [
            [[-0.1, 1.1]],
            [[np.nan, 0.5]],
            [[np.inf, 0.5]],
            [[0.2, 1.1]],
        ],
    )
    def test_raises_for_invalid_strain_fractions(self, fractions):
        with pytest.raises(ValueError):
            weekly_cases_per_strain([10], fractions)
