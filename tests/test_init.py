import influenza_case_estimation as pkg


class TestPublicAPI:
    def test_all_symbols_are_importable(self):
        for name in pkg.__all__:
            assert hasattr(pkg, name), f"{name} fehlt im Paket"

    def test_core_functions_are_callable(self):
        assert callable(pkg.weekly_influenza_activity)
        assert callable(pkg.weekly_influenza_cases)
        assert callable(pkg.weekly_cases_per_strain)

    def test_data_loaders_are_callable(self):
        assert callable(pkg.load_weighted_ili_rate)
        assert callable(pkg.load_test_positivity)
        assert callable(pkg.load_symptomatic_illnesses)
        assert callable(pkg.load_strain_fractions)
