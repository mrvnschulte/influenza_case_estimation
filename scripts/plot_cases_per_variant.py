import influenza_case_estimation as ice

activity = ice.weekly_influenza_activity(
    [0.1, 0.2],
    [0.1, 0.25],
)

print(activity)