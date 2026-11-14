import pandas as pd

from agnical.regions import assign_district


def test_points_get_the_right_state_and_outside_points_none():
    df = pd.DataFrame({"lon": [75.84, 76.98, 74.34], "lat": [30.90, 29.69, 31.55]})   # Ludhiana, Karnal, Lahore
    out = assign_district(df)
    assert out.loc[0, "state"] == "Punjab" and out.loc[0, "district"] == "Ludhiana"
    assert out.loc[1, "state"] == "Haryana" and out.loc[1, "district"] == "Karnal"
    assert pd.isna(out.loc[2, "state"])
