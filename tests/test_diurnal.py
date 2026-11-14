import numpy as np
import pandas as pd

from agnical import diurnal as D


def simulate(peak_by_year, width=0.8, days=60, seed=0):
    """Passes from three satellites at known times, fires drawn from a known daily curve."""
    rng = np.random.default_rng(seed)
    rows = []
    for y, mu in peak_by_year.items():
        for d in pd.date_range(f"{y}-10-01", periods=days):
            amount = rng.gamma(3, 200)                       # how much burned that day: cancels out
            for sat, cls, t, sens in (("Terra", "MODIS", 10.6, 1.0), ("Aqua", "MODIS", 13.8 + rng.normal(0, .1), 1.0),
                                      ("SNPP", "VIIRS_SNPP", 13.6 + rng.normal(0, .3), 4.0)):
                lam = amount * sens * np.exp(-(t - mu) ** 2 / (2 * width ** 2))
                rows.append((d, y, sat, cls, t, rng.poisson(lam)))
    p = pd.DataFrame(rows, columns=["date", "year", "sat", "cls", "t", "n"])
    return p[p["n"] > 0]


def test_fit_recovers_known_peaks_and_sensitivity():
    truth = {2015: 13.5, 2016: 13.9}
    f = D.fit(simulate(truth))
    for y, mu in truth.items():
        assert abs(f["mu"][y] - mu) < 0.15
    assert abs(f["beta"]["VIIRS_SNPP"] - np.log(4.0)) < 0.25


def test_summary_marks_unmeasurable_peak_as_beyond_view():
    p = simulate({2015: 13.6, 2016: 16.5})
    f = D.fit(p)
    draws = D.bootstrap(p, f, n_boot=15)
    s = D.summarize(f, draws, p).set_index("year")
    assert not s.loc[2015, "beyond_view"]
    assert s.loc[2016, "beyond_view"] and np.isnan(s.loc[2016, "peak_ist"])


def test_timing_breaks_need_two_indicators():
    ind = pd.DataFrame({"year": range(2012, 2026),
                        "aqua_over_snpp": [0.2] * 11 + [0.21, 0.4, 0.6],
                        "night_share": [0.003, 0.004] * 5 + [0.003] + [0.0035, 0.05, 0.07],
                        "terra_over_aqua": [0.1, 0.09] * 5 + [0.1] + [0.095, 0.01, 0.008]})
    b = D.timing_breaks(ind).set_index("year")
    assert not b.loc[2023, "timing_break"]
    assert b.loc[2024, "timing_break"] and b.loc[2025, "timing_break"]


def test_hhmm():
    assert D.hhmm(13.5) == "13:30" and D.hhmm(25.25) == "01:15"
