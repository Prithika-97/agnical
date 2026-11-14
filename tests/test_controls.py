import numpy as np
import pandas as pd

from agnical import controls as K
from agnical import diurnal as D


def _det(rng, shift_study=True):
    """Two regions seen by the same satellites. From 2024 the study region burns later:
    Aqua (drifted late) sees more, and more fires last into the night. The control does not change."""
    rows = []
    for region in (K.STUDY, K.PAKISTAN):
        for y in range(2012, 2026):
            late = shift_study and region == K.STUDY and y >= 2024
            for day in pd.date_range(f"{y}-10-10", f"{y}-11-10"):
                n = rng.poisson(60)
                rows += [(region, day, "VIIRS", "SNPP", "D", 800, 0.4)] * n
                rows += [(region, day, "MODIS", "Aqua", "D", 830 if y < 2023 else 940, 1.0)] * rng.poisson(n * (0.6 if late else 0.2))
                rows += [(region, day, "VIIRS", "SNPP", "N", 2000, 0.4)] * rng.poisson(n * (0.05 if late else 0.004))
    d = pd.DataFrame(rows, columns=["region", "date", "sensor", "satellite", "daynight", "time", "scan"])
    d["lat"], d["lon"], d["frp"] = 30.0, 75.0, 5.0
    return d


def test_diff_in_diff_separates_burning_change_from_satellite_change():
    det = _det(np.random.default_rng(0))
    did = K.diff_in_diff(det, years=(2022, 2025), n_boot=60).set_index(["indicator", "year"])
    for k in ("aqua_over_snpp", "night_share"):
        assert did.loc[(k, 2025), "ratio_lo"] > 1.5          # study moved, control did not
    assert abs(did.loc[("aqua_over_snpp", 2022), "ratio"] - 1) < 0.15   # before the change: alike


def test_no_false_alarm_when_nothing_changes():
    det = _det(np.random.default_rng(1), shift_study=False)
    did = K.diff_in_diff(det, years=(2025,), n_boot=60)
    a = did.set_index("indicator").loc["aqua_over_snpp"]
    assert abs(a["ratio"] - 1) < 0.15 and a["ratio_lo"] < a["ratio"] < a["ratio_hi"]
    reg = K.region_indicators(det)
    assert not reg["timing_break"].any()


def test_late_view_floor_matches_count_before_drift_and_exceeds_it_after():
    det = _det(np.random.default_rng(2))
    ind = D.indicators(det[det["region"] == K.STUDY])
    f = K.late_view_floor(ind).set_index("year")
    assert abs(f.loc[2018, "factor"] - 1) < 0.1
    assert f.loc[2025, "factor_lo"] > 2
    assert (f["aqua_as_snpp_lo"] <= f["aqua_as_snpp"]).all() and (f["aqua_as_snpp"] <= f["aqua_as_snpp_hi"]).all()


def test_terra_ratio_is_not_a_flag():
    assert "terra_over_aqua" not in [c for c, _ in D.FLAG_INDICATORS]
