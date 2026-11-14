import numpy as np
import pandas as pd
import pytest

from agnical import config
from agnical.burncal import anomaly, regional_series, season_stats
from agnical.grid import build_panel, period_of
from agnical.harmonize import CellRatio, GlobalRatio, build_harmonized, cross_validate, evaluate
from agnical.ingest import apply_filters, detect_sensor, drop_nrt_overlap, normalize

BBOX = (75.0, 30.0, 75.5, 30.5)


def raw_modis(n=3):
    return pd.DataFrame({
        "latitude": [30.1] * n, "longitude": [75.1] * n, "brightness": 320.0, "scan": 1, "track": 1,
        "acq_date": ["2015-10-20"] * n, "acq_time": [530] * n, "satellite": ["T", "A", "Aqua"][:n],
        "instrument": "MODIS", "confidence": [80, 10, 90][:n], "version": "6.1", "bright_t31": 300,
        "frp": 10.0, "daynight": "D", "type": [0, 0, 2][:n]})


def raw_viirs():
    return pd.DataFrame({
        "latitude": [30.1, 30.2, 35.0], "longitude": [75.1, 75.2, 75.2], "bright_ti4": 330.0,
        "scan": 0.4, "track": 0.4, "acq_date": "2015-10-20", "acq_time": 820,
        "satellite": ["N", "1", "N"], "instrument": "VIIRS", "confidence": ["n", "h", "h"],
        "version": "2", "bright_ti5": 290, "frp": 3.0, "daynight": "D", "type": 0})


def test_detect_sensor():
    assert detect_sensor(raw_modis().columns) == "MODIS"
    assert detect_sensor(raw_viirs().columns) == "VIIRS"
    with pytest.raises(ValueError):
        detect_sensor(["a", "b"])
    # Real archive downloads name the VIIRS column 'brightness' and say SNPP in 'instrument'.
    assert detect_sensor(["brightness", "instrument"], pd.Series(["SNPP", "SNPP"])) == "VIIRS"
    assert detect_sensor(["brightness", "instrument"], pd.Series(["MODIS"])) == "MODIS"


def test_filters_drop_low_conf_static_other_sats_and_outside_bbox():
    det = pd.concat([normalize(raw_modis()), normalize(raw_viirs())])
    kept = apply_filters(det, bbox=BBOX)
    m = kept[kept.sensor == "MODIS"]
    v = kept[kept.sensor == "VIIRS"]
    assert len(m) == 1 and m.iloc[0]["satellite"] == "Terra"  # low conf and type=2 removed
    assert set(v["satellite"]) == {"SNPP", "NOAA20"}            # out-of-box point removed


def test_nrt_dropped_when_standard_exists():
    a = normalize(raw_viirs(), "fire_archive_SV.csv")
    b = normalize(raw_viirs(), "fire_nrt_SV.csv")
    out = drop_nrt_overlap(pd.concat([a, b]))
    assert set(out["product"]) == {"sp"}


def test_period_of():
    d = pd.Series(pd.to_datetime(["2020-01-01", "2020-01-09", "2020-12-31"]))
    assert period_of(d).tolist() == [0, 1, 45]


def toy_detections(rng, years=range(2010, 2016)):
    rows = []
    for y in years:
        for day in pd.date_range(f"{y}-10-01", f"{y}-11-30"):
            k = rng.poisson(2)
            for _ in range(k):
                lat, lon = 30.05 + rng.random() * 0.4, 75.05 + rng.random() * 0.4
                rows.append(("MODIS", "Aqua", day, lat, lon))
                for _ in range(1 + rng.poisson(1.5)):
                    if day >= pd.Timestamp("2012-01-20"):
                        rows.append(("VIIRS", "SNPP", day, lat, lon))
    df = pd.DataFrame(rows, columns=["sensor", "satellite", "date", "lat", "lon"])
    df["daynight"], df["frp"], df["time"] = "D", 5.0, 800
    return df


def test_panel_and_validity():
    det = toy_detections(np.random.default_rng(0))
    p = build_panel(det, bbox=BBOX, modis_first_year=2010, end_date="2015-12-31", viirs_end="2015-12-31")
    assert p.groupby(["ix", "iy"]).size().nunique() == 1          # full panel per cell
    assert not p.loc[p.year == 2011, "v_valid"].any()
    assert p.loc[(p.year == 2013), "v_valid"].all()
    assert p["m_count"].sum() == (det.sensor == "MODIS").sum()


def test_cv_harmonize_and_calendar():
    det = toy_detections(np.random.default_rng(1))
    p = build_panel(det, bbox=BBOX, modis_first_year=2010, end_date="2015-12-31", viirs_end="2015-12-31")
    cv = cross_validate(p, fast=True, log=lambda *_: None)
    res = evaluate(cv)
    assert set(res["models"]) == {"global_ratio", "cell_ratio", "gbm_poisson"}
    h, _ = build_harmonized(p, fast=True)
    assert (h["viirs_equiv"] >= 0).all()
    assert set(h.loc[h.year < 2012, "source"]) == {"reconstructed"}
    s = regional_series(h, BBOX)
    st = season_stats(s, list(config.SEASONS)[0])
    assert len(st) == 6 and (st["start_doy"] <= st["end_doy"]).all()
    assert anomaly(st, 2015)["label"]


def test_ratio_baselines():
    df = pd.DataFrame({"ix": [1, 1, 2], "iy": [1, 1, 1], "m_count": [2, 2, 0], "v_count": [6, 6, 3]})
    assert GlobalRatio().fit(df).r == pytest.approx(15 / 4)
    pred = CellRatio(k=0).fit(df).predict(df)
    assert pred[0] == pytest.approx(6.0)
