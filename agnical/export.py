"""Write the small JSON / GeoJSON files the AgniCal website reads. Every number on the site
comes from these files, and every number in these files is computed here from the data."""
from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from . import config
from .grid import PERIOD_DAYS
from .harmonize import season_periods
from .regions import assign_district

KHARIF = list(config.SEASONS)[0]
RABI = list(config.SEASONS)[1]

# Independent published estimates of the peak burning hour, used only as an outside check.
SOURCES = [
    {"key": "jethva2026", "text": "Jethva, H. (2026). Timing the Flames: Geostationary Satellite Detection of Diurnally "
     "Shifting Stubble Burning in Northwestern India. Remote Sensing 18(10), 1506 (GEO-KOMPSAT-2A AMI, 2019-2025).",
     "url": "https://doi.org/10.3390/rs18101506"},
    {"key": "singh2025", "text": "Singh, N., Pradhan, R., Shukla, B. P. and Pandya, M. R. (Space Applications Centre, ISRO), "
     "communication in Current Science, 2025 (Meteosat SEVIRI), as reported by The Tribune, 10 Dec 2025.",
     "url": "https://www.tribuneindia.com/news/haryana/farmers-shift-stubble-burning-time-to-dodge-satellite-tracking-isro-study"},
    {"key": "iforest2025", "text": "iFOREST (2025). Stubble Burning Status Report 2025 (Sentinel-2 burnt area, Meteosat SEVIRI "
     "timing), as reported by The Tribune, 9 Dec 2025.",
     "url": "https://www.tribuneindia.com/news/delhi/satellites-miss-majority-of-stubble-fires-delhi-air-pollution-underestimated-report"},
]
# Burnt area mapped from Sentinel-2 by iFOREST: an outside measure that does not depend on the hour.
BURNT_AREA = {"source": "iforest2025", "state": "Punjab", "from_year": 2022, "from_km2": 31447,
              "to_year": 2025, "to_km2": 20000}

PUBLISHED = [
    {"year": 2020, "peak_ist": 13.5, "label": "Singh et al., ISRO SAC (Meteosat SEVIRI): peak ~13:30 in 2020",
     "url": SOURCES[1]["url"]},
    {"year": 2024, "peak_ist": 17.0, "label": "Singh et al., ISRO SAC (Meteosat SEVIRI): peak ~17:00 in 2024",
     "url": SOURCES[1]["url"]},
    {"year": 2025, "peak_lo": 16.0, "peak_hi": 18.0,
     "label": "Jethva 2026 (GEO-KOMPSAT-2A): fire signal spikes around and after 16:00",
     "url": SOURCES[0]["url"]},
]


def _r(x, n=1):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)


def _season_mask(doy: pd.Series, season: str) -> pd.Series:
    lo, hi = config.SEASONS[season]
    return doy.between(lo, hi + PERIOD_DAYS - 1)


def _cells_to_districts(h: pd.DataFrame, districts) -> pd.DataFrame:
    cells = h[["ix", "iy", "lat", "lon"]].drop_duplicates()
    cells = assign_district(cells, districts)
    return h.merge(cells[["ix", "iy", "district", "state"]], on=["ix", "iy"], how="left")


def season_totals(h: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    rows = []
    for s in config.SEASONS:
        sub = h[h["period"].isin(season_periods(s))]
        g = sub.groupby(by + ["year"]).agg(
            harmonized=("viirs_equiv", "sum"), snpp=("v_count", "sum"), modis=("m_count", "sum"),
            stitched=("stitched_raw", "sum"), source=("source", "first")).reset_index()
        g["season"] = s
        rows.append(g)
    return pd.concat(rows, ignore_index=True)


def write_site_data(out: Path, det, harmonized, timing, indicators, districts, meta, cv, ctrl=None) -> None:
    out.mkdir(parents=True, exist_ok=True)
    h = _cells_to_districts(harmonized, districts)
    h = h[h["state"].notna()]

    # ---------- burning time ----------
    tim = timing.merge(indicators, on="year", how="left")
    timing_json = {
        "years": [{k: _r(v, 3) if isinstance(v, float) else v for k, v in r.items()} for r in tim.to_dict("records")],
        "baseline_peak_ist": _r(timing.loc[timing["year"].between(2012, 2022), "peak_ist"].mean(), 3),
        "published": PUBLISHED, "sources": SOURCES,
        "model": meta["diurnal"],
    }
    (out / "timing.json").write_text(json.dumps(timing_json, default=float))

    # ---------- season totals: region and states ----------
    st = season_totals(h, ["state"])
    reg = season_totals(h.assign(state="Punjab + Haryana"), ["state"])
    tot = pd.concat([reg, st], ignore_index=True)
    k = metrics_band(meta)
    tot = tot.merge(timing[["year", "timing_break", "timing_confirmed"]], on="year", how="left")
    record = {}
    for (area, season), g in tot.groupby(["state", "season"]):
        q10, q90 = k.get(season, (np.nan, np.nan))
        rows = []
        for r in g.sort_values("year").itertuples():
            rec = r.source == "reconstructed"
            row = {"year": int(r.year), "source": r.source, "harmonized": _r(r.harmonized, 0),
                   "snpp": _r(r.snpp, 0) if r.snpp > 0 else None, "modis": _r(r.modis, 0),
                   "stitched": _r(r.stitched, 0)}
            if rec:
                row["lo"] = _r(r.harmonized / (1 + q90), 0)
                row["hi"] = _r(r.harmonized / (1 + q10), 0)
            if season == KHARIF:
                row["timing_break"] = bool(r.timing_break) if pd.notna(r.timing_break) else False
                row["timing_confirmed"] = bool(r.timing_confirmed) if pd.notna(r.timing_confirmed) else False
            rows.append(row)
        record.setdefault(area, {})[season] = rows
    (out / "record.json").write_text(json.dumps(record))

    # ---------- burning calendar: region, states, districts ----------
    cal = {}
    def add_cal(name, sub):
        piv = sub.groupby(["year", "period"])["viirs_equiv"].sum().unstack("period").reindex(
            columns=range(config.N_PERIODS)).fillna(0)
        src = sub.groupby("year")["source"].agg(lambda s: "observed" if (s == "observed").any() else "reconstructed")
        cal[name] = {"years": [int(y) for y in piv.index], "source": [src.get(y, "observed") for y in piv.index],
                     "z": [[int(round(v)) for v in row] for row in piv.to_numpy()]}
    add_cal("Punjab + Haryana", h)
    for s, g in h.groupby("state"):
        add_cal(s, g)
    for d, g in h.groupby("district"):
        add_cal(d, g)
    labels = [(pd.Timestamp("2021-01-01") + pd.Timedelta(days=p * PERIOD_DAYS)).strftime("%d %b")
              for p in range(config.N_PERIODS)]
    (out / "calendar.json").write_text(json.dumps({"period_labels": labels, "areas": cal}))

    # ---------- districts map ----------
    dk = season_totals(h, ["district", "state"])
    dk = dk[dk["season"] == KHARIF]
    piv = dk.pivot_table(index="district", columns="year", values="harmonized", aggfunc="sum").fillna(0)
    night = det.assign(year=det["date"].dt.year, doy=det["date"].dt.dayofyear)
    night = night[_season_mask(night["doy"], KHARIF) & (night["sensor"] == "VIIRS")]
    ns = night.groupby(["district", "year"]).apply(lambda g: (g["daynight"] != "D").mean()).unstack("year")
    gdf = districts.copy()
    ref_years = [y for y in range(2017, 2022) if y in piv.columns]
    last_full = max(y for y in piv.columns if y <= 2025)
    gdf["kharif_by_year"] = gdf["district"].map(lambda d: {int(y): int(round(v)) for y, v in piv.loc[d].items()}
                                                if d in piv.index else {})
    gdf["ref_mean_2017_21"] = gdf["district"].map(lambda d: _r(piv.loc[d, ref_years].mean(), 0) if d in piv.index else None)
    gdf[f"kharif_{last_full}"] = gdf["district"].map(lambda d: _r(piv.loc[d, last_full], 0) if d in piv.index else None)
    gdf["change_pct"] = gdf.apply(lambda r: _r((r[f"kharif_{last_full}"] / r["ref_mean_2017_21"] - 1) * 100, 0)
                                  if r["ref_mean_2017_21"] else None, axis=1)
    gdf["night_share_recent"] = gdf["district"].map(
        lambda d: _r(ns.loc[d, [y for y in (2024, 2025, 2026) if y in ns.columns]].mean(), 3) if d in ns.index else None)
    gdf["night_share_2017_21"] = gdf["district"].map(
        lambda d: _r(ns.loc[d, [y for y in ref_years if y in ns.columns]].mean(), 3) if d in ns.index else None)
    gdf["kharif_by_year"] = gdf["kharif_by_year"].map(json.dumps)
    gdf.to_file(out / "districts.geojson", driver="GeoJSON")

    # ---------- 2026 live season, one sensor (NOAA-20) for all years ----------
    vi = det[det["sensor"] == "VIIRS"]
    has26 = [s_ for s_ in ("NOAA20", "NOAA21") if ((vi["satellite"] == s_) & (vi["date"].dt.year == 2026)).any()]
    live_sat = has26[0] if has26 else "NOAA20"
    j1 = vi[vi["satellite"] == live_sat].copy()
    j1["year"] = j1["date"].dt.year
    j1["md"] = j1["date"].dt.strftime("%m-%d")
    j1 = j1[(j1["md"] >= "09-15") & (j1["md"] <= "12-15")]
    days = pd.date_range("2021-09-15", "2021-12-15").strftime("%m-%d")
    live = {"days": list(days), "series": {}, "satellite": {"NOAA20": "NOAA-20", "NOAA21": "NOAA-21"}[live_sat]}
    if ctrl is not None and len(ctrl["handover"]):
        ho = ctrl["handover"][ctrl["handover"]["satellite"] == live_sat]
        if len(ho):
            live["snpp_per_new"] = {"pooled": _r(ho["snpp"].sum() / ho["new"].sum(), 3),
                                    "lo": _r(ho["snpp_per_new"].min(), 3), "hi": _r(ho["snpp_per_new"].max(), 3),
                                    "years": [int(y) for y in ho["year"]]}
    for area, sub in [("Punjab + Haryana", j1)] + list(j1.groupby("state")):
        live["series"][area] = {}
        for y, g in sub.groupby("year"):
            dd = g[g["daynight"] == "D"].groupby("md").size().reindex(days, fill_value=0).cumsum()
            nn = g[g["daynight"] != "D"].groupby("md").size().reindex(days, fill_value=0).cumsum()
            last = g["md"].max() if y == 2026 else None
            live["series"][area][int(y)] = {"day": dd.tolist(), "night": nn.tolist(), "last": last}
    live["data_end"] = str(det["date"].max().date())
    (out / "live.json").write_text(json.dumps(live))

    # ---------- validation ----------
    m = meta["metrics"]
    reg = cv.groupby(["year", "period"])[["v_count"] + [c for c in cv.columns if c.startswith("pred_")]].sum().reset_index()
    val = {"metrics": m, "drift_test": meta["drift_test"], "feature_importance": meta["feature_importance"],
           "cv_regional": {int(y): {"period": g["period"].tolist(), "obs": g["v_count"].round(0).tolist(),
                                    "gbm": g["pred_gbm_poisson"].round(0).tolist(),
                                    "ratio": g["pred_global_ratio"].round(0).tolist()}
                           for y, g in reg.groupby("year")}}
    (out / "validation.json").write_text(json.dumps(val, default=float))

    # ---------- controls ----------
    if ctrl is not None:
        def recs(df, nd=4):
            return [{k: (_r(v, nd) if isinstance(v, (float, np.floating)) else (int(v) if isinstance(v, (np.integer,)) else
                         (bool(v) if isinstance(v, (np.bool_,)) else v))) for k, v in r.items()}
                    for r in df.to_dict("records")]
        rg = ctrl["regions"]
        keep = ["region", "year", "snpp_day", "snpp_night", "aqua_day", "terra_day", "aqua_over_snpp", "night_share",
                "terra_over_aqua", "aqua_over_snpp_x_baseline", "night_share_x_baseline",
                "terra_over_aqua_x_baseline", "n_flags", "timing_break"]
        (out / "controls.json").write_text(json.dumps({
            "study": __import__("agnical.controls", fromlist=["x"]).STUDY,
            "control": __import__("agnical.controls", fromlist=["x"]).PAKISTAN,
            "regions": recs(rg[rg["year"] >= 2012][keep]), "did": recs(ctrl["did"], 3),
            "nadir": recs(ctrl["nadir"]), "late_view": recs(ctrl["late_view"], 3),
            "viirs_window": recs(ctrl["viirs_window"], 3), "handover": recs(ctrl["handover"], 3),
            "burnt_area": BURNT_AREA}, default=float))

    # ---------- the record as one downloadable table ----------
    rows = []
    lv = ctrl["late_view"].set_index("year") if ctrl is not None else None
    for area, seasons in record.items():
        for season, rr in seasons.items():
            for r in rr:
                if r["year"] > 2025:
                    continue
                row = {"area": area, "season": season, "year": r["year"], "source": r["source"],
                       "viirs_equivalent_detections": r["harmonized"], "lo_p10": r.get("lo"), "hi_p90": r.get("hi"),
                       "raw_modis": r["modis"], "raw_suomi_npp": r["snpp"],
                       "timing_break": r.get("timing_break", ""), "timing_break_confirmed": r.get("timing_confirmed", "")}
                if lv is not None and area == "Punjab + Haryana" and season == KHARIF and r["year"] in lv.index and r["year"] >= 2023:
                    row["late_view_floor_factor"] = _r(lv.loc[r["year"], "factor"], 2)
                    row["late_view_floor_factor_lo"] = _r(lv.loc[r["year"], "factor_lo"], 2)
                    row["late_view_floor_factor_hi"] = _r(lv.loc[r["year"], "factor_hi"], 2)
                rows.append(row)
    pd.DataFrame(rows).to_csv(out / "agnical_record.csv", index=False)

    # ---------- headline numbers ----------
    kr = pd.DataFrame(record["Punjab + Haryana"][KHARIF]).set_index("year")
    ind = indicators.set_index("year")
    peak_yr = int(kr.loc[2012:2025, "snpp"].idxmax())
    summary = {
        "data_end": live["data_end"],
        "detections_in_states": meta["data"]["rows_in_states"],
        "by_satellite": meta["data"]["by_satellite"],
        "snpp_kharif": {int(y): _r(v, 0) for y, v in kr["snpp"].dropna().items()},
        "snpp_peak_year": peak_yr,
        "snpp_drop_peak_to_2025_pct": _r((kr.loc[2025, "snpp"] / kr.loc[peak_yr, "snpp"] - 1) * 100, 0),
        "night_snpp": {int(y): int(v) for y, v in ind["snpp_night"].items() if y >= 2012},
        "aqua_over_snpp": {int(y): _r(v, 3) for y, v in ind["aqua_over_snpp"].dropna().items()},
        "aqua_time": {int(y): _r(v, 2) for y, v in ind["aqua_time_ist"].dropna().items()},
        "terra_over_aqua": {int(y): _r(v, 3) for y, v in ind["terra_over_aqua"].dropna().items()},
        "baseline_peak_ist": timing_json["baseline_peak_ist"],
        "peak": {int(r.year): [_r(r.peak_ist, 2), _r(r.peak_lo, 2), _r(r.peak_hi, 2), bool(r.beyond_view)]
                 for r in timing.itertuples()},
        "timing_break_years": [int(r.year) for r in timing.itertuples() if r.timing_break == True],
        "timing_confirmed_years": [int(r.year) for r in timing.itertuples() if r.timing_confirmed == True],
        "night_snpp_ref_mean_2017_21": _r(np.mean([ind.loc[y, "snpp_night"] for y in range(2017, 2022)]), 0),
        "day_snpp_ref_mean_2017_21": _r(np.mean([ind.loc[y, "snpp_day"] for y in range(2017, 2022)]), 0),
    }
    (out / "summary.json").write_text(json.dumps(summary, default=float))


def metrics_band(meta) -> dict:
    g = meta["metrics"]["models"]["gbm_poisson"]["seasons"]
    return {s: (v["total_rel_err_q10"], v["total_rel_err_q90"]) for s, v in g.items()}
