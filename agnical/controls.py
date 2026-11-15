"""Separate a change in burning from a change in the satellites.

1. Comparison regions: Pakistani Punjab strip, and Indian states with no flag. Same sensors,
   same passes, so orbit drift or view geometry must appear there as well.
2. Difference in differences (study / control, each against its 2012-22 level); day-level bootstrap.
3. Near-nadir check: Aqua / SNPP ratio from small-pixel detections only.
4. Late-view floor: Aqua day count x 2012-22 SNPP/Aqua ratio. Lower bound on the SNPP count
   at Aqua's hour.
5. VIIRS hand-over: NOAA-20/21 to SNPP units (SNPP data delivery ended 2026-11-02).
"""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from . import config, diurnal

REGION_FILE = config.ROOT / "data" / "boundaries" / "comparison_regions.geojson"
STUDY = "Punjab + Haryana (India)"
PAKISTAN = "Punjab, Pakistan (border strip)"
BASELINE = range(2012, 2023)
SEASON_DOY = (274, 334)
NADIR_SCAN = {"MODIS": 1.3, "VIIRS": 0.42}   # pixel width (km) below which a detection is near nadir


def build_region_file(adm1_path: Path, out: Path = REGION_FILE, bbox=config.BBOX) -> gpd.GeoDataFrame:
    """Cut the download box into study and comparison regions. Run once."""
    a1 = gpd.read_file(adm1_path)
    bx = box(*bbox)
    name = a1["shapeName"]
    def part(names):
        return a1[name.isin(names)].geometry.union_all().intersection(bx)
    rows = [
        (STUDY, part(["Punjab", "Haryāna"])),
        (PAKISTAN, bx.difference(a1.geometry.union_all().buffer(0.02))),
        ("Rajasthan (northern districts, India)", part(["Rājasthān"])),
        ("Uttar Pradesh (western districts, India)", part(["Uttar Pradesh"])),
    ]
    g = gpd.GeoDataFrame({"region": [r[0] for r in rows]}, geometry=[r[1] for r in rows], crs="EPSG:4326")
    g["geometry"] = g.geometry.simplify(0.005, preserve_topology=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    g.to_file(out, driver="GeoJSON")
    return g


def assign_region(det: pd.DataFrame, regions: gpd.GeoDataFrame | None = None) -> pd.DataFrame:
    regions = gpd.read_file(REGION_FILE) if regions is None else regions
    pts = gpd.GeoDataFrame(det[["lon", "lat"]].copy(), geometry=gpd.points_from_xy(det["lon"], det["lat"]),
                           crs="EPSG:4326")
    j = gpd.sjoin(pts, regions[["region", "geometry"]], how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    return det.assign(region=j["region"].reindex(det.index).to_numpy())


def region_indicators(det: pd.DataFrame) -> pd.DataFrame:
    """Indicators and timing-break test, per region."""
    out = []
    for r, g in det.dropna(subset=["region"]).groupby("region"):
        ind = diurnal.indicators(g)
        ind = ind[ind["snpp_day"] > 0]
        if ind["year"].isin(list(BASELINE)).sum() < 5:
            continue
        b = diurnal.timing_breaks(ind)
        out.append(ind.merge(b, on="year").assign(region=r))
    return pd.concat(out, ignore_index=True)


def _daily(det: pd.DataFrame) -> pd.DataFrame:
    d = det[det["date"].dt.dayofyear.between(*SEASON_DOY)]
    day = d["daynight"] == "D"
    snpp = (d["sensor"] == "VIIRS") & (d["satellite"] == "SNPP")
    aqua = (d["sensor"] == "MODIS") & (d["satellite"] == "Aqua")
    t = pd.DataFrame({"region": d["region"], "date": d["date"],
                      "aqua_day": (aqua & day).astype(int), "snpp_day": (snpp & day).astype(int),
                      "snpp_night": (snpp & ~day).astype(int), "snpp_all": snpp.astype(int)})
    t = t.groupby(["region", "date"]).sum().reset_index()
    t["year"] = t["date"].dt.year
    return t


def _stat(t: pd.DataFrame, num: str, den: str) -> pd.Series:
    g = t.groupby("year")[[num, den]].sum()
    return g[num] / g[den].replace(0, np.nan)


def diff_in_diff(det: pd.DataFrame, study=STUDY, control=PAKISTAN, years=(2023, 2024, 2025),
                 n_boot=500, seed=0) -> pd.DataFrame:
    """Per indicator and year: (study / 2012-22 level) / (control / 2012-22 level). A value of 1 means
    both regions moved together, i.e. a satellite effect. Days resampled within year, same days in both."""
    daily = _daily(det[det["region"].isin([study, control])])
    wide = daily.pivot_table(index=["year", "date"], columns="region",
                             values=["aqua_day", "snpp_day", "snpp_night", "snpp_all"], fill_value=0)
    specs = {"aqua_over_snpp": ("aqua_day", "snpp_day"), "night_share": ("snpp_night", "snpp_all")}
    base_years = list(BASELINE)

    def calc(w):
        res = {}
        s = w.groupby(level="year").sum()
        for k, (num, den) in specs.items():
            for reg in (study, control):
                r = s[(num, reg)] / s[(den, reg)].replace(0, np.nan)
                b = s.loc[s.index.isin(base_years), (num, reg)].sum() / s.loc[s.index.isin(base_years), (den, reg)].sum()
                res[(k, reg)] = r / b
        return res

    point = calc(wide)
    rng = np.random.default_rng(seed)
    by_year = {y: g for y, g in wide.groupby(level="year")}
    draws = {(k, y): [] for k in specs for y in years}
    for _ in range(n_boot):
        w = pd.concat([g.iloc[rng.integers(0, len(g), len(g))] for g in by_year.values()])
        c = calc(w)
        for k in specs:
            for y in years:
                if y in c[(k, study)].index:
                    draws[(k, y)].append(c[(k, study)][y] / c[(k, control)][y])
    rows = []
    for k in specs:
        for y in years:
            if y not in point[(k, study)].index:
                continue
            d = np.array(draws[(k, y)], float)
            d = d[np.isfinite(d)]
            rows.append({"indicator": k, "year": int(y),
                         "study_x_baseline": float(point[(k, study)][y]),
                         "control_x_baseline": float(point[(k, control)][y]),
                         "ratio": float(point[(k, study)][y] / point[(k, control)][y]),
                         "ratio_lo": float(np.quantile(d, 0.025)) if len(d) else np.nan,
                         "ratio_hi": float(np.quantile(d, 0.975)) if len(d) else np.nan})
    return pd.DataFrame(rows)


def nadir_check(det: pd.DataFrame) -> pd.DataFrame:
    """Aqua / SNPP daytime ratio, near-nadir detections only."""
    d = det[det["date"].dt.dayofyear.between(*SEASON_DOY) & (det["daynight"] == "D")]
    d = d.assign(year=d["date"].dt.year)
    aq = d[(d["satellite"] == "Aqua") & (d["sensor"] == "MODIS")]
    sn = d[(d["satellite"] == "SNPP") & (d["sensor"] == "VIIRS")]
    out = pd.DataFrame({
        "all": aq.groupby("year").size() / sn.groupby("year").size(),
        "near_nadir": aq[aq["scan"] <= NADIR_SCAN["MODIS"]].groupby("year").size()
        / sn[sn["scan"] <= NADIR_SCAN["VIIRS"]].groupby("year").size()}).dropna().reset_index()
    for c in ("all", "near_nadir"):
        out[f"{c}_x_baseline"] = out[c] / out.loc[out["year"].isin(list(BASELINE)), c].mean()
    return out


def late_view_floor(ind: pd.DataFrame) -> pd.DataFrame:
    """Aqua day count in SNPP units, using the 2012-22 SNPP/Aqua ratio (both near 13:40 then).
    From 2023 Aqua looks later, so this is the SNPP count expected at Aqua's hour: a floor, since
    the peak is later still. Range: lowest and highest yearly ratio, 2012-22."""
    b = ind[ind["year"].isin(list(BASELINE)) & (ind["snpp_day"] > 0)]
    r_year = b["snpp_day"] / b["aqua_day"]
    r0 = b["snpp_day"].sum() / b["aqua_day"].sum()
    o = ind[ind["snpp_day"] > 0][["year", "snpp_day", "aqua_day", "aqua_time_ist", "snpp_time_ist"]].copy()
    o["aqua_as_snpp"] = o["aqua_day"] * r0
    o["aqua_as_snpp_lo"] = o["aqua_day"] * r_year.min()
    o["aqua_as_snpp_hi"] = o["aqua_day"] * r_year.max()
    for s in ("", "_lo", "_hi"):
        o[f"factor{s}"] = o[f"aqua_as_snpp{s}"] / o["snpp_day"]
    o["baseline_ratio"] = r0
    return o


def viirs_window_test(det: pd.DataFrame) -> pd.DataFrame:
    """VIIRS-only test. On days with VIIRS passes before and after 13:30 (swath overlap, ~100 min
    apart) compare the two counts; a steeper rise means burning builds up later."""
    p = diurnal.build_passes(det[det["sensor"] == "VIIRS"])
    p["late"] = p["t"] >= 13.5
    g = p.groupby(["year", "date", "late"]).agg(n=("n", "sum"), t=("t", "mean")).unstack().dropna()
    rows = []
    for y, q in g.groupby(level="year"):
        early, late = q[("n", False)].sum(), q[("n", True)].sum()
        if early < 50 or late < 50:
            continue
        dt = float((q[("t", True)] - q[("t", False)]).mean())
        rows.append({"year": int(y), "days": int(len(q)), "early": int(early), "late": int(late),
                     "gap_h": dt, "log_rise_per_h": float(np.log(late / early) / dt)})
    return pd.DataFrame(rows)


def handover(det: pd.DataFrame, new_sats=("NOAA20", "NOAA21")) -> pd.DataFrame:
    """SNPP day detections per NOAA-20/21 detection, paddy-season days with both. Multiply the newer
    count by this to continue the record in SNPP units."""
    d = det[(det["sensor"] == "VIIRS") & (det["daynight"] == "D")
            & det["date"].dt.dayofyear.between(*SEASON_DOY)]
    c = d.groupby(["date", "satellite"]).size().unstack(fill_value=0)
    rows = []
    for s in new_sats:
        if s not in c or "SNPP" not in c:
            continue
        both = c[(c[s] > 0) & (c["SNPP"] > 0)]
        for y, g in both.groupby(both.index.year):
            rows.append({"satellite": s, "year": int(y), "days": int(len(g)), "snpp": int(g["SNPP"].sum()),
                         "new": int(g[s].sum()), "snpp_per_new": float(g["SNPP"].sum() / g[s].sum())})
    return pd.DataFrame(rows)
