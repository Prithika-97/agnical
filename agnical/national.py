"""The same timing test for every state of India and every burning season.

For each state and season: the two model-free indicators (Aqua fires per Suomi NPP fire, share
of Suomi NPP fires seen on the night pass), the timing-break test against the state's own
2012-2022 range, a comparison with the rest of India (which shares the same satellites and
the same orbit drift), and the late-view floor.
"""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from . import config, controls, diurnal

STATE_FILE = config.ROOT / "data" / "boundaries" / "india_states.geojson"
REST = "Rest of India"
# Burning seasons by day of year. Post-monsoon: paddy residue (north-west), also Odisha/MP.
# Pre-monsoon: wheat residue and forest fires. Winter: sugarcane, shifting cultivation starts.
SEASONS = {
    "Post-monsoon (Oct-Nov)": (274, 334),
    "Pre-monsoon (Mar-May)": (60, 151),
    "Winter (Dec-Feb)": (335, 59),
}
MIN_FIRES = 300        # mean Suomi NPP daytime fires per season in 2012-2022 needed to test a state
MIN_NIGHT = 5          # and mean night-pass fires, or the night share is too noisy to use
ASCII = str.maketrans("āīūṇṭḍśṣ", "aiuntdss")


def build_state_file(adm1_path: Path, out: Path = STATE_FILE) -> gpd.GeoDataFrame:
    a1 = gpd.read_file(adm1_path)[["shapeName", "geometry"]].rename(columns={"shapeName": "state"})
    a1["state"] = a1["state"].str.translate(ASCII)
    a1["geometry"] = a1.geometry.simplify(0.02, preserve_topology=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    a1.to_file(out, driver="GeoJSON")
    return a1


def assign_state(det: pd.DataFrame, states: gpd.GeoDataFrame | None = None) -> pd.DataFrame:
    states = gpd.read_file(STATE_FILE) if states is None else states
    pts = gpd.GeoDataFrame(det[["lon", "lat"]].copy(), geometry=gpd.points_from_xy(det["lon"], det["lat"]),
                           crs="EPSG:4326")
    j = gpd.sjoin(pts, states[["state", "geometry"]], how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    return det.assign(state=j["state"].reindex(det.index).to_numpy())


def _in_season(doy: pd.Series, rng) -> pd.Series:
    lo, hi = rng
    return doy.between(lo, hi) if lo <= hi else ((doy >= lo) | (doy <= hi))


def _season_year(date: pd.Series, rng) -> pd.Series:
    """Seasons that cross New Year are labelled by the year they start in."""
    lo, hi = rng
    y = date.dt.year
    return y if lo <= hi else y.where(date.dt.dayofyear >= lo, y - 1)


def season_slice(det: pd.DataFrame, season: str) -> pd.DataFrame:
    rng = SEASONS[season]
    d = det[_in_season(det["date"].dt.dayofyear, rng)].copy()
    d["syear"] = _season_year(d["date"], rng)
    return d


def _counts(d: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    day = d["daynight"] == "D"
    snpp = (d["sensor"] == "VIIRS") & (d["satellite"] == "SNPP")
    aqua = (d["sensor"] == "MODIS") & (d["satellite"] == "Aqua")
    t = diurnal.to_ist_hours(d["time"])
    g = pd.DataFrame({**{k: d[k] for k in by},
                      "aqua_day": (aqua & day).astype(int), "snpp_day": (snpp & day).astype(int),
                      "snpp_night": (snpp & ~day).astype(int), "snpp_all": snpp.astype(int),
                      "aqua_t": t.where(aqua & day)}).groupby(by)
    out = g[["aqua_day", "snpp_day", "snpp_night", "snpp_all"]].sum()
    out["aqua_time_ist"] = g["aqua_t"].median()
    return out.reset_index()


def state_indicators(det: pd.DataFrame, baseline=controls.BASELINE, k=3.0) -> pd.DataFrame:
    """One row per state, season and year, with indicators, timing-break flags and the floor."""
    rows = []
    base = list(baseline)
    for season in SEASONS:
        d = season_slice(det.dropna(subset=["state"]), season)
        c = _counts(d, ["state", "syear"]).rename(columns={"syear": "year"})
        c = c[c["snpp_day"] > 0]
        c["aqua_over_snpp"] = c["aqua_day"] / c["snpp_day"]
        c["night_share"] = c["snpp_night"] / c["snpp_all"]
        for state, g in c.groupby("state"):
            b = g[g["year"].isin(base)]
            g = g.copy()
            g["season"] = season
            g["base_fires"] = b["snpp_day"].mean() if len(b) else np.nan
            g["base_night"] = b["snpp_night"].mean() if len(b) else np.nan
            testable = len(b) >= 8 and g["base_fires"].iloc[0] >= MIN_FIRES
            night_ok = testable and g["base_night"].iloc[0] >= MIN_NIGHT
            g["testable"] = testable
            flags = np.zeros(len(g), int)
            for col, ok in (("aqua_over_snpp", testable), ("night_share", night_ok)):
                if not ok:
                    g[f"{col}_x"] = np.nan
                    g[f"{col}_z"] = np.nan
                    continue
                pooled = (b["aqua_day"].sum() / b["snpp_day"].sum() if col == "aqua_over_snpp"
                          else b["snpp_night"].sum() / b["snpp_all"].sum())
                med = b[col].median()
                mad = 1.4826 * (b[col] - med).abs().median()
                scale = mad if mad > 0 else (b[col].std(ddof=1) or 1e-9)
                g[f"{col}_x"] = g[col] / pooled if pooled > 0 else np.nan
                g[f"{col}_z"] = (g[col] - med) / scale
                flags += (g[f"{col}_z"] > k).to_numpy().astype(int)
            g["n_flags"] = flags
            g["n_tests"] = int(testable) + int(night_ok)
            # late-view floor, as in controls.late_view_floor
            if testable and b["aqua_day"].sum() > 0:
                r_year = (b["snpp_day"] / b["aqua_day"].replace(0, np.nan)).dropna()
                r0 = b["snpp_day"].sum() / b["aqua_day"].sum()
                g["floor_factor"] = g["aqua_day"] * r0 / g["snpp_day"]
                g["floor_factor_lo"] = g["aqua_day"] * r_year.min() / g["snpp_day"]
                g["floor_factor_hi"] = g["aqua_day"] * r_year.max() / g["snpp_day"]
            rows.append(g)
    out = pd.concat(rows, ignore_index=True)
    out["status"] = np.select(
        [~out["testable"], (out["n_tests"] == 2) & (out["n_flags"] == 2), out["n_flags"] >= 1],
        ["too_few_fires", "timing_break", "one_indicator"], "no_break")
    return out


def versus_rest(det: pd.DataFrame, ind: pd.DataFrame, years=(2023, 2024, 2025), n_boot=300, seed=0) -> pd.DataFrame:
    """For every state and season with at least one flag in `years`: how much further did the state
    move than the rest of India, which is seen by the same drifting satellites? Day-level bootstrap."""
    hits = ind[ind["year"].isin(years) & (ind["n_flags"] >= 1)][["state", "season"]].drop_duplicates()
    out = []
    for season, hs in hits.groupby("season"):
        d = season_slice(det.dropna(subset=["state"]), season)
        d = d.assign(date=pd.to_datetime(d["syear"].astype(str) + "-01-01")
                     + pd.to_timedelta((d["date"] - pd.to_datetime(d["date"].dt.year.astype(str) + "-01-01")).dt.days
                                       + np.where(d["date"].dt.year > d["syear"], 366, 0), unit="D"))
        for state in hs["state"]:
            dd = d.assign(region=np.where(d["state"] == state, state, REST))
            r = _did(dd, state, years, n_boot, seed)
            out.append(r.assign(state=state, season=season))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def _did(d: pd.DataFrame, study: str, years, n_boot, seed) -> pd.DataFrame:
    c = _counts(d.assign(day=d["date"]), ["region", "syear", "day"]).rename(columns={"syear": "year"})
    wide = c.pivot_table(index=["year", "day"], columns="region",
                         values=["aqua_day", "snpp_day", "snpp_night", "snpp_all"], fill_value=0)
    specs = {"aqua_over_snpp": ("aqua_day", "snpp_day"), "night_share": ("snpp_night", "snpp_all")}
    base = list(controls.BASELINE)

    def calc(w):
        s = w.groupby(level="year").sum()
        res = {}
        for k, (num, den) in specs.items():
            v = []
            for reg in (study, REST):
                r = s[(num, reg)] / s[(den, reg)].replace(0, np.nan)
                bsel = s.index.isin(base)
                v.append(r / (s.loc[bsel, (num, reg)].sum() / max(s.loc[bsel, (den, reg)].sum(), 1e-9)))
            res[k] = (v[0], v[1])
        return res

    pt = calc(wide)
    rng = np.random.default_rng(seed)
    groups = [g for _, g in wide.groupby(level="year")]
    draws = {(k, y): [] for k in specs for y in years}
    for _ in range(n_boot):
        w = pd.concat([g.iloc[rng.integers(0, len(g), len(g))] for g in groups])
        cc = calc(w)
        for k in specs:
            for y in years:
                if y in cc[k][0].index:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        draws[(k, y)].append(cc[k][0][y] / cc[k][1][y])
    rows = []
    for k in specs:
        for y in years:
            if y not in pt[k][0].index:
                continue
            dr = np.array(draws[(k, y)], float)
            dr = dr[np.isfinite(dr)]
            rows.append({"indicator": k, "year": int(y), "state_x": float(pt[k][0][y]), "rest_x": float(pt[k][1][y]),
                         "ratio": float(pt[k][0][y] / pt[k][1][y]) if pt[k][1][y] else np.nan,
                         "ratio_lo": float(np.quantile(dr, 0.025)) if len(dr) > 20 else np.nan,
                         "ratio_hi": float(np.quantile(dr, 0.975)) if len(dr) > 20 else np.nan})
    return pd.DataFrame(rows)
