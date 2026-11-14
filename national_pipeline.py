"""All-India run: the timing test for every state and burning season.

    python national_pipeline.py --raw data/raw_india

Reads NASA FIRMS country downloads for India (MODIS C6.1 and VIIRS Suomi NPP at least),
writes data/processed/india_*.csv and docs/data/india.json + india_states.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from agnical import config, national
from agnical.ingest import load_firms

INDIA_BBOX = (68.0, 6.0, 98.0, 37.5)


def _r(v, nd=3):
    return None if v is None or (isinstance(v, float) and not np.isfinite(v)) else round(float(v), nd)


def export(ind: pd.DataFrame, vs: pd.DataFrame, info: dict, out: Path, years=(2023, 2024, 2025)) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    states = gpd.read_file(national.STATE_FILE)
    states["geometry"] = states.geometry.simplify(0.05, preserve_topology=True)
    (out / "india_states.json").write_text(states.to_json(drop_id=True, to_wgs84=True))
    rows = []
    for r in ind[ind["year"] >= 2012].itertuples():
        rows.append({"state": r.state, "season": r.season, "year": int(r.year), "snpp_day": int(r.snpp_day),
                     "snpp_night": int(r.snpp_night), "aqua_day": int(r.aqua_day),
                     "aqua_x": _r(r.aqua_over_snpp_x, 2), "night_x": _r(r.night_share_x, 2),
                     "n_flags": int(r.n_flags), "n_tests": int(r.n_tests), "status": r.status,
                     "base_fires": _r(r.base_fires, 0), "aqua_time": _r(r.aqua_time_ist, 2),
                     "floor": _r(getattr(r, "floor_factor", np.nan), 2),
                     "floor_lo": _r(getattr(r, "floor_factor_lo", np.nan), 2),
                     "floor_hi": _r(getattr(r, "floor_factor_hi", np.nan), 2)})
    versus = [{k: (_r(v, 2) if isinstance(v, float) else v) for k, v in r.items()} for r in vs.to_dict("records")]
    # a break is confirmed when the state also moved further than the rest of India on both indicators
    conf = {}
    if len(vs):
        ok = vs.assign(ok=vs["ratio_lo"] > 1).groupby(["state", "season", "year"])["ok"].agg(["all", "sum"])
        conf = {f"{s}|{se}|{y}": {"both": bool(a), "n": int(n)} for (s, se, y), (a, n) in ok.iterrows()}
    last = ind[ind["year"] == max(years)]
    summary = {
        "years": list(years), "seasons": list(national.SEASONS), "min_fires": national.MIN_FIRES,
        "data": info,
        "tested": {s: int(((last["season"] == s) & last["testable"]).sum()) for s in national.SEASONS},
        "breaks": {s: sorted(last[(last["season"] == s) & (last["status"] == "timing_break")]["state"]) for s in national.SEASONS},
    }
    if len(vs):
        aq = vs[vs["indicator"] == "aqua_over_snpp"]
        summary["drift_alone"] = {se: {int(y): _r(g["rest_x"].median(), 2) for y, g in d.groupby("year")}
                                  for se, d in aq.groupby("season")}
    (out / "india.json").write_text(json.dumps({"summary": summary, "rows": rows, "versus": versus, "confirmed": conf}))
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=config.ROOT / "data" / "raw_india")
    ap.add_argument("--site", type=Path, default=config.ROOT / "docs" / "data")
    ap.add_argument("--boot", type=int, default=300)
    a = ap.parse_args()
    t0 = time.time()
    print("[1/4] Reading FIRMS files for India ...")
    det, info = load_firms(a.raw, bbox=INDIA_BBOX)
    det["n"] = 1
    for f in sorted(a.raw.glob("*_counts.csv.gz")):       # pre-counted Suomi NPP (scripts/reduce_snpp.py)
        det = pd.concat([det, national.load_counted(f)], ignore_index=True)
    print(f"      {int(det['n'].sum()):,} detections ({time.time() - t0:.0f}s)")
    print("[2/4] Assigning states ...")
    det = national.assign_state(det).dropna(subset=["state"])
    print("[3/4] Indicators, timing breaks and floors per state and season ...")
    ind = national.state_indicators(det)
    vs = national.versus_rest(det, ind, n_boot=a.boot)
    config.OUT_DIR.mkdir(parents=True, exist_ok=True)
    ind.to_csv(config.OUT_DIR / "india_state_indicators.csv", index=False)
    vs.to_csv(config.OUT_DIR / "india_versus_rest.csv", index=False)
    print("[4/4] Writing site data ...")
    meta = {"detections": int(det["n"].sum()), "first": str(det["date"].min().date()), "last": str(det["date"].max().date()),
            "by_satellite": {k: int(v) for k, v in det.groupby("satellite")["n"].sum().items()}}
    s = export(ind, vs, meta, a.site)
    print(json.dumps(s["breaks"], indent=1))
    print(f"Done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
