"""AgniCal v2 pipeline: NASA FIRMS CSVs -> harmonized record + burning-time analysis -> site data.

Usage:
    python pipeline.py              # reads data/raw, writes data/processed and site/data
    python pipeline.py --fast       # smaller models and fewer bootstrap draws, for a quick look
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from agnical import config, diurnal, export
from agnical.grid import build_panel
from agnical.harmonize import build_harmonized, cross_validate, drift_test, evaluate, ratio_stability
from agnical.ingest import load_firms
from agnical.regions import keep_study_area, load_districts
from agnical.report import write_results


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", type=Path, default=config.RAW_DIR)
    ap.add_argument("--out", type=Path, default=config.OUT_DIR)
    ap.add_argument("--site", type=Path, default=config.ROOT / "docs" / "data")
    ap.add_argument("--end", default=config.DATA_END, help="last day covered by the FIRMS requests")
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--boot", type=int, default=300, help="bootstrap draws for the burning-time model")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    a.site.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print(f"[1/7] Reading FIRMS files from {a.raw} ...")
    det, info = load_firms(a.raw)
    print(f"      {info['rows_kept']:,} detections in the download box: {info['by_sensor']}")

    print("[2/7] Keeping detections inside the real borders of Punjab and Haryana ...")
    districts = load_districts()
    det = keep_study_area(det, districts)
    info["rows_in_states"] = int(len(det))
    info["by_satellite"] = {f"{k[0]} {k[1]}": int(v) for k, v in det.groupby(["sensor", "satellite"]).size().items()}
    print(f"      {len(det):,} detections inside the two states: {info['by_satellite']}")

    print("[3/7] When in the day do the fires burn? (all five satellites, day-matched) ...")
    passes = diurnal.build_passes(det)
    dfit = diurnal.fit(passes)
    n_boot = 40 if a.fast else a.boot
    draws = diurnal.bootstrap(passes, dfit, n_boot=n_boot, log=print)
    timing = diurnal.summarize(dfit, draws, passes)
    indicators = diurnal.indicators(det)
    breaks = diurnal.timing_breaks(indicators)
    timing = timing.merge(breaks, on="year", how="left")
    print(timing[["year", "peak_ist", "peak_lo", "peak_hi", "beyond_view", "n_flags", "timing_break"]]
          .round(2).tail(8).to_string(index=False))

    print("[4/7] Gridding to 0.1 degree cells x 8-day periods ...")
    panel = build_panel(det, end_date=a.end)
    print(f"      {panel[['ix', 'iy']].drop_duplicates().shape[0]:,} active cells, {len(panel):,} rows")

    print("[5/7] Validating the MODIS -> VIIRS mapping (leave one year out, 2012-2022) ...")
    cv = cross_validate(panel, fast=a.fast)
    metrics = evaluate(cv)
    drift = drift_test(panel, fast=a.fast)
    stab = ratio_stability(panel)

    print("[6/7] Building the harmonized 2003-2026 record ...")
    h, model = build_harmonized(panel, fast=a.fast)

    print("[7/7] Writing outputs ...")
    keep = ["year", "period", "ix", "iy", "lat", "lon", "m_count", "v_count", "j1_count", "j2_count",
            "model_pred", "viirs_equiv", "stitched_raw", "source"]
    hs = h[keep].copy()
    hs = hs[(hs["viirs_equiv"].fillna(0) >= 0.01) | (hs["m_count"] > 0) | (hs["v_count"] > 0)
            | (hs["j1_count"] > 0) | (hs["j2_count"] > 0)]
    for c in ["model_pred", "viirs_equiv"]:
        hs[c] = hs[c].round(3).astype("float32")
    hs.to_parquet(a.out / "harmonized.parquet", index=False)
    timing.to_csv(a.out / "burning_time.csv", index=False)
    indicators.to_csv(a.out / "timing_indicators.csv", index=False)
    meta = {"data": info, "metrics": metrics, "feature_importance": model.importance().round(4).to_dict(),
            "drift_test": drift.to_dict(orient="records"),
            "diurnal": {"sigma0_h": dfit["sigma0"], "beta": dfit["beta"], "n_days": dfit["n_days"],
                        "n_passes": dfit["n_passes"], "n_boot": n_boot,
                        "baseline_peak_ist": timing.attrs.get("baseline_peak_ist")},
            "runtime_s": round(time.time() - t0, 1)}
    (a.out / "meta.json").write_text(json.dumps(meta, indent=2, default=float))
    write_results(meta, stab, config.REPORT_DIR / "RESULTS.md")
    export.write_site_data(a.site, det=det, harmonized=hs, timing=timing, indicators=indicators,
                           districts=districts, meta=meta, cv=cv)
    print(f"Done in {round(time.time() - t0, 1)} s. Site data in {a.site}")


if __name__ == "__main__":
    main()
