"""Run the full AgniCal pipeline: FIRMS CSVs -> grid -> validated harmonization -> app data.

Usage:
    python pipeline.py                      # reads data/raw, writes data/processed
    python pipeline.py --fast               # quicker model for a first look
    python pipeline.py --raw path/to/csvs --bbox 73.8 27.6 77.6 32.6
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from agnical import config
from agnical.grid import build_panel
from agnical.harmonize import build_harmonized, cross_validate, evaluate, ratio_stability
from agnical.ingest import load_firms
from agnical.report import write_results


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", type=Path, default=config.RAW_DIR)
    ap.add_argument("--out", type=Path, default=config.OUT_DIR)
    ap.add_argument("--bbox", type=float, nargs=4, default=config.BBOX, metavar=("W", "S", "E", "N"))
    ap.add_argument("--res", type=float, default=config.GRID_RES)
    ap.add_argument("--end", default=None, help="last day covered by your FIRMS request, e.g. 2025-12-31")
    ap.add_argument("--fast", action="store_true", help="smaller model, for a quick first run")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print(f"[1/5] Reading FIRMS files from {a.raw} ...")
    det, info = load_firms(a.raw, bbox=tuple(a.bbox))
    print(f"      kept {info['rows_kept']:,} of {info['rows_raw']:,} detections: {info['by_sensor']}")
    if info["synthetic"]:
        print("      NOTE: synthetic test data detected. Results are NOT real findings.")
    if not {"MODIS", "VIIRS"} <= set(info["by_sensor"]):
        raise SystemExit("Need both MODIS and VIIRS files in the raw folder.")

    print("[2/5] Gridding to cells x 8-day periods ...")
    panel = build_panel(det, bbox=tuple(a.bbox), res=a.res, end_date=a.end)
    print(f"      {panel[['ix', 'iy']].drop_duplicates().shape[0]:,} active cells, {len(panel):,} rows")

    print("[3/5] Leave-one-year-out validation on the overlap years ...")
    cv = cross_validate(panel, fast=a.fast)
    metrics = evaluate(cv)
    stab = ratio_stability(panel)

    print("[4/5] Fitting final model and building the harmonized record ...")
    h, model = build_harmonized(panel, fast=a.fast)

    print("[5/5] Writing outputs ...")
    keep = ["year", "period", "ix", "iy", "lat", "lon", "m_count", "v_count",
            "model_pred", "viirs_equiv", "stitched_raw", "source"]
    hs = h[keep].copy()
    hs = hs[(hs["viirs_equiv"].fillna(0) >= 0.01) | (hs["m_count"] > 0) | (hs["v_count"] > 0)]
    for c in ["model_pred", "viirs_equiv"]:
        hs[c] = hs[c].round(3).astype("float32")
    hs.to_parquet(a.out / "harmonized.parquet", index=False)
    cv_reg = cv.groupby(["year", "period"])[
        ["v_count", "m_count"] + [c for c in cv.columns if c.startswith("pred_")]].sum().reset_index()
    cv_reg.to_parquet(a.out / "cv_regional.parquet", index=False)
    stab.to_parquet(a.out / "ratio_stability.parquet", index=False)
    meta = {"data": info, "bbox": list(a.bbox), "grid_res_deg": a.res, "fast_mode": a.fast,
            "metrics": metrics, "feature_importance": model.importance().round(4).to_dict(),
            "runtime_s": round(time.time() - t0, 1)}
    (a.out / "meta.json").write_text(json.dumps(meta, indent=2))
    write_results(meta, stab, config.REPORT_DIR / "RESULTS.md")
    print(f"Done in {meta['runtime_s']} s. Results: {config.REPORT_DIR / 'RESULTS.md'}")
    print("Launch the app with:  streamlit run app.py")


if __name__ == "__main__":
    main()
