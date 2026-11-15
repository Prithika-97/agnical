"""Synthetic FIRMS-like CSVs for pipeline tests. Named SYNTHETIC_*.csv; outputs are labelled as
test data. Not for results.

Usage:  python scripts/make_synthetic.py --out data/raw
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

rng = np.random.default_rng(7)
W, S, E, N = 73.8, 27.6, 77.6, 32.6
RES = 0.1


def season_curve(doy, year):
    shift = rng.normal(0, 4)
    kharif = 1.0 * np.exp(-0.5 * ((doy - (305 + shift)) / 12) ** 2)
    rabi = 0.45 * np.exp(-0.5 * ((doy - (120 + shift)) / 10) ** 2)
    return (kharif + rabi) * (0.8 + 0.4 * rng.random())


def main(out: Path, cells: int = 900):
    out.mkdir(parents=True, exist_ok=True)
    cx = rng.uniform(W, E, cells)
    cy = rng.uniform(S, N, cells)
    weight = rng.gamma(0.6, 1.0, cells)
    m_rows, v_rows = [], []
    for year in range(2003, 2026):
        days = pd.date_range(f"{year}-01-01", f"{year}-12-31")
        doy = days.dayofyear.to_numpy()
        lam = np.outer(season_curve(doy, year), weight) * 0.35
        fires = rng.poisson(lam)
        d_idx, c_idx = np.nonzero(fires)
        reps = fires[d_idx, c_idx]
        d_idx, c_idx = np.repeat(d_idx, reps), np.repeat(c_idx, reps)
        size = rng.lognormal(0, 0.8, len(d_idx))
        lat = cy[c_idx] + rng.normal(0, 0.03, len(d_idx))
        lon = cx[c_idx] + rng.normal(0, 0.03, len(d_idx))
        date = days[d_idx].strftime("%Y-%m-%d")
        for sat, hhmm, p_scale in (("Terra", 530, 2.6), ("Aqua", 820, 1.9)):
            det = rng.random(len(d_idx)) < 1 - np.exp(-size / p_scale)
            n = det.sum()
            m_rows.append(pd.DataFrame({
                "latitude": lat[det].round(4), "longitude": lon[det].round(4),
                "brightness": rng.normal(320, 8, n).round(1), "scan": 1.0, "track": 1.0,
                "acq_date": date[det], "acq_time": hhmm, "satellite": sat, "instrument": "MODIS",
                "confidence": rng.integers(30, 100, n), "version": "6.1NRT" if year == 2025 else "6.1",
                "bright_t31": 300.0, "frp": (size[det] * 12).round(1), "daynight": "D", "type": 0}))
        if year >= 2012:
            det = rng.random(len(d_idx)) < 1 - np.exp(-size / 0.6)
            if year == 2012:
                det &= pd.to_datetime(date) >= "2012-01-20"
            k = 1 + rng.poisson(size[det] * 0.8)
            idx = np.repeat(np.nonzero(det)[0], k)
            n = len(idx)
            v_rows.append(pd.DataFrame({
                "latitude": (lat[idx] + rng.normal(0, 0.004, n)).round(4),
                "longitude": (lon[idx] + rng.normal(0, 0.004, n)).round(4),
                "bright_ti4": rng.normal(335, 10, n).round(1), "scan": 0.4, "track": 0.4,
                "acq_date": date[idx], "acq_time": 820, "satellite": "N", "instrument": "VIIRS",
                "confidence": rng.choice(["l", "n", "h"], n, p=[0.1, 0.7, 0.2]), "version": "2.0",
                "bright_ti5": 295.0, "frp": (size[idx] * 4).round(1), "daynight": "D", "type": 0}))
    pd.concat(m_rows).to_csv(out / "SYNTHETIC_fire_archive_M-C61.csv", index=False)
    pd.concat(v_rows).to_csv(out / "SYNTHETIC_fire_archive_SV-C2.csv", index=False)
    print(f"Wrote synthetic MODIS and VIIRS files to {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("data/raw"))
    main(ap.parse_args().out)
