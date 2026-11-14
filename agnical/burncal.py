"""Turn the harmonized grid into a burning calendar, season statistics and anomaly flags."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .grid import PERIOD_DAYS, period_start
from .harmonize import season_periods


def in_bbox(df: pd.DataFrame, bbox) -> pd.Series:
    w, s, e, n = bbox
    return df["lon"].between(w, e) & df["lat"].between(s, n)


def regional_series(h: pd.DataFrame, bbox) -> pd.DataFrame:
    """Sum over cells in the box: one row per year x period."""
    sub = h[in_bbox(h, bbox)]
    cols = ["viirs_equiv", "stitched_raw", "m_count", "v_count", "model_pred"]
    cols = [c for c in cols if c in sub.columns]
    g = sub.groupby(["year", "period"])[cols].sum(min_count=1).reset_index()
    src = sub.groupby(["year", "period"])["source"].first().reset_index(drop=True)
    g["source"] = src.to_numpy()
    g["date"] = period_start(g["year"], g["period"]).to_numpy()
    return g


def season_stats(series: pd.DataFrame, season: str) -> pd.DataFrame:
    """Per year: season total, peak period, and the 10%-90% burning window."""
    sp = season_periods(season)
    rows = []
    for y, g in series[series["period"].isin(sp)].groupby("year"):
        g = g.sort_values("period")
        v = g["viirs_equiv"].fillna(0).to_numpy()
        total = float(v.sum())
        if total <= 0:
            continue
        cum = np.cumsum(v) / total
        peak_p = int(g["period"].iloc[int(np.argmax(v))])
        start_p = int(g["period"].iloc[int(np.searchsorted(cum, 0.10))])
        end_p = int(g["period"].iloc[min(int(np.searchsorted(cum, 0.90)), len(g) - 1)])
        rows.append({
            "year": int(y), "total": total, "source": g["source"].iloc[0],
            "peak_period": peak_p, "peak_doy": peak_p * PERIOD_DAYS + 1 + PERIOD_DAYS // 2,
            "start_doy": start_p * PERIOD_DAYS + 1, "end_doy": end_p * PERIOD_DAYS + PERIOD_DAYS,
            "stitched_total": float(g["stitched_raw"].fillna(0).sum()),
        })
    return pd.DataFrame(rows)


def anomaly(stats: pd.DataFrame, year: int, z_flag: float = config.ANOMALY_Z) -> dict | None:
    """Compare one season with all other years of the harmonized record."""
    if year not in set(stats["year"]):
        return None
    cur = stats[stats["year"] == year].iloc[0]
    ref = stats[stats["year"] != year]
    if len(ref) < 3:
        return None
    sd = ref["total"].std(ddof=1)
    z = float((cur["total"] - ref["total"].mean()) / sd) if sd > 0 else 0.0
    shift = float(cur["peak_doy"] - ref["peak_doy"].median())
    if z >= z_flag:
        label = "Unusually high"
    elif z <= -z_flag:
        label = "Unusually low"
    else:
        label = "Normal"
    return {"year": int(year), "total": float(cur["total"]), "z": z, "label": label,
            "peak_shift_days": shift, "ref_mean": float(ref["total"].mean()),
            "ref_years": f"{int(ref['year'].min())}-{int(ref['year'].max())}"}


def doy_to_label(doy: float) -> str:
    return (pd.Timestamp("2021-01-01") + pd.Timedelta(days=int(doy) - 1)).strftime("%d %b")
