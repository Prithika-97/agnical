"""Grid point detections: cell x year x 8-day period."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config

PERIOD_DAYS = 8


def period_of(dates: pd.Series) -> pd.Series:
    return ((dates.dt.dayofyear - 1) // PERIOD_DAYS).clip(upper=config.N_PERIODS - 1)


def period_start(year, period) -> pd.Series:
    year = pd.Series(year).astype(int)
    period = pd.Series(period).astype(int)
    return pd.to_datetime(year.astype(str) + "-01-01") + pd.to_timedelta(period * PERIOD_DAYS, unit="D")


def _grid_index(bbox, res):
    w, s, e, n = bbox
    ix0, iy0 = int(np.floor(w / res)), int(np.floor(s / res))
    nx = int(np.floor(e / res)) - ix0 + 1
    ny = int(np.floor(n / res)) - iy0 + 1
    return ix0, iy0, nx, ny


def build_panel(det: pd.DataFrame, bbox=config.BBOX, res=config.GRID_RES,
                modis_first_year=config.MODIS_FIRST_YEAR,
                viirs_start=config.VIIRS_START, end_date=None,
                viirs_end=config.VIIRS_END) -> pd.DataFrame:
    """One row per active cell x year x period: MODIS features, VIIRS target.
    Cells with no detection in the whole record are dropped; the rest get all periods, zeros included.
    """
    ix0, iy0, nx, ny = _grid_index(bbox, res)
    det = det.copy()
    det["year"] = det["date"].dt.year
    det = det[det["year"] >= modis_first_year]
    years = np.arange(modis_first_year, det["year"].max() + 1)
    nyr, npd = len(years), config.N_PERIODS

    yi = (det["year"] - years[0]).to_numpy()
    pi = period_of(det["date"]).to_numpy()
    xi = (np.floor(det["lon"].to_numpy() / res).astype(int) - ix0).clip(0, nx - 1)
    ji = (np.floor(det["lat"].to_numpy() / res).astype(int) - iy0).clip(0, ny - 1)

    def dense(mask, weights=None, how="sum"):
        arr = np.zeros((nyr, npd, nx, ny))
        idx = (yi[mask], pi[mask], xi[mask], ji[mask])
        w = np.ones(mask.sum()) if weights is None else weights[mask]
        (np.maximum if how == "max" else np.add).at(arr, idx, w)
        return arr

    sensor = det["sensor"].to_numpy()
    sat = det["satellite"].to_numpy()
    dn = det["daynight"].to_numpy()
    frp = det["frp"].fillna(0).to_numpy()
    is_m = sensor == "MODIS"
    is_v = (sensor == "VIIRS") & (sat == config.TARGET_SAT)
    is_j1 = (sensor == "VIIRS") & (sat == "NOAA20")
    is_j2 = (sensor == "VIIRS") & (sat == "NOAA21")

    f = {
        "m_terra_day": dense(is_m & (sat == "Terra") & (dn == "D")),
        "m_terra_night": dense(is_m & (sat == "Terra") & (dn != "D")),
        "m_aqua_day": dense(is_m & (sat == "Aqua") & (dn == "D")),
        "m_aqua_night": dense(is_m & (sat == "Aqua") & (dn != "D")),
        "m_frp_sum": dense(is_m, frp),
        "m_frp_max": dense(is_m, frp, how="max"),
        "v_count": dense(is_v),
        "v_frp_sum": dense(is_v, frp),
        "j1_count": dense(is_j1),
        "j2_count": dense(is_j2),
    }
    m = f["m_terra_day"] + f["m_terra_night"] + f["m_aqua_day"] + f["m_aqua_night"]
    f["m_count"] = m
    # Spatial context: MODIS detections in the 8 surrounding cells.
    pad = np.pad(m, ((0, 0), (0, 0), (1, 1), (1, 1)))
    box = sum(pad[:, :, 1 + dx:1 + dx + nx, 1 + dy:1 + dy + ny]
              for dx in (-1, 0, 1) for dy in (-1, 0, 1))
    f["m_nbr"] = box - m
    # Temporal context: MODIS detections in the previous and next period.
    f["m_prev"] = np.concatenate([np.zeros_like(m[:, :1]), m[:, :-1]], axis=1)
    f["m_next"] = np.concatenate([m[:, 1:], np.zeros_like(m[:, :1])], axis=1)

    active = (m.sum(axis=(0, 1)) + f["v_count"].sum(axis=(0, 1)) + f["j1_count"].sum(axis=(0, 1))
              + f["j2_count"].sum(axis=(0, 1))) > 0
    ax, ay = np.nonzero(active)
    n_act = len(ax)
    Y, P, C = np.meshgrid(np.arange(nyr), np.arange(npd), np.arange(n_act), indexing="ij")
    Y, P, C = Y.ravel(), P.ravel(), C.ravel()
    gx, gy = ax[C], ay[C]

    panel = pd.DataFrame({
        "year": years[Y].astype(np.int16),
        "period": P.astype(np.int8),
        "ix": (gx + ix0).astype(np.int32),
        "iy": (gy + iy0).astype(np.int32),
    })
    panel["lon"] = (panel["ix"] + 0.5) * res
    panel["lat"] = (panel["iy"] + 0.5) * res
    for k, arr in f.items():
        panel[k] = arr[Y, P, gx, gy].astype(np.float32)

    # Validity: a period counts only if the sensor covered all of it. Coverage ends at
    # `end_date`; by default, the last day of the month of the latest detection.
    if end_date is None:
        end_date = det["date"].max() + pd.offsets.MonthEnd(0)
    end_date = pd.Timestamp(end_date)
    keys = panel[["year", "period"]].drop_duplicates()
    ps = period_start(keys["year"], keys["period"]).to_numpy()
    year_end = pd.to_datetime(keys["year"].astype(str) + "-12-31").to_numpy()
    pe = np.minimum(ps + np.timedelta64(PERIOD_DAYS - 1, "D"), year_end)
    keys = keys.assign(
        m_valid=bool(is_m.any()) & (pe <= end_date),
        v_valid=bool(is_v.any()) & (ps >= np.datetime64(pd.Timestamp(viirs_start)))
        & (pe <= min(end_date, pd.Timestamp(viirs_end))))
    return panel.merge(keys, on=["year", "period"], how="left")


def stitched_raw(panel: pd.DataFrame) -> pd.Series:
    """Naive record: raw MODIS before VIIRS, raw VIIRS after."""
    return panel["v_count"].where(panel["v_valid"], panel["m_count"])
