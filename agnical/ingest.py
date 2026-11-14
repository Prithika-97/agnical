"""Read NASA FIRMS archive / NRT CSV files for MODIS and VIIRS into one clean table."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import config

_VIIRS_SAT = {"N": "SNPP", "SNPP": "SNPP", "SUOMI NPP": "SNPP", "S-NPP": "SNPP",
              "1": "NOAA20", "J1": "NOAA20", "N20": "NOAA20", "NOAA-20": "NOAA20",
              "2": "NOAA21", "N21": "NOAA21", "NOAA-21": "NOAA21"}
_MODIS_SAT = {"T": "Terra", "TERRA": "Terra", "A": "Aqua", "AQUA": "Aqua"}


def detect_sensor(columns, instrument=None) -> str:
    """Identify the sensor. Archive downloads use 'brightness' for both sensors, so the
    'instrument' column decides when it is present."""
    cols = {c.strip().lower() for c in columns}
    if instrument is not None and len(instrument):
        inst = str(pd.Series(instrument).dropna().astype(str).str.upper().mode().iat[0])
        return "MODIS" if "MODIS" in inst else "VIIRS"
    if "bright_ti4" in cols:
        return "VIIRS"
    if "brightness" in cols:
        return "MODIS"
    raise ValueError("Not a FIRMS MODIS or VIIRS active-fire file (no brightness column).")


def normalize(df: pd.DataFrame, source_file: str = "") -> pd.DataFrame:
    """Map one raw FIRMS file onto a common schema."""
    df = df.rename(columns=lambda c: c.strip().lower())
    sensor = detect_sensor(df.columns, df["instrument"] if "instrument" in df.columns else None)
    sat = df["satellite"].astype(str).str.strip().str.upper()
    out = pd.DataFrame({
        "lat": pd.to_numeric(df["latitude"], errors="coerce"),
        "lon": pd.to_numeric(df["longitude"], errors="coerce"),
        "date": pd.to_datetime(df["acq_date"], errors="coerce"),
        "time": pd.to_numeric(df["acq_time"], errors="coerce").fillna(0).astype(int),
        "satellite": sat.map(_MODIS_SAT if sensor == "MODIS" else _VIIRS_SAT).fillna(sat),
        "daynight": df["daynight"].astype(str).str.strip().str.upper(),
        "frp": pd.to_numeric(df["frp"], errors="coerce"),
        "confidence": df["confidence"].astype(str).str.strip().str.lower(),
        "type": pd.to_numeric(df["type"], errors="coerce") if "type" in df.columns else np.nan,
    })
    out["sensor"] = sensor
    out["product"] = "nrt" if "nrt" in Path(source_file).name.lower() else "sp"
    return out.dropna(subset=["lat", "lon", "date"])


def apply_filters(det: pd.DataFrame, bbox=config.BBOX,
                  modis_conf_min=config.MODIS_CONF_MIN,
                  viirs_drop_low=config.VIIRS_DROP_LOW,
                  viirs_sats=config.VIIRS_SATELLITES,
                  keep_types=config.KEEP_TYPES) -> pd.DataFrame:
    w, s, e, n = bbox
    keep = det["lon"].between(w, e) & det["lat"].between(s, n)
    # NRT files carry no 'type' field, so their static sources cannot be removed.
    keep &= det["type"].isna() | det["type"].isin(keep_types)
    is_m = det["sensor"] == "MODIS"
    conf_num = pd.to_numeric(det["confidence"], errors="coerce")
    keep &= ~is_m | (conf_num >= modis_conf_min)
    if viirs_drop_low:
        keep &= is_m | ~det["confidence"].isin(["l", "low"])
    keep &= is_m | det["satellite"].isin(viirs_sats)
    return det[keep]


def drop_nrt_overlap(det: pd.DataFrame) -> pd.DataFrame:
    """Where standard (sp) and near-real-time (nrt) data cover the same day for the same
    satellite, keep sp."""
    sp_days = det.loc[det["product"] == "sp", ["sensor", "satellite", "date"]].drop_duplicates()
    sp_days["_sp"] = True
    det = det.merge(sp_days, on=["sensor", "satellite", "date"], how="left")
    det = det[(det["product"] == "sp") | det["_sp"].isna()]
    return det.drop(columns="_sp")


def load_firms(raw_dir: Path = config.RAW_DIR, **filter_kw) -> tuple[pd.DataFrame, dict]:
    files = sorted(Path(raw_dir).rglob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No FIRMS CSV files found under {raw_dir}")
    frames = [normalize(pd.read_csv(f, low_memory=False), str(f)) for f in files]
    det = pd.concat(frames, ignore_index=True)
    n_raw = len(det)
    det = drop_nrt_overlap(det)
    det = apply_filters(det, **filter_kw)
    det = det.drop_duplicates(subset=["sensor", "satellite", "lat", "lon", "date", "time"])
    det = det.reset_index(drop=True)
    info = {
        "files": [f.name for f in files],
        "synthetic": any("synthetic" in f.name.lower() for f in files),
        "rows_raw": int(n_raw),
        "rows_kept": int(len(det)),
        "by_sensor": {k: int(v) for k, v in det["sensor"].value_counts().items()},
        "date_range": {s: [str(g["date"].min().date()), str(g["date"].max().date())]
                       for s, g in det.groupby("sensor")},
    }
    return det, info
