"""Real state and district boundaries for the study area (Punjab and Haryana, India).

Source: geoBoundaries gbOpen IND (https://www.geoboundaries.org).
  ADM1 states: DataMeet India community / Election Commission of India, CC BY 2.5 IN.
  ADM2 districts: Pathways Data Pvt. Ltd. / lgdirectory.gov.in, ODbL 1.0.
"""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from . import config

BOUNDARY_FILE = config.ROOT / "data" / "boundaries" / "punjab_haryana_districts.geojson"
STATES = {"Punjab": "Punjab", "Haryāna": "Haryana"}


def build_boundary_file(adm1_path: Path, adm2_path: Path, out: Path = BOUNDARY_FILE) -> gpd.GeoDataFrame:
    """One-off: cut the national files down to the districts of Punjab and Haryana."""
    a1 = gpd.read_file(adm1_path)
    a2 = gpd.read_file(adm2_path)
    st = a1[a1["shapeName"].isin(STATES)].assign(state=lambda d: d["shapeName"].map(STATES))
    pts = a2.set_geometry(a2.representative_point())
    j = gpd.sjoin(pts, st[["state", "geometry"]], predicate="within")
    d = a2.loc[j.index, ["shapeName", "geometry"]].rename(columns={"shapeName": "district"})
    d["district"] = d["district"].str.strip()
    d["state"] = j["state"].to_numpy()
    d["geometry"] = d.geometry.simplify(0.003, preserve_topology=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    d[["district", "state", "geometry"]].to_file(out, driver="GeoJSON")
    return d


def load_districts(path: Path = BOUNDARY_FILE) -> gpd.GeoDataFrame:
    return gpd.read_file(path)


def assign_district(df: pd.DataFrame, districts: gpd.GeoDataFrame | None = None,
                    lon="lon", lat="lat") -> pd.DataFrame:
    """Add 'district' and 'state' columns by point-in-polygon; rows outside both states get NaN."""
    districts = load_districts() if districts is None else districts
    pts = gpd.GeoDataFrame(df[[lon, lat]].copy(), geometry=gpd.points_from_xy(df[lon], df[lat]), crs="EPSG:4326")
    j = gpd.sjoin(pts, districts[["district", "state", "geometry"]].to_crs("EPSG:4326"),
                  how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    out = df.copy()
    out["district"] = j["district"].reindex(df.index).to_numpy()
    out["state"] = j["state"].reindex(df.index).to_numpy()
    return out


def keep_study_area(df: pd.DataFrame, districts=None) -> pd.DataFrame:
    out = assign_district(df, districts)
    return out[out["state"].notna()].reset_index(drop=True)
