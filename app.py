"""AgniCal web app.  Run:  streamlit run app.py"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from agnical import config
from agnical.burncal import anomaly, doy_to_label, in_bbox, regional_series, season_stats
from agnical.grid import PERIOD_DAYS
from agnical.harmonize import season_periods
from agnical.report import LABELS

st.set_page_config(page_title="AgniCal: burning calendar", layout="wide")
OUT = config.OUT_DIR
FIRE = "#d9480f"
MUTED = "#868e96"


@st.cache_data
def load():
    meta = json.loads((OUT / "meta.json").read_text())
    h = pd.read_parquet(OUT / "harmonized.parquet")
    cv = pd.read_parquet(OUT / "cv_regional.parquet")
    stab = pd.read_parquet(OUT / "ratio_stability.parquet")
    return meta, h, cv, stab


if not (OUT / "meta.json").exists():
    st.title("AgniCal")
    st.warning("No processed data yet. Put FIRMS CSV files in data/raw and run `python pipeline.py`.")
    st.stop()

meta, h, cv, stab = load()
metrics = meta["metrics"]

# ---------------- sidebar ----------------
st.sidebar.header("Area and season")
area = st.sidebar.selectbox("Area", list(config.AREA_PRESETS) + ["Custom box"])
if area == "Custom box":
    w0, s0, e0, n0 = meta["bbox"]
    lon = st.sidebar.slider("Longitude", w0, e0, (w0, e0), 0.1)
    lat = st.sidebar.slider("Latitude", s0, n0, (s0, n0), 0.1)
    bbox = (lon[0], lat[0], lon[1], lat[1])
else:
    bbox = config.AREA_PRESETS[area]
season = st.sidebar.radio("Season", list(config.SEASONS))

series = regional_series(h, bbox)
stats = season_stats(series, season)
if stats.empty:
    st.error("No burning detected in this area and season. Try a larger area.")
    st.stop()
year = st.sidebar.selectbox("Season to check", sorted(stats["year"], reverse=True))
st.sidebar.caption("Areas are approximate bounding boxes, not district or state boundaries.")

# ---------------- header ----------------
st.title("AgniCal")
st.markdown("**A consistent 20+ year crop-burning calendar for Punjab and Haryana, "
            "from NASA MODIS and VIIRS active-fire data.**")
if meta["data"]["synthetic"]:
    st.error("SYNTHETIC TEST DATA: this app is running on generated test files. "
             "Nothing shown here is a real finding.")

an = anomaly(stats, year)
cur = stats[stats["year"] == year].iloc[0]
c1, c2, c3, c4 = st.columns(4)
c1.metric("Peak week", doy_to_label(cur["peak_doy"]),
          f"{an['peak_shift_days']:+.0f} days vs typical" if an else None, delta_color="off",
          help="Middle of the 8-day period with the most detections.")
c2.metric("Burning window", f"{doy_to_label(cur['start_doy'])}–{doy_to_label(cur['end_doy'])}",
          help="Dates between which the middle 80% of the season's detections fall.")
c3.metric("Season total", f"{cur['total']:,.0f}",
          f"{(cur['total'] / an['ref_mean'] - 1) * 100:+.0f}% vs other years" if an else None,
          delta_color="off", help="VIIRS 375 m-equivalent fire detections in the selected area.")
if an:
    c4.metric("This season is", an["label"], f"z = {an['z']:+.1f} vs {an['ref_years']}",
              delta_color="off", help=f"Flagged when |z| exceeds {config.ANOMALY_Z}.")
if cur["source"] == "reconstructed":
    st.info(f"{year} is before VIIRS existed. Its values are reconstructed from MODIS by the "
            "harmonization model; see the uncertainty band in the trend tab.")

tab1, tab2, tab3, tab4 = st.tabs(["Burning calendar", "Long-term trend", "Map", "Method and validation"])

# ---------------- calendar ----------------
with tab1:
    piv = series.pivot(index="year", columns="period", values="viirs_equiv").reindex(
        columns=range(config.N_PERIODS)).fillna(0)
    xlab = [doy_to_label(p * PERIOD_DAYS + 1) for p in piv.columns]
    fig = go.Figure(go.Heatmap(
        z=np.log10(1 + piv.to_numpy()), x=xlab, y=piv.index.astype(str),
        customdata=piv.to_numpy(),
        colorscale=[[0, "#ffffff"], [0.08, "#fff3bf"], [0.35, "#ffa94d"], [0.65, "#e8590c"], [1, "#7f1d1d"]],
        zmin=0, zmax=max(1.0, float(np.log10(1 + piv.to_numpy().max()))),
        hovertemplate="%{y}, week of %{x}<br>%{customdata:,.0f} detections<extra></extra>",
        colorbar=dict(title="Detections<br>(log scale)", tickvals=[0, 1, 2, 3, 4],
                      ticktext=["0", "10", "100", "1k", "10k"])))
    obs_years = series.loc[series["source"] == "observed", "year"]
    if len(obs_years):
        first_obs = int(obs_years.min())
        fig.add_hline(y=first_obs - 0.5, line_dash="dash", line_color="#212529")
        fig.add_annotation(x=xlab[1], y=first_obs - 0.5, yshift=10, showarrow=False, xanchor="left",
                           text="above: reconstructed from MODIS  |  below: observed by VIIRS",
                           bgcolor="white", font=dict(size=11))
    fig.update_layout(height=620, margin=dict(l=10, r=10, t=30, b=10), xaxis=dict(nticks=12),
                      yaxis=dict(autorange="reversed", title=None),
                      title="Every 8-day period, every year, on one consistent scale")
    st.plotly_chart(fig, width="stretch")
    st.caption("Each cell is one 8-day period. Values are VIIRS 375 m-equivalent fire detections "
               "summed over the selected area.")

# ---------------- trend ----------------
with tab2:
    gbm = metrics["models"]["gbm_poisson"]["seasons"][season]
    q10, q90 = gbm["total_rel_err_q10"], gbm["total_rel_err_q90"]
    s = stats.copy()
    rec = s["source"] == "reconstructed"
    s["lo"] = np.where(rec, s["total"] / (1 + max(q90, -0.99)), np.nan)
    s["hi"] = np.where(rec, s["total"] / (1 + max(q10, -0.99)), np.nan)
    fig = go.Figure()
    fig.add_bar(x=s.loc[rec, "year"], y=s.loc[rec, "total"], name="Reconstructed from MODIS",
                marker_color="#ffa94d",
                error_y=dict(type="data", symmetric=False, array=(s["hi"] - s["total"])[rec],
                             arrayminus=(s["total"] - s["lo"])[rec], color=MUTED))
    fig.add_bar(x=s.loc[~rec, "year"], y=s.loc[~rec, "total"], name="Observed by VIIRS", marker_color=FIRE)
    fig.add_scatter(x=s["year"], y=s["stitched_total"], mode="lines+markers", name="Naive stitched record",
                    line=dict(color="#495057", dash="dot"))
    fig.update_layout(height=460, barmode="overlay", legend=dict(orientation="h", y=1.1),
                      yaxis_title="Season total (VIIRS-equivalent detections)", xaxis_title=None,
                      margin=dict(l=10, r=10, t=40, b=10))
    st.plotly_chart(fig, width="stretch")
    st.markdown(
        "The dotted line is what you get by joining raw MODIS counts to raw VIIRS counts. "
        "It jumps when VIIRS begins because VIIRS sees smaller fires, not because burning changed. "
        "AgniCal puts all years on the VIIRS scale. Error bars on reconstructed years come from "
        f"the cross-validated season-total errors (10th-90th percentile: {q10:+.0%} to {q90:+.0%}).")

# ---------------- map ----------------
with tab3:
    sp = season_periods(season)
    cells = h[in_bbox(h, bbox) & (h["year"] == year) & h["period"].isin(sp)]
    cells = cells.groupby(["lat", "lon"], as_index=False)["viirs_equiv"].sum()
    cells = cells[cells["viirs_equiv"] >= 1]
    if cells.empty:
        st.info("No detections in this area for the selected season.")
    else:
        fig = px.scatter_map(cells, lat="lat", lon="lon", color="viirs_equiv", size="viirs_equiv",
                             color_continuous_scale="YlOrRd", size_max=14, zoom=6,
                             center=dict(lat=(bbox[1] + bbox[3]) / 2, lon=(bbox[0] + bbox[2]) / 2),
                             map_style="carto-positron",
                             labels={"viirs_equiv": "Detections"})
        fig.update_layout(height=600, margin=dict(l=0, r=0, t=0, b=0))
        st.plotly_chart(fig, width="stretch")
        st.caption(f"Each dot is a {meta['grid_res_deg']}° grid cell (about 11 km), "
                   f"{season.split(' (')[0]} season {year}.")

# ---------------- method ----------------
with tab4:
    st.subheader("The problem")
    st.markdown(
        "MODIS (Terra and Aqua, about 1 km pixels) has recorded fires since 2000 and is near the end of "
        "its life. VIIRS (Suomi NPP, 375 m pixels) started in 2012 and detects many more, smaller fires. "
        "Raw counts from the two cannot be compared, so long-term trends break in 2012.")
    st.subheader("The method")
    st.markdown(
        "1. Grid all detections to 0.1° cells and 8-day periods.\n"
        "2. In the overlap years, learn how many VIIRS detections each cell-period has, given its MODIS "
        "evidence: Terra and Aqua day/night counts, fire radiative power, neighbouring cells and adjacent periods. "
        "Model: gradient-boosted trees with a Poisson loss, since the target is a count.\n"
        "3. Validate by leaving out one whole year at a time and predicting it.\n"
        "4. Apply the model to the MODIS-only years to build one VIIRS-equivalent record.")
    st.subheader("Validation: leave-one-year-out")
    rows = []
    for k, v in metrics["models"].items():
        x = v["seasons"][season]
        rows.append({"Model": LABELS.get(k, k), "Cell Poisson deviance": round(v["cell_poisson_deviance"], 3),
                     "Cell RMSE": round(v["cell_rmse"], 3), "Season total error": f"{x['total_mape_pct']:.1f}%",
                     "Bias": f"{x['total_bias_pct']:+.1f}%", "Peak error (days)": round(x["peak_err_days_mean"], 1)})
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    st.caption(f"Folds: {metrics['years'][0]}-{metrics['years'][-1]}, whole region. Lower is better everywhere.")

    ycv = st.selectbox("Held-out year to inspect", sorted(cv["year"].unique(), reverse=True))
    d = cv[cv["year"] == ycv]
    xl = [doy_to_label(p * PERIOD_DAYS + 1) for p in d["period"]]
    fig = go.Figure()
    fig.add_scatter(x=xl, y=d["v_count"], name="Observed VIIRS", line=dict(color="#212529", width=3))
    fig.add_scatter(x=xl, y=d["pred_gbm_poisson"], name="AgniCal model (from MODIS only)",
                    line=dict(color=FIRE))
    fig.add_scatter(x=xl, y=d["pred_global_ratio"], name="Regional ratio baseline",
                    line=dict(color=MUTED, dash="dot"))
    fig.update_layout(height=360, legend=dict(orientation="h", y=1.15), margin=dict(l=10, r=10, t=30, b=10),
                      yaxis_title="Detections per 8 days")
    st.plotly_chart(fig, width="stretch")

    a, b = st.columns(2)
    with a:
        fig = px.line(stab, x="year", y="ratio", markers=True,
                      title="Observed VIIRS / MODIS ratio by year")
        fig.update_traces(line_color=FIRE)
        fig.update_layout(height=320, margin=dict(l=10, r=10, t=40, b=10), yaxis_title="Ratio")
        st.plotly_chart(fig, width="stretch")
        st.caption("A stable ratio supports extending the relationship back to 2003-2011.")
    with b:
        fi = pd.Series(meta["feature_importance"]).sort_values()
        fig = px.bar(x=fi.values, y=fi.index, orientation="h", title="What the model relies on")
        fig.update_traces(marker_color=FIRE)
        fig.update_layout(height=320, margin=dict(l=10, r=10, t=40, b=10), xaxis_title="Share of gain",
                          yaxis_title=None)
        st.plotly_chart(fig, width="stretch")
    st.subheader("Limitations")
    st.markdown(
        "- Active-fire detections are not burned area; clouds, haze and short-lived fires between "
        "overpasses are missed by both sensors.\n"
        "- MODIS Terra and Aqua orbits drifted late in the mission, shifting overpass times. "
        "The ratio chart above is the check for this.\n"
        "- Reconstructed years assume the MODIS-VIIRS relationship learned in 2012 onward also held earlier.\n"
        "- Areas are bounding boxes, not administrative boundaries.")
    st.caption(f"Data: NASA FIRMS ({', '.join(meta['data']['by_sensor'])}). "
               f"Detections used: {meta['data']['rows_kept']:,}.")
