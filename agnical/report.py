"""Results summary written from the computed metrics."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

LABELS = {"global_ratio": "Baseline: one regional ratio",
          "cell_ratio": "Baseline: per-cell ratio",
          "gbm_poisson": "AgniCal model (LightGBM, Poisson)"}


def write_results(meta: dict, stability: pd.DataFrame, path: Path, timing=None, ctrl=None) -> None:
    m = meta["metrics"]
    lines = ["# AgniCal validation results", ""]
    if meta["data"]["synthetic"]:
        lines += ["> **SYNTHETIC TEST DATA.** These numbers only show the pipeline runs. "
                  "They are not findings about real fires.", ""]
    lines += [f"Leave-one-year-out cross-validation over overlap years "
              f"{m['years'][0]}-{m['years'][-1]} ({len(m['years'])} folds, {m['n_rows']:,} cell-periods).",
              "", "## Cell level (one cell, one 8-day period)", "",
              "| Model | RMSE | MAE | Poisson deviance |", "|---|---|---|---|"]
    for k, v in m["models"].items():
        lines.append(f"| {LABELS.get(k, k)} | {v['cell_rmse']:.3f} | {v['cell_mae']:.3f} | "
                     f"{v['cell_poisson_deviance']:.3f} |")
    lines += ["", "## Regional series (sum over the region, per 8-day period)", "",
              "| Model | Correlation | RMSE |", "|---|---|---|"]
    for k, v in m["models"].items():
        lines.append(f"| {LABELS.get(k, k)} | {v['regional_corr']:.3f} | {v['regional_rmse']:.1f} |")
    seasons = list(next(iter(m["models"].values()))["seasons"])
    for s in seasons:
        lines += ["", f"## {s}", "",
                  "| Model | Season total error (MAPE) | Bias | Peak timing error | Peak exact |",
                  "|---|---|---|---|---|"]
        for k, v in m["models"].items():
            x = v["seasons"][s]
            lines.append(f"| {LABELS.get(k, k)} | {x['total_mape_pct']:.1f}% | {x['total_bias_pct']:+.1f}% | "
                         f"{x['peak_err_days_mean']:.1f} days | {x['peak_exact_pct']:.0f}% of years |")
    if m.get("total_error_test"):
        lines += ["", "Paired Wilcoxon test, model vs regional ratio, absolute season-total error by year "
                  "(p above 0.05: the two are not distinguishable):", ""]
        lines += [f"- {s}: p = {p:.2f}" for s, p in m["total_error_test"].items()]
    lines += ["", "## Is the MODIS-VIIRS relationship stable over time?", "",
              "| Year | VIIRS detections | MODIS detections | Ratio |", "|---|---|---|---|"]
    for _, r in stability.iterrows():
        lines.append(f"| {int(r['year'])} | {int(r['v_count']):,} | {int(r['m_count']):,} | {r['ratio']:.2f} |")
    lines += ["", "A drifting ratio means extrapolating to 2003-2011 carries extra uncertainty; "
              "the app shows a band for reconstructed years from the cross-validation errors.", "",
              "## Feature importance (share of gain)", ""]
    for f, g in meta["feature_importance"].items():
        lines.append(f"- {f}: {g:.1%}")
    if ctrl is not None:
        lines += _controls(ctrl, timing)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _controls(ctrl: dict, timing) -> list[str]:
    L = ["", "# Timing of burning: evidence and controls", "",
         "Paddy season, 1 Oct - 30 Nov. 'x baseline' = value divided by the region's own 2012-2022 level.", "",
         "## Indicators by region (same satellites, same passes)", "",
         "| Region | Year | Suomi NPP day fires | Aqua / Suomi NPP, x baseline | Night share, x baseline | "
         "Terra / Aqua, x baseline (not used) | Flags |", "|---|---|---|---|---|---|---|"]
    rg = ctrl["regions"]
    for r in rg[rg["year"] >= 2021].itertuples():
        L.append(f"| {r.region} | {r.year} | {int(r.snpp_day):,} | {r.aqua_over_snpp_x_baseline:.2f} | "
                 f"{r.night_share_x_baseline:.2f} | {r.terra_over_aqua_x_baseline:.2f} | {int(r.n_flags)} |")
    L += ["", "Terra / Aqua fell as much across the border as in the study area, so it reflects the satellites' "
          "orbit drift and is not used as evidence.", "",
          "## Study area relative to the comparison region across the border", "",
          "| Indicator | Year | Study area x baseline | Comparison x baseline | Ratio | 95% interval |",
          "|---|---|---|---|---|---|"]
    for r in ctrl["did"].itertuples():
        L.append(f"| {r.indicator} | {r.year} | {r.study_x_baseline:.2f} | {r.control_x_baseline:.2f} | "
                 f"{r.ratio:.2f} | {r.ratio_lo:.2f} - {r.ratio_hi:.2f} |")
    L += ["", "## Near-nadir check (small pixels only)", "", "| Year | Aqua / Suomi NPP, all | near nadir | "
          "all, x baseline | near nadir, x baseline |", "|---|---|---|---|---|"]
    for r in ctrl["nadir"].tail(6).itertuples():
        L.append(f"| {r.year} | {r.all:.3f} | {r.near_nadir:.3f} | {r.all_x_baseline:.2f} | {r.near_nadir_x_baseline:.2f} |")
    L += ["", "## Late-view floor", "", "Aqua's daytime count converted to Suomi NPP units with the 2012-2022 ratio "
          f"({ctrl['late_view']['baseline_ratio'].iloc[0]:.2f}).", "",
          "| Year | Aqua looks at | Suomi NPP counted | Aqua in Suomi NPP units | Range | Factor |", "|---|---|---|---|---|---|"]
    for r in ctrl["late_view"].tail(6).itertuples():
        h = r.aqua_time_ist
        L.append(f"| {r.year} | {int(h):02d}:{int(round((h % 1) * 60)) % 60:02d} | {int(r.snpp_day):,} | "
                 f"{r.aqua_as_snpp:,.0f} | {r.aqua_as_snpp_lo:,.0f} - {r.aqua_as_snpp_hi:,.0f} | {r.factor:.2f} |")
    L += ["", "## A test that did not show the shift: VIIRS alone, inside its own window", "",
          "Days with VIIRS passes both before and after 13:30. The rise per hour through the window did not "
          "leave its 2012-2022 range, so what still burns between about 12:40 and 14:20 keeps its old shape; "
          "the change lies outside the window.", "", "| Year | Days | Early | Late | Log rise per hour |", "|---|---|---|---|---|"]
    for r in ctrl["viirs_window"].itertuples():
        L.append(f"| {r.year} | {r.days} | {r.early:,} | {r.late:,} | {r.log_rise_per_h:.2f} |")
    if len(ctrl["handover"]):
        L += ["", "## VIIRS hand-over after Suomi NPP (ended 2 Nov 2026)", "",
              "| Satellite | Year | Days both observed | Suomi NPP per detection |", "|---|---|---|---|"]
        for r in ctrl["handover"].itertuples():
            L.append(f"| {r.satellite} | {r.year} | {r.days} | {r.snpp_per_new:.3f} |")
    return L
