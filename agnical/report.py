"""Write a results summary straight from the computed metrics, so no number is typed by hand."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

LABELS = {"global_ratio": "Baseline: one regional ratio",
          "cell_ratio": "Baseline: per-cell ratio",
          "gbm_poisson": "AgniCal model (LightGBM, Poisson)"}


def write_results(meta: dict, stability: pd.DataFrame, path: Path) -> None:
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
    lines += ["", "## Is the MODIS-VIIRS relationship stable over time?", "",
              "| Year | VIIRS detections | MODIS detections | Ratio |", "|---|---|---|---|"]
    for _, r in stability.iterrows():
        lines.append(f"| {int(r['year'])} | {int(r['v_count']):,} | {int(r['m_count']):,} | {r['ratio']:.2f} |")
    lines += ["", "A drifting ratio means extrapolating to 2003-2011 carries extra uncertainty; "
              "the app shows a band for reconstructed years from the cross-validation errors.", "",
              "## Feature importance (share of gain)", ""]
    for f, g in meta["feature_importance"].items():
        lines.append(f"- {f}: {g:.1%}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
