# AgniCal validation results

Leave-one-year-out cross-validation over overlap years 2012-2022 (11 folds, 497,970 cell-periods).

## Cell level (one cell, one 8-day period)

| Model | RMSE | MAE | Poisson deviance |
|---|---|---|---|
| Baseline: one regional ratio | 4.844 | 1.102 | 7.292 |
| Baseline: per-cell ratio | 4.722 | 1.092 | 7.290 |
| AgniCal model (LightGBM, Poisson) | 3.790 | 0.870 | 0.935 |

## Regional series (sum over the region, per 8-day period)

| Model | Correlation | RMSE |
|---|---|---|
| Baseline: one regional ratio | 0.971 | 1243.4 |
| Baseline: per-cell ratio | 0.972 | 1223.7 |
| AgniCal model (LightGBM, Poisson) | 0.984 | 916.5 |

## Kharif (paddy residue, Oct-Nov)

| Model | Season total error (MAPE) | Bias | Peak timing error | Peak exact |
|---|---|---|---|---|
| Baseline: one regional ratio | 5.5% | -1.8% | 3.6 days | 55% of years |
| Baseline: per-cell ratio | 5.3% | -2.0% | 2.9 days | 64% of years |
| AgniCal model (LightGBM, Poisson) | 6.6% | -0.6% | 0.7 days | 91% of years |

## Rabi (wheat residue, Apr-May)

| Model | Season total error (MAPE) | Bias | Peak timing error | Peak exact |
|---|---|---|---|---|
| Baseline: one regional ratio | 22.5% | +11.7% | 2.9 days | 64% of years |
| Baseline: per-cell ratio | 22.8% | +12.0% | 2.9 days | 64% of years |
| AgniCal model (LightGBM, Poisson) | 14.2% | +0.3% | 2.2 days | 73% of years |

## Is the MODIS-VIIRS relationship stable over time?

| Year | VIIRS detections | MODIS detections | Ratio |
|---|---|---|---|
| 2012 | 94,223 | 22,728 | 4.15 |
| 2013 | 88,242 | 20,242 | 4.36 |
| 2014 | 100,135 | 20,489 | 4.89 |
| 2015 | 85,696 | 16,636 | 5.15 |
| 2016 | 112,679 | 23,676 | 4.76 |
| 2017 | 81,662 | 16,202 | 5.04 |
| 2018 | 83,339 | 18,960 | 4.40 |
| 2019 | 71,541 | 16,341 | 4.38 |
| 2020 | 92,205 | 20,472 | 4.50 |
| 2021 | 101,381 | 22,518 | 4.50 |
| 2022 | 68,337 | 16,082 | 4.25 |
| 2023 | 49,316 | 11,299 | 4.36 |
| 2024 | 30,429 | 8,132 | 3.74 |
| 2025 | 25,581 | 7,091 | 3.61 |
| 2026 | 0 | 0 | nan |

A drifting ratio means extrapolating to 2003-2011 carries extra uncertainty; the app shows a band for reconstructed years from the cross-validation errors.

## Feature importance (share of gain)

- m_nbr: 60.6%
- m_aqua_day: 17.6%
- m_count: 13.5%
- m_frp_sum: 3.4%
- period: 1.4%
- m_prev: 0.8%
- lon: 0.7%
- m_next: 0.7%
- lat: 0.7%
- m_frp_max: 0.6%
- m_terra_day: 0.0%
- m_terra_night: 0.0%
- m_aqua_night: 0.0%
