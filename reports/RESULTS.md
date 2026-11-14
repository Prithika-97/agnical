# AgniCal validation results

Leave-one-year-out cross-validation over overlap years 2012-2025 (14 folds, 1,146,749 cell-periods).

## Cell level (one cell, one 8-day period)

| Model | RMSE | MAE | Poisson deviance |
|---|---|---|---|
| Baseline: one regional ratio | 3.637 | 0.736 | 6.166 |
| Baseline: per-cell ratio | 3.560 | 0.735 | 6.177 |
| AgniCal model (LightGBM, Poisson) | 2.957 | 0.614 | 0.801 |

## Regional series (sum over the region, per 8-day period)

| Model | Correlation | RMSE |
|---|---|---|
| Baseline: one regional ratio | 0.968 | 1244.3 |
| Baseline: per-cell ratio | 0.969 | 1231.3 |
| AgniCal model (LightGBM, Poisson) | 0.981 | 951.4 |

## Kharif (paddy residue, Oct-Nov)

| Model | Season total error (MAPE) | Bias | Peak timing error | Peak exact |
|---|---|---|---|---|
| Baseline: one regional ratio | 15.6% | +10.4% | 3.4 days | 57% of years |
| Baseline: per-cell ratio | 15.4% | +7.6% | 3.4 days | 57% of years |
| AgniCal model (LightGBM, Poisson) | 14.9% | +7.7% | 1.1 days | 86% of years |

## Rabi (wheat residue, Apr-May)

| Model | Season total error (MAPE) | Bias | Peak timing error | Peak exact |
|---|---|---|---|---|
| Baseline: one regional ratio | 17.9% | +1.3% | 2.3 days | 71% of years |
| Baseline: per-cell ratio | 18.9% | +3.3% | 2.3 days | 71% of years |
| AgniCal model (LightGBM, Poisson) | 12.1% | -0.9% | 1.7 days | 79% of years |

## Is the MODIS-VIIRS relationship stable over time?

| Year | VIIRS detections | MODIS detections | Ratio |
|---|---|---|---|
| 2012 | 112,125 | 26,147 | 4.29 |
| 2013 | 100,546 | 22,825 | 4.41 |
| 2014 | 113,320 | 22,726 | 4.99 |
| 2015 | 95,649 | 18,244 | 5.24 |
| 2016 | 130,153 | 26,714 | 4.87 |
| 2017 | 97,571 | 19,117 | 5.10 |
| 2018 | 103,218 | 22,513 | 4.58 |
| 2019 | 83,388 | 18,545 | 4.50 |
| 2020 | 102,905 | 22,379 | 4.60 |
| 2021 | 118,376 | 25,339 | 4.67 |
| 2022 | 83,067 | 18,513 | 4.49 |
| 2023 | 59,598 | 12,977 | 4.59 |
| 2024 | 50,385 | 11,700 | 4.31 |
| 2025 | 39,121 | 9,248 | 4.23 |

A drifting ratio means extrapolating to 2003-2011 carries extra uncertainty; the app shows a band for reconstructed years from the cross-validation errors.

## Feature importance (share of gain)

- m_nbr: 62.4%
- m_aqua_day: 20.3%
- m_count: 8.4%
- m_frp_sum: 3.2%
- period: 1.6%
- lat: 1.1%
- m_prev: 1.0%
- lon: 1.0%
- m_next: 0.8%
- m_frp_max: 0.1%
- m_aqua_night: 0.1%
- m_terra_night: 0.0%
- m_terra_day: 0.0%
