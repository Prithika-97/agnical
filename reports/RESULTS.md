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

Paired Wilcoxon test, model vs regional ratio, absolute season-total error by year (p above 0.05: the two are not distinguishable):

- Kharif (paddy residue, Oct-Nov): p = 0.32
- Rabi (wheat residue, Apr-May): p = 0.03

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

# Timing of burning: evidence and controls

Paddy season, 1 Oct - 30 Nov. 'x baseline' = value divided by the region's own 2012-2022 level.

## Indicators by region (same satellites, same passes)

| Region | Year | Suomi NPP day fires | Aqua / Suomi NPP, x baseline | Night share, x baseline | Terra / Aqua, x baseline (not used) | Flags |
|---|---|---|---|---|---|---|
| Punjab + Haryana (India) | 2021 | 84,033 | 1.04 | 0.91 | 1.18 | 0 |
| Punjab + Haryana (India) | 2022 | 49,074 | 1.10 | 2.30 | 0.51 | 0 |
| Punjab + Haryana (India) | 2023 | 32,436 | 1.34 | 3.39 | 0.16 | 2 |
| Punjab + Haryana (India) | 2024 | 9,725 | 2.09 | 20.78 | 0.11 | 2 |
| Punjab + Haryana (India) | 2025 | 4,564 | 3.61 | 31.89 | 0.08 | 2 |
| Punjab, Pakistan (border strip) | 2021 | 5,695 | 1.03 | 0.79 | 1.29 | 0 |
| Punjab, Pakistan (border strip) | 2022 | 2,854 | 1.16 | 1.68 | 0.37 | 0 |
| Punjab, Pakistan (border strip) | 2023 | 2,537 | 1.47 | 1.40 | 0.15 | 1 |
| Punjab, Pakistan (border strip) | 2024 | 3,104 | 1.04 | 2.09 | 0.16 | 0 |
| Punjab, Pakistan (border strip) | 2025 | 3,110 | 0.98 | 1.48 | 0.02 | 0 |
| Rajasthan (northern districts, India) | 2021 | 866 | 1.17 | 2.16 | 0.85 | 1 |
| Rajasthan (northern districts, India) | 2022 | 629 | 1.04 | 0.35 | 0.83 | 0 |
| Rajasthan (northern districts, India) | 2023 | 612 | 1.38 | 0.98 | 0.26 | 1 |
| Rajasthan (northern districts, India) | 2024 | 911 | 1.49 | 1.54 | 0.32 | 2 |
| Rajasthan (northern districts, India) | 2025 | 970 | 1.48 | 0.79 | 0.30 | 1 |
| Uttar Pradesh (western districts, India) | 2021 | 446 | 0.82 | 0.74 | 0.62 | 0 |
| Uttar Pradesh (western districts, India) | 2022 | 213 | 0.51 | 1.19 | 2.37 | 0 |
| Uttar Pradesh (western districts, India) | 2023 | 206 | 1.05 | 2.33 | 0.30 | 1 |
| Uttar Pradesh (western districts, India) | 2024 | 267 | 0.99 | 1.45 | 0.12 | 0 |
| Uttar Pradesh (western districts, India) | 2025 | 246 | 1.95 | 1.48 | 0.13 | 1 |

Terra / Aqua fell as much across the border as in the study area, so it reflects the satellites' orbit drift and is not used as evidence.

## Study area relative to the comparison region across the border

| Indicator | Year | Study area x baseline | Comparison x baseline | Ratio | 95% interval |
|---|---|---|---|---|---|
| aqua_over_snpp | 2023 | 1.35 | 1.47 | 0.92 | 0.68 - 1.20 |
| aqua_over_snpp | 2024 | 2.10 | 1.04 | 2.02 | 1.39 - 2.93 |
| aqua_over_snpp | 2025 | 3.63 | 0.98 | 3.70 | 2.52 - 5.38 |
| night_share | 2023 | 3.65 | 1.42 | 2.56 | 1.70 - 4.21 |
| night_share | 2024 | 22.42 | 2.13 | 10.51 | 7.35 - 15.17 |
| night_share | 2025 | 34.40 | 1.51 | 22.81 | 14.78 - 37.77 |

## Near-nadir check (small pixels only)

| Year | Aqua / Suomi NPP, all | near nadir | all, x baseline | near nadir, x baseline |
|---|---|---|---|---|
| 2020 | 0.216 | 0.308 | 1.09 | 1.04 |
| 2021 | 0.208 | 0.306 | 1.04 | 1.03 |
| 2022 | 0.220 | 0.339 | 1.10 | 1.15 |
| 2023 | 0.268 | 0.370 | 1.34 | 1.25 |
| 2024 | 0.416 | 0.542 | 2.09 | 1.83 |
| 2025 | 0.719 | 0.937 | 3.61 | 3.17 |

## Late-view floor

Aqua's daytime count converted to Suomi NPP units with the 2012-2022 ratio (5.05).

| Year | Aqua looks at | Suomi NPP counted | Aqua in Suomi NPP units | Range | Factor |
|---|---|---|---|---|---|
| 2020 | 13:47 | 75,202 | 82,135 | 74,055 - 93,794 | 1.09 |
| 2021 | 13:48 | 84,011 | 88,104 | 79,436 - 100,609 | 1.05 |
| 2022 | 13:52 | 49,068 | 54,422 | 49,068 - 62,147 | 1.11 |
| 2023 | 14:07 | 32,438 | 43,848 | 39,534 - 50,071 | 1.35 |
| 2024 | 14:37 | 9,727 | 20,377 | 18,373 - 23,270 | 2.09 |
| 2025 | 15:15 | 4,571 | 16,563 | 14,934 - 18,914 | 3.62 |

## A test that did not show the shift: VIIRS alone, inside its own window

Days with VIIRS passes both before and after 13:30. The rise per hour through the window did not leave its 2012-2022 range, so what still burns between about 12:40 and 14:20 keeps its old shape; the change lies outside the window.

| Year | Days | Early | Late | Log rise per hour |
|---|---|---|---|---|
| 2012 | 23 | 11,979 | 18,041 | 0.24 |
| 2013 | 20 | 7,056 | 12,427 | 0.34 |
| 2014 | 22 | 12,346 | 19,095 | 0.26 |
| 2015 | 25 | 8,251 | 18,107 | 0.47 |
| 2016 | 25 | 10,833 | 24,006 | 0.47 |
| 2017 | 22 | 6,771 | 10,657 | 0.27 |
| 2018 | 21 | 7,381 | 16,542 | 0.48 |
| 2019 | 22 | 3,037 | 11,903 | 0.82 |
| 2020 | 23 | 11,341 | 22,091 | 0.40 |
| 2021 | 25 | 15,160 | 19,498 | 0.15 |
| 2022 | 21 | 6,371 | 12,387 | 0.40 |
| 2023 | 50 | 21,257 | 36,690 | 0.46 |
| 2024 | 47 | 5,175 | 9,113 | 0.50 |
| 2025 | 48 | 2,705 | 5,354 | 0.61 |
| 2026 | 12 | 341 | 382 | 0.07 |

## VIIRS hand-over after Suomi NPP (ended 2 Nov 2026)

| Satellite | Year | Days both observed | Suomi NPP per detection |
|---|---|---|---|
| NOAA21 | 2023 | 58 | 1.013 |
| NOAA21 | 2024 | 57 | 1.119 |
| NOAA21 | 2025 | 57 | 1.023 |
