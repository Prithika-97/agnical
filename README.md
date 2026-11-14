# AgniCal

**Punjab and Haryana's farm fires moved later in the day. The satellites that count them didn't.**

NASA Space Apps Challenge 2026 · Challenge: *Harmonization of MODIS and VIIRS Hot Spots*

Crop-residue fires in Punjab and Haryana feed north India's winter smog, and claims of
progress rest on one number: how many fires NASA's polar satellites count. Those satellites
(MODIS on Terra and Aqua; VIIRS on Suomi NPP, NOAA-20 and NOAA-21) pass at fixed hours.

That farmers now burn later in the day is **not our discovery**. Jethva (2026, *Remote
Sensing*), Singh et al. (ISRO SAC, 2025, *Current Science*) and iFOREST (2025) showed it with
geostationary satellites and burnt-area mapping. AgniCal asks what follows for the harmonized
MODIS-VIIRS record that everyone uses, and answers with NASA's polar data alone:

1. **One scale for 2003-2025**, with uncertainty. MODIS is converted to VIIRS-equivalent
   counts so the 2012 switch no longer looks like a jump in burning.
2. **A record that knows when it stops being comparable.** Two model-free indicators flag
   seasons in which burning left the satellites' viewing hours.
3. **A control for the satellites themselves.** Aqua and Terra both drifted in the same years.
   Every test is repeated on the farm fires just across the border in Pakistani Punjab, which
   the same satellites see on the same passes.
4. **A floor under the real count.** Aqua's later pass, in Suomi NPP units.

## What the data show (all computed by `pipeline.py`; see `reports/RESULTS.md`)

| | Punjab + Haryana, India | Punjab, Pakistan (border strip) |
|---|---|---|
| Suomi NPP daytime fires, 2025 vs 2012-22 average | -94% | -5% |
| Aqua fires per Suomi NPP fire, 2025, x its 2012-22 level | 3.6 | 1.0 |
| Share of fires seen on the night pass, 2025, x its 2012-22 level | 34 | 1.5 |
| Terra per Aqua fire, 2025, x its 2012-22 level (**not used**) | 0.08 | 0.02 |

- Relative to the comparison region, the Aqua ratio on the Indian side is **3.7x** in 2025
  (95% interval 2.5-5.4) and 2.0x in 2024 (1.4-2.9). The night share is 23x (15-38) and
  10.5x (7.4-15).
- **2023 is unsettled**: Aqua's ratio rose as much across the border (0.92, 0.68-1.20). Only
  the night share separates the two sides that year. 2024 and 2025 are confirmed.
- **A signal we threw out**: the morning (Terra) to afternoon (Aqua) ratio collapsed on both
  sides of the border, so it reflects orbit drift, not burning.
- **A test that found nothing**: within VIIRS's own window (two looks ~100 min apart on some
  days) the rise per hour stayed inside its 2012-2022 range. What still burns at midday keeps
  its old shape; the change lies after the window.
- **Viewing angle** does not explain it: with near-nadir detections only, the Aqua ratio in
  2025 is 3.2x its old level (3.6x with all detections).
- **Floor**: in 2025 Suomi NPP counted about 4,600 daytime fires (1 Oct - 30 Nov). At Aqua's
  hour (15:15) it would have counted about 16,600 (14,900-18,900). Against the 2016 peak the
  count says -95%; the floor says no more than about -83%. iFOREST's Sentinel-2 burnt area for
  Punjab fell 36% between 2022 and 2025.

### Harmonization, validated honestly (leave one season out, 2012-2022)

| Method | Cell deviance | Paddy total error | Wheat total error | Paddy peak period exact |
|---|---|---|---|---|
| One regional ratio | 7.29 | 5.5% | 22.5% | 55% of seasons |
| Per-cell ratio | 7.29 | 5.3% | 22.8% | 64% |
| LightGBM, Poisson | 0.93 | 6.6% | 14.2% | 91% |

For the paddy-season regional total a simple ratio is as good as the model (paired Wilcoxon
p = 0.32); all methods miss the same seasons. The model is better where a ratio cannot be:
cell-level maps, the timing of the seasonal peak, and the wheat season (p = 0.03).

Suomi NPP stopped delivering data on 2 Nov 2026. NOAA-21 counts convert to Suomi NPP units by
a factor of 1.03 (yearly range 1.01-1.12, 2023-2025).

## Run it

```bash
pip install -r requirements.txt
# put NASA FIRMS archive CSVs in data/raw/ (MODIS C6.1 and VIIRS S-NPP, NOAA-20, NOAA-21)
python pipeline.py            # about 3 minutes; writes data/processed, reports/ and docs/data
python -m http.server -d docs # then open http://localhost:8000
pytest -q
```

The website in `docs/` is static (HTML, JSON, Plotly.js) and is served as-is by GitHub Pages.
The harmonized record is one file: `docs/data/agnical_record.csv`.

## Method

| Step | Module | What it does |
|---|---|---|
| Ingest | `agnical/ingest.py` | Reads FIRMS archive and near-real-time files, keeps vegetation fires, drops low confidence, prefers standard over NRT data per satellite and day. |
| Borders | `agnical/regions.py` | Keeps detections inside Punjab and Haryana (geoBoundaries) and tags each with its district. |
| Controls | `agnical/controls.py` | Comparison regions; difference-in-differences with a day-level bootstrap; near-nadir check; late-view floor; VIIRS-only window test; VIIRS hand-over factor. |
| Time of day | `agnical/diurnal.py` | One row per satellite pass; n ~ Poisson(exp(a_day + b_sensor - (t - mu_year)^2 / 2 sigma_year^2)), with a_day removed by conditioning on each day's total; day-level bootstrap. Reports a peak only where viewing times bracket it. Timing-break test on two indicators (robust z > 3 on both). |
| One scale | `agnical/grid.py`, `agnical/harmonize.py` | 0.1 degree cells x 8-day periods; LightGBM (Poisson) predicts Suomi NPP counts from MODIS features; trained on 2012-2022 only; leave-one-season-out validation against two ratio baselines; drift test on 2023 onward. |
| Outputs | `agnical/export.py`, `agnical/report.py` | Write the JSON and CSV the website reads and `reports/RESULTS.md`. No number on the site is typed by hand. |

To run another region, change `BBOX` in `agnical/config.py` and the two boundary files built
by `regions.build_boundary_file` and `controls.build_region_file`.

## Limitations

- Active-fire detections are snapshots, not burned area; cloud, haze and short fires between passes are missed.
- Polar satellites cannot see burning after about 15:15. The floor is a lower bound, not an estimate of the total. We did not obtain geostationary fire data during the hackathon; the 16:00-18:00 peak is quoted from published studies.
- The comparison region is a narrow strip with about a tenth of the study area's fires. It controls for the satellites, not for weather or crop-calendar differences.
- Both indicators use Suomi NPP's daytime count, so they are not fully independent.
- Only one region is analysed. The timing test ran unchanged on four regions in the download box, but those are neighbours, not a global sample.
- District borders are from 2021 (before Malerkotla and Charkhi Dadri).
- 2026 near-real-time files have no fire-type field, so a few static heat sources may remain.

## Data and credits

- NASA FIRMS active fire data (MODIS C6.1; VIIRS 375 m): https://firms.modaps.eosdis.nasa.gov
- geoBoundaries gbOpen IND: states (DataMeet / Election Commission of India, CC BY 2.5 IN), districts (Pathways Data / lgdirectory.gov.in, ODbL 1.0)
- Jethva, H. (2026). Timing the Flames: Geostationary Satellite Detection of Diurnally Shifting Stubble Burning in Northwestern India. *Remote Sensing* 18(10), 1506. https://doi.org/10.3390/rs18101506
- Singh, N., Pradhan, R., Shukla, B. P., Pandya, M. R. (Space Applications Centre, ISRO), communication in *Current Science*, 2025, as reported by The Tribune, 10 Dec 2025.
- iFOREST (2025). Stubble Burning Status Report 2025, as reported by The Tribune, 9 Dec 2025.

## Use of AI

Code and text were drafted with Claude (Anthropic) and reviewed, run and checked by the author.
All figures are computed from NASA data by this pipeline. No AI-generated images or data are used.
