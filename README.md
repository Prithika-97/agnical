# AgniCal

**Punjab and Haryana's farm fires moved later in the day. The satellites that count them didn't.**

NASA Space Apps Challenge 2026 · Challenge: *Harmonization of MODIS and VIIRS Hot Spots*

Crop-residue fires in Punjab and Haryana feed north India's winter smog, and official
claims of progress rest on one number: how many fires NASA's satellites count. Those
satellites (MODIS on Terra and Aqua, VIIRS on Suomi NPP, NOAA-20 and NOAA-21) pass over at
fixed hours. AgniCal does two things:

1. **One scale for 2003-2026.** It harmonizes MODIS into VIIRS-equivalent counts, so the
   2012 switch from MODIS to VIIRS no longer looks like a jump in burning, and turns the
   record into a burning calendar for every state and district.
2. **A timing break nobody had flagged in the record.** Using only NASA's own polar-orbiter
   data, it shows that burning moved later than the satellites look after 2022, so a fall in
   counts can no longer be read as a fall in burning.

## What the data show (computed by `pipeline.py`, see `reports/RESULTS.md`)

- Suomi NPP counted about 95% fewer paddy-season fires in 2025 than in its 2016 peak.
- In the same seasons, three model-free indicators all left their 2012-2022 range:
  - **Aqua**, whose orbit drifted from ~13:50 to ~15:15 IST, now sees several times as many
    fires per Suomi NPP fire as it did when both passed at the same hour.
  - The **share of fires seen on the ~01:30 night pass** rose many-fold, and night-pass
    counts rose in absolute terms while daytime counts collapsed.
  - **Terra** (late morning) now sees almost nothing.
- A day-matched model of all five satellites puts the peak of visible burning at about
  13:40 in 2012-2022, in line with the ~13:30 that geostationary data showed for 2020; it
  moves later after 2022 and leaves the satellites' view in 2025. Independent geostationary
  studies (ISRO SAC with Meteosat; NASA with GEO-KOMPSAT-2A) report peaks of 16:00-18:00.
- Trained on 2012-2022, MODIS alone overstates what Suomi NPP counts by a growing margin in
  2023-2025, another sign that the satellites now sample a different part of the day.

## Run it

```bash
pip install -r requirements.txt
# put NASA FIRMS archive CSVs in data/raw/ (MODIS C6.1 and VIIRS S-NPP, NOAA-20, NOAA-21)
python pipeline.py            # about 5 minutes; writes data/processed, reports/ and docs/data
python -m http.server -d docs # then open http://localhost:8000
pytest -q
```

The website in `docs/` is static (HTML, JSON, Plotly.js) and is served as-is by GitHub Pages.

## Method

| Step | Module | What it does |
|---|---|---|
| Ingest | `agnical/ingest.py` | Reads FIRMS archive and near-real-time files for all five satellites, keeps vegetation fires, drops low confidence, prefers standard over NRT data per satellite and day. |
| Borders | `agnical/regions.py` | Keeps detections inside the real borders of Punjab and Haryana (geoBoundaries) and tags each with its district. |
| Time of day | `agnical/diurnal.py` | Builds one row per satellite pass; fits n ~ Poisson(exp(α_day + β_sensor − (t − μ_year)² / 2σ_year²)), removing α_day exactly by conditioning on each day's total; day-level bootstrap. Reports a peak only where viewing times bracket it. Flags timing breaks from three model-free indicators (robust z > 3 on at least two). |
| One scale | `agnical/grid.py`, `agnical/harmonize.py` | 0.1° cells × 8-day periods; LightGBM (Poisson) predicts Suomi NPP counts from MODIS features; trained on 2012-2022 only; leave-one-season-out validation against two ratio baselines; drift test on 2023 onward. |
| Outputs | `agnical/export.py` | Writes the JSON the website reads. No number on the site is typed by hand. |

## Limitations

- Active-fire detections are snapshots, not burned area; cloud, haze and short fires between passes are missed.
- Polar satellites cannot see burning after about 15:00. Where the peak is later, AgniCal reports only a lower bound. A full correction needs geostationary data (Meteosat-9 IODC, INSAT-3D/3DR, GEO-KOMPSAT-2A).
- District borders are from 2021 (before Malerkotla and Charkhi Dadri).
- 2026 near-real-time files have no fire-type field, so a few static heat sources may remain.

## Data and credits

- NASA FIRMS active fire data (MODIS C6.1; VIIRS 375 m from Suomi NPP, NOAA-20, NOAA-21): https://firms.modaps.eosdis.nasa.gov
- geoBoundaries gbOpen IND: states (DataMeet / Election Commission of India, CC BY 2.5 IN), districts (Pathways Data / lgdirectory.gov.in, ODbL 1.0)
- Published peak-hour estimates used only as outside checks: ISRO SAC study reported in The Tribune (Dec 2025); NASA / H. Jethva reported by IANS (Dec 2025).

## Use of AI

Code and text were drafted with Claude (Anthropic) and reviewed, run and checked by the author.
All figures are computed from NASA data by this pipeline. No AI-generated images or data are used.
