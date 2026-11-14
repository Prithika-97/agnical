# AgniCal: a consistent 20+ year burning calendar from MODIS and VIIRS

NASA Space Apps Challenge 2026, challenge **Harmonization of MODIS and VIIRS Hot Spots**.

Satellites have tracked active fires for over 20 years, but the record is split between
MODIS (Terra/Aqua, ~1 km pixels, since 2000) and VIIRS (Suomi NPP, 375 m pixels, since 2012).
VIIRS sees many more, smaller fires, so joining the two raw records creates a false jump in 2012.
AgniCal learns how MODIS detections map to VIIRS-equivalent detections during the overlap years,
validates that mapping one held-out year at a time, and rebuilds one consistent record.
The showcase region is Punjab and Haryana, India, where crop-residue burning drives severe
seasonal air pollution.

## What the app shows

- **Burning calendar:** every 8-day period of every year on one scale.
- **Peak week, burning window and an "unusual season" flag** for any area and season.
- **Long-term trend** with the naive stitched record for comparison, and uncertainty bands
  for reconstructed (pre-VIIRS) years.
- **Method and validation:** model vs baselines, held-out years, ratio stability, feature importance.

## Quick start

```bash
pip install -r requirements.txt
# 1. Put FIRMS archive CSVs (MODIS C6.1 and VIIRS S-NPP 375 m) in data/raw/
# 2. Build everything (about 5-10 minutes; add --fast for a quick first look)
python pipeline.py --end 2025-12-31
# 3. Open the app
streamlit run app.py
```

No data yet? `python scripts/make_synthetic.py --out data/synthetic` writes clearly labelled
synthetic files; run `python pipeline.py --raw data/synthetic --fast` to test the pipeline.
Everything built from them is flagged as SYNTHETIC in the results file and the app.

Run the tests with `pytest -q`.

## Method

1. **Ingest** (`agnical/ingest.py`): read FIRMS archive and near-real-time CSVs, keep vegetation
   fires (type 0), drop low-confidence detections, keep VIIRS Suomi NPP only, prefer standard over
   NRT data where both exist, remove duplicates.
2. **Grid** (`agnical/grid.py`): 0.1° cells x 8-day periods. MODIS features per cell-period:
   Terra/Aqua day/night counts, fire radiative power (sum, max), detections in the 8 neighbouring
   cells, and in the previous and next period.
3. **Harmonize** (`agnical/harmonize.py`): LightGBM with a Poisson objective predicts VIIRS
   detections from MODIS features. Compared against two baselines: one regional ratio, and a
   per-cell ratio shrunk toward the regional ratio.
4. **Validate:** leave-one-year-out over the overlap years. Metrics at cell level (RMSE, MAE,
   Poisson deviance), regional level (correlation, RMSE), season totals (MAPE, bias) and peak
   timing (days). `reports/RESULTS.md` is written automatically from these numbers.
5. **Calendar** (`agnical/burncal.py`): season totals, peak week, 10-90% burning window, and a
   z-score against all other years.

## Data

- NASA FIRMS active fire archive: MODIS Collection 6.1 (MCD14ML) and VIIRS S-NPP 375 m
  (VNP14IMGML). https://firms.modaps.eosdis.nasa.gov/download/

## Limitations

- Active-fire detections are not burned area. Clouds, haze and fires between overpasses are missed.
- MODIS orbits drifted late in the mission; the ratio-stability check is reported for this reason.
- Reconstructed years assume the relationship learned from 2012 onward held before 2012.
- App areas are bounding boxes, not administrative boundaries.

## Use of AI

Code was written with help from Claude (Anthropic), an AI assistant, and reviewed, run and tested
by the author. No AI-generated images or data are used; the synthetic test generator is labelled
and its outputs are never presented as results.
