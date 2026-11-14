"""Project-wide settings. Change the region here; nothing else is region-specific."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"
OUT_DIR = ROOT / "data" / "processed"
REPORT_DIR = ROOT / "reports"

# Region of interest: Punjab + Haryana, India (west, south, east, north), degrees.
BBOX = (73.8, 27.6, 77.6, 32.6)

# Approximate sub-areas offered in the app. These are bounding boxes, not
# administrative boundaries, so each one includes some neighbouring land.
AREA_PRESETS = {
    "Whole region (Punjab + Haryana box)": BBOX,
    "Punjab (approx. box)": (73.8, 29.5, 76.95, 32.6),
    "Haryana (approx. box)": (74.4, 27.6, 77.6, 30.95),
}

GRID_RES = 0.1          # degrees; ~11 km cells
N_PERIODS = 46          # 8-day periods per year, same convention as MODIS composites

MODIS_FIRST_YEAR = 2003            # first full year with both Terra and Aqua
VIIRS_START = "2012-01-20"         # Suomi NPP 375 m record starts
VIIRS_SATELLITES = ("SNPP",)       # one satellite keeps the target record consistent
MODIS_CONF_MIN = 30                # drop low-confidence MODIS detections
VIIRS_DROP_LOW = True              # drop VIIRS 'low' confidence detections
KEEP_TYPES = (0,)                  # 0 = presumed vegetation fire (archive files only)

# Burning seasons, by day of year of each 8-day period's first day.
SEASONS = {
    "Kharif (paddy residue, Oct-Nov)": (258, 349),   # 15 Sep - 15 Dec
    "Rabi (wheat residue, Apr-May)": (91, 151),       # 1 Apr - 31 May
}

ANOMALY_Z = 2.0   # |z| above this flags a season as unusual
