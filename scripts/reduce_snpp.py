"""Reduce the India SNPP download to counts per 0.1 deg cell, date and day/night.

Usage: python reduce_snpp.py
Reads DL_FIRE_SV-C2_847520.zip (or the .crdownload copy) next to this script, writes
india_snpp_counts.csv.gz. Standard library only.
"""
import csv, glob, gzip, io, os, sys, zipfile
from collections import Counter

here = os.path.dirname(os.path.abspath(__file__))
cands = glob.glob(os.path.join(here, "DL_FIRE_SV-C2_847520*.zip")) + glob.glob(os.path.join(here, "Unconfirmed*.crdownload"))
cands = [c for c in cands if os.path.getsize(c) > 100_000_000]
if not cands:
    sys.exit("Could not find the Suomi NPP download (about 160 MB) in " + here)
src = max(cands, key=os.path.getsize)
print("Reading", os.path.basename(src), f"({os.path.getsize(src) / 1e6:.0f} MB)")
try:
    z = zipfile.ZipFile(src)
except zipfile.BadZipFile:
    sys.exit("The download is not complete yet. Wait for Chrome to finish it, then run this again.")

counts, last_archive, n = Counter(), "", 0
names = sorted([m for m in z.namelist() if m.endswith(".csv")], key=lambda m: "nrt" in m)   # archive first
for name in names:
    nrt = "nrt" in name
    with z.open(name) as f:
        for row in csv.DictReader(io.TextIOWrapper(f, "utf-8")):
            n += 1
            if n % 1_000_000 == 0:
                print(f"  {n / 1e6:.0f} million rows ...")
            if row.get("confidence", "").lower().startswith("l"):
                continue                       # drop low-confidence detections
            if row.get("type", "0") not in ("0", ""):
                continue                       # keep vegetation fires only
            d = row["acq_date"]
            if nrt and d <= last_archive:
                continue                       # prefer standard data where both exist
            if not nrt and d > last_archive:
                last_archive = d
            ix = int(float(row["longitude"]) * 10 // 1)
            iy = int(float(row["latitude"]) * 10 // 1)
            counts[(d, row["daynight"], ix, iy)] += 1

out = os.path.join(here, "india_snpp_counts.csv.gz")
with gzip.open(out, "wt", newline="") as g:
    w = csv.writer(g)
    w.writerow(["date", "daynight", "ix", "iy", "n"])
    for (d, dn, ix, iy), c in sorted(counts.items()):
        w.writerow([d, dn, ix, iy, c])
print(f"Done: {n:,} rows read, {sum(counts.values()):,} kept, {len(counts):,} lines written")
print("Saved", out, f"({os.path.getsize(out) / 1e6:.1f} MB). Tell Claude it is ready.")
