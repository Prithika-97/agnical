"""Where does a fixed MODIS(Aqua) -> Suomi NPP mapping stay valid after Aqua's orbit drift?

Fit each state-season's Suomi NPP / Aqua day-fire ratio on 2012-2022 (before the drift), apply it to
2023-2025 Aqua counts and compare with what Suomi NPP counted. Reported for Punjab + Haryana and for
all other testable states. Also tests an hour-aware version (ratio adjusted for Aqua's viewing hour,
slope fitted outside the northwest on 2012-2022); it did not improve the result and is recorded as such.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
d = pd.read_csv(ROOT / "data/processed/india_state_indicators.csv")
d = d[d.testable & (d.aqua_day > 0) & (d.snpp_day > 0)].copy()
NW = ["Punjab", "Haryana"]; key = ["state", "season"]
tr = d[(d.year >= 2012) & (d.year <= 2022)]
base = tr.groupby(key).apply(lambda g: g.snpp_day.sum() / g.aqua_day.sum(), include_groups=False).rename("r0").reset_index()
te = d[(d.year >= 2023) & (d.year <= 2025)].merge(base, on=key)
te["pred"] = te.aqua_day * te.r0
te["grp"] = np.where(te.state.isin(NW), "Punjab and Haryana", "Other states")
rows = []
for (s, g, y), x in te.groupby(["season", "grp", "year"]):
    rows.append({"season": s, "group": g, "year": int(y), "n_states": int(len(x)), "snpp_observed": int(x.snpp_day.sum()),
                 "predicted_from_aqua": int(round(x.pred.sum())), "error_pct": round((x.pred.sum() / x.snpp_day.sum() - 1) * 100, 1),
                 "median_state_error_pct": round(float(((x.pred - x.snpp_day).abs() / x.snpp_day).median() * 100), 1)})
# hour-aware variant (recorded negative result)
h = d[d.aqua_time_ist.notna()].copy(); h["lr"] = np.log(h.snpp_day / h.aqua_day)
t2 = h[(h.year >= 2012) & (h.year <= 2022) & ~h.state.isin(NW)].copy()
t2["lr_m"] = t2.lr - t2.groupby(key).lr.transform("mean"); t2["h_m"] = t2.aqua_time_ist - t2.groupby(key).aqua_time_ist.transform("mean")
w = t2.snpp_day * t2.aqua_day / (t2.snpp_day + t2.aqua_day)
b = float((w * t2.h_m * t2.lr_m).sum() / (w * t2.h_m ** 2).sum())
bs = t2.groupby(key).agg(mlr=("lr", "mean"), mh=("aqua_time_ist", "mean")).reset_index()
t3 = h[(h.year >= 2023) & (h.year <= 2025) & ~h.state.isin(NW)].merge(bs, on=key)
pa = t3.aqua_day * np.exp(t3.mlr); pb = t3.aqua_day * np.exp(t3.mlr + b * (t3.aqua_time_ist - t3.mh))
ob = t3.snpp_day.sum()
out = {"rows": rows, "hour_aware_outside_nw": {"slope_per_hour": round(b, 3), "stationary_error_pct": round((pa.sum() / ob - 1) * 100, 1),
       "hour_aware_error_pct": round((pb.sum() / ob - 1) * 100, 1), "note": "The hour adjustment did not reduce the error, so it is not used."}}
(ROOT / "docs/data/harmonize_india.json").write_text(json.dumps(out, indent=1))
pd.DataFrame(rows).to_csv(ROOT / "data/processed/harmonize_india_drift.csv", index=False)
print(json.dumps(out["hour_aware_outside_nw"]))
