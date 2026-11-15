"""MODIS -> VIIRS-equivalent detections, fit on the overlap years."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import mean_poisson_deviance

from . import config
from .grid import PERIOD_DAYS

FEATURES = ["m_count", "m_terra_day", "m_terra_night", "m_aqua_day", "m_aqua_night",
            "m_frp_sum", "m_frp_max", "m_nbr", "m_prev", "m_next", "period", "lat", "lon"]
TARGET = "v_count"


def overlap(panel: pd.DataFrame, years=config.TRAIN_YEARS) -> pd.DataFrame:
    """Rows where MODIS and SNPP both observed; default 2012-2022, before the Aqua drift."""
    ov = panel[panel["m_valid"] & panel["v_valid"]]
    return ov if years is None else ov[ov["year"].isin(list(years))]


class GlobalRatio:
    """Baseline 1: one regional scale factor."""
    name = "global_ratio"

    def fit(self, df):
        self.r = df[TARGET].sum() / max(df["m_count"].sum(), 1e-9)
        return self

    def predict(self, df):
        return self.r * df["m_count"].to_numpy()


class CellRatio:
    """Baseline 2: per-cell factor, shrunk toward the regional one."""
    name = "cell_ratio"

    def __init__(self, k=20.0):
        self.k = k

    def fit(self, df):
        self.r0 = df[TARGET].sum() / max(df["m_count"].sum(), 1e-9)
        g = df.groupby(["ix", "iy"])[[TARGET, "m_count"]].sum()
        denom = (g["m_count"] + self.k).replace(0, np.nan)
        self.r = ((g[TARGET] + self.r0 * self.k) / denom).rename("r")
        return self

    def predict(self, df):
        r = df[["ix", "iy"]].merge(self.r, left_on=["ix", "iy"], right_index=True, how="left")["r"]
        return r.fillna(self.r0).to_numpy() * df["m_count"].to_numpy()


class GBMHarmonizer:
    """LightGBM, Poisson objective."""
    name = "gbm_poisson"

    def __init__(self, n_estimators=400, fast=False):
        import lightgbm as lgb
        self.model = lgb.LGBMRegressor(
            objective="poisson", n_estimators=120 if fast else n_estimators,
            learning_rate=0.08 if fast else 0.05, num_leaves=63, min_child_samples=50,
            subsample=0.8, subsample_freq=1, colsample_bytree=0.9,
            random_state=42, verbose=-1)

    def fit(self, df):
        self.model.fit(df[FEATURES], df[TARGET])
        return self

    def predict(self, df):
        return np.clip(self.model.predict(df[FEATURES]), 0, None)

    def importance(self):
        gain = self.model.booster_.feature_importance(importance_type="gain")
        s = pd.Series(gain, index=FEATURES)
        return (s / s.sum()).sort_values(ascending=False)


def make_models(fast=False):
    return [GlobalRatio(), CellRatio(), GBMHarmonizer(fast=fast)]


def cross_validate(panel, fast=False, years=None, log=print) -> pd.DataFrame:
    """Leave-one-year-out CV; returns held-out rows with predictions."""
    ov = overlap(panel)
    years = sorted(ov["year"].unique()) if years is None else years
    out = []
    for y in years:
        train, test = ov[ov["year"] != y], ov[ov["year"] == y].copy()
        for model in make_models(fast):
            test[f"pred_{model.name}"] = model.fit(train).predict(test)
        out.append(test)
        log(f"  CV fold {y}: done ({len(test):,} rows)")
    return pd.concat(out, ignore_index=True)


def season_periods(season: str) -> list[int]:
    lo, hi = config.SEASONS[season]
    return [p for p in range(config.N_PERIODS) if lo <= p * PERIOD_DAYS + 1 <= hi]


def evaluate(cv: pd.DataFrame) -> dict:
    """Metrics: cell, region, season total, peak timing."""
    names = [c.removeprefix("pred_") for c in cv.columns if c.startswith("pred_")]
    y = cv[TARGET].to_numpy()
    res = {"n_rows": int(len(cv)), "years": sorted(int(v) for v in cv["year"].unique()),
           "models": {}}
    reg = cv.groupby(["year", "period"])[[TARGET] + [f"pred_{n}" for n in names]].sum()
    for n in names:
        p = cv[f"pred_{n}"].to_numpy()
        r = reg[f"pred_{n}"]
        m = {
            "cell_rmse": float(np.sqrt(np.mean((p - y) ** 2))),
            "cell_mae": float(np.mean(np.abs(p - y))),
            "cell_poisson_deviance": float(mean_poisson_deviance(y, np.clip(p, 1e-6, None))),
            "regional_corr": float(np.corrcoef(r, reg[TARGET])[0, 1]),
            "regional_rmse": float(np.sqrt(np.mean((r - reg[TARGET]) ** 2))),
            "seasons": {},
        }
        for s in config.SEASONS:
            sp = season_periods(s)
            sub = reg[reg.index.get_level_values("period").isin(sp)]
            tot = sub.groupby("year")[[TARGET, f"pred_{n}"]].sum()
            tot = tot[tot[TARGET] > 0]
            rel = (tot[f"pred_{n}"] - tot[TARGET]) / tot[TARGET]
            pk_obs = sub[TARGET].groupby("year").idxmax().map(lambda t: t[1])
            pk_pred = sub[f"pred_{n}"].groupby("year").idxmax().map(lambda t: t[1])
            peak_err = (pk_pred - pk_obs).abs() * PERIOD_DAYS
            m["seasons"][s] = {
                "total_mape_pct": float(rel.abs().mean() * 100),
                "total_bias_pct": float(rel.mean() * 100),
                "total_rel_err_q10": float(rel.quantile(0.1)),
                "total_rel_err_q90": float(rel.quantile(0.9)),
                "peak_err_days_mean": float(peak_err.mean()),
                "peak_exact_pct": float((peak_err == 0).mean() * 100),
            }
        res["models"][n] = m
    # Is the difference in season-total error between the model and the regional ratio more than
    # chance? Paired Wilcoxon test on the absolute yearly errors.
    if {"gbm_poisson", "global_ratio"} <= set(names):
        from scipy.stats import wilcoxon
        res["total_error_test"] = {}
        for s in config.SEASONS:
            sub = reg[reg.index.get_level_values("period").isin(season_periods(s))]
            tot = sub.groupby("year").sum()
            tot = tot[tot[TARGET] > 0]
            a = ((tot["pred_gbm_poisson"] - tot[TARGET]) / tot[TARGET]).abs()
            b = ((tot["pred_global_ratio"] - tot[TARGET]) / tot[TARGET]).abs()
            try:
                res["total_error_test"][s] = float(wilcoxon(a, b).pvalue)
            except ValueError:
                res["total_error_test"][s] = float("nan")
    return res


def ratio_stability(panel: pd.DataFrame) -> pd.DataFrame:
    """Observed VIIRS/MODIS ratio per overlap year."""
    ov = overlap(panel, years=None)
    g = ov.groupby("year")[[TARGET, "m_count"]].sum()
    g["ratio"] = g[TARGET] / g["m_count"].replace(0, np.nan)
    return g.reset_index()


def drift_test(panel: pd.DataFrame, fast=False) -> pd.DataFrame:
    """Train 2012-22, predict 2023 on. Season totals, observed vs predicted."""
    ov_all = overlap(panel, years=None)
    train = overlap(panel)
    test = ov_all[~ov_all["year"].isin(list(config.TRAIN_YEARS))].copy()
    if test.empty:
        return pd.DataFrame()
    test["pred"] = GBMHarmonizer(fast=fast).fit(train).predict(test)
    rows = []
    for s in config.SEASONS:
        sp = season_periods(s)
        g = test[test["period"].isin(sp)].groupby("year")[[TARGET, "pred"]].sum()
        for y, r in g.iterrows():
            if r[TARGET] < 100:   # season not (yet) covered by Suomi NPP
                continue
            rows.append({"year": int(y), "season": s, "observed": float(r[TARGET]), "predicted": float(r["pred"]),
                         "error_pct": float((r["pred"] - r[TARGET]) / max(r[TARGET], 1) * 100)})
    return pd.DataFrame(rows)


def build_harmonized(panel: pd.DataFrame, fast=False):
    """Fit on all overlap years; VIIRS-equivalent record for every year."""
    ov = overlap(panel)
    model = GBMHarmonizer(fast=fast).fit(ov)
    h = panel[panel["m_valid"] | panel["v_valid"]].copy()
    h["model_pred"] = np.where(h["m_valid"], model.predict(h), np.nan)
    h["viirs_equiv"] = np.where(h["v_valid"], h[TARGET], h["model_pred"])
    h["source"] = np.where(h["v_valid"], "observed", "reconstructed")
    h["stitched_raw"] = np.where(h["v_valid"], h[TARGET], h["m_count"])
    return h, model
