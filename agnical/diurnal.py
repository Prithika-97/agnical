"""Time of day of crop fires from polar-orbiting satellites only.

Each pass is a snapshot at a known local time. On one day the ratio of two passes depends on
relative sensor sensitivity and the shape of the daily burning curve; the day's total cancels.
Terra and Aqua carry the same MODIS; VIIRS (SNPP, NOAA-20/21) passes near 13:30. Day-by-day
comparison gives the peak hour per season since 2003.

Model
-----
For pass p on day d, with local time t_p and sensor class k(p):
    n_dp ~ Poisson( exp(alpha_d + beta_k(p) + q_y(t_p)) )
    q_y(t) = -(t - mu_y)^2 / (2 sigma_y^2)
alpha_d (daily total) drops out by conditioning on each day's total (multinomial likelihood).
beta_MODIS = 0 for Terra and Aqua (same instrument); one beta per VIIRS satellite. mu_y: peak hour,
sigma_y: window width, with a weak hierarchical prior on sigma. Uncertainty: day-level bootstrap.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize

IST_OFFSET_H = 5.5
SENSOR_CLASSES = ["VIIRS_SNPP", "VIIRS_NOAA20", "VIIRS_NOAA21"]  # MODIS is the reference


def to_ist_hours(hhmm: pd.Series) -> pd.Series:
    return ((hhmm // 100) + (hhmm % 100) / 60 + IST_OFFSET_H) % 24


def sensor_class(det: pd.DataFrame) -> pd.Series:
    return np.where(det["sensor"] == "MODIS", "MODIS", "VIIRS_" + det["satellite"].astype(str))


def build_passes(det: pd.DataFrame, doy_range=(274, 334), gap_min=30) -> pd.DataFrame:
    """One row per daytime pass: day, sensor class, local time, count."""
    d = det[(det["daynight"] == "D")].copy()
    d["doy"] = d["date"].dt.dayofyear
    d = d[d["doy"].between(*doy_range)]
    d["cls"] = sensor_class(d)
    d["sat"] = d["satellite"].astype(str)
    d["tmin"] = (d["time"] // 100) * 60 + d["time"] % 100
    d = d.sort_values(["date", "sat", "tmin"])
    new = (d.groupby(["date", "sat"])["tmin"].diff().fillna(1e9) > gap_min).astype(int)
    d["pass_id"] = new.groupby([d["date"], d["sat"]]).cumsum()
    p = d.groupby(["date", "sat", "pass_id"]).agg(
        cls=("cls", "first"), n=("tmin", "size"), tmin=("tmin", "median")).reset_index()
    p["t"] = ((p["tmin"] / 60) + IST_OFFSET_H) % 24
    p["year"] = p["date"].dt.year
    return p[["date", "year", "sat", "cls", "t", "n"]]


def _unpack(theta, years, classes):
    ny = len(years)
    mu = dict(zip(years, theta[:ny]))
    log_sig = dict(zip(years, theta[ny:2 * ny]))
    log_sig0 = theta[2 * ny]
    beta = {"MODIS": 0.0, **dict(zip(classes, theta[2 * ny + 1:]))}
    return mu, {y: float(np.exp(v)) for y, v in log_sig.items()}, beta, float(log_sig0)


def _arrays(P, years, classes):
    day_codes, day_idx = np.unique(P["date"].astype(str).to_numpy(), return_inverse=True)
    yi = {y: i for i, y in enumerate(years)}
    ci = {c: i for i, c in enumerate(classes)}
    return {
        "day": day_idx, "nd": len(day_codes),
        "n": P["n"].to_numpy(float), "t": P["t"].to_numpy(float),
        "y": P["year"].map(yi).to_numpy(), "c": P["cls"].map(lambda c: ci.get(c, -1)).to_numpy(),
        "ny": len(years), "nc": len(classes),
    }


PRIOR_SD_LOG_SIGMA = 1.0   # weak: lets each season's window width follow its own data


def _nll_grad(theta, A):
    ny, nc = A["ny"], A["nc"]
    mu, ls, ls0, beta = theta[:ny], theta[ny:2 * ny], theta[2 * ny], theta[2 * ny + 1:]
    s2 = np.exp(2 * ls)
    dt = A["t"] - mu[A["y"]]
    b = np.where(A["c"] >= 0, beta[np.clip(A["c"], 0, None)] if nc else 0.0, 0.0)
    eta = b - dt ** 2 / (2 * s2[A["y"]])
    # per-day log-sum-exp (stable)
    mx = np.full(A["nd"], -np.inf); np.maximum.at(mx, A["day"], eta)
    ex = np.exp(eta - mx[A["day"]])
    se = np.bincount(A["day"], ex, A["nd"])
    N = np.bincount(A["day"], A["n"], A["nd"])
    lse = np.log(se) + mx
    ll = (A["n"] * eta).sum() - (N * lse).sum()
    # gradient of ll wrt eta: n_p - N_d * softmax_p
    g_eta = A["n"] - N[A["day"]] * ex / se[A["day"]]
    g_mu = np.bincount(A["y"], g_eta * dt / s2[A["y"]], ny)
    g_ls = np.bincount(A["y"], g_eta * dt ** 2 / s2[A["y"]], ny)
    g_beta = np.bincount(np.clip(A["c"], 0, None), g_eta * (A["c"] >= 0), nc) if nc else np.zeros(0)
    # priors: mu_y ~ N(14, 4^2); log sigma_y ~ N(log sigma0, sd^2); log sigma0 ~ N(log 1.5, 1)
    sd = A.get("prior_sd", PRIOR_SD_LOG_SIGMA)
    pr = ((mu - 14) / 4) ** 2 / 2 + ((ls - ls0) / sd) ** 2 / 2
    prior = pr.sum() + (ls0 - np.log(1.5)) ** 2 / 2
    g_pr_mu = (mu - 14) / 16
    g_pr_ls = (ls - ls0) / sd ** 2
    g_pr_ls0 = -((ls - ls0) / sd ** 2).sum() + (ls0 - np.log(1.5))
    f = -ll + prior
    g = np.r_[-g_mu + g_pr_mu, -g_ls + g_pr_ls, g_pr_ls0, -g_beta]
    return f, g


def fit(passes: pd.DataFrame, classes=None, x0=None, prior_sd=PRIOR_SD_LOG_SIGMA):
    P = passes[passes.groupby("date")["n"].transform("size") > 1]   # only days with >1 pass carry information
    years = sorted(P["year"].unique())
    classes = [c for c in (classes or SENSOR_CLASSES) if c in set(P["cls"])]
    A = _arrays(P, years, classes)
    A["prior_sd"] = prior_sd
    theta0 = x0 if x0 is not None else np.r_[np.full(len(years), 14.0), np.full(len(years), np.log(1.5)),
                                             np.log(1.5), np.full(len(classes), 1.5)]
    res = minimize(_nll_grad, theta0, args=(A,), jac=True, method="L-BFGS-B",
                   options={"maxiter": 2000})
    mu, sigma, beta, log_sig0 = _unpack(res.x, years, classes)
    return {"mu": mu, "sigma": sigma, "sigma0": float(np.exp(log_sig0)), "beta": beta,
            "years": years, "classes": classes,
            "theta": res.x, "success": bool(res.success), "n_days": int(A["nd"]),
            "n_passes": int(len(P)), "prior_sd": prior_sd}


def bootstrap(passes: pd.DataFrame, base_fit: dict, n_boot=100, seed=0, log=None):
    """Resample days within year, refit; draws of mu and sigma."""
    rng = np.random.default_rng(seed)
    days = passes[["date", "year"]].drop_duplicates().reset_index(drop=True)
    by_year = {y: g["date"].to_numpy() for y, g in days.groupby("year")}
    draws = []
    for b in range(n_boot):
        picks = np.concatenate([rng.choice(v, size=len(v), replace=True) for v in by_year.values()])
        pick = pd.DataFrame({"date": picks, "rep": np.arange(len(picks))})
        Pb = passes.merge(pick, on="date")
        Pb["date"] = Pb["date"].astype(str) + "_" + Pb["rep"].astype(str)   # duplicates become separate days
        f = fit(Pb, base_fit["classes"], x0=base_fit["theta"], prior_sd=base_fit.get("prior_sd", PRIOR_SD_LOG_SIGMA))
        draws.append({**{f"mu_{y}": v for y, v in f["mu"].items()},
                      **{f"sig_{y}": v for y, v in f["sigma"].items()}})
        if log and (b + 1) % 20 == 0:
            log(f"  bootstrap {b + 1}/{n_boot}")
    return pd.DataFrame(draws)


def summarize(base_fit: dict, draws: pd.DataFrame, passes: pd.DataFrame) -> pd.DataFrame:
    """Peak hour per season, 90% interval. Reported only if the interval is under 2 h and ends before
    the latest pass + 30 min; otherwise 'beyond view' and only the lower end is kept."""
    span = passes.groupby("year").agg(t_first=("t", "min"), t_last=("t", "max"), n_sat=("sat", "nunique"))
    rows = []
    for y in base_fit["years"]:
        mu_d = draws[f"mu_{y}"]
        lo, hi = mu_d.quantile(0.05), mu_d.quantile(0.95)
        t_last = span.loc[y, "t_last"] if y in span.index else np.nan
        n_sat = span.loc[y, "n_sat"] if y in span.index else 0
        if n_sat < 2:   # one satellite cannot separate time of day from day-to-day change
            continue
        ok = (hi - lo) < 2.0 and hi <= t_last + 0.5
        rows.append({
            "year": int(y),
            "peak_ist": base_fit["mu"][y] if ok else np.nan,
            "peak_lo": lo, "peak_hi": hi if ok else np.nan,
            "beyond_view": not ok,
            "width_h": base_fit["sigma"][y] if ok else np.nan,
            "latest_view_ist": t_last,
        })
    out = pd.DataFrame(rows)
    base = out[out["year"].between(2012, 2022)]["peak_ist"]
    out.attrs["baseline_peak_ist"] = float(base.mean())
    return out


FLAG_INDICATORS = [("aqua_over_snpp", 1), ("night_share", 1)]
# terra_over_aqua is reported but NOT used: Terra's own orbit drifted earlier after 2020, and
# the ratio fell just as much in the comparison region across the border (agnical/controls.py),
# so its fall is a satellite effect, not evidence about burning.
DIAGNOSTIC_INDICATORS = [("terra_over_aqua", -1)]


def timing_breaks(ind: pd.DataFrame, baseline=range(2012, 2023), k=3.0) -> pd.DataFrame:
    """Flag seasons outside the 2012-22 envelope. Per indicator: ratio to baseline mean, and distance
    from the baseline median in robust SDs (1.4826 * MAD). Break: both flag indicators outside; watch: one."""
    out = ind[["year"]].copy()
    flags = []
    for col, direction in FLAG_INDICATORS + DIAGNOSTIC_INDICATORS:
        b = ind[ind["year"].isin(list(baseline))][col].dropna()
        med, mad = b.median(), 1.4826 * (b - b.median()).abs().median()
        scale = mad if mad > 0 else (b.std(ddof=1) if b.std(ddof=1) > 0 else abs(med) * 0.05 or 1e-9)
        out[f"{col}_x_baseline"] = ind[col] / b.mean()
        z = (ind[col] - med) / scale
        out[f"{col}_z"] = z
        if (col, direction) in FLAG_INDICATORS:
            flags.append((direction * z > k).to_numpy())
    out["n_flags"] = np.sum(flags, axis=0)
    out["timing_break"] = out["n_flags"] >= len(FLAG_INDICATORS)
    out["timing_watch"] = out["n_flags"] == 1
    return out


def hhmm(h: float) -> str:
    if h is None or np.isnan(h):
        return "n/a"
    h = h % 24
    return f"{int(h):02d}:{int(round((h % 1) * 60)) % 60:02d}"


def indicators(det: pd.DataFrame, doy_range=(274, 334)) -> pd.DataFrame:
    """Per kharif season, what each satellite saw and when.

    aqua_over_snpp:  Aqua day fires / SNPP day fires.
    night_share:     SNPP night-pass (~01:30) fires / all SNPP fires.
    terra_over_aqua: Terra day fires / Aqua day fires (diagnostic only).
    """
    d = det.copy()
    d["doy"] = d["date"].dt.dayofyear
    d = d[d["doy"].between(*doy_range)]
    d["year"] = d["date"].dt.year
    d["t"] = to_ist_hours(d["time"])
    rows = []
    for y, g in d.groupby("year"):
        day = g[g["daynight"] == "D"]
        aq = day[(day["sensor"] == "MODIS") & (day["satellite"] == "Aqua")]
        te = day[(day["sensor"] == "MODIS") & (day["satellite"] == "Terra")]
        sn = g[(g["sensor"] == "VIIRS") & (g["satellite"] == "SNPP")]
        sn_day = sn[sn["daynight"] == "D"]
        rows.append({
            "year": int(y),
            "terra_day": int(len(te)), "aqua_day": int(len(aq)), "snpp_day": int(len(sn_day)),
            "snpp_night": int((sn["daynight"] != "D").sum()),
            "terra_time_ist": float(te["t"].median()) if len(te) else np.nan,
            "aqua_time_ist": float(aq["t"].median()) if len(aq) else np.nan,
            "snpp_time_ist": float(sn_day["t"].median()) if len(sn_day) else np.nan,
            "aqua_over_snpp": float(len(aq) / len(sn_day)) if len(sn_day) else np.nan,
            "night_share": float((sn["daynight"] != "D").sum() / len(sn)) if len(sn) else np.nan,
            "terra_over_aqua": float(len(te) / len(aq)) if len(aq) else np.nan,
        })
    return pd.DataFrame(rows)
