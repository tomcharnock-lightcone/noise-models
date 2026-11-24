"""Do the sensors' 1/f slopes agree with exactly 1, and what does holding them there change?

usage: python src/slope_test.py FREE.csv FIXED.csv OUT.json

FREE and FIXED are noise_models tables fitted with the slope free and held at
1, under the same window. Reported:

- population: χ² of (α − 1)/σ_α over the sensors and its p-value; the
  precision-weighted mean slope; and the intrinsic scatter τ of the slopes
  about their common mean (α_i ~ N(μ, σ_i² + τ²), μ profiled out), as a
  maximum-likelihood value and a one-sided 95% upper limit. τ bounds how far
  holding every sensor at one slope can be wrong.
- tension: every sensor whose free slope is more than 2σ from 1, how many
  that many sensors give by chance, and how far holding its slope at 1 moves
  its σ_w and f_knee, in units of the free fit's uncertainty.
- precision: the median ratio of the held fit's uncertainty to the free
  fit's, for σ_w and f_knee.
"""

from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
from scipy import optimize, stats


def _neg2lnl(tau: float, a: np.ndarray, s: np.ndarray) -> float:
    """-2 ln L of α_i ~ N(μ, s_i² + τ²), with μ at its maximum for this τ."""
    v = s**2 + tau**2
    mu = np.sum(a / v) / np.sum(1 / v)
    return float(np.sum(np.log(v) + (a - mu) ** 2 / v))


def intrinsic_scatter(a: np.ndarray, s: np.ndarray) -> tuple[float, float]:
    """ML intrinsic scatter and its one-sided 95% upper limit (Δ(−2 ln L) = 2.71)."""
    top = 10 * float(np.std(a)) + float(np.max(s))
    res = optimize.minimize_scalar(_neg2lnl, bounds=(0.0, top), args=(a, s), method="bounded",
                                   options={"xatol": 1e-7})
    tau_ml = float(res.x) if _neg2lnl(res.x, a, s) < _neg2lnl(0.0, a, s) else 0.0
    best = _neg2lnl(tau_ml, a, s)
    upper = optimize.brentq(lambda t: _neg2lnl(t, a, s) - best - 2.71, tau_ml, top)
    return tau_ml, float(upper)


def main(free_path: str, fixed_path: str, out: str) -> None:
    free = pd.read_csv(free_path).set_index("sensor")
    fixed = pd.read_csv(fixed_path).set_index("sensor")
    if free.alpha_fixed.any() or not fixed.alpha_fixed.all():
        raise ValueError("expected the free-slope table first, then the held-slope table")
    if not free.index.equals(fixed.index) or (free.window != fixed.window).any():
        raise ValueError("the two tables do not cover the same sensors under the same window")

    a, s = free.alpha.to_numpy(), free.alpha_err.to_numpy()
    pull = (free.alpha - 1) / free.alpha_err
    chi2 = float(np.sum(pull**2))
    w = 1 / s**2
    tau_ml, tau_95 = intrinsic_scatter(a, s)

    sig_shift = (fixed.sigma_white_c - free.sigma_white_c) / free.sigma_white_c_err
    fk_shift = np.log(fixed.f_knee_hz / free.f_knee_hz) / (free.f_knee_hz_err / free.f_knee_hz)
    tension = [
        {
            "sensor": sensor,
            "alpha": round(float(free.alpha[sensor]), 4),
            "alpha_err": round(float(free.alpha_err[sensor]), 4),
            "alpha_minus_one_in_errors": round(float(pull[sensor]), 2),
            "sigma_white_shift_in_errors": round(float(sig_shift[sensor]), 2),
            "f_knee_shift_in_errors": round(float(fk_shift[sensor]), 2),
        }
        for sensor in pull.abs().sort_values(ascending=False).index
        if abs(pull[sensor]) > 2
    ]
    report = {
        "window": str(free.window.iloc[0]),
        "n_sensors": len(free),
        "population": {
            "chi2_alpha_one": round(chi2, 1),
            "dof": len(free),
            "p_value": round(float(stats.chi2.sf(chi2, len(free))), 3),
            "weighted_mean_alpha": round(float(np.sum(w * a) / np.sum(w)), 4),
            "weighted_mean_alpha_err": round(float(np.sum(w) ** -0.5), 4),
            "intrinsic_scatter_ml": round(tau_ml, 4),
            "intrinsic_scatter_95_upper": round(tau_95, 4),
        },
        "tension": {
            "beyond_2_sigma_expected_by_chance": round(len(free) * 2 * float(stats.norm.sf(2)), 1),
            "beyond_3_sigma_expected_by_chance": round(len(free) * 2 * float(stats.norm.sf(3)), 2),
            "sensors_beyond_2_sigma": tension,
        },
        "precision": {
            "median_err_ratio_fixed_over_free_sigma_white": round(float(
                (fixed.sigma_white_c_err / free.sigma_white_c_err).median()), 2),
            "median_err_ratio_fixed_over_free_f_knee": round(float(
                ((fixed.f_knee_hz_err / fixed.f_knee_hz) / (free.f_knee_hz_err / free.f_knee_hz)).median()), 2),
            "median_abs_sigma_white_shift_in_errors": round(float(sig_shift.abs().median()), 2),
            "median_abs_f_knee_shift_in_errors": round(float(fk_shift.abs().median()), 2),
        },
    }
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
