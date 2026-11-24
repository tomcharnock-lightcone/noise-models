"""Fit the white + 1/f noise model to each sensor's Welch PSD.

usage: python src/fit_noise.py PSD.npz SLOPE OUT.csv

SLOPE is an option of the `slope` decision (noise_model.SLOPES): `free`
fits α with σ_w and f_knee, `one` holds α at exactly 1 (alpha_err is then 0).

One row per sensor. Uncertainties are the spread of a parametric bootstrap:
N_BOOT series are simulated from each sensor's fit, put through the same
Welch estimate and fit, and the standard deviation of their parameters taken
(in ln σ_w and ln f_knee, propagated to σ_w and f_knee; directly for α).

low_band_ratio is the mean of data / model over the lowest decade of the fit
band, and low_band_ratio_z how far it sits from the bootstrap's in units of
the bootstrap's spread: a sensor whose spectrum bends away from 1/f at low
frequency shows |z| well above 2.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

import noise_model

N_BOOT = 200
#: Fixed so a rebuild is byte-identical; each sensor draws from its own stream.
SEED = 20221212


def main(psd_path: str, slope: str, out: str) -> None:
    alpha = noise_model.SLOPES[slope]
    d = np.load(psd_path)
    f, fs = d["freq_hz"], float(d["fs_hz"])
    window, segment_s, n_samples = str(d["window"]), float(d["segment_s"]), int(d["n_samples"])
    rows = []
    for i, sensor in enumerate(d["sensors"]):
        fit = noise_model.fit_psd(f, d["psd_c2_per_hz"][i], fs, alpha)
        rng = np.random.default_rng([SEED, i])
        boot = noise_model.bootstrap(fit, n_samples, fs, window, segment_s, N_BOOT, rng)
        spread = boot.std(axis=0, ddof=1)
        rows.append({
            "sensor": str(sensor),
            "mean_c": d["mean_c"][i],
            "std_c": d["std_c"][i],
            "sigma_white_c": fit["sigma_white"],
            # σ_w = exp(ln σ_w² / 2)
            "sigma_white_c_err": fit["sigma_white"] * spread[0] / 2,
            "white_psd_c2_per_hz": 2 * fit["sigma_white"] ** 2 / fs,
            "f_knee_hz": fit["f_knee"],
            "f_knee_hz_err": fit["f_knee"] * spread[1],
            "alpha": fit["alpha"],
            "alpha_err": spread[2],
            "alpha_fixed": fit["alpha_fixed"],
            "knee_above_nyquist": fit["f_knee"] > fs / 2,
            "low_band_ratio": fit["low_band_ratio"],
            "low_band_ratio_z": (fit["low_band_ratio"] - boot[:, 3].mean()) / spread[3],
            "at_bound": fit["at_bound"],
            "window": window,
            "segment_s": segment_s,
            "fs_hz": fs,
            "fit_f_min_hz": fit["fit_f_min"],
            "fit_f_max_hz": fit["fit_f_max"],
            "n_fit_bins": fit["n_fit_bins"],
            "n_boot": N_BOOT,
        })
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.6g")


if __name__ == "__main__":
    main(*sys.argv[1:])
