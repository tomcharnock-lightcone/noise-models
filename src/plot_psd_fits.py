"""Every sensor's Welch PSD with its fitted noise model, one panel per sensor.

usage: python src/plot_psd_fits.py PSD.npz NOISE_MODELS.csv OUT.png
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import noise_model

INK, INK_2, MUTED, GRID = "#0b0b0b", "#52514e", "#a8a7a1", "#e4e3df"
MODEL, SURFACE = "#2a78d6", "#fcfcfb"


def main(psd_path: str, models_path: str, out: str) -> None:
    d = np.load(psd_path)
    models = pd.read_csv(models_path).set_index("sensor")
    f, fs = d["freq_hz"], float(d["fs_hz"])
    band = noise_model.fit_band(f, fs)
    sensors = [str(s) for s in d["sensors"]]
    ncol = 5
    nrow = int(np.ceil(len(sensors) / ncol))

    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK})
    fig, axes = plt.subplots(nrow, ncol, figsize=(13, 2.1 * nrow), sharex=True,
                             facecolor=SURFACE, constrained_layout=True)
    fm = np.geomspace(f[band][0], f[band][-1], 200)
    for ax, sensor, p in zip(axes.flat, sensors, d["psd_c2_per_hz"]):
        m = models.loc[sensor]
        ax.set_facecolor(SURFACE)
        ax.grid(True, which="major", color=GRID, lw=0.6)
        ax.loglog(f[1:], p[1:], color=MUTED, lw=0.5)
        ax.loglog(fm, noise_model.psd_model(fm, m.sigma_white_c, m.f_knee_hz, m.alpha, fs),
                  color=MODEL, lw=2)
        ax.axhline(m.white_psd_c2_per_hz, color=INK_2, lw=1, ls="--")
        if m.f_knee_hz < fs / 2:
            ax.axvline(m.f_knee_hz, color=INK_2, lw=1, ls=":")
        ax.set_title(f"{sensor}   σ_w = {m.sigma_white_c:.3f} °C   α = {m.alpha:.2f}",
                     loc="left", fontsize=8, color=INK)
    for ax in axes.flat[len(sensors):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("frequency (Hz)")
    for ax in axes[:, 0]:
        ax.set_ylabel("PSD (°C²/Hz)")

    handles = [
        plt.Line2D([], [], color=MUTED, lw=1, label="Welch estimate"),
        plt.Line2D([], [], color=MODEL, lw=2, label="white + 1/f fit (over the fitted bins)"),
        plt.Line2D([], [], color=INK_2, lw=1, ls="--", label="white floor 2σ_w²/fs"),
        plt.Line2D([], [], color=INK_2, lw=1, ls=":", label="knee (when below Nyquist)"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle(f"Chamber noise spectra, {d['window']} window, "
                 f"{float(d['segment_s']) / 3600:g} h segments", fontsize=11, color=INK)
    fig.savefig(out, dpi=110, facecolor=SURFACE, metadata={"Software": None})


if __name__ == "__main__":
    main(*sys.argv[1:])
