"""Blackman against Hann, parameter by parameter, for every sensor.

usage: python src/plot_window_comparison.py HANN.csv BLACKMAN.csv OUT.png

Each panel shows the difference (relative for σ_w and f_knee, absolute for
α) with the Hann fit's ±1σ bootstrap uncertainty shaded behind it.
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INK, INK_2, MUTED, GRID = "#0b0b0b", "#52514e", "#a8a7a1", "#e4e3df"
POINT, BAND, SURFACE = "#eb6834", "#d9d8d3", "#fcfcfb"


def main(hann_path: str, blackman_path: str, out: str) -> None:
    h = pd.read_csv(hann_path).set_index("sensor")
    b = pd.read_csv(blackman_path).set_index("sensor")
    x = np.arange(len(h))
    panels = [
        ("σ_w: Blackman / Hann − 1", (b.sigma_white_c / h.sigma_white_c - 1),
         h.sigma_white_c_err / h.sigma_white_c),
        ("f_knee: Blackman / Hann − 1", (b.f_knee_hz / h.f_knee_hz - 1), h.f_knee_hz_err / h.f_knee_hz),
        ("α: Blackman − Hann", (b.alpha - h.alpha), h.alpha_err),
    ]
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK})
    fig, axes = plt.subplots(3, 1, figsize=(11, 7), sharex=True, facecolor=SURFACE,
                             constrained_layout=True)
    for ax, (title, diff, err) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        ax.grid(True, axis="y", color=GRID, lw=0.6)
        ax.fill_between(x, -err, err, step="mid", color=BAND, lw=0)
        ax.axhline(0, color=INK_2, lw=0.8)
        ax.plot(x, diff, "o", ms=4, color=POINT, mec=SURFACE, mew=1)
        ax.set_title(title, loc="left", color=INK)
    axes[-1].set_xticks(x, h.index, rotation=90)
    axes[-1].set_xlabel("sensor")
    handles = [
        plt.Line2D([], [], color=POINT, marker="o", ls="", label="difference"),
        plt.Rectangle((0, 0), 1, 1, color=BAND, label="Hann fit ±1σ (bootstrap)"),
    ]
    fig.legend(handles=handles, loc="outside upper right", ncol=2, frameon=False)
    fig.suptitle("Does the window change the noise models?", fontsize=11, color=INK)
    fig.savefig(out, dpi=110, facecolor=SURFACE, metadata={"Software": None})


if __name__ == "__main__":
    main(*sys.argv[1:])
