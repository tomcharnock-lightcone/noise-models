"""The free slopes against 1, and what holding them at 1 does to σ_w and f_knee.

usage: python src/plot_slope_comparison.py FREE.csv FIXED.csv OUT.png

Top: each sensor's free slope with its ±1σ uncertainty, against α = 1.
Below: the held fit's σ_w and f_knee relative to the free fit's, with the
held fit's ±1σ as error bars and the free fit's ±1σ shaded behind.
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INK, INK_2, MUTED, GRID = "#0b0b0b", "#52514e", "#a8a7a1", "#e4e3df"
POINT, BAND, SURFACE = "#2a78d6", "#d9d8d3", "#fcfcfb"


def main(free_path: str, fixed_path: str, out: str) -> None:
    free = pd.read_csv(free_path).set_index("sensor")
    fixed = pd.read_csv(fixed_path).set_index("sensor")
    x = np.arange(len(free))

    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK})
    fig, axes = plt.subplots(3, 1, figsize=(11, 7.5), sharex=True, facecolor=SURFACE,
                             constrained_layout=True)
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, axis="y", color=GRID, lw=0.6)

    ax = axes[0]
    ax.axhline(1, color=INK_2, lw=0.8)
    ax.errorbar(x, free.alpha, free.alpha_err, fmt="o", ms=4, color=POINT, mec=SURFACE,
                mew=1, elinewidth=1.2, capsize=0)
    pull = (free.alpha - 1) / free.alpha_err
    for i, sensor in enumerate(free.index):
        if abs(pull[sensor]) > 2:
            ax.annotate(f"{pull[sensor]:+.1f}σ", (i, free.alpha[sensor] + free.alpha_err[sensor]),
                        xytext=(0, 3), textcoords="offset points", ha="center", va="bottom",
                        fontsize=7, color=INK_2)
    ax.set_title("free fit's slope α (±1σ), against 1", loc="left", color=INK)

    for ax, title, ratio, err_fixed, err_free in [
        (axes[1], "f_knee: held at α = 1 / free − 1",
         fixed.f_knee_hz / free.f_knee_hz - 1,
         fixed.f_knee_hz_err / free.f_knee_hz, free.f_knee_hz_err / free.f_knee_hz),
        (axes[2], "σ_w: held at α = 1 / free − 1",
         fixed.sigma_white_c / free.sigma_white_c - 1,
         fixed.sigma_white_c_err / free.sigma_white_c, free.sigma_white_c_err / free.sigma_white_c),
    ]:
        ax.fill_between(x, -err_free, err_free, step="mid", color=BAND, lw=0)
        ax.axhline(0, color=INK_2, lw=0.8)
        ax.errorbar(x, ratio, err_fixed, fmt="o", ms=4, color=POINT, mec=SURFACE, mew=1,
                    elinewidth=1.2, capsize=0)
        ax.set_title(title, loc="left", color=INK)

    axes[-1].set_xticks(x, free.index, rotation=90)
    axes[-1].set_xlabel("sensor")
    handles = [
        plt.Line2D([], [], color=POINT, marker="o", lw=1.2, label="value ±1σ (held fit, lower panels)"),
        plt.Rectangle((0, 0), 1, 1, color=BAND, label="free fit ±1σ"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    fig.suptitle(f"Holding the 1/f slope at 1 ({free.window.iloc[0]} window)", fontsize=11, color=INK)
    fig.savefig(out, dpi=110, facecolor=SURFACE, metadata={"Software": None})


if __name__ == "__main__":
    main(*sys.argv[1:])
