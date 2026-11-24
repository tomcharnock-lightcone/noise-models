"""Welch PSD of every sensor's chamber readings.

usage: python src/welch_psd.py READINGS.csv WINDOW SEGMENT_LENGTH OUT.npz

WINDOW is a scipy window name (hann, blackman); SEGMENT_LENGTH an option of
the `segment_length` decision (noise_model.SEGMENT_SECONDS).
"""

from __future__ import annotations

import sys

import numpy as np

import noise_model
import readings


def main(csv_path: str, window: str, segment_length: str, out: str) -> None:
    wide, fs = readings.load(csv_path)
    segment_s = noise_model.SEGMENT_SECONDS[segment_length]
    x = wide.to_numpy()
    f, p, info = noise_model.welch(x, fs, window, segment_s)
    with open(out, "wb") as fh:
        np.savez(
            fh,
            sensors=np.array(wide.columns, dtype=str),
            freq_hz=f,
            psd_c2_per_hz=p,
            fs_hz=fs,
            window=window,
            segment_s=float(segment_s),
            n_samples=x.shape[0],
            mean_c=x.mean(axis=0),
            std_c=x.std(axis=0, ddof=1),
            start=str(wide.index[0].isoformat()),
            end=str(wide.index[-1].isoformat()),
            **info,
        )


if __name__ == "__main__":
    main(*sys.argv[1:])
