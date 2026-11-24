"""Was the chamber's temperature constant enough that the sensors read only their own noise?

usage: python src/chamber_check.py READINGS.csv OUT.json

A drift in the chamber is common to every sensor; independent sensor noise
is not. The mean across N sensors of their mean-removed readings therefore
has a variance of Σσᵢ²/N² if the chamber held still, and more if it moved.
This compares the two, on the 10-second readings and on hourly means (where
a slow chamber drift would show first), and gives the sensors' pairwise
correlations, which a shared drift would make positive.

1/f noise makes these statistics scatter far more than for white noise, so
each comes with its null distribution: N_NULL copies of the data with every
sensor circularly shifted by its own random offset, which removes anything
the sensors share while keeping each one's own spectrum. The 2.5th and
97.5th percentiles of each statistic under that null are reported beside it.
"""

from __future__ import annotations

import json
import sys

import numpy as np

import readings

N_NULL = 200
SEED = 20221218


def _stats(x: np.ndarray) -> dict:
    x = x - x.mean(axis=0)
    n = x.shape[1]
    common = x.mean(axis=1).std(ddof=1)
    expected = np.sqrt(np.sum(x.var(axis=0, ddof=1))) / n
    r = np.corrcoef(x.T)[np.triu_indices(n, k=1)]
    return {
        "common_mode_std_c": float(common),
        "expected_if_independent_c": float(expected),
        "ratio": float(common / expected),
        "mean_pairwise_correlation": float(r.mean()),
        "max_abs_pairwise_correlation": float(np.abs(r).max()),
    }


def _hourly(x: np.ndarray, per_hour: int) -> np.ndarray:
    n = x.shape[0] // per_hour * per_hour
    return x[:n].reshape(-1, per_hour, x.shape[1]).mean(axis=1)


def _report(observed: dict, null: list[dict]) -> dict:
    out = {k: round(v, 5) for k, v in observed.items()}
    for k in ("ratio", "mean_pairwise_correlation", "max_abs_pairwise_correlation"):
        lo, hi = np.percentile([d[k] for d in null], [2.5, 97.5])
        out[f"{k}_null_95"] = [round(float(lo), 4), round(float(hi), 4)]
    return out


def main(csv_path: str, out: str) -> None:
    wide, fs = readings.load(csv_path)
    if wide.index[0] != wide.index[0].floor("1h"):
        raise ValueError("the readings do not start on the hour")
    per_hour = int(round(3600 * fs))
    x = wide.to_numpy()
    rng = np.random.default_rng(SEED)
    null_raw, null_hourly = [], []
    for _ in range(N_NULL):
        shifts = rng.integers(0, x.shape[0], size=x.shape[1])
        xs = np.column_stack([np.roll(x[:, j], s) for j, s in enumerate(shifts)])
        null_raw.append(_stats(xs))
        null_hourly.append(_stats(_hourly(xs, per_hour)))
    report = {
        "n_sensors": x.shape[1],
        "start": wide.index[0].isoformat(),
        "end": wide.index[-1].isoformat(),
        "n_null": N_NULL,
        "readings": _report(_stats(x), null_raw),
        "hourly_means": _report(_stats(_hourly(x, per_hour)), null_hourly),
    }
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
