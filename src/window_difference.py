"""How far the Hann and Blackman noise models differ, sensor by sensor.

usage: python src/window_difference.py HANN.csv BLACKMAN.csv [HANN.csv BLACKMAN.csv ...] OUT.json

Each pair of tables is fitted the same way but for the window; the report has
one section per pair, named by its slope (slope_free or slope_one). For σ_w
and f_knee the difference is relative (Blackman / Hann − 1); for a fitted α
it is absolute. Each is also given in units of the Hann fit's bootstrap
uncertainty: much below 1 means the choice of window moves a parameter by
far less than the chamber run can measure it.
"""

from __future__ import annotations

import json
import sys

import pandas as pd


def compare(h: pd.DataFrame, b: pd.DataFrame) -> dict:
    if not h.index.equals(b.index):
        raise ValueError("the two tables do not list the same sensors")
    if (h.window == b.window).any() or not h.alpha_fixed.equals(b.alpha_fixed):
        raise ValueError("a pair must differ in window and nothing else")
    params = [("sigma_white", "sigma_white_c", True), ("f_knee", "f_knee_hz", True)]
    if not h.alpha_fixed.any():
        params.append(("alpha", "alpha", False))
    report = {}
    for name, col, relative in params:
        diff = b[col] - h[col]
        size = (diff / h[col]).abs() if relative else diff.abs()
        in_err = (diff / h[f"{col}_err"]).abs()
        report[name] = {
            ("max_abs_relative_difference" if relative else "max_abs_difference"):
                round(float(size.max()), 4),
            "sensor_of_max": size.idxmax(),
            "max_abs_difference_in_errors": round(float(in_err.max()), 3),
            "median_abs_difference_in_errors": round(float(in_err.median()), 3),
        }
    return report


def main(*args: str) -> None:
    *tables, out = args
    if not tables or len(tables) % 2:
        raise SystemExit(__doc__)
    report = {}
    for hann_path, blackman_path in zip(tables[::2], tables[1::2]):
        h = pd.read_csv(hann_path).set_index("sensor")
        b = pd.read_csv(blackman_path).set_index("sensor")
        key = "slope_one" if h.alpha_fixed.all() else "slope_free"
        report[key] = compare(h, b)
        report["n_sensors"] = len(h)
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
