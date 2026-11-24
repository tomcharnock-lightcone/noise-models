"""How far the Hann and Blackman noise models differ, sensor by sensor.

usage: python src/window_difference.py HANN.csv BLACKMAN.csv OUT.json

For σ_w and f_knee the difference is relative (Blackman / Hann − 1); for α
it is absolute. Each is also given in units of the Hann fit's bootstrap
uncertainty: much below 1 means the choice of window moves a parameter by
far less than the chamber run can measure it.
"""

from __future__ import annotations

import json
import sys

import pandas as pd


def main(hann_path: str, blackman_path: str, out: str) -> None:
    h = pd.read_csv(hann_path).set_index("sensor")
    b = pd.read_csv(blackman_path).set_index("sensor")
    if not h.index.equals(b.index):
        raise ValueError("the two tables do not list the same sensors")
    report = {}
    for name, col, relative in [("sigma_white", "sigma_white_c", True),
                                ("f_knee", "f_knee_hz", True),
                                ("alpha", "alpha", False)]:
        diff = b[col] - h[col]
        size = (diff / h[col]).abs() if relative else diff.abs()
        in_err = (diff / h[f"{col}_err"]).abs()
        worst = size.idxmax()
        report[name] = {
            ("max_abs_relative_difference" if relative else "max_abs_difference"):
                round(float(size.max()), 4),
            "sensor_of_max": worst,
            "max_abs_difference_in_errors": round(float(in_err.max()), 3),
            "median_abs_difference_in_errors": round(float(in_err.median()), 3),
        }
    report["n_sensors"] = len(h)
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
