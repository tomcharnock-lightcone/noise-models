"""Load the chamber readings as one column per sensor, refusing anything irregular."""

from __future__ import annotations

import numpy as np
import pandas as pd


def load(path: str) -> tuple[pd.DataFrame, float]:
    """Read a long `time,sensor,temp_c` CSV into a (time × sensor) frame.

    Returns the frame and the sampling rate in Hz. Raises ValueError on a
    missing or duplicated reading, a gap, or an uneven cadence: the Welch
    estimate assumes an unbroken, evenly sampled series for every sensor.
    """
    df = pd.read_csv(path, dtype={"sensor": str, "temp_c": float})
    df["time"] = pd.to_datetime(df["time"], utc=True)
    if df["temp_c"].isna().any():
        raise ValueError(f"{df['temp_c'].isna().sum()} readings are empty")
    dup = df.duplicated(["time", "sensor"]).sum()
    if dup:
        raise ValueError(f"{dup} (time, sensor) pairs appear more than once")
    wide = df.pivot(index="time", columns="sensor", values="temp_c").sort_index()
    if wide.isna().any().any():
        missing = wide.isna().sum()
        raise ValueError(f"readings missing: {missing[missing > 0].to_dict()}")
    dt = np.diff(wide.index.tz_convert(None).to_numpy()).astype("timedelta64[ms]").astype(np.int64)
    if len(np.unique(dt)) != 1:
        raise ValueError(f"uneven cadence: intervals of {sorted(set(dt.tolist()))} ms")
    return wide, 1000.0 / dt[0]
