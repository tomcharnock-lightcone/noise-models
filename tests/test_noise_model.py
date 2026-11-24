"""The fit recovers what it is given, and the loader refuses irregular readings."""

import numpy as np
import pytest

import noise_model
import readings

FS = 0.1
WEEK = 60_480


def test_simulated_white_noise_has_the_requested_variance():
    rng = np.random.default_rng(0)
    # A knee far below the lowest frequency leaves white noise.
    x = noise_model.simulate(rng, WEEK, FS, 0.2, 1e-12, 1.0, n_series=4)
    assert x.shape == (WEEK, 4)
    assert np.std(x) == pytest.approx(0.2, rel=0.01)


@pytest.mark.parametrize("window", ["hann", "blackman"])
@pytest.mark.parametrize("sigma, f_knee, alpha", [(0.15, 5e-3, 1.0), (0.08, 7e-2, 1.1)])
def test_fit_recovers_simulated_parameters(window, sigma, f_knee, alpha):
    rng = np.random.default_rng(1)
    x = noise_model.simulate(rng, WEEK, FS, sigma, f_knee, alpha, n_series=20)
    f, p, info = noise_model.welch(x, FS, window, noise_model.SEGMENT_SECONDS["one_day"])
    assert info["n_segments"] == 13
    fits = [noise_model.fit_psd(f, pi, FS) for pi in p]
    assert all(r["at_bound"] == "" for r in fits)
    est = np.median([[r["sigma_white"], r["f_knee"], r["alpha"]] for r in fits], axis=0)
    assert est[0] == pytest.approx(sigma, rel=0.03)
    assert est[1] == pytest.approx(f_knee, rel=0.15)
    assert est[2] == pytest.approx(alpha, abs=0.03)


def test_bootstrap_spread_matches_scatter_between_realisations():
    rng = np.random.default_rng(2)
    x = noise_model.simulate(rng, WEEK, FS, 0.2, 2e-3, 1.0, n_series=60)
    f, p, _ = noise_model.welch(x, FS, "hann", 86_400)
    theta = np.array([noise_model.fit_psd(f, pi, FS)["theta"] for pi in p])
    boot = noise_model.bootstrap(noise_model.fit_psd(f, p[0], FS), WEEK, FS, "hann", 86_400, 60, rng)
    ratio = boot[:, :3].std(axis=0) / theta.std(axis=0)
    # Two 60-draw spreads: each is good to ~10%, so their ratio to ~30% at 2σ.
    assert np.all((ratio > 0.65) & (ratio < 1.5)), ratio


def test_welch_dof_is_two_per_segment_without_overlap():
    w = np.hanning(100)
    assert noise_model.welch_dof(w, step=100, n_segments=7) == pytest.approx(14)


def _write(tmp_path, rows):
    path = tmp_path / "r.csv"
    path.write_text("time,sensor,temp_c\n" + "\n".join(rows) + "\n")
    return str(path)


def test_loader_reads_a_regular_series(tmp_path):
    rows = [f"2022-12-12T00:00:{s:02d}Z,{sid},20.0" for s in (0, 10, 20) for sid in ("s01", "s02")]
    wide, fs = readings.load(_write(tmp_path, rows))
    assert list(wide.columns) == ["s01", "s02"]
    assert fs == pytest.approx(0.1)


@pytest.mark.parametrize("rows, message", [
    (["2022-12-12T00:00:00Z,s01,20.0", "2022-12-12T00:00:10Z,s01,20.0",
      "2022-12-12T00:00:30Z,s01,20.0"], "uneven cadence"),
    (["2022-12-12T00:00:00Z,s01,20.0", "2022-12-12T00:00:00Z,s01,20.1"], "more than once"),
    (["2022-12-12T00:00:00Z,s01,20.0", "2022-12-12T00:00:00Z,s02,20.0",
      "2022-12-12T00:00:10Z,s01,20.0"], "missing"),
    (["2022-12-12T00:00:00Z,s01,"], "empty"),
])
def test_loader_refuses_irregular_readings(tmp_path, rows, message):
    with pytest.raises(ValueError, match=message):
        readings.load(_write(tmp_path, rows))
