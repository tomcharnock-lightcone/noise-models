"""The per-sensor noise model, shared by the recipes and by anyone reusing the fits.

Each sensor's noise is modelled as white noise plus a 1/f-type ("flicker")
component, through its one-sided power spectral density in °C²/Hz:

    P(f) = (2 σ_w² / fs) · [1 + (f_knee / f)^α]

σ_w is the standard deviation, in °C, of the white part per sample at sampling
rate fs (Hz); 2 σ_w² / fs is the white floor, which does not depend on fs.
f_knee is the frequency at which the 1/f part equals the white floor, and α is
the slope of the 1/f part.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize, signal

#: Segment lengths for the Welch estimate, by the option ids of the
#: `segment_length` decision, in seconds.
SEGMENT_SECONDS = {"one_day": 86_400, "six_hours": 21_600}

#: Lowest Welch bin the fit uses. Removing each segment's mean zeroes the DC
#: bin, and the window spreads that hole over its main lobe: two bins either
#: side of DC for Hann, three for Blackman. The fit starts past both.
FIRST_FIT_BIN = 3


def psd_model(f, sigma_white, f_knee, alpha, fs):
    """One-sided noise PSD in °C²/Hz at frequencies f (Hz)."""
    return 2.0 * sigma_white**2 / fs * (1.0 + (f_knee / f) ** alpha)


def welch_dof(window: np.ndarray, step: int, n_segments: int) -> float:
    """Equivalent degrees of freedom of a Welch average of overlapping segments.

    Percival & Walden (1993), eq. 292b: each segment's periodogram is χ²₂, and
    overlap correlates neighbouring segments by the window's overlap integral.
    """
    h = window / np.sqrt(np.sum(window**2))
    n = len(h)
    s = 0.0
    for m in range(1, n_segments):
        shift = m * step
        if shift >= n:
            break
        s += (1 - m / n_segments) * np.sum(h[: n - shift] * h[shift:]) ** 2
    return 2 * n_segments / (1 + 2 * s)


def welch(x: np.ndarray, fs: float, window_name: str, segment_s: float):
    """Welch PSD of each column of x (time along axis 0), with its statistics.

    Segments overlap by half; each segment's mean is removed before windowing.
    Returns (f, psd[n_columns, n_f], info) where info holds the segmenting and
    the estimate's equivalent degrees of freedom.
    """
    nperseg = int(round(segment_s * fs))
    if nperseg > x.shape[0]:
        raise ValueError(f"segment of {nperseg} samples is longer than the series ({x.shape[0]})")
    noverlap = nperseg // 2
    step = nperseg - noverlap
    n_segments = 1 + (x.shape[0] - nperseg) // step
    win = signal.get_window(window_name, nperseg)
    f, p = signal.welch(
        x, fs=fs, window=win, nperseg=nperseg, noverlap=noverlap,
        detrend="constant", scaling="density", return_onesided=True, axis=0,
    )
    info = {
        "nperseg": nperseg,
        "noverlap": noverlap,
        "n_segments": n_segments,
        "dof": welch_dof(win, step, n_segments),
    }
    return f, p.T, info


def fit_band(f: np.ndarray, fs: float) -> np.ndarray:
    """Mask of the bins the fit uses: from FIRST_FIT_BIN up to, not including, Nyquist.

    scipy does not double the Nyquist bin of a one-sided density, and its
    periodogram is χ²₁ rather than χ²₂, so it is left out.
    """
    idx = np.arange(len(f))
    return (idx >= FIRST_FIT_BIN) & (f < fs / 2)


def _log_model_and_jacobian(theta, f, fs):
    """ln P(f) and its derivatives with respect to θ = (ln σ_w², ln f_knee, α)."""
    ln_s2, ln_fk, alpha = theta
    ln_ratio = ln_fk - np.log(f)
    u = np.exp(np.clip(alpha * ln_ratio, -700, 700))
    g = u / (1.0 + u)
    ln_m = np.log(2.0 / fs) + ln_s2 + np.log1p(u)
    jac = np.stack([np.ones_like(f), alpha * g, g * ln_ratio])
    return ln_m, jac


def simulate(rng: np.random.Generator, n_samples: int, fs: float, sigma_white: float,
             f_knee: float, alpha: float, n_series: int = 1) -> np.ndarray:
    """Draw Gaussian noise series with the model's PSD, shape (n_samples, n_series).

    Built in the Fourier domain over twice the length and cut to the first
    half, so the series is not periodic over its own length. The mean (DC
    term) is zero.
    """
    n = 2 * n_samples
    f = np.fft.rfftfreq(n, 1 / fs)
    p = np.zeros_like(f)
    p[1:] = psd_model(f[1:], sigma_white, f_knee, alpha, fs)
    # E|X_k|² = n fs P / 2 for a one-sided P; real and imaginary parts share it.
    amp = np.sqrt(n * fs * p / 4)
    x = amp * (rng.standard_normal((n_series, len(f))) + 1j * rng.standard_normal((n_series, len(f))))
    x[:, -1] = np.sqrt(n * fs * p[-1] / 2) * rng.standard_normal(n_series)
    return np.fft.irfft(x, n=n, axis=1)[:, :n_samples].T


#: Bounds of the fit in θ = (ln σ_w², ln f_knee, α); ln f_knee's are set per fit.
ALPHA_BOUNDS = (0.05, 4.0)


def _nll(theta, f, p, fs):
    ln_m, jac = _log_model_and_jacobian(theta, f, fs)
    r = p * np.exp(-ln_m)
    return np.sum(ln_m + r), jac @ (1.0 - r), jac


def fit_psd(f: np.ndarray, p: np.ndarray, fs: float) -> dict:
    """Fit the white + 1/f model to a Welch PSD by maximum Whittle likelihood.

    A Welch estimate is P̂ = P · χ²_ν / ν at each bin, so the negative log
    likelihood is (ν/2) Σ [ln P + P̂/P] up to a constant, whose minimum does
    not depend on ν. L-BFGS-B finds it, then Fisher scoring polishes it.
    Uncertainties are not computed here: neighbouring bins and segments are
    correlated, so the independent-bin Fisher matrix understates them, and
    `bootstrap` measures them instead.
    """
    band = fit_band(f, fs)
    fb, pb = f[band], p[band]

    # Start: white floor from the top decade's median, knee where the
    # spectrum last sits above twice that floor.
    floor = np.median(pb[fb > fb[-1] / 10])
    above = np.nonzero(pb > 2 * floor)[0]
    fk0 = fb[above[-1]] if len(above) else fb[0]
    theta0 = np.array([np.log(floor * fs / 2), np.log(fk0), 1.0])

    lo = np.array([-np.inf, np.log(fb[0] / 1e3), ALPHA_BOUNDS[0]])
    hi = np.array([np.inf, np.log(fs * 1e3), ALPHA_BOUNDS[1]])
    res = optimize.minimize(lambda t: _nll(t, fb, pb, fs)[:2], theta0, jac=True,
                            method="L-BFGS-B", bounds=list(zip(lo, hi)))
    theta = np.clip(res.x, lo, hi)

    # Fisher scoring: the expected Hessian of Σ[ln P + P̂/P] is Σ ∂lnP ∂lnPᵀ.
    # The Newton decrement is twice the expected drop still to come; one unit
    # of the likelihood is ~ν/2 ≈ 10 of these, so 1e-9 is far below anything
    # that moves a parameter.
    converged = False
    for _ in range(50):
        value, grad, jac = _nll(theta, fb, pb, fs)
        step = -np.linalg.solve(jac @ jac.T, grad)
        decrement = -grad @ step
        if decrement < 1e-9:
            converged = True
            break
        t = 1.0
        while t > 1e-6:
            trial = np.clip(theta + t * step, lo, hi)
            if _nll(trial, fb, pb, fs)[0] < value:
                break
            t /= 2
        else:
            # No step lowers it: rounding, or a parameter held at its bound.
            pinned = np.isclose(theta, lo) | np.isclose(theta, hi)
            converged = decrement < 1e-6 or bool(pinned.any())
            break
        theta = trial
    if not converged:
        raise RuntimeError(f"PSD fit did not converge (L-BFGS-B: {res.message})")

    ln_m, _ = _log_model_and_jacobian(theta, fb, fs)
    at_bound = [name for name, v, a, b in zip(("ln_sigma2", "ln_f_knee", "alpha"), theta, lo, hi)
                if np.isclose(v, a) or np.isclose(v, b)]
    return {
        "theta": theta,
        "sigma_white": float(np.exp(theta[0] / 2)),
        "f_knee": float(np.exp(theta[1])),
        "alpha": float(theta[2]),
        "fit_f_min": float(fb[0]),
        "fit_f_max": float(fb[-1]),
        "n_fit_bins": int(band.sum()),
        # Over the whole band the mean of P̂/P is 1 by construction (the
        # σ_w equation of the fit says so). Over the lowest decade it is a
        # test of the 1/f shape where it is most stretched.
        "low_band_ratio": float(np.mean((pb / np.exp(ln_m))[fb < 10 * fb[0]])),
        "at_bound": ";".join(at_bound),
    }


def bootstrap(fit: dict, n_samples: int, fs: float, window_name: str, segment_s: float,
              n_boot: int, rng: np.random.Generator) -> np.ndarray:
    """Parametric bootstrap: refit n_boot series simulated from a fit.

    Returns (ln σ_w², ln f_knee, α, low_band_ratio) for each draw, shape
    (n_boot, 4). The spread is the fit's uncertainty, with the window, the overlap and the
    correlation between bins all accounted for, because the simulated series
    go through the same Welch estimate and fit as the data.
    """
    x = simulate(rng, n_samples, fs, fit["sigma_white"], fit["f_knee"], fit["alpha"], n_boot)
    f, p, _ = welch(x, fs, window_name, segment_s)
    fits = [fit_psd(f, pi, fs) for pi in p]
    return np.array([[*r["theta"], r["low_band_ratio"]] for r in fits])
