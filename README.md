# noise-models

Per-sensor noise models for the 40 temperature sensors `s01` to `s40`, from
one week in a climate chamber at constant temperature. Each sensor's noise is
white plus 1/f, given by its power spectral density

    P(f) = (2 σ_w² / fs) · [1 + (f_knee / f)^α]      (°C²/Hz, one-sided)

- Models, one row per sensor with 1σ uncertainties:
  `results/hann/noise_models.csv` and `results/blackman/noise_models.csv`
  (Welch window; the two agree within their uncertainties).
- Model code and a noise simulator: `src/noise_model.py`.
- Method, checks and caveats: [`index.md`](index.md). Specification: `astra.yaml`.

The data and results are stored with git-annex in this repository's Git LFS
storage. In a clone, run `lc fetch` to get the results and
`lc fetch chamber_run` to get the raw readings.
