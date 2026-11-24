# noise-models

Per-sensor noise models for the 40 temperature sensors `s01` to `s40`, from
one week in a climate chamber at constant temperature. Each sensor's noise is
white plus 1/f, given by its power spectral density

    P(f) = (2 σ_w² / fs) · [1 + (f_knee / f)^α]      (°C²/Hz, one-sided)

- Models, one row per sensor with 1σ uncertainties, in
  `results/<universe>/noise_models.csv`. The universes cross the Welch window
  (Hann or Blackman, which agree within their uncertainties) with the 1/f
  slope α (fitted per sensor, or held at 1: `*_slope_one`).
- **Default: `results/hann_slope_one/noise_models.csv`** (α = 1), except
  sensor s38, whose data reject α = 1. For s38 use its row in
  `results/hann/noise_models.csv`.
- Model code and a noise simulator: `src/noise_model.py`.
- Method, checks and caveats: [`index.md`](index.md). Specification: `astra.yaml`.

The data and results are stored with git-annex in this repository's Git LFS
storage. In a clone, run `lc fetch` to get the results and
`lc fetch chamber_run` to get the raw readings.
