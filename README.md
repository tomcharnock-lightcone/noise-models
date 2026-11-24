# noise-models

Per-sensor noise models for the 40 temperature sensors `s01` to `s40`, from
one week in a climate chamber at constant temperature. Each sensor's noise is
white plus 1/f, given by its power spectral density

    P(f) = (2 σ_w² / fs) · [1 + (f_knee / f)^α]      (°C²/Hz, one-sided)

- Models, one row per sensor with 1σ uncertainties, in
  `results/<universe>/noise_models.csv`. The universes cross the Welch window
  (Hann or Blackman, which agree within their uncertainties) with the 1/f
  slope α (fitted per sensor, or held at 1: `*_slope_one`).
- **Default, for every sensor: `results/hann_slope_one/noise_models.csv`**
  (α held at 1; the Blackman table is interchangeable).

**Caveat: s38.** Its free fit gives α = 1.043 ± 0.013, 3.3σ above 1, the only
sensor past 3σ. Among 40 sensors that happens by chance about one time in
ten, so the default keeps s38 at α = 1 like the rest. Holding its slope at 1
puts its σ_w 10% below the free fit and its f_knee 41% above, both about 3σ
of the free fit's uncertainty. If that matters for a project, s38's
free-slope row is in `results/hann/noise_models.csv`, and
`results/slope_test.json` has the comparison.
- Model code and a noise simulator: `src/noise_model.py`.
- Method, checks and caveats: [`index.md`](index.md). Specification: `astra.yaml`.

The data and results are stored with git-annex in this repository's Git LFS
storage. In a clone, run `lc fetch` to get the results and
`lc fetch chamber_run` to get the raw readings.
