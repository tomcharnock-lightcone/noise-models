# noise-models

Noise models for the 40 temperature sensors `s01` to `s40`, measured in a
climate chamber held at constant temperature. The sensors keep these ids
outside the chamber, so other projects look a sensor's model up by its id.
The analysis is specified in `astra.yaml` (ASTRA 0.1) and run with
lightcone-cli (`lc`). This report quotes no measured value: the numbers live
in the results it names.

## Data

`data/chamber-run.csv` is one week (2022-12-12 to 2022-12-18, UTC) of
readings, one every 10 seconds per sensor, 60,480 per sensor with none
missing or duplicated, in steps of 0.01 °C. Columns `time`, `sensor`,
`temp_c`. Its git-annex key is its SHA-256, the one the delivery's
`SHA256SUMS` lists, and its bytes are in this repository's Git LFS storage:
`lc fetch chamber_run` brings them into a clone.

## The noise model

Every sensor's noise is coloured: on top of a white floor it rises toward
low frequency roughly as 1/f, so it cannot be summarised by one standard
deviation (a 1/f series' standard deviation keeps growing with the length of
the record). Each sensor is therefore described by its one-sided power
spectral density, in °C²/Hz,

$$
P(f) = \frac{2\sigma_w^2}{f_s}\left[1 + \left(\frac{f_\mathrm{knee}}{f}\right)^{\alpha}\right]
$$

- `sigma_white_c`, σ_w: the white noise per reading, in °C, at the chamber
  run's 10-second cadence (f_s = 0.1 Hz). `white_psd_c2_per_hz` is the same
  level as 2σ_w²/f_s, which carries over to other cadences provided the
  sensor's readout does not change.
- `f_knee_hz`: the frequency at which the 1/f part equals the white floor.
  Below it the 1/f part dominates.
- `alpha`: the slope of the 1/f part. Held at exactly 1 (the `slope`
  decision's `one`), the model is P(f) = W + K/f, with W = 2σ_w²/f_s the
  white floor and K = W · f_knee the 1/f amplitude.

`src/noise_model.py` holds the model (`psd_model`) and a simulator that
draws noise series from it (`simulate`), for projects that want to inject
realistic noise rather than only weight by it.

## Methods

**Spectrum.** `src/welch_psd.py` estimates each sensor's PSD with Welch's
method: one-day segments overlapping by half (13 over the week), each
segment's mean removed before windowing.

**The window is a decision, built both ways.** Steep low-frequency power can
leak through a window's sidelobes onto the white floor. Blackman's sidelobes
are far lower than Hann's, but its main lobe is wider. Which matters more
here was not known, so both are built.

**So is the 1/f slope.** α is either fitted per sensor with σ_w and f_knee
(`free`), or held at exactly 1 (`one`), the textbook 1/f. Holding it narrows
the other two parameters' uncertainties, but is only right if the data allow
α = 1, which `src/slope_test.py` tests. With both decisions there are four
universes:

| Universe | Window | Slope | Results |
|---|---|---|---|
| `hann` | Hann | fitted | `results/hann/` |
| `blackman` | Blackman | fitted | `results/blackman/` |
| `hann_slope_one` | Hann | held at 1 | `results/hann_slope_one/` |
| `blackman_slope_one` | Blackman | held at 1 | `results/blackman_slope_one/` |

The segment length is a second decision, `segment_length`. One day is built;
six hours is recorded in `astra.yaml` as an alternative and not built,
because it would start the fit above several sensors' knees.

**Fit.** `src/fit_noise.py` fits the model to each sensor's spectrum by
maximum Whittle likelihood, the likelihood of a Welch estimate, from the
third frequency bin above zero up to Nyquist, not including it. Removing
each segment's mean leaves a hole at zero frequency that the window spreads
over its main lobe, two bins for Hann and three for Blackman, so the fit
starts past both.

**Uncertainties** come from a parametric bootstrap: 200 series are simulated
from each sensor's fit and put through the same Welch estimate and fit, and
the spread of their parameters is the uncertainty. The textbook Fisher-matrix
errors were tried first and rejected: on simulated noise they understated the
scatter of α by a factor of about 1.6, because neighbouring frequency bins and
overlapping segments are correlated. The bootstrap's spread matches the
scatter between independent simulated realisations, a check that
`tests/test_noise_model.py` repeats.

**Is it the sensors, or the chamber?** If the chamber's temperature had
wandered, all 40 sensors would have recorded the same wander, and it would be
mistaken for sensor noise. `src/chamber_check.py` compares the spread of the
across-sensor mean with what independent noise alone would give, and looks
at the sensors' pairwise correlations, on the readings and on hourly means.
Because 1/f noise makes these statistics scatter widely, each is judged
against a null made by shifting every sensor's series in time by its own
random offset. That breaks anything the sensors share while keeping each
one's own spectrum.

## Results

**The noise models** are `results/<universe>/noise_models.csv`, one row per
sensor in each of the four universes. Beside each parameter is its 1σ
uncertainty (`_err`), and `knee_above_nyquist`, `low_band_ratio_z` and
`at_bound` flag the fits to read with care (below). Every column is
described under `noise_models` in `astra.yaml`.

```{figure} results/hann/psd_fits.png
:label: fig-psd-fits
Each sensor's Welch spectrum (Hann window), its fitted white + 1/f model,
the white floor (dashed) and the knee (dotted). The model is drawn over the
bins it was fitted to; the bins left of it are excluded (see Methods).
Each other universe has the same figure, `results/<universe>/psd_fits.png`.
```

Every sensor's spectrum is a white floor plus a 1/f rise, with α close to 1
throughout, and the fitted model follows the data across the whole band
([](#fig-psd-fits)). The sensors differ widely, both in their white level
and in how far up the 1/f part reaches.

**The window does not matter here.** For every sensor and every parameter,
the Hann and Blackman fits differ by well under the fits' own statistical
uncertainty, whether the slope is fitted or held at 1
([](#fig-window-comparison), `results/window_difference.json`, finding
`window_choice_immaterial` in `astra.yaml`). The spectra are close
enough to 1/f that Hann's higher sidelobes leak too little power to move the
white floor. Both universes are kept, and either set of models can be used.

```{figure} results/window_comparison.png
:label: fig-window-comparison
Blackman against Hann, sensor by sensor, for each parameter, against the
Hann fit's ±1σ bootstrap uncertainty (grey).
```

**The slope can be held at 1, and should be by default, except for s38.**
Fitted freely, the slopes cluster around 1. Holding them there is the
`slope` decision's option `one`, built in the universes `hann_slope_one`
and `blackman_slope_one`. `results/slope_test.json` (finding
`slope_one_allowed`) asks whether the data allow it:

- Taken together, the 40 free slopes are consistent with exactly 1, and
  their precision-weighted mean sits on it. The slopes' intrinsic scatter
  about their common value is bounded tightly enough that holding every
  sensor at one slope costs at most a few percent in the 1/f power, even a
  decade below the measured band.
- Three sensors' free slopes lie beyond 2σ of 1, against about two expected
  by chance among 40. One, s38, lies beyond 3σ, and holding its slope moves
  its σ_w and f_knee by about 3σ of their free uncertainty. Under the held
  model its spectrum, and s37's, also misfit the lowest frequencies
  (`low_band_ratio_z` past 2). Without s38, the other 39 agree with 1 as
  well as chance allows (`population.without_most_discrepant`).
- Holding the slope narrows the uncertainties of σ_w and f_knee by about a
  third, and changes the model only slightly within the measured band.

```{figure} results/slope_comparison.png
:label: fig-slope-comparison
Top: each sensor's free slope against 1. Below: what holding it at 1 does to
f_knee and σ_w, against the free fit's ±1σ (grey). Hann window.
```

**The chamber held still.** Neither the readings nor their hourly means show
anything shared between sensors beyond what independent noise gives: every
statistic in `results/chamber_check.json` lies inside its null's 95% range
(finding `chamber_held_still`). The noise models therefore describe the
sensors, not the chamber.

### Fits to read with care

- `knee_above_nyquist` is true for a few sensors: their white floor lies
  above every frequency the 10-second readings reach, so their σ_w is the
  model extrapolated past the data. It is still constrained, through the
  curvature of the spectrum near Nyquist, but less tightly, and its `_err`
  says by how much.
- `low_band_ratio_z` tests the 1/f shape over the lowest decade fitted. With
  the slope fitted, s15 and s31 sit just past |z| = 2 under both windows,
  about as many as chance gives among 40. They are worth a second look if the
  chamber run is ever repeated. With the slope held at 1 it is s37 and s38
  instead, the two sensors whose free slopes lie furthest above 1.
- `at_bound` is empty for every sensor: no fit ran into a parameter limit.

### Using the models elsewhere

- **By default use the slope held at 1**, `results/hann_slope_one/noise_models.csv`
  (`results/blackman_slope_one/noise_models.csv` is interchangeable with it),
  and look a sensor up by its id. **For s38 use its row in the free-slope
  table**, `results/hann/noise_models.csv`: its data reject α = 1, and holding
  it there shifts its σ_w and f_knee well outside their uncertainties. s37 is
  borderline. Its free and held models agree within about 2σ.
- A held-slope model's uncertainties are only as good as the assumption that
  α is exactly 1. The data allow the sensors' slopes to scatter only a
  little (`slope_test.json`, `intrinsic_scatter_95_upper`), and that bounds
  the error well inside the measured band. Below the band, at periods longer
  than about a day, the gap grows with the distance extrapolated. A use case
  that extrapolates far should compare the free-slope model.
- For noise per reading at another cadence f_s′, the white part's standard
  deviation is √(`white_psd_c2_per_hz` · f_s′ / 2). This holds only if the
  sensor's readout (its internal averaging) is the same as in the chamber.
- The 1/f part means averages beat down slower than 1/√N: over long
  intervals the scatter of a mean is dominated by the 1/f part, not σ_w.
  Integrate the model over the frequencies the use case is sensitive to, or
  simulate from it with `noise_model.simulate`.
- The chamber run covers one week at one temperature, about 20 °C. The models
  say nothing about how the noise changes with temperature or over longer
  than a week.

## Reproducing

`lc status` reports where each result stands, `lc fetch` brings results'
bytes into a clone, and `lc materialize <cluster>` remakes whatever is
stale. `uv run pytest` runs the model's tests.
