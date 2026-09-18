# tweezer-experiment (experiment_toolkit)

Offline experiment planning toolkit for optical-tweezer transport, extracted
from the research repository `monte-carlo-for-aodslm` (branch
`movement_effects_v1`, commit `1fed632`, physics of the 2026-09-16
calibration study). Phase 1 scope: **mechanical motion, survival and
excitation of Cs-like atoms** — input validation, protocol compilation,
forward prediction, constrained control search, and waveform export with
read-back verification. Spin/DD/RB physics and real-device control are
explicit non-goals of phase 1.

This directory is standalone: copy it out, build a wheel, install, run. It
never imports research scripts, never scans historical outputs, and reads
paper data only through the optional `adapters.paper` module with
user-supplied paths.

## Install

```bash
python -m pip install .            # runtime deps: numpy, PyYAML
python -m pip install .[test]      # + pytest
tweezer-experiment selftest        # 5 fast unit anchors, exit 0
```

Python >= 3.10, numpy-only backend (no numba/scipy required; scipy is used
only for `cubic` CSV interpolation when present).

## Five commands

```bash
tweezer-experiment validate examples/static_hold.yaml        # F1: errors/capabilities
tweezer-experiment config   examples/static_hold.yaml        # effective config + sources
tweezer-experiment simulate examples/short_transfer.yaml --out runs/a
tweezer-experiment compare  examples/static_hold.yaml examples/short_transfer.yaml  # same-seed diff
tweezer-experiment optimize examples/search_pickup.yaml --out runs/s
tweezer-experiment export   examples/export_demo.yaml runs/e # F5 + read-back verify
tweezer-experiment report   runs/a                           # reprint a stored result
```

Python API (CLI calls the same layer):

```python
from tweezer_experiment import (validate_experiment, compile_protocol, simulate,
                                search_controls, export_controls, ExperimentBundle)
bundle = ExperimentBundle.load("examples/short_transfer.yaml")
report = validate_experiment(bundle)      # ValidationReport
result = simulate(bundle)                 # SimulationResult (arrays, no files forced)
search = search_controls(bundle)          # SearchResult (requires search: section)
manifest = export_controls(bundle, "runs/e")  # physical table (+RF if calibrated)
```

Importing the package has no side effects (no runs, no cwd change, no global
RNG state) — verified by test.

## Input contract (units and sources)

One YAML with sections `device, preparation, noise, protocol, numerics`
(required) and `calibration, search` (optional). **Field names carry their
unit** (`temperature_uK`, `waist_um`, `duration_us`, `psd_per_hz`); values
are converted to SI at parse time. Non-SI combinations and unknown fields
fail with a stable `E_*` code and the dotted field path (see
`errors.py`). Note YAML 1.1 quirks: write `2.0e+7` (not `2.0e7`) and quote
bare `off`.

| Section | Key fields |
| --- | --- |
| `device` | `species` (Cs133 preset or `mass_u`), per-beam `wavelength_nm/waist_um/depth_uK/center_um`, `lensing {mode: off\|velocity_shift, vs_m_s}` |
| `preparation` | `temperature_uK`, site disorder sigmas (`slm_depth_rel_sigma`, `aod_depth_rel_sigma`, `aod_waist_rel_sigma`, `alignment_sigma_um`, `shot_jitter_alignment_sigma_nm`), sampler `site_bound_harmonic` |
| `noise` | `model: none \| diffusion \| local_mean \| legacy_constant_power`; per-beam `*_rin_psd_per_hz` / `*_pointing_psd_m2_per_hz` (one-sided per Hz), `recoil_per_joule`; `legacy_power_w` only for the legacy model |
| `protocol` | segments (`static`, `transfer` with `path` family + `ramp_fraction` + `order: depth_then_move\|move_then_depth`, `csv` with `interpolation`), optional `observe: aod\|slm` per segment, `sequence` with nested `repeat`, optional `final_hold_us` |
| `numerics` | `dt_us`, `noise_dt_us` (integer multiple of dt), `shots`, `seed`, `backend`, `save_trajectories` |
| `calibration` | `clock_hz`, `quantization_bits`, monotone `channels` (`aod_x`, `aod_y`, `aod_depth`) with `input/input_unit/output/unit`, `limits` |
| `search` | `variables` (`segment.<id>.duration_us/.ramp_fraction/.depth_to_uK/.path`), `objective` (maximize_capture or minimize_duration+bound), pool sizes/seeds, `max_candidates` |

Every numeric field may carry a companion `<field>__source:
{tag: measured|derived|assumed|fitted|synthetic, note: ...}`. The
capability classifier in `validate` reports which measured inputs are
missing, so a run with all-assumed parameters is labeled an idealized
simulation, never a calibrated prediction.

## Physics included (and its provenance)

- 3D Gaussian tweezer fields with optional acousto lensing z-shifts
  (`s = 2 z_R v_cmd / v_s`), wavelengths/waists/depths fully parametric.
- Site-bound initial-state preparation (drawn site depths, harmonic
  proposal, E<0 rejection, retry at the same site).
- Velocity-Verlet with segment-boundary switch ledger and single-beam
  binding judges at observed checkpoints (absorbing by default;
  `absorb=False` for recapture-allowed diagnosis).
- Stochastic force diffusion (research model of 2026-09-16): one-sided PSDs
  per beam, intensity kicks through each beam's force vector, pointing
  kicks through dF/dx, recoil proportional to local depth, noise step
  independent of the mechanical step. `local_mean` is the deterministic
  diagnostic of the same channels; `legacy_constant_power` is the
  historical fixed-power heating, kept only as an explicitly named legacy
  mode.
- Move families: `paper_manual` (3u^2-2u^3, endpoint acceleration jumps —
  reported, not smoothed), `minimum_jerk`, `constant_jerk` (4-segment
  symmetric), `adiabatic_sine`, `linear` (+`smoothstep` alias).

Known limits (honest): no spin-dependent forces, no common-path
batch noise (per-atom independent streams only — MC intervals do not
describe day-to-day experimental batches), no vacuum/imaging loss model,
`cubic` CSV interpolation requires scipy, one species preset is a mass
only (Rb/Yb mechanisms are not validated).

## Examples (`examples/`)

| file | demonstrates |
| --- | --- |
| `static_hold.yaml` | smallest prediction; diffusion heating baseline |
| `short_transfer.yaml` | A07: 350/410 µs asymmetric pickup/dropoff, 3.1 µm, 82 µs wait, 13 repeats |
| `long_transport.yaml` | AOD-only 375 µm diagonal transport with far-end observation |
| `csv_trajectory.yaml` + `.csv` | user CSV import, derived velocities, extra observation |
| `search_pickup.yaml` | constrained 8-candidate search with separate validation pool |
| `export_demo.yaml` | synthetic calibration, RF table, limits, read-back verification |

All example numbers are synthetic interface validation, not experimental
recommendations.

## Tests and CI

```bash
python -m pytest tests/ -q        # 48 tests: units, schemas, physics
                                  # benchmarks (A05), semantics, search,
                                  # export round-trip, CLI parity, purity
```

`.github/workflows/ci.yml` runs the same suite on push (GitHub Actions).
Physics benchmarks include independent expectations: finite-difference
force checks, oscillation frequency vs the analytic trap frequency,
second-order Verlet convergence, diffusion ledger realized-vs-expected,
bitwise zero-noise limits, exact pickup/dropoff time-mirror equality.

## Provenance of the extraction

Physics kernels were extracted from
`simulation/level2_joint_transfer/src/continuous_transfer/{engine,noisy_survival,
initialization,protocol,paper_waveforms}.py` (engine/protocol semantics
unchanged; wavelengths and beam parameters parametrized; checkpoints
generalized from the matched 2R+1 layout to arbitrary observed segments).
The audit trail and acceptance report live in `acceptance/`.
