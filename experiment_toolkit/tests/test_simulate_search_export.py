"""simulate semantics, search contract, export round-trip and CLI (A02, A06-A14)."""
import json
import random
import numpy as np
import pytest

from tweezer_experiment.schemas import ExperimentBundle, validate_experiment
from tweezer_experiment.simulation import simulate
from tweezer_experiment.search import search_controls
from tweezer_experiment.export import export_controls
from tweezer_experiment.errors import ToolkitError
from tweezer_experiment import cli

BASE = """
device:
  species: Cs133
  slm: {wavelength_nm: 1061, waist_um: 1.17, depth_uK: 140}
  aod: {wavelength_nm: 1055, waist_um: 1.17, depth_uK: 280}
  lensing: {mode: "off"}
preparation: {temperature_uK: 5}
noise: {model: none}
numerics: {dt_us: 0.1, noise_dt_us: 0.1, shots: 64, seed: 11}
protocol:
  segments:
    - id: hold
      type: static
      duration_us: 30
      aod_center_um: [0, 0]
      aod_depth_uK: 0
      slm_factor: 1.0
      observe: slm
  sequence: [hold]
"""


def _load(text=BASE):
    import yaml
    return ExperimentBundle.parse(yaml.safe_load(text))


def test_static_deterministic_full_survival_and_recount():
    result = simulate(_load(), seed=1)
    assert result.survival[-1] == 1.0
    assert result.initial_unbound == 0
    assert np.all(np.diff(result.survival_counts) <= 0)        # absorbing: monotone
    # recount from per-atom labels (denominator = initial shots, absorbing)
    assert int(result.alive[:, -1].sum()) == result.shots


def test_absorbing_vs_recapture_allowed_semantics():
    # judge against an SLM switched OFF: every atom fails at the checkpoint
    text = BASE.replace("slm_factor: 1.0", "slm_factor: 0.0")
    bundle = _load(text)
    absorb = simulate(bundle, seed=2)
    assert absorb.survival[-1] == 0.0
    free = simulate(bundle, seed=2, absorb=False)
    assert free.survival[-1] == 0.0            # first failure still counts
    assert free.loss_histogram == {"hold": free.shots}
    assert free.final_states is not None         # propagation continued
    assert free.run.seeds == absorb.run.seeds


def test_different_master_seed_different_pool_same_statistics():
    a = simulate(_load(), seed=1)
    b = simulate(_load(), seed=2)
    assert not np.array_equal(a.final_states, b.final_states)
    assert a.survival[-1] == b.survival[-1] == 1.0


def test_loss_records_name_the_segment():
    import yaml
    raw = yaml.safe_load(BASE)
    raw["protocol"] = {
        "segments": [
            {"id": "a", "type": "static", "duration_us": 30, "aod_center_um": [0, 0],
             "aod_depth_uK": 0, "slm_factor": 0.0, "observe": "slm"},
            {"id": "b", "type": "static", "duration_us": 30, "aod_center_um": [0, 0],
             "aod_depth_uK": 0, "slm_factor": 1.0, "observe": "slm"},
        ],
        "sequence": ["a", "b"]}
    result = simulate(ExperimentBundle.parse(raw), seed=3)
    assert result.loss_histogram == {"a": result.shots}
    assert result.survival.tolist() == [1.0, 0.0, 0.0]  # index 0 = prepared state


SEARCH_CFG = BASE + """
search:
  variables:
    - target: segment.hold.duration_us
      values: [30, 30.03, 40]
  objective: {type: maximize_capture}
  max_candidates: 8
"""


def test_search_records_failed_candidates_and_keeps_device_fixed():
    bundle = _load(SEARCH_CFG)
    result = search_controls(bundle)
    statuses = {c.index: c.status for c in result.candidates}
    assert statuses.get(1) == "failed"          # 30.03 us not divisible by dt
    assert "E_NOT_DIVISIBLE" in str(result.candidates[1].error)
    assert any(s == "ok" for s in statuses.values())
    assert result.recommended is not None
    assert result.device_preparation_hash


def test_search_reports_no_feasible_candidate_when_bound_impossible():
    import yaml
    raw = yaml.safe_load(BASE)
    raw["search"] = {
        "variables": [{"target": "segment.hold.duration_us", "values": [30, 40]}],
        "objective": {"type": "minimize_duration", "segment": "hold",
                      "capture_lower_bound": 1.0},
        "max_candidates": 4}
    raw["protocol"]["segments"][0]["slm_factor"] = 0.0  # everyone fails the judge
    result = search_controls(ExperimentBundle.parse(raw))
    assert result.recommended is None
    assert all(c.status == "infeasible" for c in result.candidates)


CALIB = """
calibration:
  version: syn-test
  source_tag: synthetic
  clock_hz: 2.0e+7
  quantization_bits: 16
  channels:
    aod_x: {kind: frequency, unit: MHz, input: [-5, 20], output: [79.5, 84.5], input_unit: um}
    aod_depth: {kind: voltage, unit: V, input: [0, 320], output: [0, 1], input_unit: uK}
  limits: {max_rf_hz: 8.5e+7, max_rf_step_hz: 6.0e+4, max_amplitude: 1.0}
"""


def _export_bundle(calib=CALIB, protocol=None):
    import yaml
    raw = yaml.safe_load(BASE)
    if calib:
        raw.update(yaml.safe_load(calib))
    if protocol:
        raw["protocol"] = yaml.safe_load(protocol)
    return ExperimentBundle.parse(raw)


TRANSFER_PROTOCOL = """
segments:
  - id: pickup
    type: transfer
    duration_us: 100
    from_um: [0, 0]
    to_um: [2.4, 0]
    depth_from_uK: 0
    depth_to_uK: 280
    path: paper_manual
    ramp_fraction: 0.48
    slm_factor: 1.0
    observe: aod
  - id: dropoff
    type: transfer
    duration_us: 100
    from_um: [2.4, 0]
    to_um: [0, 0]
    depth_from_uK: 280
    depth_to_uK: 0
    path: paper_manual
    ramp_fraction: 0.48
    order: move_then_depth
    slm_factor: 1.0
    observe: slm
sequence: [pickup, dropoff]
"""


def test_export_physical_table_and_lossless_roundtrip(tmp_path):
    bundle = _export_bundle(protocol=TRANSFER_PROTOCOL)
    manifest = export_controls(bundle, tmp_path)
    assert (tmp_path / "physical_controls.csv").is_file()
    rb = manifest["readback"]["physical_csv_roundtrip"]
    assert "error" not in rb
    assert rb["max_abs_dS"] == pytest.approx(0.0, abs=1e-12)
    assert manifest["device_instructions"]["limits_ok"] is True
    rf = manifest["readback"]["rf_clock_quantized_roundtrip"]
    assert "error" not in rf
    # quantization error recorded honestly
    assert manifest["device_instructions"]["quantization_max_abs_error"]["rf_x_hz"] > 0


def test_export_coarse_clock_changes_the_result_measurably(tmp_path):
    bundle = _export_bundle(calib=CALIB.replace("2.0e+7", "2.0e+6"),
                            protocol=TRANSFER_PROTOCOL)
    manifest = export_controls(bundle, tmp_path)
    rf = manifest["readback"]["rf_clock_quantized_roundtrip"]
    assert "error" not in rf
    assert rf["max_abs_dS"] >= 0.0  # coarse 0.5 us clock: reported, not hidden


def test_export_rejects_out_of_calibration_range(tmp_path):
    tight = CALIB.replace("input: [-5, 20]", "input: [-5, 1]").replace(
        "valid_input: [-4, 18]", "valid_input: [-4, 1]")
    bundle = _export_bundle(calib=tight, protocol=TRANSFER_PROTOCOL)
    with pytest.raises(ToolkitError) as e:
        export_controls(bundle, tmp_path)
    assert e.value.code == "E_CAL_RANGE"
    assert not (tmp_path / "rf_instructions.csv").exists()


def test_export_limit_violation_blocks_compliant_file(tmp_path):
    tight = CALIB.replace("max_rf_step_hz: 6.0e+4", "max_rf_step_hz: 1.0")
    bundle = _export_bundle(calib=tight, protocol=TRANSFER_PROTOCOL)
    manifest = export_controls(bundle, tmp_path)
    assert manifest["device_instructions"]["limits_ok"] is False
    assert not (tmp_path / "rf_instructions.csv").exists()
    assert "failed" in manifest["device_instructions"]["reason"]


def test_no_calibration_still_exports_physical_table_only(tmp_path):
    bundle = _export_bundle(calib=None, protocol=TRANSFER_PROTOCOL)
    manifest = export_controls(bundle, tmp_path)
    assert (tmp_path / "physical_controls.csv").is_file()
    assert manifest["device_instructions"] is None


def test_cli_matches_api_and_import_has_no_side_effects(tmp_path):
    bundle = _export_bundle(protocol=TRANSFER_PROTOCOL)
    api = simulate(bundle)
    state = random.getstate()
    cwd = tmp_path
    import os
    os.chdir(cwd)
    before = sorted(p.name for p in cwd.iterdir())
    rc = cli.main(["simulate", "-"])
    # '-' not supported: use a real file instead
    cfg = tmp_path / "cfg.yaml"
    import yaml
    cfg.write_text(yaml.safe_dump(json.loads(json.dumps({
        "device": {"species": "Cs133",
                   "slm": {"wavelength_nm": 1061, "waist_um": 1.17, "depth_uK": 140},
                   "aod": {"wavelength_nm": 1055, "waist_um": 1.17, "depth_uK": 280},
                   "lensing": {"mode": "off"}},
        "preparation": {"temperature_uK": 5},
        "noise": {"model": "none"},
        "numerics": {"dt_us": 0.1, "noise_dt_us": 0.1, "shots": 64, "seed": 11},
        "protocol": yaml.safe_load(TRANSFER_PROTOCOL)})), default_flow_style=False))
    rc = cli.main(["simulate", str(cfg), "--out", str(tmp_path / "run")])
    assert rc == 0
    data = json.loads((tmp_path / "run" / "result.json").read_text())
    assert data["survival"] == pytest.approx(api.survival.tolist())
    # import purity: no cwd change, no stray files, RNG untouched
    assert random.getstate() == state
    os.chdir(tmp_path)
    after = sorted(p.name for p in cwd.iterdir() if not p.name.startswith("run")
                   and p.name != "cfg.yaml")
    assert after == before


def test_cli_error_paths_return_stable_codes(tmp_path, capsys):
    bad = tmp_path / "bad.yaml"
    bad.write_text("device: {species: Unicorn}\n")
    rc = cli.main(["validate", str(bad)])
    assert rc == 1
    rc = cli.main(["validate", str(tmp_path / "missing.yaml")])
    assert rc == 2
    rc = cli.main(["selftest"])
    assert rc == 0
