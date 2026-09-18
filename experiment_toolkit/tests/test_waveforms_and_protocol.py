"""Waveform families and protocol compilation tests (F2, A08)."""
import numpy as np
import pytest

from tweezer_experiment.waveforms import FAMILIES, transfer_controls
from tweezer_experiment.errors import ToolkitError
from tweezer_experiment.schemas import ExperimentBundle
from tweezer_experiment.protocol_compile import compile_protocol, apply_override
from tweezer_experiment import units as U

BASE = """
device:
  species: Cs133
  slm: {wavelength_nm: 1061, waist_um: 1.17, depth_uK: 140}
  aod: {wavelength_nm: 1055, waist_um: 1.17, depth_uK: 280}
  lensing: {mode: "off"}
preparation: {temperature_uK: 5}
noise: {model: none}
numerics: {dt_us: 0.1, noise_dt_us: 0.1, shots: 16, seed: 3}
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


@pytest.mark.parametrize("name", sorted(FAMILIES))
def test_family_endpoints_and_derivative_consistency(name):
    fn = FAMILIES[name]
    u = np.linspace(0.0, 1.0, 2001)
    q, dq = fn(u)
    assert q[0] == pytest.approx(0.0, abs=1e-12)
    assert q[-1] == pytest.approx(1.0, abs=1e-9)
    # analytic derivative must match finite differences of the position
    fd = np.gradient(q, u)
    assert np.max(np.abs(fd - dq)) < 2e-3


def test_constant_jerk_zero_endpoint_acceleration():
    fn = FAMILIES["constant_jerk"]
    u = np.linspace(0, 1, 4001)
    q, dq = fn(u)
    accel = np.gradient(dq, u)
    assert abs(accel[0]) < 5e-3 and abs(accel[-1]) < 5e-3
    assert dq.max() == pytest.approx(2.0, abs=1e-9)   # peak = 2 * average
    assert q[2000] == pytest.approx(0.5, abs=1e-9)


def test_paper_manual_has_endpoint_acceleration_jump_and_it_is_reported():
    bundle = _load(BASE.replace(
        """  segments:
    - id: hold
      type: static
      duration_us: 30
      aod_center_um: [0, 0]
      aod_depth_uK: 0
      slm_factor: 1.0
      observe: slm
  sequence: [hold]""",
        """  segments:
    - id: move
      type: transfer
      duration_us: 100
      from_um: [0, 0]
      to_um: [2.4, 0]
      depth_from_uK: 0
      depth_to_uK: 280
      path: paper_manual
      ramp_fraction: 0.48
  sequence: [move]"""))
    compiled = compile_protocol(bundle)
    notes = [a for a in compiled.audit["accel_notes"] if a["segment"] == "move"]
    assert notes, "paper_manual acceleration jump must be disclosed in the audit"
    # the jump itself: velocity derivative discontinuity at move start
    c = compiled.controls
    ramp_end_idx = int(round(48 / 0.1))
    v_before = (c[ramp_end_idx, 1] - c[ramp_end_idx - 1, 1]) / (c[ramp_end_idx, 0] - c[ramp_end_idx - 1, 0])
    v_after = (c[ramp_end_idx + 1, 1] - c[ramp_end_idx, 1]) / (c[ramp_end_idx + 1, 0] - c[ramp_end_idx, 0])
    assert v_before == pytest.approx(0.0, abs=1e-12)
    assert v_after > 0.0  # acceleration starts at 6*d/tau^2, not zero


def test_dropoff_is_exact_time_mirror_of_pickup():
    dur, rf, dist = 400e-6, 0.48, 2.4e-6
    d0, d1 = 0.0, 280 * 1.380649e-29
    pickup = transfer_controls("paper_manual", rf, dur, (0, 0), (dist, 0), d0, d1,
                               order="depth_then_move")
    drop = transfer_controls("paper_manual", rf, dur, (dist, 0), (0, 0), d1, d0,
                             order="move_then_depth")
    t = np.linspace(0, dur, 800)
    x1, _, _, _, dep1 = pickup(dur - t)   # dropoff = exact time reverse of pickup
    x2, _, _, _, dep2 = drop(t)
    assert np.allclose(x2, x1, atol=1e-15, rtol=0)          # position mirror
    assert np.allclose(dep2, dep1, atol=1e-30, rtol=0)      # depth mirror


def test_checkpoint_count_follows_protocol_definition():
    text = BASE.replace("sequence: [hold]", """sequence:
    - repeat: 13
      items: [hold]""")
    bundle = _load(text)
    compiled = compile_protocol(bundle)
    assert compiled.n_checkpoints == 13          # not 2*rounds+1
    assert len(compiled.segment_names) == 13


def test_final_hold_appends_checkpoint():
    text = BASE.replace("sequence: [hold]", "sequence: [hold]\n  final_hold_us: 50")
    compiled = compile_protocol(_load(text))
    assert compiled.n_checkpoints == 2
    assert compiled.checkpoints[-1][0] == "after_final_hold"


def test_boundary_switch_events_are_recorded_not_smoothed():
    text = """device:
  species: Cs133
  slm: {wavelength_nm: 1061, waist_um: 1.17, depth_uK: 140}
  aod: {wavelength_nm: 1055, waist_um: 1.17, depth_uK: 280}
  lensing: {mode: "off"}
preparation: {temperature_uK: 5}
noise: {model: none}
numerics: {dt_us: 0.1, noise_dt_us: 0.1, shots: 16, seed: 3}
protocol:
  segments:
    - id: slm_on
      type: static
      duration_us: 20
      aod_center_um: [0, 0]
      aod_depth_uK: 0
      slm_factor: 1.0
    - id: slm_off
      type: static
      duration_us: 20
      aod_center_um: [0, 0]
      aod_depth_uK: 0
      slm_factor: 0.0
  sequence: [slm_on, slm_off]
"""
    compiled = compile_protocol(_load(text))
    slm_events = [e for e in compiled.events if e["field"] == "slm_factor"]
    assert len(slm_events) == 1
    assert slm_events[0]["left"] == pytest.approx(1.0)
    assert slm_events[0]["right"] == pytest.approx(0.0)
    # duplicated boundary node is preserved in the table
    t = compiled.controls[:, 0]
    dup = np.where(np.diff(t) == 0)[0]
    assert len(dup) == 1


def test_non_divisible_duration_error():
    text = BASE.replace("duration_us: 30", "duration_us: 30.07")
    report_text = text
    from tweezer_experiment.schemas import validate_experiment
    report = validate_experiment(_load(report_text) and __import__("yaml").safe_load(report_text))
    assert not report.ok
    assert any("E_NOT_DIVISIBLE" in e for e in report.errors)


def test_apply_override_changes_compiled_duration_only():
    bundle = _load(BASE)
    spec2 = apply_override(bundle.protocol, "segment.hold.duration_us", 60)
    c1 = compile_protocol(bundle)
    c2 = compile_protocol(bundle, spec2)
    assert c2.duration_s == pytest.approx(2 * c1.duration_s)
    with pytest.raises(ToolkitError):
        apply_override(bundle.protocol, "segment.nosuch.duration_us", 60)
