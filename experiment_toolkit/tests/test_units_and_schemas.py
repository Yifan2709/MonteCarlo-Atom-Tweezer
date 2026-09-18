"""Units, schema errors and contract tests (A03, A04)."""
import numpy as np
import pytest

from tweezer_experiment import units as U
from tweezer_experiment.errors import ToolkitError
from tweezer_experiment.schemas import ExperimentBundle, validate_experiment


def test_unit_anchors_against_independent_values():
    # 1 uK of energy must equal k_B * 1e-6 J exactly.
    assert U.temperature_uK_to_joule(1.0) == pytest.approx(1.380649e-29, rel=1e-12)
    assert U.joule_to_temperature_uK(U.temperature_uK_to_joule(140.0)) == pytest.approx(140.0)
    assert U.us_to_s(400.0) == pytest.approx(4.0e-4)
    assert U.um_to_m(2.4) == pytest.approx(2.4e-6)
    assert U.hz_to_rad_per_s(24611.0) == pytest.approx(2 * np.pi * 24611.0)
    assert U.rad_per_s_to_hz(U.hz_to_rad_per_s(1000.0)) == pytest.approx(1000.0)
    assert U.mass_u_to_kg(132.90545196) == pytest.approx(2.2069e-25, rel=1e-3)
    # arrays are supported (CSV columns)
    assert np.allclose(U.us_to_s(np.array([0.0, 10.0, 20.0])), [0.0, 1e-5, 2e-5])


def test_to_si_rejects_unknown_units():
    with pytest.raises(ToolkitError) as e:
        U.to_si(1.0, "furlongs", "length")
    assert e.value.code == "E_UNIT"


BASE = """
device:
  species: Cs133
  slm: {wavelength_nm: 1061, waist_um: 1.17, depth_uK: 140}
  aod: {wavelength_nm: 1055, waist_um: 1.17, depth_uK: 280}
  lensing: {mode: "off"}
preparation: {temperature_uK: 15}
noise: {model: none}
numerics: {dt_us: 0.1, noise_dt_us: 0.1, shots: 32, seed: 7}
protocol:
  segments:
    - id: hold
      type: static
      duration_us: 50
      aod_center_um: [0, 0]
      aod_depth_uK: 0
      slm_factor: 1.0
      observe: slm
  sequence: [hold]
"""


def _load(text=BASE):
    import yaml
    return ExperimentBundle.parse(yaml.safe_load(text))


def test_unknown_field_reports_path_and_code():
    import yaml
    bad = yaml.safe_load(BASE)
    bad["device"]["quantum_flux"] = 3
    with pytest.raises(ToolkitError) as e:
        ExperimentBundle.parse(bad)
    assert e.value.code == "E_UNKNOWN_FIELD"
    assert "quantum_flux" in e.value.field


def test_missing_required_section():
    import yaml
    bad = yaml.safe_load(BASE)
    del bad["noise"]
    with pytest.raises(ToolkitError) as e:
        ExperimentBundle.parse(bad)
    assert e.value.code == "E_MISSING_FIELD"


def test_nan_rejected_with_field_path():
    import yaml
    bad = yaml.safe_load(BASE)
    bad["preparation"]["temperature_uK"] = float("nan")
    with pytest.raises(ToolkitError) as e:
        ExperimentBundle.parse(bad)
    assert e.value.code == "E_NOT_FINITE"
    assert e.value.field == "preparation.temperature_uK"


def test_unknown_noise_model_is_model_error():
    import yaml
    bad = yaml.safe_load(BASE)
    bad["noise"]["model"] = "spin_boson"
    with pytest.raises(ToolkitError) as e:
        ExperimentBundle.parse(bad)
    assert e.value.code == "E_MODEL_UNKNOWN"


def test_non_divisible_duration_reports_nearest():
    import yaml
    bad = yaml.safe_load(BASE)
    bad["protocol"]["segments"][0]["duration_us"] = 50.03
    report = validate_experiment(bad)
    assert not report.ok
    assert any("E_NOT_DIVISIBLE" in err and "50.03" in err for err in report.errors)


def test_yaml_off_boolean_is_coerced():
    bundle = _load()
    assert bundle.device.lensing_mode == "off"


def test_noise_dt_smaller_than_dt_is_cross_field_error():
    import yaml
    bad = yaml.safe_load(BASE)
    bad["numerics"]["noise_dt_us"] = 0.01
    report = validate_experiment(bad)
    assert not report.ok
    assert any("noise_dt" in err for err in report.errors)


def test_none_model_with_amplitudes_is_conflict():
    import yaml
    bad = yaml.safe_load(BASE)
    bad["noise"]["aod_rin_psd_per_hz"] = 1e-8
    report = validate_experiment(bad)
    assert not report.ok
    assert any("model=none" in err for err in report.errors)


def test_effective_config_shows_sources_and_overrides():
    import yaml
    raw = yaml.safe_load(BASE)
    raw["device"]["slm"]["depth_uK__source"] = {"tag": "measured", "note": "lab 2026"}
    raw["preparation"]["temperature_uK__source"] = {"tag": "assumed"}
    bundle = ExperimentBundle.parse(raw)
    report = validate_experiment(bundle)
    assert report.ok
    cfg = report.effective_config
    assert cfg["device"]["slm"]["nominal_depth_j"] == pytest.approx(140 * 1.380649e-29)
    assert cfg["sources"]["device.slm.depth_uK"]["tag"] == "measured"
    assert report.capability["source_tag_counts"]["measured"] == 1
    # missing-measurement list reacts to the tag
    assert "device.slm.depth_uK" not in report.capability["missing_measurements"]


def test_parameter_reaches_physics_core():
    """A04: changing wavelength/waist/duration actually changes the numbers."""
    from tweezer_experiment.physics.field import BeamContext, beam_field
    b1 = _load()
    ctx1 = BeamContext.from_device(b1.device)
    # trap curvatures: fx slope at origin = 4D/(w^2) * ... via finite diff
    d = 140 * 1.380649e-29
    eps = 1e-9
    f_plus = beam_field(eps, 0, 0, d, ctx1.slm_waist_m, ctx1.slm_zr_m)[1]
    f_minus = beam_field(-eps, 0, 0, d, ctx1.slm_waist_m, ctx1.slm_zr_m)[1]
    k1 = (f_minus - f_plus) / (2 * eps)
    import yaml
    b2 = ExperimentBundle.parse(yaml.safe_load(
        BASE.replace("waist_um: 1.17", "waist_um: 1.0", 1)))
    ctx2 = BeamContext.from_device(b2.device)
    f_plus = beam_field(eps, 0, 0, d, ctx2.slm_waist_m, ctx2.slm_zr_m)[1]
    f_minus = beam_field(-eps, 0, 0, d, ctx2.slm_waist_m, ctx2.slm_zr_m)[1]
    k2 = (f_minus - f_plus) / (2 * eps)
    assert k2 > k1 * 1.3  # tighter waist -> stiffer trap, actually computed
    # duration change reaches the compiled grid
    from tweezer_experiment.protocol_compile import compile_protocol
    c1 = compile_protocol(b1)
    b3 = ExperimentBundle.parse(yaml.safe_load(BASE.replace("duration_us: 50", "duration_us: 100")))
    c3 = compile_protocol(b3)
    assert c3.duration_s == pytest.approx(2 * c1.duration_s)
    assert len(c3.controls) > len(c1.controls)


def test_search_grid_bounded_and_duplicate_targets_rejected():
    import yaml
    raw = yaml.safe_load(BASE)
    raw["protocol"]["segments"] = [
        {"id": "mv", "type": "transfer", "duration_us": 60, "from_um": [0, 0],
         "to_um": [2.0, 0], "depth_from_uK": 0, "depth_to_uK": 280,
         "path": "paper_manual", "ramp_fraction": 0.48, "slm_factor": 1.0,
         "observe": "aod"}]
    raw["protocol"]["sequence"] = ["mv"]
    raw["search"] = {
        "variables": [
            {"target": "segment.mv.duration_us", "grid": [40, 60, 8]},
            {"target": "segment.mv.ramp_fraction", "values": [0.4, 0.48]},
        ],
        "objective": {"type": "maximize_capture"},
        "max_candidates": 2,
    }
    with pytest.raises(ToolkitError) as e:
        ExperimentBundle.parse(raw)
    assert e.value.code == "E_BAD_VALUE"          # 16 candidates > max 2
    raw["search"]["max_candidates"] = 64
    raw["search"]["variables"].append({"target": "segment.mv.duration_us", "values": [50]})
    with pytest.raises(ToolkitError) as e:
        ExperimentBundle.parse(raw)
    assert e.value.code == "E_CONFLICT"           # duplicate target
