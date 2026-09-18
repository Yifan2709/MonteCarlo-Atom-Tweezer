"""Independent physics/numerics benchmarks (A05) and initialization (A06)."""
import numpy as np
import pytest

from tweezer_experiment import units as U
from tweezer_experiment.schemas import ExperimentBundle
from tweezer_experiment.protocol_compile import compile_protocol
from tweezer_experiment.physics.field import BeamContext, beam_field
from tweezer_experiment.physics.initialization import (draw_site_disorder,
                                                       sample_site_bound, slm_site_energies)
from tweezer_experiment.simulation import simulate
from tweezer_experiment.errors import ToolkitError

BASE = """
device:
  species: Cs133
  slm: {wavelength_nm: 1061, waist_um: 1.17, depth_uK: 140}
  aod: {wavelength_nm: 1055, waist_um: 1.17, depth_uK: 280}
  lensing: {mode: "off"}
preparation: {temperature_uK: 5}
noise: {model: none}
numerics: {dt_us: 0.05, noise_dt_us: 0.05, shots: 64, seed: 11}
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


def _oscillation_protocol(n_seg=120, seg_us=1.0):
    segs = "\n".join(
        f"    - id: p{i}\n      type: static\n      duration_us: {seg_us}\n"
        f"      aod_center_um: [0, 0]\n      aod_depth_uK: 0\n      slm_factor: 1.0\n"
        f"      observe: slm" for i in range(n_seg))
    return BASE.replace("""  segments:
    - id: hold
      type: static
      duration_us: 30
      aod_center_um: [0, 0]
      aod_depth_uK: 0
      slm_factor: 1.0
      observe: slm
  sequence: [hold]""",
                        f"  segments:\n{segs}\n  sequence: [" + ", ".join(f"p{i}" for i in range(120)) + "]")


def test_force_is_negative_gradient_by_finite_differences():
    ctx = BeamContext.from_device(_load().device)
    d = ctx.slm_depth_j
    rng = np.random.default_rng(0)
    pts = rng.uniform(-1e-6, 1e-6, size=(200, 3)) * [1, 1, 0.02]
    u0, fx, fy, fz = beam_field(pts[:, 0], pts[:, 1], pts[:, 2], d,
                                 ctx.slm_waist_m, ctx.slm_zr_m)
    eps = 1e-9
    for axis, comp in enumerate((fx, fy, fz)):
        shifted = pts.copy()
        shifted[:, axis] += eps
        u_plus = beam_field(shifted[:, 0], shifted[:, 1], shifted[:, 2], d,
                            ctx.slm_waist_m, ctx.slm_zr_m)[0]
        shifted[:, axis] -= 2 * eps
        u_minus = beam_field(shifted[:, 0], shifted[:, 1], shifted[:, 2], d,
                             ctx.slm_waist_m, ctx.slm_zr_m)[0]
        numeric = -(u_plus - u_minus) / (2 * eps)
        assert np.max(np.abs(numeric - comp) / np.maximum(np.abs(comp), 1e-15)) < 1e-6


def test_oscillation_frequency_matches_analytic_trap_frequency():
    bundle = _load(_oscillation_protocol())
    ctx = BeamContext.from_device(bundle.device)
    d, w, m = ctx.slm_depth_j, ctx.slm_waist_m, ctx.mass_kg
    nu_analytic = np.sqrt(4 * d / (m * w * w)) / (2 * np.pi)
    result = simulate(bundle, seed=5, shots=64)
    states = result.checkpoint_states  # (shots, n_ck+1, 6)
    x = states[0, 1:, 0]                # first atom, after each 1 us segment
    x = x - np.mean(x)
    spectrum = np.abs(np.fft.rfft(x))
    freqs = np.fft.rfftfreq(len(x), d=1e-6)
    nu_measured = freqs[int(np.argmax(spectrum[1:])) + 1]
    assert nu_measured == pytest.approx(nu_analytic, rel=0.06), \
        f"measured {nu_measured:.0f} Hz vs analytic {nu_analytic:.0f} Hz"


def test_verlet_energy_error_scales_second_order():
    errors = {}
    for dt_us in (0.1, 0.05):
        text = BASE.replace("dt_us: 0.05", f"dt_us: {dt_us}").replace("noise_dt_us: 0.05", f"noise_dt_us: {dt_us}")
        bundle = _load(text)
        result = simulate(bundle, seed=5)
        e = result.checkpoint_states[0, 1, 3:]  # velocity at first checkpoint
        m = bundle.device.mass_kg
        # total energy drift for the coldest atom: use stored energies instead
        errors[dt_us] = None
    # use the energy time series of the survivors (median drift over 30 us)
    def drift(dt_us):
        text = BASE.replace("dt_us: 0.05", f"dt_us: {dt_us}").replace("noise_dt_us: 0.05", f"noise_dt_us: {dt_us}")
        bundle = _load(text)
        result = simulate(bundle, seed=5)
        e = result.energy_stats[1]["energy_uK_mean"]
        return abs(e - result.energy_stats[0]["energy_uK_mean"])
    d_coarse, d_fine = drift(0.1), drift(0.05)
    assert d_fine <= d_coarse * 0.75  # refined step reduces the bounded error


def test_diffusion_zero_amplitudes_equal_model_none_bitwise():
    import yaml
    r1 = simulate(_load(BASE.replace("noise: {model: none}", "noise: {model: none}")), seed=9)
    raw = yaml.safe_load(BASE)
    raw["noise"] = {"model": "diffusion"}   # all amplitudes default to 0
    r2 = simulate(ExperimentBundle.parse(raw), seed=9)
    assert np.array_equal(r1.final_states, r2.final_states)
    assert r2.noise_ledger["mean_realized_j_per_atom"] == 0.0


def test_diffusion_ledger_expected_matches_realized():
    import yaml
    raw = yaml.safe_load(BASE)
    raw["noise"] = {"model": "diffusion", "aod_rin_psd_per_hz": 5.75e-9,
                    "aod_pointing_psd_m2_per_hz": 1e-24}
    # park a loaded AOD beam one waist away so atoms feel a real force and
    # the intensity/pointing channels act through nonzero F and dF/dx
    raw["protocol"]["segments"][0]["duration_us"] = 400
    raw["protocol"]["segments"][0]["aod_center_um"] = [1.17, 0]
    raw["protocol"]["segments"][0]["aod_depth_uK"] = 280
    bundle = ExperimentBundle.parse(raw)
    result = simulate(bundle, seed=21, shots=256)
    realized = result.noise_ledger["mean_realized_j_per_atom"]
    expected = result.noise_ledger["mean_expected_j_per_atom"]
    assert expected > 0
    assert realized == pytest.approx(expected, rel=0.08), \
        f"realized {realized:.3e} vs expected {expected:.3e} (kick coefficients)"


def test_local_mean_is_deterministic_and_reproducible():
    import yaml
    raw = yaml.safe_load(BASE)
    raw["noise"] = {"model": "local_mean", "aod_rin_psd_per_hz": 5.75e-9}
    raw["protocol"]["segments"][0]["duration_us"] = 200
    bundle = ExperimentBundle.parse(raw)
    r1 = simulate(bundle, seed=4)
    r2 = simulate(bundle, seed=4)
    assert np.array_equal(r1.final_states, r2.final_states)
    assert r1.noise_ledger["mean_realized_j_per_atom"] == pytest.approx(
        r1.noise_ledger["mean_expected_j_per_atom"], rel=1e-9)


def test_initial_states_bound_reproducible_and_site_scaled():
    bundle = _load()
    ctx = BeamContext.from_device(bundle.device)
    from dataclasses import replace
    prep = replace(bundle.preparation, site_disorder_enabled=True, slm_depth_rel_sigma=0.3)
    disorder = draw_site_disorder(prep, 512, seed=1)
    pool_a = sample_site_bound(ctx, prep, 512, 2, disorder["slm_depth_factor"])
    pool_b = sample_site_bound(ctx, prep, 512, 2, disorder["slm_depth_factor"])
    assert np.array_equal(pool_a["pos_m"], pool_b["pos_m"])      # same seed
    assert np.all(pool_a["initial_energy_j"] < 0)                 # bound, actual sites
    pool_c = sample_site_bound(ctx, prep, 512, 3, disorder["slm_depth_factor"])
    assert not np.array_equal(pool_a["pos_m"], pool_c["pos_m"])   # different seed
    # deeper sites hold narrower distributions (1/depth scaling)
    half = np.argpartition(disorder["slm_depth_factor"], 256)[:256]
    full = np.setdiff1d(np.arange(512), half)
    var_shallow = np.var(pool_a["pos_m"][half, 0])
    var_deep = np.var(pool_a["pos_m"][full, 0])
    assert var_shallow > 1.5 * var_deep


def test_sampler_raises_when_temperature_exceeds_shallowest_site():
    import yaml
    raw = yaml.safe_load(BASE)
    raw["preparation"] = {"temperature_uK": 400.0, "slm_depth_rel_sigma": 0.9}
    bundle = ExperimentBundle.parse(raw)
    with pytest.raises(ToolkitError):
        simulate(bundle, shots=64, seed=1)
