import numpy as np
import pytest
import yaml

from continuous_transfer.cli import ROOT, case_inputs
from continuous_transfer.initialization import sample_site_bound_harmonic, slm_energies


def settings():
    return yaml.safe_load((ROOT / "configs/continuous_level2c_345.yaml").read_text(encoding="utf-8"))


@pytest.mark.parametrize("arm", ["baseline_1d", "matched"])
def test_every_prepared_atom_is_bound_in_its_actual_site(arm):
    _, physics, pos, vel, draws, _, sampler = case_inputs(settings(), 40, 4096, 23003, arm)
    assert sampler.startswith("site_bound_harmonic_")
    # Check independently against the scalar Gaussian kernel used by dynamics.
    from continuous_transfer.engine import field
    w = physics["slm_waist_m"]
    zr = np.pi * w * w / 1061e-9
    potential = np.array([field(*p, physics["slm_depth_j"] * d, w, zr, 0, 0)[0]
                          for p, d in zip(pos, draws.slm_depth_factor)])
    energy = .5 * physics["mass_kg"] * (vel ** 2).sum(axis=1) + potential
    assert np.all(energy < 0)
    if arm == "baseline_1d":
        assert np.all(pos[:, 1:] == 0)
        assert np.all(vel[:, 1:] == 0)


def test_nominal_depth_legacy_exposes_the_archived_initialization_failure():
    _, physics, pos, vel, draws, _, _ = case_inputs(
        settings() | {"initial_sampling": "legacy_nominal"}, 25, 4096, 23003, "matched")
    assert np.sum(slm_energies(pos, vel, physics, draws.slm_depth_factor) >= 0) == 50


def test_local_harmonic_widths_and_site_identity_survive_rejection():
    _, physics, *_ = case_inputs(settings(), 5, 2, 23003, "matched")
    factors = np.r_[np.full(12000, .5), np.full(12000, 1.5)]
    low = sample_site_bound_harmonic(91, 3, physics, factors)
    positions = low["pos_m"]
    # Harmonic position variance scales inversely with local depth; velocities do not.
    ratio = positions[:12000, 0].var() / positions[12000:, 0].var()
    assert ratio == pytest.approx(3, rel=.06)
    assert low["vel_m_per_s"][:12000, 0].var() / low["vel_m_per_s"][12000:, 0].var() == pytest.approx(1, rel=.06)
    hot = sample_site_bound_harmonic(91, 80, physics, factors[:200])
    assert len(hot["pos_m"]) == 200
    assert np.all(hot["initial_energy_j"] < 0)
    assert np.any(hot["attempts_per_site"] > 1)
    repeated = sample_site_bound_harmonic(91, 80, physics, factors[:200])
    np.testing.assert_array_equal(hot["pos_m"], repeated["pos_m"])
    np.testing.assert_array_equal(hot["vel_m_per_s"], repeated["vel_m_per_s"])


def test_initial_temperature_and_heating_reference_can_be_varied_independently():
    fixed_heat = settings() | {"parametric_reference_temperature_uK": 25}
    a = case_inputs(fixed_heat, 5, 256, 23003, "matched")
    b = case_inputs(fixed_heat, 25, 256, 23003, "matched")
    assert a[5] == b[5]
    assert not np.array_equal(a[2], b[2])
    c = case_inputs(settings() | {"parametric_reference_temperature_uK": 50}, 5, 256, 23003, "matched")
    np.testing.assert_array_equal(a[2], c[2])
    np.testing.assert_array_equal(a[3], c[3])
    assert c[5]["parametric"] == pytest.approx(2 * a[5]["parametric"])
    for name in ["recoil", "pointing"]:
        assert c[5][name] == a[5][name]
