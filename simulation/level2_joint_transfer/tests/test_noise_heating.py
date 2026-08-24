"""平均强度确定性噪声加热的测试：机制、账本、注入接口与配置。"""
from __future__ import annotations

import math
import pickle
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import yaml

from level0_static_trap.constants import microkelvin_to_joule
from level0_static_trap.potentials import gaussian_force, gaussian_potential
from level2_joint_transfer.config import load_config
from level2_joint_transfer.level2_simulation import (evaluate_waveform_on_states,
                                                     evaluate_waveform_parallel,
                                                     physics_from_config)
from level2_joint_transfer.noise_heating import (MeanHeating, build_mean_heating,
                                                 heating_from_config,
                                                 parametric_rate_from_rin,
                                                 pointing_power_from_psd,
                                                 recoil_energy_j,
                                                 recoil_scattering_rate_cs,
                                                 scaled_heating,
                                                 velocity_verlet_time_dependent_heating)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG = PROJECT_ROOT / "configs" / "level2_joint_transfer.yaml"
UK_PER_JOULE = 1.0 / microkelvin_to_joule(1.0)


@pytest.fixture(scope="module")
def cfg():
    return load_config(CONFIG)


@pytest.fixture(scope="module")
def physics(cfg):
    return physics_from_config(cfg)


# ------------------------------------------------ 导出函数与物理常数
def test_recoil_energy_and_scattering_rate_cs_1061nm(cfg):
    mass = cfg.atom.mass_kg
    e_r = recoil_energy_j(1061e-9, mass)
    assert abs(e_r / 1.380649e-23 * 1e9 - 64.0) < 3.0  # E_r/k_B ≈ 64 nK
    rate = recoil_scattering_rate_cs(1061e-9, cfg.aod.initial_depth_j)
    assert 2.0 < rate < 4.0  # 280 μK 两能级近似 ≈ 2.8 s^-1


def test_parametric_and_pointing_formulas(cfg):
    nu0 = 25.46e3
    rate = parametric_rate_from_rin(nu0, 1.0e-8)
    assert abs(rate / (0.25 * math.pi ** 2 * nu0 ** 2 * 1.0e-8) - 1.0) < 1e-12
    assert 14.0 < rate < 18.0  # (π²/4)ν₀²S_RIN ≈ 16 s^-1
    omega0 = 2.0 * math.pi * nu0
    power = pointing_power_from_psd(cfg.atom.mass_kg, nu0, 1.0e-18)
    analytic = cfg.atom.mass_kg * omega0 ** 4 * 1.0e-18 / 8.0
    assert abs(power / analytic - 1.0) < 1e-12
    assert abs(power * UK_PER_JOULE - 1.31e6) / 1.31e6 < 0.05  # 1 nm/√Hz 共振 ≈ 1.3 K/s


# ------------------------------------------------ 机制：常数功率线性增长
def test_constant_power_grows_energy_linearly(physics):
    mass = physics["mass_kg"]
    depth, waist = physics["slm_depth_j"], physics["slm_waist_m"]
    dt, total = 0.05e-6, 500e-6
    times = np.arange(round(total / dt) + 1) * dt
    power_w = microkelvin_to_joule(100.0)  # 100 μK/s
    heating = build_mean_heating(recoil_power_w=power_w)
    force = lambda x, t: gaussian_force(x, depth, waist, 0.0)
    e_osc = lambda x, v, t: max(0.0, 0.5 * mass * v * v
                                + float(gaussian_potential(x, depth, waist, 0.0)) + depth)
    _t, x, v, _a, work = velocity_verlet_time_dependent_heating(
        force, mass, 0.2e-6, 0.0, times, heating, e_osc)
    energy = 0.5 * mass * v ** 2 + gaussian_potential(x, depth, waist, 0.0)
    expected = microkelvin_to_joule(100.0) * times
    assert np.max(np.abs(energy - energy[0] - expected)) / depth < 1e-6
    assert abs(work["recoil_J"][-1] - expected[-1]) / expected[-1] < 1e-12
    assert np.allclose(work["total_J"], work["recoil_J"], rtol=1e-15)


# ------------------------------------------------ 机制：参量通道指数增长
def test_parametric_channel_grows_oscillation_energy_exponentially(physics):
    mass = physics["mass_kg"]
    depth, waist = physics["slm_depth_j"], physics["slm_waist_m"]
    dt, total = 0.05e-6, 500e-6
    times = np.arange(round(total / dt) + 1) * dt
    rate = 1000.0
    heating = build_mean_heating(parametric_rate_s=rate)
    force = lambda x, t: gaussian_force(x, depth, waist, 0.0)

    def e_osc(x, v, t):
        return max(0.0, 0.5 * mass * v * v + float(gaussian_potential(x, depth, waist, 0.0)) + depth)

    _t, x, v, _a, work = velocity_verlet_time_dependent_heating(
        force, mass, 0.3e-6, 0.0, times, heating, e_osc)
    e_osc_series = 0.5 * mass * v ** 2 + gaussian_potential(x, depth, waist, 0.0) + depth
    growth = e_osc_series[-1] / e_osc_series[0]
    assert abs(growth / math.exp(rate * total) - 1.0) < 1e-3
    assert np.allclose(work["parametric_J"], work["total_J"], rtol=1e-15)


# ------------------------------------------------ 机制：自定义注入函数
def test_custom_time_dependent_callables_are_respected(physics):
    mass = physics["mass_kg"]
    depth, waist = physics["slm_depth_j"], physics["slm_waist_m"]
    dt, total = 0.05e-6, 200e-6
    times = np.arange(round(total / dt) + 1) * dt
    p0 = microkelvin_to_joule(50.0)
    modulation = 2.0 * math.pi / total

    def timed_power(x, v, t, e_osc):
        return p0 * (1.0 + 0.5 * math.sin(modulation * t))

    heating = MeanHeating(recoil=timed_power)
    force = lambda x, t: gaussian_force(x, depth, waist, 0.0)
    e_osc = lambda x, v, t: 0.0
    _t, x, v, _a, work = velocity_verlet_time_dependent_heating(
        force, mass, 0.0, 0.0, times, heating, e_osc)
    expected = p0 * (times[-1] + 0.5 * (1.0 - math.cos(modulation * times[-1])) / modulation)
    assert abs(work["recoil_J"][-1] - expected) / expected < 1e-12


# ------------------------------------------------ 全零通道 = 与关闭时逐位一致
def test_zero_heating_channels_are_bit_identical_to_disabled(cfg, physics):
    states = {"x0_m": np.array([2.2e-6, 2.5e-6]), "v0_m_per_s": np.array([0.01, -0.02]),
              "initial_aod_energy_J": np.array([-1.0e-28, -1.1e-28]),
              "initial_excitation_J": np.array([2e-28, 3e-28]), "shots": 2}
    spec = {"type": "sequential", "duration_s": 20e-6, "move_fraction": 0.52}
    off = evaluate_waveform_on_states(dict(spec), states, physics, 0.05e-6, 5e-6)
    zero = evaluate_waveform_on_states(dict(spec), states, physics, 0.05e-6, 5e-6,
                                       heating=MeanHeating())
    for a, b in zip(off["records"], zero["records"]):
        for key in ("final_x_um", "final_v_m_s", "final_slm_energy_uK", "captured"):
            assert a[key] == b[key]


# ------------------------------------------------ 功-能账本含噪声通道
def test_work_energy_residual_absorbs_noise_channel(cfg, physics):
    states = {"x0_m": np.array([2.3e-6, 2.45e-6]), "v0_m_per_s": np.array([0.0, 0.03]),
              "initial_aod_energy_J": np.array([-1.0e-28, -1.1e-28]),
              "initial_excitation_J": np.array([2e-28, 3e-28]), "shots": 2}
    spec = dict({"type": "sequential", "duration_s": 20e-6, "move_fraction": 0.52})
    heating = build_mean_heating(recoil_power_w=microkelvin_to_joule(2000.0),
                                 parametric_rate_s=500.0,
                                 pointing_power_w=microkelvin_to_joule(1000.0))
    result = evaluate_waveform_on_states(spec, states, physics, 0.05e-6, 5e-6,
                                         work_energy=True, heating=heating)
    depth = physics["slm_depth_j"]
    for record in result["records"]:
        assert record["noise_total_energy_uK"] > 0
        assert (record["noise_recoil_energy_uK"] > 0
                and record["noise_parametric_energy_uK"] > 0
                and record["noise_pointing_energy_uK"] > 0)
        assert record["max_abs_work_energy_residual_over_depth"] < 1e-3
    base = evaluate_waveform_on_states(spec, states, physics, 0.05e-6, 5e-6)
    for noisy, clean in zip(result["records"], base["records"]):
        assert noisy["final_slm_energy_uK"] > clean["final_slm_energy_uK"]  # 加热使末能量升高


# ------------------------------------------------ 随时间越发不稳定
def test_heating_grows_with_duration_and_destabilizes_capture(cfg, physics):
    # 注意：记录中 *_uK 字段沿用既有约定，数值单位实为开尔文（见 README）。
    states = {"x0_m": np.array([2.3e-6, 2.5e-6]), "v0_m_per_s": np.array([0.0, 0.01]),
              "initial_aod_energy_J": np.array([-1.0e-28, -1.1e-28]),
              "initial_excitation_J": np.array([2e-28, 3e-28]), "shots": 2}
    heating = build_mean_heating(parametric_rate_s=20000.0)
    injected = []
    for duration_us in (20.0, 60.0):
        spec = {"type": "sequential", "duration_s": duration_us * 1e-6, "move_fraction": 0.52}
        result = evaluate_waveform_on_states(spec, states, physics, 0.05e-6, 5e-6, heating=heating)
        injected.append(np.mean([r["noise_total_energy_uK"] for r in result["records"]]))
    assert injected[1] > 3.0 * injected[0]  # 指数型累积（纯指数预期 ~4.7×），时间越长越不稳定
    assert injected[1] > 1.0e-4  # 60 μs 末注入 > 100 μK（单位：K，即 1e-4）


# ------------------------------------------------ 配置路径
def test_heating_from_config_disabled_and_enabled(cfg):
    assert heating_from_config(cfg) is None
    enabled = replace(cfg, noise_mean_heating=replace(cfg.noise_mean_heating, enabled=True))
    heating = heating_from_config(enabled)
    assert heating is not None
    described = heating.describe()
    assert described["model"] == "mean_rate_deterministic"
    assert described["provenance"] == "assumed_sensitivity_only"
    derived = described["channels"]
    assert 0.2 < derived["recoil"]["power_uK_per_s"] < 0.6  # 2×64nK×2.8/s ≈ 0.36 μK/s
    assert 12.0 < derived["parametric"]["rate_per_s"] < 20.0
    # 1 pm/√Hz @25.5 kHz → m·ω₀⁴·S_x/8 ≈ 1.3 μK/s
    assert 0.8 < derived["pointing"]["power_uK_per_s"] < 2.0
    # 显式速率优先于推导
    explicit = replace(cfg, noise_mean_heating=replace(cfg.noise_mean_heating, enabled=True,
                                                       parametric_rate_s=123.0))
    assert heating_from_config(explicit).describe()["channels"]["parametric"]["rate_per_s"] == 123.0


def test_scaled_heating_multiplies_all_channels():
    base = build_mean_heating(recoil_power_w=1.0, parametric_rate_s=10.0, pointing_power_w=2.0)
    scaled = scaled_heating(base, 100.0)
    powers = scaled.channel_powers(0.0, 0.0, 0.0, 5.0)
    assert powers["recoil"] == 100.0
    assert powers["parametric"] == 100.0 * 10.0 * 5.0
    assert powers["pointing"] == 200.0
    clone = pickle.loads(pickle.dumps(scaled))  # 敏感性扫描需跨进程
    assert clone.channel_powers(0.0, 0.0, 0.0, 5.0) == powers


def test_negative_noise_parameter_is_rejected(cfg, tmp_path):
    raw = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    raw["noise_mean_heating"]["parametric_rate_s"] = -1.0
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="必须非负"):
        load_config(path)


# ------------------------------------------------ 并行与 pickle
def test_model_pickle_roundtrip():
    heating = build_mean_heating(recoil_power_w=1e-27, parametric_rate_s=5.0, pointing_power_w=2e-27)
    clone = pickle.loads(pickle.dumps(heating))
    assert clone.channel_powers(0.0, 0.0, 0.0, 3e-28) == heating.channel_powers(0.0, 0.0, 0.0, 3e-28)


def test_parallel_matches_serial_with_heating(cfg, physics):
    from level2_joint_transfer.level2_simulation import sample_initial_states
    states = sample_initial_states(cfg, 90210, 8)
    spec = {"type": "sequential", "duration_s": 20e-6, "move_fraction": 0.52}
    heating = build_mean_heating(recoil_power_w=microkelvin_to_joule(500.0),
                                 parametric_rate_s=300.0)
    serial = evaluate_waveform_on_states(dict(spec), states, physics, 0.05e-6, 5e-6, heating=heating)
    parallel = evaluate_waveform_parallel(dict(spec), states, physics, 0.05e-6, 5e-6, heating=heating)
    for a, b in zip(serial["records"], parallel["records"]):
        assert a["final_x_um"] == b["final_x_um"]
        assert a["noise_total_energy_uK"] == b["noise_total_energy_uK"]
        assert a["captured"] == b["captured"]
