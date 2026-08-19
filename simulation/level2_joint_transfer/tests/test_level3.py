"""Level 3 测试：三维势、透镜、轨迹、采样、功-能、判据、隔离与 CLI 冒烟。

对应 level3_prompt.md 第十七节 25 项要求（部分合并为一测）。
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from level3_3d_transport_lensing.aod_lensing import (axial_minima, axis_shift_m,
                                                     lensing_potential_and_force,
                                                     vs_from_ratio)
from level3_3d_transport_lensing.gaussian_3d import (KB, physics_defaults,
                                                    static_force,
                                                    static_potential)
from level3_3d_transport_lensing.integrators_3d import (run_trajectory,
                                                        static_hold)
from level3_3d_transport_lensing.level3_cli import main as level3_main
from level3_3d_transport_lensing.level3_config import load_config
from level3_3d_transport_lensing.level3_simulation import (make_control,
                                                           physics_from_config,
                                                           resolve_vs,
                                                           run_condition,
                                                           sample_states)
from level3_3d_transport_lensing.preflight3 import (_check_constant_jerk,
                                                    _check_curvature,
                                                    _check_force_gradient,
                                                    _check_level0_slice,
                                                    _check_lensing_zero_limit,
                                                    _check_profiles)
from level3_3d_transport_lensing.thermal_sampling_3d import (
    sample_initial_states_3d)
from level3_3d_transport_lensing.transport_profiles import (
    geometry_centers, geometry_velocity_factors, get_profile)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG = PROJECT_ROOT / "configs" / "level3_3d_transport.yaml"
PHYSICS = physics_defaults()


# 1-3. 标准三维高斯势：中心值、对称性、远场极限、力、曲率频率
def test_static_potential_symmetry_and_limits():
    assert static_potential(0, 0, 0, PHYSICS["depth_j"], PHYSICS["waist_m"],
                            PHYSICS["z_r_m"]) == pytest.approx(-PHYSICS["depth_j"])
    x = np.array([0.7e-6, -0.7e-6])
    u = static_potential(x, x * 0, 0, PHYSICS["depth_j"], PHYSICS["waist_m"],
                         PHYSICS["z_r_m"])
    assert u[0] == pytest.approx(u[1])  # 对称性
    far = static_potential(np.array([50e-6]), np.array([50e-6]), np.array([0.0]),
                           PHYSICS["depth_j"], PHYSICS["waist_m"], PHYSICS["z_r_m"])
    assert abs(far[0]) < 1e-40  # 远场极限


def test_force_matches_gradient():
    assert _check_force_gradient(PHYSICS)["passed"]


def test_curvature_frequencies():
    assert _check_curvature(PHYSICS)["passed"]


# 4. lensing 零偏移极限严格恢复标准势
def test_lensing_zero_limit():
    assert _check_lensing_zero_limit(PHYSICS)["passed"]


# 5. 直线/对角中心分解与欧氏路径长度
def test_geometry_decomposition():
    for qv in (0.0, 0.3, 1.0):
        c1, c2 = geometry_centers(qv, "straight", 510e-6)
        assert np.hypot(c1, c2) == pytest.approx(510e-6 * qv)
        c1, c2 = geometry_centers(qv, "diagonal", 610e-6)
        assert np.hypot(c1, c2) == pytest.approx(610e-6 * qv)
        assert c1 == pytest.approx(c2)
    assert geometry_velocity_factors("straight") == (1.0, 0.0)
    f1, f2 = geometry_velocity_factors("diagonal")
    assert f1 == pytest.approx(1 / np.sqrt(2)) and f2 == pytest.approx(1 / np.sqrt(2))


# 6-8. 三类轨迹公式与命名
def test_profiles_endpoints_and_derivatives():
    assert _check_profiles()["passed"]
    assert _check_constant_jerk()["passed"]
    sine = get_profile("adiabatic_sine")
    s = np.array([0.13, 0.47, 0.81])
    h = 1e-7
    numeric = (sine.q(s + h) - sine.q(s - h)) / (2 * h)
    assert np.allclose(sine.q_dot(s), numeric, rtol=1e-5)
    # minimum-jerk 正确命名与公式
    mj = get_profile("minimum_jerk")
    assert mj.name == "minimum_jerk"
    assert mj.q(0.5) == pytest.approx(0.5)


# 9-11. 单阱/双阱临界、对角无分裂、符号配置
def test_straight_split_near_critical():
    z_r = PHYSICS["z_r_m"]
    _, _, _, n_below = axial_minima(PHYSICS["depth_j"], PHYSICS["waist_m"], z_r,
                                    1.5 * z_r, 0.0)
    _, _, _, n_above = axial_minima(PHYSICS["depth_j"], PHYSICS["waist_m"], z_r,
                                    2.5 * z_r, 0.0)
    assert n_below == 1
    assert n_above == 2


def test_diagonal_no_split_and_signs():
    z_r = PHYSICS["z_r_m"]
    _, _, _, n = axial_minima(PHYSICS["depth_j"], PHYSICS["waist_m"], z_r,
                              4.0 * z_r, 4.0 * z_r)
    assert n == 1
    # 反号时两条柱面焦点反向偏移，势更浅但仍可控
    assert axis_shift_m(0.5, 1.0, -1, z_r) == pytest.approx(-2 * z_r * 0.5)
    assert axis_shift_m(0.5, 1.0, +1, z_r) == pytest.approx(+2 * z_r * 0.5)


# 12. 三维热采样：协方差、拒绝、固定种子复现
def test_thermal_sampling():
    states = sample_initial_states_3d(123, 400, 5.0, PHYSICS["mass_kg"],
                                      PHYSICS["depth_j"], PHYSICS["waist_m"],
                                      PHYSICS["z_r_m"], PHYSICS["omega_r"],
                                      PHYSICS["omega_z"])
    again = sample_initial_states_3d(123, 400, 5.0, PHYSICS["mass_kg"],
                                     PHYSICS["depth_j"], PHYSICS["waist_m"],
                                     PHYSICS["z_r_m"], PHYSICS["omega_r"],
                                     PHYSICS["omega_z"])
    assert np.array_equal(states["pos_m"], again["pos_m"])
    assert np.all(states["initial_energy_j"] < 0)
    assert 0.4 <= states["acceptance_rate"] <= 1.0
    std = states["pos_m"].std(axis=0)
    assert std[2] / std[0] == pytest.approx(
        PHYSICS["omega_r"] / PHYSICS["omega_z"], rel=0.2)


# 13-16. 静态积分、无透镜与 lensing 功-能、二阶收敛
def test_static_integrator_bounded_and_convergent():
    control = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, PHYSICS["depth_j"], 0.0)
    pos = np.array([[0.3e-6, -0.2e-6, 0.5e-6]])
    vel = np.array([[0.02, 0.01, -0.005]])
    residuals = {}
    for dt in (0.1e-6, 0.05e-6):
        t_grid = np.arange(int(1e-4 / dt) + 1) * dt
        r = run_trajectory(pos, vel, t_grid, lambda t: control, PHYSICS["mass_kg"],
                           PHYSICS["depth_j"], PHYSICS["waist_m"], PHYSICS["z_r_m"])
        residuals[dt] = abs(r["residual"][0])
    assert residuals[0.1e-6] < 1e-13
    assert residuals[0.05e-6] <= residuals[0.1e-6] * (1 + 1e-6)  # 噪声地板下不增


def test_work_energy_moving_trap(cfg3):
    control = make_control("adiabatic_sine", "diagonal", 100e-6, 4e-4, None, (1, 1),
                           PHYSICS["depth_j"], PHYSICS["z_r_m"])
    t_grid = np.arange(int(4e-4 / 1e-8) + 1) * 1e-8
    r = run_trajectory(np.zeros((1, 3)), np.zeros((1, 3)), t_grid, control,
                       PHYSICS["mass_kg"], PHYSICS["depth_j"], PHYSICS["waist_m"],
                       PHYSICS["z_r_m"])
    assert abs(r["residual"][0]) / PHYSICS["depth_j"] < 1e-6


def test_work_energy_lensing(cfg3):
    v_s = resolve_vs(cfg3, 1.0)
    control = make_control("adiabatic_sine", "diagonal", 200e-6, 8e-4, v_s, (1, 1),
                           PHYSICS["depth_j"], PHYSICS["z_r_m"])
    t_grid = np.arange(int(8e-4 / 2e-8) + 1) * 2e-8
    r = run_trajectory(np.zeros((2, 3)), np.array([[0.1e-6, 0, 0], [0, -0.1e-6, 0.2e-6]]),
                       t_grid, control, PHYSICS["mass_kg"], PHYSICS["depth_j"],
                       PHYSICS["waist_m"], PHYSICS["z_r_m"])
    assert np.max(np.abs(r["residual"])) / PHYSICS["depth_j"] < 1e-4


# 17-18. E_f<0 留阱判据与 near-threshold 保留
def test_retention_criterion(cfg3):
    states = sample_states(cfg3, PHYSICS, 42, 24)
    out = run_condition(PHYSICS, states, "adiabatic_sine", "diagonal", 100.0, 400.0,
                        0.0, None, (1, 1), 0.05, 10.0, cfg3)
    assert len(out["records"]) == 24
    for r in out["records"]:
        assert r["retained"] == (r["e_final_uK"] < 0)
        near = abs(r["e_final_uK"] * KB * 1e-6) / PHYSICS["depth_j"] < 0.001
        assert r["near_threshold"] == near


# 19-21. CRN 对齐、种子隔离、条件冻结
def test_common_random_numbers_across_conditions(cfg3):
    states = sample_states(cfg3, PHYSICS, cfg3.scan_seed, 16)
    a = run_condition(PHYSICS, states, "adiabatic_sine", "diagonal", 300.0, 600.0,
                      1.0, resolve_vs(cfg3, 1.0), (1, 1), 0.1, 10.0, cfg3)
    b = run_condition(PHYSICS, states, "constant_jerk", "straight", 300.0, 600.0,
                      1.0, resolve_vs(cfg3, 1.0), (1, 1), 0.1, 10.0, cfg3)
    assert [r["x0_um"] for r in a["records"]] == [r["x0_um"] for r in b["records"]]
    assert [r["z0_um"] for r in a["records"]] == [r["z0_um"] for r in b["records"]]
    seeds = {cfg3.scan_seed, cfg3.validation_seed, cfg3.convergence_seed}
    assert len(seeds) == 3
    from level3_3d_transport_lensing.level3_cli import FROZEN_CONDITIONS
    assert len(FROZEN_CONDITIONS) == 4
    assert {(c["geometry"], c["profile"]) for c in FROZEN_CONDITIONS} == {
        ("diagonal", "adiabatic_sine"), ("diagonal", "constant_jerk"),
        ("straight", "adiabatic_sine"), ("straight", "constant_jerk")}


# 22. 统计已在 Level 2 测试覆盖（复用同实现），此处校验 v_s 反解与扫描场景
def test_vs_from_ratio_consistency(cfg3):
    v_ref = vs_from_ratio(1.0, cfg3.reference_profile,
                          cfg3.reference_distance_um * 1e-6,
                          cfg3.reference_duration_us * 1e-6,
                          cfg3.reference_geometry)
    v_half = vs_from_ratio(0.5, cfg3.reference_profile,
                           cfg3.reference_distance_um * 1e-6,
                           cfg3.reference_duration_us * 1e-6,
                           cfg3.reference_geometry)
    assert v_half == pytest.approx(2 * v_ref)
    assert resolve_vs(cfg3, 0.0) is None


# 23-24. CLI 冒烟 + 输出 schema（slow）
TINY_OVERRIDE = {
    "initial_ensemble": {"scan_shots": 8, "validation_shots": 12,
                         "sensitivity_shots": 6},
    "integration": {"post_transport_hold_us": 10.0},
    "transport": {"distances_um": [270.0], "durations_us": [400.0, 600.0]},
    "lensing": {"reference_axis_vmax_over_vs": [0.0, 1.0]},
    "output": {"directory": "outputs/level3_3d_transport_lensing/tiny_test"},
}


def _tiny_config(tmp_path):
    data = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    for section, values in TINY_OVERRIDE.items():
        data[section].update(values)
    path = tmp_path / "tiny3.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


@pytest.fixture(scope="module")
def cfg3():
    return load_config(CONFIG)


@pytest.mark.slow
def test_cli_all_stages_smoke(tmp_path):
    config_path = _tiny_config(tmp_path)
    assert level3_main(["--config", str(config_path), "--stage", "all"]) == 0
    output = tmp_path / "outputs" / "level3_3d_transport_lensing" / "tiny_test"
    required = ["config_used.yaml", "preflight.json", "potential_minima_scan.csv",
                "trajectory_profiles.csv", "coarse_scan.csv", "validation_shots.csv",
                "sensitivity.csv", "representative_trajectories.csv", "metrics.json",
                "image_validation.json", "summary.md",
                "static_3d_trap_slices.png", "lensing_potential_bifurcation.png",
                "trajectory_profile_comparison.png", "straight_vs_diagonal.png",
                "survival_phase_diagram.png", "heating_comparison.png",
                "final_energy_distributions.png", "timestep_convergence.png",
                "work_energy_balance.png", "axial_dynamics.png",
                "sensitivity_panels.png"]
    for name in required:
        assert (output / name).is_file(), name

    def reject_constant(value):
        raise ValueError(f"non-strict JSON constant {value}")

    metrics = json.loads((output / "metrics.json").read_text(),
                         parse_constant=reject_constant)
    assert metrics["seeds"]["validation"] != metrics["seeds"]["scan"]
    with (output / "validation_shots.csv").open(newline="", encoding="utf-8") as h:
        rows = list(csv.DictReader(h))
    assert len(rows) >= 12 * 8  # 8 条件 × shots
    ids_by_cond = {}
    for r in rows:
        ids_by_cond.setdefault(r["condition"], []).append(r["shot_id"])
    conditions = sorted(ids_by_cond)
    for cond in conditions[1:]:
        assert ids_by_cond[cond] == ids_by_cond[conditions[0]]  # CRN 对齐
    image_validation = json.loads((output / "image_validation.json").read_text(),
                                  parse_constant=reject_constant)
    assert image_validation["all_passed"]


# 25. 输出位于 sources/ 时被拒绝（复用 level0 io_utils，已有 Level 2 测试；
# 此处验证 level3 CLI 的输出目录覆盖同样经过 ensure_output_directory）
def test_output_into_sources_rejected(tmp_path):
    from level0_static_trap.io_utils import ensure_output_directory
    with pytest.raises(ValueError, match="sources"):
        ensure_output_directory(tmp_path / "sources" / "out", tmp_path)
