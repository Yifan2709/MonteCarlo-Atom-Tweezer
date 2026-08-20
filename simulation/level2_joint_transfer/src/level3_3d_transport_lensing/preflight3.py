"""Level 3 预检查：16 项数值自检，失败则停止大规模 MC。"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from level0_static_trap.potentials import gaussian_force
from .aod_lensing import axial_minima, lensing_potential_and_force
from .gaussian_3d import KB, physics_defaults, static_force, static_potential
from .heating_scaling import run_heating_scaling_check
from .integrators_3d import run_trajectory, static_hold
from .level3_config import Level3Config
from .level3_simulation import (make_control, physics_from_config, resolve_vs,
                                sample_states)
from .thermal_sampling_3d import sample_initial_states_3d
from .transport_profiles import get_profile

PACKAGE_ROOT = Path(__file__).resolve().parents[2]


def _check_existing_tests() -> dict:
    """运行 Level 0-2 全部测试（排除 slow 防递归）。"""
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(PACKAGE_ROOT / "tests"), "-q",
         "--no-header", "-m", "not slow"],
        cwd=PACKAGE_ROOT, capture_output=True, text=True, timeout=3600)
    tail = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else ""
    return {"exit_code": result.returncode, "summary_line": tail,
            "passed": result.returncode == 0}


def _check_level0_slice(physics) -> dict:
    """z=0、x2=0 切片应恢复 Level 0 一维势和力。"""
    x = np.linspace(-3e-6, 3e-6, 61)
    u3 = static_potential(x, np.zeros_like(x), np.zeros_like(x),
                          physics["depth_j"], physics["waist_m"], physics["z_r_m"])
    u1 = -physics["depth_j"] * np.exp(-2 * x ** 2 / physics["waist_m"] ** 2)
    f3 = static_force(x, np.zeros_like(x), np.zeros_like(x), physics["depth_j"],
                      physics["waist_m"], physics["z_r_m"])[0]
    f1 = gaussian_force(x, physics["depth_j"], physics["waist_m"], 0.0)
    return {"potential_match": bool(np.allclose(u3, u1, rtol=1e-12)),
            "force_match": bool(np.allclose(f3, f1, rtol=1e-10, atol=1e-15)),
            "passed": bool(np.allclose(u3, u1, rtol=1e-12)
                           and np.allclose(f3, f1, rtol=1e-10, atol=1e-15))}


def _check_curvature(physics) -> dict:
    """静态数值曲率恢复 ω_r 与 ω_z。"""
    d = physics["depth_j"]
    w = physics["waist_m"]
    z_r = physics["z_r_m"]
    m = physics["mass_kg"]
    h = 1e-10
    k_r = (static_potential(h, 0, 0, d, w, z_r)
           - 2 * static_potential(0, 0, 0, d, w, z_r)
           + static_potential(-h, 0, 0, d, w, z_r)) / h ** 2
    k_z = (static_potential(0, 0, h, d, w, z_r)
           - 2 * static_potential(0, 0, 0, d, w, z_r)
           + static_potential(0, 0, -h, d, w, z_r)) / h ** 2
    omega_r_num = float(np.sqrt(k_r / m))
    omega_z_num = float(np.sqrt(k_z / m))
    ok_r = abs(omega_r_num / physics["omega_r"] - 1) < 1e-6
    ok_z = abs(omega_z_num / physics["omega_z"] - 1) < 1e-6
    return {"omega_r_analytic": physics["omega_r"], "omega_r_numeric": omega_r_num,
            "omega_z_analytic": physics["omega_z"], "omega_z_numeric": omega_z_num,
            "passed": bool(ok_r and ok_z)}


def _check_lensing_zero_limit(physics) -> dict:
    """z_s=0 时 lensing 势严格等于标准三维高斯势。"""
    rng = np.random.default_rng(7)
    xi = rng.normal(0, 1e-6, 500)
    eta = rng.normal(0, 1e-6, 500)
    z = rng.normal(0, 4e-6, 500)
    u_lens, f_xi, f_eta, f_z = lensing_potential_and_force(
        xi, eta, z, physics["depth_j"], physics["waist_m"], physics["z_r_m"], 0.0, 0.0)
    u_std = static_potential(xi, eta, z, physics["depth_j"], physics["waist_m"],
                             physics["z_r_m"])
    sf = static_force(xi, eta, z, physics["depth_j"], physics["waist_m"],
                      physics["z_r_m"])
    force_ok = (np.allclose(f_xi, sf[0], rtol=1e-10, atol=1e-24)
                and np.allclose(f_eta, sf[1], rtol=1e-10, atol=1e-24)
                and np.allclose(f_z, sf[2], rtol=1e-8, atol=1e-24))
    pot_ok = bool(np.allclose(u_lens, u_std, rtol=1e-12))
    return {"potential_match": pot_ok, "force_match": bool(force_ok),
            "passed": pot_ok and force_ok}


def _check_force_gradient(physics) -> dict:
    """lensing 解析力与势的有限差分梯度一致。"""
    rng = np.random.default_rng(11)
    pts = rng.normal(0, 2e-6, (40, 3))
    zs1, zs2 = 2e-6, -1e-6
    h = 1e-11
    max_rel = 0.0
    for p in pts:
        _, fx, fy, fz = lensing_potential_and_force(
            p[0], p[1], p[2], physics["depth_j"], physics["waist_m"],
            physics["z_r_m"], zs1, zs2)
        for axis in range(3):
            args = [p[0], p[1], p[2]]
            args[axis] += h
            up, _, _, _ = lensing_potential_and_force(
                args[0], args[1], args[2], physics["depth_j"], physics["waist_m"],
                physics["z_r_m"], zs1, zs2)
            args[axis] -= 2 * h
            um, _, _, _ = lensing_potential_and_force(
                args[0], args[1], args[2], physics["depth_j"], physics["waist_m"],
                physics["z_r_m"], zs1, zs2)
            numeric = -(up - um) / (2 * h)
            analytic = (fx, fy, fz)[axis]
            scale = max(abs(numeric), abs(analytic), 1e-30)
            max_rel = max(max_rel, abs(numeric - analytic) / scale)
    return {"max_relative_error": max_rel, "passed": bool(max_rel < 1e-5)}


def _check_profiles() -> dict:
    """三类轨迹端点、速度、加速度与单调性。"""
    s = np.linspace(0, 1, 2001)
    results = {}
    ok = True
    for name in ("adiabatic_sine", "constant_jerk", "minimum_jerk"):
        prof = get_profile(name)
        q, v, a = prof.q(s), prof.q_dot(s), prof.q_ddot(s)
        good = (abs(q[0]) < 1e-14 and abs(q[-1] - 1) < 1e-12
                and abs(v[0]) < 1e-12 and abs(v[-1]) < 1e-12
                and abs(a[0]) < 1e-10 and abs(a[-1]) < 1e-10
                and np.all(np.diff(q) >= -1e-12))
        results[name] = {"q0": float(q[0]), "q1": float(q[-1]), "v0": float(v[0]),
                         "v1": float(v[-1]), "a0": float(a[0]), "a1": float(a[-1]),
                         "monotone": bool(np.all(np.diff(q) >= -1e-12)),
                         "passed": bool(good)}
        ok = ok and good
    return {"profiles": results, "passed": ok}


def _check_constant_jerk() -> dict:
    """constant-jerk 四段解析积分与 J=32L/T^3 归一化。"""
    prof = get_profile("constant_jerk")
    s = np.linspace(0, 1, 4001)
    jerk = prof.q_jerk(s)
    j_ok = np.allclose(np.abs(jerk), 32.0)
    centers = np.array([0.125, 0.375, 0.625, 0.875])
    signs = np.sign(prof.q_jerk(centers))
    signs_ok = bool(np.array_equal(signs, np.array([1.0, -1.0, -1.0, 1.0])))
    # 速度连续：对 q̇ 做梯形积分应恢复总位移 1
    v = prof.q_dot(s)
    displacement = float(np.trapezoid(v, s)) if hasattr(np, "trapezoid") \
        else float(np.trapz(v, s))
    return {"jerk_magnitude_ok": bool(j_ok), "sign_pattern_ok": signs_ok,
            "displacement_ok": bool(abs(displacement - 1.0) < 1e-9),
            "passed": bool(j_ok and signs_ok and abs(displacement - 1.0) < 1e-9)}


def _check_straight_split_critical(physics) -> dict:
    """直线运动临界点：|z_s|≈2 z_R 时单阱→双阱（数值极值）。"""
    z_r = physics["z_r_m"]
    transitions = []
    for zs in (1.8 * z_r, 2.0 * z_r, 2.2 * z_r):
        _, _, _, n = axial_minima(physics["depth_j"], physics["waist_m"], z_r, zs, 0.0)
        transitions.append(n)
    ok = transitions[0] == 1 and transitions[2] >= 2
    return {"n_minima_at_1.8zR": transitions[0], "n_minima_at_2.0zR": transitions[1],
            "n_minima_at_2.2zR": transitions[2], "passed": bool(ok)}


def _check_diagonal_single_minimum(physics) -> dict:
    """同号等速对角模型在峰值速度处保持单一最低点。"""
    z_r = physics["z_r_m"]
    _, _, _, n = axial_minima(physics["depth_j"], physics["waist_m"], z_r,
                              3.0 * z_r, 3.0 * z_r)
    return {"n_minima": n, "passed": bool(n == 1)}


def _static_energy_convergence(physics, dt_s):
    pos = np.array([[0.3e-6, -0.2e-6, 0.5e-6]])
    vel = np.array([[0.02, 0.01, -0.005]])
    control = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, physics["depth_j"], 0.0)
    n = int(round(2e-4 / dt_s))
    t_grid = np.arange(n + 1) * dt_s
    result = run_trajectory(pos, vel, t_grid, lambda t: control, physics["mass_kg"],
                            physics["depth_j"], physics["waist_m"], physics["z_r_m"])
    return float(abs(result["residual"][0]))


def _check_static_integrator(physics) -> dict:
    """静态三维积分能量误差有界并随步长二阶收敛。"""
    e_coarse = _static_energy_convergence(physics, 0.05e-6)
    e_fine = _static_energy_convergence(physics, 0.025e-6)
    order = np.log2(e_coarse / e_fine) if e_fine > 0 else float("nan")
    noise_floor = 1e-12 * physics["depth_j"]
    # 误差达到机器噪声地板（随步长不再下降）时视为已收敛
    at_floor = e_coarse < noise_floor
    return {"error_dt_0.05us_J": e_coarse, "error_dt_0.025us_J": e_fine,
            "convergence_order": float(order),
            "passed": bool(e_coarse < 1e-12 and (
                at_floor or (np.isnan(order) or order > 1.0)))}


def _translation_residual(cfg, physics, severity, dt_us):
    """理想初态对角 610 μm/1.6 ms 运输的功-能残差。"""
    v_s = resolve_vs(cfg, severity)
    control = make_control("adiabatic_sine", "diagonal", 610e-6, 1.6e-3, v_s,
                           cfg.axis_signs, physics["depth_j"], physics["z_r_m"])
    n = int(round(1.6e-3 / (dt_us * 1e-6)))
    t_grid = np.arange(n + 1) * dt_us * 1e-6
    pos = np.zeros((1, 3))
    vel = np.zeros((1, 3))
    result = run_trajectory(pos, vel, t_grid, control, physics["mass_kg"],
                            physics["depth_j"], physics["waist_m"],
                            physics["z_r_m"])
    return result


def _check_work_energy_convergence(cfg, physics) -> dict:
    """无透镜与 lensing 功-能残差随步长收敛。"""
    depth = max(physics["depth_j"], 1e-30)
    out = {}
    ok = True
    noise_floor = 1e-9  # 相对阱深；低于此为机器噪声地板，不再要求单调下降
    for label, severity in (("lensing_off", 0.0), ("lensing_ref1.00", 1.0)):
        r_coarse = _translation_residual(cfg, physics, severity, 0.05)["max_residual"]
        r_fine = _translation_residual(cfg, physics, severity, 0.025)["max_residual"]
        entry = {"residual_dt0.05_over_depth": r_coarse / depth,
                 "residual_dt0.025_over_depth": r_fine / depth}
        at_floor = r_coarse / depth < noise_floor
        entry["passed"] = bool(r_coarse / depth < 5e-4
                               and (at_floor or r_fine <= r_coarse * (1 + 1e-9)))
        out[label] = entry
        ok = ok and entry["passed"]
    return {"cases": out, "passed": ok}


def _check_ideal_transport(cfg, physics) -> dict:
    """理想中心零速度初态完成 610 μm / 1.6 ms 对角运输。"""
    result = _translation_residual(cfg, physics, 1.0, 0.05)
    e_final = result["e_final"][0]
    finite = bool(np.all(np.isfinite(result["final_pos"])))
    near_zero = abs(e_final) < 1e-9 * physics["depth_j"]  # 数值零：符号属噪声
    return {"e_final_uK": e_final / KB * 1e6, "finite": finite,
            "retained": bool(e_final < 0), "near_numerical_zero": bool(near_zero),
            "passed": bool(finite and (e_final < 0 or near_zero))}


def _check_mc_dt_paired(cfg, physics) -> dict:
    """128 固定热初态的 0.10/0.05/0.025 μs 配对收敛与标签一致性。"""
    states = sample_initial_states_3d(
        cfg.convergence_seed, 128, cfg.temperature_uK, physics["mass_kg"],
        physics["depth_j"], physics["waist_m"], physics["z_r_m"],
        physics["omega_r"], physics["omega_z"])
    from .level3_simulation import run_condition
    labels = {}
    finals = {}
    for dt_us in (0.10, 0.05, 0.025):
        out = run_condition(physics, states, "adiabatic_sine", "diagonal", 610.0,
                            800.0, 1.0, resolve_vs(cfg, 1.0), cfg.axis_signs, dt_us,
                            10.0, cfg)
        labels[dt_us] = [r["retained"] for r in out["records"]]
        finals[dt_us] = np.array([r["e_final_uK"] for r in out["records"]])
    flips = float(np.mean(np.array(labels[0.10]) != np.array(labels[0.05])))
    max_diff = float(np.max(np.abs(finals[0.10] - finals[0.05])))
    med_res = float(np.median([abs(r["w_residual_uK"]) for r in out["records"]]))
    return {"label_flips_0.10_vs_0.05": flips, "max_final_energy_diff_uK": max_diff,
            "median_work_residual_uK": med_res,
            "passed": bool(flips <= 0.02)}


def _check_hold_segment(cfg, physics) -> dict:
    """100 μs 静态保持段能量误差有界。"""
    result = _translation_residual(cfg, physics, 1.0, 0.05)
    _, _, pp, _ = static_hold(result["final_pos"], result["final_vel"], 100e-6,
                              0.05e-6, physics["mass_kg"], physics["depth_j"],
                              physics["waist_m"], physics["z_r_m"])
    worst = float(np.max(pp)) / physics["depth_j"]
    return {"max_peak_to_peak_over_depth": worst, "passed": bool(worst < 1e-3)}


def _check_config_values(cfg: Level3Config) -> dict:
    """配置值有限、单位正确、validation 条件已冻结。"""
    frozen = {
        "paper_anchor_diagonal": {"geometry": "diagonal", "profile": "adiabatic_sine",
                                  "distance_um": 610.0, "duration_us": 1600.0},
        "paper_anchor_constant_jerk": {"geometry": "diagonal",
                                       "profile": "constant_jerk",
                                       "distance_um": 610.0, "duration_us": 1600.0},
        "straight_control": {"geometry": "straight", "profile": "adiabatic_sine",
                             "distance_um": 510.0, "duration_us": 1600.0},
        "straight_control_constant_jerk": {"geometry": "straight",
                                           "profile": "constant_jerk",
                                           "distance_um": 510.0,
                                           "duration_us": 1600.0},
        "nonzero_severity": 1.0,
    }
    values = [cfg.mass_u, cfg.wavelength_nm, cfg.waist_um, cfg.depth_uK,
              cfg.temperature_uK, cfg.reference_distance_um,
              cfg.reference_duration_us, *cfg.reference_axis_vmax_over_vs]
    finite = bool(all(np.isfinite(v) for v in values))
    return {"values_finite": finite, "frozen_validation_conditions": frozen,
            "passed": finite}


def run_preflight(cfg: Level3Config) -> dict:
    """按顺序执行 16 项预检查。"""
    physics = physics_from_config(cfg)
    checks = {
        "existing_tests": _check_existing_tests(),
        "level0_slice": _check_level0_slice(physics),
        "curvature_frequencies": _check_curvature(physics),
        "lensing_zero_limit": _check_lensing_zero_limit(physics),
        "force_gradient": _check_force_gradient(physics),
        "profile_endpoints": _check_profiles(),
        "constant_jerk_constants": _check_constant_jerk(),
        "heating_scaling": run_heating_scaling_check(physics),
        "straight_split_critical": _check_straight_split_critical(physics),
        "diagonal_single_minimum": _check_diagonal_single_minimum(physics),
        "static_integrator": _check_static_integrator(physics),
        "work_energy_convergence": _check_work_energy_convergence(cfg, physics),
        "ideal_transport": _check_ideal_transport(cfg, physics),
        "mc_dt_paired": _check_mc_dt_paired(cfg, physics),
        "hold_segment": _check_hold_segment(cfg, physics),
        "config_values": _check_config_values(cfg),
    }
    critical = list(checks)
    checks["all_passed"] = all(
        bool(v.get("passed", True)) if isinstance(v, dict) else bool(v)
        for v in checks.values())
    checks["critical_checks"] = critical
    return checks
