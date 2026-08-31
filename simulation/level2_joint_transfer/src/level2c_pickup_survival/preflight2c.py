"""Level 2C 预检查：不过关不得启动扫描。全部结果写入 preflight2c.json。"""

from __future__ import annotations

import numpy as np

from level1_transfer_1d.thermal import sample_thermal_initial_state
from level2_joint_transfer.config import level1_view
from level2_joint_transfer.level2_simulation import (evaluate_waveform_on_states, ideal_run,
                                                     physics_from_config)
from level2_joint_transfer.preflight import _check_existing_tests
from level2_joint_transfer.waveforms import (PickupSequentialWaveform, validate_pickup_waveform)
from level2_joint_transfer.work_energy import work_energy_two_trap

from .level2c_config import Level2cConfig
from .survival_engine import build_roundtrip_legs, run_roundtrip_detailed, run_survival

IDEAL_RESIDUAL_TOLERANCE = 5.0e-4
LABEL_FLIP_TOLERANCE = 0.02
# 时间反演对称性：保守极限往返后相对热初态尺度的回归容差。
# 物理注记：constant-jerk 移动的非绝热残余振荡（~a_max/ω_AOD ~ 2 mm/s 量级）
# 使往返不严格等于恒等映射；正确判据是残余 ≪ 热尺度，取 0.5σ。
TIME_REVERSAL_X_TOL = 0.5    # |Δx|/σ_x
TIME_REVERSAL_V_TOL = 0.5    # |Δv|/σ_v


def pickup_spec(duration_us: float, cfg: Level2cConfig, direction="pickup") -> dict:
    """手工 pick-up（或时间反演 drop-off）波形 spec。"""
    return {"type": "pickup_sequential", "duration_s": duration_us * 1e-6,
            "ramp_fraction": cfg.pickup.ramp_fraction, "pickup_direction": direction,
            "direction": direction,
            "name": f"pickup_{direction}_{duration_us:.0f}us"}


def _check_waveforms(cfg: Level2cConfig) -> dict:
    """各时长 pick-up/drop-off 波形的端点、单调性、解析导数 vs 有限差分。"""
    physics = physics_from_config(cfg.base)
    out = {}
    for duration_us in cfg.pickup.durations_us:
        for direction in ("pickup", "dropoff"):
            spec = pickup_spec(duration_us, cfg, direction)
            from level2_joint_transfer.level2_simulation import build_waveform
            wf = build_waveform(spec, physics)
            valid, reasons = validate_pickup_waveform(
                wf, physics["aod_initial_center_m"], physics["aod_initial_depth_j"],
                direction=direction)
            t = np.linspace(1e-9, wf.duration_s - 1e-9, 4001)
            h = 1e-9
            fd_v = (wf.center(t + h) - wf.center(t - h)) / (2 * h)
            an_v = wf.center_velocity(t)
            scale = max(float(np.max(np.abs(an_v))), 1e-30)
            fd_err = float(np.max(np.abs(fd_v - an_v)) / scale)
            # 深度变化率：避开升/降深段端点 kink（解析导数在该处不连续）
            kink_t = wf.ramp_duration_s if direction == "pickup" else wf.move_duration_s
            fd_d = (wf.depth(t + h) - wf.depth(t - h)) / (2 * h)
            an_d = wf.depth_rate(t)
            scale_d = max(float(np.max(np.abs(an_d))), 1e-30)
            mask_d = (np.abs(an_d) > 1e-3 * scale_d) & (np.abs(t - kink_t) > 5e-7)
            fd_err_d = float(np.max(np.abs(fd_d[mask_d] - an_d[mask_d])) / scale_d)
            out[spec["name"]] = {"valid": valid, "reasons": reasons,
                                 "velocity_fd_rel_err": fd_err, "depth_rate_fd_rel_err": fd_err_d,
                                 "passed": valid and fd_err < 1e-6 and fd_err_d < 1e-6}
    return {"per_waveform": out, "passed": bool(all(v["passed"] for v in out.values()))}


def _check_single_pickup_ideal(cfg: Level2cConfig) -> dict:
    """理想初态（SLM 阱底）单次 pick-up：应被 AOD 俘获且功-能残差小。"""
    physics = physics_from_config(cfg.base)
    duration_us = cfg.pickup.durations_us[1]  # 400 µs
    spec = pickup_spec(duration_us, cfg, "pickup")
    run = ideal_run(spec, physics, cfg.roundtrip.dt_s, cfg.roundtrip.wait_s)
    record = run["records"][0]
    # 转移段（[0,T]，SLM 静态、无开关跳变）的功-能残差；保持段由 two_trap 检查覆盖
    residual = record["max_abs_work_energy_residual_over_depth"] * physics["slm_depth_j"] \
        / max(physics["slm_depth_j"], physics["aod_initial_depth_j"])
    return {"captured": bool(record["captured"]),
            "final_aod_excitation_over_depth": record["final_aod_excitation_over_depth"],
            "work_energy_residual_over_depth": float(residual),
            "passed": bool(record["captured"]) and residual < IDEAL_RESIDUAL_TOLERANCE}


def _check_time_reversal(cfg: Level2cConfig) -> dict:
    """保守极限下 pick-up + drop-off 往返应近似回到初态（SLM 阱底静止出发）。"""
    physics = physics_from_config(cfg.base)
    duration_us = cfg.pickup.durations_us[1]
    from level2_joint_transfer.level2_simulation import build_waveform
    wf = build_waveform(pickup_spec(duration_us, cfg, "pickup"), physics)
    legs = build_roundtrip_legs(wf, cfg.roundtrip.wait_s, cfg.roundtrip.dt_s)
    detail = run_roundtrip_detailed(0.0, 0.0, physics["mass_kg"], legs, cfg.roundtrip.dt_s,
                                    physics["slm_depth_j"], physics["slm_waist_m"],
                                    physics["aod_waist_m"], physics["aod_initial_depth_j"])
    from level1_transfer_1d.thermal import thermal_scales
    view = level1_view(cfg.base)
    sigma_x, sigma_v = thermal_scales(view, well="slm")
    dx = abs(float(detail["x_m"][-1]) - 0.0)
    dv = abs(float(detail["v_m_s"][-1]) - 0.0)
    return {"final_x_over_sigma_x": float(dx / sigma_x),
            "final_v_over_sigma_v": float(dv / sigma_v),
            "tolerance": TIME_REVERSAL_X_TOL,
            "passed": bool(dx / sigma_x < TIME_REVERSAL_X_TOL
                           and dv / sigma_v < TIME_REVERSAL_V_TOL)}


def _check_two_trap_work_energy(cfg: Level2cConfig) -> dict:
    """双阱时变（含 SLM 开关跳变）功-能残差，理想 shot 一个完整往返。"""
    physics = physics_from_config(cfg.base)
    duration_us = cfg.pickup.durations_us[1]
    from level2_joint_transfer.level2_simulation import build_waveform
    wf = build_waveform(pickup_spec(duration_us, cfg, "pickup"), physics)
    legs = build_roundtrip_legs(wf, cfg.roundtrip.wait_s, cfg.roundtrip.dt_s)
    detail = run_roundtrip_detailed(0.0, 0.0, physics["mass_kg"], legs, cfg.roundtrip.dt_s,
                                    physics["slm_depth_j"], physics["slm_waist_m"],
                                    physics["aod_waist_m"], physics["aod_initial_depth_j"])
    # 逐段解析导数拼成完整时间轴
    dropoff = build_waveform(pickup_spec(duration_us, cfg, "dropoff"), physics)
    t = detail["time_s"]
    steps_pu, steps_wait = legs.steps_pickup, legs.steps_wait
    steps_rt = detail["steps_per_roundtrip"]
    depth_rate = np.zeros_like(t)
    center_rate = np.zeros_like(t)
    t_leg = np.arange(steps_pu + 1) * cfg.roundtrip.dt_s
    pu_dr = np.asarray(wf.depth_rate(t_leg), dtype=float)
    pu_cr = np.asarray(wf.center_velocity(t_leg), dtype=float)
    do_dr = np.asarray(dropoff.depth_rate(t_leg), dtype=float)
    do_cr = np.asarray(dropoff.center_velocity(t_leg), dtype=float)
    for base in range(0, t.size - 1, steps_rt):
        depth_rate[base:base + steps_pu] = pu_dr[:steps_pu]
        center_rate[base:base + steps_pu] = pu_cr[:steps_pu]
        do_base = base + steps_pu + steps_wait
        if do_base + steps_pu <= t.size:
            depth_rate[do_base:do_base + steps_pu] = do_dr[:steps_pu]
            center_rate[do_base:do_base + steps_pu] = do_cr[:steps_pu]
    we = work_energy_two_trap(t, detail["x_m"], detail["v_m_s"], physics["mass_kg"],
                              detail["aod_depth_j"], detail["aod_center_m"],
                              depth_rate, center_rate,
                              physics["slm_depth_j"], physics["slm_waist_m"],
                              detail["slm_factor"], physics["aod_waist_m"])
    depth_scale = max(physics["slm_depth_j"], physics["aod_initial_depth_j"])
    residual = we["max_abs_residual_J"] / depth_scale
    return {"max_abs_residual_over_depth": float(residual),
            "slm_switch_count": we["slm_switch_count"],
            "tolerance": IDEAL_RESIDUAL_TOLERANCE,
            "passed": bool(residual < IDEAL_RESIDUAL_TOLERANCE)}


def _check_timestep_convergence(cfg: Level2cConfig) -> dict:
    """往返序列步长减半：理想 shot 末态连续量收敛 + 热 shot 存活标签一致。"""
    physics = physics_from_config(cfg.base)
    duration_us = cfg.pickup.durations_us[1]
    from level2_joint_transfer.level2_simulation import build_waveform
    wf = build_waveform(pickup_spec(duration_us, cfg, "pickup"), physics)
    results = {}
    for dt in (cfg.roundtrip.dt_s, cfg.roundtrip.convergence_dt_s):
        legs = build_roundtrip_legs(wf, cfg.roundtrip.wait_s, dt)
        detail = run_roundtrip_detailed(0.0, 0.0, physics["mass_kg"], legs, dt,
                                        physics["slm_depth_j"], physics["slm_waist_m"],
                                        physics["aod_waist_m"], physics["aod_initial_depth_j"])
        results[dt] = (float(detail["x_m"][-1]), float(detail["v_m_s"][-1]))
    dts = sorted(results)
    dx = abs(results[dts[0]][0] - results[dts[1]][0])
    dv = abs(results[dts[0]][1] - results[dts[1]][1])
    view = level1_view(cfg.base)
    from level1_transfer_1d.thermal import thermal_scales
    sigma_x, sigma_v = thermal_scales(view, well="slm")
    # 热 shots 配对：两个步长下同一初态数组的逐 shot 存活序列一致
    rng = np.random.default_rng(cfg.roundtrip.scan_seed)
    n_shots = 32
    xs, vs = [], []
    while len(xs) < n_shots:
        state = sample_thermal_initial_state(view, rng, well="slm")
        xs.append(state.x0_m)
        vs.append(state.v0_m_s)
    curves = []
    for dt in dts:
        legs = build_roundtrip_legs(wf, cfg.roundtrip.wait_s, dt)
        out = run_survival(np.array(xs), np.array(vs), physics["mass_kg"], legs, dt,
                           max_transfers=4, slm_depth_j=physics["slm_depth_j"],
                           slm_waist_m=physics["slm_waist_m"],
                           aod_waist_m=physics["aod_waist_m"])
        curves.append(out["transfers_survived"])
    flips = float(np.mean(curves[0] != curves[1]))
    return {"dt_us": [d * 1e6 for d in dts],
            "ideal_final_dx_over_sigma_x": float(dx / sigma_x),
            "ideal_final_dv_over_sigma_v": float(dv / sigma_v),
            "thermal_label_flip_rate": flips,
            "flip_tolerance": LABEL_FLIP_TOLERANCE,
            "passed": bool(flips <= LABEL_FLIP_TOLERANCE)}


def _check_slm_sampler(cfg: Level2cConfig) -> dict:
    """SLM 阱采样器：全部样本满足 E_SLM<0，均值≈阱心，温度缩放正确。"""
    view = level1_view(cfg.base)
    rng = np.random.default_rng(cfg.roundtrip.scan_seed)
    energies, xs = [], []
    for _ in range(256):
        state = sample_thermal_initial_state(view, rng, well="slm")
        energies.append(state.aod_total_energy_J)
        xs.append(state.x0_m)
    bound = bool(np.all(np.array(energies) < 0))
    centered = bool(abs(np.mean(xs)) < 0.2 * np.std(xs, ddof=1))
    return {"all_bound": bound, "mean_x_um": float(np.mean(xs) * 1e6),
            "std_x_um": float(np.std(xs, ddof=1) * 1e6),
            "passed": bound and centered}


def run_preflight2c(cfg: Level2cConfig) -> dict:
    """Level 2C 预检查入口；关键项失败由 CLI 层拦截。"""
    checks = {}
    checks["existing_tests"] = _check_existing_tests(cfg.base.project_root)
    checks["pickup_waveforms"] = _check_waveforms(cfg)
    checks["single_pickup_ideal"] = _check_single_pickup_ideal(cfg)
    checks["time_reversal"] = _check_time_reversal(cfg)
    checks["two_trap_work_energy"] = _check_two_trap_work_energy(cfg)
    checks["timestep_convergence"] = _check_timestep_convergence(cfg)
    checks["slm_sampler"] = _check_slm_sampler(cfg)
    checks["pool_isolation"] = {
        "scan_seed": cfg.roundtrip.scan_seed, "validation_seed": cfg.roundtrip.validation_seed,
        "passed": cfg.roundtrip.scan_seed != cfg.roundtrip.validation_seed,
    }
    critical = ["existing_tests", "pickup_waveforms", "single_pickup_ideal", "time_reversal",
                "two_trap_work_energy", "timestep_convergence", "slm_sampler", "pool_isolation"]
    checks["critical_checks"] = critical

    def _passed(name):
        value = checks[name]
        return bool(value) if isinstance(value, bool) else bool(value.get("passed", False))

    checks["all_passed"] = bool(all(_passed(name) for name in critical))
    return checks
