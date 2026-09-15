"""Level 5 命令行入口：preflight / channel / reference_rb / interleaved_rb /
array_extension / sensitivity / all。

执行顺序：加载上游 → 生成/复用冻结轨迹池 → preflight（失败停止正式 RB）
→ 单通道表征 → RB development → 冻结设计 → 47-site 正式 RB/IRB →
195-site 扩展 → 敏感性 → 报告与图片。
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import yaml

from level0_static_trap.io_utils import (ensure_output_directory, save_json,
                                         save_markdown, save_yaml,
                                         validate_pngs)
from level3_3d_transport_lensing.gaussian_3d import KB
from level4_roundtrip_coherence import protocol_segments as ps
from level4_roundtrip_coherence import upstream_artifacts as l4_upstream
from level4_roundtrip_coherence.level4_cli import _write_rows_csv
from level4_roundtrip_coherence.level4_config import load_config \
    as load_level4_config

from . import level5_visualization as viz
from .array_ensemble import (build_site_ensemble, common_local_decomposition,
                             site_resolved_summary)
from .channel_tomography import (CARDINAL_INPUTS, _reference_inputs,
                                 affine_fit,
                                 average_state_fidelity_from_affine,
                                 characterize_from_states, ptm_bootstrap)
from .clifford_group import CliffordGroup
from .dd_finite_pulses import validate_transformed_frame
from .interleaved_rb import generate_irb_sequences, run_interleaved_rb
from .level5_config import Level5Config, config_as_dict, load_config
from .level5_statistics import hierarchical_bootstrap, seed_for, wilson_interval
from .noise_psd import (DetuningNoiseModel, validate_psd_synthesis,
                        welch_psd)
from .preflight5 import Preflight5
from .pulse_library import make_bare_pulse, make_scrofulous_pulses, \
    scrofulous_angles
from .qubit_hamiltonian import (INPUT_STATES, equatorial_rotation_unitary,
                                return_probability_z, statevector_to_bloch)
from .reference_rb import generate_reference_sequences, run_reference_rb
from .rb_fitting import (early_window_slope, fit_exponential_rb,
                         fit_irb_decay, instantaneous_ratios, irb_estimate,
                         model_selection, safe_fit_summary)
from .spin_channels import MoveChannelEngine, StaticGateEngine
from .trajectory_replay import (TrajectoryPool,
                                generate_trajectory_pool,
                                replay_phase_for_pool)

TWO_PI = 2.0 * np.pi
HBAR = 1.054571817e-34
US = 1e-6


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Level 5 有限脉冲单比特动力学与多站点 RB/IRB")
    parser.add_argument("--config", required=True)
    parser.add_argument("--stage", choices=("preflight", "channel",
                                            "reference_rb",
                                            "interleaved_rb",
                                            "array_extension",
                                            "sensitivity", "all"),
                        default="all")
    parser.add_argument("--output-dir", help="覆盖输出目录")
    parser.add_argument("--fast-preflight", action="store_true",
                        help="开发模式：preflight 跳过 Level 0-4 测试重跑"
                             "（正式运行禁止使用）")
    return parser


# --------------------------------------------------------------- 上游/池
def _discover_level4_output(cfg: Level5Config) -> Path:
    """定位 Level 4 demo 输出目录。"""
    if cfg.upstream_level4_output != "auto_discover":
        p = Path(cfg.upstream_level4_output)
        if not p.is_absolute():
            p = cfg.project_root / p
        return p.resolve()
    hits = sorted((cfg.project_root / "outputs").glob(
        "level4_end_to_end_roundtrip/*"))
    hits = [h for h in hits if h.is_dir() and (h / "metrics.json").is_file()]
    if not hits:
        raise FileNotFoundError("未找到 Level 4 输出目录（需要 metrics.json）")
    return hits[0].resolve()


def _build_upstream_manifest(cfg: Level5Config, l4_dir: Path,
                             pools: dict) -> dict:
    """构造 Level 5 的上游 manifest（绝对路径 + SHA-256）。"""
    import hashlib
    from datetime import datetime, timezone

    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()

    files = {}
    for key, rel in (("level4_metrics", "metrics.json"),
                     ("level4_frozen_protocols",
                      "frozen_validation_protocols.yaml"),
                     ("level4_config_used", "config_used.yaml"),
                     ("level4_phase_records", "phase_records.csv"),
                     ("level4_repeated_round_metrics",
                      "repeated_round_metrics.csv"),
                     ("level4_representative_trajectories",
                      "representative_trajectories.csv")):
        path = l4_dir / rel
        if path.is_file():
            files[key] = {"absolute_path": str(path),
                          "sha256": sha256(path),
                          "role": f"Level 4 冻结产物 {rel}"}
    l2 = cfg.project_root / "outputs/level2_joint_transfer/demo"
    l3 = cfg.project_root / "outputs/level3_3d_transport_lensing/demo"
    for base, tag in ((l2, "level2"), (l3, "level3")):
        for rel in ("best_waveform.yaml", "metrics.json"):
            path = base / rel
            if path.is_file():
                files[f"{tag}_{rel.split('.')[0]}"] = {
                    "absolute_path": str(path), "sha256": sha256(path),
                    "role": f"{tag} 冻结产物 {rel}"}
    for name, pool in pools.items():
        files[f"pool_{name}"] = {
            "absolute_path": str(pool.get("_path", "")),
            "sha256": pool.get("_sha256", ""),
            "role": f"Level 5 冻结轨迹池 {name}"}
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "require_sha256_manifest": cfg.require_sha256_manifest,
        "files": files,
        "trajectory_mode": cfg.trajectory_mode,
        "key_parameters": {
            "drive_rabi_frequency_kHz": cfg.drive_rabi_frequency_khz,
            "eta_slm": cfg.eta_slm, "eta_aod": cfg.eta_aod,
            "dd_transport_mode": cfg.transport_mode,
            "pulse_family": cfg.default_family,
            "nominal_noise_enabled": cfg.nominal_noise_enabled,
        },
        "provenance": {
            "drive_rabi_frequency_kHz": "paper_anchor (Fig. 4a)",
            "eta": "paper_theory_scale + level4 assumption (1.3e-4)",
            "durations": "level4_frozen (Extended Data Fig. 10e 锚点)",
            "noise": "assumed（nominal 关闭）",
        },
    }


def _ensure_pools(cfg: Level5Config, cfg4, protocols, output_dir: Path,
                  log=print) -> dict:
    """生成或复用冻结轨迹池（NPZ 缓存）。"""
    pool_dir = output_dir / "frozen_pools"
    pool_dir.mkdir(parents=True, exist_ok=True)
    specs = [
        ("protocol_a", protocols["protocol_a"], cfg.pool_full_trajectories),
        ("static_idle", protocols["static_idle"], cfg.pool_full_trajectories),
        ("protocol_b", protocols["protocol_b"],
         cfg.pool_secondary_trajectories),
        ("transfer_only", ps.build_transfer_only(
            cfg4, l4_upstream.build_manifest(cfg4)[1]
            ["nominal_vs_m_per_s"]), cfg.pool_secondary_trajectories),
        ("transport_only", ps.build_transport_only(
            cfg4, l4_upstream.build_manifest(cfg4)[1]
            ["nominal_vs_m_per_s"]), cfg.pool_secondary_trajectories),
        ("protocol_a_no_lensing", protocols["protocol_a_no_lensing"],
         cfg.pool_secondary_trajectories),
    ]
    pools = {}
    for name, timeline, n_traj in specs:
        path = pool_dir / f"{name}.npz"
        if path.is_file():
            pool = TrajectoryPool.load(path)
            if pool.n_traj >= n_traj:
                pool = _shrink_pool(pool, n_traj)
                pools[name] = pool
                pools[name]._meta_extra = {"path": str(path)}
                log(f"[pools] 复用 {name}（{pool.n_traj} 条）")
                continue
        pool = generate_trajectory_pool(
            name, timeline, cfg4, cfg.pool_seed, n_traj,
            dt_classical_us=cfg.pool_dt_classical_us,
            record_stride_us=cfg.pool_record_stride_us,
            chunk_size=cfg.pool_chunk_shots, log=log)
        pool.save(path)
        pools[name] = pool
        log(f"[pools] 生成 {name}（{pool.n_traj} 条，"
            f"sha256={pool.sha256()[:12]}）")
    return pools


def _shrink_pool(pool: TrajectoryPool, n: int) -> TrajectoryPool:
    """把池截断到前 n 条（复用缓存时保持规模一致）。"""
    if pool.n_traj <= n:
        return pool
    import dataclasses
    return dataclasses.replace(
        pool, u_slm=pool.u_slm[:n], u_aod=pool.u_aod[:n],
        roundtrip_success=pool.roundtrip_success[:n],
        final_retained=pool.final_retained[:n],
        first_failure_stage=pool.first_failure_stage[:n],
        u_slm_final=pool.u_slm_final[:n], u_aod_final=pool.u_aod_final[:n],
        level4_a_slm=pool.level4_a_slm[:n],
        level4_a_aod=pool.level4_a_aod[:n],
        meta={**pool.meta, "n_traj": n, "shrunk_from": pool.n_traj})


# --------------------------------------------------------- channel 表征
def _propagate_six_states(pool: TrajectoryPool, omega: float, dd_mode: str,
                          cfg: Level5Config, n_traj: int, family="bare",
                          eta_slm=None, eta_aod=None, rabi_scale=None,
                          delta_site=None, amp_error=0.0,
                          free_noise_traces=None, envelope=None,
                          edge_us=None):
    """对 6 个输入态批量传播通道并返回 [n_traj, 6, 3] Bloch 输出。"""
    eta_slm = cfg.eta_slm if eta_slm is None else eta_slm
    eta_aod = cfg.eta_aod if eta_aod is None else eta_aod
    envelope = cfg.envelope if envelope is None else envelope
    edge_us = cfg.cosine_edge_us if edge_us is None else edge_us
    n = min(n_traj, pool.n_traj)
    engine = MoveChannelEngine(pool, omega, dd_mode=dd_mode,
                               pulse_family=family, envelope=envelope,
                               edge_us=edge_us)
    bloch = np.empty((n, 6, 3))
    for s, name in enumerate(CARDINAL_INPUTS):
        psi = np.tile(INPUT_STATES[name], (n, 1)).astype(complex)
        engine.propagate(psi, np.arange(n), eta_slm, eta_aod,
                         dt_free_us=cfg.dt_free_us,
                         dt_pulse_us=cfg.dt_pulse_us,
                         rabi_scale=rabi_scale, delta_site=delta_site,
                         amp_error=amp_error,
                         free_noise_traces=free_noise_traces)
        bloch[:, s, :] = statevector_to_bloch(psi)
    surv = pool.roundtrip_success[:n].astype(bool)
    return bloch, surv, engine


def stage_channel(cfg: Level5Config, pools: dict, output_dir: Path,
                  plot_records: list, log=print) -> dict:
    """Stage B：多条件单通道 PTM 表征。"""
    omega = TWO_PI * cfg.drive_rabi_frequency_khz * 1e3
    group = CliffordGroup()
    conditions = []
    for dd in ("none", "spin_echo", "xy4", "xy16"):
        conditions.append(("static_idle", dd, "nominal"))
    conditions += [("transfer_only", "xy4", "nominal"),
                   ("transport_only", "none", "nominal"),
                   ("transport_only", "xy4_per_long_move", "nominal"),
                   ("transport_only", "xy16", "nominal"),
                   ("protocol_a", "none", "nominal"),
                   ("protocol_a", "spin_echo", "nominal"),
                   ("protocol_a", "xy4_per_long_move", "nominal"),
                   ("protocol_a", "xy16", "nominal"),
                   ("protocol_b", "none", "nominal"),
                   ("protocol_b", "xy4_per_long_move", "nominal"),
                   ("protocol_a_no_lensing", "xy4_per_long_move",
                    "lensing_off"),
                   ("protocol_a", "xy4_per_long_move", "noise_quasistatic"),
                   ("protocol_a", "xy4_per_long_move", "noise_ou")]
    rows = []
    metrics = {}
    bloch_axes = {}
    for pool_name, dd, variant in conditions:
        pool = pools[pool_name]
        n = cfg.channel_trajectories if pool_name in ("protocol_a",
                                                      "static_idle") \
            else cfg.channel_dev_trajectories * 2
        n = min(n, pool.n_traj)
        noise_traces = None
        surv_override = None
        if variant == "noise_quasistatic":
            rng = np.random.default_rng(seed_for("noise_calibration_pool"))
            qs = rng.normal(0, TWO_PI * 20.0, size=n)
            # quasistatic：每轨迹一个常数失谐 → 传给 engine 的 noise_qs
            # （此处复用 free_noise_traces 接口：常数列）
            noise_traces = np.repeat(qs[:, None], pool.times_s.size, axis=1)
        elif variant == "noise_ou":
            rng = np.random.default_rng(seed_for(
                "noise_calibration_pool", 1))
            noise_traces = _ou_traces(rng, n, pool.times_s.size,
                                      pool.sample_interval_s,
                                      TWO_PI * cfg.ou_rms_hz
                                      if cfg.ou_rms_hz > 0 else TWO_PI * 10.0,
                                      cfg.ou_tau_us * US)
        started = time.time()
        eta_s, eta_a = cfg.eta_slm, cfg.eta_aod
        bloch, surv, engine = _propagate_six_states(
            pool, omega, dd, cfg, n, family="bare",
            eta_slm=eta_s, eta_aod=eta_a,
            free_noise_traces=noise_traces)
        result = characterize_from_states(bloch, surv)
        boot = ptm_bootstrap(bloch, surv, n_boot=cfg.bootstrap_samples,
                             seed=seed_for("channel_validation_pool"))
        result["bootstrap"] = boot
        result["condition"] = f"{pool_name}|{dd}|{variant}"
        result["n_trajectories"] = int(n)
        result["dd_pulse_count"] = engine.n_pulses
        result["free_vs_pulse_time"] = {
            "free_s": engine.phase_diagnostics(
                np.arange(min(8, n)), eta_s, eta_a)[
                "free_evolution_time_s"],
            "pulse_s": engine.phase_diagnostics(
                np.arange(min(8, n)), eta_s, eta_a)["pulse_total_time_s"]}
        metrics[result["condition"]] = result
        m = np.asarray(result["ptm_m"])
        rows.append({"condition": result["condition"], "pool": pool_name,
                     "dd": dd, "variant": variant,
                     "n_traj": n, "survival": result.get("survival"),
                     "f_avg_conditional": result["f_avg"],
                     "f_avg_survival_weighted":
                         result.get("f_avg_survival_weighted"),
                     "unitarity": result["unitarity"]["unitarity"],
                     "coherent_angle_rad":
                         result["coherent"]["coherent_angle_rad"],
                     "f_avg_ci_low": boot["f_avg_ci_low"],
                     "f_avg_ci_high": boot["f_avg_ci_high"],
                     "m_xx": m[0, 0], "m_yy": m[1, 1], "m_zz": m[2, 2],
                     "ptm_m": result["ptm_m"], "affine_c": result["affine_c"],
                     "n_dd_pulses": engine.n_pulses})
        bloch_axes[result["condition"]] = {
            axis: (m @ np.eye(3)[{"x": 0, "y": 1, "z": 2}[axis]]).tolist()
            for axis in ("x", "y", "z")}
        log(f"[channel] {result['condition']}: F_avg="
            f"{result['f_avg']:.6f} survival={result.get('survival'):.4f}"
            f"（{time.time() - started:.1f}s）")
    _write_rows_csv(output_dir / "channel_ptm.csv",
                    [{**r, "ptm_m": json.dumps(r["ptm_m"]),
                      "affine_c": json.dumps(r["affine_c"])}
                     for r in rows])
    save_json(output_dir / "channel_metrics.json",
              {"conditions": metrics,
               "conventions_note": "conditional PTM 未把损失归一化；"
                                   "survival 与 survival-weighted 并列报告",
               "assumed_noise_variants": ["noise_quasistatic (20 Hz, assumed)",
                                          "noise_ou (10 Hz, τ=200 μs, "
                                          "assumed)"]})
    plot_records.append(viz.plot_channel_ptm(
        [{k: r[k] for k in ("condition", "ptm_m")} for r in rows],
        output_dir / "channel_ptm.png"))
    plot_records.append(viz.plot_channel_bloch_sphere(
        dict(list(bloch_axes.items())[:6]),
        output_dir / "channel_bloch_sphere.png"))
    return {"rows": rows, "metrics": metrics}


def _ou_traces(rng, n, n_samples, dt, sigma_rad_s, tau_s):
    """向量化 OU 轨迹 [n, n_samples]（rad/s 失谐）。"""
    decay = np.exp(-dt / tau_s)
    scale = sigma_rad_s * np.sqrt(max(0.0, 1.0 - decay ** 2))
    out = np.empty((n, n_samples), dtype=np.float32)
    out[:, 0] = sigma_rad_s * rng.standard_normal(n)
    noise = rng.standard_normal((n, n_samples - 1)).astype(np.float32)
    for i in range(1, n_samples):
        out[:, i] = out[:, i - 1] * decay + scale * noise[:, i - 1]
    return out


# ----------------------------------------------------------- reference RB
def _run_reference_stage(cfg, pools, group, n_sites, n_seq_per_length,
                         pool_tag, seed, log):
    """执行一组 reference RB（返回逐序列行与汇总）。"""
    static_pool = pools["static_idle"]
    omega = TWO_PI * cfg.drive_rabi_frequency_khz * 1e3
    gate_engine = StaticGateEngine(static_pool, omega,
                                   family=cfg.default_family,
                                   envelope=cfg.envelope, group=group)
    seqs = generate_reference_sequences(group, cfg.rb_lengths,
                                        n_seq_per_length, seed)
    result = run_reference_rb(seqs, gate_engine, group, static_pool,
                              cfg.eta_slm, cfg.eta_aod,
                              seed + 1, n_sites=n_sites)
    rows = []
    probs = result["return_probabilities"]
    durations = result["durations_s"]
    for i, (length, seq_idx) in enumerate(result["sequence_ids"]):
        for site in range(n_sites):
            rows.append({
                "pool_tag": pool_tag, "length": int(length),
                "sequence_index": int(seq_idx), "site": int(site),
                "n_gates": int(result["n_gates_per_sequence"][i]),
                "return_probability": float(probs[i, site]),
                "duration_s": float(durations[i]),
                "k_static_trace": int(result["k_trace"][i]),
                "pulse_family": cfg.default_family,
                "dd": "none",
                "survival": 1.0,
                "joint_return": float(probs[i, site]),
                "fit_group": f"ref|{pool_tag}",
                "seed": int(seed),
            })
    return result, rows


def _fit_reference_rb(cfg, rows, log):
    """reference RB 拟合（多模型 + 诊断）。"""
    lengths = np.array(sorted({r["length"] for r in rows}))
    by_len = {n: [r["return_probability"] for r in rows
                  if r["length"] == n] for n in lengths}
    means = np.array([np.mean(by_len[n]) for n in lengths])
    stds = np.array([np.std(by_len[n]) for n in lengths])
    sigma = np.array([max(np.std(by_len[n]), 1e-4) / np.sqrt(len(by_len[n]))
                      for n in lengths])
    selection = model_selection(lengths, means, sigma)
    exp_fit = selection["exponential"]
    p_hat = exp_fit.get("params", {}).get("p") if exp_fit.get("ok") else None
    fit_curves = {}
    if exp_fit.get("ok"):
        pars = exp_fit["params"]
        fit_curves["exponential"] = ([pars["A"], pars["p"], pars["B"]],
                                     viz.C_A)
    if selection.get("stretched", {}).get("ok"):
        pars = selection["stretched"]["params"]
        fit_curves["stretched"] = ([pars["A"], pars["tau"], pars["beta"],
                                    pars["B"]], viz.C_B)
    residuals = (means - (exp_fit["params"]["A"]
                          * exp_fit["params"]["p"] ** lengths
                          + exp_fit["params"]["B"])) if exp_fit.get("ok") \
        else np.zeros_like(means)
    # hierarchical bootstrap：length × (sequence×site) 两层（无独立噪声层）
    matrix = np.array([[r["return_probability"] for r in rows
                        if r["length"] == n] for n in lengths])
    n_rep = matrix.shape[1]
    hier = hierarchical_bootstrap(
        np.repeat(matrix[:, :, None], 1, axis=2), len(lengths),
        max(1, n_rep), 1,
        n_boot=min(cfg.statistics_bootstrap_samples, 500), seed=4242)
    return {"lengths": lengths, "means": means, "stds": stds,
            "selection": selection, "p": p_hat,
            "fit_curves": {k: v[0] for k, v in fit_curves.items()},
            "residuals": residuals,
            "sequence_std": stds,
            "early_slope": early_window_slope(lengths, means),
            "hierarchical": hier,
            "assumptions": {
                "gate_independence": "有限脉冲门均为同一 Clifford 分布"
                                     "（twirling 假设按构造满足到门无关"
                                     "误差阶）",
                "markovianity": "quasistatic 光移失谐 shot-to-shot 独立"
                                "时近似成立；长序列内相关噪声会偏离单指数",
                "single_exponential_preferred": selection.get("preferred"),
            }}


def stage_reference_rb(cfg, pools, output_dir, plot_records, log=print):
    """Stage C/D：reference RB（dev 8 序列 → formal 72 序列 × 47 sites）。"""
    group = CliffordGroup()
    all_rows = []
    summaries = {}
    for tag, n_seq, n_sites, seed in (
            ("dev", cfg.rb_dev_sequences_per_length, cfg.development_sites,
             seed_for("reference_rb_pool", 0)),
            ("formal_47", cfg.rb_sequences_per_length, cfg.validation_sites,
             cfg.rb_sequence_seed)):
        result, rows = _run_reference_stage(cfg, pools, group, n_sites,
                                            n_seq, tag, seed, log)
        all_rows.extend(rows)
        summaries[tag] = _fit_reference_rb(cfg, rows, log)
        p = summaries[tag]["p"]
        log(f"[reference_rb:{tag}] p="
            f"{p if p is None else round(p, 6)}")
    _write_rows_csv(output_dir / "reference_rb_sequence_results.csv",
                    all_rows)
    # 图
    formal = summaries["formal_47"]
    scatter = [[r["return_probability"] for r in all_rows
                if r["pool_tag"] == "formal_47" and r["site"] == 0
                and r["length"] == n] for n in formal["lengths"]]
    plot_records.append(viz.plot_reference_rb_decay(
        {"lengths": formal["lengths"], "means": formal["means"],
         "stds": formal["stds"], "scatter_by_length": scatter,
         "fits": {k: (v, color) for k, (v, color) in {
             "exponential": (formal["fit_curves"].get(
                 "exponential"), viz.C_A),
             "stretched": (formal["fit_curves"].get("stretched"),
                           viz.C_B)}.items()}},
        output_dir / "reference_rb_decay.png"))
    plot_records.append(viz.plot_rb_residuals(
        {"lengths": formal["lengths"], "residuals": formal["residuals"],
         "sequence_std": formal["sequence_std"], "irb_residuals": None},
        output_dir / "rb_residuals.png"))
    return {"rows": all_rows, "summaries": summaries}


# ------------------------------------------------------- interleaved RB
def _run_irb_stage(cfg, pools, group, n_sites, n_seq, pool_tag, seed, log,
                   site_params=None, move_pool_name="protocol_a"):
    """执行一组 IRB（返回逐序列行与原始结果）。"""
    move_pool = pools[move_pool_name]
    static_pool = pools["static_idle"]
    omega = TWO_PI * cfg.drive_rabi_frequency_khz * 1e3
    move_engine = MoveChannelEngine(move_pool, omega,
                                    dd_mode=cfg.transport_mode,
                                    pulse_family="bare",
                                    envelope=cfg.envelope)
    gate_engine = StaticGateEngine(static_pool, omega,
                                   family=cfg.default_family,
                                   envelope=cfg.envelope, group=group)
    seqs = generate_irb_sequences(group, cfg.irb_total_cliffords,
                                  cfg.irb_interleaved_counts, n_seq, seed)
    rabi_scale = None
    delta_site = None
    eta_scale = None
    if site_params is not None:
        rabi_scale = site_params["params"]["rabi_scale"]
        delta_site = TWO_PI * site_params["params"]["static_detuning_hz"]
        eta_scale = site_params["params"]["eta_scale"]
    started = time.time()
    result = run_interleaved_rb(
        seqs, move_engine, gate_engine, move_pool, static_pool,
        cfg.eta_slm, cfg.eta_aod, seed + 2, n_sites=n_sites,
        rabi_scale=rabi_scale, delta_site=delta_site,
        dt_free_us=cfg.dt_free_us, dt_pulse_us=cfg.dt_pulse_us,
        site_eta_scale=eta_scale)
    rows = []
    probs = result["return_probabilities"]
    surv = result["survival"]
    for i, (m_val, seq_idx) in enumerate(result["sequence_ids"]):
        for site in range(n_sites):
            rows.append({
                "pool_tag": pool_tag, "M": int(m_val),
                "sequence_index": int(seq_idx), "site": int(site),
                "string_id": int(result["string_ids"][i]),
                "n_gates": int(cfg.irb_total_cliffords + 1),
                "survival": float(surv[i, site]),
                "conditional_return": float(probs[i, site]),
                "joint_return": float(probs[i, site] * surv[i, site]),
                "loss_detected_return": float(probs[i, site]
                                              * surv[i, site]),
                "k_static_trace": int(result["k_static"][
                    result["string_ids"][i]]),
                "pulse_family": cfg.default_family,
                "dd": cfg.transport_mode,
                "operation": cfg.irb_operation,
                "n_moves": int(m_val),
                "fit_group": f"irb|{pool_tag}",
                "seed": int(seed),
                "clifford_count": int(cfg.irb_total_cliffords),
            })
    log(f"[irb:{pool_tag}] 完成（{time.time() - started:.0f}s）")
    return result, rows


def _summarize_irb(rows):
    """IRB 汇总：survival / conditional / joint 三口径曲线与拟合。"""
    m_vals = np.array(sorted({r["M"] for r in rows}))
    out = {}
    for key in ("survival", "conditional_return", "joint_return"):
        per_m = [np.array([r[key] for r in rows if r["M"] == m])
                 for m in m_vals]
        out[key] = {"means": [float(x.mean()) for x in per_m],
                    "stds": [float(x.std()) for x in per_m],
                    "n": [int(x.size) for x in per_m]}
    cond_fit = fit_exponential_rb(m_vals, np.array(
        out["conditional_return"]["means"]))
    joint_fit = fit_exponential_rb(m_vals, np.array(
        out["joint_return"]["means"]))
    surv_fit = fit_irb_decay(m_vals, np.array(out["survival"]["means"]),
                             "clipped_boltzmann_x_depol")
    out["fits"] = {"conditional_exponential": safe_fit_summary(cond_fit),
                   "joint_exponential": safe_fit_summary(joint_fit),
                   "survival_clipped_boltzmann_x_depol":
                       safe_fit_summary(surv_fit)}
    out["instantaneous_ratios"] = {
        "conditional": instantaneous_ratios(
            out["conditional_return"]["means"],
            out["survival"]["means"]),
        "survival": instantaneous_ratios(out["survival"]["means"])}
    out["early_slope"] = {
        key: early_window_slope(m_vals, out[key]["means"])
        for key in ("survival", "conditional_return", "joint_return")}
    return out, m_vals


def stage_interleaved_rb(cfg, pools, output_dir, plot_records, log=print):
    """Stage C/D：IRB development 与 47-site formal（含 site ensemble）。"""
    group = CliffordGroup()
    all_rows = []
    raw = {}
    dev_result, dev_rows = _run_irb_stage(
        cfg, pools, group, cfg.development_sites,
        cfg.irb_dev_sequences_per_count, "dev",
        seed_for("interleaved_rb_pool", 0), log)
    all_rows.extend(dev_rows)
    ens_47 = build_site_ensemble(
        {"trap_depth_relative_std": cfg.trap_depth_relative_std,
         "rabi_relative_std": cfg.rabi_relative_std,
         "static_detuning_std_Hz": cfg.static_detuning_std_Hz,
         "eta_relative_std": cfg.eta_relative_std,
         "spatial_noise_correlation": cfg.spatial_noise_correlation},
        cfg.validation_sites, cfg.site_seed)
    formal_result, formal_rows = _run_irb_stage(
        cfg, pools, group, cfg.validation_sites,
        cfg.irb_sequences_per_count, "formal_47", cfg.irb_sequence_seed,
        log, site_params=ens_47)
    all_rows.extend(formal_rows)
    _write_rows_csv(output_dir / "interleaved_rb_sequence_results.csv",
                    all_rows)
    summary, m_vals = _summarize_irb(formal_rows)
    raw["formal_47"] = {"summary": summary, "m_values": m_vals}
    # IRB estimate（适用性检查；reference p 从 reference CSV 读取）
    cond_p = (summary["fits"]["conditional_exponential"].get("params")
              or {}).get("p")
    ref_fit = _ref_fit_for_metrics(output_dir, cfg)
    p_ref = ref_fit.get("p") if isinstance(ref_fit, dict) else None
    if cond_p is not None and p_ref is not None:
        irb = irb_estimate(float(p_ref), float(cond_p),
                           design="fixed_cliffords_variable_moves")
        irb["p_ref"] = p_ref
        irb["p_per_move_fit"] = cond_p
    else:
        irb = {"interpretable": False,
               "reason": "reference/interleaved 拟合不可用"}
    exp_ok = summary["fits"]["conditional_exponential"].get(
        "status") == "ok"
    ratios_ok = summary["instantaneous_ratios"]["conditional"].get(
        "stable")
    irb["assumption_checks"] = {
        "conditional_exponential_fit_ok": bool(exp_ok),
        "instantaneous_ratios_stable": ratios_ok,
        "survival_tail_risk_set":
            summary["instantaneous_ratios"]["survival"].get("min_risk_set"),
    }
    if not (exp_ok and ratios_ok):
        irb["interpretable"] = False
        irb["reason"] = (irb.get("reason", "") + "；"
                         "conditional 衰减非单指数或相邻比值不稳定 → "
                         "IRB estimate not interpretable under standard "
                         "assumptions")
    summary["irb_estimate"] = irb
    summary["assumptions"] = {
        "gate_independence": "frozen 轨迹 draw i.i.d.、门与通道交替 → "
                             "Markov/twirl 假设到一阶",
        "markovianity": "move 通道重复独立（无跨 move 长相关噪声）",
        "loss_convention": "joint = survival × conditional；"
                           "loss-detected 口径数值相同",
        "verdict": ("IRB estimate interpretable under standard assumptions"
                    if irb.get("interpretable") else
                    "IRB estimate not interpretable under standard "
                    "assumptions"),
    }
    save_json(output_dir / "rb_fit_metrics.json", {
        "reference_rb": ref_fit,
        "interleaved_rb": summary,
        "loss_conventions": ["conditional spin return", "atomic survival",
                             "survival-weighted (joint)",
                             "loss-detected"],
        "early_window": summary["early_slope"],
        "instantaneous_ratios": summary["instantaneous_ratios"]})
    # 图
    plot_payload = {**summary, "m_values": m_vals}
    plot_records.append(viz.plot_interleaved_rb_decay(
        plot_payload, output_dir / "interleaved_rb_decay.png"))
    plot_records.append(viz.plot_survival_vs_spin_return(
        plot_payload, output_dir / "survival_vs_spin_return.png"))
    means = np.array(summary["conditional_return"]["means"])
    resid = np.zeros_like(means)
    if cond_p:
        pred = 0.5 + 0.5 * float(cond_p) ** m_vals
        resid = means - pred
    plot_records.append(viz.plot_rb_residuals(
        {"lengths": m_vals, "residuals": resid,
         "sequence_std": np.array(summary["conditional_return"]["stds"]),
         "irb_residuals": resid, "irb_m": m_vals},
        output_dir / "rb_residuals.png"))
    return {"rows": all_rows, "summary": summary, "m_values": m_vals,
            "raw": raw, "formal_rows": formal_rows,
            "site_ensemble": ens_47["summary"],
            "site_ensemble_full": ens_47}


def _ref_fit_for_metrics(output_dir: Path, cfg):
    """读取 reference RB 拟合（若无则从 CSV 重算）。"""
    csv_path = output_dir / "reference_rb_sequence_results.csv"
    if not csv_path.is_file():
        return {"status": "not_run"}
    import csv as _csv
    rows = list(_csv.DictReader(csv_path.open(encoding="utf-8")))
    rows = [r for r in rows if r.get("pool_tag") == "formal_47"]
    if not rows:
        return {"status": "no_formal_rows"}
    lengths = np.array(sorted({int(r["length"]) for r in rows}))
    means = np.array([np.mean([float(r["return_probability"])
                               for r in rows
                               if int(r["length"]) == n])
                      for n in lengths])
    sigma = np.array([max(np.std([float(r["return_probability"])
                                  for r in rows
                                  if int(r["length"]) == n]), 1e-4)
                      / np.sqrt(sum(1 for r in rows
                                    if int(r["length"]) == n))
                      for n in lengths])
    selection = model_selection(lengths, means, sigma)
    p = (selection["exponential"].get("params") or {}).get("p")
    return {"status": "ok", "p": p, "lengths": lengths.tolist(),
            "means": means.tolist(), "selection": {
                k: safe_fit_summary(v) for k, v in selection.items()
                if isinstance(v, dict)}}


# --------------------------------------------------------- array 扩展
def _site_rows_from(rows, n_sites, ens, tag):
    """IRB 行 → site-resolved 指标行。"""
    probs = np.array([r["conditional_return"] for r in rows]).reshape(
        -1, n_sites)
    surv = np.array([r["survival"] for r in rows]).reshape(-1, n_sites)
    site_rows = []
    for s in range(n_sites):
        site_rows.append({
            "pool_tag": tag, "site": s,
            "survival": float(surv[:, s].mean()),
            "conditional_return": float(probs[:, s].mean()),
            "joint_return": float((probs[:, s] * surv[:, s]).mean()),
            "rabi_scale": float(ens["params"]["rabi_scale"][s]),
            "depth_scale": float(ens["params"]["depth_scale"][s]),
            "eta_scale": float(ens["params"]["eta_scale"][s]),
            "static_detuning_hz":
                float(ens["params"]["static_detuning_hz"][s]),
        })
    stats_cond = site_resolved_summary(probs, ens["params"])
    stats_surv = site_resolved_summary(surv, ens["params"])
    site_data = {"survival": surv.mean(axis=0).tolist(),
                 "conditional_return": probs.mean(axis=0).tolist(),
                 "joint_return": (probs * surv).mean(axis=0).tolist(),
                 "params": {k: np.asarray(v).tolist()
                            for k, v in ens["params"].items()
                            if isinstance(v, np.ndarray)},
                 "n_sites": n_sites}
    return site_rows, stats_cond, stats_surv, site_data


def stage_array_extension(cfg, pools, output_dir, plot_records, log=print,
                          formal_rows_47=None, ens_47=None):
    """Stage E：195-site 扩展（冻结模型，仅改 site 数）+ site 统计。"""
    group = CliffordGroup()
    summaries = {}
    site_data = {}
    arr_rows = []
    if formal_rows_47 and ens_47 is not None:
        rows_47, stats_cond, stats_surv, sd47 = _site_rows_from(
            formal_rows_47, cfg.validation_sites, ens_47, "sites_47")
        _write_rows_csv(output_dir / "site_resolved_metrics_sites_47.csv",
                        rows_47)
        summaries["sites_47"] = {"site_stats": {
            "conditional": stats_cond, "survival": stats_surv}}
        site_data["sites_47"] = sd47
    ens_195 = build_site_ensemble(
        {"trap_depth_relative_std": cfg.trap_depth_relative_std,
         "rabi_relative_std": cfg.rabi_relative_std,
         "static_detuning_std_Hz": cfg.static_detuning_std_Hz,
         "eta_relative_std": cfg.eta_relative_std,
         "spatial_noise_correlation": cfg.spatial_noise_correlation},
        cfg.extension_sites, cfg.site_seed + 1)
    result, rows_195 = _run_irb_stage(
        cfg, pools, group, cfg.extension_sites,
        cfg.irb_sequences_per_count, "sites_195", cfg.irb_sequence_seed,
        log, site_params=ens_195)
    summary_195, m_vals = _summarize_irb(rows_195)
    rows_s195, stats_cond, stats_surv, sd195 = _site_rows_from(
        rows_195, cfg.extension_sites, ens_195, "sites_195")
    _write_rows_csv(output_dir / "site_resolved_metrics_sites_195.csv",
                    rows_s195)
    _write_rows_csv(output_dir / "site_resolved_metrics.csv",
                    rows_s195)  # 兼容 prompt 命名（195-site 主表）
    summaries["sites_195"] = {"summary": summary_195,
                              "site_ensemble": ens_195["summary"],
                              "site_stats": {"conditional": stats_cond,
                                             "survival": stats_surv},
                              "m_values": m_vals.tolist()}
    site_data["sites_195"] = sd195
    common_local = common_local_decomposition(
        np.asarray(sd195["conditional_return"])[None, :])
    summaries["sites_195"]["common_local_note"] = common_local
    log(f"[array:sites_195] conditional mean="
        f"{np.mean(sd195['conditional_return']):.5f}, between-site std="
        f"{stats_cond['between_site_std']:.5f}")
    for tag in summaries:
        st = summaries[tag]["site_stats"]["conditional"]
        surv_m = summaries[tag]["site_stats"]["survival"]["array_average"]
        arr_rows.append({
            "pool_tag": tag, "n_sites": len(
                site_data[tag]["conditional_return"]),
            "array_conditional_return": st["array_average"],
            "between_site_std": st["between_site_std"],
            "within_site_std": st["within_site_std"],
            "worst_decile": st["worst_decile"],
            "best_decile": st["best_decile"],
            "array_survival": surv_m,
        })
    _write_rows_csv(output_dir / "array_summary.csv", arr_rows)
    plot_records.append(viz.plot_site_resolved_array(
        site_data.get("sites_47", sd195), output_dir
        / "site_resolved_array.png"))
    return {"summaries": summaries, "rows_195": rows_195,
            "site_data": site_data}


# ------------------------------------------------------------ 敏感性
def stage_sensitivity(cfg, pools, output_dir, plot_records, log=print):
    """敏感性/消融（冻结主结果后运行）。"""
    omega = TWO_PI * cfg.drive_rabi_frequency_khz * 1e3
    pool = pools["protocol_a"]
    n = min(cfg.sensitivity_trajectories, pool.n_traj)
    rows = []

    def run_factor(factor, value, **kwargs):
        variant = dict(dd_mode=cfg.transport_mode, family="bare",
                       eta_slm=cfg.eta_slm, eta_aod=cfg.eta_aod,
                       rabi_scale=None, delta_site=None, amp_error=0.0,
                       free_noise_traces=None, envelope=cfg.envelope,
                       edge_us=None)
        variant.update(kwargs)
        bloch, surv, engine = _propagate_six_states(
            pool, omega, variant["dd_mode"], cfg, n,
            family=variant["family"], eta_slm=variant["eta_slm"],
            eta_aod=variant["eta_aod"], rabi_scale=variant["rabi_scale"],
            delta_site=variant["delta_site"],
            amp_error=variant["amp_error"],
            free_noise_traces=variant["free_noise_traces"],
            envelope=variant["envelope"], edge_us=variant["edge_us"])
        result = characterize_from_states(bloch, surv)
        rows.append({
            "factor": factor, "value": str(value), "n_traj": n,
            "f_avg_conditional": result["f_avg"],
            "survival": result.get("survival"),
            "f_avg_survival_weighted":
                result.get("f_avg_survival_weighted"),
            "unitarity": result["unitarity"]["unitarity"],
            "coherent_angle_rad":
                result["coherent"]["coherent_angle_rad"],
            "dd": variant["dd_mode"],
            "note": "quantum-only 因子不改变经典 survival（frozen replay）",
        })
        log(f"[sensitivity] {factor}={value}: "
            f"F={result['f_avg']:.6f}")

    run_factor("baseline", "nominal")
    for scale in (0.98, 1.02):
        run_factor("rabi_scale", scale, rabi_scale=np.full(n, scale))
    for det in (-50.0, 50.0):
        run_factor("delta_drive_hz", det,
                   delta_site=np.full(n, TWO_PI * det))
    for eta in (cfg.eta_sensitivity if cfg.eta_sensitivity
                else (1.5e-4,)):
        run_factor("eta", eta, eta_slm=eta, eta_aod=eta)
    rng = np.random.default_rng(seed_for("sensitivity_pool"))
    for sigma in (10.0, 50.0):
        qs = rng.normal(0, TWO_PI * sigma, size=(n, pool.times_s.size))
        run_factor("quasistatic_detuning_rms_hz", sigma,
                   free_noise_traces=qs.astype(np.float32))
    for rms in (5.0, 20.0):
        ou = _ou_traces(rng, n, pool.times_s.size,
                        pool.sample_interval_s, TWO_PI * rms,
                        cfg.ou_tau_us * US)
        run_factor("ou_rms_hz_tau200us", rms, free_noise_traces=ou)
    for amp in (0.01, 0.03):
        run_factor("pulse_amp_error", amp, amp_error=amp)
    run_factor("envelope", "cosine_edge_1us", envelope="cosine_edge",
               edge_us=1.0)
    for dd in ("none", "spin_echo", "xy4", "xy16",
               "xy4_per_long_move"):
        run_factor("dd_mode", dd, dd_mode=dd)
    # paper 场景：11.4% 阱深不均匀（量子侧线性缩放 η_SLM）
    rng2 = np.random.default_rng(seed_for("sensitivity_pool", 1))
    depth = 1.0 + rng2.normal(0, cfg.paper_scenario_trap_depth_std, n)
    run_factor("paper_trap_depth_std_114pct", "site_depth_spread",
               eta_slm=cfg.eta_slm * depth.mean(),
               eta_aod=cfg.eta_aod)
    _write_rows_csv(output_dir / "sensitivity.csv", rows)
    plot_records.append(viz.plot_sensitivity_panels(
        rows, output_dir / "sensitivity_panels.png"))
    return {"rows": rows}


# ------------------------------------------------------ 静态产物与诊断图
def _write_pulse_library_csv(cfg, output_dir):
    """pulse_library.csv：bare/SCROFULOUS/DD 脉冲规格。"""
    from .pulse_library import PulseSpec, VirtualZ
    omega = TWO_PI * cfg.drive_rabi_frequency_khz * 1e3
    rows = []
    for angle, label in ((0.5 * np.pi, "pi/2"), (np.pi, "pi")):
        for family in ("bare", "scrofulous"):
            pulses = ([make_bare_pulse(0.0, angle, omega, cfg.envelope,
                                       cfg.cosine_edge_us)]
                      if family == "bare" else
                      make_scrofulous_pulses(0.0, angle, omega,
                                             cfg.envelope,
                                             cfg.cosine_edge_us))
            for i, p in enumerate(pulses):
                row = p.as_row()
                row["role"] = f"{family}_{label}_comp{i}"
                rows.append(row)
    from .dd_finite_pulses import dd_pulse_events
    for mode in ("spin_echo", "xy4", "xy8", "xy16"):
        for p in dd_pulse_events(mode, 0.0, 100e-6, omega,
                                 envelope=cfg.envelope):
            row = p.as_row()
            row["role"] = f"dd_{mode}"
            rows.append(row)
    rows.append(VirtualZ(np.pi).as_row())
    _write_rows_csv(output_dir / "pulse_library.csv", rows)
    return rows


def _write_clifford_table_csv(cfg, group, output_dir):
    """clifford_table.csv（含面积统计与论文 2.02π 比较）。"""
    rows = group.as_rows(TWO_PI * cfg.drive_rabi_frequency_khz * 1e3,
                         family=cfg.default_family)
    _write_rows_csv(output_dir / "clifford_table.csv", rows)
    from collections import Counter
    counts = Counter(r["n_physical_pulses"] for r in rows)
    # SCROFULOUS 每旋转面积 = 2θc+π；分解里每个物理 π/2 旋转一个 composite
    areas = []
    for i in range(24):
        n_rot = sum(1 for k, v in group.decompose(i) if k == "P")
        ang = scrofulous_angles(0.5 * np.pi)
        areas.append(n_rot * ang.total_area / np.pi)
    return {"rows": rows, "pulse_count_histogram": dict(counts),
            "avg_scrofulous_area_pi_per_clifford": float(np.mean(areas)),
            "paper_anchor_avg_area_pi": 2.02}


def _write_noise_artifacts(cfg, output_dir, plot_records):
    """noise_psd.csv / noise_validation.json + 验证图。"""
    spec = [{"model": "white", "params": {"level": 1e-4}},
            {"model": "one_over_f", "params": {"level": 2e-4, "alpha": 1.0}},
            {"model": "lorentzian", "params": {"sigma": 10.0,
                                               "tau_s": 2e-4}},
            {"model": "line", "params": {"f_hz": 60.0, "width_hz": 2.0,
                                         "amp": 1e-3, "harmonics": 3}}]
    validation = validate_psd_synthesis(spec, n_samples=20000, fs=1e6,
                                        seed=11)
    rows = []
    for f, tgt, m, lo, hi in zip(validation["freqs"], validation["target"],
                                 validation["welch_mean"],
                                 validation["welch_lo"],
                                 validation["welch_hi"]):
        rows.append({"freq_hz": f, "target_psd": tgt, "welch_mean": m,
                     "welch_lo": lo, "welch_hi": hi})
    _write_rows_csv(output_dir / "noise_psd.csv", rows)
    noise_model = DetuningNoiseModel(
        cfg.detuning_model if cfg.nominal_noise_enabled else "none",
        rms=TWO_PI * cfg.quasistatic_detuning_rms_hz,
        tau_s=None if cfg.ou_tau_us <= 0 else cfg.ou_tau_us * US,
        psd_spec=cfg.psd_spec, seed=cfg.noise_seed)
    save_json(output_dir / "noise_validation.json", {
        "psd_validation": validation,
        "nominal_model": noise_model.describe(),
        "nominal_enabled": cfg.nominal_noise_enabled,
        "assumed_parameters": ["quasistatic_rms", "ou_rms", "ou_tau",
                               "psd levels", "60 Hz line strengths"],
        "note": "nominal 噪声关闭；所有非零噪声参数均为 assumed 敏感性"})
    plot_records.append(viz.plot_noise_psd_validation(
        validation, output_dir / "noise_psd_validation.png"))
    return validation


def _diagnostic_plots(cfg, pools, group, output_dir, plot_records,
                      l4_metrics_path):
    """脉冲响应、SCROFULOUS 鲁棒性、Ramsey/echo/DD、Level4 对照、收敛图。"""
    from .quantum_propagators import apply_su2, su2_propagator
    omega = TWO_PI * cfg.drive_rabi_frequency_khz * 1e3
    # 图1：脉冲响应
    t_us = np.linspace(0, 120, 1201)
    dt = (t_us[1] - t_us[0]) * 1e-6
    pops = [0.0]
    psi = np.array([[1.0, 0.0]], dtype=complex)
    for _ in range(t_us.size - 1):
        apply_su2(psi, *su2_propagator(np.zeros(1), np.full(1, omega),
                                       0.0, dt))
        pops.append(1.0 - abs(psi[0, 0]) ** 2)
    pops = np.array(pops)
    det_axis = np.linspace(-5000, 5000, 201)
    final_pop = []
    for det in det_axis:
        psi = np.array([[1.0, 0.0]], dtype=complex)
        apply_su2(psi, *su2_propagator(np.zeros(1), np.full(1, omega), 0.0,
                                       0.5 * np.pi / omega))
        apply_su2(psi, *su2_propagator(np.full(1, TWO_PI * det),
                                       np.zeros(1), 0.0, np.pi / omega))
        final_pop.append(1.0 - abs(psi[0, 0]) ** 2)
    envelopes = {}
    from .pulse_library import envelope_value
    for env, edge in (("square", 0.0), ("cosine_edge", 2.0)):
        t_env = np.linspace(0, 40, 401)
        envelopes[f"{env}({edge}μs)"] = (
            t_env, [envelope_value(env, t / 40, edge, t, 40.0)
                    for t in t_env])
    plot_records.append(viz.plot_pulse_response(
        {"rabi_t_us": t_us, "rabi_population": pops,
         "detuning_scan_hz": det_axis,
         "detuning_final_population": np.array(final_pop),
         "envelopes": envelopes},
        output_dir / "pulse_response.png"))
    # 图2：SCROFULOUS 鲁棒性
    scan = {"targets": ["pi/2", "pi"], "eps_axis": np.linspace(-0.1, 0.1, 41),
            "infidelity": {}, "det_axis": np.linspace(-2000, 2000, 81),
            "det_infidelity": {}}
    for label, theta in (("pi/2", 0.5 * np.pi), ("pi", np.pi)):
        scan["infidelity"][label] = {}
        for family in ("bare", "scrofulous"):
            infs = []
            for eps in scan["eps_axis"]:
                if family == "bare":
                    u = equatorial_rotation_unitary(0.0, theta * (1 + eps))
                else:
                    ang = scrofulous_angles(theta, 0.0)
                    p1, p2, _ = ang.pulse_phases
                    u = (equatorial_rotation_unitary(p1, ang.theta_c
                                                     * (1 + eps))
                         @ equatorial_rotation_unitary(p2, np.pi * (1 + eps))
                         @ equatorial_rotation_unitary(p1, ang.theta_c
                                                       * (1 + eps)))
                tgt = equatorial_rotation_unitary(0.0, theta)
                prod = u @ tgt.conj().T
                ph = prod.flat[int(np.argmax(np.abs(prod)))]
                ph = ph / abs(ph)
                infs.append(float(np.abs(prod / ph - np.eye(2)).max() ** 2
                                  / 2))
            scan["infidelity"][label][family] = np.array(infs)
    for family in ("bare", "scrofulous"):
        infs = []
        theta = 0.5 * np.pi
        for det in scan["det_axis"]:
            ang = scrofulous_angles(theta, 0.0)
            p1, p2, _ = ang.pulse_phases
            if family == "scrofulous":
                u = (equatorial_rotation_unitary(p1, ang.theta_c)
                     @ _detuning_sandwich(p2, np.pi, TWO_PI * det)
                     @ equatorial_rotation_unitary(p1, ang.theta_c))
            else:
                u = _detuning_sandwich(0.0, theta, TWO_PI * det)
            tgt = equatorial_rotation_unitary(0.0, theta)
            prod = u @ tgt.conj().T
            ph = prod.flat[int(np.argmax(np.abs(prod)))]
            ph = ph / abs(ph)
            infs.append(float(np.abs(prod / ph - np.eye(2)).max() ** 2 / 2))
        scan["det_infidelity"][family] = np.array(infs)
    plot_records.append(viz.plot_scrofulous_robustness(
        scan, output_dir / "scrofulous_robustness.png"))
    # 图4：Ramsey/echo/XY4/XY16（static pool 上做变时长演化）
    static = pools["static_idle"]
    coherence = {}
    rng = np.random.default_rng(seed_for("channel_validation_pool", 5))
    durations_ms = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 5.5])
    for dd in ("none", "spin_echo", "xy4", "xy16"):
        contrast = []
        for dur_ms in durations_ms:
            n_sub_pool = min(512, static.n_traj)
            sub = _truncate_pool(static, dur_ms * 1e-3)
            from .spin_channels import MoveChannelEngine as MCE
            eng = MCE(sub, omega, dd_mode=dd, pulse_family="bare")
            psi = np.tile(INPUT_STATES["x_plus"], (n_sub_pool, 1)
                          ).astype(complex)
            k = rng.integers(0, static.n_traj, size=n_sub_pool)
            eng.propagate(psi, k, cfg.eta_slm, cfg.eta_aod,
                          dt_free_us=cfg.dt_free_us,
                          dt_pulse_us=cfg.dt_pulse_us)
            bloch = statevector_to_bloch(psi)
            contrast.append(float(np.abs(np.mean(bloch[:, 0]
                                                 + 1j * bloch[:, 1]))))
        coherence[dd] = (durations_ms.tolist(), contrast)
    plot_records.append(viz.plot_ramsey_echo_dd(coherence,
                                                output_dir
                                                / "ramsey_echo_dd.png"))
    # 图12：与 Level 4 对比（repeated moves 1..10）
    pool = pools["protocol_a"]
    l4 = json.loads(Path(l4_metrics_path).read_text(encoding="utf-8"))
    l4_curve = {}
    csv_path = l4_metrics_path.parent / "repeated_round_metrics.csv"
    if csv_path.is_file():
        import csv as _csv
        rr = [r for r in _csv.DictReader(csv_path.open(encoding="utf-8"))
              if r["protocol"] == "protocol_a"]
        rounds = [int(r["round"]) for r in rr]
        for key, col in (("contrast_none", "contrast_none"),
                         ("contrast_xy4", "contrast_xy4"),
                         ("contrast_spin_echo", "contrast_spin_echo")):
            l4_curve[key] = (rounds, [float(r[col]) for r in rr])
    l5_curve = {}
    rng = np.random.default_rng(seed_for("channel_validation_pool", 7))
    max_rounds = min(10, len(l4_curve.get("contrast_xy4", ([], []))[0])
                     or 10)
    for dd in ("none", "xy4_per_long_move"):
        rounds_l, contrasts = [], []
        for rnd in range(1, max_rounds + 1):
            n_sub = 512
            psi = np.tile(INPUT_STATES["x_plus"], (n_sub, 1)).astype(complex)
            surv = np.ones(n_sub, dtype=bool)
            for _ in range(rnd):
                k = rng.integers(0, pool.n_traj, size=n_sub)
                eng = MoveChannelEngine(pool, omega, dd_mode=dd,
                                        pulse_family="bare")
                eng.propagate(psi, k, cfg.eta_slm, cfg.eta_aod,
                              dt_free_us=cfg.dt_free_us,
                              dt_pulse_us=cfg.dt_pulse_us)
                surv &= pool.roundtrip_success[k]
            bloch = statevector_to_bloch(psi)
            contrasts.append(float(np.abs(np.mean(
                bloch[surv, 0] + 1j * bloch[surv, 1])))
                if surv.any() else 0.0)
            rounds_l.append(rnd)
        l5_curve[dd] = (rounds_l, contrasts)
    plot_records.append(viz.plot_coherence_vs_round(
        l4_curve, l5_curve, output_dir / "coherence_vs_round.png"))
    return {"l4_curve": {k: {"rounds": v[0], "contrast": v[1]}
                         for k, v in l4_curve.items()},
            "l5_curve": {k: {"rounds": v[0], "contrast": v[1]}
                         for k, v in l5_curve.items()}}


def _detuning_sandwich(phase, angle, delta_total):
    """含总失谐的单旋转（子步数值传播，作图用）。"""
    from .quantum_propagators import su2_propagator
    n_sub = 64
    omega = 2 * np.pi * 24.611e3
    dt = (angle / omega) / n_sub
    u = np.eye(2, dtype=complex)
    for _ in range(n_sub):
        u00, u01, u10, u11 = su2_propagator(np.array([delta_total]),
                                            np.array([omega]), phase, dt)
        step = np.array([[u00[0], u01[0]], [u10[0], u11[0]]],
                        dtype=complex)
        u = step @ u
    return u


def _truncate_pool(pool: TrajectoryPool, duration_s: float):
    """把池截断到前 duration_s（静态相干性扫描用）。"""
    import dataclasses
    n_t = int(np.searchsorted(pool.times_s, duration_s, side="right"))
    n_t = max(2, min(n_t, pool.times_s.size))
    return dataclasses.replace(pool, times_s=pool.times_s[:n_t],
                               u_slm=pool.u_slm[:, :n_t],
                               u_aod=pool.u_aod[:, :n_t],
                               total_duration_s=float(pool.times_s[n_t - 1]),
                               segment_boundaries_s=np.array(
                                   [0.0, float(pool.times_s[n_t - 1])]),
                               segment_names=["static_idle_trunc"])


def _write_rb_sequences_jsonl(cfg, group, output_dir):
    """rb_sequences.jsonl：reference 与 IRB 序列（可重放）。"""
    ref = generate_reference_sequences(group, cfg.rb_lengths,
                                       cfg.rb_sequences_per_length,
                                       cfg.rb_sequence_seed)
    irb = generate_irb_sequences(group, cfg.irb_total_cliffords,
                                 cfg.irb_interleaved_counts,
                                 cfg.irb_sequences_per_count,
                                 cfg.irb_sequence_seed)
    with (output_dir / "rb_sequences.jsonl").open("w",
                                                  encoding="utf-8") as fh:
        for s in ref:
            fh.write(json.dumps({"kind": "reference_rb", **s},
                                ensure_ascii=False) + "\n")
        for s in irb:
            fh.write(json.dumps({"kind": "interleaved_rb",
                                 "operation": cfg.irb_operation, **s},
                                ensure_ascii=False) + "\n")
    return len(ref) + len(irb)


def _frozen_design(cfg) -> dict:
    """frozen_rb_design.yaml 内容（打开正式结果前冻结）。"""
    return {
        "frozen_at": "before_formal_rb",
        "reference_rb": {"lengths": list(cfg.rb_lengths),
                         "sequences_per_length":
                             cfg.rb_sequences_per_length,
                         "seed": cfg.rb_sequence_seed},
        "interleaved_rb": {"N": cfg.irb_total_cliffords,
                           "M_values": list(cfg.irb_interleaved_counts),
                           "sequences_per_count":
                               cfg.irb_sequences_per_count,
                           "operation": cfg.irb_operation,
                           "seed": cfg.irb_sequence_seed},
        "sites": {"validation": cfg.validation_sites,
                  "extension": cfg.extension_sites, "seed": cfg.site_seed},
        "fit_windows": {"reference": "全部预注册长度",
                        "irb": "全部预注册 M"},
        "noise_realizations": {"reference": cfg.rb_noise_realizations,
                               "interleaved": cfg.irb_noise_realizations},
        "seed_pools": "reference_rb_pool / interleaved_rb_pool / "
                      "site_validation_pool 与 calibration 池分离",
        "note": "设计先于数据冻结；不得依据 formal 结果调参重跑"}


def _build_summary(cfg, preflight, channel, ref, irb, arr, sens,
                   clifford_stats, coherence_cmp) -> str:
    """生成 summary.md。"""
    lines = [
        "# Level 5 有限脉冲单比特动力学、噪声谱与多站点 RB/IRB", "",
        "在 Level 4 冻结的经典轨迹与协议之上，把内部态升级为旋转框架下"
        "有限时长微波脉冲驱动的 ^133Cs 钟态量子比特，完成单通道 PTM 表征、"
        "reference RB、transport-interleaved RB 与多站点统计。", "",
        f"- Rabi 锚点 Ω₀ = 2π×{cfg.drive_rabi_frequency_khz} kHz"
        f"（论文 Fig. 4a）；η_SLM=η_AOD={cfg.eta_slm}；"
        f"nominal 噪声{'开' if cfg.nominal_noise_enabled else '关'}。",
        f"- 脉冲族 {cfg.default_family}（Clifford 门）/ bare（DD π 脉冲）；"
        f"DD transport 模式 {cfg.transport_mode}。",
        f"- Clifford 群 24 元素；平均 SCROFULOUS 面积 "
        f"{clifford_stats['avg_scrofulous_area_pi_per_clifford']:.3f}π"
        f"（论文平均锚点 2.02π，分解不同不作归一）。",
    ]
    if preflight:
        lines.append(f"- preflight：{len(preflight['critical_checks'])} 项"
                     f"{'全部通过' if preflight['all_passed'] else '存在失败'}"
                     f"（{preflight.get('elapsed_s_total', 0):.0f} s）。")
    if channel:
        lines += ["", "## 单通道表征（survival-conditioned PTM）", "",
                  "| 条件 | DD | survival | F_avg(cond) | 95% CI | F_avg"
                  "(joint) |", "|---|---|---:|---:|---|---:|"]
        for c in channel["rows"]:
            ci = f"[{c['f_avg_ci_low']:.5f}, {c['f_avg_ci_high']:.5f}]" \
                if c.get("f_avg_ci_low") else "—"
            lines.append(
                f"| {c['pool']} | {c['dd']} | {c['survival']:.4f} | "
                f"{c['f_avg_conditional']:.6f} | {ci} | "
                f"{c['f_avg_survival_weighted']:.6f} |")
    if ref:
        formal = ref["summaries"].get("formal_47")
        if formal:
            lines += ["", "## Reference RB（47 sites × 72 序列）", "",
                      f"- 单指数 p = {formal['p']}；"
                      f"首选模型 {formal['selection'].get('preferred')}；"
                      f"早期斜率 "
                      f"{formal['early_slope'].get('log_slope')}。"]
    if irb:
        s = irb["summary"]
        est = s.get("irb_estimate", {})
        lines += ["", "## Interleaved RB（N=80，47 sites × 72 strings）", ""]
        lines.append(f"- IRB 判定：{s['assumptions']['verdict']}")
        if est.get("interpretable"):
            lines.append(f"- 标准 IRB 每移动 F_avg = "
                         f"{est['f_avg_interleaved']:.6f}"
                         f"（p_ref={est.get('p_ref'):.6f}, "
                         f"p_comb={est.get('p_interleaved_combined'):.6f}）")
        else:
            lines.append(f"- 不输出标准 IRB 保真度：{est.get('reason', '假设未满足')}。")
        lines.append(f"- survival（M=80）= "
                     f"{s['survival']['means'][-1]:.4f}；conditional（M=80"
                     f"）= {s['conditional_return']['means'][-1]:.4f}；"
                     f"joint（M=80）= {s['joint_return']['means'][-1]:.4f}。")
        lines.append(f"- 相邻比值稳定性（conditional）："
                     f"{s['instantaneous_ratios']['conditional'].get('stable')}"
                     f"；尾部 risk set = "
                     f"{s['instantaneous_ratios']['survival'].get('min_risk_set')}"
                     f"。")
    if arr:
        for tag, st in arr["summaries"].items():
            sc = st["site_stats"]["conditional"]
            lines += ["", f"## {tag}", "",
                      f"- array conditional = {sc['array_average']:.5f}；"
                      f"between-site std = {sc['between_site_std']:.5f}；"
                      f"within-site std = {sc['within_site_std']:.5f}；"
                      f"worst decile = {sc['worst_decile']:.5f}。",
                      f"- 站点参数相关性：{sc['parameter_correlations']}。"]
    if sens:
        lines += ["", "## 敏感性（one-factor-at-a-time，assumed）", ""]
        for r in sens["rows"]:
            lines.append(f"- {r['factor']} = {r['value']}: "
                         f"F_avg(cond)={r['f_avg_conditional']:.6f}")
    if coherence_cmp:
        lines += ["", "## 与 Level 4 半经典模型对照", ""]
        for k, v in coherence_cmp.get("l4_curve", {}).items():
            lines.append(f"- Level 4 {k}: 10 轮对比度 "
                         f"{v['contrast'][-1]:.4f}。")
        for k, v in coherence_cmp.get("l5_curve", {}).items():
            lines.append(f"- Level 5 {k}: {len(v['rounds'])} 轮有限脉冲"
                         f"对比度 {v['contrast'][-1]:.4f}。")
    lines += ["", "## 与论文的关系与边界", "",
              "- 本级所有数值为模拟量：survival 是经典轨迹标签的经验均值，"
              "conditional spin channel 是冻结轨迹上的幺正系综平均；两者"
              "都**不是**论文实验 IRB fidelity。",
              "- 论文的 transformed Clifford frame XY4 未启用（配置"
              "transformed_clifford_frame=false）；本实现的 xy4_per_long_move "
              "是不同的、更简单的模式。",
              "- nominal 噪声关闭；噪声敏感性参数全部标记 assumed。"
              "paper-coherence anchor 仅作对照，不用于校准。",
              "- 47/195 站点是条件独立单原子 + 可配置共同/局域经典噪声；"
              "不包含原子间相互作用或并行门串扰。",
              "- 经典轨迹池冻结重放：站点阱深不均匀只通过势能线性缩放进入"
              "量子通道，经典 survival 未按站点重新生成。",
              "- 未包含：Rydberg 相互作用、两比特门、测量/SPAM 误差、"
              "多能级泄漏与自发 Raman 散射。", "",
              "图片完成程序化检查；未进行人工视觉审阅。", ""]
    return "\n".join(lines)


# ------------------------------------------------------------------ main
def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    cfg = load_config(Path(args.config))
    output_dir = ensure_output_directory(
        Path(args.output_dir) if args.output_dir else
        cfg.project_root / cfg.output_directory, cfg.project_root)
    l4_dir = _discover_level4_output(cfg)
    cfg4 = load_level4_config(cfg.project_root / cfg.level4_config_path)
    _, resolved = l4_upstream.build_manifest(cfg4)
    protocols = {
        "protocol_a": ps.build_protocol_a(cfg4,
                                          resolved["nominal_vs_m_per_s"]),
        "protocol_a_no_lensing": ps.build_protocol_a(cfg4, None),
        "protocol_b": ps.build_protocol_b(cfg4, resolved["level2_spec"],
                                          resolved["nominal_vs_m_per_s"]),
        "static_idle": ps.build_static_idle(cfg4),
    }
    save_yaml(output_dir / "config_used.yaml", config_as_dict(cfg))
    frozen_path = output_dir / "frozen_rb_design.yaml"
    if not frozen_path.is_file():
        save_yaml(frozen_path, _frozen_design(cfg))
    group = CliffordGroup()
    pools = _ensure_pools(cfg, cfg4, protocols, output_dir)
    manifest = _build_upstream_manifest(cfg, l4_dir, {
        name: {"_path": output_dir / "frozen_pools" / f"{name}.npz",
               "_sha256": pool.sha256()}
        for name, pool in pools.items()})
    save_json(output_dir / "upstream_manifest.json", manifest)
    stage = args.stage
    results = {"preflight": None, "channel": None, "reference_rb": None,
               "interleaved_rb": None, "array_extension": None,
               "sensitivity": None}
    plot_records = []
    clifford_stats = _write_clifford_table_csv(cfg, group, output_dir)
    _write_pulse_library_csv(cfg, output_dir)
    _write_noise_artifacts(cfg, output_dir, plot_records)
    _write_rb_sequences_jsonl(cfg, group, output_dir)
    started_all = time.time()
    if stage in ("preflight", "all"):
        pf = Preflight5(cfg, manifest, resolved, group=group,
                        pool=pools["protocol_a"],
                        static_pool=pools["static_idle"],
                        timeline=protocols["protocol_a"], cfg4=cfg4,
                        run_upstream_tests=(cfg.run_upstream_tests_in_preflight
                                            and not args.fast_preflight),
                        log=print)
        checks = pf.run()
        save_json(output_dir / "preflight.json", checks)
        results["preflight"] = checks
        if not checks["all_passed"]:
            failed = [n for n in checks["critical_checks"]
                      if not checks["checks"][n]["passed"]]
            raise RuntimeError(f"preflight 未通过，停止正式 RB：{failed}")
    if stage in ("channel", "all"):
        results["channel"] = stage_channel(cfg, pools, output_dir,
                                           plot_records)
    if stage in ("reference_rb", "all"):
        results["reference_rb"] = stage_reference_rb(cfg, pools, output_dir,
                                                     plot_records)
    if stage in ("interleaved_rb", "all"):
        results["interleaved_rb"] = stage_interleaved_rb(cfg, pools,
                                                         output_dir,
                                                         plot_records)
    if stage in ("array_extension", "all"):
        irb = results["interleaved_rb"]
        results["array_extension"] = stage_array_extension(
            cfg, pools, output_dir, plot_records,
            formal_rows_47=(irb["formal_rows"] if irb else None),
            ens_47=(irb["site_ensemble_full"] if irb else None))
    if stage in ("sensitivity", "all"):
        results["sensitivity"] = stage_sensitivity(cfg, pools, output_dir,
                                                   plot_records)
    # 诊断图（不依赖阶段）
    coherence_cmp = _diagnostic_plots(cfg, pools, group, output_dir,
                                      plot_records,
                                      l4_dir / "metrics.json")
    pf_json = output_dir / "preflight.json"
    if pf_json.is_file() and results["preflight"] is None:
        results["preflight"] = json.loads(pf_json.read_text(
            encoding="utf-8"))
    if results["preflight"] and results["preflight"].get(
            "dt_convergence", {}).get("rows"):
        plot_records.append(viz.plot_timestep_convergence(
            results["preflight"]["dt_convergence"]["rows"],
            output_dir / "timestep_convergence.png"))
    image_validation = validate_pngs([r for r in plot_records if r])
    save_json(output_dir / "image_validation.json", image_validation)
    metrics = {
        "level5": {"name": cfg.name, "stage": stage,
                   "elapsed_s_total": round(time.time() - started_all, 1)},
        "upstream": {"manifest_files": len(manifest["files"]),
                     "trajectory_mode": cfg.trajectory_mode,
                     "level4_output_dir": str(l4_dir)},
        "preflight": {"all_passed": (results["preflight"] or {}).get(
            "all_passed")},
        "pulse_clifford_noise": {
            "rabi_khz": cfg.drive_rabi_frequency_khz,
            "family": cfg.default_family, "dd_transport": cfg.transport_mode,
            "eta_slm": cfg.eta_slm, "eta_aod": cfg.eta_aod,
            "nominal_noise_enabled": cfg.nominal_noise_enabled,
            "clifford_pulse_histogram":
                clifford_stats["pulse_count_histogram"],
            "avg_scrofulous_area_pi_per_clifford":
                clifford_stats["avg_scrofulous_area_pi_per_clifford"],
            "paper_anchor_area_pi": 2.02},
        "channel": None if not results["channel"] else {
            "conditions": {c: {"survival": v.get("survival"),
                               "f_avg": v["f_avg"],
                               "f_avg_survival_weighted":
                                   v.get("f_avg_survival_weighted"),
                               "bootstrap": v["bootstrap"]}
                           for c, v in results["channel"]["metrics"]
                           .items()}},
        "reference_rb": None if not results["reference_rb"] else {
            "formal_p": results["reference_rb"]["summaries"][
                "formal_47"]["p"],
            "preferred_model": results["reference_rb"]["summaries"][
                "formal_47"]["selection"].get("preferred")},
        "interleaved_rb": None if not results["interleaved_rb"] else {
            "verdict": results["interleaved_rb"]["summary"]["assumptions"][
                "verdict"],
            "irb_estimate": results["interleaved_rb"]["summary"][
                "irb_estimate"],
            "survival_M80": results["interleaved_rb"]["summary"][
                "survival"]["means"][-1],
            "conditional_M80": results["interleaved_rb"]["summary"][
                "conditional_return"]["means"][-1],
            "joint_M80": results["interleaved_rb"]["summary"][
                "joint_return"]["means"][-1]},
        "array": None if not results["array_extension"] else {
            tag: st["site_stats"]["conditional"]
            for tag, st in results["array_extension"]["summaries"].items()},
        "paper_coherence_anchor": {
            "enabled": cfg.anchor_enabled,
            "ensemble_t2star_ms": cfg.anchor_ensemble_t2star_ms,
            "site_t2star_ms": cfg.anchor_site_t2star_ms,
            "xy16_t2_s": cfg.anchor_xy16_t2_s,
            "used_for_calibration": cfg.anchor_use_for_calibration},
        "coherence_level4_vs_level5": coherence_cmp,
        "convergence_ablation_sensitivity": {
            "dt_convergence_rows": (results["preflight"] or {}).get(
                "dt_convergence", {}).get("rows", []),
            "sensitivity_rows": (results["sensitivity"] or {}).get(
                "rows", [])},
        "provenance": manifest.get("provenance", {}),
        "image_validation_all_passed": image_validation["all_passed"],
        "visual_review_performed": False,
        "level6_recommendations": [
            "Rydberg 两比特门与阻塞丢失", "多体串扰与并行门",
            "测量/SPAM 模型与读出误差", "泄漏态（|F,mF⟩≠0）与 Raman 散射",
            "逻辑层电路与纠错周期"],
    }
    save_json(output_dir / "metrics.json", metrics)
    save_markdown(output_dir / "summary.md", _build_summary(
        cfg, results["preflight"], results["channel"],
        results["reference_rb"], results["interleaved_rb"],
        results["array_extension"], results["sensitivity"],
        clifford_stats, coherence_cmp))
    print(f"Level 5 stage '{stage}' completed: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
