"""Level 5 配置加载与校验（沿用 Level 0-4 dataclass 风格）。"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class Level5Config:
    """解析后的 Level 5 配置。"""

    name: str
    # 上游
    upstream_level4_output: str
    require_sha256_manifest: bool
    trajectory_mode: str
    joint_validation_shots: int
    # 量子比特
    drive_rabi_frequency_khz: float
    drive_detuning_hz: float
    eta_slm: float
    eta_aod: float
    eta_sensitivity: tuple
    rotating_wave_approximation: bool
    include_state_dependent_force: bool
    include_bloch_siegert: bool
    # 量子积分
    method: str
    max_dt_us: float
    convergence_dt_us: tuple
    dt_free_us: float
    dt_pulse_us: float
    norm_tolerance: float
    density_matrix_min_eigenvalue_tolerance: float
    n_sub_pulse: int
    # 脉冲
    default_family: str
    envelope: str
    cosine_edge_us: float
    include_amplitude_error: bool
    amplitude_error_rms: float
    include_phase_error: bool
    # DD
    dd_modes: tuple
    transport_mode: str
    transformed_clifford_frame: bool
    # 噪声
    nominal_noise_enabled: bool
    detuning_model: str
    amplitude_model: str
    psd_model: str
    psd_spec: list
    quasistatic_detuning_rms_hz: float
    ou_rms_hz: float
    ou_tau_us: float
    common_fraction: float
    noise_seed: int
    # paper anchor
    anchor_enabled: bool
    anchor_trap_depth_uK: float
    anchor_temperature_uK: float
    anchor_ensemble_t2star_ms: float
    anchor_site_t2star_ms: float
    anchor_xy16_t2_s: float
    anchor_use_for_calibration: bool
    # channel 表征
    channel_trajectories: int
    channel_dev_trajectories: int
    bootstrap_samples: int
    input_states: tuple
    # reference RB
    rb_lengths: tuple
    rb_sequences_per_length: int
    rb_dev_sequences_per_length: int
    rb_sequence_seed: int
    # IRB
    irb_total_cliffords: int
    irb_interleaved_counts: tuple
    irb_sequences_per_count: int
    irb_dev_sequences_per_count: int
    irb_operation: str
    irb_sequence_seed: int
    irb_noise_realizations: int
    rb_noise_realizations: int
    # 阵列
    development_sites: int
    validation_sites: int
    extension_sites: int
    site_seed: int
    trap_depth_relative_std: float
    rabi_relative_std: float
    static_detuning_std_Hz: float
    eta_relative_std: float
    spatial_noise_correlation: float
    paper_scenario_trap_depth_std: float
    # 轨迹池
    pool_seed: int
    pool_chunk_shots: int
    pool_dt_classical_us: float
    pool_record_stride_us: float
    pool_full_trajectories: int
    pool_secondary_trajectories: int
    # 统计
    statistics_bootstrap_samples: int
    confidence_level: float
    # 输出
    output_directory: str
    save_sequence_resolved_data: bool
    save_site_resolved_data: bool
    save_noise_realizations_for_replay: bool
    # 敏感性
    sensitivity_trajectories: int
    # 其他
    project_root: Path
    config_path: Path
    level4_config_path: str
    run_upstream_tests_in_preflight: bool


def _root(path: Path) -> Path:
    for candidate in (path.parent, *path.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate.resolve()
    cwd = Path.cwd().resolve()
    if (cwd / "pyproject.toml").is_file():
        return cwd
    return path.parent.resolve()


def load_config(path: Path) -> Level5Config:
    """读取并校验 Level 5 YAML 配置。"""
    path = Path(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    up = data["upstream"]
    q = data["qubit"]
    qi = data["quantum_integration"]
    pu = data["pulses"]
    dd = data["dd"]
    no = data["noise"]
    an = data.get("paper_coherence_anchor", {})
    ch = data.get("channel_characterization", {})
    rb = data["reference_rb"]
    ir = data["interleaved_rb"]
    ar = data["array"]
    st = data.get("statistics", {})
    ou = data["output"]
    po = data.get("trajectory_pool", {})
    se = data.get("sensitivity", {})
    if up.get("trajectory_mode") not in ("frozen_replay", "joint"):
        raise ValueError("trajectory_mode 仅支持 frozen_replay/joint")
    if ir["operation"] not in ("one_way_transport", "aod_roundtrip",
                              "full_pickup_transport_dropoff_roundtrip",
                              "static_idle_matched_duration"):
        raise ValueError(f"未知 interleaved operation {ir['operation']}")
    cfg = Level5Config(
        name=data["level5"]["name"],
        upstream_level4_output=str(up.get("level4_output", "auto_discover")),
        require_sha256_manifest=bool(up.get("require_sha256_manifest", True)),
        trajectory_mode=str(up["trajectory_mode"]),
        joint_validation_shots=int(up.get("joint_validation_shots", 128)),
        drive_rabi_frequency_khz=float(q["drive_rabi_frequency_kHz"]),
        drive_detuning_hz=float(q.get("drive_detuning_Hz", 0.0)),
        eta_slm=float(q["eta_slm"]),
        eta_aod=float(q["eta_aod"]),
        eta_sensitivity=tuple(float(v) for v in q.get("eta_sensitivity", ())),
        rotating_wave_approximation=bool(q.get(
            "rotating_wave_approximation", True)),
        include_state_dependent_force=bool(q.get(
            "include_state_dependent_force", False)),
        include_bloch_siegert=bool(q.get("include_bloch_siegert", False)),
        method=str(qi.get("method", "analytic_su2_midpoint")),
        max_dt_us=float(qi.get("max_dt_us", 0.1)),
        convergence_dt_us=tuple(float(v) for v in qi.get(
            "convergence_dt_us", (0.10, 0.05, 0.025))),
        dt_free_us=float(qi.get("dt_free_us", 0.25)),
        dt_pulse_us=float(qi.get("dt_pulse_us", 0.1)),
        norm_tolerance=float(qi.get("norm_tolerance", 1e-10)),
        density_matrix_min_eigenvalue_tolerance=float(qi.get(
            "density_matrix_min_eigenvalue_tolerance", -1e-10)),
        n_sub_pulse=int(qi.get("n_sub_pulse", 4)),
        default_family=str(pu.get("default_family", "scrofulous")),
        envelope=str(pu.get("envelope", "square")),
        cosine_edge_us=float(pu.get("cosine_edge_us", 0.0)),
        include_amplitude_error=bool(pu.get("include_amplitude_error", False)),
        amplitude_error_rms=float(pu.get("amplitude_error_rms", 0.0)),
        include_phase_error=bool(pu.get("include_phase_error", False)),
        dd_modes=tuple(str(v) for v in dd.get("modes", ())),
        transport_mode=str(dd.get("transport_mode", "xy4_per_long_move")),
        transformed_clifford_frame=bool(dd.get(
            "transformed_clifford_frame", False)),
        nominal_noise_enabled=bool(no.get("nominal_enabled", False)),
        detuning_model=str(no.get("detuning_model", "none")),
        amplitude_model=str(no.get("amplitude_model", "none")),
        psd_model=str(no.get("psd_model", "none")),
        psd_spec=list(no.get("psd_spec", [])),
        quasistatic_detuning_rms_hz=float(no.get(
            "quasistatic_detuning_rms_hz", 0.0)),
        ou_rms_hz=float(no.get("ou_rms_hz", 0.0)),
        ou_tau_us=float(no.get("ou_tau_us", 200.0)),
        common_fraction=float(no.get("common_fraction", 0.0)),
        noise_seed=int(no.get("noise_seed", 51999)),
        anchor_enabled=bool(an.get("enabled", True)),
        anchor_trap_depth_uK=float(an.get("trap_depth_uK", 55.0)),
        anchor_temperature_uK=float(an.get("temperature_uK", 4.3)),
        anchor_ensemble_t2star_ms=float(an.get("ensemble_t2star_ms", 14.0)),
        anchor_site_t2star_ms=float(an.get("site_t2star_ms", 25.5)),
        anchor_xy16_t2_s=float(an.get("xy16_t2_s", 12.6)),
        anchor_use_for_calibration=bool(an.get("use_for_calibration",
                                               False)),
        channel_trajectories=int(ch.get("classical_trajectories", 2000)),
        channel_dev_trajectories=int(ch.get("dev_trajectories", 256)),
        bootstrap_samples=int(ch.get("bootstrap_samples", 2000)),
        input_states=tuple(str(v) for v in ch.get(
            "input_states", ("z_plus", "z_minus", "x_plus", "x_minus",
                             "y_plus", "y_minus"))),
        rb_lengths=tuple(int(v) for v in rb["lengths"]),
        rb_sequences_per_length=int(rb["sequences_per_length"]),
        rb_dev_sequences_per_length=int(rb.get(
            "dev_sequences_per_length", 8)),
        rb_sequence_seed=int(rb["sequence_seed"]),
        irb_total_cliffords=int(ir["total_cliffords"]),
        irb_interleaved_counts=tuple(int(v) for v in ir["interleaved_counts"]),
        irb_sequences_per_count=int(ir["sequences_per_count"]),
        irb_dev_sequences_per_count=int(ir.get(
            "dev_sequences_per_count", 8)),
        irb_operation=str(ir["operation"]),
        irb_sequence_seed=int(ir["sequence_seed"]),
        irb_noise_realizations=int(ir.get("noise_realizations", 1)),
        rb_noise_realizations=int(rb.get("noise_realizations", 1)),
        development_sites=int(ar["development_sites"]),
        validation_sites=int(ar["validation_sites"]),
        extension_sites=int(ar["extension_sites"]),
        site_seed=int(ar["site_seed"]),
        trap_depth_relative_std=float(ar.get("trap_depth_relative_std", 0.0)),
        rabi_relative_std=float(ar.get("rabi_relative_std", 0.0)),
        static_detuning_std_Hz=float(ar.get("static_detuning_std_Hz", 0.0)),
        eta_relative_std=float(ar.get("eta_relative_std", 0.0)),
        spatial_noise_correlation=float(ar.get(
            "spatial_noise_correlation", 0.0)),
        paper_scenario_trap_depth_std=float(ar.get(
            "paper_scenario_trap_depth_relative_std", 0.114)),
        pool_seed=int(po.get("seed", 50900)),
        pool_chunk_shots=int(po.get("chunk_shots", 250)),
        pool_dt_classical_us=float(po.get("dt_classical_us", 0.1)),
        pool_record_stride_us=float(po.get("record_stride_us", 0.5)),
        pool_full_trajectories=int(po.get("full_pool_trajectories", 2000)),
        pool_secondary_trajectories=int(po.get(
            "secondary_pool_trajectories", 1000)),
        statistics_bootstrap_samples=int(st.get("bootstrap_samples", 2000)),
        confidence_level=float(st.get("confidence_level", 0.95)),
        output_directory=str(ou["directory"]),
        save_sequence_resolved_data=bool(ou.get(
            "save_sequence_resolved_data", True)),
        save_site_resolved_data=bool(ou.get("save_site_resolved_data",
                                            True)),
        save_noise_realizations_for_replay=bool(ou.get(
            "save_noise_realizations_for_replay", True)),
        sensitivity_trajectories=int(se.get("trajectories", 500)),
        project_root=_root(path), config_path=path,
        level4_config_path=str(up.get("level4_config",
                                      "configs/level4_roundtrip.yaml")),
        run_upstream_tests_in_preflight=bool(
            up.get("run_upstream_tests_in_preflight", True)),
    )
    return cfg


def config_as_dict(cfg: Level5Config) -> dict:
    """配置转字典（写入 config_used.yaml）。"""
    return {
        "level5": {"name": cfg.name},
        "upstream": {
            "level4_output": cfg.upstream_level4_output,
            "require_sha256_manifest": cfg.require_sha256_manifest,
            "trajectory_mode": cfg.trajectory_mode,
            "joint_validation_shots": cfg.joint_validation_shots,
            "level4_config": cfg.level4_config_path,
            "run_upstream_tests_in_preflight":
                cfg.run_upstream_tests_in_preflight},
        "qubit": {
            "basis": "cs133_clock_states",
            "drive_rabi_frequency_kHz": cfg.drive_rabi_frequency_khz,
            "drive_detuning_Hz": cfg.drive_detuning_hz,
            "eta_slm": cfg.eta_slm, "eta_aod": cfg.eta_aod,
            "eta_sensitivity": list(cfg.eta_sensitivity),
            "rotating_wave_approximation":
                cfg.rotating_wave_approximation,
            "include_state_dependent_force":
                cfg.include_state_dependent_force,
            "include_bloch_siegert": cfg.include_bloch_siegert},
        "quantum_integration": {
            "method": cfg.method, "max_dt_us": cfg.max_dt_us,
            "convergence_dt_us": list(cfg.convergence_dt_us),
            "dt_free_us": cfg.dt_free_us, "dt_pulse_us": cfg.dt_pulse_us,
            "norm_tolerance": cfg.norm_tolerance,
            "density_matrix_min_eigenvalue_tolerance":
                cfg.density_matrix_min_eigenvalue_tolerance,
            "n_sub_pulse": cfg.n_sub_pulse},
        "pulses": {
            "default_family": cfg.default_family, "envelope": cfg.envelope,
            "cosine_edge_us": cfg.cosine_edge_us,
            "include_amplitude_error": cfg.include_amplitude_error,
            "amplitude_error_rms": cfg.amplitude_error_rms,
            "include_phase_error": cfg.include_phase_error},
        "dd": {"modes": list(cfg.dd_modes),
               "transport_mode": cfg.transport_mode,
               "transformed_clifford_frame":
                   cfg.transformed_clifford_frame,
               "pulses_are_finite": True},
        "noise": {
            "nominal_enabled": cfg.nominal_noise_enabled,
            "detuning_model": cfg.detuning_model,
            "amplitude_model": cfg.amplitude_model,
            "psd_model": cfg.psd_model, "psd_spec": cfg.psd_spec,
            "quasistatic_detuning_rms_hz": cfg.quasistatic_detuning_rms_hz,
            "ou_rms_hz": cfg.ou_rms_hz, "ou_tau_us": cfg.ou_tau_us,
            "common_fraction": cfg.common_fraction,
            "noise_seed": cfg.noise_seed,
            "provenance": "assumed_unless_calibrated"},
        "paper_coherence_anchor": {
            "enabled": cfg.anchor_enabled,
            "trap_depth_uK": cfg.anchor_trap_depth_uK,
            "temperature_uK": cfg.anchor_temperature_uK,
            "ensemble_t2star_ms": cfg.anchor_ensemble_t2star_ms,
            "site_t2star_ms": cfg.anchor_site_t2star_ms,
            "xy16_t2_s": cfg.anchor_xy16_t2_s,
            "use_for_calibration": cfg.anchor_use_for_calibration,
            "note": "独立 paper-anchor 对照配置；与 Level 4 transport 协议"
                    "参数不同，不用于拟合 nominal 自由参数"},
        "channel_characterization": {
            "classical_trajectories": cfg.channel_trajectories,
            "dev_trajectories": cfg.channel_dev_trajectories,
            "bootstrap_samples": cfg.bootstrap_samples,
            "input_states": list(cfg.input_states)},
        "reference_rb": {
            "lengths": list(cfg.rb_lengths),
            "sequences_per_length": cfg.rb_sequences_per_length,
            "dev_sequences_per_length": cfg.rb_dev_sequences_per_length,
            "sequence_seed": cfg.rb_sequence_seed,
            "noise_realizations": cfg.rb_noise_realizations},
        "interleaved_rb": {
            "total_cliffords": cfg.irb_total_cliffords,
            "interleaved_counts": list(cfg.irb_interleaved_counts),
            "sequences_per_count": cfg.irb_sequences_per_count,
            "dev_sequences_per_count": cfg.irb_dev_sequences_per_count,
            "operation": cfg.irb_operation,
            "sequence_seed": cfg.irb_sequence_seed,
            "noise_realizations": cfg.irb_noise_realizations},
        "array": {
            "development_sites": cfg.development_sites,
            "validation_sites": cfg.validation_sites,
            "extension_sites": cfg.extension_sites,
            "site_seed": cfg.site_seed,
            "trap_depth_relative_std": cfg.trap_depth_relative_std,
            "rabi_relative_std": cfg.rabi_relative_std,
            "static_detuning_std_Hz": cfg.static_detuning_std_Hz,
            "eta_relative_std": cfg.eta_relative_std,
            "spatial_noise_correlation": cfg.spatial_noise_correlation,
            "paper_scenario_trap_depth_relative_std":
                cfg.paper_scenario_trap_depth_std},
        "trajectory_pool": {
            "seed": cfg.pool_seed, "chunk_shots": cfg.pool_chunk_shots,
            "dt_classical_us": cfg.pool_dt_classical_us,
            "record_stride_us": cfg.pool_record_stride_us,
            "full_pool_trajectories": cfg.pool_full_trajectories,
            "secondary_pool_trajectories":
                cfg.pool_secondary_trajectories},
        "statistics": {
            "bootstrap_samples": cfg.statistics_bootstrap_samples,
            "hierarchical_bootstrap": True,
            "confidence_level": cfg.confidence_level},
        "sensitivity": {"trajectories": cfg.sensitivity_trajectories},
        "output": {
            "directory": cfg.output_directory,
            "save_sequence_resolved_data":
                cfg.save_sequence_resolved_data,
            "save_site_resolved_data": cfg.save_site_resolved_data,
            "save_noise_realizations_for_replay":
                cfg.save_noise_realizations_for_replay},
    }
