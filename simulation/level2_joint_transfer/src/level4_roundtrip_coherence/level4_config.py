"""Level 4 配置加载与校验（沿用 Level 0-3 的 dataclass 配置风格）。"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class Level4Config:
    """解析后的 Level 4 配置。"""

    name: str
    mass_u: float
    # SLM 阱（静态三维高斯）
    slm_wavelength_nm: float
    slm_waist_um: float
    slm_initial_depth_uK: float
    slm_transport_depth_uK: float
    slm_center_um: tuple
    slm_protocol_b_depth_uK: float
    # AOD 阱（Level 3 crossed-AOD 透镜势）
    aod_wavelength_nm: float
    aod_waist_um: float
    aod_transport_depth_uK: float
    split_distance_um: float
    # 论文时序锚点（Extended Data Fig. 10e）
    initial_hold_us: float
    pickup_us: float
    split_us: float
    long_distance_um: float
    outbound_us: float
    remote_operation_hold_us: float
    inbound_us: float
    merge_us: float
    dropoff_us: float
    final_hold_us: float
    long_geometry: str
    long_profile: str
    split_profile: str
    # lensing（默认继承 Level 3 冻结场景）
    inherit_lensing_from_level3: bool
    nominal_vs_m_per_s: float | None
    lensing_axis_signs: tuple
    sensitivity_vs_m_per_s: tuple
    # 初态与样本池
    temperature_uK: float
    single_round_validation_shots: int
    repeated_round_shots: int
    long_tail_shots: int
    sensitivity_shots: int
    single_round_seed: int
    repeated_round_seed: int
    long_tail_seed: int
    sensitivity_seed: int
    development_seed: int
    # 积分
    validation_dt_us: float
    repeated_dt_us: float
    convergence_dt_us: float
    max_rounds: int
    long_tail_max_rounds: int
    checkpoint_rounds: tuple
    # 相干性
    coherence_enabled: bool
    eta_slm: float
    eta_aod: float
    eta_sensitivity: tuple
    dd_modes: tuple
    dd_pulse_offset_us: float
    custom_pulse_times_us: tuple
    # 技术噪声（仅敏感性）
    noise_enabled_nominal: bool
    noise_model: str
    noise_rms_grid: tuple
    noise_tau_us_grid: tuple
    noise_seed: int
    noise_cross_correlation: float
    # 分类
    energy_threshold_J: float
    near_threshold_fraction_of_depth: float
    ambiguous_basin_tolerance: float
    absorbing_protocol_success: bool
    # 上游冻结
    upstream_level2_waveform: str
    upstream_level2_metrics: str
    upstream_level3_metrics: str
    upstream_level3_config: str
    require_sha256_manifest: bool
    allow_reconstruction_without_reoptimization: bool
    # 输出
    output_directory: str
    save_all_checkpoint_states: bool
    save_representative_trajectories: bool
    representative_shots_per_failure_stage: int
    bootstrap_samples: int
    bootstrap_seed: int
    project_root: Path
    config_path: Path


def _root(path: Path) -> Path:
    for candidate in (path.parent, *path.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate.resolve()
    return path.parent.resolve()


def load_config(path: Path) -> Level4Config:
    """读取并校验 Level 4 YAML 配置。"""
    path = Path(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    atom = data["atom"]
    slm = data["slm"]
    aod = data["aod"]
    anchor = data["protocol_anchor"]
    lens = data.get("lensing", {})
    ens = data["initial_ensemble"]
    integ = data["integration"]
    coh = data["coherence"]
    noise = data["intensity_noise"]
    cls = data["classification"]
    up = data["upstream"]
    out = data["output"]
    seeds = {ens["single_round_seed"], ens["repeated_round_seed"],
             ens["long_tail_seed"], ens["sensitivity_seed"],
             ens.get("development_seed", 40999), noise["noise_seed"]}
    if len(seeds) != 6:
        raise ValueError("development/single_round/repeated/long_tail/sensitivity/noise 种子必须互不相同")
    if anchor["long_geometry"] not in ("diagonal", "straight"):
        raise ValueError("long_geometry 仅支持 diagonal/straight")
    cfg = Level4Config(
        name=data["level4"]["name"],
        mass_u=float(atom["mass_u"]),
        slm_wavelength_nm=float(slm["wavelength_nm"]),
        slm_waist_um=float(slm["waist_um"]),
        slm_initial_depth_uK=float(slm["initial_depth_uK"]),
        slm_transport_depth_uK=float(slm["transport_depth_uK"]),
        slm_center_um=tuple(float(v) for v in slm.get("center_um", (0.0, 0.0))),
        slm_protocol_b_depth_uK=float(slm.get("protocol_b_depth_uK", 140.0)),
        aod_wavelength_nm=float(aod["wavelength_nm"]),
        aod_waist_um=float(aod["waist_um"]),
        aod_transport_depth_uK=float(aod["transport_depth_uK"]),
        split_distance_um=float(aod["split_distance_um"]),
        initial_hold_us=float(anchor["initial_hold_us"]),
        pickup_us=float(anchor["pickup_us"]),
        split_us=float(anchor["split_us"]),
        long_distance_um=float(anchor["long_distance_um"]),
        outbound_us=float(anchor["outbound_us"]),
        remote_operation_hold_us=float(anchor["remote_operation_hold_us"]),
        inbound_us=float(anchor["inbound_us"]),
        merge_us=float(anchor["merge_us"]),
        dropoff_us=float(anchor["dropoff_us"]),
        final_hold_us=float(anchor["final_hold_us"]),
        long_geometry=str(anchor["long_geometry"]),
        long_profile=str(anchor["long_profile"]),
        split_profile=str(anchor.get("split_profile", "constant_jerk")),
        inherit_lensing_from_level3=bool(lens.get("inherit_from_level3", True)),
        nominal_vs_m_per_s=lens.get("nominal_vs_m_per_s"),
        lensing_axis_signs=tuple(int(v) for v in lens.get("axis_signs", (1, 1))),
        sensitivity_vs_m_per_s=tuple(
            float(v) for v in lens.get("sensitivity_vs_m_per_s", ())),
        temperature_uK=float(ens["temperature_uK"]),
        single_round_validation_shots=int(ens["single_round_validation_shots"]),
        repeated_round_shots=int(ens["repeated_round_shots"]),
        long_tail_shots=int(ens["long_tail_shots"]),
        sensitivity_shots=int(ens.get("sensitivity_shots", 500)),
        single_round_seed=int(ens["single_round_seed"]),
        repeated_round_seed=int(ens["repeated_round_seed"]),
        long_tail_seed=int(ens["long_tail_seed"]),
        sensitivity_seed=int(ens["sensitivity_seed"]),
        development_seed=int(ens.get("development_seed", 40999)),
        validation_dt_us=float(integ["validation_dt_us"]),
        repeated_dt_us=float(integ["repeated_dt_us"]),
        convergence_dt_us=float(integ["convergence_dt_us"]),
        max_rounds=int(integ["max_rounds"]),
        long_tail_max_rounds=int(integ["long_tail_max_rounds"]),
        checkpoint_rounds=tuple(int(v) for v in integ.get(
            "checkpoint_rounds", (1, 2, 5, 10, 20, 40))),
        coherence_enabled=bool(coh["enabled"]),
        eta_slm=float(coh["eta_slm"]),
        eta_aod=float(coh["eta_aod"]),
        eta_sensitivity=tuple(float(v) for v in coh.get("eta_sensitivity", ())),
        dd_modes=tuple(str(v) for v in coh["dd_modes"]),
        dd_pulse_offset_us=float(coh.get("dd_pulse_offset_us", 0.025)),
        custom_pulse_times_us=tuple(
            float(v) for v in coh.get("custom_pulse_times_us", ())),
        noise_enabled_nominal=bool(noise["enabled_nominal"]),
        noise_model=str(noise["model"]),
        noise_rms_grid=tuple(float(v) for v in noise.get("rms_grid", ())),
        noise_tau_us_grid=tuple(
            None if v is None else float(v)
            for v in noise.get("tau_us_grid", ())),
        noise_seed=int(noise["noise_seed"]),
        noise_cross_correlation=float(noise.get("cross_correlation", 0.0)),
        energy_threshold_J=float(cls["energy_threshold_J"]),
        near_threshold_fraction_of_depth=float(cls["near_threshold_fraction_of_depth"]),
        ambiguous_basin_tolerance=float(cls["ambiguous_basin_tolerance"]),
        absorbing_protocol_success=bool(cls.get("absorbing_protocol_success", True)),
        upstream_level2_waveform=str(up.get("level2_best_waveform", "auto_discover")),
        upstream_level2_metrics=str(up.get("level2_metrics", "auto_discover")),
        upstream_level3_metrics=str(up.get("level3_metrics", "auto_discover")),
        upstream_level3_config=str(up.get("level3_config", "auto_discover")),
        require_sha256_manifest=bool(up.get("require_sha256_manifest", True)),
        allow_reconstruction_without_reoptimization=bool(
            up.get("allow_reconstruction_without_reoptimization", True)),
        output_directory=str(out["directory"]),
        save_all_checkpoint_states=bool(out.get("save_all_checkpoint_states", True)),
        save_representative_trajectories=bool(
            out.get("save_representative_trajectories", True)),
        representative_shots_per_failure_stage=int(
            out.get("representative_shots_per_failure_stage", 3)),
        bootstrap_samples=int(out.get("bootstrap_samples", 10000)),
        bootstrap_seed=int(out.get("bootstrap_seed", 41020)),
        project_root=_root(path), config_path=path,
    )
    return cfg


def config_as_dict(cfg: Level4Config) -> dict:
    """配置转字典（用于 config_used.yaml）。"""
    return {
        "level4": {"name": cfg.name},
        "atom": {"isotope": "Cs-133", "mass_u": cfg.mass_u},
        "slm": {
            "wavelength_nm": cfg.slm_wavelength_nm, "waist_um": cfg.slm_waist_um,
            "initial_depth_uK": cfg.slm_initial_depth_uK,
            "transport_depth_uK": cfg.slm_transport_depth_uK,
            "center_um": list(cfg.slm_center_um),
            "protocol_b_depth_uK": cfg.slm_protocol_b_depth_uK,
        },
        "aod": {
            "wavelength_nm": cfg.aod_wavelength_nm, "waist_um": cfg.aod_waist_um,
            "transport_depth_uK": cfg.aod_transport_depth_uK,
            "split_distance_um": cfg.split_distance_um,
        },
        "protocol_anchor": {
            "initial_hold_us": cfg.initial_hold_us, "pickup_us": cfg.pickup_us,
            "split_us": cfg.split_us, "long_distance_um": cfg.long_distance_um,
            "outbound_us": cfg.outbound_us,
            "remote_operation_hold_us": cfg.remote_operation_hold_us,
            "inbound_us": cfg.inbound_us, "merge_us": cfg.merge_us,
            "dropoff_us": cfg.dropoff_us, "final_hold_us": cfg.final_hold_us,
            "long_geometry": cfg.long_geometry, "long_profile": cfg.long_profile,
            "split_profile": cfg.split_profile,
        },
        "lensing": {
            "inherit_from_level3": cfg.inherit_lensing_from_level3,
            "nominal_vs_m_per_s": cfg.nominal_vs_m_per_s,
            "axis_signs": list(cfg.lensing_axis_signs),
            "sensitivity_vs_m_per_s": list(cfg.sensitivity_vs_m_per_s),
        },
        "initial_ensemble": {
            "temperature_uK": cfg.temperature_uK,
            "single_round_validation_shots": cfg.single_round_validation_shots,
            "repeated_round_shots": cfg.repeated_round_shots,
            "long_tail_shots": cfg.long_tail_shots,
            "sensitivity_shots": cfg.sensitivity_shots,
            "single_round_seed": cfg.single_round_seed,
            "repeated_round_seed": cfg.repeated_round_seed,
            "long_tail_seed": cfg.long_tail_seed,
            "sensitivity_seed": cfg.sensitivity_seed,
            "development_seed": cfg.development_seed,
        },
        "integration": {
            "validation_dt_us": cfg.validation_dt_us,
            "repeated_dt_us": cfg.repeated_dt_us,
            "convergence_dt_us": cfg.convergence_dt_us,
            "max_rounds": cfg.max_rounds,
            "long_tail_max_rounds": cfg.long_tail_max_rounds,
            "checkpoint_rounds": list(cfg.checkpoint_rounds),
        },
        "coherence": {
            "enabled": cfg.coherence_enabled,
            "eta_slm": cfg.eta_slm, "eta_aod": cfg.eta_aod,
            "eta_sensitivity": list(cfg.eta_sensitivity),
            "dd_modes": list(cfg.dd_modes),
            "dd_pulse_offset_us": cfg.dd_pulse_offset_us,
            "custom_pulse_times_us": list(cfg.custom_pulse_times_us),
        },
        "intensity_noise": {
            "enabled_nominal": cfg.noise_enabled_nominal,
            "model": cfg.noise_model,
            "rms_grid": list(cfg.noise_rms_grid),
            "tau_us_grid": list(cfg.noise_tau_us_grid),
            "cross_correlation": cfg.noise_cross_correlation,
            "noise_seed": cfg.noise_seed,
            "provenance": "assumed_sensitivity_only",
        },
        "classification": {
            "energy_threshold_J": cfg.energy_threshold_J,
            "near_threshold_fraction_of_depth":
                cfg.near_threshold_fraction_of_depth,
            "ambiguous_basin_tolerance": cfg.ambiguous_basin_tolerance,
            "absorbing_protocol_success": cfg.absorbing_protocol_success,
        },
        "upstream": {
            "level2_best_waveform": cfg.upstream_level2_waveform,
            "level2_metrics": cfg.upstream_level2_metrics,
            "level3_metrics": cfg.upstream_level3_metrics,
            "level3_config": cfg.upstream_level3_config,
            "require_sha256_manifest": cfg.require_sha256_manifest,
            "allow_reconstruction_without_reoptimization":
                cfg.allow_reconstruction_without_reoptimization,
        },
        "output": {
            "directory": cfg.output_directory,
            "save_all_checkpoint_states": cfg.save_all_checkpoint_states,
            "save_representative_trajectories":
                cfg.save_representative_trajectories,
            "representative_shots_per_failure_stage":
                cfg.representative_shots_per_failure_stage,
            "bootstrap_samples": cfg.bootstrap_samples,
            "bootstrap_seed": cfg.bootstrap_seed,
        },
    }
