"""Level 3 配置加载与校验（沿用 Level 0-2 的配置风格）。"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass
class Level3Config:
    """解析后的 Level 3 配置。"""
    name: str
    mass_u: float
    wavelength_nm: float
    waist_um: float
    depth_uK: float
    waist_source: str
    temperature_uK: float
    scan_shots: int
    validation_shots: int
    sensitivity_shots: int
    scan_seed: int
    validation_seed: int
    convergence_seed: int
    scan_dt_us: float
    validation_dt_us: float
    convergence_dt_us: float
    post_transport_hold_us: float
    geometries: tuple
    profiles: tuple
    distances_um: tuple
    durations_us: tuple
    calibration_mode: str
    axis_signs: tuple
    reference_geometry: str
    reference_profile: str
    reference_distance_um: float
    reference_duration_us: float
    reference_axis_vmax_over_vs: tuple
    absolute_critical_velocity_m_per_s: float | None
    commanded_depth_mode: str
    enable_exploratory_lensing_compensation: bool
    max_power_multiplier: float
    retained_energy_threshold_J: float
    near_threshold_fraction_of_depth: float
    output_directory: str
    project_root: Path
    config_path: Path


def _root(path: Path) -> Path:
    for candidate in (path.parent, *path.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate.resolve()
    return path.parent.resolve()


def load_config(path: Path) -> Level3Config:
    """读取并校验 Level 3 YAML 配置。"""
    path = Path(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    atom = data["atom"]
    trap = data["aod_trap"]
    ens = data["initial_ensemble"]
    integ = data["integration"]
    tr = data["transport"]
    lens = data["lensing"]
    depth = data["depth_control"]
    cls = data["classification"]
    out = data["output"]
    seeds = {ens["scan_seed"], ens["validation_seed"], ens["convergence_seed"]}
    if len(seeds) != 3:
        raise ValueError("scan/validation/convergence 种子必须互不相同")
    mode = lens["calibration_mode"]
    if mode not in ("lensing_off", "absolute_critical_velocity",
                    "normalized_reference_scan"):
        raise ValueError(f"未知 lensing 校准模式 {mode}")
    if mode == "absolute_critical_velocity" and not lens.get(
            "absolute_critical_velocity_m_per_s"):
        raise ValueError("absolute_critical_velocity 模式必须给出 v_s")
    cfg = Level3Config(
        name=data["level3"]["name"], mass_u=float(atom["mass_u"]),
        wavelength_nm=float(trap["wavelength_nm"]), waist_um=float(trap["waist_um"]),
        depth_uK=float(trap["depth_uK"]), waist_source=str(trap["waist_source"]),
        temperature_uK=float(ens["temperature_uK"]),
        scan_shots=int(ens["scan_shots"]), validation_shots=int(ens["validation_shots"]),
        sensitivity_shots=int(ens["sensitivity_shots"]),
        scan_seed=int(ens["scan_seed"]), validation_seed=int(ens["validation_seed"]),
        convergence_seed=int(ens["convergence_seed"]),
        scan_dt_us=float(integ["scan_dt_us"]),
        validation_dt_us=float(integ["validation_dt_us"]),
        convergence_dt_us=float(integ["convergence_dt_us"]),
        post_transport_hold_us=float(integ["post_transport_hold_us"]),
        geometries=tuple(tr["geometries"]), profiles=tuple(tr["profiles"]),
        distances_um=tuple(float(v) for v in tr["distances_um"]),
        durations_us=tuple(float(v) for v in tr["durations_us"]),
        calibration_mode=mode, axis_signs=tuple(int(v) for v in lens["axis_signs"]),
        reference_geometry=str(lens["reference_geometry"]),
        reference_profile=str(lens["reference_profile"]),
        reference_distance_um=float(lens["reference_distance_um"]),
        reference_duration_us=float(lens["reference_duration_us"]),
        reference_axis_vmax_over_vs=tuple(
            float(v) for v in lens["reference_axis_vmax_over_vs"]),
        absolute_critical_velocity_m_per_s=lens.get(
            "absolute_critical_velocity_m_per_s"),
        commanded_depth_mode=str(depth["commanded_depth_mode"]),
        enable_exploratory_lensing_compensation=bool(
            depth.get("enable_exploratory_lensing_compensation", False)),
        max_power_multiplier=float(depth.get("max_power_multiplier", 2.0)),
        retained_energy_threshold_J=float(cls["retained_energy_threshold_J"]),
        near_threshold_fraction_of_depth=float(
            cls["near_threshold_fraction_of_depth"]),
        output_directory=str(out["directory"]),
        project_root=_root(path), config_path=path,
    )
    return cfg


def config_as_dict(cfg: Level3Config) -> dict:
    """配置转字典（用于 config_used.yaml）。"""
    return {
        "level3": {"name": cfg.name},
        "atom": {"isotope": "Cs-133", "mass_u": cfg.mass_u},
        "aod_trap": {"wavelength_nm": cfg.wavelength_nm, "waist_um": cfg.waist_um,
                     "depth_uK": cfg.depth_uK, "waist_source": cfg.waist_source},
        "initial_ensemble": {
            "temperature_uK": cfg.temperature_uK, "scan_shots": cfg.scan_shots,
            "validation_shots": cfg.validation_shots,
            "sensitivity_shots": cfg.sensitivity_shots, "scan_seed": cfg.scan_seed,
            "validation_seed": cfg.validation_seed,
            "convergence_seed": cfg.convergence_seed},
        "integration": {
            "scan_dt_us": cfg.scan_dt_us, "validation_dt_us": cfg.validation_dt_us,
            "convergence_dt_us": cfg.convergence_dt_us,
            "post_transport_hold_us": cfg.post_transport_hold_us},
        "transport": {
            "geometries": list(cfg.geometries), "profiles": list(cfg.profiles),
            "distances_um": list(cfg.distances_um),
            "durations_us": list(cfg.durations_us)},
        "lensing": {
            "calibration_mode": cfg.calibration_mode,
            "axis_signs": list(cfg.axis_signs),
            "reference_geometry": cfg.reference_geometry,
            "reference_profile": cfg.reference_profile,
            "reference_distance_um": cfg.reference_distance_um,
            "reference_duration_us": cfg.reference_duration_us,
            "reference_axis_vmax_over_vs": list(cfg.reference_axis_vmax_over_vs),
            "absolute_critical_velocity_m_per_s":
                cfg.absolute_critical_velocity_m_per_s},
        "depth_control": {
            "commanded_depth_mode": cfg.commanded_depth_mode,
            "enable_exploratory_lensing_compensation":
                cfg.enable_exploratory_lensing_compensation,
            "max_power_multiplier": cfg.max_power_multiplier},
        "classification": {
            "retained_energy_threshold_J": cfg.retained_energy_threshold_J,
            "near_threshold_fraction_of_depth":
                cfg.near_threshold_fraction_of_depth},
        "output": {"directory": cfg.output_directory},
    }
