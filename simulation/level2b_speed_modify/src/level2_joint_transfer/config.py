"""Level 2 严格 YAML 配置加载与 SI 单位视图。"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import math
from types import SimpleNamespace
from typing import Any, Mapping

import yaml

from level0_static_trap.constants import atomic_mass_to_kg, micrometer_to_meter, microkelvin_to_joule


@dataclass(frozen=True)
class Level2Info:
    name: str


@dataclass(frozen=True)
class AtomConfig:
    isotope: str
    mass_u: float

    @property
    def mass_kg(self):
        return atomic_mass_to_kg(self.mass_u)


@dataclass(frozen=True)
class SlmTrap:
    depth_uK: float
    waist_um: float
    center_um: float
    waist_source: str

    @property
    def depth_j(self):
        return microkelvin_to_joule(self.depth_uK)

    @property
    def waist_m(self):
        return micrometer_to_meter(self.waist_um)

    @property
    def center_m(self):
        return micrometer_to_meter(self.center_um)


@dataclass(frozen=True)
class AodTrap:
    initial_depth_uK: float
    waist_um: float
    initial_center_um: float
    waist_source: str

    @property
    def initial_depth_j(self):
        return microkelvin_to_joule(self.initial_depth_uK)

    @property
    def waist_m(self):
        return micrometer_to_meter(self.waist_um)

    @property
    def initial_center_m(self):
        return micrometer_to_meter(self.initial_center_um)


@dataclass(frozen=True)
class InitialEnsemble:
    temperature_uK: float
    sampler: str
    optimization_shots: int
    selection_shots: int
    validation_shots: int
    optimization_seed: int
    selection_seed: int
    validation_seed: int
    level1_reproduction_seed: int


@dataclass(frozen=True)
class IntegrationConfig:
    optimization_dt_us: float
    validation_dt_us: float
    convergence_dt_us: float
    post_transfer_hold_us: float

    @property
    def optimization_dt_s(self):
        return self.optimization_dt_us * 1e-6

    @property
    def validation_dt_s(self):
        return self.validation_dt_us * 1e-6

    @property
    def convergence_dt_s(self):
        return self.convergence_dt_us * 1e-6

    @property
    def post_transfer_hold_s(self):
        return self.post_transfer_hold_us * 1e-6


@dataclass(frozen=True)
class WaveformConfig:
    sequential_baseline_duration_us: float
    compressed_baseline_duration_us: float
    optimized_duration_us: float
    level1_move_fraction: float
    level1_ramp_fraction: float
    min_segment_fraction: float
    ramp_power_bounds: tuple

    @property
    def sequential_baseline_duration_s(self):
        return self.sequential_baseline_duration_us * 1e-6

    @property
    def compressed_baseline_duration_s(self):
        return self.compressed_baseline_duration_us * 1e-6

    @property
    def optimized_duration_s(self):
        return self.optimized_duration_us * 1e-6

    def as_constraint_dict(self):
        """供 validate_waveform 使用的约束字典。"""
        return {
            "min_segment_fraction": self.min_segment_fraction,
            "ramp_power_bounds": list(self.ramp_power_bounds),
        }


@dataclass(frozen=True)
class OptimizationConfig:
    method: str
    lhs_seed: int
    max_candidates: int
    finalists: int
    refine_steps_per_finalist: int
    refine_perturbation_sigma: float
    enable_monotone_spline_refinement: bool
    spline_control_points: int
    spline_candidates: int
    capture_weight: float
    excitation_weight: float
    smoothness_weight: float
    margin_quantile: float


@dataclass(frozen=True)
class ClassificationConfig:
    capture_energy_threshold_J: float
    near_threshold_fraction_of_slm_depth: float


@dataclass(frozen=True)
class RobustnessConfig:
    temperatures_uK: tuple
    final_alignment_offsets_um: tuple
    shots_per_condition: int
    timestep_subset_shots: int
    bootstrap_samples: int
    bootstrap_seed: int


@dataclass(frozen=True)
class NoiseMeanHeatingConfig:
    """三类噪声的平均强度确定性加热（常数功率、无随机性）。

    显式速率优先；为 null 时按物理常数计算：反冲取标称 AOD 深度的
    Γ_sc（Cs D2 两能级近似），参量/指向以 SLM 解析阱频为参考频率。
    参量通道为线性化常数功率 P = Γ_par·k_B·T_ref，
    T_ref 缺省取 initial_ensemble.temperature_uK。
    所有幅度参数均为 assumed_sensitivity_only。
    """
    enabled: bool
    trap_wavelength_nm: float
    recoil_scattering_rate_s: float | None
    parametric_rate_s: float | None
    parametric_reference_temperature_uK: float | None
    parametric_rin_psd_per_hz: float
    pointing_heating_rate_uK_per_s: float | None
    pointing_position_psd_m2_per_hz: float
    reference_nu0_hz: float | None


@dataclass(frozen=True)
class OutputConfig:
    directory: str
    resolved_directory: Path
    save_candidate_table: bool
    save_representative_trajectories: bool
    representative_shots_per_outcome: int


@dataclass(frozen=True)
class Level2Config:
    level2: Level2Info
    atom: AtomConfig
    slm: SlmTrap
    aod: AodTrap
    initial_ensemble: InitialEnsemble
    integration: IntegrationConfig
    waveforms: WaveformConfig
    optimization: OptimizationConfig
    classification: ClassificationConfig
    robustness: RobustnessConfig
    noise_mean_heating: NoiseMeanHeatingConfig
    output: OutputConfig
    config_path: Path
    project_root: Path


_KEYS = {
    "level2": {"name"},
    "atom": {"isotope", "mass_u"},
    "traps": {"slm", "aod"},
    "initial_ensemble": {
        "temperature_uK", "sampler", "optimization_shots", "selection_shots", "validation_shots",
        "optimization_seed", "selection_seed", "validation_seed", "level1_reproduction_seed",
    },
    "integration": {"optimization_dt_us", "validation_dt_us", "convergence_dt_us", "post_transfer_hold_us"},
    "waveforms": {
        "sequential_baseline_duration_us", "compressed_baseline_duration_us", "optimized_duration_us",
        "level1_move_fraction", "level1_ramp_fraction", "min_segment_fraction", "ramp_power_bounds",
    },
    "optimization": {
        "method", "lhs_seed", "max_candidates", "finalists", "refine_steps_per_finalist",
        "refine_perturbation_sigma", "enable_monotone_spline_refinement", "spline_control_points",
        "spline_candidates", "capture_weight", "excitation_weight", "smoothness_weight", "margin_quantile",
    },
    "classification": {"capture_energy_threshold_J", "near_threshold_fraction_of_slm_depth"},
    "robustness": {
        "temperatures_uK", "final_alignment_offsets_um", "shots_per_condition",
        "timestep_subset_shots", "bootstrap_samples", "bootstrap_seed",
    },
    "noise_mean_heating": {
        "enabled", "trap_wavelength_nm", "recoil_scattering_rate_s", "parametric_rate_s",
        "parametric_reference_temperature_uK", "parametric_rin_psd_per_hz",
        "pointing_heating_rate_uK_per_s",
        "pointing_position_psd_m2_per_hz", "reference_nu0_hz",
    },
    "output": {
        "directory", "save_candidate_table", "save_representative_trajectories",
        "representative_shots_per_outcome",
    },
}


def _number(value, label) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise ValueError(f"配置项 {label} 必须是有限数值，实际为 {value!r}")
    return float(value)


def _integer(value, label) -> int:
    number = _number(value, label)
    if not math.isclose(number, round(number), abs_tol=1e-12):
        raise ValueError(f"配置项 {label} 必须是整数，实际为 {number}")
    return int(round(number))


def _root(path: Path) -> Path:
    for candidate in (path.parent, *path.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate.resolve()
    return path.parent.resolve()


def load_config(path) -> Level2Config:
    """读取并校验 Level 2 YAML，返回结构化配置。"""
    config_path = Path(path).expanduser().resolve()
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(raw, Mapping) or set(raw) != set(_KEYS):
        raise ValueError("配置顶层段必须与 schema 完全一致")
    for section, keys in _KEYS.items():
        if not isinstance(raw[section], Mapping) or set(raw[section]) != keys:
            raise ValueError(f"配置段 {section!r} 字段与 schema 不一致")

    n = lambda s, k: _number(raw[s][k], f"{s}.{k}")
    i = lambda s, k: _integer(raw[s][k], f"{s}.{k}")
    traps = raw["traps"]
    ens = raw["initial_ensemble"]
    opt = raw["optimization"]
    rob = raw["robustness"]
    nse = raw["noise_mean_heating"]

    def opt_number(section, key, label):
        value = section[key]
        if value is None:
            return None
        return _number(value, label)

    out = raw["output"]
    root = _root(config_path)
    config = Level2Config(
        level2=Level2Info(str(raw["level2"]["name"])),
        atom=AtomConfig(str(raw["atom"]["isotope"]), n("atom", "mass_u")),
        slm=SlmTrap(_number(traps["slm"]["depth_uK"], "traps.slm.depth_uK"),
                    _number(traps["slm"]["waist_um"], "traps.slm.waist_um"),
                    _number(traps["slm"]["center_um"], "traps.slm.center_um"),
                    str(traps["slm"]["waist_source"])),
        aod=AodTrap(_number(traps["aod"]["initial_depth_uK"], "traps.aod.initial_depth_uK"),
                    _number(traps["aod"]["waist_um"], "traps.aod.waist_um"),
                    _number(traps["aod"]["initial_center_um"], "traps.aod.initial_center_um"),
                    str(traps["aod"]["waist_source"])),
        initial_ensemble=InitialEnsemble(
            n("initial_ensemble", "temperature_uK"), str(ens["sampler"]),
            i("initial_ensemble", "optimization_shots"), i("initial_ensemble", "selection_shots"),
            i("initial_ensemble", "validation_shots"), i("initial_ensemble", "optimization_seed"),
            i("initial_ensemble", "selection_seed"), i("initial_ensemble", "validation_seed"),
            i("initial_ensemble", "level1_reproduction_seed")),
        integration=IntegrationConfig(
            n("integration", "optimization_dt_us"), n("integration", "validation_dt_us"),
            n("integration", "convergence_dt_us"), n("integration", "post_transfer_hold_us")),
        waveforms=WaveformConfig(
            n("waveforms", "sequential_baseline_duration_us"), n("waveforms", "compressed_baseline_duration_us"),
            n("waveforms", "optimized_duration_us"), n("waveforms", "level1_move_fraction"),
            n("waveforms", "level1_ramp_fraction"), n("waveforms", "min_segment_fraction"),
            tuple(_number(v, "waveforms.ramp_power_bounds") for v in raw["waveforms"]["ramp_power_bounds"])),
        optimization=OptimizationConfig(
            str(opt["method"]), i("optimization", "lhs_seed"), i("optimization", "max_candidates"),
            i("optimization", "finalists"), i("optimization", "refine_steps_per_finalist"),
            n("optimization", "refine_perturbation_sigma"),
            bool(opt["enable_monotone_spline_refinement"]), i("optimization", "spline_control_points"),
            i("optimization", "spline_candidates"), n("optimization", "capture_weight"),
            n("optimization", "excitation_weight"), n("optimization", "smoothness_weight"),
            n("optimization", "margin_quantile")),
        classification=ClassificationConfig(
            n("classification", "capture_energy_threshold_J"),
            n("classification", "near_threshold_fraction_of_slm_depth")),
        robustness=RobustnessConfig(
            tuple(_number(v, "robustness.temperatures_uK") for v in rob["temperatures_uK"]),
            tuple(_number(v, "robustness.final_alignment_offsets_um") for v in rob["final_alignment_offsets_um"]),
            i("robustness", "shots_per_condition"), i("robustness", "timestep_subset_shots"),
            i("robustness", "bootstrap_samples"), i("robustness", "bootstrap_seed")),
        noise_mean_heating=NoiseMeanHeatingConfig(
            bool(nse["enabled"]), n("noise_mean_heating", "trap_wavelength_nm"),
            opt_number(nse, "recoil_scattering_rate_s", "noise_mean_heating.recoil_scattering_rate_s"),
            opt_number(nse, "parametric_rate_s", "noise_mean_heating.parametric_rate_s"),
            opt_number(nse, "parametric_reference_temperature_uK",
                       "noise_mean_heating.parametric_reference_temperature_uK"),
            n("noise_mean_heating", "parametric_rin_psd_per_hz"),
            opt_number(nse, "pointing_heating_rate_uK_per_s",
                       "noise_mean_heating.pointing_heating_rate_uK_per_s"),
            n("noise_mean_heating", "pointing_position_psd_m2_per_hz"),
            opt_number(nse, "reference_nu0_hz", "noise_mean_heating.reference_nu0_hz")),
        output=OutputConfig(
            str(out["directory"]), (root / str(out["directory"])).resolve(),
            bool(out["save_candidate_table"]), bool(out["save_representative_trajectories"]),
            i("output", "representative_shots_per_outcome")),
        config_path=config_path, project_root=root,
    )
    _validate(config)
    return config


def _validate(c: Level2Config) -> None:
    if c.atom.mass_u <= 0 or c.slm.depth_uK <= 0 or c.aod.initial_depth_uK <= 0:
        raise ValueError("质量与阱深必须为正")
    if c.slm.waist_um <= 0 or c.aod.waist_um <= 0:
        raise ValueError("束腰必须为正")
    if c.aod.initial_center_um <= 0:
        raise ValueError("AOD 初始中心必须为正（以 SLM 为原点）")
    if c.initial_ensemble.temperature_uK <= 0:
        raise ValueError("初态温度必须为正")
    if c.initial_ensemble.sampler != "level1_harmonic_rejection":
        raise ValueError("初态采样器必须复用 level1_harmonic_rejection")
    seeds = {c.initial_ensemble.optimization_seed, c.initial_ensemble.selection_seed,
             c.initial_ensemble.validation_seed}
    if len(seeds) != 3:
        raise ValueError("optimization/selection/validation 种子必须互不相同")
    if c.initial_ensemble.optimization_shots <= 0:
        raise ValueError("样本数必须为正")
    if c.initial_ensemble.selection_shots <= c.initial_ensemble.optimization_shots:
        raise ValueError("selection 池应大于 optimization 池")
    if c.initial_ensemble.validation_shots <= c.initial_ensemble.selection_shots:
        raise ValueError("validation 池应大于 selection 池")
    for label in ("optimization_dt_us", "validation_dt_us", "convergence_dt_us", "post_transfer_hold_us"):
        if getattr(c.integration, label) <= 0:
            raise ValueError(f"integration.{label} 必须为正")
    if c.integration.convergence_dt_us >= c.integration.validation_dt_us:
        raise ValueError("convergence 步长必须小于 validation 步长")
    w = c.waveforms
    if w.compressed_baseline_duration_us >= w.sequential_baseline_duration_us:
        raise ValueError("压缩基线时长必须短于顺序基线")
    if w.optimized_duration_us <= 0 or w.sequential_baseline_duration_us <= 0:
        raise ValueError("波形时长必须为正")
    if not math.isclose(w.level1_move_fraction + w.level1_ramp_fraction, 1.0, abs_tol=1e-12):
        raise ValueError("Level 1 移动/降深配比之和必须为 1")
    if not (0.0 < w.min_segment_fraction < 0.5):
        raise ValueError("min_segment_fraction 必须在 (0, 0.5)")
    if len(w.ramp_power_bounds) != 2 or w.ramp_power_bounds[0] >= w.ramp_power_bounds[1]:
        raise ValueError("ramp_power_bounds 必须是递增二元组")
    if c.optimization.margin_quantile <= 0 or c.optimization.margin_quantile >= 1:
        raise ValueError("margin_quantile 必须在 (0,1)")
    if not (0.5 <= w.ramp_power_bounds[0] and w.ramp_power_bounds[1] <= 4.0):
        raise ValueError("ramp_power_bounds 应位于 [0.5, 4.0]")
    if c.robustness.shots_per_condition <= 0:
        raise ValueError("鲁棒性每条件样本数必须为正")
    noise = c.noise_mean_heating
    if noise.trap_wavelength_nm <= 0:
        raise ValueError("noise_mean_heating.trap_wavelength_nm 必须为正")
    for label in ("recoil_scattering_rate_s", "parametric_rate_s",
                  "parametric_reference_temperature_uK",
                  "parametric_rin_psd_per_hz", "pointing_heating_rate_uK_per_s",
                  "pointing_position_psd_m2_per_hz", "reference_nu0_hz"):
        value = getattr(noise, label)
        if value is not None and value < 0:
            raise ValueError(f"noise_mean_heating.{label} 必须非负")
    # 步长合法性：最大阱频（双阱叠加的保守估计）乘最大步长需小于 0.05
    omega_max = math.sqrt(4 * (c.slm.depth_j + c.aod.initial_depth_j) /
                          (c.atom.mass_kg * min(c.slm.waist_m, c.aod.waist_m) ** 2))
    for label in ("optimization_dt_s", "validation_dt_s", "convergence_dt_s"):
        dt = getattr(c.integration, label)
        if omega_max * dt >= 0.05:
            raise ValueError(f"{label}={dt*1e6} us 未满足 omega_max*dt<0.05（{omega_max*dt:.4g}）")


def with_output_override(c: Level2Config, directory) -> Level2Config:
    """应用 CLI 输出目录覆盖；相对路径按项目根解析。"""
    text = str(directory)
    path = Path(text).expanduser()
    resolved = (path if path.is_absolute() else c.project_root / path).resolve()
    return replace(c, output=replace(c.output, directory=text, resolved_directory=resolved))


def config_as_dict(c: Level2Config) -> dict:
    """导出含实际覆盖值与全部种子的 YAML 字典。"""

    def clean(value: Any) -> Any:
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, tuple):
            return list(value)
        return value

    data = {}
    for field in ("level2", "atom", "initial_ensemble", "integration",
                  "waveforms", "optimization", "classification", "robustness",
                  "noise_mean_heating"):
        data[field] = {k: clean(v) for k, v in vars(getattr(c, field)).items()}
    data["traps"] = {
        "slm": {k: clean(v) for k, v in vars(c.slm).items()},
        "aod": {k: clean(v) for k, v in vars(c.aod).items()},
    }
    data["output"] = {
        "directory": c.output.directory,
        "resolved_directory": str(c.output.resolved_directory),
        "save_candidate_table": c.output.save_candidate_table,
        "save_representative_trajectories": c.output.save_representative_trajectories,
        "representative_shots_per_outcome": c.output.representative_shots_per_outcome,
    }
    data["config_path"] = str(c.config_path)
    data["project_root"] = str(c.project_root)
    return data


def level1_view(c: Level2Config, temperature_uK=None, aod_initial_center_m=None):
    """构造与 Level 1 采样器接口兼容的视图对象，保证逐调用复用同一实现。"""
    temperature = c.initial_ensemble.temperature_uK if temperature_uK is None else temperature_uK
    center = c.aod.initial_center_m if aod_initial_center_m is None else aod_initial_center_m
    return SimpleNamespace(
        atom=SimpleNamespace(mass_kg=c.atom.mass_kg),
        aod=SimpleNamespace(depth_j=c.aod.initial_depth_j, waist_m=c.aod.waist_m, initial_center_m=center),
        thermal=SimpleNamespace(temperature_uK=temperature),
    )
