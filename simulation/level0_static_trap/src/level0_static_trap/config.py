"""YAML 配置加载、严格校验和 SI 单位视图。"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import math
from typing import Any, Mapping

import yaml

from .constants import atomic_mass_to_kg, micrometer_to_meter, microkelvin_to_joule


@dataclass(frozen=True)
class ModelConfig:
    coordinate: str


@dataclass(frozen=True)
class AtomConfig:
    isotope: str
    mass_u: float

    @property
    def mass_kg(self) -> float:
        return atomic_mass_to_kg(self.mass_u)


@dataclass(frozen=True)
class TrapConfig:
    name: str
    depth_uK: float
    waist_um: float
    center_um: float
    depth_source: str
    waist_source: str

    @property
    def depth_j(self) -> float:
        return microkelvin_to_joule(self.depth_uK)

    @property
    def waist_m(self) -> float:
        return micrometer_to_meter(self.waist_um)

    @property
    def center_m(self) -> float:
        return micrometer_to_meter(self.center_um)

    @property
    def waist_is_assumed(self) -> bool:
        return "assumed" in self.waist_source.lower()


@dataclass(frozen=True)
class GridConfig:
    min_offset_waist: float
    max_offset_waist: float
    num_points: int
    curvature_step_waist: float


@dataclass(frozen=True)
class TrajectoryConfig:
    initial_displacement_waist: float
    initial_velocity_m_per_s: float
    duration_us: float
    dt_us: float
    frequency_estimator: str
    minimum_complete_cycles: int

    @property
    def duration_s(self) -> float:
        return self.duration_us * 1.0e-6

    @property
    def dt_s(self) -> float:
        return self.dt_us * 1.0e-6


@dataclass(frozen=True)
class ValidationConfig:
    max_omega_dt: float
    curvature_relative_tolerance: float
    harmonic_frequency_relative_tolerance: float
    trajectory_frequency_relative_tolerance: float
    max_energy_error_over_depth: float


@dataclass(frozen=True)
class OutputConfig:
    directory: str
    resolved_directory: Path
    save_png: bool
    save_csv: bool
    save_json: bool


@dataclass(frozen=True)
class SimulationConfig:
    model: ModelConfig
    atom: AtomConfig
    slm: TrapConfig
    aod: TrapConfig
    grid: GridConfig
    trajectory: TrajectoryConfig
    validation: ValidationConfig
    output: OutputConfig
    config_path: Path
    project_root: Path


_SECTION_KEYS = {
    "model": {"coordinate"},
    "atom": {"isotope", "mass_u"},
    "slm": {"name", "depth_uK", "waist_um", "center_um", "depth_source", "waist_source"},
    "aod": {"name", "depth_uK", "waist_um", "center_um", "depth_source", "waist_source"},
    "grid": {"min_offset_waist", "max_offset_waist", "num_points", "curvature_step_waist"},
    "trajectory": {
        "initial_displacement_waist", "initial_velocity_m_per_s", "duration_us", "dt_us",
        "frequency_estimator", "minimum_complete_cycles",
    },
    "validation": {
        "max_omega_dt", "curvature_relative_tolerance", "harmonic_frequency_relative_tolerance",
        "trajectory_frequency_relative_tolerance", "max_energy_error_over_depth",
    },
    "output": {"directory", "save_png", "save_csv", "save_json"},
}


def _finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise ValueError(f"配置项 {label} 必须是有限数值，实际为 {value!r}")
    return float(value)


def _strict_section(raw: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    if name not in raw or not isinstance(raw[name], Mapping):
        raise ValueError(f"缺少配置段 {name!r} 或其不是映射")
    keys = set(raw[name])
    expected = _SECTION_KEYS[name]
    unknown = keys - expected
    missing = expected - keys
    if unknown:
        raise ValueError(f"配置段 {name!r} 含未知字段: {sorted(unknown)}")
    if missing:
        raise ValueError(f"配置段 {name!r} 缺少字段: {sorted(missing)}")
    return raw[name]


def _project_root(config_path: Path) -> Path:
    for candidate in (config_path.parent, *config_path.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate.resolve()
    return config_path.parent.resolve()


def _resolve_output(directory: str, project_root: Path) -> Path:
    path = Path(directory).expanduser()
    return (path if path.is_absolute() else project_root / path).resolve()


def _validate_config(config: SimulationConfig) -> None:
    if config.model.coordinate != "radial_at_focus":
        raise ValueError("model.coordinate 必须为 radial_at_focus（焦平面径向坐标）")
    if not config.atom.isotope.strip():
        raise ValueError("atom.isotope 不得为空")
    if config.atom.mass_u <= 0.0:
        raise ValueError(f"原子质量必须大于零，实际为 {config.atom.mass_u}")
    for trap in (config.slm, config.aod):
        if trap.depth_uK <= 0.0 or trap.waist_um <= 0.0:
            raise ValueError(
                f"{trap.name} 的深度和束腰必须大于零，实际为 depth={trap.depth_uK} μK, waist={trap.waist_um} μm"
            )
    if config.grid.num_points < 5:
        raise ValueError(f"grid.num_points 至少为 5，实际为 {config.grid.num_points}")
    if config.grid.min_offset_waist >= config.grid.max_offset_waist:
        raise ValueError("网格无量纲偏移下限必须小于上限")
    if config.grid.curvature_step_waist <= 0.0:
        raise ValueError("曲率差分步长必须为正")
    trajectory = config.trajectory
    if trajectory.duration_us <= 0.0 or trajectory.dt_us <= 0.0:
        raise ValueError(
            f"轨迹时长和步长必须大于零，实际为 duration={trajectory.duration_us} μs, dt={trajectory.dt_us} μs"
        )
    if trajectory.frequency_estimator != "interpolated_zero_crossings":
        raise ValueError("frequency_estimator 目前仅支持 interpolated_zero_crossings")
    if trajectory.minimum_complete_cycles < 1:
        raise ValueError("minimum_complete_cycles 至少为 1")
    ratio = trajectory.duration_s / trajectory.dt_s
    if not math.isclose(ratio, round(ratio), rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError(f"duration_us/dt_us 必须为整数以构造等间隔时间数组，实际为 {ratio:.12g}")
    for label, value in vars(config.validation).items():
        if value <= 0.0:
            raise ValueError(f"validation.{label} 必须大于零，实际为 {value}")

    mass = config.atom.mass_kg
    dt = trajectory.dt_s
    for trap in (config.slm, config.aod):
        omega = math.sqrt(4.0 * trap.depth_j / (mass * trap.waist_m**2))
        omega_dt = omega * dt
        if omega_dt >= config.validation.max_omega_dt:
            raise ValueError(
                f"{trap.name} 时间步长过大：ωΔt={omega_dt:.6g}，必须小于阈值 {config.validation.max_omega_dt:.6g}"
            )
        cycles = trajectory.duration_s * omega / (2.0 * math.pi)
        if cycles < trajectory.minimum_complete_cycles:
            raise ValueError(
                f"{trap.name} 预计完整周期数不足：{cycles:.3f}，要求至少 {trajectory.minimum_complete_cycles}"
            )
        q0 = trajectory.initial_displacement_waist
        initial_potential = -trap.depth_j * math.exp(-2.0 * q0 * q0)
        initial_energy = 0.5 * mass * trajectory.initial_velocity_m_per_s**2 + initial_potential
        if initial_energy >= 0.0:
            raise ValueError(
                f"{trap.name} 初始状态不束缚：E0={initial_energy:.6e} J，必须满足 E0<0；停止轨迹频率分析"
            )


def load_config(path: str | Path) -> SimulationConfig:
    """读取 YAML，拒绝未知字段并返回经过完整校验的结构化配置。"""

    config_path = Path(path).expanduser().resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, Mapping):
        raise ValueError("YAML 顶层必须是映射")
    unknown_sections = set(raw) - set(_SECTION_KEYS)
    missing_sections = set(_SECTION_KEYS) - set(raw)
    if unknown_sections:
        raise ValueError(f"配置含未知顶层段: {sorted(unknown_sections)}")
    if missing_sections:
        raise ValueError(f"配置缺少顶层段: {sorted(missing_sections)}")
    sections = {name: _strict_section(raw, name) for name in _SECTION_KEYS}

    def numeric(section: str, key: str) -> float:
        return _finite_number(sections[section][key], f"{section}.{key}")

    root = _project_root(config_path)
    output_directory = str(sections["output"]["directory"])
    config = SimulationConfig(
        model=ModelConfig(coordinate=str(sections["model"]["coordinate"])),
        atom=AtomConfig(isotope=str(sections["atom"]["isotope"]), mass_u=numeric("atom", "mass_u")),
        slm=TrapConfig(
            name=str(sections["slm"]["name"]), depth_uK=numeric("slm", "depth_uK"),
            waist_um=numeric("slm", "waist_um"), center_um=numeric("slm", "center_um"),
            depth_source=str(sections["slm"]["depth_source"]),
            waist_source=str(sections["slm"]["waist_source"]),
        ),
        aod=TrapConfig(
            name=str(sections["aod"]["name"]), depth_uK=numeric("aod", "depth_uK"),
            waist_um=numeric("aod", "waist_um"), center_um=numeric("aod", "center_um"),
            depth_source=str(sections["aod"]["depth_source"]),
            waist_source=str(sections["aod"]["waist_source"]),
        ),
        grid=GridConfig(
            min_offset_waist=numeric("grid", "min_offset_waist"),
            max_offset_waist=numeric("grid", "max_offset_waist"),
            num_points=int(numeric("grid", "num_points")),
            curvature_step_waist=numeric("grid", "curvature_step_waist"),
        ),
        trajectory=TrajectoryConfig(
            initial_displacement_waist=numeric("trajectory", "initial_displacement_waist"),
            initial_velocity_m_per_s=numeric("trajectory", "initial_velocity_m_per_s"),
            duration_us=numeric("trajectory", "duration_us"), dt_us=numeric("trajectory", "dt_us"),
            frequency_estimator=str(sections["trajectory"]["frequency_estimator"]),
            minimum_complete_cycles=int(numeric("trajectory", "minimum_complete_cycles")),
        ),
        validation=ValidationConfig(**{
            key: numeric("validation", key) for key in _SECTION_KEYS["validation"]
        }),
        output=OutputConfig(
            directory=output_directory,
            resolved_directory=_resolve_output(output_directory, root),
            save_png=bool(sections["output"]["save_png"]),
            save_csv=bool(sections["output"]["save_csv"]),
            save_json=bool(sections["output"]["save_json"]),
        ),
        config_path=config_path,
        project_root=root,
    )
    _validate_config(config)
    return config


def with_output_override(config: SimulationConfig, directory: str | Path) -> SimulationConfig:
    """应用 CLI 输出目录覆盖；相对路径始终以项目根目录解析。"""

    text = str(directory)
    output = replace(config.output, directory=text, resolved_directory=_resolve_output(text, config.project_root))
    return replace(config, output=output)


def config_as_dict(config: SimulationConfig) -> dict[str, Any]:
    """生成保留输入单位、来源、假设标签和实际输出覆盖值的 YAML 字典。"""

    return {
        "model": {"coordinate": config.model.coordinate},
        "atom": {"isotope": config.atom.isotope, "mass_u": config.atom.mass_u},
        "slm": {
            "name": config.slm.name, "depth_uK": config.slm.depth_uK, "waist_um": config.slm.waist_um,
            "center_um": config.slm.center_um, "depth_source": config.slm.depth_source,
            "waist_source": config.slm.waist_source, "waist_assumption": config.slm.waist_is_assumed,
        },
        "aod": {
            "name": config.aod.name, "depth_uK": config.aod.depth_uK, "waist_um": config.aod.waist_um,
            "center_um": config.aod.center_um, "depth_source": config.aod.depth_source,
            "waist_source": config.aod.waist_source, "waist_assumption": config.aod.waist_is_assumed,
        },
        "grid": vars(config.grid),
        "trajectory": vars(config.trajectory),
        "validation": vars(config.validation),
        "output": {
            "directory": config.output.directory, "resolved_directory": str(config.output.resolved_directory),
            "save_png": config.output.save_png, "save_csv": config.output.save_csv,
            "save_json": config.output.save_json,
        },
    }

