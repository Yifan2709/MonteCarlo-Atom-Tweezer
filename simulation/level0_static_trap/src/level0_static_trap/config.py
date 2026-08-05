"""配置加载与校验。

读取 YAML，转换为 SI，校验所有数值条件，并解析输出路径（拒绝写入 ``sources/``）。
失败时抛出 :class:`ConfigError`，信息用中文，并给出实际值与阈值。
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .constants import is_finite_number, um_to_m, u_to_kg, uk_to_j


class ConfigError(ValueError):
    """配置校验失败。"""


# 默认配置文件中允许的顶层键集合；遇到未知键会给出明确警告。
_REQUIRED_TOP_KEYS = {
    "model",
    "atom",
    "slm",
    "aod",
    "grid",
    "trajectory",
    "validation",
    "output",
}


@dataclass
class TrapConfig:
    """单势阱（SLM 或 AOD）的已解析配置（SI 单位，保留来源/假设标签）。"""

    name: str
    depth_uK: float                 # 输入势阱深度 (μK)
    waist_um: float                 # 输入束腰 (μm)
    center_um: float                # 输入中心位置 (μm)
    depth_source: str
    waist_source: str
    # SI 换算结果：
    depth_J: float = 0.0
    waist_m: float = 0.0
    center_m: float = 0.0
    waist_assumed: bool = False     # 束腰是否为假设值

    def to_summary(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "depth_uK": self.depth_uK,
            "waist_um": self.waist_um,
            "center_um": self.center_um,
            "depth_source": self.depth_source,
            "waist_source": self.waist_source,
            "depth_J": self.depth_J,
            "waist_m": self.waist_m,
            "center_m": self.center_m,
            "waist_assumed": self.waist_assumed,
        }


@dataclass
class Config:
    """整个 Level 0 运行的已解析配置。"""

    coordinate: str
    isotope: str
    mass_u: float
    mass_kg: float
    slm: TrapConfig
    aod: TrapConfig
    # 网格
    min_offset_waist: float
    max_offset_waist: float
    num_points: int
    curvature_step_waist: float
    # 轨迹
    initial_displacement_waist: float
    initial_velocity_m_per_s: float
    duration_s: float
    dt_s: float
    frequency_estimator: str
    minimum_complete_cycles: int
    # 校验
    max_omega_dt: float
    curvature_relative_tolerance: float
    harmonic_frequency_relative_tolerance: float
    trajectory_frequency_relative_tolerance: float
    max_energy_error_over_depth: float
    # 输出
    output_dir: str
    save_png: bool = True
    save_csv: bool = True
    save_json: bool = True
    save_md: bool = True
    # 原始输入（用于 config_used.yaml）
    raw: dict[str, Any] = field(default_factory=dict)
    # CLI 覆盖记录
    output_dir_override: str | None = None
    # 运行时未知键警告
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# 校验辅助
# ---------------------------------------------------------------------------

def _check_finite(section: str, key: str, value: Any) -> float:
    if not is_finite_number(value):
        raise ConfigError(f"配置项 [{section}].{key} 不是有限数值，实际值：{value!r}")
    return float(value)


def _check_positive(section: str, key: str, value: Any) -> float:
    v = _check_finite(section, key, value)
    if v <= 0.0:
        raise ConfigError(
            f"配置项 [{section}].{key} 必须为正数，实际值：{v:g}（阈值：> 0）"
        )
    return v


def _check_nonnegative(section: str, key: str, value: Any) -> float:
    v = _check_finite(section, key, value)
    if v < 0.0:
        raise ConfigError(
            f"配置项 [{section}].{key} 不能为负，实际值：{v:g}（阈值：>= 0）"
        )
    return v


def _parse_trap(section_name: str, raw: dict[str, Any], warnings: list[str]) -> TrapConfig:
    expected_keys = {
        "name",
        "depth_uK",
        "waist_um",
        "center_um",
        "depth_source",
        "waist_source",
    }
    extra = set(raw.keys()) - expected_keys
    for k in sorted(extra):
        warnings.append(f"[{section_name}] 含未知字段 {k!r}，已忽略")

    if "name" not in raw:
        raise ConfigError(f"[{section_name}].name 缺失")
    name = str(raw["name"])
    if not name:
        raise ConfigError(f"[{section_name}].name 不能为空")

    depth_uK = _check_positive(section_name, "depth_uK", raw.get("depth_uK", None))
    waist_um = _check_positive(section_name, "waist_um", raw.get("waist_um", None))
    center_um = _check_finite(section_name, "center_um", raw.get("center_um", 0.0))
    depth_source = str(raw.get("depth_source", "unspecified"))
    waist_source = str(raw.get("waist_source", "unspecified"))

    waist_assumed = waist_source.lower().startswith("assumed")

    return TrapConfig(
        name=name,
        depth_uK=depth_uK,
        waist_um=waist_um,
        center_um=center_um,
        depth_source=depth_source,
        waist_source=waist_source,
        depth_J=float(uk_to_j(depth_uK)),
        waist_m=float(um_to_m(waist_um)),
        center_m=float(um_to_m(center_um)),
        waist_assumed=waist_assumed,
    )


# ---------------------------------------------------------------------------
# 输出路径解析
# ---------------------------------------------------------------------------

def _resolve_output_dir(raw_dir: str, override: str | None) -> str:
    """解析输出目录：规范化为绝对路径，拒绝写入 ``sources/`` 内。

    解析规则：相对于 *当前工作目录* 解析相对路径，再做 realpath 规范化。
    若 realpath 落在任意名为 ``sources`` 的目录之内，拒绝。
    """
    chosen = override if override else raw_dir
    if not chosen:
        raise ConfigError("输出目录为空（output.directory 与 --output-dir 均未提供）")
    abs_dir = os.path.abspath(chosen)
    real_dir = os.path.realpath(abs_dir)

    # 防御性检查：拒绝任何解析后位于 sources/ 内的路径。
    parts = Path(real_dir).parts
    if "sources" in parts:
        raise ConfigError(
            f"输出路径解析到只读参考目录 sources/ 内（{real_dir}），已拒绝"
        )
    return real_dir


# ---------------------------------------------------------------------------
# 顶层加载
# ---------------------------------------------------------------------------

def load_config_from_dict(
    data: dict[str, Any], output_dir_override: str | None = None
) -> Config:
    """从已解析的字典构造 :class:`Config`，进行完整校验。"""
    warnings: list[str] = []

    if not isinstance(data, dict):
        raise ConfigError(f"配置顶层必须是映射，实际类型：{type(data).__name__}")

    missing = _REQUIRED_TOP_KEYS - set(data.keys())
    if missing:
        raise ConfigError(f"配置缺少顶层键：{sorted(missing)}")

    unknown_top = set(data.keys()) - _REQUIRED_TOP_KEYS
    for k in sorted(unknown_top):
        warnings.append(f"顶层含未知键 {k!r}，已忽略")

    # ---- model ----
    model = data["model"]
    if not isinstance(model, dict):
        raise ConfigError("[model] 必须是映射")
    coordinate = str(model.get("coordinate", "radial_at_focus"))

    # ---- atom ----
    atom = data["atom"]
    if not isinstance(atom, dict):
        raise ConfigError("[atom] 必须是映射")
    isotope = str(atom.get("isotope", "unspecified"))
    mass_u = _check_positive("atom", "mass_u", atom.get("mass_u", None))
    mass_kg = u_to_kg(mass_u)

    # ---- slm / aod ----
    slm = _parse_trap("slm", data["slm"], warnings)
    aod = _parse_trap("aod", data["aod"], warnings)

    # ---- grid ----
    grid = data["grid"]
    if not isinstance(grid, dict):
        raise ConfigError("[grid] 必须是映射")
    min_off = _check_finite("grid", "min_offset_waist", grid.get("min_offset_waist", None))
    max_off = _check_finite("grid", "max_offset_waist", grid.get("max_offset_waist", None))
    if min_off >= max_off:
        raise ConfigError(
            f"[grid] 偏移范围无序：min={min_off:g} 必须 < max={max_off:g}"
        )
    num_points_raw = grid.get("num_points", None)
    if not isinstance(num_points_raw, int) or isinstance(num_points_raw, bool):
        raise ConfigError(
            f"[grid].num_points 必须是整数，实际值：{num_points_raw!r}"
        )
    if num_points_raw < 5:
        raise ConfigError(
            f"[grid].num_points 至少为 5，实际值：{num_points_raw}（阈值：>= 5）"
        )
    num_points = int(num_points_raw)
    curvature_step = _check_positive(
        "grid", "curvature_step_waist", grid.get("curvature_step_waist", None)
    )

    # ---- trajectory ----
    traj = data["trajectory"]
    if not isinstance(traj, dict):
        raise ConfigError("[trajectory] 必须是映射")
    init_disp = _check_finite(
        "trajectory", "initial_displacement_waist", traj.get("initial_displacement_waist", None)
    )
    init_vel = _check_finite(
        "trajectory", "initial_velocity_m_per_s", traj.get("initial_velocity_m_per_s", 0.0)
    )
    duration_us = _check_positive(
        "trajectory", "duration_us", traj.get("duration_us", None)
    )
    dt_us = _check_positive("trajectory", "dt_us", traj.get("dt_us", None))
    freq_estimator = str(traj.get("frequency_estimator", "interpolated_zero_crossings"))
    min_cycles_raw = traj.get("minimum_complete_cycles", None)
    if not isinstance(min_cycles_raw, int) or isinstance(min_cycles_raw, bool):
        raise ConfigError(
            f"[trajectory].minimum_complete_cycles 必须是整数，实际值：{min_cycles_raw!r}"
        )
    if min_cycles_raw < 1:
        raise ConfigError(
            f"[trajectory].minimum_complete_cycles 至少为 1，实际值：{min_cycles_raw}"
        )
    min_cycles = int(min_cycles_raw)

    duration_s = duration_us * 1.0e-6
    dt_s = dt_us * 1.0e-6

    # 时间数组等间隔且严格递增的隐式检查：dt > 0 已保证；duration 与 dt 的
    # 整数关系不强制（采样起点一致即可），但需保证至少有 2 个采样点。
    if duration_s <= dt_s:
        raise ConfigError(
            f"[trajectory] duration={duration_s:g}s 必须 > dt={dt_s:g}s"
        )

    # ---- validation ----
    vali = data["validation"]
    if not isinstance(vali, dict):
        raise ConfigError("[validation] 必须是映射")
    max_omega_dt = _check_positive("validation", "max_omega_dt", vali.get("max_omega_dt", None))
    curv_tol = _check_positive(
        "validation", "curvature_relative_tolerance", vali.get("curvature_relative_tolerance", None)
    )
    harm_tol = _check_positive(
        "validation",
        "harmonic_frequency_relative_tolerance",
        vali.get("harmonic_frequency_relative_tolerance", None),
    )
    traj_tol = _check_positive(
        "validation",
        "trajectory_frequency_relative_tolerance",
        vali.get("trajectory_frequency_relative_tolerance", None),
    )
    energy_tol = _check_positive(
        "validation",
        "max_energy_error_over_depth",
        vali.get("max_energy_error_over_depth", None),
    )

    # ---- output ----
    out = data["output"]
    if not isinstance(out, dict):
        raise ConfigError("[output] 必须是映射")
    out_dir_raw = str(out.get("directory", ""))
    save_png = bool(out.get("save_png", True))
    save_csv = bool(out.get("save_csv", True))
    save_json = bool(out.get("save_json", True))
    save_md = bool(out.get("save_md", True))
    output_dir = _resolve_output_dir(out_dir_raw, output_dir_override)

    # ---- 物理量级校验：omega*dt < max_omega_dt（用解析小振幅频率）----
    # 用两个势阱中更深（曲率更大、频率更高）的那个作为保守估计。
    from .potentials import trap_frequency_analytic  # 局部导入避免循环依赖

    f_slm = trap_frequency_analytic(slm.depth_J, slm.waist_m, mass_kg)
    f_aod = trap_frequency_analytic(aod.depth_J, aod.waist_m, mass_kg)
    f_max = max(f_slm, f_aod)
    omega_max = 2.0 * math.pi * f_max
    omega_dt = omega_max * dt_s
    if omega_dt >= max_omega_dt:
        raise ConfigError(
            f"[validation] omega*dt = {omega_dt:g} 超过阈值 {max_omega_dt:g}"
            f"（最高解析频率 f0={f_max:g} Hz，dt={dt_s:g} s）。"
            f"请减小 dt 或提高 max_omega_dt。"
        )

    return Config(
        coordinate=coordinate,
        isotope=isotope,
        mass_u=mass_u,
        mass_kg=mass_kg,
        slm=slm,
        aod=aod,
        min_offset_waist=min_off,
        max_offset_waist=max_off,
        num_points=num_points,
        curvature_step_waist=curvature_step,
        initial_displacement_waist=init_disp,
        initial_velocity_m_per_s=init_vel,
        duration_s=duration_s,
        dt_s=dt_s,
        frequency_estimator=freq_estimator,
        minimum_complete_cycles=min_cycles,
        max_omega_dt=max_omega_dt,
        curvature_relative_tolerance=curv_tol,
        harmonic_frequency_relative_tolerance=harm_tol,
        trajectory_frequency_relative_tolerance=traj_tol,
        max_energy_error_over_depth=energy_tol,
        output_dir=output_dir,
        save_png=save_png,
        save_csv=save_csv,
        save_json=save_json,
        save_md=save_md,
        raw=data,
        output_dir_override=output_dir_override,
        warnings=warnings,
    )


def load_config(path: str | os.PathLike, output_dir_override: str | None = None) -> Config:
    """从 YAML 文件加载并校验配置。"""
    path = os.fspath(path)
    if not os.path.isfile(path):
        raise ConfigError(f"配置文件不存在：{path}")
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
    except yaml.YAMLError as exc:  # pragma: no cover - 异常路径
        raise ConfigError(f"YAML 解析失败：{exc}") from exc
    if data is None:
        raise ConfigError("配置文件为空")
    return load_config_from_dict(data, output_dir_override=output_dir_override)


def config_to_used_yaml(cfg: Config) -> dict[str, Any]:
    """生成写入 ``config_used.yaml`` 的可序列化字典（含 CLI 覆盖标记）。"""
    used = _to_plain(cfg.raw)
    # 标注实际使用的输出目录与 CLI 覆盖
    used.setdefault("output", {})["resolved_directory"] = cfg.output_dir
    if cfg.output_dir_override is not None:
        used["output"]["cli_output_dir_override"] = cfg.output_dir_override
    used.setdefault("runtime", {})
    used["runtime"]["warnings"] = list(cfg.warnings)
    used["runtime"]["mass_kg"] = cfg.mass_kg
    used["runtime"]["omega_dt"] = _omega_dt(cfg)
    return used


def _omega_dt(cfg: Config) -> float:
    from .potentials import trap_frequency_analytic

    f_slm = trap_frequency_analytic(cfg.slm.depth_J, cfg.slm.waist_m, cfg.mass_kg)
    f_aod = trap_frequency_analytic(cfg.aod.depth_J, cfg.aod.waist_m, cfg.mass_kg)
    return 2.0 * math.pi * max(f_slm, f_aod) * cfg.dt_s


def _to_plain(obj: Any) -> Any:
    """把含 numpy/Path 等的对象转为纯 Python 可序列化结构。"""
    import numpy as np

    if isinstance(obj, dict):
        return {str(k): _to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_plain(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return v if math.isfinite(v) else None
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    return obj
