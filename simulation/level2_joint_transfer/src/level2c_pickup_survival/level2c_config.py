"""Level 2C pick-up 存活率复现：level2c 配置段解析（复用 Level 2 配置体系）。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import yaml

from level2_joint_transfer.config import Level2Config, load_config


@dataclass(frozen=True)
class PickupTrajectoryConfig:
    """手工 pick-up 轨迹参数（论文 Fig. 6b 结构）。"""
    durations_us: tuple           # 手工轨迹三档单程时长
    ramp_fraction: float          # 二次升深段占比（论文 48%）
    ml_control_points: int        # ML 轨迹控制点数（论文 14）
    ml_duration_us: float         # ML 单程时长（400 µs）

    @property
    def move_fraction(self):
        return 1.0 - self.ramp_fraction


@dataclass(frozen=True)
class RoundtripConfig:
    """往返序列与 survival(n) 扫描参数。"""
    wait_us: float                # pick-up 后等待段（SLM 关闭）
    max_one_way_transfers: int    # 最大单程转移次数（论文 ML 目标 60）
    scan_shots: int               # 扫描池样本数
    scan_seed: int
    validation_shots: int         # 独立复核池（不复用于扫描决策）
    validation_seed: int
    dt_us: float
    convergence_dt_us: float

    @property
    def wait_s(self):
        return self.wait_us * 1e-6

    @property
    def dt_s(self):
        return self.dt_us * 1e-6

    @property
    def convergence_dt_s(self):
        return self.convergence_dt_us * 1e-6


@dataclass(frozen=True)
class DisorderConfig:
    """位点间差异与实验不完美性系综（逐 shot 抽样；幅度为标定假设）。

    - alignment_sigma_um：AOD–SLM 残余相对对准误差（论文先对 PSF 再按存活率优化）；
    - aod_depth_rel_sigma：RF 功率→AOD 深度标定误差（相对）；
    - slm_depth_rel_sigma：SLM 深度逐位点差异（相对）；
    - aod_waist_rel_sigma：1055/1061 nm 两光路束腰失配（相对，把 w_AOD 从
      假设值升级为带不确定性的参数）；
    - shot_jitter_alignment_sigma_um：逐转移重抽的对准抖动（shot-to-shot）。
    """
    enabled: bool
    alignment_sigma_um: float
    aod_depth_rel_sigma: float
    slm_depth_rel_sigma: float
    aod_waist_rel_sigma: float
    shot_jitter_alignment_sigma_um: float


@dataclass(frozen=True)
class MLOptimizationConfig:
    """ML 轨迹黑盒优化预算（替代论文在线 ML 优化器；见报告说明）。

    优化目标：固定 transfers_for_objective 次连续单程转移的存活率。
    搜索空间：深度/位置各 ml_control_points 个控制点（端点固定），
    内部点 LHS 采样 + 局部精修；只允许训练池，冻结后一次性复核。
    """
    candidates: int
    finalists: int
    refine_steps: int
    perturbation_sigma: float
    shots_per_candidate: int
    candidate_seed: int
    transfers_for_objective: int
    dt_us: float


@dataclass(frozen=True)
class Level2cConfig:
    """Level 2C 聚合配置：base 为 Level 2 结构化配置。"""
    base: Level2Config
    pickup: PickupTrajectoryConfig
    roundtrip: RoundtripConfig
    disorder: DisorderConfig
    ml_optimization: MLOptimizationConfig
    reference_csv: Path
    acceptance_tolerance_band: float
    section: dict


_KEYS = {
    "pickup": {"durations_us", "ramp_fraction", "ml_control_points", "ml_duration_us"},
    "roundtrip": {"wait_us", "max_one_way_transfers", "scan_shots", "scan_seed",
                  "validation_shots", "validation_seed", "dt_us", "convergence_dt_us"},
    "disorder": {"enabled", "alignment_sigma_um", "aod_depth_rel_sigma", "slm_depth_rel_sigma",
                 "aod_waist_rel_sigma", "shot_jitter_alignment_sigma_um"},
    "ml_optimization": {"candidates", "finalists", "refine_steps", "perturbation_sigma",
                        "shots_per_candidate", "candidate_seed", "transfers_for_objective",
                        "dt_us"},
    "acceptance": {"tolerance_band"},
    "reference": {"csv"},
}


def _num(v, label):
    import math
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(float(v)):
        raise ValueError(f"level2c 配置项 {label} 必须是有限数值，实际为 {v!r}")
    return float(v)


def _int(v, label):
    import math
    number = _num(v, label)
    if not math.isclose(number, round(number), abs_tol=1e-12):
        raise ValueError(f"level2c 配置项 {label} 必须是整数，实际为 {number}")
    return int(round(number))


def load_level2c_config(path) -> Level2cConfig:
    """读取含 level2c 段的 YAML；Level 2 段走原有严格校验。"""
    config_path = Path(path).expanduser().resolve()
    base = load_config(config_path, extra_sections=("level2c",))
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))["level2c"]
    if not isinstance(raw, Mapping) or set(raw) != set(_KEYS):
        raise ValueError("level2c 配置段字段与 schema 不一致")
    for section, keys in _KEYS.items():
        if not isinstance(raw[section], Mapping) or set(raw[section]) != keys:
            raise ValueError(f"level2c.{section} 字段与 schema 不一致")

    pk, rt, ds, ml, ac, rf = (raw[s] for s in
                              ("pickup", "roundtrip", "disorder", "ml_optimization",
                               "acceptance", "reference"))
    cfg = Level2cConfig(
        base=base,
        pickup=PickupTrajectoryConfig(
            tuple(_num(v, "pickup.durations_us") for v in pk["durations_us"]),
            _num(pk["ramp_fraction"], "pickup.ramp_fraction"),
            _int(pk["ml_control_points"], "pickup.ml_control_points"),
            _num(pk["ml_duration_us"], "pickup.ml_duration_us")),
        roundtrip=RoundtripConfig(
            _num(rt["wait_us"], "roundtrip.wait_us"),
            _int(rt["max_one_way_transfers"], "roundtrip.max_one_way_transfers"),
            _int(rt["scan_shots"], "roundtrip.scan_shots"), _int(rt["scan_seed"], "roundtrip.scan_seed"),
            _int(rt["validation_shots"], "roundtrip.validation_shots"),
            _int(rt["validation_seed"], "roundtrip.validation_seed"),
            _num(rt["dt_us"], "roundtrip.dt_us"),
            _num(rt["convergence_dt_us"], "roundtrip.convergence_dt_us")),
        disorder=DisorderConfig(
            bool(ds["enabled"]), _num(ds["alignment_sigma_um"], "disorder.alignment_sigma_um"),
            _num(ds["aod_depth_rel_sigma"], "disorder.aod_depth_rel_sigma"),
            _num(ds["slm_depth_rel_sigma"], "disorder.slm_depth_rel_sigma"),
            _num(ds["aod_waist_rel_sigma"], "disorder.aod_waist_rel_sigma"),
            _num(ds["shot_jitter_alignment_sigma_um"], "disorder.shot_jitter_alignment_sigma_um")),
        ml_optimization=MLOptimizationConfig(
            _int(ml["candidates"], "ml_optimization.candidates"), _int(ml["finalists"], "ml_optimization.finalists"),
            _int(ml["refine_steps"], "ml_optimization.refine_steps"),
            _num(ml["perturbation_sigma"], "ml_optimization.perturbation_sigma"),
            _int(ml["shots_per_candidate"], "ml_optimization.shots_per_candidate"),
            _int(ml["candidate_seed"], "ml_optimization.candidate_seed"),
            _int(ml["transfers_for_objective"], "ml_optimization.transfers_for_objective"),
            _num(ml["dt_us"], "ml_optimization.dt_us")),
        reference_csv=(base.project_root / str(rf["csv"])).resolve(),
        acceptance_tolerance_band=_num(ac["tolerance_band"], "acceptance.tolerance_band"),
        section=raw,
    )
    _validate(cfg)
    return cfg


def _validate(c: Level2cConfig) -> None:
    if not (0.0 < c.pickup.ramp_fraction < 1.0):
        raise ValueError("level2c.pickup.ramp_fraction 必须在 (0,1)")
    if any(d <= 0 for d in c.pickup.durations_us):
        raise ValueError("level2c.pickup.durations_us 必须为正")
    if c.pickup.ml_control_points < 4:
        raise ValueError("level2c.pickup.ml_control_points 至少为 4")
    rt = c.roundtrip
    if rt.wait_us <= 0 or rt.max_one_way_transfers < 2:
        raise ValueError("level2c.roundtrip 等待与转移次数非法")
    if rt.scan_seed == rt.validation_seed:
        raise ValueError("level2c 扫描池与复核池种子必须不同")
    if rt.convergence_dt_us >= rt.dt_us:
        raise ValueError("level2c 收敛步长必须小于扫描步长")
    if c.disorder.enabled:
        for label in ("alignment_sigma_um", "aod_depth_rel_sigma", "slm_depth_rel_sigma",
                      "aod_waist_rel_sigma", "shot_jitter_alignment_sigma_um"):
            if getattr(c.disorder, label) < 0:
                raise ValueError(f"level2c.disorder.{label} 必须非负")
    if not (0.0 < c.acceptance_tolerance_band < 1.0):
        raise ValueError("level2c.acceptance.tolerance_band 必须在 (0,1)")
    if not c.reference_csv.is_file():
        raise ValueError(f"找不到论文参考曲线 {c.reference_csv}（先生成 reference 数据）")
    ml = c.ml_optimization
    if ml.candidates <= 0 or ml.finalists <= 0 or ml.shots_per_candidate <= 0:
        raise ValueError("level2c.ml_optimization 预算必须为正")
    if ml.transfers_for_objective < 2 or ml.transfers_for_objective > rt.max_one_way_transfers:
        raise ValueError("ML 优化目标转移次数必须在 [2, max_one_way_transfers]")
