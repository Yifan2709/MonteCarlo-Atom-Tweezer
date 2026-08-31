"""Level 2B 配置加载：本级 YAML + Level 2 基线配置（物理与采样保持一致）。"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import yaml

from level2_joint_transfer.config import load_config as load_level2_config

REQUIRED_SECTIONS = {
    "level2b", "initial_ensemble", "speed_scan", "follow_loss", "integration",
    "comparison", "noise_speed_validation", "max_speed_probe", "output",
}


def load_level2b_config(path: str | Path) -> dict:
    """读取并校验 Level 2B YAML，返回 (本级配置 dict, Level 2 配置对象)。"""
    config_path = Path(path).expanduser().resolve()
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(raw, Mapping) or not REQUIRED_SECTIONS.issubset(raw):
        missing = REQUIRED_SECTIONS - set(raw if isinstance(raw, Mapping) else {})
        raise ValueError(f"Level 2B 配置缺少段落：{sorted(missing)}")
    bcfg = _to_plain(raw)
    level2_name = Path(str(bcfg["level2b"].get(
        "level2_baseline_config", "configs/level2_joint_transfer.yaml"))).name
    candidates = [
        Path(str(bcfg["level2b"]["level2_baseline_config"])),
        config_path.parent / level2_name,
        Path(__file__).resolve().parents[2] / "configs" / level2_name,
    ]
    level2_yaml = next((p for p in candidates if p.is_file()), None)
    if level2_yaml is None:
        raise ValueError(f"找不到 Level 2 基线配置，尝试过：{[str(p) for p in candidates]}")
    cfg2 = load_level2_config(level2_yaml)
    if int(bcfg["initial_ensemble"]["scan_seed"]) == int(bcfg["initial_ensemble"]["boundary_seed"]):
        raise ValueError("scan_seed 与 boundary_seed 必须互不相同")
    return bcfg, cfg2


def _to_plain(value: Any) -> Any:
    """递归转换为可 JSON 序列化的普通容器。"""
    if isinstance(value, Mapping):
        return {str(k): _to_plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_plain(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def config_with_output(bcfg: dict, directory: str | Path) -> dict:
    """返回输出目录被 CLI 覆盖后的配置副本。"""
    import copy
    updated = copy.deepcopy(bcfg)
    updated["output"]["directory"] = str(directory)
    return updated
