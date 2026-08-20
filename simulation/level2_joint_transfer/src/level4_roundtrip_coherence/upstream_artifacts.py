"""发现、hash 并冻结 Level 2/3 上游产物，生成 upstream_manifest.json。"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import yaml

from level0_static_trap.io_utils import save_json

UPSTREAM_KEYS = ("level2_best_waveform", "level2_metrics", "level2_config_used",
                 "level3_metrics", "level3_config_used")


def _sha256(path: Path) -> str:
    """文件 SHA-256。"""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _record(path: Path, role: str, note: str) -> dict:
    """单个上游文件的 manifest 记录。"""
    stat = path.stat()
    return {
        "role": role,
        "absolute_path": str(path.resolve()),
        "sha256": _sha256(path),
        "size_bytes": int(stat.st_size),
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "note": note,
    }


def _discover(project_root: Path) -> dict[str, Path]:
    """在工作区中搜索 Level 2/3 冻结产物（不假定固定相对路径）。"""
    found: dict[str, Path] = {}
    candidates = {
        "level2_best_waveform": "outputs/level2_joint_transfer/demo/best_waveform.yaml",
        "level2_metrics": "outputs/level2_joint_transfer/demo/metrics.json",
        "level2_config_used": "outputs/level2_joint_transfer/demo/config_used.yaml",
        "level3_metrics": "outputs/level3_3d_transport_lensing/demo/metrics.json",
        "level3_config_used": "outputs/level3_3d_transport_lensing/demo/config_used.yaml",
    }
    for key, rel in candidates.items():
        direct = project_root / rel
        if direct.is_file():
            found[key] = direct
            continue
        hits = sorted(project_root.glob(f"**/{Path(rel).name}"))
        hits = [h for h in hits if ".venv" not in h.parts and ".git" not in h.parts]
        if hits:
            found[key] = hits[0]
    return found


def load_level3_frozen(cfg) -> dict:
    """从 Level 3 metrics.json 提取冻结的 lensing 场景与验证条件。"""
    path = _discover(cfg.project_root).get("level3_metrics")
    if path is None:
        raise FileNotFoundError("未找到 Level 3 metrics.json；请先完成并冻结 Level 3")
    data = json.loads(path.read_text(encoding="utf-8"))
    vs_table = data["lensing"]["fixed_vs_m_per_s"]
    nominal_vs = float(vs_table["1.0"]) if "1.0" in vs_table else None
    return {
        "metrics_path": path,
        "calibration_mode": data["lensing"]["calibration_mode"],
        "fixed_vs_m_per_s": {k: float(v) for k, v in vs_table.items()},
        "nominal_vs_m_per_s": nominal_vs,
        "frozen_conditions": data.get("frozen_validation_conditions", []),
        "preflight_all_passed": bool(data.get("preflight", {}).get("all_passed", False)),
    }


def load_level2_waveform_spec(cfg) -> dict:
    """读取 Level 2 冻结最佳波形 spec（含 provenance）。"""
    path = _discover(cfg.project_root).get("level2_best_waveform")
    if path is None:
        raise FileNotFoundError("未找到 Level 2 best_waveform.yaml；请先完成并冻结 Level 2")
    spec = yaml.safe_load(path.read_text(encoding="utf-8"))
    if spec.get("type") != "overlapped":
        raise ValueError(f"Level 4 仅适配 overlapped 冻结波形，得到 {spec.get('type')!r}")
    return {"path": path, "spec": spec,
            "duration_us": float(spec["duration_us"]),
            "provenance": spec.get("provenance", {})}


def build_manifest(cfg) -> tuple[dict, dict]:
    """生成 upstream manifest 与解析后的上游参数。

    返回 (manifest_dict, resolved)，resolved 含 nominal v_s、Level 2 波形 spec、
    各锚点参数及其来源标注（paper_anchor / level2_frozen / level3_frozen / level4_assumption）。
    """
    discovered = _discover(cfg.project_root)
    manifest_files = {}
    roles = {
        "level2_best_waveform": "Level 2 最终 transfer 波形（Protocol B pickup/dropoff）",
        "level2_metrics": "Level 2 配置/验证摘要与随机种子记录",
        "level2_config_used": "Level 2 解析配置（代码版本与冻结条件）",
        "level3_metrics": "Level 3 三维 lensing 配置与冻结验证摘要",
        "level3_config_used": "Level 3 解析配置",
    }
    missing = []
    for key in UPSTREAM_KEYS:
        path = discovered.get(key)
        if path is None:
            missing.append(key)
            continue
        manifest_files[key] = _record(path, key, roles[key])
    if missing:
        raise FileNotFoundError(f"上游冻结产物缺失: {missing}；停止主 validation")

    level3 = load_level3_frozen(cfg)
    level2 = load_level2_waveform_spec(cfg)
    nominal_vs = (cfg.nominal_vs_m_per_s if cfg.nominal_vs_m_per_s is not None
                  else level3["nominal_vs_m_per_s"])
    if nominal_vs is None:
        raise ValueError("无法确定 nominal v_s（Level 3 冻结表缺失 severity 1.0）")

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "require_sha256_manifest": cfg.require_sha256_manifest,
        "reconstructed_from_config": False,
        "files": manifest_files,
        "level3_lensing": {
            "calibration_mode": level3["calibration_mode"],
            "fixed_vs_m_per_s": level3["fixed_vs_m_per_s"],
            "nominal_vs_m_per_s": nominal_vs,
            "scenario_note": "v_s 沿用 Level 3 normalized_reference_scan 冻结值，"
                             "无量纲场景参数，非论文测量值",
            "level3_preflight_all_passed": level3["preflight_all_passed"],
        },
        "level2_waveform": {
            "type": level2["spec"]["type"],
            "name": level2["spec"].get("name"),
            "duration_us": level2["duration_us"],
            "params": {k: level2["spec"][k] for k in
                       ("move_start_fraction", "move_end_fraction",
                        "ramp_start_fraction", "ramp_end_fraction", "ramp_power")
                       if k in level2["spec"]},
            "provenance": level2["provenance"],
        },
        "parameter_provenance": {
            "durations_us": "paper_anchor (Extended Data Fig. 10e 组合操作锚点)",
            "depths_uK": "paper_anchor + level4_assumption (ramp 形状 smootherstep 为模型假设)",
            "split_distance_um": "level2_frozen (2.4 μm 初始 AOD-SLM 分离)",
            "long_profile": f"level3_frozen ({cfg.long_profile})",
            "split_profile": f"level3_frozen ({cfg.split_profile})",
            "vs_m_per_s": "level3_frozen (severity 1.0 场景)",
            "eta_slm/eta_aod": "level4_assumption (nominal 1.3e-4，敏感性对照 1.5e-4 等)",
            "protocol_b_slm_depth_uK": "level2_frozen (恒定 140 μK)",
        },
    }
    resolved = {
        "nominal_vs_m_per_s": nominal_vs,
        "level2_spec": dict(level2["spec"]),
        "level2_duration_us": level2["duration_us"],
        "level3_fixed_vs": level3["fixed_vs_m_per_s"],
    }
    return manifest, resolved


def save_manifest(path: Path, manifest: dict) -> None:
    """保存 upstream_manifest.json（严格 JSON，禁止 NaN）。"""
    save_json(path, manifest)
