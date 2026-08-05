"""CLI 冒烟测试（spec §9 第 13–15 项）。

通过 subprocess 调用 ``python -m level0_static_trap.cli``，验证退出码为 0、
仅生成所选模式对应的文件、metrics.json 可严格解析且不含 NaN/Inf，所有 PNG
通过程序化图像检查（image_validation.json 中记录）。
"""

from __future__ import annotations

import io
import json
import os
import re
import subprocess
import sys

import pytest

EXPECTED_BOTH = {
    "config_used.yaml",
    "static_grid.csv",
    "slm_trajectory.csv",
    "aod_trajectory.csv",
    "metrics.json",
    "image_validation.json",
    "summary.md",
    "potential_comparison.png",
    "force_comparison.png",
    "slm_trajectory.png",
    "aod_trajectory.png",
    "slm_energy.png",
    "aod_energy.png",
    "frequency_comparison.png",
}


def _run_cli(config, trap, out_dir, project_root):
    env = dict(os.environ)
    proc = subprocess.run(
        [sys.executable, "-m", "level0_static_trap.cli",
         "--config", str(config), "--trap", trap, "--output-dir", str(out_dir)],
        cwd=str(project_root),
        env=env,
        capture_output=True,
        text=True,
    )
    return proc


def _assert_no_nan_inf_in_json(obj):
    """递归检查 JSON 对象中不含 NaN/Infinity（已由 allow_nan=False 保证，这里双保险）。"""
    if isinstance(obj, dict):
        for v in obj.values():
            _assert_no_nan_inf_in_json(v)
    elif isinstance(obj, list):
        for v in obj:
            _assert_no_nan_inf_in_json(v)
    elif isinstance(obj, float):
        import math
        assert math.isfinite(obj), f"JSON 含非有限值：{obj}"


def test_cli_both_exit_zero_and_expected_files(project_root, tmp_path):
    proc = _run_cli(project_root / "configs" / "default.yaml", "both", tmp_path / "both", project_root)
    assert proc.returncode == 0, f"CLI 失败：\nstdout={proc.stdout}\nstderr={proc.stderr}"
    files = set(os.listdir(tmp_path / "both"))
    assert files == EXPECTED_BOTH, f"文件集合不符：多={files-EXPECTED_BOTH}，缺={EXPECTED_BOTH-files}"


def test_cli_slm_only_generates_slm_files(project_root, tmp_path):
    proc = _run_cli(project_root / "configs" / "default.yaml", "slm", tmp_path / "slm", project_root)
    assert proc.returncode == 0, f"CLI 失败：\nstderr={proc.stderr}"
    files = set(os.listdir(tmp_path / "slm"))
    # SLM 模式：必须有 slm_trajectory.csv，不得伪造 aod_trajectory.csv
    assert "slm_trajectory.csv" in files
    assert "aod_trajectory.csv" not in files
    # 单势阱模式仍生成势/力比较图、slm 轨迹与能量图、频率图
    assert "potential_comparison.png" in files
    assert "slm_trajectory.png" in files
    assert "slm_energy.png" in files
    assert "frequency_comparison.png" in files
    assert "aod_trajectory.png" not in files
    assert "aod_energy.png" not in files


def test_cli_aod_only_generates_aod_files(project_root, tmp_path):
    proc = _run_cli(project_root / "configs" / "default.yaml", "aod", tmp_path / "aod", project_root)
    assert proc.returncode == 0, f"CLI 失败：\nstderr={proc.stderr}"
    files = set(os.listdir(tmp_path / "aod"))
    assert "aod_trajectory.csv" in files
    assert "slm_trajectory.csv" not in files
    assert "aod_trajectory.png" in files
    assert "slm_trajectory.png" not in files


def test_metrics_json_strict_parse_no_nan_inf(project_root, tmp_path):
    proc = _run_cli(project_root / "configs" / "default.yaml", "both", tmp_path / "both2", project_root)
    assert proc.returncode == 0, proc.stderr
    with open(tmp_path / "both2" / "metrics.json", "r", encoding="utf-8") as fh:
        # json.load 默认允许 NaN/Infinity；先用 raw 文本检查
        text = fh.read()
    assert not re.search(r":\s*NaN", text, re.IGNORECASE)
    assert not re.search(r":\s*Infinity", text, re.IGNORECASE)
    obj = json.loads(text)  # 应可严格解析
    _assert_no_nan_inf_in_json(obj)
    # 结构：trap_choice 与 traps
    assert obj["trap_choice"] == "both"
    assert set(obj["traps"].keys()) == {"slm", "aod"}
    for name, m in obj["traps"].items():
        assert "frequencies" in m
        assert "error_classification" in m
        assert "energy" in m
        assert "validation" in m


def test_all_pngs_pass_programmatic_check(project_root, tmp_path):
    proc = _run_cli(project_root / "configs" / "default.yaml", "both", tmp_path / "both3", project_root)
    assert proc.returncode == 0, proc.stderr
    with open(tmp_path / "both3" / "image_validation.json", "r", encoding="utf-8") as fh:
        iv = json.load(fh)
    assert iv["all_pass"] is True
    for name, info in iv["images"].items():
        assert info["checks"]["overall_pass"] is True, f"{name} 程序化检查未通过：{info['checks']}"
        assert info["checks"]["size_above_lower_bound"] is True
        assert info["checks"]["reopen_ok"] is True


def test_cli_invalid_trap_choice_rejected(project_root, tmp_path):
    proc = _run_cli(project_root / "configs" / "default.yaml", "bogus", tmp_path / "x", project_root)
    assert proc.returncode != 0


def test_cli_oversized_dt_returns_nonzero(project_root, tmp_path):
    """构造一个 dt 过大的配置，CLI 应在配置校验阶段失败退出。"""
    import yaml

    with open(project_root / "configs" / "default.yaml", "r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    raw["trajectory"]["dt_us"] = 50.0  # 远超 omega*dt 阈值
    bad_cfg = tmp_path / "bad.yaml"
    with open(bad_cfg, "w", encoding="utf-8") as fh:
        yaml.safe_dump(raw, fh)
    proc = _run_cli(bad_cfg, "both", tmp_path / "badout", project_root)
    assert proc.returncode != 0
    assert "omega" in proc.stderr or "配置错误" in proc.stderr
