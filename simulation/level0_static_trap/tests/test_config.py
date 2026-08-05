"""配置校验测试（spec §9 第 12 项）。"""

from __future__ import annotations

import math

import pytest

from level0_static_trap.config import ConfigError, load_config_from_dict


def _set(raw, path, value):
    """在嵌套字典中按点路径赋值。"""
    keys = path.split(".")
    d = raw
    for k in keys[:-1]:
        d = d[k]
    d[keys[-1]] = value


def test_default_config_loads(cfg):
    assert cfg.mass_kg > 0
    assert cfg.slm.depth_J > 0
    assert cfg.aod.depth_J > 0
    assert cfg.output_dir.endswith("out")


def test_rejects_non_positive_mass(raw_config, tmp_path):
    _set(raw_config, "atom.mass_u", 0.0)
    with pytest.raises(ConfigError, match="mass_u"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_rejects_non_positive_depth(raw_config, tmp_path):
    _set(raw_config, "slm.depth_uK", -10.0)
    with pytest.raises(ConfigError, match="depth_uK"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_rejects_non_positive_waist(raw_config, tmp_path):
    _set(raw_config, "aod.waist_um", 0.0)
    with pytest.raises(ConfigError, match="waist_um"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_rejects_non_finite_value(raw_config, tmp_path):
    _set(raw_config, "atom.mass_u", float("nan"))
    with pytest.raises(ConfigError, match="有限"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_rejects_unordered_grid(raw_config, tmp_path):
    _set(raw_config, "grid.min_offset_waist", 4.0)
    _set(raw_config, "grid.max_offset_waist", -4.0)
    with pytest.raises(ConfigError, match="偏移范围无序"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_rejects_too_few_grid_points(raw_config, tmp_path):
    _set(raw_config, "grid.num_points", 3)
    with pytest.raises(ConfigError, match="num_points"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_rejects_non_integer_grid_points(raw_config, tmp_path):
    _set(raw_config, "grid.num_points", 4001.5)
    with pytest.raises(ConfigError, match="num_points.*整数"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_rejects_oversized_dt(raw_config, tmp_path):
    # omega*dt 超过阈值：把 dt 调大
    _set(raw_config, "trajectory.dt_us", 50.0)  # 默认 0.05 -> 50，显著超阈值
    with pytest.raises(ConfigError, match="omega"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_rejects_non_positive_duration(raw_config, tmp_path):
    _set(raw_config, "trajectory.duration_us", 0.0)
    with pytest.raises(ConfigError, match="duration_us"):
        load_config_from_dict(raw_config, output_dir_override=str(tmp_path))


def test_aod_waist_marked_assumed(cfg):
    assert cfg.aod.waist_assumed is True
    assert cfg.slm.waist_assumed is False


def test_unknown_top_key_warns(raw_config, tmp_path):
    raw_config["bogus_section"] = {"a": 1}
    cfg = load_config_from_dict(raw_config, output_dir_override=str(tmp_path / "o"))
    assert any("bogus_section" in w for w in cfg.warnings)


def test_unknown_trap_field_warns(raw_config, tmp_path):
    raw_config["slm"]["extra_field"] = 123
    cfg = load_config_from_dict(raw_config, output_dir_override=str(tmp_path / "o"))
    assert any("extra_field" in w for w in cfg.warnings)
