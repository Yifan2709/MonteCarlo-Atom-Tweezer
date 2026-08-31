from copy import deepcopy
from pathlib import Path

import pytest

from level0_static_trap.config import load_config, with_output_override
from level0_static_trap.io_utils import ensure_output_directory


def test_default_config_loads_and_resolves_against_project_root():
    config_path = Path(__file__).resolve().parents[1] / "configs" / "default.yaml"
    config = load_config(config_path)
    assert config.atom.mass_kg > 0.0
    assert config.slm.depth_j > 0.0
    assert config.aod.waist_is_assumed
    assert config.output.resolved_directory == (config.project_root / config.output.directory).resolve()


def test_cli_output_override_is_recorded_and_root_relative():
    config_path = Path(__file__).resolve().parents[1] / "configs" / "default.yaml"
    config = with_output_override(load_config(config_path), "outputs/custom")
    assert config.output.directory == "outputs/custom"
    assert config.output.resolved_directory == (config.project_root / "outputs/custom").resolve()


@pytest.mark.parametrize(
    "mutator, message",
    [
        (lambda data: data["slm"].update(depth_uK=0.0), "深度和束腰必须大于零"),
        (lambda data: data["aod"].update(waist_um=float("nan")), "有限数值"),
        (lambda data: data["atom"].update(mass_u=-1.0), "原子质量必须大于零"),
        (lambda data: data["grid"].update(num_points=4), "至少为 5"),
        (lambda data: data["grid"].update(min_offset_waist=2.0, max_offset_waist=-2.0), "下限必须小于上限"),
        (lambda data: data["trajectory"].update(initial_velocity_m_per_s=1.0), "初始状态不束缚"),
        (lambda data: data["trajectory"].update(dt_us=1.0), "时间步长过大"),
    ],
)
def test_invalid_physical_and_numerical_config_is_rejected(default_mapping, write_config, mutator, message):
    data = deepcopy(default_mapping)
    mutator(data)
    with pytest.raises(ValueError, match=message):
        load_config(write_config(data))


def test_unknown_key_is_rejected(default_mapping, write_config):
    data = deepcopy(default_mapping)
    data["slm"]["mystery"] = 12
    with pytest.raises(ValueError, match="未知字段"):
        load_config(write_config(data))


def test_too_few_cycles_is_rejected(default_mapping, write_config):
    data = deepcopy(default_mapping)
    data["trajectory"].update(duration_us=20.0, dt_us=0.05)
    with pytest.raises(ValueError, match="预计完整周期数不足"):
        load_config(write_config(data))


def test_output_inside_sources_is_rejected(tmp_path):
    root = tmp_path / "project"
    with pytest.raises(ValueError, match="sources"):
        ensure_output_directory(root / "sources" / "bad", root)

