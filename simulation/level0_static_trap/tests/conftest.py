"""共享 pytest fixtures。"""

from __future__ import annotations

import copy
import os
from pathlib import Path

import pytest

from level0_static_trap.config import load_config_from_dict

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "configs" / "default.yaml"


def _base_raw():
    """返回默认配置字典的深拷贝（供测试修改）。"""
    import yaml

    with open(DEFAULT_CONFIG_PATH, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


@pytest.fixture(scope="session")
def default_raw():
    return _base_raw()


@pytest.fixture
def raw_config(default_raw):
    """每次测试一份独立拷贝。"""
    return copy.deepcopy(default_raw)


@pytest.fixture
def cfg(default_raw, tmp_path):
    """加载默认配置，输出目录指向临时目录。"""
    raw = copy.deepcopy(default_raw)
    return load_config_from_dict(raw, output_dir_override=str(tmp_path / "out"))


@pytest.fixture
def project_root():
    return PROJECT_ROOT
