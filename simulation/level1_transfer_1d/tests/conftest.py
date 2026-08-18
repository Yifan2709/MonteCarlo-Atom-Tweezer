from __future__ import annotations

from pathlib import Path

import pytest
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "default.yaml"


@pytest.fixture
def default_mapping() -> dict:
    with DEFAULT_CONFIG.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


@pytest.fixture
def write_config(tmp_path):
    def _write(mapping: dict) -> Path:
        path = tmp_path / "config.yaml"
        path.write_text(yaml.safe_dump(mapping, allow_unicode=True, sort_keys=False), encoding="utf-8")
        return path

    return _write

