from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import pytest
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("selection", ["slm", "aod", "both"])
def test_cli_modes_generate_only_selected_trap_and_valid_images(tmp_path, default_mapping, selection):
    data = default_mapping
    data["trajectory"].update(duration_us=120.0, dt_us=0.1, minimum_complete_cycles=2)
    config_path = tmp_path / f"{selection}.yaml"
    config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    output = tmp_path / selection
    completed = subprocess.run(
        [
            sys.executable, "-m", "level0_static_trap.cli", "--config", str(config_path),
            "--trap", selection, "--output-dir", str(output),
        ],
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    expected_common = {
        "config_used.yaml", "static_grid.csv", "metrics.json", "image_validation.json", "summary.md",
        "potential_comparison.png", "force_comparison.png", "frequency_comparison.png",
    }
    assert expected_common <= {path.name for path in output.iterdir()}
    for key in ("slm", "aod"):
        selected = selection == "both" or selection == key
        assert (output / f"{key}_trajectory.csv").exists() is selected
        assert (output / f"{key}_trajectory.png").exists() is selected
        assert (output / f"{key}_energy.png").exists() is selected

    metrics_text = (output / "metrics.json").read_text(encoding="utf-8")
    assert "NaN" not in metrics_text and "Infinity" not in metrics_text
    metrics = json.loads(metrics_text, parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
    assert metrics["overall_validation_passed"]
    assert set(metrics["trap_results"]) == ({"slm", "aod"} if selection == "both" else {selection})
    image_validation = json.loads((output / "image_validation.json").read_text(encoding="utf-8"))
    assert image_validation["all_passed"]
    assert all(item["passed"] for item in image_validation["images"].values())


def test_config_used_contains_cli_output_override(tmp_path, default_mapping):
    default_mapping["trajectory"].update(duration_us=120.0, dt_us=0.1, minimum_complete_cycles=2)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(default_mapping, sort_keys=False), encoding="utf-8")
    output = tmp_path / "override"
    completed = subprocess.run(
        [sys.executable, "-m", "level0_static_trap.cli", "--config", str(config_path), "--trap", "slm", "--output-dir", str(output)],
        cwd=PROJECT_ROOT, text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 0, completed.stderr
    used = yaml.safe_load((output / "config_used.yaml").read_text(encoding="utf-8"))
    assert used["output"]["directory"] == str(output)
    assert Path(used["output"]["resolved_directory"]) == output.resolve()

