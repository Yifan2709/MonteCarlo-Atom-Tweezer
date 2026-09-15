from pathlib import Path
import pytest
from simulation_workflow.cli import LEVELS, PROJECT_ROOT, command_for, inventory


def test_all_registered_level_files_exist():
    assert all(row["runnable_files_present"] for row in inventory())


@pytest.mark.parametrize("level", ["2c", "3", "4", "5", "full"])
def test_dispatch_uses_correct_module_and_fresh_output(level):
    command = command_for(level, "preflight")
    assert command[2] == LEVELS[level]["module"]
    assert Path(command[4]) == PROJECT_ROOT / "configs" / LEVELS[level]["config"]
    assert Path(command[-1]).name.startswith("run_")
    assert Path(command[-1]).parent == PROJECT_ROOT / "outputs" / LEVELS[level]["output"]


def test_invalid_cross_level_stage_is_rejected():
    with pytest.raises(ValueError):
        command_for("3", "interleaved_rb")
