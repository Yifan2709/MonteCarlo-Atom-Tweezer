"""Level 0 静态高斯光镊命令行入口。"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .config import config_as_dict, load_config, with_output_override
from .constants import BOLTZMANN_CONSTANT
from .io_utils import ensure_output_directory, save_csv, save_json, save_markdown, save_yaml, validate_pngs
from .simulation import build_summary, run_simulation
from .visualization import (
    plot_energy_evolution,
    plot_force_comparison,
    plot_frequency_comparison,
    plot_potential_comparison,
    plot_single_trajectory,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="一维单原子静态高斯光镊 Level 0 数值检查")
    parser.add_argument("--config", required=True, help="YAML 配置路径")
    parser.add_argument("--trap", choices=("slm", "aod", "both"), default="both")
    parser.add_argument("--output-dir", help="覆盖输出目录；相对路径按项目根目录解析")
    return parser


def _write_outputs(config, simulation: dict) -> Path:
    output = ensure_output_directory(config.output.resolved_directory, config.project_root)
    save_yaml(output / "config_used.yaml", config_as_dict(config))
    grid = simulation["grid_m"]
    if config.output.save_csv:
        static_columns = {"x_m": grid, "x_um": grid * 1e6}
        for key, result in simulation["traps"].items():
            static_columns[f"{key}_potential_J"] = result["arrays"]["potential_grid_J"]
            static_columns[f"{key}_potential_uK"] = result["arrays"]["potential_grid_J"] / BOLTZMANN_CONSTANT * 1e6
            static_columns[f"{key}_force_N"] = result["arrays"]["force_grid_N"]
        save_csv(output / "static_grid.csv", static_columns)
        for key, result in simulation["traps"].items():
            arrays = result["arrays"]
            center = result["config"].center_m
            trajectory_columns = {
                "time_s": arrays["time_s"],
                "time_us": arrays["time_s"] * 1e6,
                "position_m": arrays["position_m"],
                "displacement_from_center_m": arrays["position_m"] - center,
                "position_um": arrays["position_m"] * 1e6,
                "velocity_m_per_s": arrays["velocity_m_per_s"],
                "acceleration_m_per_s2": arrays["acceleration_m_per_s2"],
                "kinetic_energy_J": arrays["kinetic_energy_J"],
                "potential_energy_J": arrays["potential_energy_J"],
                "total_energy_J": arrays["total_energy_J"],
                "energy_error_J": arrays["energy_error_J"],
                "energy_error_over_depth": arrays["energy_error_over_depth"],
                "abs_energy_error_over_depth": abs(arrays["energy_error_over_depth"]),
            }
            save_csv(output / f"{key}_trajectory.csv", trajectory_columns)

    records: list[dict] = []
    if config.output.save_png:
        records.append(plot_potential_comparison(grid, simulation["traps"], output / "potential_comparison.png"))
        records.append(plot_force_comparison(grid, simulation["traps"], output / "force_comparison.png"))
        for key, result in simulation["traps"].items():
            records.append(plot_single_trajectory(key, result, output / f"{key}_trajectory.png"))
            records.append(plot_energy_evolution(key, result, output / f"{key}_energy.png"))
        records.append(plot_frequency_comparison(simulation["traps"], output / "frequency_comparison.png"))
    image_validation = validate_pngs(records)

    metrics = {
        "model": {
            "coordinate": simulation["coordinate"],
            "coordinate_definition": "one radial coordinate in the Gaussian-beam focal plane, not the optical axis",
            "scope_and_limitations": simulation["scope_and_limitations"],
        },
        "selection": simulation["selection"],
        "trap_results": {key: result["metrics"] for key, result in simulation["traps"].items()},
        "image_validation_all_passed": image_validation["all_passed"],
        "visual_review_performed": image_validation["visual_review_performed"],
        "overall_validation_passed": image_validation["all_passed"] and all(
            result["metrics"]["validation"]["all_passed"] for result in simulation["traps"].values()
        ),
    }
    if config.output.save_json:
        save_json(output / "metrics.json", metrics)
        save_json(output / "image_validation.json", image_validation)
    save_markdown(output / "summary.md", build_summary(simulation, image_validation))
    return output


def main(argv: Sequence[str] | None = None) -> int:
    """加载配置、运行模拟并写出经验证的全部结果。"""

    args = _parser().parse_args(argv)
    config = load_config(Path(args.config))
    if args.output_dir:
        config = with_output_override(config, args.output_dir)
    simulation = run_simulation(config, args.trap)
    output = _write_outputs(config, simulation)
    print(f"Level 0 simulation completed: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

