"""冻结波形后的鲁棒性检查：温度、末端对准误差与时间步长敏感性。"""
from __future__ import annotations

import numpy as np

from .config import Level2Config
from .level2_simulation import (evaluate_waveform_on_states, evaluate_waveform_parallel,
                                physics_from_config, sample_initial_states)
from .noise_heating import MeanHeating
from .statistics import summarize_capture


def temperature_scan(cfg: Level2Config, specs: dict, heating: MeanHeating | None = None) -> list[dict]:
    """在多个初态温度下比较三类主波形（独立样本，含 Wilson 区间）。"""
    physics = physics_from_config(cfg)
    rows = []
    for index, temperature in enumerate(cfg.robustness.temperatures_uK):
        states = sample_initial_states(
            cfg, cfg.initial_ensemble.validation_seed + 100 + index,
            cfg.robustness.shots_per_condition, temperature_uK=temperature)
        for name, spec in specs.items():
            evaluation = evaluate_waveform_parallel(
                dict(spec, name=name), states, physics, cfg.integration.validation_dt_s,
                cfg.integration.post_transfer_hold_s, heating=heating)
            summary = summarize_capture(evaluation["captured_flags"], name)
            margins = [r["capture_margin"] for r in evaluation["records"] if r["captured"]]
            rows.append({
                "scan": "temperature", "condition": temperature, "waveform": name,
                "shots": summary["shots"], "captured": summary["captured"],
                "capture_fraction": summary["capture_fraction"],
                "wilson_low": summary["wilson_low"], "wilson_high": summary["wilson_high"],
                "mean_capture_margin": float(np.mean(margins)) if margins else None,
                "label_flips_vs_reference": "", "final_energy_max_abs_diff_uK": "",
            })
    return rows


def alignment_scan(cfg: Level2Config, specs: dict, heating: MeanHeating | None = None) -> list[dict]:
    """末端对准偏差敏感性：AOD 中心整体平移 offset（初态跟随实际初始中心采样）。"""
    physics = physics_from_config(cfg)
    rows = []
    for index, offset_um in enumerate(cfg.robustness.final_alignment_offsets_um):
        offset_m = offset_um * 1e-6
        states = sample_initial_states(
            cfg, cfg.initial_ensemble.validation_seed + 200 + index,
            cfg.robustness.shots_per_condition,
            aod_initial_center_m=physics["aod_initial_center_m"] + offset_m)
        for name, spec in specs.items():
            evaluation = evaluate_waveform_parallel(
                dict(spec, name=name), states, physics, cfg.integration.validation_dt_s,
                cfg.integration.post_transfer_hold_s, alignment_offset_m=offset_m, heating=heating)
            summary = summarize_capture(evaluation["captured_flags"], name)
            margins = [r["capture_margin"] for r in evaluation["records"] if r["captured"]]
            rows.append({
                "scan": "alignment", "condition": offset_um, "waveform": name,
                "shots": summary["shots"], "captured": summary["captured"],
                "capture_fraction": summary["capture_fraction"],
                "wilson_low": summary["wilson_low"], "wilson_high": summary["wilson_high"],
                "mean_capture_margin": float(np.mean(margins)) if margins else None,
                "label_flips_vs_reference": "", "final_energy_max_abs_diff_uK": "",
            })
    return rows


def timestep_scan(cfg: Level2Config, specs: dict, validation_states: dict,
                  heating: MeanHeating | None = None) -> list[dict]:
    """步长敏感性：固定 validation 子集，比较 0.10/0.05/0.025 us。

    连续量收敛用与最细步长的末态能量差表示；分类跳变用与最细步长的标签翻转率。
    """
    physics = physics_from_config(cfg)
    subset = {key: value[:cfg.robustness.timestep_subset_shots] for key, value in validation_states.items()
              if isinstance(value, np.ndarray)}
    rows = []
    for name, spec in specs.items():
        reference = evaluate_waveform_on_states(
            dict(spec, name=name), subset, physics, cfg.integration.convergence_dt_s,
            cfg.integration.post_transfer_hold_s, heating=heating)
        reference_energy = np.array([r["final_slm_energy_uK"] for r in reference["records"]])
        reference_flags = np.array(reference["captured_flags"])
        for dt in (cfg.integration.optimization_dt_s, cfg.integration.validation_dt_s,
                   cfg.integration.convergence_dt_s):
            evaluation = evaluate_waveform_on_states(
                dict(spec, name=name), subset, physics, dt, cfg.integration.post_transfer_hold_s,
                heating=heating)
            flags = np.array(evaluation["captured_flags"])
            energy = np.array([r["final_slm_energy_uK"] for r in evaluation["records"]])
            rows.append({
                "scan": "timestep", "condition": dt * 1e6, "waveform": name,
                "shots": len(flags), "captured": int(flags.sum()),
                "capture_fraction": float(flags.mean()),
                "wilson_low": "", "wilson_high": "",
                "mean_capture_margin": None,
                "label_flips_vs_reference": float(np.mean(flags != reference_flags)),
                "final_energy_max_abs_diff_uK": float(np.max(np.abs(energy - reference_energy))),
            })
    return rows
