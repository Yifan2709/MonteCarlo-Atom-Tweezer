"""实验 5（最高速度测算）的测试：条件覆盖、同时长剖面 spec、端到端冒烟。"""
from __future__ import annotations

import copy
import math
from pathlib import Path

import numpy as np
import pytest

from level2b_speed_limit.level2b_config import load_level2b_config
from level2b_speed_limit.max_speed_probe import (adiabatic_speed_m_per_s,
                                                 condition_physics,
                                                 profile_spec_set,
                                                 run_max_speed_stage)
from level2b_speed_limit.speed_engine import physics_from_config
from level2b_speed_limit.speed_waveforms import build_speed_waveform, speed_calibers

PROJECT_ROOT = Path(__file__).resolve().parents[1]
LEVEL2B_CONFIG = PROJECT_ROOT / "configs" / "level2b_speed_limit.yaml"


@pytest.fixture(scope="module")
def setup():
    bcfg, cfg2 = load_level2b_config(LEVEL2B_CONFIG)
    return bcfg, cfg2, physics_from_config(cfg2)


def test_condition_physics_overrides_and_adiabatic_scaling(setup):
    _bcfg, _cfg2, physics = setup
    doubled = condition_physics(physics, depth_scale=2.0)
    assert math.isclose(doubled["aod_initial_depth_j"],
                        2.0 * physics["aod_initial_depth_j"], rel_tol=1e-15)
    assert doubled["aod_waist_m"] == physics["aod_waist_m"]
    narrowed = condition_physics(physics, waist_um=0.90)
    assert math.isclose(narrowed["aod_waist_m"], 0.90e-6, rel_tol=1e-15)
    assert narrowed["aod_initial_depth_j"] == physics["aod_initial_depth_j"]
    nominal = condition_physics(physics)
    assert nominal == physics  # 无覆盖时逐键一致
    # v_ad = w·ω，ω ∝ √D/w → v_ad ∝ √D：深度 ×2 → ×√2，束腰不改变 v_ad
    assert math.isclose(adiabatic_speed_m_per_s(doubled)
                        / adiabatic_speed_m_per_s(nominal), math.sqrt(2.0),
                        rel_tol=1e-12)
    assert math.isclose(adiabatic_speed_m_per_s(narrowed),
                        adiabatic_speed_m_per_s(nominal), rel_tol=1e-12)


def test_same_duration_profiles_share_duration_and_displacement(setup):
    _bcfg, _cfg2, physics = setup
    duration = 40e-6
    accel_by_profile = {"smootherstep": None, "trapezoid_30": 0.30,
                        "trapezoid_15": 0.15, "triangular_50": 0.50}
    specs = profile_spec_set(duration, accel_by_profile)
    assert set(specs) == set(accel_by_profile)
    displacements, durations, ratios = [], [], {}
    for name, spec in specs.items():
        waveform = build_speed_waveform(spec, physics)
        cal = speed_calibers(waveform)
        assert math.isclose(waveform.duration_s, duration, rel_tol=1e-12)
        displacements.append(abs(cal["total_displacement_m"]))
        durations.append(cal["move_duration_s"])
        ratios[name] = cal["v_peak_m_per_s"] / cal["v_avg_m_per_s"]
    # 同时长、同位移（= 同 v_avg），仅速度形状不同。
    # move_duration 口径差：smootherstep 取速度非零区间跨度、梯形取窗口宽度，
    # 相对差 ~3e-4，故用 1e-3 容差。
    assert np.allclose(displacements, displacements[0], rtol=1e-9)
    assert np.allclose(durations, durations[0], rtol=1e-3)
    # 峰值/平均比的理论值受 v_avg 口径影响有 ~3e-4 偏差
    assert math.isclose(ratios["smootherstep"], 1.875, rel_tol=1e-3)
    assert math.isclose(ratios["triangular_50"], 2.0, rel_tol=1e-3)
    assert math.isclose(ratios["trapezoid_30"], 1.0 / 0.70, rel_tol=1e-3)
    assert math.isclose(ratios["trapezoid_15"], 1.0 / 0.85, rel_tol=1e-3)


def test_run_max_speed_stage_smoke(setup, tmp_path):
    bcfg, cfg2, physics = setup
    probe_cfg = copy.deepcopy(bcfg)
    probe_cfg["max_speed_probe"] = {
        "temperature_uK": 20.0,
        "v95_structure": "frozen_shape_scaled",
        "v95_durations_us": [90, 45, 32],
        "profile_durations_us": [44],
        "shots_per_point": 8,
        "hardware_conditions": [{"name": "depth_x2.0", "aod_depth_scale": 2.0}],
    }
    summary = run_max_speed_stage(cfg2, probe_cfg, physics, tmp_path, num_procs=1)
    # 5a：nominal + 1 个硬件条件，各 3 点
    hardware_csv = tmp_path / "max_speed_hardware.csv"
    rows = hardware_csv.read_text(encoding="utf-8").strip().splitlines()
    assert len(rows) == 1 + 2 * 3
    assert {e["condition"] for e in summary["v95_table"]} == {"nominal", "depth_x2.0"}
    # 8 shots 下 v95 可能未覆盖，但表结构与旗标档案必须完整
    for entry in summary["v95_table"]:
        assert entry["boundary_covered"] in (True, False)
        assert entry["adiabatic_speed_mm_per_s"] > 0
    # 5b：1 个时长 × 4 剖面，配对差只对非参考剖面存在
    profiles_csv = tmp_path / "max_speed_profiles.csv"
    prows = profiles_csv.read_text(encoding="utf-8").strip().splitlines()
    assert len(prows) == 1 + 4
    names = {p["profile"] for p in summary["profiles"]}
    assert names == {"smootherstep", "trapezoid_30", "trapezoid_15", "triangular_50"}
    for item in summary["profiles"]:
        assert 0.0 <= item["capture_rate"] <= 1.0
        if item["profile"] == "trapezoid_30":
            assert item["pair_vs_reference_difference"] is None
        else:
            assert item["pair_vs_reference_difference"] is not None
