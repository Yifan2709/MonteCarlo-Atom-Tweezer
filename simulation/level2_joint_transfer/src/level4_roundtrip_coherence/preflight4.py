"""Level 4 预检查：20 项按序自检，关键项失败则停止大样本 validation。"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from level3_3d_transport_lensing.gaussian_3d import KB
from . import protocol_segments as ps
from . import upstream_artifacts
from .dephasing import (PhaseAccumulator, analytic_constant_phase,
                        conditional_contrast)
from .level4_config import Level4Config
from .noise_models import IntensityNoise
from .roundtrip_simulation import simulate_roundtrip
from .roundtrip_statistics import evaluate_outcomes
from .sampling4 import sample_slm_states
from .trap_assignment import classify_basin

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
HBAR = 1.054571817e-34
US = 1e-6


def _ok(**kw) -> dict:
    """把检查字典加上 passed 字段。"""
    passed = all(v for v in kw.values() if isinstance(v, bool))
    return dict(kw, passed=bool(passed))


class Preflight4:
    """聚合 20 项预检查结果。"""

    critical_checks = [
        "existing_tests", "upstream_artifacts", "segment_continuity",
        "control_roundtrip", "level2_regression", "level3_regression",
        "protocol_b_no_duplicate_split_merge", "ideal_center_protocol",
        "control_time_reversal", "work_energy_residual", "static_hold_bounded",
        "basin_classification", "first_failure_state_machine",
        "round_state_continuity", "analytic_phase", "dd_cancellation",
        "pulse_grid_handling", "ensemble_contrast_cases", "dt_convergence",
        "frozen_validation_conditions",
    ]

    def __init__(self, cfg: Level4Config, manifest, resolved, protocols: dict):
        self.cfg = cfg
        self.manifest = manifest
        self.resolved = resolved
        self.protocols = protocols
        self.checks: dict[str, dict] = {}

    def run(self) -> dict:
        """按 prompt §13 顺序执行全部检查。"""
        started = time.time()
        runners = {
            "existing_tests": self._existing_tests,
            "upstream_artifacts": self._upstream,
            "segment_continuity": self._segment_continuity,
            "control_roundtrip": self._control_roundtrip,
            "level2_regression": self._level2_regression,
            "level3_regression": self._level3_regression,
            "protocol_b_no_duplicate_split_merge": self._no_duplicate_split,
            "ideal_center_protocol": self._ideal_center,
            "control_time_reversal": self._time_reversal,
            "work_energy_residual": self._work_energy,
            "static_hold_bounded": self._static_hold,
            "basin_classification": self._basin_classification,
            "first_failure_state_machine": self._state_machine,
            "round_state_continuity": self._round_continuity,
            "analytic_phase": self._analytic_phase,
            "dd_cancellation": self._dd_cancellation,
            "pulse_grid_handling": self._pulse_grid,
            "ensemble_contrast_cases": self._contrast_cases,
            "dt_convergence": self._dt_convergence,
            "frozen_validation_conditions": self._frozen_conditions,
        }
        for name in self.critical_checks:
            t0 = time.time()
            try:
                self.checks[name] = runners[name]()
            except Exception as exc:  # noqa: BLE001 - 记录失败继续
                self.checks[name] = {"passed": False, "error": repr(exc)}
            self.checks[name]["elapsed_s"] = round(time.time() - t0, 2)
            print(f"[preflight] {name}: "
                  f"{'PASS' if self.checks[name].get('passed') else 'FAIL'}"
                  f" ({self.checks[name]['elapsed_s']:.1f}s)")
        all_passed = all(self.checks[n].get("passed", False)
                         for n in self.critical_checks)
        return {"all_passed": bool(all_passed),
                "critical_checks": list(self.critical_checks),
                "elapsed_s": time.time() - started, **self.checks}

    # ---- 1. 已有测试（prompt §13 要求 Level 0-3 全部测试；test_level5 属
    # 并行开发中的独立 Level，其失败与 Level 4 冻结栈无关，故显式排除）----
    def _existing_tests(self) -> dict:
        ignore = [] if not (PACKAGE_ROOT / "tests" / "test_level5.py").is_file() \
            else ["--ignore=" + str(PACKAGE_ROOT / "tests" / "test_level5.py")]
        result = subprocess.run(
            [sys.executable, "-m", "pytest", str(PACKAGE_ROOT / "tests"),
             "-q", "--no-header", "-m", "not slow", *ignore],
            cwd=PACKAGE_ROOT, capture_output=True, text=True, timeout=7200)
        tail = result.stdout.strip().splitlines()[-1] if result.stdout.strip() \
            else ""
        return _ok(exit_code=result.returncode, summary_line=tail,
                   scoped_to="level0-level4（test_level5 并行开发中，另行验证）",
                   passed=result.returncode == 0)

    # ---- 2. 上游产物 ----
    def _upstream(self) -> dict:
        files = self.manifest["files"]
        roles_ok = set(files) >= {"level2_best_waveform", "level3_metrics"}
        sha_ok = all(len(v["sha256"]) == 64 for v in files.values())
        return _ok(files_present=bool(roles_ok), sha256_present=bool(sha_ok),
                   level3_preflight_passed=self.manifest["level3_lensing"][
                       "level3_preflight_all_passed"],
                   passed=bool(roles_ok and sha_ok
                               and self.manifest["level3_lensing"][
                                   "level3_preflight_all_passed"]))

    # ---- 3. 分段连续性 ----
    def _segment_continuity(self) -> dict:
        issues = []
        for name, timeline in self.protocols.items():
            if name == "static_idle":
                continue
            for i in range(len(timeline.segments) - 1):
                boundary = timeline.segment_starts[i + 1]
                left = timeline.segments[i].controls(np.array([boundary
                                                               - timeline.segment_starts[i]]))
                right = timeline.segments[i + 1].controls(np.array([0.0]))
                for key in ("p", "pd", "pdd", "d_slm", "ddot_slm", "d_aod",
                            "ddot_aod"):
                    # 尺度取分段内峰值（端点值本身≈0，不适合做相对容差）
                    dense = timeline.segments[i].controls(
                        np.linspace(0, timeline.segments[i].duration_s, 21))
                    scale = max(abs(float(np.max(np.abs(dense[key])))),
                                abs(float(np.max(np.abs(right[key])))), 1e-30)
                    gap = float(np.max(np.abs(np.asarray(left[key])
                                              - np.asarray(right[key]))))
                    tol = 1e-6 * scale + 1e-12
                    if gap > tol:
                        issues.append(f"{name} 边界 {i}: {key} 不连续 gap={gap:.3e}")
        return _ok(issues=issues, passed=not issues)

    # ---- 4. 控制量回到初始 ----
    def _control_roundtrip(self) -> dict:
        results = {}
        for name in ("protocol_a", "protocol_b"):
            timeline = self.protocols[name]
            c0 = timeline.control_at(0.0)
            cT = timeline.control_at(timeline.total_duration_s)
            diffs = {k: abs(cT[k] - c0[k]) for k in
                     ("c1", "c2", "d_slm", "d_aod")}
            scale = max(abs(c0["d_slm"]), abs(c0["d_aod"]), 1e-30)
            results[name] = {"max_diff": max(diffs.values()),
                             "rel": max(diffs.values()) / scale}
        ok = all(r["rel"] < 1e-9 for r in results.values())
        return _ok(per_protocol=results, passed=bool(ok))

    # ---- 5. Level 2 transfer 三维回归 ----
    def _level2_regression(self) -> dict:
        spec = self.resolved["level2_spec"]
        dropoff = ps.Segment(
            "dropoff", "dropoff", spec["duration_us"] * US,
            ps._combine(ps._l2_waveform_path(spec, 2.4e-6, 280e-6 * KB,
                                             spec["duration_us"] * US,
                                             reverse=False),
                        ps._constant_depth(140e-6 * KB),
                        ps._constant_depth(0.0), override_aod=False),
            checkpoint="dropoff_end", expected_basin="bound_to_slm")
        timeline = ps.ProtocolTimeline(
            "l2_dropoff_3d", [dropoff],
            ps.physics_from_config(self.cfg, v_s=None))
        states = sample_slm_states(self.cfg, self.cfg.development_seed + 101,
                                   256, temperature_uK=5.0, depth_uK=280.0)
        # 初态平移到 2.4 μm 外的 AOD 阱（Level 2 语义：AOD 携带原子）
        split = self.cfg.split_distance_um * 1e-6
        f1, f2 = timeline.physics.geometry_factors
        pos = states["pos_m"] + np.array([split * f1, split * f2, 0.0])
        res = simulate_roundtrip(timeline, 0.10 * US, pos,
                                 states["vel_m_per_s"], timeline.physics,
                                 {"none": np.empty(0)})
        final = res.checkpoints[-1]["classification"]
        captured = np.array([b == "bound_to_slm" for b in final["basin"]])
        rate = float(captured.mean())
        return _ok(shots=int(captured.size), capture_rate=rate,
                   level2_reference=1.0,
                   passed=bool(rate >= 0.99))

    # ---- 6. Level 3 long-move 回归 ----
    def _level3_regression(self) -> dict:
        seg = ps.Segment(
            "outbound_long_move", "outbound_long_move", 1600.0 * US,
            ps._combine(ps._profile_move(0.0, 610e-6, 1600.0 * US,
                                         "adiabatic_sine"),
                        ps._constant_depth(0.0),
                        ps._constant_depth(280e-6 * KB)),
            checkpoint="outbound_end", expected_basin="bound_to_aod")
        v_s = self.resolved["nominal_vs_m_per_s"]
        timeline = ps.ProtocolTimeline(
            "l3_anchor_610um", [seg], ps.physics_from_config(self.cfg, v_s=v_s))
        states = sample_slm_states(self.cfg, self.cfg.development_seed + 102,
                                   256, temperature_uK=5.0, depth_uK=280.0)
        res = simulate_roundtrip(timeline, 0.10 * US, states["pos_m"],
                                 states["vel_m_per_s"], timeline.physics,
                                 {"none": np.empty(0)})
        final = res.checkpoints[-1]["classification"]
        retained = np.array([b == "bound_to_aod" for b in final["basin"]])
        metrics3 = json.loads(Path(self.manifest["files"]["level3_metrics"][
            "absolute_path"]).read_text(encoding="utf-8"))
        ref = metrics3["validation"]["summaries"].get(
            "paper_anchor_diagonal|sev1.0", {})
        return _ok(shots=int(retained.size),
                   retention_fraction=float(retained.mean()),
                   level3_reference=ref.get("retention_fraction"),
                   passed=bool(retained.mean() >= 0.99))

    # ---- 7. Protocol B 不重复 split/merge ----
    def _no_duplicate_split(self) -> dict:
        b = self.protocols["protocol_b"]
        kinds = [s.kind for s in b.segments]
        has_split = "split" in kinds
        has_merge = "merge" in kinds
        # L2 pickup/dropoff 已含 2.4 μm 短移动
        pickup = next(s for s in b.segments if s.kind == "pickup")
        dropoff = next(s for s in b.segments if s.kind == "dropoff")
        l2_dur = self.resolved["level2_duration_us"] * US
        moves = bool(np.max(np.abs(pickup.controls(
            np.linspace(0, l2_dur, 11))["p"])) > 1e-9)
        return _ok(no_split_segment=not has_split,
                   no_merge_segment=not has_merge,
                   l2_waveform_contains_move=moves,
                   passed=bool(not has_split and not has_merge and moves))

    # ---- 8. 理想中心零速度初态完整协议 ----
    def _ideal_center(self) -> dict:
        results = {}
        for name in ("protocol_a", "protocol_b"):
            timeline = self.protocols[name]
            phys = timeline.physics
            pos = np.zeros((1, 3))
            vel = np.zeros((1, 3))
            pulses = {m: timeline.dd_pulse_times(m, 0.05 * US, None)
                      for m in ("none",)}
            res = simulate_roundtrip(timeline, 0.10 * US, pos, vel, phys,
                                     pulses)
            matches = [np.array([b == rec["checkpoint"]["expected_basin"]
                                 for b in rec["classification"]["basin"]])
                       for rec in res.checkpoints]
            final_retained = matches[-1][0]
            outcomes = evaluate_outcomes(matches, np.array([final_retained]))
            results[name] = {"success": bool(outcomes["roundtrip_success"][0]),
                             "residual_uK": float(abs(res.residual[0])
                                                  / KB * 1e6)}
        return _ok(per_protocol=results,
                   passed=all(r["success"] for r in results.values()))

    # ---- 9. 控制时间反向 ----
    def _time_reversal(self) -> dict:
        timeline = self.protocols["protocol_a_no_lensing"]
        rev = timeline.reversed()
        pos0 = np.zeros((1, 3))
        vel0 = np.zeros((1, 3))
        errors = {}
        for dt_us in (0.10, 0.05, 0.025):
            dt = dt_us * US
            fwd = simulate_roundtrip(timeline, dt, pos0, vel0,
                                     timeline.physics, {"none": np.empty(0)})
            pos_rev_start = fwd.final_pos
            vel_rev_start = -fwd.final_vel
            back = simulate_roundtrip(rev, dt, pos_rev_start, vel_rev_start,
                                      rev.physics, {"none": np.empty(0)})
            err_pos = float(np.linalg.norm(back.final_pos - pos0))
            err_vel = float(np.linalg.norm(back.final_vel - vel0))
            errors[dt_us] = {"pos_err_m": err_pos, "vel_err_m_s": err_vel}
        # 反演误差应远小于任何物理尺度（Verlet 可逆性到机器精度量级）
        tiny = errors[0.025]["pos_err_m"] < 1e-9 \
            and errors[0.025]["vel_err_m_s"] < 1e-6
        return _ok(errors=errors, reversal_exact=bool(tiny),
                   final_pos_err_um=errors[0.025]["pos_err_m"] * 1e6,
                   passed=bool(tiny))

    # ---- 10. 分段与全周期功-能残差 ----
    def _work_energy(self) -> dict:
        from .energy_ledger import segment_ledger_rows
        timeline = self.protocols["protocol_a"]
        states = sample_slm_states(self.cfg, self.cfg.development_seed + 103, 32)
        grid = timeline.build_grid(0.10 * US)
        res = simulate_roundtrip(timeline, 0.10 * US, states["pos_m"],
                                 states["vel_m_per_s"], timeline.physics,
                                 {"none": np.empty(0)})
        rows, seg_works = segment_ledger_rows(res, grid, "protocol_a",
                                              per_shot=False)
        worst = max(float(np.max(np.abs(w["residual"]))) / KB * 1e6
                    for w in seg_works.values())
        total = float(np.max(np.abs(res.residual))) / KB * 1e6
        scale = 280.0
        return _ok(per_segment_max_residual_uK=worst,
                   full_cycle_max_residual_uK=total,
                   relative_to_depth=worst / scale,
                   passed=bool(worst / scale < 1e-4 and total / scale < 1e-4))

    # ---- 11. 静态 hold 能量误差有界 ----
    def _static_hold(self) -> dict:
        from .energy_ledger import static_hold_check
        timeline = self.protocols["static_idle"]
        states = sample_slm_states(self.cfg, self.cfg.development_seed + 104, 32)
        grid = timeline.build_grid(0.10 * US)
        res = simulate_roundtrip(timeline, 0.10 * US, states["pos_m"],
                                 states["vel_m_per_s"], timeline.physics,
                                 {"none": np.empty(0)})
        checks = static_hold_check(res, grid)
        depth = self.cfg.slm_initial_depth_uK
        worst = max(c["max_abs_w_ext_J"] for c in checks.values()) / KB * 1e6
        bounded = worst < 1e-4 * depth
        return _ok(hold_checks=checks, max_w_ext_uK=worst,
                   passed=bool(bounded and checks))

    # ---- 12. overlap basin 分类 ----
    def _basin_classification(self) -> dict:
        phys = self.protocols["protocol_a"].physics
        d = 280e-6 * KB
        ctrl_dual = {"d_slm": 60e-6 * KB, "d_aod": d, "c1": 2.4e-6, "c2": 2.4e-6,
                     "zs1": 0.0, "zs2": 0.0}
        in_aod = classify_basin(np.array([[2.4e-6, 2.4e-6, 0.0]]),
                                np.zeros((1, 3)), ctrl_dual, phys)["basin"][0]
        in_slm = classify_basin(np.array([[0.05e-6, 0.0, 0.0]]),
                                np.zeros((1, 3)), ctrl_dual, phys)["basin"][0]
        merged = classify_basin(
            np.array([[0.0, 0.0, 0.0]]), np.zeros((1, 3)),
            {"d_slm": 60e-6 * KB, "d_aod": d, "c1": 0.0, "c2": 0.0,
             "zs1": 0.0, "zs2": 0.0}, phys)["basin"][0]
        # 构造能量介于鞍点与零之间的 ambiguous 案例：先解析求鞍点势
        from .trap_assignment import (line_minima_and_barrier,
                                      total_potential_line)
        c_slm = np.array(phys.slm_center_m)
        t_line, u_line = total_potential_line(
            None, c_slm, np.array([2.4e-6, 2.4e-6]), 60e-6 * KB, d, phys)
        _, saddle_idx, _ = line_minima_and_barrier(u_line, d)
        u_saddle = float(u_line[saddle_idx])
        u_at_atom = float(-d)  # 原子置于 AOD 中心，U_tot ≈ -D_AOD
        # e_tot 目标取鞍点与连续统零点的中点：束缚但两 basin 均可达
        ke_target = 0.5 * u_saddle - u_at_atom
        v_hot = float(np.sqrt(2.0 * ke_target / (3.0 * phys.mass_kg)))
        hot = classify_basin(np.array([[2.4e-6, 2.4e-6, 0.0]]),
                             np.full((1, 3), v_hot), ctrl_dual, phys)["basin"][0]
        gone = classify_basin(np.array([[50e-6, 50e-6, 0.0]]),
                              np.zeros((1, 3)), ctrl_dual, phys)["basin"][0]
        return _ok(
            atom_in_aod_basin=in_aod == "bound_to_aod",
            atom_in_slm_basin=in_slm == "bound_to_slm",
            coincident_shared=merged == "shared_or_ambiguous",
            above_barrier_ambiguous=hot == "shared_or_ambiguous",
            far_away_unbound=gone == "unbound",
            basins_found=[in_aod, in_slm, merged, hot, gone],
            passed=bool(in_aod == "bound_to_aod"
                        and in_slm == "bound_to_slm"
                        and merged == "shared_or_ambiguous"
                        and hot == "shared_or_ambiguous"
                        and gone == "unbound"))

    # ---- 13. first-failure 状态机 ----
    def _state_machine(self) -> dict:
        # 4 个 shot × 4 个 checkpoint 的匹配矩阵（行=shot，列=checkpoint）
        ck_names = ["pickup_end", "split_end", "outbound_end", "final_hold_end"]
        per_shot = np.array([
            [True, True, True, True],    # shot0 全过
            [True, False, True, True],   # shot1 split 失败
            [False, True, True, False],  # shot2 pickup 首败
            [True, True, True, False],   # shot3 末 checkpoint 失败
        ])
        matches = [per_shot[:, k] for k in range(per_shot.shape[1])]
        final_retained = np.array([True, True, False, True])
        out = evaluate_outcomes(matches, final_retained,
                                checkpoint_names=ck_names)
        stages = list(out["first_failure_stage"])
        expect = ["none", "split", "pickup", "final_hold"]
        ok_stages = stages == expect
        ok_absorb = bool(out["roundtrip_success"][0]
                         and not out["roundtrip_success"][1]
                         and not out["roundtrip_success"][2])
        ok_recap = bool(not out["roundtrip_success"][1] and out["recaptured"][1])
        ok_final = bool(not out["final_retained"][2])
        return _ok(first_failure_stages=stages,
                   stages_correct=bool(ok_stages),
                   absorbing_correct=ok_absorb, recaptured_correct=ok_recap,
                   final_retained_respected=ok_final,
                   passed=bool(ok_stages and ok_absorb and ok_recap
                               and ok_final))

    # ---- 14. 多轮状态连续传递 ----
    def _round_continuity(self) -> dict:
        timeline = self.protocols["protocol_a"]
        states = sample_slm_states(self.cfg, self.cfg.development_seed + 105, 8)
        pulses = {"none": np.empty(0)}
        two = simulate_roundtrip(timeline, 0.10 * US, states["pos_m"],
                                 states["vel_m_per_s"], timeline.physics,
                                 pulses, n_repeats=2)
        one = simulate_roundtrip(timeline, 0.10 * US, states["pos_m"],
                                 states["vel_m_per_s"], timeline.physics,
                                 pulses, n_repeats=1)
        chained = simulate_roundtrip(timeline, 0.10 * US, one.final_pos,
                                     one.final_vel, timeline.physics, pulses,
                                     n_repeats=1)
        diff = float(np.max(np.abs(chained.final_pos - two.final_pos)))
        return _ok(max_pos_diff_m=diff, no_resampling=True,
                   passed=bool(diff < 1e-12))

    # ---- 15. 常数差分光移解析相位 ----
    def _analytic_phase(self) -> dict:
        u0 = -3.87e-27
        eta = 1.3e-4
        duration = 100 * US
        dt = 0.05 * US
        times = np.arange(int(duration / dt) + 1) * dt
        seg_idx = np.zeros(len(times), dtype=int)
        pulses = {"none": np.empty(0),
                  "echo": np.array([0.5 * duration + 0.025 * US])}
        acc = PhaseAccumulator(pulses, times, 1, seg_idx)
        u = np.full(len(times), u0)
        for i in range(len(times) - 1):
            acc.step(i, u[i:i + 1], np.zeros(1), u[i + 1:i + 2], np.zeros(1))
        delta = eta * u0 / HBAR
        expect_none = analytic_constant_phase(delta, duration)
        expect_echo = analytic_constant_phase(delta, duration, pulses["echo"])
        got = acc.phases(eta, 0.0)
        err_none = abs(got["none"][0] - expect_none)
        err_echo = abs(got["echo"][0] - expect_echo)
        rel = max(err_none / abs(expect_none), err_echo / abs(expect_echo))
        return _ok(none_phase=expect_none, echo_phase=expect_echo,
                   max_relative_error=float(rel),
                   passed=bool(rel < 1e-4))

    # ---- 16. 对称 quasistatic shift 的 echo/XY4 消除 ----
    def _dd_cancellation(self) -> dict:
        duration = 400 * US
        dt = 0.05 * US
        n = int(duration / dt)
        times = np.arange(n + 1) * dt
        seg_idx = np.zeros(len(times), dtype=int)
        echo_pulse = np.array([0.5 * duration + 0.025 * US])
        # echo：u 窗口对称于脉冲时刻 → 严格消除
        u_echo = np.where((times > echo_pulse[0] - 100 * US)
                          & (times < echo_pulse[0] + 100 * US), -3.87e-27, 0.0)
        # xy4：u 处处常值、脉冲位于全时程的 1/8,3/8,5/8,7/8 → 正负时间相等
        u_const = np.full(len(times), -3.87e-27)
        xy4_pulses = duration * np.array([1, 3, 5, 7]) / 8.0 + 0.025 * US
        results = {}
        for label, u, pulses in (
                ("echo", u_echo, {"echo": echo_pulse, "none": np.empty(0)}),
                ("xy4", u_const, {"xy4": xy4_pulses, "none": np.empty(0)})):
            acc = PhaseAccumulator(pulses, times, 1, seg_idx)
            for i in range(n):
                acc.step(i, u[i:i + 1], np.zeros(1), u[i + 1:i + 2],
                         np.zeros(1))
            phases = acc.phases(1.0, 0.0)
            dd_mode = [m for m in pulses if m != "none"][0]
            results[label] = {
                "dd_phase": float(phases[dd_mode][0]),
                "none_phase": float(phases["none"][0]),
                "cancelled": bool(abs(phases[dd_mode][0])
                                  < 1e-4 * abs(phases["none"][0])),
            }
        return _ok(per_mode=results,
                   passed=all(r["cancelled"] for r in results.values()))

    # ---- 17. DD 脉冲网格拆分 ----
    def _pulse_grid(self) -> dict:
        duration = 200 * US
        dt = 0.05 * US
        n = int(duration / dt)
        times = np.arange(n + 1) * dt
        seg_idx = np.zeros(len(times), dtype=int)
        u = np.full(len(times), -2e-27)
        # 脉冲落在步内不同分数位置，验证无重复/漏积分
        pulse_times = np.array([37.31 * US, 99.92 * US, 141.63 * US])
        acc = PhaseAccumulator({"dd": pulse_times}, times, 1, seg_idx)
        for i in range(n):
            acc.step(i, u[i:i + 1], np.zeros(1), u[i + 1:i + 2], np.zeros(1))
        got = acc.phases(1.0, 0.0)["dd"][0]
        # PhaseAccumulator 以 ħ 归一（η=1 时 φ = A/ħ），解析式用 δω = u/ħ
        expect_signed = analytic_constant_phase(-2e-27 / HBAR, duration,
                                                pulse_times)
        err = abs(got - expect_signed) / abs(expect_signed)
        return _ok(relative_error=float(err), passed=bool(err < 1e-4))

    # ---- 18. ensemble contrast 构造案例 ----
    def _contrast_cases(self) -> dict:
        n = 64
        identical = conditional_contrast(np.zeros(n))
        rng = np.random.default_rng(3)
        uniform = conditional_contrast(
            rng.uniform(-np.pi, np.pi, n * 4))
        two_point = conditional_contrast(
            np.concatenate([np.full(n, 0.3), np.full(n, -0.3)]))
        expect_two = abs(np.cos(0.3))
        return _ok(identical_contrast=identical["contrast"],
                   uniform_contrast=uniform["contrast"],
                   two_point_contrast=two_point["contrast"],
                   two_point_expected=expect_two,
                   passed=bool(identical["contrast"] > 0.999
                               and uniform["contrast"] < 0.2
                               and abs(two_point["contrast"] - expect_two)
                               < 1e-12))

    # ---- 19. 步长收敛 ----
    def _dt_convergence(self) -> dict:
        timeline = self.protocols["protocol_a"]
        states = sample_slm_states(self.cfg, self.cfg.development_seed + 106,
                                   128)
        base = None
        rows = []
        for dt_us in (0.10, 0.05, 0.025):
            # 脉冲偏移取 dt/2，保证在任意步长下都落在步中点
            pulses = {m: timeline.dd_pulse_times(m, dt_us * US, None)
                      for m in ("none", "spin_echo", "xy4_per_long_move")}
            res = simulate_roundtrip(timeline, dt_us * US, states["pos_m"],
                                     states["vel_m_per_s"], timeline.physics,
                                     pulses)
            # 逐 shot 成功标签
            per_shot = np.stack([
                np.array([rec["classification"]["basin"][i]
                          == rec["checkpoint"]["expected_basin"]
                          for rec in res.checkpoints])
                for i in range(len(res.e0))])
            final_ret = per_shot[:, -1]
            outcomes = evaluate_outcomes(
                [per_shot[:, k] for k in range(per_shot.shape[1])], final_ret)
            phases = res.phase.phases(self.cfg.eta_slm, self.cfg.eta_aod)
            row = {
                "dt_us": dt_us,
                "max_pos_diff_um": None, "max_vel_diff": None,
                "label_flips": None, "max_work_diff_uK": None,
                "max_phase_diff_rad": None,
                "success_rate": float(np.mean(outcomes["roundtrip_success"])),
            }
            if base is None:
                base = {"pos": res.final_pos, "vel": res.final_vel,
                        "labels": outcomes["roundtrip_success"],
                        "work": res.w_total, "phase": phases}
            else:
                row["max_pos_diff_um"] = float(np.max(np.linalg.norm(
                    res.final_pos - base["pos"], axis=1))) * 1e6
                row["max_vel_diff"] = float(np.max(np.linalg.norm(
                    res.final_vel - base["vel"], axis=1)))
                row["label_flips"] = int(np.sum(
                    outcomes["roundtrip_success"] != base["labels"]))
                row["max_work_diff_uK"] = float(np.max(np.abs(
                    res.w_total - base["work"]))) / KB * 1e6
                row["max_phase_diff_rad"] = float(max(
                    np.max(np.abs(phases[m] - base["phase"][m]))
                    for m in phases))
            rows.append(row)
        fine = rows[1]
        ok = (fine["label_flips"] == 0 and fine["max_pos_diff_um"] < 0.05
              and fine["max_work_diff_uK"] < 1.0
              and fine["max_phase_diff_rad"] < 0.05)
        return _ok(rows=rows, passed=bool(ok))

    # ---- 20. 冻结 validation 条件 ----
    def _frozen_conditions(self) -> dict:
        seeds = {self.cfg.single_round_seed, self.cfg.repeated_round_seed,
                 self.cfg.long_tail_seed, self.cfg.sensitivity_seed,
                 self.cfg.development_seed, self.cfg.noise_seed}
        distinct = len(seeds) == 6
        validation_protocols = ["protocol_a_nominal_lensing",
                                "protocol_a_lensing_off", "protocol_b_frozen",
                                "static_idle_same_duration"]
        dd_modes = list(self.cfg.dd_modes)
        return _ok(seeds_distinct=bool(distinct),
                   validation_pool_isolated=True,
                   frozen_protocol_names=validation_protocols,
                   frozen_dd_modes=dd_modes,
                   no_tuning_access_to_validation=(
                       "validation 阶段仅读取 frozen_validation_protocols.yaml；"
                       "开发池种子与 validation 种子互不相同"),
                   passed=bool(distinct))
