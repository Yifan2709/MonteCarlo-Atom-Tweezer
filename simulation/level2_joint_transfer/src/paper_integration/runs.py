"""正式运行编排：preflight / scan / validation / integration / all。

每个 run 写 run_manifest.json（解析后配置、种子、耗时、版本、合同哈希）
与逐阶段 JSON 结果。B0 读取既有归档做回归核对，不重跑 Cs 物理也不改写
历史结论。
"""
from __future__ import annotations

import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import (contract_data, integration, rb_circuit, rb_entangled,
               rb_transport, species, yb_code, yb_gate, yb_leakage,
               yb_transport, yb_vls)

OUTPUT_ROOT = Path(__file__).resolve().parents[2] / "outputs" / "paper_integration"


def _manifest(stage: str, out_dir: Path, extra: dict) -> dict:
    contract = json.loads((contract_data.CONTRACT_DIR
                           / "acceptance_contract.json").read_text("utf-8"))
    return {
        "stage": stage, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0], "platform": platform.platform(),
        "numpy": np.__version__,
        "contract_sha256": contract.get("contract_sha256"),
        "config": species.species_summary(),
        "output_dir": str(out_dir), **extra,
    }


def _dump(out_dir: Path, name: str, obj) -> None:
    def _default(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return str(o)
    (out_dir / name).write_text(
        json.dumps(obj, ensure_ascii=False, indent=1, default=_default),
        encoding="utf-8")


def stage_preflight(out_dir: Path) -> dict:
    t0 = time.time()
    results: dict = {"b0": {}, "m1": {}}
    # B0：读取既有归档核对初始化修复与历史偏差保留
    repo = Path(__file__).resolve().parents[2]   # level2_joint_transfer 根
    init_status = repo / ("outputs/continuous_level2c_345/"
                          "run_20260915_manual600_initialization/RESULT_STATUS.md")
    scan_json = repo / ("outputs/continuous_level2c_345/"
                        "run_20260915_manual600_initialization/"
                        "initfix_T25_matched_manual_600us_N4096_dt0.0125_initSite_Tref25.json")
    b0 = {"init_status_file_exists": init_status.exists(),
          "known_deviation_still_reported": None,
          "initial_unbound_actual_site": None,
          "s60": None}
    if scan_json.exists():
        d = json.loads(scan_json.read_text("utf-8"))
        b0["initial_unbound_actual_site"] = d.get("initial_unbound_actual_site")
        b0["s60"] = d.get("survival", [None])[-1]
        b0["known_deviation_still_reported"] = True
    results["b0"] = b0
    # M1：账本与键控随机流自检
    from .events import Ledger, TrueEvent, DetectionRecord, stream_seed
    led = Ledger()
    led.add_true(TrueEvent(0, 0, 1e-3, "leakage", "scatter"))
    led.add_detection(DetectionRecord(0, 0, 2e-3, "erasure_check", True))
    s1 = stream_seed(7, 0, 0, "ch")
    s2 = stream_seed(7, 0, 0, "ch")
    s3 = stream_seed(7, 1, 0, "ch")
    results["m1"] = {"ledger_summary": led.summary(),
                     "stream_stable": s1 == s2, "stream_keyed": s1 != s3}
    # 修正（验证 D8）：B0 增加本版本代码的初态束缚实检（不再只读历史归档）
    from continuous_transfer.initialization import sample_site_bound_harmonic
    KBc = 1.380649e-23
    cs_physics = {"slm_depth_j": 140e-6 * KBc, "slm_waist_m": 1.17e-6,
                  "slm_center_m": 0.0,
                  "mass_kg": 132.90545196 * 1.66053906660e-27}
    init = sample_site_bound_harmonic(90210, 25.0, cs_physics,
                                       np.ones(256))
    e = init["initial_energy_j"]
    results["b0"]["current_code_init_check"] = {
        "sampler": init["sampler"], "n": 256,
        "all_bound": bool(np.all(e < 0)),
        "mean_energy_uK": float(np.mean(e) / KBc / 1e-6)}
    ok = (b0["initial_unbound_actual_site"] == 0 and s1 == s2 and s1 != s3
          and b0["known_deviation_still_reported"]
          and results["b0"]["current_code_init_check"]["all_bound"])
    results["preflight_ok"] = bool(ok)
    _dump(out_dir, "preflight.json", results)
    _dump(out_dir, "run_manifest.json",
          _manifest("preflight", out_dir, {"ok": ok, "elapsed_s": time.time() - t0}))
    return results


def stage_scan_rb(out_dir: Path, shots: int, seed: int) -> dict:
    t0 = time.time()
    durations = [60e-6, 100e-6, 150e-6, 200e-6, 250e-6, 300e-6, 400e-6, 500e-6, 800e-6]
    rows = rb_transport.scan_transports(110e-6, durations, shots, seed, 0.05e-6)
    crit = rb_transport.evaluate_structure_criteria(rows)
    # R2：锚点标定 + 形状扫描
    cal = rb_entangled.calibrate_c_phi(110e-6, 300e-6, shots=shots, seed=seed)
    scan2 = rb_entangled.scan_entangled(
        110e-6, [100e-6, 200e-6, 300e-6, 600e-6], shots, seed + 1,
        cal["c_phi_rad_per_quanta"])
    crit2 = rb_entangled.evaluate_structure(scan2)
    dd_off = rb_entangled.entangled_transport(
        110e-6, 1000e-6, shots, seed + 2, cal["c_phi_rad_per_quanta"], dd=False)
    dd_on = [r for r in scan2 if r["time_s"] == 600e-6][0]
    dd_on_1ms = rb_entangled.entangled_transport(
        110e-6, 1000e-6, shots, seed + 3, cal["c_phi_rad_per_quanta"], dd=True)
    dd_on = dd_on_1ms  # DD 对照在静态失相干可分辨时长（σ=2T/T2*，T=1ms→0.5 rad）
    crit2["C3_dd_on_greater_than_off"] = bool(
        dd_on["fidelity_raw_loss_as_one"] > dd_off["fidelity_raw_loss_as_one"])
    out = {"r1_rows": rows, "r1_criteria": crit,
           "r2_calibration": cal, "r2_scan": scan2, "r2_criteria": crit2,
           "r2_dd_off": dd_off["fidelity_raw_loss_as_one"],
           "r2_dd_on": dd_on["fidelity_raw_loss_as_one"]}
    _dump(out_dir, "scan_rb2022.json", out)
    _dump(out_dir, "run_manifest.json",
          _manifest("scan_rb2022", out_dir,
                    {"shots": shots, "seed": seed,
                     "elapsed_s": time.time() - t0}))
    return out


def stage_scan_yb(out_dir: Path, shots: int, seed: int) -> dict:
    t0 = time.time()
    times = [0.4e-3, 0.5e-3, 0.6e-3, 0.7e-3, 0.89e-3, 1.5e-3, 3.0e-3]
    rows = yb_transport.scan_moves(
        100e-6, times, shots, seed, 0.1e-6,
        taus=[0.0, 8.0], trajs=["low_peak_v", "gamma:1.875"])
    crit = yb_transport.evaluate_structure_criteria(rows)
    # 绝热窗口复核（无透镜存活>0.98 的区间上判 C2/C3）
    out = {"y1_rows": rows, "y1_criteria": crit,
           "note": "τ_eff=8 为 assumed 上界扫描；绝对曲线不验收"}
    _dump(out_dir, "scan_yb2026.json", out)
    _dump(out_dir, "run_manifest.json",
          _manifest("scan_yb2026", out_dir,
                    {"shots": shots, "seed": seed,
                     "elapsed_s": time.time() - t0}))
    return out


def stage_validation(out_dir: Path, trials: int, seed: int,
                     which: str) -> dict:
    t0 = time.time()
    out: dict = {}
    if which in ("rb2022", "all"):
        cz = rb_circuit.single_cz_layer_fidelity(200000, seed)
        st = rb_circuit.steane_circuit(trials, seed + 1,
                                       transport_layer_loss=0.005)
        out["r3"] = {"cz_layer": cz, "steane": st,
                     "anchors": rb_circuit.evaluate_anchors(st)}
    if which in ("yb2026", "all"):
        out["y2"] = yb_vls.evaluate_anchors()
        out["y2_contrast"] = {c: yb_vls.ramsey_contrast_vs_trips(c, 30, 4000,
                                                                seed + 2)
                              for c in ("static_slm", "optimal", "parallel")}
        rt = yb_leakage.simulate_roundtrips(4096, 5, seed + 3)
        out["y3"] = {"result": rt, "anchors": yb_leakage.evaluate_anchors(rt),
                     "confusion": yb_leakage.confusion_matrix()}
        z = yb_gate.run_zone_mc(200000, seed + 4)
        sc = {"AR": yb_gate.scan_scaling("AR", n_trials=100000),
              "TO": yb_gate.scan_scaling("TO", n_trials=100000)}
        out["y4"] = {"zone": z, "scaling": sc,
                     "anchors": yb_gate.evaluate_anchors(z, sc)}
        prep = {m: yb_code.measure_prep_fidelity(trials, seed + 5, m)
                for m in ("raw", "flag", "erasure_postselect")}
        storage = [yb_code.storage_experiment(trials, seed + 6 + i, t)
                   for i, t in enumerate((0.0, 0.25, 0.5, 1.0))]
        out["y5"] = {"prep": prep, "storage": storage}
        tele_ad = yb_code.teleportation(trials, seed + 10, policy="adaptive")
        tele_rd = yb_code.teleportation(trials, seed + 11, policy="random")
        out["y6"] = {"adaptive": tele_ad, "random": tele_rd}
    _dump(out_dir, f"validation_{which}.json", out)
    _dump(out_dir, "run_manifest.json",
          _manifest(f"validation_{which}", out_dir,
                    {"trials": trials, "seed": seed,
                     "elapsed_s": time.time() - t0}))
    return out


def stage_integration(out_dir: Path) -> dict:
    t0 = time.time()
    yb = integration.yb_perturbation(move_time_s=0.4e-3, distance_m=100e-6,
                                     shots=192, seed=4242)
    rb = integration.rb_perturbation(shots=192, seed=777)
    out = {"yb_chain": yb, "rb_chain": rb,
           "all_propagated": bool(yb["all_propagated"] and rb["all_propagated"])}
    _dump(out_dir, "integration.json", out)
    _dump(out_dir, "run_manifest.json",
          _manifest("integration", out_dir,
                    {"elapsed_s": time.time() - t0}))
    return out
