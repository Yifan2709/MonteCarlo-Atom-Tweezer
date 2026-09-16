"""Y3：亚稳态泄漏、光散射/光电离、机械失阱与擦除测量（Zhang 2026 图2/ED3）。

通道（每往返 RT = 2 单程）：
- 稳态泄漏：2×交接 0.1% + 移动 0.1% = 0.3%/RT（论文散射锚点），其中
  r_e,l=0.72 可被擦除检查检出；
- 首往返瞬态（paper_derived 锚点：RT1 总损失 1.1(4)%、约一半可擦除）：
  额外机械 0.55% + 额外泄漏 0.55%（瞬态近阈值原子/交接激发，一次性）；
- 加热驱动机械损失：能量线性累积 + 近阈值初始份额（0.55%，assumed）。

检测：擦除检查 FP/FN=1.4%/1.4%；终端成像 FP=0.1%/FN=0.5%。
真实事件进 Ledger.true_events；观测进 detections（策略只可读后者）。
擦除占比口径 = 检出且对应真实泄漏丢失 / 该 RT 总丢失（FP 单独报告）。

论文锚点（结构判据）：RT1 损失 1.1(4)%（≈一半可擦除）；RT5 擦除占比
降至 ~10%（机械份额随加热上升）；混淆矩阵与 FP/FN 一致。
"""
from __future__ import annotations

import numpy as np

from .events import (DetectionRecord, Ledger, TrueEvent, channel_rng,
                     erasure_check_detection)
from .provenance import REGISTRIES

# 瞬态通道（paper_derived：由 RT1 总量 1.1% 与"约一半可擦除"拆分）
RT1_EXTRA_MECH = 0.0034
RT1_EXTRA_LEAK = 0.0046
NEAR_THRESHOLD_FRAC = 0.002   # 初始近阈值份额（assumed）


DIFFUSION_SIGMA_UK_PER_RT = 14.0   # assumed：能量随机行走展宽（扫描项）


def simulate_roundtrips(n_trials: int, n_roundtrips: int, seed: int,
                        heating_per_trip_uK: float = 4.0,
                        leak_per_rt: float | None = None,
                        storage_depth_uK: float | None = None,
                        erasure_checks: tuple[int, ...] = (1, 2, 3, 4, 5),
                        rt1_extra_mech_scale: float = 1.0,
                        ) -> dict:
    yb = REGISTRIES["yb171_2026_native"]
    if leak_per_rt is None:
        leak_per_rt = (yb.get("scatter_prob_per_handoff").value * 2
                       + yb.get("scatter_prob_per_move").value)
    if storage_depth_uK is None:
        storage_depth_uK = yb.get("storage_depth_uK").value
    fp = yb.get("erasure_check_fp").value
    fn = yb.get("erasure_check_fn").value
    re_frac = yb.get("leakage_erasure_fraction").value
    ledger = Ledger()
    rng = np.random.default_rng(seed)
    # 逐 trial 能量（μK）：主体 gamma 分布 + 显式近阈值份额
    energy = rng.gamma(6.0, 2.5, n_trials)
    near_thr = rng.random(n_trials) < NEAR_THRESHOLD_FRAC
    energy = np.where(near_thr, storage_depth_uK * (0.96 + 0.03 * rng.random(n_trials)), energy)
    alive = np.ones(n_trials, bool)
    leak_state = np.zeros(n_trials, bool)      # 已泄漏（待判丢失）
    records = []
    for rt in range(1, n_roundtrips + 1):
        n_loss_leak = 0
        n_loss_mech = 0
        # 1) 稳态泄漏
        new_leak = alive & ~leak_state & (rng.random(n_trials) < leak_per_rt)
        if rt == 1:
            new_leak |= alive & ~leak_state & (rng.random(n_trials) < RT1_EXTRA_LEAK)
        leak_state |= new_leak
        for tr in np.nonzero(new_leak)[0]:
            ledger.add_true(TrueEvent(tr, tr, rt * 2.5e-3, "leakage", "scatter"))
        # 2) 机械：加热推进（漂移+扩散，assumed σ=14 μK/RT）+ 阈值判定
        step = (2 * heating_per_trip_uK
                + DIFFUSION_SIGMA_UK_PER_RT * rng.standard_normal(n_trials))
        energy = energy + step * (alive & ~leak_state)
        mech = alive & ~leak_state & (energy >= storage_depth_uK)
        if rt == 1:
            mech |= alive & ~leak_state & (
                rng.random(n_trials) < RT1_EXTRA_MECH * rt1_extra_mech_scale)
        lost = alive & (leak_state | mech)
        n_loss_leak = int(np.sum(alive & leak_state))
        n_loss_mech = int(np.sum(mech))
        for tr in np.nonzero(lost)[0]:
            ledger.add_true(TrueEvent(
                tr, tr, rt * 2.5e-3,
                "mechanical_loss" if mech[tr] else "leakage", "transport"))
        alive = alive & ~lost
        leak_state = np.zeros(n_trials, bool)
        # 3) 擦除检查（声明的 RT 上；真实泄漏→检出率 1−FN；否则 FP）
        n_tp = n_fp = 0
        lost_this_rt = set(np.nonzero(lost)[0].tolist())
        if rt in erasure_checks:
            for tr in range(n_trials):
                # 检查只针对当前 RT 新丢失（was_leak）与仍存活原子（FP）
                was_leak = tr in lost_this_rt and _loss_kind(ledger, tr) == "leakage"
                rec = erasure_check_detection(
                    tr, tr, rt * 2.5e-3, was_leak, fp, fn,
                    channel_rng(seed, tr, tr, f"erase_check_rt{rt}"))
                ledger.add_detection(rec)
                if rec.detected:
                    if was_leak:
                        n_tp += 1
                    else:
                        n_fp += 1
        records.append({"roundtrip": rt, "alive": int(alive.sum()),
                         "loss_leak": n_loss_leak, "loss_mech": n_loss_mech,
                         "erasure_tp": n_tp, "erasure_fp": n_fp})
    frac = []
    for r in records:
        tot = r["loss_leak"] + r["loss_mech"]
        # 可擦除检出 = 泄漏中 r_e,l 份额的 TP 检出（口径：TP/总丢失）
        frac.append(r["erasure_tp"] / tot if tot > 0 else float("nan"))
    return {
        "n_trials": n_trials, "n_roundtrips": n_roundtrips,
        "per_roundtrip": records,
        "first_roundtrip_loss_frac": (records[0]["loss_leak"]
                                      + records[0]["loss_mech"]) / n_trials,
        "erasure_fraction_by_rt": frac,
        "ledger_summary": ledger.summary(),
        "params": {"leak_per_rt": leak_per_rt, "rt1_extra_mech": RT1_EXTRA_MECH,
                   "rt1_extra_leak": RT1_EXTRA_LEAK,
                   "heating_per_trip_uK": heating_per_trip_uK,
                   "near_threshold_frac": NEAR_THRESHOLD_FRAC},
    }


def _loss_kind(ledger: Ledger, trial_id: int) -> str:
    kinds = [e.kind for e in ledger.true_events if e.trial_id == trial_id
             and e.kind in ("leakage", "mechanical_loss")]
    return kinds[-1] if kinds else ""


def confusion_matrix(n_checks: int = 20000, seed: int = 5) -> dict:
    yb = REGISTRIES["yb171_2026_native"]
    fp, fn = yb.get("erasure_check_fp").value, yb.get("erasure_check_fn").value
    rng = np.random.default_rng(seed)
    has_leak = rng.random(n_checks) < 0.02
    det = np.where(has_leak, rng.random(n_checks) >= fn,
                   rng.random(n_checks) < fp)
    return {
        "P_detect_given_leak": float(det[has_leak].mean()),
        "P_detect_given_clean": float(det[~has_leak].mean()),
        "paper_fp": fp, "paper_fn": fn,
        "latency_s": yb.get("detection_latency_s").value,
    }


def evaluate_anchors(result: dict) -> dict:
    yb = REGISTRIES["yb171_2026_native"]
    f = result["erasure_fraction_by_rt"]
    return {
        "rt1_loss": result["first_roundtrip_loss_frac"],
        "rt1_loss_paper": yb.get("roundtrip1_loss_prob").value,
        "rt1_erasure_frac": f[0],
        "rt1_erasure_frac_paper": yb.get("roundtrip1_erasure_frac").value,
        "rt5_erasure_frac": f[4] if len(f) >= 5 else float("nan"),
        "rt5_erasure_frac_paper": yb.get("roundtrip5_erasure_frac").value,
        "erasure_frac_declines": bool(
            np.isfinite(f[0]) and len(f) >= 5 and np.isfinite(f[4])
            and f[4] < f[0]),
        "rt1_transient_elevated": bool(
            result["per_roundtrip"][0]["loss_leak"]
            + result["per_roundtrip"][0]["loss_mech"]
            > result["per_roundtrip"][1]["loss_leak"]
            + result["per_roundtrip"][1]["loss_mech"]),
    }
