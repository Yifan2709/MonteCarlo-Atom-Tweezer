"""分段功-能账本：把 RoundtripResult 的 checkpoint 快照组装成 ledger 行。

显式时变阶段满足 ΔE = W_ext + R_W；静态 hold 段应只有有界数值能量误差。
"""
from __future__ import annotations

import numpy as np

from level3_3d_transport_lensing.gaussian_3d import KB

TARGET_WELL = {
    "initial_slm_hold": "slm", "pickup": "aod", "split": "aod",
    "outbound_long_move": "aod", "remote_operation_hold": "aod",
    "inbound_long_move": "aod", "merge": "aod", "dropoff": "slm",
    "final_slm_hold": "slm", "final_aod_hold": "aod", "static_idle": "slm",
}


def segment_ledger_rows(result, grid, protocol_name, shot_ids=None,
                        per_shot=True, uK=True):
    """按 (segment, shot) 生成功-能账本行。

    需要 result.checkpoints 的累积功快照与 grid 的分段边界；首段起点取
    e0 与零功。返回行列表；per_shot=False 时返回按 shot 聚合的均值行。
    """
    n_shots = len(result.e0)
    ids = np.arange(n_shots) if shot_ids is None else np.asarray(shot_ids)
    seg_names = grid["segment_names"]
    rows = []
    marks = result.segment_marks
    seg_works = {}
    for k in range(len(marks) - 1):
        m0, m1 = marks[k], marks[k + 1]
        seg_idx = m0["segment_index"]
        seg = seg_names[seg_idx]
        e_start, e_end = m0["e_total"], m1["e_total"]
        w_move = m1["cum_w_move"] - m0["cum_w_move"]
        w_shift = m1["cum_w_shift"] - m0["cum_w_shift"]
        w_aod = m1["cum_w_aod_depth"] - m0["cum_w_aod_depth"]
        w_slm = m1["cum_w_slm_depth"] - m0["cum_w_slm_depth"]
        w_ext = w_move + w_shift + w_aod + w_slm
        delta_e = e_end - e_start
        residual = delta_e - w_ext
        target = TARGET_WELL.get(seg, "aod")
        d_target = (grid["d_aod"] if target == "aod" else grid["d_slm"])
        bottom = d_target[m1["step"]]
        excitation = e_end - bottom
        seg_works[seg] = {"w_move": w_move, "w_shift": w_shift,
                          "w_aod_depth": w_aod, "w_slm_depth": w_slm,
                          "delta_e": delta_e, "residual": residual}
        if per_shot:
            for i in ids:
                rows.append(_ledger_row(
                    protocol_name, seg, i, e_start[i], e_end[i],
                    w_move[i], w_shift[i], w_aod[i], w_slm[i], w_ext[i],
                    delta_e[i], residual[i], excitation[i], uK))
    return rows, seg_works


def _ledger_row(protocol, segment, shot, e_start, e_end, w_move, w_shift,
                w_aod, w_slm, w_ext, delta_e, residual, excitation, uK=True):
    """单行账本（μK 或 J）。"""
    scale = 1.0 / KB * 1e6 if uK else 1.0
    unit = "_uK" if uK else "_J"
    row = {"protocol": protocol, "segment": segment, "shot_id": int(shot)}
    for name, value in (("e_start", e_start), ("e_end", e_end),
                        ("w_move", w_move), ("w_shift", w_shift),
                        ("w_aod_depth_work", w_aod),
                        ("w_slm_depth_work", w_slm), ("w_ext_total", w_ext),
                        ("delta_e", delta_e), ("residual", residual),
                        ("excitation_vs_target_bottom", excitation)):
        row[name + unit] = float(value * scale)
    return row


def segment_ledger_summary(seg_works, uK=True):
    """每段的跨 shot 聚合（均值与最大绝对残差）。"""
    rows = []
    scale = 1.0 / KB * 1e6 if uK else 1.0
    for seg, works in seg_works.items():
        if not isinstance(works.get("delta_e"), np.ndarray):
            continue
        rows.append({
            "segment": seg,
            "mean_delta_e_uK": float(np.mean(works["delta_e"]) * scale),
            "mean_w_ext_uK": float(np.mean(
                works["w_move"] + works["w_shift"] + works["w_aod_depth"]
                + works["w_slm_depth"]) * scale),
            "max_abs_residual_uK": float(np.max(np.abs(works["residual"])) * scale),
            "mean_abs_residual_uK": float(np.mean(np.abs(works["residual"])) * scale),
            "mean_w_move_uK": float(np.mean(works["w_move"]) * scale),
            "mean_w_shift_uK": float(np.mean(works["w_shift"]) * scale),
            "mean_w_aod_depth_uK": float(np.mean(works["w_aod_depth"]) * scale),
            "mean_w_slm_depth_uK": float(np.mean(works["w_slm_depth"]) * scale),
        })
    return rows


def static_hold_check(result, grid, hold_kinds=("hold", "final_hold",
                                                "remote_hold")):
    """静态 hold 段的功-能检查：外部功应≈0，残差即数值误差。"""
    seg_names = grid["segment_names"]
    seg_kinds = grid["segment_kinds"]
    checks = {}
    ck_seg = [rec["checkpoint"]["segment_index"] for rec in result.checkpoints]
    cum = ([np.zeros(len(result.e0))]
           + [rec["cum_work"] for rec in result.checkpoints])
    for k, seg_idx in enumerate(ck_seg):
        seg = seg_names[seg_idx]
        if seg_kinds[seg_idx] not in hold_kinds:
            continue
        w_ext = cum[k + 1] - cum[k]
        checks[seg] = {
            "max_abs_w_ext_over_depth": float(np.max(np.abs(w_ext))
                                              / max(abs(grid["d_slm"][0]),
                                                    abs(grid["d_aod"][0]),
                                                    1e-30)),
            "max_abs_w_ext_J": float(np.max(np.abs(w_ext))),
        }
    return checks
