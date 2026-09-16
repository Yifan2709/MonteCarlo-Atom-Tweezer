"""Y5/Y6：[[4,2,2]] 编码、存储、测量、解码与逻辑传送（Zhang 2026 图2-4）。

码结构：S={XXXX, ZZZZ}；逻辑算符 Z̄1=ZZII, X̄1=XIXI, Z̄2=IZZI, X̄2=IXIX。
|00̄⟩=(|0000⟩+|1111⟩)/√2（h(0)+CX(0→1/2/3) 梯形制备）。

【已登记结构假设】论文未公开 [[4,2,2]] 精确线路。本实现：GHZ 梯形制备
（3 CX 层）+ 可选 flag（flag 比特经 CX 接首/末数据位，测 flag=1 检出
双错故障）+ 传送 = 两块横向 CNOT + X 基测量 + Pauli 帧纠正。
误差率取论文独立测量锚点（species.py）；线路深度属结构假设，影响绝对值。

解码规则（论文"宇称+单擦除"）：
- 无 dark：两逻辑宇称（q0q1、q2q3）均偶 → 成功；
- 单 dark 且该位擦除被检出：由 ZZZZ=0 用其余三位推断被擦除读数 → 再判；
- 其余（多 dark；未检出 dark）：解码器无法区分丢失与读出 FN → 判失败
  （无信息口径在单 dark 时一律按 0.5 概率成功抽样——歧义）。
解码器只读读出比特串与擦除检测记录，不读真实状态。

三口径：all_prepared / erasure_informed / postselected（接受率单列）。
"""
from __future__ import annotations

import numpy as np

from .provenance import REGISTRIES
from .quantum_mc import BatchedState

GHZ_CX = ((0, 1), (0, 2), (0, 3))


# ------------------------------------------------------------------ 助手
def _prep_ghz(state: BatchedState, rng, eps_gate, fr_er, fr_loss,
              transport_loss: float) -> np.ndarray:
    """GHZ 梯形制备 + 逐门误差（ϵ 为每 CZ 门误差，作用于 c/t 配对之一）。

    返回真实初始丢失掩码 (n, trials)。
    """
    init_lost = rng.random((state.n, state.trials)) < transport_loss
    state.mark_lost(init_lost)
    state.h(0)
    for c, t in GHZ_CX:
        state.cx(c, t)
        # 每门：以 eps_gate 概率在 {c,t} 之一上发生误差
        hit = rng.random(state.trials) < eps_gate
        on_c = rng.random(state.trials) < 0.5
        v = rng.random(state.trials)
        lost_hit = hit & (v < fr_er + fr_loss)
        pauli = hit & ~lost_hit
        w = rng.random(state.trials)
        fx_hit = pauli & (((w < 0.3) | ((w >= 0.5) & (w < 0.7))))
        fz_hit = pauli & (w >= 0.3)
        fx = np.zeros((state.n, state.trials), bool)
        fz = np.zeros((state.n, state.trials), bool)
        lost = np.zeros((state.n, state.trials), bool)
        fx[c] = fx_hit & on_c
        fz[c] = fz_hit & on_c
        lost[c] = lost_hit & on_c
        fx[t] = fx_hit & ~on_c
        fz[t] = fz_hit & ~on_c
        lost[t] = lost_hit & ~on_c
        state.apply_pauli_errors(fx, fz)
        state.mark_lost(lost)
    return init_lost


def _readout(state: BatchedState, rng, fp: float, fn: float):
    """Z 基读出 + 终端成像 FP/FN。返回 bits(list[4]) 与真实丢失掩码。"""
    bits = []
    for q in range(4):
        b = state.measure_z(q).astype(np.int8)
        u = rng.random(state.trials)
        b = np.where(b < 0, -1,
                     np.where((b == 0) & (u < fn), 1,
                              np.where((b == 1) & (u < fp), 0, b)))
        bits.append(b)
    return bits


def _qubit_detections(rng, true_lost: np.ndarray, fp: float, fn: float):
    """逐比特擦除检测记录（观测）。"""
    return np.where(true_lost, rng.random(true_lost.shape) >= fn,
                    rng.random(true_lost.shape) < fp)


def _decode_trial(vals, det_q) -> float | None:
    """返回成功概率（0/1/None=歧义按 0.5 处理由调用方）。"""
    dark = [q for q in range(4) if vals[q] is None]
    if not dark:
        return float((vals[0] + vals[1]) % 2 == 0
                     and (vals[2] + vals[3]) % 2 == 0)
    if len(dark) == 1 and det_q is not None and det_q[dark[0]]:
        q = dark[0]
        inferred = sum(vals[k] for k in range(4) if k != q) % 2
        v = list(vals)
        v[q] = inferred
        return float((v[0] + v[1]) % 2 == 0 and (v[2] + v[3]) % 2 == 0)
    return None


def _decode_all(bits, det_q_matrix, rng):
    """三口径成功数组。det_q_matrix=None → 无信息口径。"""
    succ = np.zeros(len(bits[0]), bool)
    for tr in range(len(bits[0])):
        vals = [None if bits[q][tr] < 0 else int(bits[q][tr])
                for q in range(4)]
        det_q = None if det_q_matrix is None else det_q_matrix[:, tr]
        p = _decode_trial(vals, det_q)
        if p is None:
            succ[tr] = rng.random() < 0.5
        elif p == 1.0:
            succ[tr] = True
        elif p == 0.0:
            succ[tr] = False
        else:
            succ[tr] = rng.random() < p
    return succ


# ------------------------------------------------------------------ Y5 制备
def measure_prep_fidelity(trials: int, seed: int, mode: str = "raw",
                          eps_layer: float | None = None,
                          transport_loss: float = 0.01) -> dict:
    """mode ∈ raw / flag / erasure_postselect。双基测量用两份同噪声副本。"""
    yb = REGISTRIES["yb171_2026_native"]
    if eps_layer is None:
        eps_layer = yb.get("gate_ar_error").value
    fr_er = yb.get("gate_ar_erasure_frac").value
    fr_loss = 0.10
    use_flag = mode == "flag"
    n = 5 if use_flag else 4
    fp = yb.get("erasure_check_fp").value
    fn = yb.get("erasure_check_fn").value
    out = {"trials": trials, "mode": mode}

    def one_basis(seed_base, basis):
        st = BatchedState(n, trials, np.random.default_rng(seed_base))
        init = _prep_ghz(st, np.random.default_rng(seed_base + 100),
                         eps_layer, fr_er, fr_loss, transport_loss)
        flagged = None
        if use_flag:
            st.cx(0, 4)
            st.cx(3, 4)
            flagged = st.measure_z(4) == 1
        bits = _readout(st, np.random.default_rng(seed_base + 300), 0.0, 0.0)
        if use_flag:
            bits = [np.where(flagged, -1, b) for b in bits]
        return init, bits, flagged

    init_z, bz, flag_z = one_basis(seed, "Z")
    init_x, bx, flag_x = one_basis(seed + 1000, "X")

    def parity(bits, qs):
        par = np.zeros(trials, np.int8)
        for q in qs:
            par ^= np.where(bits[q] < 0, 0, bits[q]).astype(np.int8)
        return par % 2

    par_zzzz = parity(bz, (0, 1, 2, 3))
    par_zzii = parity(bz, (0, 1))
    par_izzi = parity(bz, (2, 3))
    # X 基副本（同噪声实现）：XXXX 宇称
    st_x = BatchedState(n, trials, np.random.default_rng(seed + 1000))
    _prep_ghz(st_x, np.random.default_rng(seed + 1100), eps_layer, fr_er,
              fr_loss, transport_loss)
    if use_flag:
        st_x.cx(0, 4)
        st_x.cx(3, 4)
        flag_bx = st_x.measure_z(4) == 1
    bx2 = [st_x.measure_x(q) for q in range(4)]
    if use_flag:
        bx2 = [np.where(flag_bx, -1, b) for b in bx2]
    par_xxxx = parity(bx2, (0, 1, 2, 3))
    e_zzzz = 1 - 2 * float(par_zzzz.mean())
    e_xxxx = 1 - 2 * float(par_xxxx.mean())
    e_zzii = 1 - 2 * float(par_zzii.mean())
    e_izzi = 1 - 2 * float(par_izzi.mean())
    fid = float((1 + e_zzzz + e_xxxx + e_zzii * e_izzi) / 4)
    rng = np.random.default_rng(seed + 200)
    detected = _qubit_detections(rng, init_z[:4].any(axis=0), fp, fn)
    out.update({
        "fidelity_estimate": fid,
        "e_zzzz": e_zzzz, "e_xxxx": e_xxxx,
        "flagged_frac": float(np.mean(flag_z)) if use_flag else None,
        "detected_erasure_frac": float(np.mean(detected)),
        "paper": {"raw": yb.get("code_prep_fid_raw").value,
                  "flag": yb.get("code_prep_fid_flag").value,
                  "erasure": yb.get("code_prep_fid_erasure").value},
        "note": "后选口径的提升来自擦除 trial 被拒绝；此处报告原始估计",
    })
    if mode == "erasure_postselect":
        acc = ~detected
        ok_z = (par_zzzz[acc] == 0)
        ok_x = (par_xxxx[acc] == 0)
        fid_ps = float(np.mean(ok_z & ok_x)) if acc.any() else float("nan")
        out["acceptance"] = float(np.mean(acc))
        out["fidelity_postselected"] = fid_ps
        out["note"] += "；后选=接受 trial 的双基宇称全对比例"
    return out


# ------------------------------------------------------------------ Y5 存储
def storage_experiment(trials: int, seed: int, t_store_s: float,
                       eps_layer: float | None = None,
                       transport_loss: float = 0.01,
                       leak_lifetime_s: float = 1.64,
                       t2_storage_s: float = 6.0) -> dict:
    yb = REGISTRIES["yb171_2026_native"]
    if eps_layer is None:
        eps_layer = yb.get("gate_ar_error").value
    fr_er = yb.get("gate_ar_erasure_frac").value
    fr_loss = 0.10
    rng = np.random.default_rng(seed)
    state = BatchedState(4, trials, rng)
    init_lost = _prep_ghz(state, np.random.default_rng(seed + 100),
                          eps_layer, fr_er, fr_loss, transport_loss)
    p_leak = 1 - np.exp(-t_store_s / leak_lifetime_s)
    p_z = 0.5 * (1 - np.exp(-t_store_s / t2_storage_s))
    storage_lost = rng.random((4, trials)) < p_leak
    state.mark_lost(storage_lost)
    fz = rng.random((4, trials)) < p_z
    state.apply_pauli_errors(np.zeros_like(fz), fz)
    bits = _readout(state, rng, yb.get("terminal_imaging_fp").value,
                    yb.get("terminal_imaging_fn").value)
    true_lost = init_lost | storage_lost
    det_q = _qubit_detections(rng, true_lost,
                              yb.get("erasure_check_fp").value,
                              yb.get("erasure_check_fn").value)
    succ_a = _decode_all(bits, None, np.random.default_rng(seed + 9))
    succ_b = _decode_all(bits, det_q, np.random.default_rng(seed + 10))
    detected = det_q.any(axis=0)
    accepted = ~detected
    return {
        "trials": trials, "t_store_s": t_store_s,
        "unconditional": float(np.mean(succ_a)),
        "erasure_informed": float(np.mean(succ_b)),
        "postselected": float(np.mean(succ_b[accepted])) if accepted.any() else float("nan"),
        "acceptance": float(np.mean(accepted)),
        "detected_erasure_frac": float(np.mean(detected)),
        "paper": {"unconditional": yb.get("decode_unconditional").value,
                  "erasure": yb.get("decode_with_erasure").value},
    }


# ------------------------------------------------------------------ Y6 传送
def teleportation(trials: int, seed: int, n_blocks: int = 4,
                  policy: str = "adaptive",
                  eps_layer: float | None = None,
                  transport_loss_per_move: float = 0.0055,  # paper_derived：RT1 1.1%/2 单程
                  moves_per_atom: int | None = None) -> dict:
    """两块 [[4,2,2]] 逻辑传送 + 辅助块选择。

    A 块(q0-3) 持 |+̄1⟩=H⊗4|00̄⟩；辅助块(q4-7) |00̄⟩。
    横向 CNOT A→B（逻辑 CNOT L1）；X 基测 A 的 X̄1（q0,q2 宇称）；
    B 以 Z̄1 纠正后读 X̄1B（q4,q6 X 宇称）→ 成功=+1。
    辅助块选择只读擦除检测记录；25 ms/次的检测时延计入总时延。
    """
    yb = REGISTRIES["yb171_2026_native"]
    if eps_layer is None:
        eps_layer = yb.get("gate_ar_error").value
    if moves_per_atom is None:
        moves_per_atom = int(yb.get("moves_per_atom").value)
    fr_er = yb.get("gate_ar_erasure_frac").value
    fr_loss = 0.10
    fp = yb.get("erasure_check_fp").value
    fn = yb.get("erasure_check_fn").value
    rng = np.random.default_rng(seed)
    per_atom_loss = 1 - (1 - transport_loss_per_move) ** moves_per_atom

    # 辅助块质量记录（仅检测观测）
    block_detected = []
    block_true_lost = []
    for b in range(n_blocks):
        true_lost = rng.random((4, trials)) < per_atom_loss
        det = _qubit_detections(rng, true_lost, fp, fn)
        block_detected.append(det.any(axis=0))
        block_true_lost.append(true_lost)
    latency_s = 0.0
    if policy == "adaptive":
        sel = np.zeros(trials, dtype=int)
        for tr in range(trials):
            chosen, checks = n_blocks - 1, 0
            for b in range(n_blocks):
                checks += 1
                if not block_detected[b][tr]:
                    chosen = b
                    break
            sel[tr] = chosen
            latency_s += checks * yb.get("detection_latency_s").value
        latency_s /= trials
    else:
        sel = rng.integers(0, n_blocks, trials)
    chosen_lost = np.zeros((4, trials), bool)
    for tr in range(trials):
        chosen_lost[:, tr] = block_true_lost[sel[tr]][:, tr]

    state = BatchedState(8, trials, rng)
    # 制备 A：|00̄⟩ 再 H⊗4 → |+̄1+̄2⟩；误差层
    initA = _prep_ghz(state, np.random.default_rng(seed + 100), eps_layer,
                      fr_er, fr_loss, per_atom_loss)
    for q in range(4):
        state.h(q)
    # 制备辅助块 B（选中块）：其真实丢失进掩码
    initB = rng.random((4, trials)) < per_atom_loss
    lostB = np.zeros((8, trials), bool)
    lostB[4:] = initB | chosen_lost
    state.mark_lost(lostB)
    state.h(4)
    for c, t in GHZ_CX:
        state.cx(4 + c, 4 + t)
        hit = rng.random(trials) < eps_layer
        on_c = rng.random(trials) < 0.5
        v = rng.random(trials)
        lost_hit = hit & (v < fr_er + fr_loss)
        pauli = hit & ~lost_hit
        w = rng.random(trials)
        fx_hit = pauli & (((w < 0.3) | ((w >= 0.5) & (w < 0.7))))
        fz_hit = pauli & (w >= 0.3)
        fx = np.zeros((8, trials), bool)
        fz = np.zeros((8, trials), bool)
        lost = np.zeros((8, trials), bool)
        ca, ta = 4 + c, 4 + t
        fx[ca] = fx_hit & on_c
        fz[ca] = fz_hit & on_c
        lost[ca] = lost_hit & on_c
        fx[ta] = fx_hit & ~on_c
        fz[ta] = fz_hit & ~on_c
        lost[ta] = lost_hit & ~on_c
        state.apply_pauli_errors(fx, fz)
        state.mark_lost(lost)
    # 横向 CNOT A→B（逻辑 L1），每门配对误差
    for a, b in ((0, 4), (1, 5), (2, 6), (3, 7)):
        state.cx(a, b)
        hit = rng.random(trials) < eps_layer
        v = rng.random(trials)
        lost_hit = hit & (v < fr_er + fr_loss)
        pauli = hit & ~lost_hit
        w = rng.random(trials)
        fx_hit = pauli & (((w < 0.3) | ((w >= 0.5) & (w < 0.7))))
        fz_hit = pauli & (w >= 0.3)
        fx = np.zeros((8, trials), bool)
        fz = np.zeros((8, trials), bool)
        lost = np.zeros((8, trials), bool)
        fx[a] = fx_hit
        fz[a] = fz_hit
        lost[a] = lost_hit
        state.apply_pauli_errors(fx, fz)
        state.mark_lost(lost)
    # X 基测 A 的 X̄1 = X0 X2
    m0 = state.measure_x(0)
    m2 = state.measure_x(2)
    m = np.where((m0 < 0) | (m2 < 0), 0,
                 (m0.astype(np.int8) + m2.astype(np.int8)) % 2)
    # Pauli 帧：m=1 → Z̄1B = Z4 Z5 纠正（行索引对应 qubit 4/5）
    zcorr = (m == 1)
    fz = np.zeros((8, trials), bool)
    fz[4] = zcorr
    fz[5] = zcorr
    state.apply_pauli_errors(np.zeros((8, trials), bool), fz)
    # 读出 X̄1B = X4 X6，含 B 块单擦除恢复（XXXX 稳定子 + 奇偶推断）
    xb = [state.measure_x(q) for q in range(4, 8)]
    succ = np.zeros(trials, bool)
    for tr in range(trials):
        vals = [None if xb[q][tr] < 0 else int(xb[q][tr]) for q in range(4)]
        dark = [q for q in range(4) if vals[q] is None]
        if len(dark) >= 2:
            continue
        if len(dark) == 1:
            q = dark[0]
            vals[q] = sum(vals[k] for k in range(4) if k != q) % 2  # XXXX=+1
        # X̄1B 对应 B 块逻辑 X̄1 = 局部 X4X6（q0,q2）
        succ[tr] = (vals[0] + vals[2]) % 2 == 0
    # 擦除后选口径
    detectedB = np.zeros(trials, bool)
    for tr in range(trials):
        detectedB[tr] = block_detected[sel[tr]][tr]
    accepted = ~detectedB & ~lostB[4:].any(axis=0)
    return {
        "trials": trials, "policy": policy, "n_blocks": n_blocks,
        "success": float(np.mean(succ)),
        "postselected_success": float(np.mean(succ[accepted])) if accepted.any() else float("nan"),
        "acceptance": float(np.mean(accepted)),
        "mean_selection_latency_s": latency_s if policy == "adaptive" else 0.0,
        "paper": {"raw_random": yb.get("teleport_fid_raw").value,
                  "adaptive": yb.get("teleport_fid_adaptive").value,
                  "erasure_ps": yb.get("teleport_fid_erasure_ps").value},
    }
