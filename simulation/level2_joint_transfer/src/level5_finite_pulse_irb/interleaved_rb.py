"""Transport-interleaved RB：固定总门数 N，前 M 个门间各插入一次移动通道。

论文 Methods 设计：N=80 个 Clifford 固定，M 次 move 插在前 M 个门之后，
剩余 N−M 个门补齐使总门数相同，最后施加总 inverse。interleaved
operation 的含义显式配置（one_way_transport / aod_roundtrip /
full_pickup_transport_dropoff_roundtrip / static_idle_matched_duration），
不同 channel 不混成一个“move”。reference（M=0）与 transported 使用相同
Clifford strings、初态与配对的轨迹/噪声 draw（同一 string 的前 M 个
move draw 在不同 M 之间完全一致）。保存 atomic survival、conditional
spin return、survival-weighted return 与 loss-detected 口径。
"""
from __future__ import annotations

import numpy as np

from .qubit_hamiltonian import return_probability_z
from .spin_channels import MoveChannelEngine, StaticGateEngine

OPERATION_KINDS = ("one_way_transport", "aod_roundtrip",
                   "full_pickup_transport_dropoff_roundtrip",
                   "static_idle_matched_duration")


def generate_irb_sequences(group, total_cliffords, interleaved_counts,
                           sequences_per_count, seed):
    """生成共享 Clifford strings 的 IRB 设计（所有 M 用同一组 strings）。"""
    rng = np.random.default_rng(int(seed))
    strings = [[int(c) for c in group.sample_uniform(rng, total_cliffords)]
               for _ in range(sequences_per_count)]
    seqs = []
    for m in interleaved_counts:
        for j in range(sequences_per_count):
            cliffords = strings[j]
            inv = group.inverse[cliffords[-1]]
            for c in reversed(cliffords[:-1]):
                inv = int(group.mult[group.inverse[c], inv])
            seqs.append({"M": int(m), "sequence_index": int(j),
                         "cliffords": cliffords, "inverse": int(inv),
                         "string_id": int(j), "seed": int(seed)})
    return seqs


def run_interleaved_rb(seqs, move_engine: MoveChannelEngine,
                       gate_engine: StaticGateEngine, move_pool,
                       static_pool, eta_slm, eta_aod, seed,
                       n_sites=1, rabi_scale=None, delta_site=None,
                       delta_drive=0.0, noise_qs_fn=None, amp_error=0.0,
                       n_sub_pulse=4, dt_free_us=0.25, dt_pulse_us=0.1,
                       site_eta_scale=None, input_state="z_plus",
                       move_survival_override=None):
    """执行 IRB：所有 M 批量推进（同 strings 前缀配对）。

    noise_qs_fn(slot) 返回 [batch] quasistatic 噪声（None 关闭）。
    move_survival_override: 可选 [n_traj] 数组替换 move_pool 标签
    （static_idle_matched_duration 口径用）。
    """
    n_seq = len(seqs)
    rng = np.random.default_rng(int(seed))
    n_cliff = len(seqs[0]["cliffords"])
    max_moves = max(s["M"] for s in seqs)
    n_strings = len({s["string_id"] for s in seqs})
    # 配对 draw：k_moves[string, slot, site]（不同 M 的同 string 共享前缀）
    k_moves = rng.integers(0, move_pool.n_traj,
                           size=(n_strings, max_moves, n_sites))
    k_static = rng.integers(0, static_pool.n_traj, size=n_strings)
    surv_pool = (move_pool.roundtrip_success
                 if move_survival_override is None
                 else np.asarray(move_survival_override, dtype=bool))
    eta_s_base = np.broadcast_to(np.asarray(eta_slm, dtype=float), (n_sites,))
    eta_a_base = np.broadcast_to(np.asarray(eta_aod, dtype=float), (n_sites,))
    if site_eta_scale is not None:
        eta_s_site = eta_s_base * np.asarray(site_eta_scale, float)
        eta_a_site = eta_a_base * np.asarray(site_eta_scale, float)
    else:
        eta_s_site = np.broadcast_to(eta_s_base, (n_sites,)).copy()
        eta_a_site = np.broadcast_to(eta_a_base, (n_sites,)).copy()
    rabi_site = (np.ones(n_sites) if rabi_scale is None
                 else np.asarray(rabi_scale, float))
    dsite = None if delta_site is None else np.asarray(delta_site, float)

    flat_psi = np.zeros((n_seq, n_sites, 2), dtype=complex)
    flat_psi[..., 0 if input_state == "z_plus" else 1] = 1.0
    flat_psi = flat_psi.reshape(-1, 2)
    flat_trace = np.repeat(k_static[[s["string_id"] for s in seqs]],
                           n_sites)
    flat_tstart = np.zeros(n_seq * n_sites)
    survival = np.ones(n_seq * n_sites)
    flat_rabi = np.tile(rabi_site, n_seq)
    flat_dsite = None if dsite is None else np.tile(dsite, n_seq)

    for slot in range(n_cliff):
        flat_cliff = np.repeat(
            np.array([s["cliffords"][slot] for s in seqs], dtype=int),
            n_sites)
        noise_qs = None if noise_qs_fn is None else noise_qs_fn(slot)
        gate_engine.apply(flat_psi, flat_cliff, flat_trace, flat_tstart,
                          np.tile(eta_s_site, n_seq),
                          np.tile(eta_a_site, n_seq),
                          delta_site=flat_dsite,
                          noise_qs=noise_qs, rabi_scale=flat_rabi,
                          amp_error=amp_error, n_sub_pulse=n_sub_pulse,
                          delta_drive=delta_drive)
        for i, s in enumerate(seqs):
            flat_tstart[i * n_sites:(i + 1) * n_sites] += \
                gate_engine.clifford_duration_s(s["cliffords"][slot])
        if slot < max_moves:
            rows = [i for i, s in enumerate(seqs) if s["M"] > slot]
            if rows:
                k_sel = np.stack([k_moves[seqs[i]["string_id"], slot]
                                  for i in rows])  # [n_rows, n_sites]
                surv_sel = surv_pool[k_sel].astype(float)  # [n_rows, n_sites]
                idx = (np.array(rows)[:, None]
                       + np.zeros(n_sites, dtype=int)[None, :]).ravel()
                k_flat = k_sel.ravel()
                psi_move = flat_psi[idx].copy()
                eta_s_flat = np.tile(eta_s_site, len(rows))
                eta_a_flat = np.tile(eta_a_site, len(rows))
                noise_move = (None if noise_qs_fn is None
                              else np.asarray(noise_qs_fn(slot),
                                              float)[idx])
                move_engine.propagate(
                    psi_move, k_flat, eta_s_flat, eta_a_flat,
                    delta_drive=delta_drive,
                    delta_site=None if dsite is None else np.tile(dsite,
                                                                 len(rows)),
                    noise_qs=noise_move, dt_free_us=dt_free_us,
                    dt_pulse_us=dt_pulse_us,
                    rabi_scale=None if rabi_scale is None
                    else np.tile(rabi_site, len(rows)))
                flat_psi[idx] = psi_move
                survival[idx] *= surv_sel.ravel()
    flat_cliff = np.repeat(
        np.array([s["inverse"] for s in seqs], dtype=int), n_sites)
    noise_qs = None if noise_qs_fn is None else noise_qs_fn(-1)
    gate_engine.apply(flat_psi, flat_cliff, flat_trace, flat_tstart,
                      np.tile(eta_s_site, n_seq), np.tile(eta_a_site, n_seq),
                      delta_site=flat_dsite,
                      noise_qs=noise_qs, rabi_scale=flat_rabi,
                      amp_error=amp_error, n_sub_pulse=n_sub_pulse,
                      delta_drive=delta_drive)
    probs = return_probability_z(flat_psi).reshape(n_seq, n_sites)
    surv = survival.reshape(n_seq, n_sites)
    return {"return_probabilities": probs,
            "survival": surv,
            "conditional_return": probs,
            "joint_return": probs * surv,
            "loss_detected_return": probs * surv,
            "k_static": k_static,
            "k_moves": k_moves,
            "m_values": sorted({s["M"] for s in seqs}),
            "n_moves_per_sequence": np.array([s["M"] for s in seqs]),
            "string_ids": np.array([s["string_id"] for s in seqs]),
            "sequence_ids": [(s["M"], s["sequence_index"]) for s in seqs],
            "clifford_lists": [s["cliffords"] + [s["inverse"]]
                               for s in seqs],
            "operation": None}
