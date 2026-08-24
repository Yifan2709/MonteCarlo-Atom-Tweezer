"""Reference RB：静态 idle 条件下的标准单比特随机基准。

对每个 sequence length n 均匀抽 n 个 Clifford + 精确 inverse，用有限
脉冲（bare/SCROFULOUS）在静态池轨迹回放上传播，测回到初态的概率。
逐序列保存（含 Clifford 列表与 seed），不同 M 的 IRB 与 reference 使用
相同的门引擎与可配对的噪声 draw。
"""
from __future__ import annotations

import numpy as np

from .qubit_hamiltonian import return_probability_z
from .spin_channels import StaticGateEngine


def generate_reference_sequences(group, lengths, sequences_per_length,
                                 seed):
    """生成 RB 序列（均匀 Clifford + 精确 inverse）。"""
    rng = np.random.default_rng(int(seed))
    seqs = []
    for n in lengths:
        for j in range(sequences_per_length):
            cliffords = [int(c) for c in group.sample_uniform(rng, n)]
            inv = group.inverse[cliffords[-1]]
            for c in reversed(cliffords[:-1]):
                inv = int(group.mult[group.inverse[c], inv])
            seqs.append({"length": int(n), "sequence_index": int(j),
                         "cliffords": cliffords, "inverse": int(inv),
                         "seed": int(seed)})
    return seqs


def run_reference_rb(seqs, gate_engine: StaticGateEngine, group,
                     static_pool, eta_slm, eta_aod, seed,
                     n_sites=1, rabi_scale=None, delta_site=None,
                     delta_drive=0.0, noise_qs=None, amp_error=0.0,
                     n_sub_pulse=4, input_state="z_plus",
                     site_eta_scale=None):
    """执行 reference RB：批量传播所有序列 × 站点。

    返回逐 (sequence, site) 的返回概率、trace draw 索引与序列时长。
    noise_qs 需为 [n_seq*n_sites] 或 None（quasistatic per shot）。
    """
    n_seq = len(seqs)
    n_traj = static_pool.n_traj
    rng = np.random.default_rng(int(seed))
    k_trace = rng.integers(0, n_traj, size=n_seq)
    if site_eta_scale is not None and np.ndim(site_eta_scale):
        eta_slm = np.tile(np.asarray(eta_slm, float) * site_eta_scale,
                          n_seq)
        eta_aod = np.tile(np.asarray(eta_aod, float) * site_eta_scale,
                          n_seq)
    else:
        eta_slm = float(eta_slm)
        eta_aod = float(eta_aod)
    # 批量布局：[seq, site]；无门槽位填 identity（空分解 → 无操作）
    durations = np.array([
        sum(gate_engine.clifford_duration_s(c)
            for c in seqs[i]["cliffords"] + [seqs[i]["inverse"]])
        for i in range(n_seq)])
    psi = np.zeros((n_seq, n_sites, 2), dtype=complex)
    psi[..., 0 if input_state == "z_plus" else 1] = 1.0
    max_gates = max(len(s["cliffords"]) + 1 for s in seqs)
    flat_psi = psi.reshape(-1, 2)
    flat_cliff = np.zeros(n_seq * n_sites, dtype=int)
    flat_trace = np.repeat(k_trace, n_sites)
    flat_tstart = np.zeros(n_seq * n_sites)
    flat_rabi = (None if rabi_scale is None
                 else np.tile(np.asarray(rabi_scale, float), n_seq))
    flat_dsite = (None if delta_site is None
                  else np.tile(np.asarray(delta_site, float), n_seq))
    flat_noise = (None if noise_qs is None
                  else np.asarray(noise_qs, float).reshape(-1))
    for g in range(max_gates):
        have_gate = False
        flat_cliff[:] = 0
        for i, s in enumerate(seqs):
            full = s["cliffords"] + [s["inverse"]]
            if g < len(full):
                flat_cliff[i * n_sites:(i + 1) * n_sites] = full[g]
                have_gate = True
        if not have_gate:
            break
        gate_engine.apply(flat_psi, flat_cliff, flat_trace, flat_tstart,
                          eta_slm, eta_aod, delta_site=flat_dsite,
                          noise_qs=flat_noise, rabi_scale=flat_rabi,
                          amp_error=amp_error, n_sub_pulse=n_sub_pulse,
                          delta_drive=delta_drive)
        for i, s in enumerate(seqs):
            full = s["cliffords"] + [s["inverse"]]
            if g < len(full):
                flat_tstart[i * n_sites:(i + 1) * n_sites] += \
                    gate_engine.clifford_duration_s(full[g])
    probs = return_probability_z(flat_psi).reshape(n_seq, n_sites)
    return {"return_probabilities": probs, "k_trace": k_trace,
            "durations_s": durations,
            "sequence_ids": [(s["length"], s["sequence_index"])
                             for s in seqs],
            "clifford_lists": [s["cliffords"] + [s["inverse"]]
                               for s in seqs],
            "n_gates_per_sequence": np.array(
                [len(s["cliffords"]) + 1 for s in seqs])}
