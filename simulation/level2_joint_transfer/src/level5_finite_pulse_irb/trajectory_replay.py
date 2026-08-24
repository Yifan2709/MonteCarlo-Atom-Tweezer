"""Level 4 经典轨迹的冻结池、回放插值与 joint（联合经典-自旋）模式。

两种耦合模式：
A. frozen replay——把 Level 4 引擎生成的经典轨迹的 u_SLM/u_AOD 时序
   冻结为 NPZ 池（含 survival 标签与 Level 4 相位积分），量子动力学在
   插值轨迹上重复回放；默认用于大规模 RB。
B. joint propagation——复刻 Level 4 velocity-Verlet 步进，在同一时间网格
   上同时推进经典位置与量子态（预检查/小规模验证），用于校验 replay
   的相位与量子结果一致性；超过阈值应停止主 RB。

回放在脉冲边界、DD 翻转点与协议分段边界显式插入时间点；Level 4 的
survival 标签保持不变。
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from level3_3d_transport_lensing.aod_lensing import \
    lensing_potential_and_force
from level3_3d_transport_lensing.gaussian_3d import static_potential
from level4_roundtrip_coherence.level4_cli import _outcome_rows, _pulses_for
from level4_roundtrip_coherence.roundtrip_simulation import simulate_roundtrip
from level4_roundtrip_coherence.sampling4 import sample_slm_states

US = 1e-6
HBAR = 1.054571817e-34


@dataclass
class TrajectoryPool:
    """一组冻结经典轨迹的差分光移时序与标签。"""

    name: str
    times_s: np.ndarray            # [n_t] 记录时刻（含 0 与 T）
    u_slm: np.ndarray              # [n_traj, n_t] J（沿轨迹瞬时 SLM 势能）
    u_aod: np.ndarray              # [n_traj, n_t] J
    roundtrip_success: np.ndarray  # [n_traj] bool（Level 4 标签，不改）
    final_retained: np.ndarray     # [n_traj] bool
    first_failure_stage: np.ndarray  # [n_traj] int 代码
    stage_names: list              # 代码 → 阶段名
    u_slm_final: np.ndarray        # [n_traj] 末端静止势（J）
    u_aod_final: np.ndarray        # [n_traj]
    level4_a_slm: np.ndarray       # [n_traj] Level 4 相位积分 A_SLM（DD none）
    level4_a_aod: np.ndarray       # [n_traj]
    segment_boundaries_s: np.ndarray  # [n_seg+1] 分段边界
    segment_names: list
    total_duration_s: float
    meta: dict = field(default_factory=dict)

    @property
    def n_traj(self) -> int:
        """轨迹条数。"""
        return int(self.u_slm.shape[0])

    @property
    def sample_interval_s(self) -> float:
        """记录采样间隔。"""
        return float(self.times_s[1] - self.times_s[0])

    def sample_index(self, t: float) -> np.ndarray:
        """时刻 → 左侧样本索引。"""
        return np.searchsorted(self.times_s, t, side="right") - 1

    def u_at(self, u: np.ndarray, k_idx: np.ndarray, t: np.ndarray):
        """按轨迹索引与时刻线性插值势能（k_idx [B], t 广播到 [B]。）。"""
        t = np.broadcast_to(np.asarray(t, dtype=float), np.shape(k_idx))
        pos = np.clip(t / self.sample_interval_s, 0.0,
                      self.times_s.size - 1.000001)
        i0 = np.floor(pos).astype(np.int64)
        frac = pos - i0
        rows = u[k_idx]                       # [B, n_t]
        v0 = np.take_along_axis(rows, i0[:, None], axis=-1)
        v1 = np.take_along_axis(rows, np.minimum(i0 + 1,
                                                 self.times_s.size - 1)
                                [:, None], axis=-1)
        return v0 + (v1 - v0) * frac[:, None]

    def phase_integrals(self, k_idx, eta_slm, eta_aod, weights_t):
        """给定样本权重（[B, n_t] 或 None=均匀梯形）累计 A_SLM/A_AOD。"""
        if weights_t is None:
            raise ValueError("必须显式给定权重矩阵")
        a_slm = np.einsum("bt,bt->b", self.u_slm[k_idx].astype(np.float64),
                          weights_t)
        a_aod = np.einsum("bt,bt->b", self.u_aod[k_idx].astype(np.float64),
                          weights_t)
        return (eta_slm * a_slm + eta_aod * a_aod) / HBAR

    def save(self, path: Path) -> None:
        """保存 NPZ（float32 压缩主数组）。"""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path, times_s=self.times_s, u_slm=self.u_slm.astype(np.float32),
            u_aod=self.u_aod.astype(np.float32),
            roundtrip_success=self.roundtrip_success,
            final_retained=self.final_retained,
            first_failure_stage=self.first_failure_stage,
            u_slm_final=self.u_slm_final, u_aod_final=self.u_aod_final,
            level4_a_slm=self.level4_a_slm, level4_a_aod=self.level4_a_aod,
            segment_boundaries_s=self.segment_boundaries_s,
            total_duration_s=self.total_duration_s,
            meta_json=json.dumps(
                {"name": self.name, "stage_names": self.stage_names,
                 "segment_names": self.segment_names, **self.meta}))

    @classmethod
    def load(cls, path: Path) -> "TrajectoryPool":
        """读取 NPZ。"""
        data = np.load(Path(path), allow_pickle=False)
        meta = json.loads(str(data["meta_json"]))
        name = meta.pop("name")
        stage_names = meta.pop("stage_names")
        segment_names = meta.pop("segment_names")
        return cls(name=name, times_s=data["times_s"],
                   u_slm=data["u_slm"].astype(np.float64),
                   u_aod=data["u_aod"].astype(np.float64),
                   roundtrip_success=data["roundtrip_success"].astype(bool),
                   final_retained=data["final_retained"].astype(bool),
                   first_failure_stage=data["first_failure_stage"],
                   u_slm_final=data["u_slm_final"], u_aod_final=data["u_aod_final"],
                   level4_a_slm=data["level4_a_slm"],
                   level4_a_aod=data["level4_a_aod"],
                   segment_boundaries_s=data["segment_boundaries_s"],
                   total_duration_s=float(data["total_duration_s"]),
                   segment_names=segment_names, stage_names=stage_names,
                   meta=meta)

    def sha256(self) -> str:
        """内容哈希（u 数组 + 时间轴 + 标签）。"""
        digest = hashlib.sha256()
        for arr in (self.times_s, self.u_slm, self.u_aod,
                    self.roundtrip_success, self.final_retained,
                    self.first_failure_stage):
            digest.update(np.ascontiguousarray(arr).tobytes())
        digest.update(np.float64(self.total_duration_s).tobytes())
        return digest.hexdigest()


def generate_trajectory_pool(name: str, timeline, cfg4, seed: int,
                             n_traj: int, dt_classical_us: float = 0.1,
                             record_stride_us: float = 0.5,
                             chunk_size: int = 250, log=print) -> TrajectoryPool:
    """用 Level 4 引擎生成冻结轨迹池（分块控制内存）。"""
    dt = dt_classical_us * US
    stride = max(1, int(round(record_stride_us / dt_classical_us)))
    n_steps = int(round(timeline.total_duration_s / dt))
    # 丢弃 t=0 样本后的记录数：步号 {1, 1+stride, ...} ≤ n_steps
    n_t_expected = (n_steps - 1) // stride + 1
    all_u_slm = np.empty((n_traj, n_t_expected), dtype=np.float32)
    all_u_aod = np.empty((n_traj, n_t_expected), dtype=np.float32)
    success = np.zeros(n_traj, dtype=bool)
    retained = np.zeros(n_traj, dtype=bool)
    first_stage = np.zeros(n_traj, dtype=np.int16)
    a_slm = np.zeros(n_traj)
    a_aod = np.zeros(n_traj)
    stage_names = ["none"]
    states = sample_slm_states(cfg4, int(seed), n_traj)
    pulses = _pulses_for(timeline, dt, modes=("none",))
    done = 0
    while done < n_traj:
        n_chunk = min(chunk_size, n_traj - done)
        sl = slice(done, done + n_chunk)
        res = simulate_roundtrip(
            timeline, dt, states["pos_m"][sl], states["vel_m_per_s"][sl],
            timeline.physics, pulses, keep_ids=list(range(n_chunk)),
            record_stride=stride,
            ambiguous_tol=cfg4.ambiguous_basin_tolerance)
        out, matches, names = _outcome_rows(res, n_chunk)
        traj = res.trajectory
        t_raw = np.asarray(traj["time_s"])
        # Level 4 引擎记录约定：初始样本 t=0 之后是 1+stride·k 步
        # （首间隔不均匀）→ 丢弃 t=0 样本，得到均匀网格 [dt, T']
        t_rec = t_raw[1:]
        pos_rec = np.asarray(traj["position_m"])[1:]  # [n_rec, n_chunk, 3]
        controls = timeline.controls_on_grid(t_rec)
        col = {k: np.asarray(controls[k])[:, None] for k in
               ("d_slm", "d_aod", "c1", "c2", "zs1", "zs2")}
        u_slm = static_potential(
            pos_rec[..., 0] - timeline.physics.slm_center_m[0],
            pos_rec[..., 1] - timeline.physics.slm_center_m[1],
            pos_rec[..., 2], col["d_slm"],
            timeline.physics.slm_waist_m, timeline.physics.slm_zr_m)
        u_aod = lensing_potential_and_force(
            pos_rec[..., 0] - col["c1"],
            pos_rec[..., 1] - col["c2"], pos_rec[..., 2],
            col["d_aod"], timeline.physics.aod_waist_m,
            timeline.physics.aod_zr_m, col["zs1"], col["zs2"])[0]
        all_u_slm[sl] = u_slm.T.astype(np.float32)
        all_u_aod[sl] = u_aod.T.astype(np.float32)
        success[sl] = out["roundtrip_success"]
        retained[sl] = out["final_retained"]
        for i, st in enumerate(out["first_failure_stage"]):
            if str(st) not in stage_names:
                stage_names.append(str(st))
            first_stage[done + i] = stage_names.index(str(st))
        a_slm[sl] = res.phase.state["none"]["a_slm"]
        a_aod[sl] = res.phase.state["none"]["a_aod"]
        done += n_chunk
        log(f"[pool:{name}] {done}/{n_traj} 轨迹完成")
    grid = timeline.build_grid(dt)
    seg_bounds = [0.0] + [float(timeline.segment_starts[i] + seg.duration_s)
                          for i, seg in enumerate(timeline.segments)]
    total_covered = float(t_rec[-1])
    pool = TrajectoryPool(
        name=name, times_s=t_rec, u_slm=all_u_slm, u_aod=all_u_aod,
        roundtrip_success=success, final_retained=retained,
        first_failure_stage=first_stage, stage_names=stage_names,
        u_slm_final=all_u_slm[:, -1].astype(np.float64),
        u_aod_final=all_u_aod[:, -1].astype(np.float64),
        level4_a_slm=a_slm, level4_a_aod=a_aod,
        segment_boundaries_s=np.asarray(seg_bounds),
        total_duration_s=total_covered,
        segment_names=list(grid["segment_names"]),
        meta={"seed": int(seed), "n_traj": int(n_traj),
              "dt_classical_us": dt_classical_us,
              "record_stride_us": record_stride_us,
              "temperature_uK": states["temperature_uK"],
              "sampler_note": states.get("sampler_note", ""),
              "engine": "level4_roundtrip_coherence.roundtrip_simulation."
                        "simulate_roundtrip"})
    return pool


def time_info(res) -> str:
    """结果摘要字符串（兼容保留）。"""
    return f"residual_max={float(np.max(np.abs(res.residual))):.2e}J"


# --------------------------------------------------------------- joint 模式
def joint_classical_spin(timeline, cfg4, seed, n_shots, eta_slm, eta_aod,
                         dd_pulse_times_s=(), dt_us=0.1,
                         return_quantum=False, quantum_events=None):
    """联合经典-自旋推进（模式 B）：复刻 Level 4 velocity-Verlet 步进。

    每步用与 Level 4 PhaseAccumulator 相同的梯形中点差分光移推进量子
    z-相位；瞬时 DD 翻转按 Level 4 的步内分数精确拆分。返回最终相位
    Φ（toggling 积分）与可选的量子态。该函数用于校验 frozen replay，
    不用于大规模 RB。
    """
    from .qubit_hamiltonian import light_shift_detuning
    dt = dt_us * US
    states = sample_slm_states(cfg4, int(seed), n_shots)
    grid = timeline.build_grid(dt)
    times = grid["times"]
    ctrl = {k: np.asarray(grid[k]) for k in
            ("d_slm", "ddot_slm", "c1", "c2", "c1d", "c2d", "zs1", "zs2",
             "zs1d", "zs2d", "d_aod", "ddot_aod")}
    phys = timeline.physics
    mass = phys.mass_kg
    from level3_3d_transport_lensing.gaussian_3d import static_force
    pos = states["pos_m"].astype(float).copy()
    vel = states["vel_m_per_s"].astype(float).copy()
    phi = np.zeros(n_shots)
    psi = (np.zeros((n_shots, 2), dtype=complex)
           if return_quantum else None)
    if return_quantum:
        psi[:, 0] = 1.0
    pulses_sorted = np.sort(np.asarray(dd_pulse_times_s, dtype=float))
    pulse_idx = _pulses_to_step_fractions(pulses_sorted, times)

    def evaluate(i, p):
        xi_s = p[:, 0] - phys.slm_center_m[0]
        eta_s = p[:, 1] - phys.slm_center_m[1]
        z = p[:, 2]
        u_a, fxa, fya, fza = lensing_potential_and_force(
            p[:, 0] - ctrl["c1"][i], p[:, 1] - ctrl["c2"][i], z,
            ctrl["d_aod"][i], phys.aod_waist_m, phys.aod_zr_m,
            ctrl["zs1"][i], ctrl["zs2"][i])
        fxs, fys, fzs = static_force(xi_s, eta_s, z, ctrl["d_slm"][i],
                                     phys.slm_waist_m, phys.slm_zr_m)
        u_s = static_potential(xi_s, eta_s, z, ctrl["d_slm"][i],
                               phys.slm_waist_m, phys.slm_zr_m)
        acc = np.stack([fxa + fxs, fya + fys, fza + fzs], axis=1) / mass
        return acc, u_s, u_a

    acc, u_s, u_a = evaluate(0, pos)
    y = 1.0
    next_pulse = 0
    for i in range(grid["n_steps"]):
        pos = pos + vel * dt + 0.5 * acc * dt ** 2
        acc_new, u_s_new, u_a_new = evaluate(i + 1, pos)
        vel = vel + 0.5 * (acc + acc_new) * dt
        delta_mid = light_shift_detuning(
            0.5 * (u_s + u_s_new), 0.5 * (u_a + u_a_new), eta_slm, eta_aod)
        total = delta_mid * dt
        extra = [(s, f) for (s, f) in pulse_idx if s == i]
        if not extra:
            phi += y * total
            if return_quantum:
                _apply_z(psi, y * total)
        else:
            left = np.zeros(n_shots)
            for _, frac in extra:
                part = delta_mid * (frac * dt)
                phi += y * (part - left)
                if return_quantum:
                    _apply_z(psi, y * (part - left))
                left = part
                y = -y
                if return_quantum:
                    _apply_x_pi(psi)
            phi += y * (total - left)
            if return_quantum:
                _apply_z(psi, y * (total - left))
        acc, u_s, u_a = acc_new, u_s_new, u_a_new
    out = {"phi": phi, "final_pos": pos, "final_vel": vel,
           "initial_states": {k: v for k, v in states.items()
                              if isinstance(v, (int, float, str))}}
    if return_quantum:
        out["psi"] = psi
    return out


def _pulses_to_step_fractions(pulses, times):
    """脉冲时刻 → (步, 步内分数)（与 Level 4 dephasing 相同语义）。"""
    out = []
    for p in pulses:
        idx = int(np.searchsorted(times, p, side="left")) - 1
        if idx < 0 or idx >= len(times) - 1:
            raise ValueError(f"脉冲 {p:.3e} s 落在时间网格之外")
        out.append((idx, (p - times[idx]) / (times[idx + 1] - times[idx])))
    return out


def _apply_z(psi, angle):
    """批量 z 相位。"""
    half = 0.5 * np.asarray(angle, dtype=float)
    psi[:, 0] *= np.exp(-1j * half)
    psi[:, 1] *= np.exp(1j * half)


def _apply_x_pi(psi):
    """瞬时 X π 翻转。"""
    psi[:, [0, 1]] = psi[:, [1, 0]]


def replay_phase_for_pool(pool: TrajectoryPool, k_idx, eta_slm, eta_aod,
                          dd_pulse_times_s=()):
    """在冻结池上用中点矩形积分复算 toggling 相位（一致性检查用）。

    自由段相位 = Σ_sample Δ(t_mid)·Δt_record；瞬时脉冲处 y 翻转。
    """
    times = pool.times_s
    dt = pool.sample_interval_s
    pulses = np.sort(np.asarray(dd_pulse_times_s, dtype=float))
    mids = 0.5 * (times[:-1] + times[1:])
    y = np.ones(mids.size)
    if pulses.size:
        flips = np.searchsorted(mids, pulses, side="left")
        for f in flips:
            y[f:] *= -1.0
    u_slm = pool.u_slm[k_idx].astype(np.float64)
    u_aod = pool.u_aod[k_idx].astype(np.float64)
    delta = (eta_slm * u_slm + eta_aod * u_aod) / HBAR
    # 中点失谐 = 相邻样本平均（线性插值），矩形权重 dt
    delta_mid = 0.5 * (delta[:, :-1] + delta[:, 1:])
    return delta_mid @ (y * dt)
