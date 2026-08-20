"""Level 4 单轮/多轮向量化三维往返推进引擎。

总势 U_tot = U_SLM(三维高斯，静止中心) + U_AOD(Level 3 crossed-AOD 透镜势)。
velocity-Verlet 推进，外部功分解为 SLM 深度功、AOD 深度功、横向移动功与
lensing 焦点偏移功四通道，梯形求积；相位按 dephasing.PhaseAccumulator
与轨迹同步累积。短距离 transfer 与长运输使用同一套三维势实现。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from level3_3d_transport_lensing.aod_lensing import lensing_potential_and_force
from level3_3d_transport_lensing.gaussian_3d import static_force
from .dephasing import PhaseAccumulator
from .trap_assignment import classify_basin

CTRL_KEYS = ("d_slm", "ddot_slm", "c1", "c2", "c1d", "c2d",
             "zs1", "zs2", "zs1d", "zs2d", "d_aod", "ddot_aod")


def repeat_grid(grid: dict, n_repeats: int) -> dict:
    """把单轮控制网格重复 n 轮（多轮状态机共用同一条时间轴）。

    拼接时去掉重复的轮次边界样本：ctrl 数组长度必须与 times 一致
    （n_repeats·n_steps+1），否则边界之后控制量整体错位一步。
    """
    n_steps = grid["n_steps"]
    out = dict(grid)
    for key in CTRL_KEYS + ("segment_index",):
        base = np.asarray(grid[key])
        # n-1 个去掉末样本的整段 + 完整末段，保证总长 n·n_steps+1
        out[key] = np.concatenate([base[:-1]] * (n_repeats - 1) + [base])
    out["times"] = np.arange(n_steps * n_repeats + 1) * grid["dt"]
    out["n_steps"] = n_steps * n_repeats
    out["n_repeats"] = n_repeats
    seg_names = grid["segment_names"]
    out["segment_names"] = seg_names * n_repeats
    out["segment_kinds"] = grid["segment_kinds"] * n_repeats
    out["segment_steps"] = [r * n_steps + s for r in range(n_repeats)
                            for s in grid["segment_steps"]]
    out["round_index"] = (np.arange(n_steps * n_repeats + 1) // n_steps)
    checkpoints = []
    for r in range(n_repeats):
        for ck in grid["checkpoints"]:
            checkpoints.append(dict(ck, step=ck["step"] + r * n_steps,
                                    time_s=(ck["step"] + r * n_steps) * grid["dt"],
                                    round=r + 1))
    out["checkpoints"] = checkpoints
    return out


class _Engine:
    """一批初态在一组（可重复）控制网格上的推进器。"""

    def __init__(self, grid: dict, phys, dd_pulse_sets: dict, noise=None,
                 keep_ids=(), record_stride=0, ambiguous_tol=1e-3,
                 phase_snap_steps=()):
        self.grid = grid
        self.phys = phys
        self.noise = noise
        self.times = grid["times"]
        self.dt = grid["dt"]
        self.n_steps = grid["n_steps"]
        self.ctrl = {key: np.asarray(grid[key]) for key in CTRL_KEYS}
        self.pulse_sets = dd_pulse_sets
        self.ambiguous_tol = ambiguous_tol
        self.keep_ids = list(keep_ids)
        self.record_stride = int(record_stride)
        self.phase_snap_steps = dict(phase_snap_steps)  # step -> 标签
        self.phase_snaps: dict = {}
        self.traj = {"time_s": [], "position_m": [], "velocity_m_per_s": [],
                     "w_total_j": [], "energy_j": []}

    def run(self, pos0, vel0):
        """执行推进，返回 RoundtripResult。"""
        n_shots = len(pos0)
        self.phase = PhaseAccumulator(self.pulse_sets, self.times, n_shots,
                                      self.grid["segment_index"])
        phys = self.phys
        mass = phys.mass_kg
        w_s, zr_s = phys.slm_waist_m, phys.slm_zr_m
        w_a, zr_a = phys.aod_waist_m, phys.aod_zr_m
        slm_x, slm_y = phys.slm_center_m
        dt = self.dt
        times = self.times
        ctrl = self.ctrl
        pos = pos0.astype(float).copy()
        vel = vel0.astype(float).copy()
        if self.noise is not None:
            if self.noise.n_shots != n_shots:
                raise ValueError("噪声 shot 数与初态数不一致")
        mult_s = np.ones(n_shots) if self.noise is None else None
        mult_a = np.ones(n_shots) if self.noise is None else None

        def controls_at(i):
            return {key: float(ctrl[key][i]) for key in CTRL_KEYS}

        def evaluate(i, eps_s, eps_a, eps_s_old=None, eps_a_old=None):
            c = ctrl
            xi_s = pos[:, 0] - slm_x
            eta_s = pos[:, 1] - slm_y
            z = pos[:, 2]
            m_s = 1.0 if eps_s is None else 1.0 + eps_s
            m_a = 1.0 if eps_a is None else 1.0 + eps_a
            d_s = float(c["d_slm"][i]) * m_s
            d_a = float(c["d_aod"][i]) * m_a
            xi_a = pos[:, 0] - float(c["c1"][i])
            eta_a = pos[:, 1] - float(c["c2"][i])
            zs1, zs2 = float(c["zs1"][i]), float(c["zs2"][i])
            u_a, fxa, fya, fza = lensing_potential_and_force(
                xi_a, eta_a, z, d_a, w_a, zr_a, zs1, zs2)
            fxs, fys, fzs = static_force(xi_s, eta_s, z, d_s, w_s, zr_s)
            acc = np.stack([fxa + fxs, fya + fys, fza + fzs], axis=1) / mass
            # ---- 功率分解 ----
            move = fxa * float(c["c1d"][i]) + fya * float(c["c2d"][i])
            dz1 = (z - zs1) / zr_a
            dz2 = (z - zs2) / zr_a
            a1 = 1.0 + dz1 ** 2
            a2 = 1.0 + dz2 ** 2
            g1 = np.exp(-2.0 * xi_a ** 2 / (w_a ** 2 * a1)) / np.sqrt(a1)
            g2 = np.exp(-2.0 * eta_a ** 2 / (w_a ** 2 * a2)) / np.sqrt(a2)
            dg1_dz = g1 * dz1 / zr_a * (-1.0 / a1 + 4.0 * xi_a ** 2
                                        / (w_a ** 2 * a1 ** 2))
            dg2_dz = g2 * dz2 / zr_a * (-1.0 / a2 + 4.0 * eta_a ** 2
                                        / (w_a ** 2 * a2 ** 2))
            zs1d, zs2d = float(c["zs1d"][i]), float(c["zs2d"][i])
            shift = d_a * (dg1_dz * g2 * zs1d + g1 * dg2_dz * zs2d)
            g1g2 = g1 * g2
            ddot_a = float(c["ddot_aod"][i]) * (1.0 if eps_a is None
                                                else 1.0 + eps_a)
            if eps_a is not None and eps_a_old is not None:
                ddot_a = ddot_a + float(c["d_aod"][i]) * (eps_a - eps_a_old) / dt
            aod_depth_power = -ddot_a * g1g2
            den_s = 1.0 + (z / zr_s) ** 2
            g_slm = np.exp(-2.0 * (xi_s ** 2 + eta_s ** 2)
                           / (w_s ** 2 * den_s)) / den_s
            ddot_s = float(c["ddot_slm"][i]) * (1.0 if eps_s is None
                                                else 1.0 + eps_s)
            if eps_s is not None and eps_s_old is not None:
                ddot_s = ddot_s + float(c["d_slm"][i]) * (eps_s - eps_s_old) / dt
            slm_depth_power = -ddot_s * g_slm
            u_s = -d_s * g_slm
            return acc, u_s, u_a, (move, shift, aod_depth_power, slm_depth_power)

        # ---- 初值 ----
        if self.noise is not None:
            mult_s, mult_a = self.noise.depth_multipliers()
        else:
            mult_s = mult_a = None
        acc, u_s, u_a, powers = evaluate(0, mult_s, mult_a)
        e0 = 0.5 * mass * np.sum(vel ** 2, axis=1) + u_s + u_a
        w_move = np.zeros(n_shots)
        w_shift = np.zeros(n_shots)
        w_aod_depth = np.zeros(n_shots)
        w_slm_depth = np.zeros(n_shots)
        max_speed = float(np.max(np.linalg.norm(vel, axis=1))) if n_shots else 0.0
        max_radial = np.zeros(n_shots)
        max_axial = np.zeros(n_shots)
        checkpoint_records = []
        ck_by_step = {ck["step"]: ck for ck in self.grid["checkpoints"]}
        seg_marks = []
        seg_mark_steps = set(int(s) for s in self.grid["segment_steps"])
        seg_mark_steps.discard(0)
        self._record(0, pos, vel, e0, np.zeros(n_shots))
        checkpoint_steps = set(ck_by_step)
        for i in range(self.n_steps):
            ip1 = i + 1
            pos = pos + vel * dt + 0.5 * acc * dt ** 2
            if self.noise is not None:
                eps_s_old, eps_a_old = mult_s, mult_a
                self.noise.advance(dt)
                mult_s, mult_a = self.noise.depth_multipliers()
            else:
                eps_s_old = eps_a_old = None
            acc_new, u_s_new, u_a_new, powers_new = evaluate(
                ip1, mult_s, mult_a, eps_s_old, eps_a_old)
            vel = vel + 0.5 * (acc + acc_new) * dt
            acc = acc_new
            w_move += 0.5 * (powers[0] + powers_new[0]) * dt
            w_shift += 0.5 * (powers[1] + powers_new[1]) * dt
            w_aod_depth += 0.5 * (powers[2] + powers_new[2]) * dt
            w_slm_depth += 0.5 * (powers[3] + powers_new[3]) * dt
            self.phase.step(i, u_s, u_a, u_s_new, u_a_new)
            if ip1 in self.phase_snap_steps:
                self.phase_snaps[self.phase_snap_steps[ip1]] = \
                    self.phase.snapshot()
            u_s, u_a = u_s_new, u_a_new
            powers = powers_new
            speed = np.linalg.norm(vel, axis=1)
            max_speed = max(max_speed, float(np.max(speed)) if n_shots else 0.0)
            radial = np.hypot(pos[:, 0] - float(ctrl["c1"][ip1]),
                              pos[:, 1] - float(ctrl["c2"][ip1]))
            axial = np.abs(pos[:, 2] - 0.5 * (float(ctrl["zs1"][ip1])
                                              + float(ctrl["zs2"][ip1])))
            np.maximum(max_radial, radial, out=max_radial)
            np.maximum(max_axial, axial, out=max_axial)
            e_now = 0.5 * mass * np.sum(vel ** 2, axis=1) + u_s + u_a
            w_tot_now = w_move + w_shift + w_aod_depth + w_slm_depth
            if self.record_stride and (i % self.record_stride == 0):
                self._record(ip1, pos, vel, e_now, w_tot_now)
            if ip1 in seg_mark_steps or ip1 in checkpoint_steps:
                e_now_mark = e_now
                seg_marks.append({
                    "step": ip1,
                    "segment_index": int(self.grid["segment_index"][ip1 - 1]),
                    "e_total": e_now_mark.copy(),
                    "cum_w_move": w_move.copy(),
                    "cum_w_shift": w_shift.copy(),
                    "cum_w_aod_depth": w_aod_depth.copy(),
                    "cum_w_slm_depth": w_slm_depth.copy(),
                })
            if ip1 in checkpoint_steps:
                ck = ck_by_step[ip1]
                ctrl_now = controls_at(ip1)
                classification = classify_basin(
                    pos, vel, ctrl_now, phys,
                    ambiguous_tol=self.ambiguous_tol)
                checkpoint_records.append({
                    "checkpoint": ck, "classification": classification,
                    "pos": pos.copy(), "vel": vel.copy(),
                    "e_total": e_now.copy(),
                    "cum_work": w_tot_now.copy(),
                    "cum_w_move": w_move.copy(),
                    "cum_w_shift": w_shift.copy(),
                    "cum_w_aod_depth": w_aod_depth.copy(),
                    "cum_w_slm_depth": w_slm_depth.copy(),
                })
        e_final = 0.5 * mass * np.sum(vel ** 2, axis=1) + u_s + u_a
        w_total = w_move + w_shift + w_aod_depth + w_slm_depth
        residual = e_final - e0 - w_total
        seg_marks.append({
            "step": self.n_steps,
            "segment_index": int(self.grid["segment_index"][-1]),
            "e_total": e_final.copy(),
            "cum_w_move": w_move.copy(),
            "cum_w_shift": w_shift.copy(),
            "cum_w_aod_depth": w_aod_depth.copy(),
            "cum_w_slm_depth": w_slm_depth.copy(),
        })
        initial_mark = {
            "step": 0, "segment_index": 0,
            "e_total": e0.copy(),
            "cum_w_move": np.zeros(n_shots), "cum_w_shift": np.zeros(n_shots),
            "cum_w_aod_depth": np.zeros(n_shots),
            "cum_w_slm_depth": np.zeros(n_shots),
        }
        return RoundtripResult(
            final_pos=pos, final_vel=vel, e0=e0, e_final=e_final,
            w_move=w_move, w_shift=w_shift, w_aod_depth=w_aod_depth,
            w_slm_depth=w_slm_depth, w_total=w_total, residual=residual,
            max_speed=max_speed, max_radial=max_radial, max_axial=max_axial,
            checkpoints=checkpoint_records, phase=self.phase,
            trajectory=dict(self.traj) if self.traj["time_s"] else None,
            segment_marks=[initial_mark] + sorted(
                seg_marks, key=lambda m: m["step"]),
            phase_snaps=self.phase_snaps)

    def _record(self, step, pos, vel, energy, w_total):
        if not self.record_stride or not self.keep_ids:
            return
        ids = self.keep_ids
        self.traj["time_s"].append(self.times[step])
        self.traj["position_m"].append(pos[ids].copy())
        self.traj["velocity_m_per_s"].append(vel[ids].copy())
        self.traj["w_total_j"].append(w_total[ids].copy())
        self.traj["energy_j"].append(energy[ids].copy())


@dataclass
class RoundtripResult:
    """一轮（或多轮合并）推进结果。"""

    final_pos: np.ndarray
    final_vel: np.ndarray
    e0: np.ndarray
    e_final: np.ndarray
    w_move: np.ndarray
    w_shift: np.ndarray
    w_aod_depth: np.ndarray
    w_slm_depth: np.ndarray
    w_total: np.ndarray
    residual: np.ndarray
    max_speed: float
    max_radial: np.ndarray
    max_axial: np.ndarray
    checkpoints: list
    phase: PhaseAccumulator
    trajectory: dict | None
    segment_marks: list
    phase_snaps: dict


def simulate_roundtrip(timeline, dt, pos0, vel0, phys, dd_pulse_sets,
                       noise=None, keep_ids=(), record_stride=0,
                       ambiguous_tol=1e-3, n_repeats=1, phase_snap_steps=()):
    """构建网格并推进 n_repeats 轮（控制按轮重复，状态连续传递）。"""
    grid = timeline.build_grid(dt)
    if n_repeats > 1:
        grid = repeat_grid(grid, n_repeats)
    if noise is not None:
        noise.n_shots = len(pos0)
    engine = _Engine(grid, phys, dd_pulse_sets, noise=noise,
                     keep_ids=keep_ids, record_stride=record_stride,
                     ambiguous_tol=ambiguous_tol,
                     phase_snap_steps=phase_snap_steps)
    return engine.run(pos0, vel0)


def integrate_static(pos, vel, ctrl, phys, duration_s, dt):
    """冻结控制的静态推进（诊断副本归属判定用）。"""
    n = int(round(duration_s / dt))
    mass = phys.mass_kg
    p = pos.astype(float).copy()
    v = vel.astype(float).copy()
    d_s = float(ctrl["d_slm"])
    d_a = float(ctrl["d_aod"])
    c1, c2 = float(ctrl["c1"]), float(ctrl["c2"])

    def accel(p):
        xi_s = p[:, 0] - phys.slm_center_m[0]
        eta_s = p[:, 1] - phys.slm_center_m[1]
        u_a, fxa, fya, fza = lensing_potential_and_force(
            p[:, 0] - c1, p[:, 1] - c2, p[:, 2], d_a, phys.aod_waist_m,
            phys.aod_zr_m, 0.0, 0.0)
        fxs, fys, fzs = static_force(xi_s, eta_s, p[:, 2], d_s,
                                     phys.slm_waist_m, phys.slm_zr_m)
        return np.stack([fxa + fxs, fya + fys, fza + fzs], axis=1) / mass

    a = accel(p)
    for _ in range(n):
        p = p + v * dt + 0.5 * a * dt ** 2
        a_new = accel(p)
        v = v + 0.5 * (a + a_new) * dt
        a = a_new
    return p, v
