"""紧凑三维动力学：速度 Verlet + 高斯阱（可含 crossed-AOD 声透镜）。

复用 level3_3d_transport_lensing 的解析力函数（同一物理与单位约定）：
- 无透镜：static_force（标准三维高斯）
- 含透镜：lensing_potential_and_force，轴间焦点偏移 zs_i = τ_eff·zR·ẋ_i
  （level3 表达 zs=2σ·zR·ċ/v_s 的等价参数化，τ_eff=2σ/v_s）

损失判据：总能量（动能+全部阱势，阱底为 0、无穷远为 0）E ≥ 0 判失阱；
与 continuous_transfer 的吸收性损失同口径。能量在检查点（交接/移动边界）评估。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from level3_3d_transport_lensing.aod_lensing import lensing_potential_and_force
from level3_3d_transport_lensing.gaussian_3d import static_force, static_potential

KB = 1.380649e-23


@dataclass
class TrapPhysics:
    mass_kg: float
    depth_j: float
    waist_m: float
    zr_m: float
    # 静态第二阱（SLM 存储阱），depth_j=0 表示关闭
    slm_depth_j: float = 0.0
    slm_offset_m: float = 0.0        # AOD 与 SLM 中心失配（逐原子可不同）
    lensing_tau_eff_s_per_m: float = 0.0   # 0 → 无透镜

    @property
    def lensing_on(self) -> bool:
        return self.lensing_tau_eff_s_per_m > 0.0


@dataclass
class Segment:
    """一段控制：移动轨迹或深度 ramp。cx/cy 为 AOD 中心，depth_aod 为深度数组。"""
    name: str
    t_s: np.ndarray
    cx_m: np.ndarray
    cy_m: np.ndarray
    depth_aod_j: np.ndarray
    vx_m_s: np.ndarray = field(default_factory=lambda: np.zeros(1))
    vy_m_s: np.ndarray = field(default_factory=lambda: np.zeros(1))
    slm_depth_j: np.ndarray | None = None   # 时变 SLM 深度；None 用 ph 默认
    judge: bool = True               # 段末是否判定失阱


def total_potential(x, y, z, center, depth_j, ph: TrapPhysics,
                    lens_zs=(0.0, 0.0), slm_depth_j: float | None = None):
    """AOD 阱（可选透镜）+ SLM 静态阱的总势，零点在阱底/无穷远一致。"""
    cx, cy = center
    if ph.lensing_on:
        u, _, _, _ = lensing_potential_and_force(
            x - cx, y - cy, z, depth_j, ph.waist_m, ph.zr_m, *lens_zs)
    else:
        u = static_potential(x - cx, y - cy, z, depth_j, ph.waist_m, ph.zr_m)
    sd = ph.slm_depth_j if slm_depth_j is None else slm_depth_j
    if sd > 0:
        u = u + static_potential(x - ph.slm_offset_m, y, z, sd,
                                 ph.waist_m, ph.zr_m)
    return u


def _forces(x, y, z, center, depth_j, ph, lens_zs, slm_depth_j=None):
    cx, cy = center
    if ph.lensing_on:
        _, fx, fy, fz = lensing_potential_and_force(
            x - cx, y - cy, z, depth_j, ph.waist_m, ph.zr_m, *lens_zs)
    else:
        fx, fy, fz = static_force(x - cx, y - cy, z, depth_j, ph.waist_m,
                                  ph.zr_m)
    sd = ph.slm_depth_j if slm_depth_j is None else slm_depth_j
    if sd > 0:
        sx, sy, sz = static_force(x - ph.slm_offset_m, y, z, sd,
                                  ph.waist_m, ph.zr_m)
        fx, fy, fz = fx + sx, fy + sy, fz + sz
    return fx, fy, fz


def sample_bound_thermal_3d(seed: int, temperature_uK: float, n: int,
                            ph: TrapPhysics, max_attempts: int = 200):
    """三维谐振子高斯提议 + 全高斯势束缚筛选（E<0，阱中心坐标系）。"""
    rng = np.random.default_rng(seed)
    thermal = KB * temperature_uK * 1e-6
    sigma_r = np.sqrt(thermal * ph.waist_m**2 / (4 * ph.depth_j))
    sigma_z = np.sqrt(thermal * ph.zr_m**2 / (2 * ph.depth_j))
    sigma_v = np.sqrt(thermal / ph.mass_kg)
    pos = np.zeros((n, 3))
    vel = np.zeros((n, 3))
    pending = np.arange(n)
    for _ in range(max_attempts):
        k = len(pending)
        pos[pending, 0] = sigma_r * rng.standard_normal(k)
        pos[pending, 1] = sigma_r * rng.standard_normal(k)
        pos[pending, 2] = sigma_z * rng.standard_normal(k)
        vel[pending] = sigma_v * rng.standard_normal((k, 3))
        e = (0.5 * ph.mass_kg * np.sum(vel**2, axis=1)
             + total_potential(pos[:, 0], pos[:, 1], pos[:, 2], (0.0, 0.0),
                               ph.depth_j, ph))
        pending = pending[(e[pending] >= 0) | ~np.isfinite(e[pending])]
        if not len(pending):
            break
    if len(pending):
        raise RuntimeError(f"{len(pending)} 个初态无法束缚")
    return pos, vel


def propagate(pos0, vel0, segments: list[Segment], ph, dt_s: float,
              offsets_m: np.ndarray | None = None):
    """ph 可为单个 TrapPhysics 或与 segments 等长的列表（逐段切换，
    如交接后关闭 SLM）。"""
    if not isinstance(ph, (list, tuple)):
        ph = [ph] * len(segments)
    if len(ph) != len(segments):
        raise ValueError("ph 列表长度必须等于段数")
    """沿段序列传播。返回逐检查点存活与末态。

    offsets_m：逐原子 AOD-SLM 中心失配（叠加在 cx 上模拟对准误差）。
    透镜焦点偏移由控制段给出的轴速度驱动：zs_x=τ_eff·zR·vx。
    """
    pos = np.array(pos0, dtype=float)
    vel = np.array(vel0, dtype=float)
    n0 = len(pos)
    offs = np.zeros(n0) if offsets_m is None else np.asarray(offsets_m, float)
    alive = np.ones(n0, dtype=bool)
    checkpoints: list[dict] = []
    stage_of_atom = np.full(n0, "", dtype=object)

    for seg, ph in zip(segments, ph):
        steps = len(seg.t_s) - 1
        if steps <= 0:
            continue
        cx = seg.cx_m[:, None] + offs[None, :]      # (steps+1, atoms)
        cy = np.broadcast_to(seg.cy_m[:, None], cx.shape)
        depth = seg.depth_aod_j
        if len(seg.vx_m_s) != len(seg.t_s):
            vx = np.zeros_like(seg.t_s)
        else:
            vx = seg.vx_m_s
        if len(seg.vy_m_s) != len(seg.t_s):
            vy = np.zeros_like(seg.t_s)
        else:
            vy = seg.vy_m_s
        zs_x = ph.lensing_tau_eff_s_per_m * ph.zr_m * vx if ph.lensing_on \
            else np.zeros_like(vx)
        zs_y = ph.lensing_tau_eff_s_per_m * ph.zr_m * vy if ph.lensing_on \
            else np.zeros_like(vy)

        slm_d = seg.slm_depth_j

        def accel_at(i, x, y, z):
            sd = None if slm_d is None else slm_d[i]
            fx, fy, fz = _forces(x, y, z, (cx[i], cy[i]), depth[i], ph,
                                 (zs_x[i], zs_y[i]), sd)
            return fx / ph.mass_kg, fy / ph.mass_kg, fz / ph.mass_kg

        ax, ay, az = accel_at(0, pos[:, 0], pos[:, 1], pos[:, 2])
        for i in range(steps):
            dt_s = float(seg.t_s[i + 1] - seg.t_s[i])
            if dt_s <= 0:
                raise ValueError("Segment times must be strictly increasing")
            pos[:, 0] += vel[:, 0] * dt_s + 0.5 * ax * dt_s**2
            pos[:, 1] += vel[:, 1] * dt_s + 0.5 * ay * dt_s**2
            pos[:, 2] += vel[:, 2] * dt_s + 0.5 * az * dt_s**2
            axn, ayn, azn = accel_at(i + 1, pos[:, 0], pos[:, 1], pos[:, 2])
            vel[:, 0] += 0.5 * (ax + axn) * dt_s
            vel[:, 1] += 0.5 * (ay + ayn) * dt_s
            vel[:, 2] += 0.5 * (az + azn) * dt_s
            ax, ay, az = axn, ayn, azn
        if seg.judge:
            sd_end = None if slm_d is None else slm_d[-1]
            u = total_potential(pos[:, 0], pos[:, 1], pos[:, 2],
                                (cx[-1], cy[-1]), depth[-1], ph,
                                (zs_x[-1], zs_y[-1]), sd_end)
            # A single uniformly translating trap must be judged in its rest frame.
            # Two traps moving relative to one another have no common rest frame.
            moving = abs(vx[-1]) + abs(vy[-1]) > 1e-12
            if moving and (ph.slm_depth_j if sd_end is None else sd_end) > 0:
                raise ValueError("Judge multi-trap handoffs only at rest")
            relative_vel = vel - np.array([vx[-1], vy[-1], 0.])
            e = 0.5 * ph.mass_kg * np.sum(relative_vel**2, axis=1) + u
            lost_now = alive & (e >= 0)
            stage_of_atom[lost_now] = seg.name
            alive &= e < 0
            checkpoints.append({
                "segment": seg.name, "t_s": float(seg.t_s[-1]),
                "alive": int(alive.sum()),
                "mean_energy_bottom_uK": float(
                    np.mean(e[alive]) / KB / 1e-6) if alive.any() else float("nan"),
            })
    ph_last = ph[-1] if isinstance(ph, (list, tuple)) else ph
    u = total_potential(pos[:, 0], pos[:, 1], pos[:, 2],
                       (segments[-1].cx_m[-1] + offs, segments[-1].cy_m[-1]),
                       segments[-1].depth_aod_j[-1], ph_last, (zs_x[-1], zs_y[-1]),
                       None if slm_d is None else slm_d[-1])
    relative_vel = vel - np.array([vx[-1], vy[-1], 0.])
    e_final = 0.5 * ph_last.mass_kg * np.sum(relative_vel**2, axis=1) + u
    return {
        "n0": n0,
        "alive": alive,
        "final_pos": pos, "final_vel": vel,
        "final_energy_j": e_final,
        "loss_stage": stage_of_atom,
        "checkpoints": checkpoints,
        "survival_final": float(alive.sum()) / n0,
    }


def excitation_quanta(ph: TrapPhysics, omega_r_rad_s: float,
                      energies_j: np.ndarray) -> np.ndarray:
    """能量 → 谐振子量子数估计 N = E/(ħω_r)（R 论文 ΔN 口径）。"""
    hbar = 1.054571817e-34
    return energies_j / (hbar * omega_r_rad_s)
