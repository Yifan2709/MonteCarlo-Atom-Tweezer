"""Level 2C：往返重复转移的存活率引擎（逐 shot 连续模拟 n 次单程转移）。

序列（论文 Fig. 6a）：pick-up & split → 100 µs 等待（SLM 关闭，清除残留原子）
→ SLM 重开 → merge & drop-off（pick-up 轨迹的时间反演）。一次 pick-up 或
drop-off 计为一次“单程转移”；逐 shot 中途丢失即终止并记失败（不换初态）。

实现要点：
- 所有 shots 用 NumPy 数组同步推进（velocity-Verlet），丢失的 shot 在
  判定节点被移出活跃池（早停），不向池内补充新初态；
- 判定节点：等待段末（奇数次，AOD 单阱束缚）与 drop-off 末（偶数次，
  SLM 单阱束缚）。静态段内机械能守恒、噪声注入单调增能，故端点检查
  等价于“段内任何时刻逃出”；
- 位点差异/抖动以逐 shot 参数数组进入力计算（对准偏移、深度标定、束腰）；
- 加热为平均强度确定性功率（复用 level2 的 MeanHeating 三通道常数功率）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from level0_static_trap.potentials import gaussian_potential
from level2_joint_transfer.waveforms import MLCubicWaveform, PickupSequentialWaveform

K_B_UK = 1.380649e-29  # J per µK


def time_reverse_waveform(waveform):
    """构造波形的时间反演（merge & drop-off = pick-up 的时间反演）。

    PickupSequentialWaveform：同参数 direction='dropoff'；
    MLCubicWaveform：控制点值序列反转（三次插值的严格时间反演）。
    """
    if isinstance(waveform, PickupSequentialWaveform):
        direction = "dropoff" if waveform.direction == "pickup" else "pickup"
        return PickupSequentialWaveform(waveform.duration_s, waveform.ramp_fraction,
                                        waveform.final_center_m, waveform.final_depth_j,
                                        direction=direction)
    if isinstance(waveform, MLCubicWaveform):
        return MLCubicWaveform(waveform.duration_s, waveform.control_s,
                               waveform.center_values_m[::-1].copy(),
                               waveform.depth_values_j[::-1].copy(),
                               waveform.initial_center_m, waveform.initial_depth_j)
    raise ValueError(f"不支持时间反演的波形类型 {type(waveform).__name__!r}")


@dataclass(frozen=True)
class RoundtripLegs:
    """一个往返的分段表（在各自步数网格上预计算的势场程序）。"""

    center_pickup: np.ndarray      # pick-up 段 AOD 中心（含端点）
    depth_pickup: np.ndarray
    center_wait: np.ndarray        # 等待段（常数，SLM 关闭）
    depth_wait: np.ndarray
    center_dropoff: np.ndarray     # drop-off 段（时间反演）
    depth_dropoff: np.ndarray
    slm_factor: np.ndarray         # 整段拼接后的 SLM 开关（1 开 / 0 关），含段边界重复点
    steps_pickup: int
    steps_wait: int
    steps_dropoff: int


def build_roundtrip_legs(pickup_waveform, wait_s: float, dt_s: float) -> RoundtripLegs:
    """把 pick-up 波形 + 等待段 + 时间反演 drop-off 拼成一个往返程序。"""
    t_pu = pickup_waveform.duration_s
    steps_pickup = int(round(t_pu / dt_s))
    steps_wait = int(round(wait_s / dt_s))
    if not math.isclose(t_pu / dt_s, steps_pickup, abs_tol=1e-9) or \
            not math.isclose(wait_s / dt_s, steps_wait, abs_tol=1e-9):
        raise ValueError("T_pickup 与 wait 必须能被 dt 整除")
    dropoff = time_reverse_waveform(pickup_waveform)
    t1 = np.arange(steps_pickup + 1) * dt_s
    t2 = np.arange(steps_wait + 1) * dt_s
    center_pickup = np.asarray(pickup_waveform.center(t1), dtype=float)
    depth_pickup = np.asarray(pickup_waveform.depth(t1), dtype=float)
    center_wait = np.full(steps_wait + 1, center_pickup[-1])
    depth_wait = np.full(steps_wait + 1, depth_pickup[-1])
    center_dropoff = np.asarray(dropoff.center(t1), dtype=float)
    depth_dropoff = np.asarray(dropoff.depth(t1), dtype=float)
    slm_factor = np.concatenate([np.ones(steps_pickup + 1), np.zeros(steps_wait + 1),
                                 np.ones(steps_pickup + 1)])
    return RoundtripLegs(center_pickup, depth_pickup, center_wait, depth_wait,
                         center_dropoff, depth_dropoff, slm_factor,
                         steps_pickup, steps_wait, steps_pickup)


@dataclass(frozen=True)
class DisorderDraws:
    """逐 shot 的位点差异抽样（静态）与逐转移对准抖动种子。"""

    slm_depth_factor: np.ndarray      # D_SLM 逐 shot 乘子
    aod_depth_factor: np.ndarray      # AOD 深度标定逐 shot 乘子
    aod_waist_factor: np.ndarray      # AOD 束腰逐 shot 乘子
    alignment_offset_m: np.ndarray    # AOD–SLM 静态相对对准偏移
    jitter_seed: int                  # 逐转移抖动抽种子（引擎内按转移序号派生）

    @classmethod
    def disabled(cls, shots: int):
        return cls(np.ones(shots), np.ones(shots), np.ones(shots), np.zeros(shots), 0)


def draw_disorder(cfg_disorder, shots: int, seed: int) -> DisorderDraws:
    """按配置抽样位点差异；未启用时返回全 1/全 0（严格无效应）。"""
    if not cfg_disorder.enabled:
        return DisorderDraws.disabled(shots)
    rng = np.random.default_rng(seed)
    return DisorderDraws(
        slm_depth_factor=1.0 + cfg_disorder.slm_depth_rel_sigma * rng.standard_normal(shots),
        aod_depth_factor=1.0 + cfg_disorder.aod_depth_rel_sigma * rng.standard_normal(shots),
        aod_waist_factor=1.0 + cfg_disorder.aod_waist_rel_sigma * rng.standard_normal(shots),
        alignment_offset_m=cfg_disorder.alignment_sigma_um * 1e-6 * rng.standard_normal(shots),
        jitter_seed=int(rng.integers(0, 2**31 - 1)),
    )


def _leg_steps(legs: RoundtripLegs):
    """返回 [(center, depth, slm_on, judge)] 三段；judge 为该段末的单程转移判定。"""
    return [
        (legs.center_pickup, legs.depth_pickup, True, None),
        (legs.center_wait, legs.depth_wait, False, "pickup_wait"),
        (legs.center_dropoff, legs.depth_dropoff, True, "dropoff"),
    ]


def run_survival(x0_m, v0_m_s, mass_kg, legs: RoundtripLegs, dt_s, max_transfers: int,
                 slm_depth_j, slm_waist_m, aod_waist_m,
                 heating=None, disorder: DisorderDraws | None = None,
                 jitter_sigma_m: float = 0.0, record_transfers=()):
    """逐 shot 连续模拟最多 max_transfers 次单程转移（奇=pick-up，偶=drop-off）。

    返回逐 shot 的存活转移数、丢失位置/方式、逐转移存活计数曲线 S(n)、
    每次判定时的阱能量（相对阱底），以及噪声累计注入能量。
    record_transfers：额外保存这些转移边界上的 (x, v) 快照（诊断用）。
    """
    x0_m = np.asarray(x0_m, dtype=float)
    v0_m_s = np.asarray(v0_m_s, dtype=float)
    n0 = x0_m.size
    if disorder is None:
        disorder = DisorderDraws.disabled(n0)
    jitter_rng = np.random.default_rng(disorder.jitter_seed) if jitter_sigma_m > 0 else None

    x = x0_m.copy()
    v = v0_m_s.copy()
    shot_ids = np.arange(n0)
    noise_total = np.zeros(n0)          # 逐 shot 累计注入能量（J）
    transfers_survived = np.zeros(n0, dtype=int)
    lost_leg = np.array([""] * n0, dtype=object)
    lost_transfer = np.zeros(n0, dtype=int)
    judgement_energy_over_depth: list[list[float]] = [[] for _ in range(n0)]
    snapshots: dict[int, dict[str, np.ndarray]] = {}
    survival_counts = np.zeros(max_transfers + 1, dtype=int)
    survival_counts[0] = n0

    slm_depths = slm_depth_j * disorder.slm_depth_factor
    aod_waists = aod_waist_m * disorder.aod_waist_factor

    def aod_force(xarr, center, depth_nominal):
        q = (xarr - center) / aod_waists
        return -(4.0 * depth_nominal * disorder.aod_depth_factor / aod_waists) * q * np.exp(-2.0 * q**2)

    def slm_force(xarr):
        return -(4.0 * slm_depths / slm_waist_m) * (xarr / slm_waist_m) \
            * np.exp(-2.0 * (xarr / slm_waist_m) ** 2)

    # 每个往返含 2 次单程转移；ceil 覆盖奇数上限
    roundtrips = (max_transfers + 1) // 2
    transfer_index = 0
    for rt in range(roundtrips):
        if x.size == 0 or transfer_index >= max_transfers:
            break
        for center_arr, depth_arr, slm_on, judge in _leg_steps(legs):
            if transfer_index >= max_transfers:
                break
            steps = center_arr.size - 1
            # 逐段（=逐半次转移）按当前活跃池尺寸重抽抖动，避免中途收缩后尺寸错位
            jitter = jitter_sigma_m * jitter_rng.standard_normal(x.size) \
                if jitter_rng is not None else 0.0
            offsets = disorder.alignment_offset_m + (jitter if np.ndim(jitter) else 0.0)
            accel = (slm_force(x) if slm_on else 0.0) \
                + aod_force(x, center_arr[0] + offsets, depth_arr[0])
            accel = accel / mass_kg
            for i in range(steps):
                x = x + v * dt_s + 0.5 * accel * dt_s**2
                accel_new = (slm_force(x) if slm_on else 0.0) \
                    + aod_force(x, center_arr[i + 1] + offsets, depth_arr[i + 1])
                accel_new = accel_new / mass_kg
                v = v + 0.5 * (accel + accel_new) * dt_s
                accel = accel_new
                if heating is not None:
                    powers = heating.channel_powers(x, v, 0.0, np.zeros_like(x))
                    delta = powers["total"] * dt_s
                    if np.any(delta != 0.0):
                        v = np.copysign(np.sqrt(np.maximum(0.0, v**2 + 2.0 * delta / mass_kg)), v)
                        noise_total = noise_total + delta
            if judge is None:
                continue
            transfer_index += 1
            # --- 判定 ---
            if judge == "pickup_wait":
                depth_end = depth_arr[-1] * disorder.aod_depth_factor
                center_end = center_arr[-1] + offsets
                energy = 0.5 * mass_kg * v**2 \
                    + gaussian_potential(x, depth_end, aod_waists, center_end)
                depth_scale = depth_end
            else:
                energy = 0.5 * mass_kg * v**2 + gaussian_potential(x, slm_depths, slm_waist_m, 0.0)
                depth_scale = slm_depths
            bound = energy < 0.0
            for pos, sid in enumerate(shot_ids):
                judgement_energy_over_depth[sid].append(
                    float(energy[pos] / depth_scale[pos] + 1.0))
            transfers_survived[shot_ids[bound]] = transfer_index
            survival_counts[transfer_index] = int(bound.sum())
            if transfer_index in record_transfers:
                snapshots[transfer_index] = {"x_m": x.copy(), "v_m_s": v.copy(),
                                             "shot_ids": shot_ids.copy()}
            if not bound.all():
                lost_now = shot_ids[~bound]
                lost_leg[lost_now] = judge
                lost_transfer[lost_now] = transfer_index
                keep = bound
                x, v = x[keep], v[keep]
                shot_ids = shot_ids[keep]
                disorder = DisorderDraws(disorder.slm_depth_factor[keep],
                                         disorder.aod_depth_factor[keep],
                                         disorder.aod_waist_factor[keep],
                                         disorder.alignment_offset_m[keep],
                                         disorder.jitter_seed)
                slm_depths = slm_depth_j * disorder.slm_depth_factor
                aod_waists = aod_waist_m * disorder.aod_waist_factor
            if x.size == 0:
                break

    # 未完成的转移数（提前团灭时 survival_counts 尾部已为 0）
    survival_fraction = survival_counts / n0
    return {
        "shots": n0,
        "max_transfers": int(max_transfers),
        "transfers_survived": transfers_survived,
        "lost_transfer": lost_transfer,
        "lost_leg": lost_leg.tolist(),
        "survival_counts": survival_counts,
        "survival_fraction": survival_fraction,
        "judgement_energy_over_depth": judgement_energy_over_depth,
        "noise_total_energy_uK": noise_total / K_B_UK,
        "snapshots": snapshots,
        "final_active": int(x.size),
    }


def run_roundtrip_detailed(x0_m, v0_m_s, mass_kg, legs: RoundtripLegs, dt_s,
                           slm_depth_j, slm_waist_m, aod_waist_m, aod_depth_nominal_j,
                           roundtrips: int = 1):
    """单 shot 全量记录版：返回整段 (t, x, v) 与每步势场参数，供功-能与时间反演检查。

    保守模型（无噪声无 disorder）。slm_factor 跳变处的外部功跳变由
    generalized_work_energy 另行核算。
    """
    t_leg_end = [legs.steps_pickup, legs.steps_pickup + legs.steps_wait,
                 legs.steps_pickup + legs.steps_wait + legs.steps_dropoff]
    steps_rt = t_leg_end[-1]
    total_steps = steps_rt * roundtrips
    times = np.arange(total_steps + 1) * dt_s
    x = np.empty(total_steps + 1)
    v = np.empty(total_steps + 1)
    aod_center = np.empty(total_steps + 1)
    aod_depth = np.empty(total_steps + 1)
    slm_factor = np.empty(total_steps + 1)

    legs_arrays = [(legs.center_pickup, legs.depth_pickup, 1.0),
                   (legs.center_wait, legs.depth_wait, 0.0),
                   (legs.center_dropoff, legs.depth_dropoff, 1.0)]

    x[0], v[0] = x0_m, v0_m_s
    idx = 0

    def force(xv, center, depth, factor):
        f = -(4.0 * depth / aod_waist_m) * ((xv - center) / aod_waist_m) \
            * math.exp(-2.0 * ((xv - center) / aod_waist_m) ** 2)
        if factor:
            f += -(4.0 * slm_depth_j / slm_waist_m) * (xv / slm_waist_m) \
                * math.exp(-2.0 * (xv / slm_waist_m) ** 2)
        return f

    for _ in range(roundtrips):
        for center_arr, depth_arr, factor in legs_arrays:
            steps = center_arr.size - 1
            a_now = force(x[idx], center_arr[0], depth_arr[0], factor) / mass_kg
            for i in range(steps):
                aod_center[idx] = center_arr[i]
                aod_depth[idx] = depth_arr[i]
                slm_factor[idx] = factor
                x[idx + 1] = x[idx] + v[idx] * dt_s + 0.5 * a_now * dt_s**2
                a_next = force(x[idx + 1], center_arr[i + 1], depth_arr[i + 1], factor) / mass_kg
                v[idx + 1] = v[idx] + 0.5 * (a_now + a_next) * dt_s
                a_now = a_next
                idx += 1
    aod_center[total_steps] = legs_arrays[-1][0][-1]
    aod_depth[total_steps] = legs_arrays[-1][1][-1]
    slm_factor[total_steps] = 1.0
    return {"time_s": times, "x_m": x, "v_m_s": v, "aod_center_m": aod_center,
            "aod_depth_j": aod_depth, "slm_factor": slm_factor,
            "steps_per_roundtrip": steps_rt,
            "leg_end_steps": [t_leg_end[0], t_leg_end[1], t_leg_end[2]]}
