"""积分器测试（spec §9 第 8、10、11 项）。"""

from __future__ import annotations

import math

import numpy as np
import pytest

from level0_static_trap.analysis import (
    estimate_frequency_from_trajectory,
    estimate_nonlinear_frequency_by_quadrature,
)
from level0_static_trap.constants import u_to_kg
from level0_static_trap.integrators import IntegrationError, velocity_verlet
from level0_static_trap.potentials import gaussian_force, gaussian_potential, trap_frequency_analytic

U0 = 1.0e-27
W = 1.0e-6
MASS = u_to_kg(132.90545196)


def _force(x):
    return gaussian_force(x, U0, W, 0.0)


# 第 8 项：小振幅 (0.01w) 轨迹频率接近小振幅频率
def test_small_amplitude_trajectory_frequency_matches_analytic():
    f0 = trap_frequency_analytic(U0, W, MASS)
    amp = 0.01 * W
    period = 1.0 / f0
    dt = period / 2000.0
    n = int(40 * period / dt)
    times = np.linspace(0, n * dt, n + 1)
    _, x, v, a = velocity_verlet(_force, MASS, amp, 0.0, times)
    res = estimate_frequency_from_trajectory(times, x, minimum_complete_cycles=20)
    # 采样估计 + 极小非简谐（0.01w）+ 积分误差，放宽到 5e-3
    assert abs(res.frequency_hz - f0) / f0 < 5e-3


# 第 9 项：默认 0.10w 有限振幅——不把物理非简谐频移判为积分器失败
def test_default_amplitude_trajectory_matches_nonlinear_not_small_amp():
    f0 = trap_frequency_analytic(U0, W, MASS)
    f_nl = estimate_nonlinear_frequency_by_quadrature(U0, W, 0.0, MASS, amplitude=0.10 * W)
    amp = 0.10 * W
    period = 1.0 / f_nl
    dt = period / 2000.0
    n = int(40 * period / dt)
    times = np.linspace(0, n * dt, n + 1)
    _, x, v, a = velocity_verlet(_force, MASS, amp, 0.0, times)
    res = estimate_frequency_from_trajectory(times, x, minimum_complete_cycles=20)
    # 与非线性基准一致（轨迹积分误差，严格数值容差）
    assert abs(res.frequency_hz - f_nl) / f_nl < 1e-3
    # 与小振幅频率的差异（非简谐软化，约 0.7%，落在 1–3%）
    diff_small = abs(res.frequency_hz - f0) / f0
    assert 1e-3 < diff_small < 3e-2


# 第 10 项：位置/相位误差随步长二阶收敛（Richardson 自收敛）
def test_velocity_verlet_second_order_self_convergence():
    # 取固定末时刻，比较 dt、dt/2、dt/4 三种步长的位移
    f0 = trap_frequency_analytic(U0, W, MASS)
    amp = 0.05 * W
    period = 1.0 / f0
    T_end = 5.0 * period
    errs = []
    coarse_n = 200
    for k in (1, 2, 4):
        dt = T_end / (coarse_n * k)
        n = int(round(T_end / dt))
        times = np.linspace(0, n * dt, n + 1)
        _, x, _, _ = velocity_verlet(_force, MASS, amp, 0.0, times)
        errs.append(float(x[-1]))
    e1 = abs(errs[0] - errs[2])  # 粗 vs 细
    e2 = abs(errs[1] - errs[2])  # 中 vs 细
    ratio = e1 / e2 if e2 > 0 else math.inf
    # 二阶方法：粗/中误差比应接近 4
    assert 3.0 < ratio < 5.5, f"二阶自收敛比 {ratio} 偏离 4"


# 第 11 项：固定物理时长减半步长，能量误差按二阶缩放或两者均低于阈值
def test_energy_error_scales_with_step():
    amp = 0.10 * W
    f0 = trap_frequency_analytic(U0, W, MASS)
    period = 1.0 / f0
    T_end = 30.0 * period
    pkpk = []
    for k in (1, 2):
        dt = T_end / (4000 * k)
        n = int(round(T_end / dt))
        times = np.linspace(0, n * dt, n + 1)
        _, x, v, _ = velocity_verlet(_force, MASS, amp, 0.0, times)
        ke = 0.5 * MASS * v * v
        pe = gaussian_potential(x, U0, W, 0.0)
        total = ke + pe
        e0 = float(total[0])
        dE_over = (total - e0) / U0
        pkpk.append(float(np.max(dE_over) - np.min(dE_over)))
    # 要么二阶缩放（pkpk[0]/pkpk[1] 接近 4），要么两者都低于阈值
    scaled = pkpk[1] * 4.0
    below_threshold = (pkpk[0] < 1e-4) and (pkpk[1] < 1e-4)
    ratio_ok = pkpk[0] > scaled * 0.5  # 粗误差至少接近细误差的 4 倍的一半
    assert below_threshold or ratio_ok, (
        f"能量误差未二阶缩放且未低于阈值：pkpk(dt)={pkpk[0]:.3e}, pkpk(dt/2)={pkpk[1]:.3e}"
    )


def test_velocity_verlet_rejects_non_uniform_times():
    times = np.array([0.0, 0.1, 0.3, 0.4])  # 不等间隔
    with pytest.raises(IntegrationError, match="等间隔"):
        velocity_verlet(_force, MASS, 0.0, 0.0, times)


def test_velocity_verlet_rejects_short_times():
    with pytest.raises(IntegrationError):
        velocity_verlet(_force, MASS, 0.0, 0.0, np.array([0.0]))


def test_velocity_verlet_returns_expected_shapes():
    times = np.linspace(0, 1e-5, 11)
    t, x, v, a = velocity_verlet(_force, MASS, 1e-7, 0.0, times)
    assert t.shape == x.shape == v.shape == a.shape == (11,)


def test_velocity_verlet_harmonic_exact():
    # 谐振子参考：F = -k x，VV 在一个周期内能量近似守恒
    k = 4.0 * U0 / W**2  # 用高斯阱中心曲率作为谐振子刚度
    times = np.linspace(0, 1e-4, 4001)
    _, x, v, _ = velocity_verlet(lambda xx: -k * xx, MASS, 0.05 * W, 0.0, times)
    ke = 0.5 * MASS * v * v
    pe = 0.5 * k * x * x
    total = ke + pe
    e0 = float(total[0])
    rel = np.max(np.abs((total - e0) / e0))
    assert rel < 1e-4
