"""频率与曲率分析测试（spec §9 第 5–9 项）。"""

from __future__ import annotations

import math

import numpy as np
import pytest

from level0_static_trap.analysis import (
    calculate_relative_error,
    estimate_frequency_from_trajectory,
    estimate_nonlinear_frequency_by_quadrature,
    find_equilibrium_position,
    numerical_curvature_pair,
    numerical_trap_frequency,
)
from level0_static_trap.constants import u_to_kg
from level0_static_trap.potentials import (
    gaussian_force,
    gaussian_potential,
    trap_frequency_analytic,
)

U0 = 1.0e-27
W = 1.0e-6
MASS = u_to_kg(132.90545196)
XC_OFF = 2.5e-6  # 非零中心


def _pot(x):
    return gaussian_potential(x, U0, W, XC_OFF)


def _force(x):
    return gaussian_force(x, U0, W, XC_OFF)


# 第 5 项：非零中心时数值平衡位置接近配置中心
def test_equilibrium_recovers_nonzero_center():
    res = find_equilibrium_position(_pot, _force, XC_OFF, W)
    assert res.x_eq == pytest.approx(XC_OFF, abs=1e-15)
    assert res.force_residual == pytest.approx(0.0, abs=1e-21)
    assert res.u_eq == pytest.approx(-U0, rel=1e-12, abs=0)


# 第 6 项：h 与 h/2 的数值曲率收敛到解析曲率
def test_numerical_curvature_converges_to_analytic():
    curv = numerical_curvature_pair(_pot, U0, W, XC_OFF, step_waist=0.002)
    err_h = abs(curv.numerical_h - curv.analytic) / curv.analytic
    err_h2 = abs(curv.numerical_h2 - curv.analytic) / curv.analytic
    # h/2 误差远小于 h（O(h^2) 二阶差分）
    assert err_h2 < err_h
    assert err_h2 < 1e-4  # 容差：数值离散化误差（严格）
    assert curv.convergence_diff < abs(curv.numerical_h - curv.analytic)


# 第 7 项：数值曲率频率接近解析小振幅频率
def test_numerical_curvature_frequency_matches_analytic():
    curv = numerical_curvature_pair(_pot, U0, W, XC_OFF, step_waist=0.002)
    f_num = numerical_trap_frequency(curv.numerical_h2, MASS)
    f_ana = trap_frequency_analytic(U0, W, MASS)
    assert abs(f_num - f_ana) / f_ana < 1e-4  # 严格数值容差


def test_nonlinear_frequency_recovers_harmonic_in_small_amp_limit():
    # 谐振子极限：振幅很小 -> f_nonlinear -> f_analytic
    f_ana = trap_frequency_analytic(U0, W, MASS)
    f_nl_small = estimate_nonlinear_frequency_by_quadrature(
        U0, W, 0.0, MASS, amplitude=1e-3 * W
    )
    assert abs(f_nl_small - f_ana) / f_ana < 1e-3


def test_trajectory_frequency_matches_nonlinear_at_default_amplitude():
    # 第 9 项（数值部分）：默认 0.10w 有限振幅下，轨迹频率应与非线性积分基准一致，
    # 而不是与小振幅频率一致（差异即非简谐频移）。
    f_ana = trap_frequency_analytic(U0, W, MASS)
    f_nl = estimate_nonlinear_frequency_by_quadrature(U0, W, 0.0, MASS, amplitude=0.10 * W)
    # 非简谐软化约 0.7%
    assert (f_ana - f_nl) / f_ana > 1e-3
    assert (f_ana - f_nl) / f_ana < 3e-2

    # 用解析周期生成一条无离散误差的轨迹，验证估计器能恢复非线性频率
    period = 1.0 / f_nl
    duration = 50 * period
    n = 200_000
    t = np.linspace(0, duration, n)
    # 真实轨迹非简谐；这里用同频率纯正弦近似，仅检查估计器在已知频率下的回收
    disp = (0.10 * W) * np.sin(2 * math.pi * f_nl * t)
    res = estimate_frequency_from_trajectory(t, disp, minimum_complete_cycles=20)
    assert abs(res.frequency_hz - f_nl) / f_nl < 1e-3


def test_calculate_relative_error_handles_zero_reference():
    assert calculate_relative_error(1.0, 1.0) == pytest.approx(0.0)
    assert math.isnan(calculate_relative_error(1.0, 0.0))
    assert math.isinf(calculate_relative_error(2.0, 0.0, "return_signed_inf"))
    with pytest.raises(ValueError):
        calculate_relative_error(1.0, 0.0, "raise")
