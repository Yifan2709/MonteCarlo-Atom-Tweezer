"""势函数测试（spec §9 第 1–4 项）。"""

from __future__ import annotations

import math

import numpy as np
import pytest

from level0_static_trap.potentials import (
    gaussian_curvature_analytic,
    gaussian_force,
    gaussian_potential,
    trap_frequency_analytic,
)

U0 = 1.0e-27   # J
W = 1.0e-6     # m
XC = 0.0


def test_center_potential_is_minus_depth():
    # 第 1 项：中心势能 = -U0
    u = float(gaussian_potential(XC, U0, W, XC))
    assert u == pytest.approx(-U0, rel=1e-15, abs=0)


def test_center_force_is_zero():
    # 第 2 项：中心受力为零
    f = float(gaussian_force(XC, U0, W, XC))
    assert abs(f) == pytest.approx(0.0, abs=1e-30)


def test_force_points_inward_on_both_sides():
    # 第 3 项：中心两侧受力指向中心
    dx = 0.1 * W
    f_right = float(gaussian_force(XC + dx, U0, W, XC))  # 在中心右侧
    f_left = float(gaussian_force(XC - dx, U0, W, XC))   # 在中心左侧
    assert f_right < 0.0, "右侧力应指向左（负方向）"
    assert f_left > 0.0, "左侧力应指向右（正方向）"


def test_force_matches_finite_difference_gradient():
    # 第 4 项：gaussian_force 与势能有限差分梯度一致
    # 用两个步长验证：h 与 h/2 的有限差分都应随步长缩小而趋近解析值（O(h^2)），
    # 从而确认解析力确实是势能的导数，而非仅“巧合接近”。
    x = np.linspace(-1.5 * W, 1.5 * W, 401)
    f_analytic = gaussian_force(x, U0, W, XC)
    mask = np.abs(x) > 0.05 * W  # 排除梯度近零处相对误差放大

    def fd_err(h):
        u_plus = gaussian_potential(x + h, U0, W, XC)
        u_minus = gaussian_potential(x - h, U0, W, XC)
        grad_fd = -(u_plus - u_minus) / (2.0 * h)  # F = -dU/dx
        return np.max(np.abs((grad_fd[mask] - f_analytic[mask]) / f_analytic[mask]))

    err_h = fd_err(1e-3 * W)
    err_h2 = fd_err(0.5e-3 * W)
    # h/2 的误差应约为 h 的 1/4（中心差分二阶收敛），且本身足够小
    assert err_h2 < err_h
    assert err_h / err_h2 > 3.0  # 接近 4
    assert err_h2 < 1e-5


def test_potential_and_force_accept_arrays_and_scalars():
    xs_scalar = 0.3 * W
    xs_array = np.linspace(-W, W, 11)
    assert np.isfinite(gaussian_potential(xs_scalar, U0, W, XC))
    arr_u = gaussian_potential(xs_array, U0, W, XC)
    arr_f = gaussian_force(xs_array, U0, W, XC)
    assert arr_u.shape == xs_array.shape
    assert arr_f.shape == xs_array.shape


def test_analytic_curvature_and_frequency_signs():
    k = gaussian_curvature_analytic(U0, W)
    assert k == pytest.approx(4.0 * U0 / W**2, rel=1e-15)
    f = trap_frequency_analytic(U0, W, 1.0e-25)
    assert f > 0.0
    # f0 = sqrt(k/m)/(2pi)
    assert f == pytest.approx(math.sqrt(k / 1.0e-25) / (2 * math.pi), rel=1e-15)


def test_potential_rejects_non_positive_params():
    with pytest.raises(ValueError):
        gaussian_potential(0.0, -1.0, W, XC)
    with pytest.raises(ValueError):
        gaussian_potential(0.0, U0, -W, XC)
    with pytest.raises(ValueError):
        gaussian_force(0.0, U0, 0.0, XC)
