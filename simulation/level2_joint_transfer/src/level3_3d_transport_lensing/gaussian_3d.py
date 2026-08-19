"""Level 3 三维物理：标准三维高斯光镊势、力、Rayleigh 长度与频率。"""
from __future__ import annotations

import numpy as np

KB = 1.380649e-23
ATOM_MASS_U = 132.90545196
U_TO_KG = 1.66053906660e-27
WAVELENGTH_M = 1055e-9
WAIST_M = 1.17e-6
DEPTH_J = 280e-6 * KB


def physics_defaults(depth_j: float = DEPTH_J, waist_m: float = WAIST_M,
                      wavelength_m: float = WAVELENGTH_M) -> dict:
    """返回默认物理常数字典（质量、z_R、解析频率等）。"""
    mass_kg = ATOM_MASS_U * U_TO_KG
    z_r = np.pi * waist_m ** 2 / wavelength_m
    return {
        "mass_kg": mass_kg,
        "depth_j": depth_j,
        "depth_uK": depth_j / KB / 1e-6,
        "waist_m": waist_m,
        "wavelength_m": wavelength_m,
        "z_r_m": z_r,
        "omega_r": float(np.sqrt(4.0 * depth_j / (mass_kg * waist_m ** 2))),
        "omega_z": float(np.sqrt(2.0 * depth_j / (mass_kg * z_r ** 2))),
    }


def static_potential(xi, eta, z, depth_j, waist_m, z_r_m):
    """标准三维高斯势 U0 = -D/(1+(z/zR)^2) * exp(-2(xi^2+eta^2)/(w0^2(1+(z/zR)^2)))。"""
    den = 1.0 + (np.asarray(z) / z_r_m) ** 2
    return -depth_j / den * np.exp(-2.0 * (np.asarray(xi) ** 2 + np.asarray(eta) ** 2)
                                   / (waist_m ** 2 * den))


def static_force(xi, eta, z, depth_j, waist_m, z_r_m):
    """标准三维高斯势的解析力，返回 (F_xi, F_eta, F_z)。"""
    z = np.asarray(z, dtype=float)
    den = 1.0 + (z / z_r_m) ** 2
    prefac = depth_j / den * np.exp(-2.0 * (np.asarray(xi) ** 2 + np.asarray(eta) ** 2)
                                    / (waist_m ** 2 * den))
    f_xi = -prefac * 4.0 * np.asarray(xi) / (waist_m ** 2 * den)
    f_eta = -prefac * 4.0 * np.asarray(eta) / (waist_m ** 2 * den)
    # F_z = -D·E·z/(z_R²·A²)·[2 - 4ρ²/(w²A)]，A = den，E = exp 项
    rho2 = np.asarray(xi) ** 2 + np.asarray(eta) ** 2
    f_z = -prefac * z / z_r_m ** 2 / den * (2.0 - 4.0 * rho2 / (waist_m ** 2 * den))
    return f_xi, f_eta, f_z
