"""Exact Bluvstein cubic and Zhang S2; derivatives are one-sided at endpoints.

Cubic and gamma=1.5625 have acceleration jumps at stationary joins. Their
ordinary interior jerk does not include the resulting distributional impulses.
"""
import numpy as np
from numpy.polynomial import Polynomial


def gamma_family_coeffs(gamma):
    """c2,c3,c4,c5. Gamma is midpoint speed, not necessarily global peak."""
    if not np.isfinite(gamma):
        raise ValueError("gamma must be finite")
    return 15-8*gamma, -50+32*gamma, 60-40*gamma, -24+16*gamma


def polynomial(kind):
    if kind in ("constant_jerk", "rb_cubic"):
        return Polynomial([0, 0, 3, -2])
    if kind in ("low_peak_v", "zero_jerk"):
        gamma = 1.5625
    elif kind in ("min_jerk", "minimum_jerk"):
        gamma = 1.875
    elif kind.startswith("gamma:"):
        gamma = float(kind.split(":")[1])
    else:
        raise ValueError(kind)
    return Polynomial([0, 0, *gamma_family_coeffs(gamma)])


class Trajectory:
    def __init__(self, t_s, x_m, v_m_s, name, a_m_s2=None, j_m_s3=None):
        self.t_s, self.x_m, self.v_m_s, self.name = t_s, x_m, v_m_s, name
        self.a_m_s2, self.j_m_s3 = a_m_s2, j_m_s3

    @property
    def duration_s(self):
        return float(self.t_s[-1])


def poly_trajectory(kind, distance_m, duration_s, dt_s):
    if duration_s <= 0 or dt_s <= 0:
        raise ValueError("positive times required")
    t = np.linspace(0, duration_s, max(1, int(np.ceil(duration_s/dt_s)))+1)
    s = t/duration_s
    if kind == "linear":
        return Trajectory(t, distance_m*s, np.full_like(s, distance_m/duration_s), kind)
    if kind == "constant_accel":
        q = np.where(s < .5, 2*s*s, 1-2*(1-s)**2)
        v = np.where(s < .5, 4*s, 4*(1-s))
        a = np.where(s < .5, 4., -4.)
        return Trajectory(t, distance_m*q, distance_m/duration_s*v, kind,
                          distance_m/duration_s**2*a, np.zeros_like(s))
    p = polynomial(kind)
    return Trajectory(t, distance_m*p(s), distance_m/duration_s*p.deriv()(s), kind,
                      distance_m/duration_s**2*p.deriv(2)(s),
                      distance_m/duration_s**3*p.deriv(3)(s))


MIN_JERK = "gamma:1.875"
ZERO_JERK = "low_peak_v"
