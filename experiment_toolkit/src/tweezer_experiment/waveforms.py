"""Analytic move-profile families with consistent derivatives (F2).

Every family maps normalized time u in [0,1] to normalized position q(u)
and velocity q'(u) (so physical velocity is q' * distance / duration).
Segment compilers use these both for controls and for exported velocity
columns, so the physical table always matches the simulated path.

``paper_manual`` reproduces the paper Fig. 6b single-cubic x/d = 3u^2-2u^3:
zero endpoint velocity but an acceleration jump at the motion endpoints.
That jump is intentional and reported by the compiler; callers must not
smooth it silently to pass hardware checks.

``constant_jerk`` is the four-segment +-J symmetric S-curve (zero endpoint
velocity AND acceleration), with normalized jerk J = 32 d / tau^3:
    t in [0, tau/4]:      v = J t^2 / 2,            x = J t^3 / 6
    t in [tau/4, tau/2]:  v = J tau^2/32 + J tau s/4 - J s^2/2,  s = t - tau/4
                          x = J tau^3/384 + J tau^2 s/32 + J tau s^2/8 - J s^3/6
    t > tau/2 mirrored:   x(t) = d - x(tau - t),  v(t) = v(tau - t)
"""
from __future__ import annotations

import numpy as np

from .errors import ToolkitError, E_MODEL_UNKNOWN, E_BAD_VALUE


def _paper_cubic(u):
    return 3.0 * u**2 - 2.0 * u**3, 6.0 * u * (1.0 - u)


def _minimum_jerk(u):
    return (10.0 * u**3 - 15.0 * u**4 + 6.0 * u**5,
            30.0 * u**2 - 60.0 * u**3 + 30.0 * u**4)


def _adiabatic_sine(u):
    return u - np.sin(2.0 * np.pi * u) / (2.0 * np.pi), 1.0 - np.cos(2.0 * np.pi * u)


def _linear(u):
    return u, np.ones_like(u)


def _constant_jerk(u):
    u = np.asarray(u, dtype=float)
    q = np.empty_like(u)
    dq = np.empty_like(u)
    first = u <= 0.25
    second = (u > 0.25) & (u <= 0.5)
    rest = u > 0.5
    t = u[first]
    q[first] = 16.0 * t**3 / 3.0          # J t^3/6 with J=32
    dq[first] = 16.0 * t**2               # J t^2/2
    s = u[second] - 0.25
    q[second] = 1.0 / 12.0 + s + 4.0 * s**2 - 16.0 * s**3 / 3.0
    dq[second] = 1.0 + 8.0 * s - 16.0 * s**2
    if rest.any():
        qr, vr = _constant_jerk(1.0 - u[rest])
        q[rest] = 1.0 - qr
        dq[rest] = vr
    return q, dq


FAMILIES = {
    "paper_manual": _paper_cubic,
    "minimum_jerk": _minimum_jerk,
    "constant_jerk": _constant_jerk,
    "adiabatic_sine": _adiabatic_sine,
    "linear": _linear,
    # smoothstep is numerically identical to the paper single cubic; kept as
    # an accepted alias so historical names parse, mapped to the same math.
    "smoothstep": _paper_cubic,
}

FAMILY_NOTES = {
    "paper_manual": "single cubic 3u^2-2u^3; endpoint acceleration jumps (paper Fig. 6b)",
    "constant_jerk": "four-segment +-J symmetric S-curve; zero endpoint acceleration",
}


def family_fn(name: str):
    fn = FAMILIES.get(name)
    if fn is None:
        raise ToolkitError(E_MODEL_UNKNOWN, f"unknown path family '{name}'", "path")
    return fn


def transfer_controls(family: str, ramp_fraction: float, duration_s: float,
                      from_m, to_m, depth_from_j: float, depth_to_j: float,
                      order: str = "depth_then_move"):
    """Two-phase transfer controls over local time in [0, duration].

    ``depth_then_move`` (pickup orientation): depth ramps quadratically in
    time over the first ``ramp_fraction`` of the duration while the center
    parks at ``from``; the center then moves along the straight line to
    ``to`` with the chosen family over the remainder.

    ``move_then_depth`` (dropoff orientation): the center moves first over
    the first ``1 - ramp_fraction`` of the duration, then the depth ramps
    quadratically while parked at ``to``. This is the exact time reverse of
    the pickup semantics at equal duration and ramp fraction.

    The quadratic depth-ramp shape is fixed by the paper pickup semantics
    and does not follow the move family.
    """
    fn = family_fn(family)
    ramp_s = ramp_fraction * duration_s
    move_s = duration_s - ramp_s
    if ramp_s <= 0.0 or move_s <= 0.0:
        raise ToolkitError(E_BAD_VALUE,
                           "ramp_fraction must be > 0 and leave positive time for both phases",
                           "ramp_fraction")
    fx, fy = from_m
    tx, ty = to_m

    def controls(local_t):
        t = np.asarray(local_t, dtype=float)
        if order == "depth_then_move":
            s = np.clip(t / ramp_s, 0.0, 1.0)
            depth = depth_from_j + (depth_to_j - depth_from_j) * s * s
            u = np.clip((t - ramp_s) / move_s, 0.0, 1.0)
            moving = (t >= ramp_s) & (t <= duration_s)
        else:
            # exact time mirror of the pickup ramp: quadratic with zero
            # slope at the END of the ramp window (paper dropoff = reversed
            # pickup gives D*(1-s)^2, reproduced by this parameterization)
            s = np.clip((t - move_s) / ramp_s, 0.0, 1.0)
            depth = depth_from_j + (depth_to_j - depth_from_j) * (1.0 - (1.0 - s) ** 2)
            u = np.clip(t / move_s, 0.0, 1.0)
            moving = (t >= 0.0) & (t <= move_s)
        q, dq = fn(u)
        x = fx + (tx - fx) * q
        y = fy + (ty - fy) * q
        vx = (tx - fx) * np.where(moving, dq, 0.0) / move_s
        vy = (ty - fy) * np.where(moving, dq, 0.0) / move_s
        return x, y, vx, vy, depth

    return controls
