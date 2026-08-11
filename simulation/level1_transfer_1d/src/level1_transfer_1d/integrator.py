"""Velocity-Verlet for an explicitly time-dependent force F(x, t)."""
import numpy as np
def velocity_verlet_time_dependent(force_function, mass, x_initial, v_initial, times):
    time = np.asarray(times, dtype=float)
    if time.ndim != 1 or time.size < 2 or not np.all(np.isfinite(time)): raise ValueError("invalid times")
    steps = np.diff(time)
    if np.any(steps <= 0) or not np.allclose(steps, steps[0], rtol=1e-11, atol=max(1e-18, steps[0] * 1e-13)): raise ValueError("times must be uniformly increasing")
    if mass <= 0 or not np.isfinite(mass): raise ValueError("invalid mass")
    dt = float(steps[0]); x = np.empty_like(time); v = np.empty_like(time); a = np.empty_like(time)
    x[0], v[0] = x_initial, v_initial; a[0] = float(force_function(x[0], time[0])) / mass
    for i in range(time.size - 1):
        x[i + 1] = x[i] + v[i] * dt + 0.5 * a[i] * dt**2
        a[i + 1] = float(force_function(x[i + 1], time[i + 1])) / mass
        v[i + 1] = v[i] + 0.5 * (a[i] + a[i + 1]) * dt
    return time, x, v, a
