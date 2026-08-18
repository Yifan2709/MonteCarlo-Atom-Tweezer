"""Two-phase merge/drop-off waveform."""
import numpy as np
def _out(value, original): return float(value) if np.asarray(original).ndim == 0 else value
def aod_center(t, config):
    original = t; t = np.clip(np.asarray(t, dtype=float), 0.0, config.transfer.duration_s)
    s = np.minimum(t / config.transfer.move_duration_s, 1.0); h = 3*s**2 - 2*s**3
    return _out(config.aod.initial_center_m * (1-h), original)
def aod_depth(t, config):
    original = t; t = np.clip(np.asarray(t, dtype=float), 0.0, config.transfer.duration_s)
    s = np.clip((t-config.transfer.move_duration_s)/(config.transfer.duration_s-config.transfer.move_duration_s), 0.0, 1.0)
    return _out(config.aod.depth_j * (1-s)**2, original)
