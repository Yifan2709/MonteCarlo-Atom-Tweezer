"""Unit conventions and conversions (A03 contract).

Canonical rule of this project: **configuration fields declare their unit in
the field name** (temperature_uK, waist_um, duration_us, psd_per_hz) and are
converted to strict SI at parse time; every value returned by the public API
is SI unless its name carries an explicit unit suffix. There is exactly one
conversion site per unit (this module) and every factor is tested against an
independently written expectation.
"""
from __future__ import annotations

import numpy as np

from .errors import ToolkitError, E_UNIT

K_B = 1.380649e-23            # J/K, exact (SI 2019)
ATOMIC_MASS_U_KG = 1.66053906660e-27
H_PLANCK = 6.62607015e-34     # J*s, exact

#: microkelvin expressed in joules of energy (E = k_B * T).
U_K_IN_J = K_B * 1e-6


def temperature_uK_to_joule(t_uK):
    return np.asarray(t_uK, dtype=float) * U_K_IN_J


def joule_to_temperature_uK(e_j):
    return np.asarray(e_j, dtype=float) / U_K_IN_J


def um_to_m(x):
    return np.asarray(x, dtype=float) * 1e-6


def m_to_um(x):
    return np.asarray(x, dtype=float) * 1e6


def us_to_s(t):
    return np.asarray(t, dtype=float) * 1e-6


def s_to_us(t):
    return np.asarray(t, dtype=float) * 1e6


def nm_to_m(x):
    return np.asarray(x, dtype=float) * 1e-9


def hz_to_rad_per_s(f):
    return 2.0 * float(f) * 3.141592653589793


def rad_per_s_to_hz(w):
    return float(w) / (2.0 * 3.141592653589793)


def mass_u_to_kg(m_u):
    return float(m_u) * ATOMIC_MASS_U_KG


#: Supported textual units for the few free-form "{value, unit}" slots
#: (hardware calibration tables). Anything else is E_UNIT.
_LENGTH_UNITS = {"m": 1.0, "um": 1e-6, "nm": 1e-9}
_FREQUENCY_UNITS = {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6}
_VOLTAGE_UNITS = {"V": 1.0, "mV": 1e-3}
_TIME_UNITS = {"s": 1.0, "us": 1e-6, "ms": 1e-3, "ns": 1e-9}
_UNIT_FAMILIES = {"length": _LENGTH_UNITS, "frequency": _FREQUENCY_UNITS,
                  "voltage": _VOLTAGE_UNITS, "time": _TIME_UNITS}


def to_si(value: float, unit: str, family: str) -> float:
    """Convert `value` in textual `unit` of a `family` to SI; else E_UNIT."""
    table = _UNIT_FAMILIES.get(family)
    if table is None or unit not in table:
        raise ToolkitError(E_UNIT, f"unsupported {family} unit '{unit}' "
                                   f"(supported: {sorted(table or {})})", field="unit")
    return float(value) * table[unit]
