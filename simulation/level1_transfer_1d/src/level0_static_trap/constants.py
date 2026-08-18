"""物理常数和集中式 SI 单位换算。"""

from scipy.constants import Boltzmann, physical_constants

BOLTZMANN_CONSTANT = float(Boltzmann)
ATOMIC_MASS_CONSTANT = float(physical_constants["atomic mass constant"][0])
MICROKELVIN_TO_KELVIN = 1.0e-6
MICROMETER_TO_METER = 1.0e-6


def microkelvin_to_joule(value: float) -> float:
    """将以微开尔文表示的势阱深度换算为焦耳。"""

    return float(value) * MICROKELVIN_TO_KELVIN * BOLTZMANN_CONSTANT


def micrometer_to_meter(value: float) -> float:
    """将微米换算为米。"""

    return float(value) * MICROMETER_TO_METER


def atomic_mass_to_kg(value: float) -> float:
    """将原子质量单位换算为千克。"""

    return float(value) * ATOMIC_MASS_CONSTANT

