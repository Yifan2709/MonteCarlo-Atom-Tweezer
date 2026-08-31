"""Bound thermal initial states in the exact static AOD Gaussian."""
from dataclasses import dataclass
import math
from level0_static_trap.constants import BOLTZMANN_CONSTANT, MICROKELVIN_TO_KELVIN
from level0_static_trap.potentials import gaussian_potential
@dataclass(frozen=True)
class ThermalInitialState:
    x0_m: float
    v0_m_s: float
    aod_total_energy_J: float
    aod_excitation_energy_J: float
    aod_energy_over_depth: float
def _well_params(config, well):
    """按采样阱返回 (depth_j, waist_m, center_m)；默认 AOD，可选 SLM。"""
    if well == "aod":
        return config.aod.depth_j, config.aod.waist_m, config.aod.initial_center_m
    if well == "slm":
        return config.slm.depth_j, config.slm.waist_m, config.slm.center_m
    raise ValueError(f"未知采样阱 {well!r}（可选 'aod' / 'slm'）")
def thermal_scales(config, well="aod"):
    depth_j, waist_m, _ = _well_params(config, well)
    omega = math.sqrt(4*depth_j/(config.atom.mass_kg*waist_m**2)); kbt = BOLTZMANN_CONSTANT*config.thermal.temperature_uK*MICROKELVIN_TO_KELVIN
    return math.sqrt(kbt/(config.atom.mass_kg*omega**2)), math.sqrt(kbt/config.atom.mass_kg)
def sample_thermal_initial_state(config, rng, max_attempts=100000, return_attempts=False, well="aod"):
    """简谐提议 + 束缚拒绝采样；well 选择采样阱（'aod' 为 Level 1 原行为）。

    well='slm' 用于 Level 2C 的 SLM→AOD pick-up 方向：在 SLM 单阱
    （140 μK、1.17 μm、中心 0）中按 E_SLM(x0,v0)<0 拒绝。
    返回的 dataclass 字段名沿用 aod_*（历史原因），其数值含义为
    “采样阱”的总能量/激发能。
    """
    depth_j, waist_m, center_m = _well_params(config, well)
    sigma_x, sigma_v = thermal_scales(config, well=well)
    for attempts in range(1, max_attempts+1):
        x = center_m + rng.normal(0, sigma_x); v = rng.normal(0, sigma_v)
        energy = 0.5*config.atom.mass_kg*v**2 + float(gaussian_potential(x, depth_j, waist_m, center_m))
        if energy < 0:
            excitation = energy + depth_j
            state = ThermalInitialState(x, v, energy, excitation, energy/depth_j)
            return (state, attempts) if return_attempts else state
    raise RuntimeError("could not sample bound state")


@dataclass(frozen=True)
class ThermalInitialState2D:
    """二维横向热初态（Level 2C 条件触发的 2D 模型用）。"""
    x0_m: float
    y0_m: float
    vx0_m_s: float
    vy0_m_s: float
    total_energy_J: float
    excitation_energy_J: float
    energy_over_depth: float


def sample_thermal_initial_state_2d(config, rng, max_attempts=100000, well="slm"):
    """二维高斯阱中的热初态采样：两轴各自简谐提议 + 总能量 E(x,y,vx,vy)<0 拒绝。

    阱为圆对称高斯 U(x,y)=−D·exp(−2((x−cx)²+y²)/w²)；提议分布两轴独立同分布，
    温度取 config.thermal.temperature_uK。
    """
    depth_j, waist_m, center_m = _well_params(config, well)
    sigma_x, sigma_v = thermal_scales(config, well=well)
    for _attempt in range(1, max_attempts + 1):
        x = center_m + rng.normal(0, sigma_x)
        y = rng.normal(0, sigma_x)
        vx = rng.normal(0, sigma_v)
        vy = rng.normal(0, sigma_v)
        potential = -depth_j * math.exp(-2.0 * ((x - center_m) ** 2 + y ** 2) / waist_m ** 2)
        energy = 0.5 * config.atom.mass_kg * (vx ** 2 + vy ** 2) + potential
        if energy < 0:
            excitation = energy + depth_j
            return ThermalInitialState2D(x, y, vx, vy, energy, excitation, energy / depth_j)
    raise RuntimeError("could not sample bound 2D state")
