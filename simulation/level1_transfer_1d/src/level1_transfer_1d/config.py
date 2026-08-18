"""Strict Level 1 configuration loading and SI-unit views."""
from __future__ import annotations
from dataclasses import asdict, dataclass, replace
from pathlib import Path
import math
from typing import Any, Mapping
import yaml
from level0_static_trap.constants import atomic_mass_to_kg, micrometer_to_meter, microkelvin_to_joule

@dataclass(frozen=True)
class ModelConfig:
    coordinate: str
    description: str
@dataclass(frozen=True)
class AtomConfig:
    species: str
    mass_u: float
    @property
    def mass_kg(self): return atomic_mass_to_kg(self.mass_u)
@dataclass(frozen=True)
class SlmConfig:
    center_um: float
    depth_uK: float
    waist_um: float
    depth_source: str
    waist_source: str
    @property
    def center_m(self): return micrometer_to_meter(self.center_um)
    @property
    def depth_j(self): return microkelvin_to_joule(self.depth_uK)
    @property
    def waist_m(self): return micrometer_to_meter(self.waist_um)
@dataclass(frozen=True)
class AodConfig:
    initial_center_um: float
    depth_uK: float
    waist_um: float
    depth_source: str
    waist_source: str
    waist_is_assumed: bool
    @property
    def initial_center_m(self): return micrometer_to_meter(self.initial_center_um)
    @property
    def depth_j(self): return microkelvin_to_joule(self.depth_uK)
    @property
    def waist_m(self): return micrometer_to_meter(self.waist_um)
@dataclass(frozen=True)
class TransferConfig:
    duration_us: float
    duration_source: str
    move_fraction: float
    ramp_fraction: float
    @property
    def duration_s(self): return self.duration_us * 1e-6
    @property
    def move_duration_s(self): return self.duration_s * self.move_fraction
@dataclass(frozen=True)
class ThermalConfig:
    temperature_uK: float
    source: str
@dataclass(frozen=True)
class TrajectoryConfig:
    dt_us: float
    @property
    def dt_s(self): return self.dt_us * 1e-6
@dataclass(frozen=True)
class MonteCarloConfig:
    shots: int
    seed: int
@dataclass(frozen=True)
class OutputConfig:
    directory: str
    resolved_directory: Path
@dataclass(frozen=True)
class SimulationConfig:
    model: ModelConfig
    atom: AtomConfig
    slm: SlmConfig
    aod: AodConfig
    transfer: TransferConfig
    thermal: ThermalConfig
    trajectory: TrajectoryConfig
    monte_carlo: MonteCarloConfig
    output: OutputConfig
    config_path: Path
    project_root: Path

_KEYS = {
    "model": {"coordinate", "description"}, "atom": {"species", "mass_u"},
    "slm": {"center_um", "depth_uK", "waist_um", "depth_source", "waist_source"},
    "aod": {"initial_center_um", "depth_uK", "waist_um", "depth_source", "waist_source", "waist_is_assumed"},
    "transfer": {"duration_us", "duration_source", "move_fraction", "ramp_fraction"},
    "thermal": {"temperature_uK", "source"}, "trajectory": {"dt_us"},
    "monte_carlo": {"shots", "seed"}, "output": {"directory"},
}
def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise ValueError(f"{label} must be a finite number")
    return float(value)
def _root(path: Path) -> Path:
    return next((p.resolve() for p in (path.parent, *path.parents) if (p / "pyproject.toml").is_file()), path.parent.resolve())
def load_config(path: str | Path) -> SimulationConfig:
    config_path = Path(path).expanduser().resolve()
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(raw, Mapping) or set(raw) != set(_KEYS): raise ValueError("configuration sections must exactly match schema")
    for section, keys in _KEYS.items():
        if not isinstance(raw[section], Mapping) or set(raw[section]) != keys: raise ValueError(f"invalid fields in {section}")
    n = lambda s, k: _number(raw[s][k], f"{s}.{k}")
    root = _root(config_path)
    c = SimulationConfig(
        ModelConfig(str(raw["model"]["coordinate"]), str(raw["model"]["description"])),
        AtomConfig(str(raw["atom"]["species"]), n("atom", "mass_u")),
        SlmConfig(n("slm", "center_um"), n("slm", "depth_uK"), n("slm", "waist_um"), str(raw["slm"]["depth_source"]), str(raw["slm"]["waist_source"])),
        AodConfig(n("aod", "initial_center_um"), n("aod", "depth_uK"), n("aod", "waist_um"), str(raw["aod"]["depth_source"]), str(raw["aod"]["waist_source"]), bool(raw["aod"]["waist_is_assumed"])),
        TransferConfig(n("transfer", "duration_us"), str(raw["transfer"]["duration_source"]), n("transfer", "move_fraction"), n("transfer", "ramp_fraction")),
        ThermalConfig(n("thermal", "temperature_uK"), str(raw["thermal"]["source"])),
        TrajectoryConfig(n("trajectory", "dt_us")), MonteCarloConfig(int(n("monte_carlo", "shots")), int(n("monte_carlo", "seed"))),
        OutputConfig(str(raw["output"]["directory"]), (root / str(raw["output"]["directory"])).resolve()), config_path, root)
    _validate(c)
    return c
def _validate(c):
    if c.model.coordinate != "radial_1d" or c.model.description != "classical_aod_to_slm_transfer": raise ValueError("unsupported model")
    values = (c.atom.mass_u, c.slm.depth_uK, c.slm.waist_um, c.aod.depth_uK, c.aod.waist_um, c.transfer.duration_us, c.thermal.temperature_uK, c.trajectory.dt_us, c.monte_carlo.shots)
    if any(v <= 0 for v in values): raise ValueError("physical scales and shots must be positive")
    if not math.isclose(c.transfer.move_fraction + c.transfer.ramp_fraction, 1.0, abs_tol=1e-12): raise ValueError("transfer fractions must sum to one")
    if not math.isclose(c.transfer.duration_us / c.trajectory.dt_us, round(c.transfer.duration_us / c.trajectory.dt_us), abs_tol=1e-9): raise ValueError("duration/dt must be integral")
    omega_max = math.sqrt(4 * (c.aod.depth_j + c.slm.depth_j) / (c.atom.mass_kg * min(c.aod.waist_m, c.slm.waist_m) ** 2))
    if omega_max * c.trajectory.dt_s >= 0.05: raise ValueError("transfer timestep fails omega_max * dt < 0.05")
def with_output_override(c, directory: str):
    p = Path(directory).expanduser(); resolved = (p if p.is_absolute() else c.project_root / p).resolve()
    return replace(c, output=OutputConfig(directory, resolved))
def config_as_dict(c):
    data = asdict(c); data["output"]["resolved_directory"] = str(c.output.resolved_directory); data["config_path"] = str(c.config_path); data["project_root"] = str(c.project_root)
    return data
