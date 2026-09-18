"""tweezer-experiment: offline experiment planning for tweezer transport.

Public API (stable semantics; see README for the mapping table):
    validate_experiment(...) -> ValidationReport
    compile_protocol(...)    -> CompiledProtocol
    simulate(...)            -> SimulationResult
    search_controls(...)     -> SearchResult
    export_controls(...)     -> dict (export manifest)

Importing this package has no side effects: no simulation starts, no
historical outputs are read, no working-directory or global RNG state is
touched.
"""
from .errors import ToolkitError
from .schemas import (ExperimentBundle, DeviceConfig, PreparationConfig, NoiseConfig,
                      ProtocolSpec, NumericsConfig, HardwareCalibration, SearchSpec,
                      validate_experiment, effective_config)
from .protocol_compile import compile_protocol, apply_override
from .simulation import simulate, save_result
from .search import search_controls
from .export import export_controls

__version__ = "1.0.0"

__all__ = ["ExperimentBundle", "DeviceConfig", "PreparationConfig", "NoiseConfig",
           "ProtocolSpec", "NumericsConfig", "HardwareCalibration", "SearchSpec",
           "validate_experiment", "effective_config", "compile_protocol", "apply_override",
           "simulate", "save_result", "search_controls", "export_controls",
           "ToolkitError", "__version__"]
