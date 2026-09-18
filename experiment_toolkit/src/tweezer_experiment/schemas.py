"""Versioned input schemas and strict validation (F1).

Seven configuration objects: DeviceConfig, PreparationConfig, NoiseConfig,
ProtocolSpec, NumericsConfig, HardwareCalibration, SearchSpec. Parsing is
strict: unknown fields, non-finite numbers, wrong types, bad units, out-of-
domain values and inconsistent fields all raise ToolkitError with a stable
code and the dotted field path. Field names carry their unit suffix; values
are converted to SI here and stay SI everywhere else.

Source tags: any field may declare a companion ``<field>__source`` mapping
``{tag: measured|derived|assumed|fitted|synthetic, note: str}``. Tags feed
the capability classification: a simulation whose device parameters are all
*assumed* is runnable and idealized, never a calibrated prediction.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np
import yaml

from .errors import (ToolkitError, E_UNKNOWN_FIELD, E_MISSING_FIELD, E_BAD_TYPE,
                     E_NOT_FINITE, E_BAD_VALUE, E_MODEL_UNKNOWN, E_CONFLICT,
                     E_LENGTH_MISMATCH, E_CAL_RANGE)
from . import units as U
from .provenance import SOURCE_TAGS

SCHEMA_VERSION = "1"

SPECIES_PRESETS = {"Cs133": 132.90545196, "Rb87": 86.90918053, "Yb171": 170.9363313}
PATH_FAMILIES = ("paper_manual", "minimum_jerk", "constant_jerk", "adiabatic_sine",
                 "linear", "smoothstep")
NOISE_MODELS = ("none", "diffusion", "local_mean", "legacy_constant_power")
INTERPOLATIONS = ("linear", "cubic", "hold")


# ---------------------------------------------------------------- helpers --

def _mapping(value: Any, path: str) -> dict:
    if not isinstance(value, dict):
        raise ToolkitError(E_BAD_TYPE, f"expected mapping, got {type(value).__name__}", path)
    return value


def _sources_of(raw: dict, key: str, path: str) -> dict | None:
    """Parse an optional ``<key>__source`` companion field."""
    tag_key = f"{key}__source"
    if tag_key not in raw:
        return None
    src = _mapping(raw.pop(tag_key), f"{path}.{tag_key}")
    tag = src.get("tag")
    if tag not in SOURCE_TAGS:
        raise ToolkitError(E_BAD_VALUE, f"source tag must be one of {SOURCE_TAGS}", f"{path}.{tag_key}.tag")
    return {"tag": tag, "note": str(src.get("note", ""))}


def parse_fields(raw: dict, path: str, spec: dict[str, Callable[[Any, str], Any]],
                 optional: dict[str, Any] | None = None) -> dict:
    """Strictly map `raw` through `spec` (required) and `optional` (defaults).

    Unknown keys raise E_UNKNOWN_FIELD; `__source` companions are collected
    separately so they never collide with real fields.
    """
    raw = dict(raw)
    sources: dict[str, dict] = {}
    out: dict[str, Any] = {}
    known = set(spec) | set(optional or {}) | {f"{k}__source" for k in list(spec) + list(optional or {})}
    for key in list(raw):
        if key.endswith("__source"):
            base = key.removesuffix("__source")
            if base not in known:
                raise ToolkitError(E_UNKNOWN_FIELD, f"unknown field '{key}'", f"{path}.{key}")
            sources[base] = _sources_of(raw, base, path) or {}
            continue
    for key in raw:
        if key.endswith("__source"):
            continue
        if key not in spec and key not in (optional or {}):
            raise ToolkitError(E_UNKNOWN_FIELD, f"unknown field '{key}'", f"{path}.{key}")
    for key, fn in spec.items():
        if key not in raw:
            raise ToolkitError(E_MISSING_FIELD, f"missing required field '{key}'", f"{path}.{key}")
        out[key] = fn(raw[key], f"{path}.{key}")
    for key, default in (optional or {}).items():
        if key in raw:
            out[key] = raw[key] if default is None else _reparse(raw[key], default, f"{path}.{key}")
        else:
            out[key] = None if default is None else default
    if sources:
        out["__sources"] = sources
    return out


def _reparse(value: Any, default: Any, path: str) -> Any:
    """Re-parse an explicitly provided optional value against its default."""
    if isinstance(default, bool):
        return _bool(value, path)
    if isinstance(default, int) and not isinstance(default, bool):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ToolkitError(E_BAD_TYPE, f"expected integer, got {type(value).__name__}", path)
        number = _number(value, path)
        if number != int(number):
            raise ToolkitError(E_BAD_VALUE, f"expected integer, got {value}", path)
        return int(number)
    if isinstance(default, float):
        return _number(value, path)
    if isinstance(default, str):
        return _str(value, path)
    if isinstance(default, list):
        return _vector(value, path)
    if isinstance(default, dict):
        return _mapping(value, path)
    return value


def _number(value: Any, path: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ToolkitError(E_BAD_TYPE, f"expected number, got {type(value).__name__}", path)
    v = float(value)
    if not np.isfinite(v):
        raise ToolkitError(E_NOT_FINITE, "value must be finite", path)
    return v


def _positive(value: Any, path: str) -> float:
    v = _number(value, path)
    if v <= 0:
        raise ToolkitError(E_BAD_VALUE, f"value must be > 0, got {v}", path)
    return v


def _nonneg(value: Any, path: str) -> float:
    v = _number(value, path)
    if v < 0:
        raise ToolkitError(E_BAD_VALUE, f"value must be >= 0, got {v}", path)
    return v


def _bool(value: Any, path: str) -> bool:
    if not isinstance(value, bool):
        raise ToolkitError(E_BAD_TYPE, f"expected bool, got {type(value).__name__}", path)
    return value


def _str(value: Any, path: str) -> str:
    if not isinstance(value, str) or not value:
        raise ToolkitError(E_BAD_TYPE, f"expected nonempty string", path)
    return value


def _enum(options: tuple[str, ...]) -> Callable:
    def parse(value, path):
        s = _str(value, path)
        if s not in options:
            raise ToolkitError(E_MODEL_UNKNOWN, f"'{s}' not in {list(options)}", path)
        return s
    return parse


def _vector(value: Any, path: str, length: int | None = None) -> list[float]:
    if not isinstance(value, (list, tuple)):
        raise ToolkitError(E_BAD_TYPE, f"expected list of numbers", path)
    if length is not None and len(value) != length:
        raise ToolkitError(E_LENGTH_MISMATCH, f"expected length {length}, got {len(value)}", path)
    return [_number(v, f"{path}[{i}]") for i, v in enumerate(value)]


def _vector2(value, path):
    v = _vector(value, path, 2)
    return v


# ------------------------------------------------------------------ device --

@dataclass(frozen=True)
class BeamConfig:
    wavelength_m: float
    waist_m: float
    nominal_depth_j: float
    center_m: tuple[float, float, float]


@dataclass(frozen=True)
class DeviceConfig:
    mass_kg: float
    species: str
    slm: BeamConfig
    aod: BeamConfig
    lensing_mode: str          # "off" | "velocity_shift"
    lensing_vs_m_s: float | None
    sources: dict = field(default_factory=dict, compare=False)

    @staticmethod
    def parse(raw: Any) -> "DeviceConfig":
        raw = _mapping(raw, "device")

        beam_sources: dict = {}

        def beam(key: str) -> BeamConfig:
            b = _mapping(raw[key], f"device.{key}")
            parsed = parse_fields(b, f"device.{key}", {
                "wavelength_nm": lambda v, p: U.nm_to_m(_positive(v, p)),
                "waist_um": lambda v, p: U.um_to_m(_positive(v, p)),
                "depth_uK": lambda v, p: U.temperature_uK_to_joule(_positive(v, p)),
            }, {"center_um": None})
            parsed.setdefault("__sources", {})
            center = parsed.get("center_um") or [0.0, 0.0, 0.0]
            center = _vector(center, f"device.{key}.center_um", 3)
            for src_key, src_val in parsed["__sources"].items():
                beam_sources[f"{key}.{src_key}"] = src_val
            return BeamConfig(parsed["wavelength_nm"], parsed["waist_um"],
                              parsed["depth_uK"], tuple(U.um_to_m(c) for c in center))

        head_allowed = {"species", "mass_u", "slm", "aod", "lensing"}
        for key in raw:
            if key.endswith("__source"):
                continue
            if key not in head_allowed:
                raise ToolkitError(E_UNKNOWN_FIELD, f"unknown field '{key}'", f"device.{key}")
        if "species" not in raw:
            raise ToolkitError(E_MISSING_FIELD, "missing required field 'species'", "device.species")
        head = {"species": _str(raw["species"], "device.species")}
        if head["species"] not in SPECIES_PRESETS and "mass_u" not in raw:
            raise ToolkitError(E_BAD_VALUE,
                               f"unknown species {head['species']!r}; presets: {list(SPECIES_PRESETS)} "
                               f"or provide device.mass_u", "device.species")
        mass = (U.mass_u_to_kg(_positive(raw["mass_u"], "device.mass_u"))
                if "mass_u" in raw else U.mass_u_to_kg(SPECIES_PRESETS[head["species"]]))
        lens = dict(_mapping(raw.get("lensing", {"mode": "off"}), "device.lensing"))
        # YAML 1.1 parses bare `off` as boolean False; coerce before parsing.
        if lens.get("mode") is False:
            lens["mode"] = "off"
        if lens.get("mode") is True:
            raise ToolkitError(E_BAD_VALUE, "lensing.mode: use 'velocity_shift', not 'on'",
                               "device.lensing.mode")
        lensing = parse_fields(lens, "device.lensing", {}, {
            "mode": "off", "vs_m_s": None})
        mode = _enum(("off", "velocity_shift"))(lensing["mode"], "device.lensing.mode")
        vs: float | None = None
        if mode == "velocity_shift":
            vs = _positive(lensing["vs_m_s"], "device.lensing.vs_m_s") if lensing["vs_m_s"] is not None else None
            if vs is None:
                raise ToolkitError(E_MISSING_FIELD,
                                   "velocity_shift lensing requires vs_m_s (scan speed)", "device.lensing.vs_m_s")
        return DeviceConfig(mass, head["species"], beam("slm"), beam("aod"),
                            mode, vs, beam_sources)


# ------------------------------------------------------------- preparation --

@dataclass(frozen=True)
class PreparationConfig:
    temperature_uK: float
    sampler: str
    site_disorder_enabled: bool
    slm_depth_rel_sigma: float
    aod_depth_rel_sigma: float
    aod_waist_rel_sigma: float
    alignment_sigma_m: float
    shot_jitter_sigma_m: float
    sources: dict = field(default_factory=dict, compare=False)

    @staticmethod
    def parse(raw: Any) -> "PreparationConfig":
        raw = _mapping(raw, "preparation")
        p = parse_fields(raw, "preparation", {
            "temperature_uK": _positive,
        }, {
            "sampler": "site_bound_harmonic",
            "site_disorder": True,
            "slm_depth_rel_sigma": 0.0,
            "aod_depth_rel_sigma": 0.0,
            "aod_waist_rel_sigma": 0.0,
            "alignment_sigma_um": 0.0,
            "shot_jitter_alignment_sigma_nm": 0.0,
        })
        disorder = p["site_disorder"]
        sigmas = {k: _nonneg(p[k], f"preparation.{k}") for k in
                  ("slm_depth_rel_sigma", "aod_depth_rel_sigma", "aod_waist_rel_sigma")}
        return PreparationConfig(
            p["temperature_uK"], _enum(("site_bound_harmonic",))(p["sampler"], "preparation.sampler"),
            _bool(disorder, "preparation.site_disorder"),
            sigmas["slm_depth_rel_sigma"], sigmas["aod_depth_rel_sigma"],
            sigmas["aod_waist_rel_sigma"],
            U.um_to_m(_nonneg(p["alignment_sigma_um"], "preparation.alignment_sigma_um")),
            U.nm_to_m(_nonneg(p["shot_jitter_alignment_sigma_nm"],
                              "preparation.shot_jitter_alignment_sigma_nm")),
            p.get("__sources", {}))


# ------------------------------------------------------------------- noise --

@dataclass(frozen=True)
class NoiseConfig:
    model: str
    aod_rin_psd_per_hz: float
    aod_pointing_psd_m2_per_hz: float
    slm_rin_psd_per_hz: float
    slm_pointing_psd_m2_per_hz: float
    recoil_per_joule: float
    legacy_power_w: float | None
    sources: dict = field(default_factory=dict, compare=False)

    @staticmethod
    def parse(raw: Any) -> "NoiseConfig":
        raw = _mapping(raw, "noise")
        p = parse_fields(raw, "noise", {
            "model": _enum(NOISE_MODELS),
        }, {
            "aod_rin_psd_per_hz": 0.0, "aod_pointing_psd_m2_per_hz": 0.0,
            "slm_rin_psd_per_hz": 0.0, "slm_pointing_psd_m2_per_hz": 0.0,
            "recoil_per_joule": 0.0, "legacy_power_w": None,
        })
        for key in ("aod_rin_psd_per_hz", "aod_pointing_psd_m2_per_hz",
                    "slm_rin_psd_per_hz", "slm_pointing_psd_m2_per_hz",
                    "recoil_per_joule"):
            p[key] = _nonneg(p[key], f"noise.{key}")
        legacy = p["legacy_power_w"]
        if p["model"] == "legacy_constant_power":
            if legacy is None:
                raise ToolkitError(E_MISSING_FIELD,
                                   "legacy_constant_power model requires legacy_power_w", "noise.legacy_power_w")
            legacy = _nonneg(legacy, "noise.legacy_power_w")
        elif legacy is not None:
            raise ToolkitError(E_CONFLICT,
                               "legacy_power_w is only meaningful with model=legacy_constant_power",
                               "noise.legacy_power_w")
        return NoiseConfig(p["model"], p["aod_rin_psd_per_hz"], p["aod_pointing_psd_m2_per_hz"],
                           p["slm_rin_psd_per_hz"], p["slm_pointing_psd_m2_per_hz"],
                           p["recoil_per_joule"], legacy, p.get("__sources", {}))


# ---------------------------------------------------------------- protocol --

@dataclass(frozen=True)
class Segment:
    id: str
    kind: str                    # static | transfer | csv
    duration_s: float | None     # None for csv (inferred from table)
    static: dict | None          # {center_m: (x,y), depth_j, slm_factor}
    transfer: dict | None        # {from_m, to_m, depth_from_j, depth_to_j, family, ramp_fraction, slm_factor}
    csv: dict | None             # {file, interpolation}
    observe: str | None          # None | "aod" | "slm"


@dataclass(frozen=True)
class SequenceItem:
    """Either one segment id or a repeated group of items."""
    segment: str | None = None
    repeat: int | None = None
    items: tuple["SequenceItem", ...] = ()


@dataclass(frozen=True)
class ProtocolSpec:
    segments: dict[str, Segment]
    sequence: tuple[SequenceItem, ...]
    final_hold_s: float | None
    final_observe: str | None
    sources: dict = field(default_factory=dict, compare=False)

    @staticmethod
    def parse(raw: Any) -> "ProtocolSpec":
        raw = _mapping(raw, "protocol")
        p = parse_fields(raw, "protocol", {
            "segments": lambda v, path: v,
            "sequence": lambda v, path: v,
        }, {"final_hold_us": None, "final_observe": None})
        segs_raw = p["segments"]
        if not isinstance(segs_raw, list) or not segs_raw:
            raise ToolkitError(E_BAD_TYPE, "protocol.segments must be a nonempty list", "protocol.segments")
        segments: dict[str, Segment] = {}
        for i, s in enumerate(segs_raw):
            seg = _parse_segment(s, f"protocol.segments[{i}]")
            if seg.id in segments:
                raise ToolkitError(E_CONFLICT, f"duplicate segment id '{seg.id}'", f"protocol.segments[{i}].id")
            segments[seg.id] = seg

        def parse_item(item: Any, path: str) -> SequenceItem:
            if isinstance(item, str):
                if item not in segments:
                    raise ToolkitError(E_BAD_VALUE, f"unknown segment id '{item}'", path)
                return SequenceItem(segment=item)
            m = _mapping(item, path)
            if "repeat" in m:
                count = m["repeat"]
                if isinstance(count, bool) or not isinstance(count, int) or count < 1:
                    raise ToolkitError(E_BAD_VALUE, "repeat must be a positive integer", f"{path}.repeat")
                inner = m.get("items")
                if not isinstance(inner, list) or not inner:
                    raise ToolkitError(E_MISSING_FIELD, "repeat requires items", f"{path}.items")
                return SequenceItem(repeat=count,
                                    items=tuple(parse_item(x, f"{path}.items[{j}]") for j, x in enumerate(inner)))
            raise ToolkitError(E_BAD_TYPE, "sequence items are segment ids or {repeat, items}", path)

        seq_raw = p["sequence"]
        if not isinstance(seq_raw, list) or not seq_raw:
            raise ToolkitError(E_BAD_TYPE, "protocol.sequence must be a nonempty list", "protocol.sequence")
        sequence = tuple(parse_item(x, f"protocol.sequence[{j}]") for j, x in enumerate(seq_raw))
        referenced = _collect_ids(sequence)
        unused = set(segments) - referenced
        if unused:
            raise ToolkitError(E_CONFLICT, f"segments defined but never used: {sorted(unused)}", "protocol.sequence")
        final_hold = U.us_to_s(_nonneg(p["final_hold_us"], "protocol.final_hold_us")) if p["final_hold_us"] is not None else None
        if p["final_observe"] is not None:
            _enum(("aod", "slm"))(p["final_observe"], "protocol.final_observe")
        return ProtocolSpec(segments, sequence, final_hold, p["final_observe"],
                            p.get("__sources", {}))


def _collect_ids(items: tuple[SequenceItem, ...]) -> set[str]:
    out: set[str] = set()
    for it in items:
        if it.segment:
            out.add(it.segment)
        else:
            out |= _collect_ids(it.items)
    return out


_SEGMENT_BODY_KEYS = {
    "static": {"duration_us", "aod_center_um", "aod_depth_uK", "slm_factor"},
    "transfer": {"duration_us", "from_um", "to_um", "depth_from_uK", "depth_to_uK",
                 "path", "ramp_fraction", "slm_factor", "order"},
    "csv": {"file", "interpolation", "velocity_columns"},
}


def _parse_segment(s: Any, path: str) -> Segment:
    s = _mapping(s, path)
    for key in s:
        if key.endswith("__source"):
            continue
        if key not in ("id", "type", "observe") and (
                "type" not in s or key not in _SEGMENT_BODY_KEYS.get(s["type"], set())):
            raise ToolkitError(E_UNKNOWN_FIELD, f"unknown field '{key}' for this segment type",
                               f"{path}.{key}")
    if "id" not in s or "type" not in s:
        raise ToolkitError(E_MISSING_FIELD, "segment needs 'id' and 'type'", path)
    p = {"id": _str(s["id"], f"{path}.id"),
         "type": _enum(("static", "transfer", "csv"))(s["type"], f"{path}.type"),
         "observe": s.get("observe")}
    observe = p["observe"]
    if observe is not None:
        _enum(("aod", "slm"))(observe, f"{path}.observe")
    body_source = {k: v for k, v in s.items() if k not in ("id", "type", "observe")}
    kind = p["type"]
    if kind == "static":
        body = parse_fields(body_source, path, {
            "duration_us": lambda v, pp: U.us_to_s(_positive(v, pp)),
            "aod_center_um": _vector2,
            "aod_depth_uK": _nonneg,
            "slm_factor": lambda v, pp: _in_range(v, pp, 0.0, 1.0),
        })
        return Segment(p["id"], kind, body["duration_us"],
                       {"center_m": tuple(U.um_to_m(c) for c in body["aod_center_um"]),
                        "depth_j": U.temperature_uK_to_joule(body["aod_depth_uK"]),
                        "slm_factor": body["slm_factor"]}, None, None, observe)
    if kind == "transfer":
        body = parse_fields(body_source, path, {
            "duration_us": lambda v, pp: U.us_to_s(_positive(v, pp)),
            "from_um": _vector2, "to_um": _vector2,
            "depth_from_uK": _nonneg, "depth_to_uK": _nonneg,
        }, {"path": None, "ramp_fraction": 0.48, "slm_factor": 1.0, "order": "depth_then_move"})
        family = body["path"] or "paper_manual"
        if isinstance(family, str):
            _enum(PATH_FAMILIES)(family, f"{path}.path")
        else:  # {"family": ..., optional override}
            fam = _mapping(family, f"{path}.path")
            family = _enum(PATH_FAMILIES)(fam.get("family"), f"{path}.path.family")
        ramp = _in_range(body["ramp_fraction"], f"{path}.ramp_fraction", 0.0, 0.95)
        order = _enum(("depth_then_move", "move_then_depth"))(body["order"], f"{path}.order")
        return Segment(p["id"], kind, body["duration_us"], None,
                       {"from_m": tuple(U.um_to_m(c) for c in body["from_um"]),
                        "to_m": tuple(U.um_to_m(c) for c in body["to_um"]),
                        "depth_from_j": U.temperature_uK_to_joule(body["depth_from_uK"]),
                        "depth_to_j": U.temperature_uK_to_joule(body["depth_to_uK"]),
                        "family": family, "ramp_fraction": ramp,
                        "order": order,
                        "slm_factor": _in_range(body["slm_factor"], f"{path}.slm_factor", 0.0, 1.0)},
                       None, observe)
    body = parse_fields(body_source, path, {
        "file": _str,
    }, {"interpolation": "linear", "velocity_columns": False})
    _enum(INTERPOLATIONS)(body["interpolation"], f"{path}.interpolation")
    return Segment(p["id"], kind, None, None, None,
                   {"file": body["file"], "interpolation": body["interpolation"],
                    "velocity_columns": _bool(body["velocity_columns"], f"{path}.velocity_columns")},
                   observe)


def _in_range(value: Any, path: str, lo: float, hi: float) -> float:
    v = _number(value, path)
    if not (lo <= v <= hi):
        raise ToolkitError(E_BAD_VALUE, f"value must be in [{lo}, {hi}], got {v}", path)
    return v


# ---------------------------------------------------------------- numerics --

@dataclass(frozen=True)
class NumericsConfig:
    dt_s: float
    noise_dt_s: float
    shots: int
    seed: int
    backend: str
    save_trajectories: str
    sources: dict = field(default_factory=dict, compare=False)

    @staticmethod
    def parse(raw: Any) -> "NumericsConfig":
        raw = _mapping(raw, "numerics")
        p = parse_fields(raw, "numerics", {
            "dt_us": lambda v, pp: U.us_to_s(_positive(v, pp)),
            "shots": lambda v, pp: _shots(v, pp),
            "seed": lambda v, pp: _seed(v, pp),
        }, {
            "noise_dt_us": 0.05, "backend": "numpy", "save_trajectories": "checkpoints",
        })
        backend = _enum(("numpy",))(p["backend"], "numerics.backend")
        save = _enum(("none", "checkpoints", "all"))(p["save_trajectories"], "numerics.save_trajectories")
        noise_dt = U.us_to_s(_positive(p["noise_dt_us"], "numerics.noise_dt_us"))
        return NumericsConfig(p["dt_us"], noise_dt, p["shots"], p["seed"], backend, save,
                              p.get("__sources", {}))


def _shots(value: Any, path: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ToolkitError(E_BAD_VALUE, "shots must be a positive integer", path)
    return value


def _seed(value: Any, path: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not (0 <= value < 2**63):
        raise ToolkitError(E_BAD_VALUE, "seed must be an integer in [0, 2^63)", path)
    return value


# ------------------------------------------------------------- calibration --

@dataclass(frozen=True)
class CalibrationChannel:
    kind: str                     # frequency | voltage
    unit: str                     # output unit (converted to SI)
    input_unit: str               # "um" (position) or "uK" (depth)
    input_axis: tuple[float, ...]   # strictly increasing, in input_unit
    output_si: tuple[float, ...]    # strictly increasing (Hz or V)
    valid_input: tuple[float, float]


@dataclass(frozen=True)
class HardwareCalibration:
    version: str
    source_tag: str
    note: str
    clock_hz: float
    quantization_bits: int | None
    channels: dict[str, CalibrationChannel]
    limits: dict
    sources: dict = field(default_factory=dict, compare=False)

    @staticmethod
    def parse(raw: Any) -> "HardwareCalibration":
        raw = _mapping(raw, "calibration")
        p = parse_fields(raw, "calibration", {
            "version": _str, "clock_hz": _positive,
        }, {"source_tag": "synthetic", "note": "", "quantization_bits": None,
            "channels": None, "limits": None})
        if p["source_tag"] not in SOURCE_TAGS:
            raise ToolkitError(E_BAD_VALUE, f"source_tag must be one of {SOURCE_TAGS}", "calibration.source_tag")
        chans_raw = p["channels"] or {}
        channels: dict[str, CalibrationChannel] = {}
        for name, c in _mapping(chans_raw, "calibration.channels").items():
            c = _mapping(c, f"calibration.channels.{name}")
            cp = parse_fields(c, f"calibration.channels.{name}", {
                "kind": _enum(("frequency", "voltage")),
                "unit": _str,
                "input_unit": _enum(("um", "uK")),
                "input": lambda v, pp: _vector(v, pp),
                "output": lambda v, pp: _vector(v, pp),
            }, {"valid_input": None})
            input_axis = tuple(cp["input"])
            output = tuple(cp["output"])
            if len(input_axis) != len(output) or len(input_axis) < 2:
                raise ToolkitError(E_LENGTH_MISMATCH, "input/output tables need equal length >= 2",
                                   f"calibration.channels.{name}")
            out_si = [U.to_si(v, cp["unit"], cp["kind"]) for v in output]
            for arr, label in ((input_axis, "input"), (out_si, "output")):
                if any(b <= a for a, b in zip(arr, arr[1:])):
                    raise ToolkitError(E_CAL_RANGE, f"{label} must be strictly increasing",
                                       f"calibration.channels.{name}.{label}")
            valid = tuple(cp["valid_input"]) if cp["valid_input"] else (input_axis[0], input_axis[-1])
            valid = tuple(_vector(list(valid), f"calibration.channels.{name}.valid_input", 2))
            if not (valid[0] < valid[1]):
                raise ToolkitError(E_CAL_RANGE, "valid range must satisfy min < max",
                                   f"calibration.channels.{name}.valid_input")
            channels[name] = CalibrationChannel(cp["kind"], cp["unit"], cp["input_unit"],
                                                input_axis, tuple(out_si), valid)
        limits = p["limits"] or {}
        limits = {k: _number(v, f"calibration.limits.{k}") for k, v in _mapping(limits, "calibration.limits").items()}
        bits = p["quantization_bits"]
        if bits is not None and (isinstance(bits, bool) or not isinstance(bits, int) or not 1 <= bits <= 32):
            raise ToolkitError(E_BAD_VALUE, "quantization_bits must be int in [1,32]", "calibration.quantization_bits")
        return HardwareCalibration(p["version"], p["source_tag"], p["note"], p["clock_hz"],
                                   bits, channels, limits, p.get("__sources", {}))


# ------------------------------------------------------------------ search --

@dataclass(frozen=True)
class SearchVariable:
    target: str        # e.g. "segment.pickup.duration_us"
    values: tuple[float, ...] | tuple[str, ...]
    kind: str          # "number" | "choice"


@dataclass(frozen=True)
class SearchSpec:
    variables: tuple[SearchVariable, ...]
    objective: dict
    screening_shots: int
    screening_seed: int
    validation_shots: int
    validation_seed: int
    max_candidates: int
    sources: dict = field(default_factory=dict, compare=False)

    @staticmethod
    def parse(raw: Any) -> "SearchSpec":
        raw = _mapping(raw, "search")
        p = parse_fields(raw, "search", {
            "variables": lambda v, pp: v,
            "objective": lambda v, pp: _mapping(v, pp),
        }, {
            "screening_shots": 256, "screening_seed": None,
            "validation_shots": 1024, "validation_seed": None,
            "max_candidates": 64,
        })
        vars_raw = p["variables"]
        if not isinstance(vars_raw, list) or not vars_raw:
            raise ToolkitError(E_BAD_TYPE, "search.variables must be a nonempty list", "search.variables")
        variables = []
        for i, v in enumerate(vars_raw):
            v = _mapping(v, f"search.variables[{i}]")
            vp = parse_fields(v, f"search.variables[{i}]", {"target": _str}, {"values": None, "grid": None})
            target = vp["target"]
            if not target.startswith("segment.") or target.count(".") != 2:
                raise ToolkitError(E_BAD_VALUE,
                                   "target must be 'segment.<id>.<field>'; supported fields: "
                                   "duration_us, ramp_fraction, depth_to_uK, path", f"search.variables[{i}].target")
            fieldname = target.split(".")[2]
            if fieldname not in ("duration_us", "ramp_fraction", "depth_to_uK", "path"):
                raise ToolkitError(E_BAD_VALUE, f"unsupported search field '{fieldname}'", f"search.variables[{i}].target")
            values: tuple
            if vp["values"] is not None:
                vals = vp["values"]
                if fieldname == "path":
                    values = tuple(_enum(PATH_FAMILIES)(x, f"search.variables[{i}].values") for x in vals)
                    variables.append(SearchVariable(target, values, "choice"))
                    continue
                values = tuple(_in_range(x, f"search.variables[{i}].values", 0.0, float("inf")) if fieldname != "ramp_fraction"
                               else _in_range(x, f"search.variables[{i}].values", 0.0, 0.95) for x in vals)
                variables.append(SearchVariable(target, values, "number"))
            elif vp["grid"] is not None:
                g = _vector(vp["grid"], f"search.variables[{i}].grid", 3)
                lo, hi, n = g
                if fieldname == "path":
                    raise ToolkitError(E_BAD_VALUE, "grid not supported for path choices", f"search.variables[{i}].grid")
                if not (0 <= lo <= hi) or n < 2 or n != int(n):
                    raise ToolkitError(E_BAD_VALUE, "grid must be [lo, hi, n], 0<=lo<=hi, n>=2", f"search.variables[{i}].grid")
                values = tuple(np.linspace(lo, hi, int(n)).tolist())
                variables.append(SearchVariable(target, values, "number"))
            else:
                raise ToolkitError(E_MISSING_FIELD, "variable needs values or grid", f"search.variables[{i}]")
        seen_targets = set()
        for v in variables:
            if v.target in seen_targets:
                raise ToolkitError(E_CONFLICT, f"duplicate variable target '{v.target}' "
                                               f"(later overrides would silently win)",
                                   f"search.variables")
            seen_targets.add(v.target)
        objective = p["objective"]
        otype = _enum(("maximize_capture", "minimize_duration"))(objective.get("type"), "search.objective.type")
        if otype == "maximize_capture":
            _enum(("final",))(objective.get("at", "final"), "search.objective.at")
        else:
            _str(objective.get("segment"), "search.objective.segment")
            _in_range(objective.get("capture_lower_bound"), "search.objective.capture_lower_bound", 0.0, 1.0)
        for key in ("screening_shots", "validation_shots"):
            p[key] = _shots(p[key], f"search.{key}")
        p["max_candidates"] = _shots(p["max_candidates"], "search.max_candidates")
        total = 1
        for v in variables:
            total *= len(v.values)
        if total > p["max_candidates"]:
            raise ToolkitError(E_BAD_VALUE,
                               f"grid expands to {total} candidates > max_candidates {p['max_candidates']}",
                               "search.max_candidates")
        screening_seed = p["screening_seed"] if p["screening_seed"] is not None else 424242
        validation_seed = p["validation_seed"] if p["validation_seed"] is not None else 424243
        return SearchSpec(tuple(variables), objective, p["screening_shots"], _seed(screening_seed, "search.screening_seed"),
                          p["validation_shots"], _seed(validation_seed, "search.validation_seed"),
                          p["max_candidates"], p.get("__sources", {}))


# --------------------------------------------------------------- experiment --

@dataclass(frozen=True)
class ExperimentBundle:
    device: DeviceConfig
    preparation: PreparationConfig
    noise: NoiseConfig
    protocol: ProtocolSpec
    numerics: NumericsConfig
    calibration: HardwareCalibration | None
    search: SearchSpec | None

    @staticmethod
    def parse(raw: Any) -> "ExperimentBundle":
        raw = _mapping(raw, "experiment")
        for key in ("device", "preparation", "noise", "protocol", "numerics"):
            if key not in raw:
                raise ToolkitError(E_MISSING_FIELD, f"missing required section '{key}'", key)
        known = {"device", "preparation", "noise", "protocol", "numerics", "calibration", "search", "schema_version"}
        for key in raw:
            if key not in known:
                raise ToolkitError(E_UNKNOWN_FIELD, f"unknown section '{key}'", key)
        if "schema_version" in raw and str(raw["schema_version"]) != SCHEMA_VERSION:
            raise ToolkitError(E_CONFLICT, f"schema_version {raw['schema_version']} != supported {SCHEMA_VERSION}",
                               "schema_version")
        return ExperimentBundle(
            DeviceConfig.parse(raw["device"]),
            PreparationConfig.parse(raw["preparation"]),
            NoiseConfig.parse(raw["noise"]),
            ProtocolSpec.parse(raw["protocol"]),
            NumericsConfig.parse(raw["numerics"]),
            HardwareCalibration.parse(raw["calibration"]) if "calibration" in raw else None,
            SearchSpec.parse(raw["search"]) if "search" in raw else None)

    @staticmethod
    def load(path: str | Path) -> "ExperimentBundle":
        text = Path(path).read_text(encoding="utf-8")
        try:
            raw = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            raise ToolkitError(E_BAD_TYPE, f"YAML parse failure: {exc}", str(path)) from exc
        return ExperimentBundle.parse(raw)


# --------------------------------------------------------- validation (F1) --

#: Experimental inputs that cannot be supplied by offline software. Used for
#: the capability classification: missing entries keep a run "idealized".
REQUIRED_MEASUREMENTS = (
    "preparation.temperature_uK", "device.slm.waist_um", "device.slm.depth_uK",
    "device.aod.waist_um", "device.aod.depth_uK", "noise.aod_rin_psd_per_hz",
    "noise.aod_pointing_psd_m2_per_hz", "calibration.channels.aod_x",
    "calibration.channels.aod_depth", "calibration.clock_hz",
)


@dataclass
class ValidationReport:
    ok: bool
    errors: list[str]
    warnings: list[str]
    capability: dict
    effective_config: dict

    def __bool__(self) -> bool:
        return self.ok


def validate_experiment(experiment: dict | ExperimentBundle,
                        calibration: HardwareCalibration | None = None,
                        search: SearchSpec | None = None) -> ValidationReport:
    """Validate an experiment mapping (or bundle) and classify capabilities.

    Returns ok=False with error strings (never raises) so CLI callers get a
    complete field-level report in one pass; parse errors carry codes/paths.
    """
    warnings: list[str] = []
    try:
        if isinstance(experiment, ExperimentBundle):
            bundle = experiment
        else:
            raw = dict(experiment)
            if calibration is not None:
                raw["calibration"] = calibration
            if search is not None:
                raw["search"] = search
            bundle = ExperimentBundle.parse(raw)
    except ToolkitError as exc:
        return ValidationReport(False, [f"{exc}"], warnings, {}, {})

    errors: list[str] = []
    cap = _capability(bundle, warnings)

    # Cross-field checks that need the compiled view are deferred; structural
    # ones live here so users get fast feedback without a compile.
    if bundle.numerics.noise_dt_s < bundle.numerics.dt_s:
        errors.append("numerics.noise_dt_us must be >= dt_us (multiple of it)")
    if bundle.noise.model == "none" and (bundle.noise.aod_rin_psd_per_hz or bundle.noise.recoil_per_joule
                                         or bundle.noise.aod_pointing_psd_m2_per_hz):
        errors.append("noise.model=none but nonzero noise amplitudes given (noise.*)")
    csv_segments = [s for s in bundle.protocol.segments.values() if s.kind == "csv"]
    for seg in csv_segments:
        if not Path(seg.csv["file"]).is_file():
            errors.append(f"protocol csv file missing: {seg.csv['file']} (segment {seg.id})")
    if errors:
        return ValidationReport(False, errors, warnings, cap, {})

    from .protocol_compile import compile_protocol  # local import: avoid cycle
    try:
        compiled = compile_protocol(bundle)
    except ToolkitError as exc:
        return ValidationReport(False, [f"{exc}"], warnings, cap, {})
    cap["total_duration_s"] = compiled.duration_s
    cap["checkpoint_count"] = len(compiled.checkpoints)
    cap["segments"] = len(compiled.segment_names)
    return ValidationReport(True, [], warnings, cap, effective_config(bundle, compiled))


def effective_config(bundle: ExperimentBundle, compiled=None) -> dict:
    """Merged effective configuration with defaults, overrides and sources.

    Excludes per-segment compile artifacts; used by the `config` CLI command
    to prove which values reach the computation and where they came from.
    """
    from .provenance import json_safe
    cfg = {
        "schema_version": SCHEMA_VERSION,
        "device": {"mass_kg": bundle.device.mass_kg, "species": bundle.device.species,
                   "slm": {"wavelength_m": bundle.device.slm.wavelength_m,
                           "waist_m": bundle.device.slm.waist_m,
                           "nominal_depth_j": bundle.device.slm.nominal_depth_j,
                           "center_m": list(bundle.device.slm.center_m)},
                   "aod": {"wavelength_m": bundle.device.aod.wavelength_m,
                           "waist_m": bundle.device.aod.waist_m,
                           "nominal_depth_j": bundle.device.aod.nominal_depth_j,
                           "center_m": list(bundle.device.aod.center_m)},
                   "lensing_mode": bundle.device.lensing_mode,
                   "lensing_vs_m_s": bundle.device.lensing_vs_m_s},
        "preparation": {"temperature_uK": bundle.preparation.temperature_uK,
                        "sampler": bundle.preparation.sampler,
                        "site_disorder": {
                            "enabled": bundle.preparation.site_disorder_enabled,
                            "slm_depth_rel_sigma": bundle.preparation.slm_depth_rel_sigma,
                            "aod_depth_rel_sigma": bundle.preparation.aod_depth_rel_sigma,
                            "aod_waist_rel_sigma": bundle.preparation.aod_waist_rel_sigma,
                            "alignment_sigma_m": bundle.preparation.alignment_sigma_m,
                            "shot_jitter_sigma_m": bundle.preparation.shot_jitter_sigma_m}},
        "noise": {"model": bundle.noise.model,
                  "aod_rin_psd_per_hz": bundle.noise.aod_rin_psd_per_hz,
                  "aod_pointing_psd_m2_per_hz": bundle.noise.aod_pointing_psd_m2_per_hz,
                  "slm_rin_psd_per_hz": bundle.noise.slm_rin_psd_per_hz,
                  "slm_pointing_psd_m2_per_hz": bundle.noise.slm_pointing_psd_m2_per_hz,
                  "recoil_per_joule": bundle.noise.recoil_per_joule,
                  "legacy_power_w": bundle.noise.legacy_power_w},
        "numerics": {"dt_s": bundle.numerics.dt_s, "noise_dt_s": bundle.numerics.noise_dt_s,
                     "shots": bundle.numerics.shots, "seed": bundle.numerics.seed,
                     "backend": bundle.numerics.backend,
                     "save_trajectories": bundle.numerics.save_trajectories},
        "protocol_summary": {"segments": sorted(bundle.protocol.segments),
                             "final_hold_s": bundle.protocol.final_hold_s},
        "sources": _merged_sources(bundle),
    }
    if compiled is not None:
        cfg["protocol_summary"] = {
            "segments": sorted(bundle.protocol.segments),
            "total_duration_s": compiled.duration_s,
            "checkpoints": [(name, trap) for name, _, trap in compiled.checkpoints],
            "final_hold_s": bundle.protocol.final_hold_s}
    if bundle.calibration is not None:
        cfg["calibration"] = {"version": bundle.calibration.version,
                              "source_tag": bundle.calibration.source_tag,
                              "clock_hz": bundle.calibration.clock_hz}
    if bundle.search is not None:
        cfg["search"] = {"variables": [v.target for v in bundle.search.variables],
                         "objective": bundle.search.objective}
    return json_safe(cfg)


def _merged_sources(bundle: ExperimentBundle) -> dict:
    out: dict = {}
    for section, cfg in (("device", bundle.device), ("preparation", bundle.preparation),
                         ("noise", bundle.noise), ("numerics", bundle.numerics)):
        for key, src in (getattr(cfg, "sources", {}) or {}).items():
            if src:
                out[f"{section}.{key}"] = src
    return out


def _capability(bundle: ExperimentBundle, warnings: list[str]) -> dict:
    """Classify what the run can honestly claim (F1)."""
    sources = _merged_sources(bundle)
    measured = {k for k, v in sources.items() if v.get("tag") == "measured"}
    missing_measurements = [m for m in REQUIRED_MEASUREMENTS if m not in measured]
    has_calibration = bundle.calibration is not None
    if has_calibration:
        missing_measurements = [m for m in missing_measurements if not m.startswith("calibration.")]
    # parameters lacking any explicit source default to *assumed*
    warnings.append("parameters without an explicit __source tag are treated as assumed")
    capability = {
        "runnable_idealized_simulation": True,
        "noise_model": bundle.noise.model,
        "has_hardware_calibration": has_calibration,
        "calibrated_prediction": len(missing_measurements) == 0 and has_calibration,
        "missing_measurements": missing_measurements,
        "source_tag_counts": {tag: sum(1 for v in sources.values() if v.get("tag") == tag)
                              for tag in SOURCE_TAGS},
        "protocol_units_note": "one-way transfer counts are metadata; checkpoints follow segments",
    }
    if has_calibration and bundle.calibration.source_tag == "synthetic":
        warnings.append("calibration is tagged synthetic: RF tables are interface demos, "
                        "not device-ready instructions")
    return capability
