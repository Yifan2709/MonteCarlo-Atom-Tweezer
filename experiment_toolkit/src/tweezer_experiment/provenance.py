"""Parameter provenance tags and run records (F1 contract).

Every physical parameter carries a source tag so reports can separate
measured device facts from assumptions. Run records capture everything
needed to decide whether a cached result may be reused.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import platform
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

#: Closed set of source tags (see PROMPT F1: measured / derived / assumed /
#: fitted / synthetic). A value tagged *fitted* is a calibration output and
#: must never be presented as an independently measured device parameter.
SOURCE_TAGS = ("measured", "derived", "assumed", "fitted", "synthetic")


@dataclass(frozen=True)
class Sourced:
    """A scalar with a source tag and optional note."""
    value: float
    source: str
    note: str = ""

    def __post_init__(self):
        if self.source not in SOURCE_TAGS:
            raise ValueError(f"unknown source tag {self.source!r}")


def stable_hash(obj: Any) -> str:
    """Order-stable SHA-256 of a JSON-able object (first 16 hex chars)."""
    payload = json.dumps(obj, sort_keys=True, default=_default, allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _default(o):
    if dataclasses.is_dataclass(o) and not isinstance(o, type):
        return asdict(o)
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"not JSON-able: {type(o)}")


@dataclass
class RunRecord:
    """Identity and environment of one real run."""
    run_id: str
    schema_version: str
    model_version: str
    started_utc: str
    inputs_hash: str
    calibration_hash: str
    seeds: dict
    backend: str
    versions: dict
    elapsed_s: float = 0.0
    status: str = "ok"

    @staticmethod
    def begin(inputs_hash: str, calibration_hash: str, seeds: dict,
              backend: str = "numpy", schema_version: str = "1",
              model_version: str = "1.0.0") -> "RunRecord":
        return RunRecord(
            run_id=f"run_{time.strftime('%Y%m%dT%H%M%S')}_{inputs_hash[:8]}",
            schema_version=schema_version, model_version=model_version,
            started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            inputs_hash=inputs_hash, calibration_hash=calibration_hash,
            seeds=dict(seeds), backend=backend,
            versions={"python": platform.python_version(),
                      "numpy": __import__("numpy").__version__,
                      "toolkit": model_version})


def json_safe(obj: Any) -> Any:
    """Convert to JSON-safe structure; non-finite floats become null."""
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (bool, int, str)) or obj is None:
        return obj
    if isinstance(obj, float):
        return obj if obj == obj and abs(obj) != float("inf") else None
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return json_safe(asdict(obj))
    if hasattr(obj, "item"):
        return json_safe(obj.item())
    return str(obj)


def write_json(path: Path, obj: Any) -> None:
    Path(path).write_text(json.dumps(json_safe(obj), indent=2, ensure_ascii=False),
                          encoding="utf-8")
