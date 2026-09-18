"""Hash-keyed run store and resumable search journal (A15).

Caching rule (explicit, never silent): a stored result may be reused only
when the current input hash (effective configuration incl. device,
preparation, noise, protocol), model version, seed and shot count all match
the stored record. Any change -> the store reports the mismatch reason and
the caller recomputes. The search journal appends one JSON line per finished
candidate so an interrupted scan resumes without redoing or mixing pools.
"""
from __future__ import annotations

import json
from pathlib import Path

from .provenance import stable_hash, write_json, json_safe
from .schemas import ExperimentBundle, effective_config
from .simulation import MODEL_VERSION


def run_signature(bundle: ExperimentBundle, *, seed: int | None = None,
                  shots: int | None = None) -> dict:
    """Everything a cached physics result depends on (code version included)."""
    return {"model_version": MODEL_VERSION,
            "inputs_hash": stable_hash(effective_config(bundle)),
            "seed": seed if seed is not None else bundle.numerics.seed,
            "shots": shots if shots is not None else bundle.numerics.shots}


class RunStore:
    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def save(self, bundle: ExperimentBundle, result, *, seed=None, shots=None) -> Path:
        seed = seed if seed is not None else int(result.run.seeds.get("master", 0))
        shots = shots if shots is not None else int(result.shots)
        sig = run_signature(bundle, seed=seed, shots=shots)
        path = self.root / f"{sig['inputs_hash']}_{sig['seed']}_{sig['shots']}.json"
        write_json(path, {"signature": sig, "result": result.to_dict()})
        return path

    def load(self, bundle: ExperimentBundle, *, seed=None, shots=None) -> dict:
        """Return {"hit": stored} on exact match, else {"hit": None, "reason": ...}."""
        sig = run_signature(bundle, seed=seed, shots=shots)
        for path in self.root.glob(f"{sig['inputs_hash']}_{sig['seed']}_*.json"):
            stored = json.loads(path.read_text(encoding="utf-8"))
            s = stored["signature"]
            mismatches = [k for k in ("model_version", "inputs_hash", "seed", "shots")
                          if s.get(k) != sig[k]]
            if not mismatches:
                return {"hit": stored["result"], "path": str(path)}
            return {"hit": None, "reason": f"stale cache: {mismatches} differ ({path.name})",
                    "stored": s}
        return {"hit": None, "reason": "no stored run for this signature"}


class SearchJournal:
    """Append-only candidate journal for resumable searches."""

    def __init__(self, path: str | Path, header: dict):
        self.path = Path(path)
        self.header = header
        self.records: dict[int, dict] = {}
        if self.path.exists():
            for line in self.path.read_text(encoding="utf-8").splitlines():
                if not line.strip():
                    continue
                entry = json.loads(line)
                if "index" not in entry:      # header line
                    continue
                if entry.get("header") != header:
                    self.records = {}      # different objective/device: never mix
                    self.path.unlink(missing_ok=True)
                    break
                self.records[entry["index"]] = entry["record"]
        if not self.path.exists():
            self.path.write_text(json.dumps({"header": header}) + "\n", encoding="utf-8")

    def has(self, index: int) -> bool:
        return index in self.records

    def get(self, index: int) -> dict:
        return self.records[index]

    def append(self, index: int, record: dict) -> None:
        self.records[index] = record
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({"header": self.header, "index": index,
                                     "record": json_safe(record)}) + "\n")
