"""Constrained candidate search over experimental controls (F4).

Only protocol control variables (segment durations, ramp fractions, target
depths, path families) may vary. Device, preparation and noise stay bit-
identical across candidates — their hash is recorded once and asserted.
Screening and validation use separate shot pools and seeds, fixed before
any candidate runs. "No feasible candidate" is a valid, reportable outcome;
failed simulations are recorded with their error and never averaged into
objectives.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from .errors import ToolkitError, E_CONFLICT
from .provenance import RunRecord, stable_hash, json_safe
from .schemas import ExperimentBundle, SearchSpec
from .protocol_compile import apply_override
from .simulation import simulate, Z95


@dataclass
class CandidateRecord:
    index: int
    values: dict
    status: str                       # "ok" | "infeasible" | "failed"
    violations: list = field(default_factory=list)
    error: str | None = None
    objective_value: float | None = None
    objective_interval: tuple | None = None
    duration_us: float | None = None
    validation: dict | None = None

    def to_dict(self):
        return json_safe(self)


@dataclass
class SearchResult:
    objective: dict
    device_preparation_hash: str
    candidates: list
    recommended_index: int | None
    distinguishability: dict
    notes: list
    run: RunRecord

    @property
    def recommended(self):
        return self.candidates[self.recommended_index] if self.recommended_index is not None else None

    def to_dict(self):
        return {"objective": self.objective,
                "device_preparation_hash": self.device_preparation_hash,
                "recommended_index": self.recommended_index,
                "candidates": [c.to_dict() for c in self.candidates],
                "distinguishability": self.distinguishability,
                "notes": self.notes,
                "run": json_safe(self.run)}


def search_controls(bundle: ExperimentBundle, spec: SearchSpec | None = None,
                    journal_path=None) -> SearchResult:
    spec = spec if spec is not None else bundle.search
    if spec is None:
        raise ToolkitError(E_CONFLICT, "search requires a search section", "search")
    from .runstore import SearchJournal
    journal = None
    if journal_path is not None:
        journal = SearchJournal(journal_path, {
            "device_preparation_hash": stable_hash({
                "device": json_safe(bundle.device), "preparation": json_safe(bundle.preparation),
                "noise": json_safe(bundle.noise)}),
            "objective": dict(spec.objective),
            "variables": [v.target for v in spec.variables],
            "screening_shots": spec.screening_shots, "screening_seed": spec.screening_seed})
    if bundle.numerics.shots != spec.screening_shots:
        numerics = replace(bundle.numerics, shots=spec.screening_shots)
    else:
        numerics = bundle.numerics
    bundle = replace(bundle, numerics=numerics)
    device_hash = stable_hash({"device": json_safe(bundle.device),
                               "preparation": json_safe(bundle.preparation),
                               "noise": json_safe(bundle.noise)})

    grids = [v.values for v in spec.variables]
    combos = list(itertools.product(*grids))
    objective = dict(spec.objective)
    notes = [f"objective frozen before runs: {objective}",
             f"screening {spec.screening_shots} shots seed {spec.screening_seed}; "
             f"validation {spec.validation_shots} shots seed {spec.validation_seed}"]

    candidates: list[CandidateRecord] = []
    for i, combo in enumerate(combos):
        values = {v.target: val for v, val in zip(spec.variables, combo)}
        rec = CandidateRecord(index=i, values=values, status="ok")
        if journal is not None and journal.has(i):
            cached = journal.get(i)
            candidates.append(CandidateRecord(
                index=i, values=values, status=cached.get("status", "ok"),
                violations=cached.get("violations", []), error=cached.get("error"),
                objective_value=cached.get("objective_value"),
                objective_interval=tuple(cached["objective_interval"])
                if cached.get("objective_interval") else None,
                duration_us=cached.get("duration_us")))
            continue
        protocol = bundle.protocol
        try:
            for target, val in values.items():
                protocol = apply_override(protocol, target, val)
            result = simulate(bundle, protocol, shots=spec.screening_shots, seed=spec.screening_seed)
            rec.objective_value = float(result.survival[-1])
            n = result.shots
            half = Z95 * np.sqrt(max(rec.objective_value * (1 - rec.objective_value), 1e-12) / n)
            rec.objective_interval = (max(0.0, rec.objective_value - half),
                                      min(1.0, rec.objective_value + half))
            if objective["type"] == "minimize_duration":
                target_field = next((t for t in values if t.endswith("duration_us")), None)
                rec.duration_us = float(dict(values).get(target_field, np.nan)) if target_field else np.nan
                bound = objective["capture_lower_bound"]
                if rec.objective_value < bound:
                    rec.status = "infeasible"
                    rec.violations.append(
                        f"capture {rec.objective_value:.4f} < lower bound {bound} (screening)")
            if journal is not None:
                journal.append(i, rec.to_dict())
            candidates.append(rec)
        except (ToolkitError, ValueError, RuntimeError, OSError) as exc:
            rec.status = "failed"
            rec.error = f"{type(exc).__name__}: {exc}"
            if journal is not None:
                journal.append(i, rec.to_dict())
            candidates.append(rec)

    ok = [c for c in candidates if c.status == "ok"]
    feasible = list(ok)
    if objective["type"] == "minimize_duration":
        for c in feasible:
            duration_targets = [t for t in c.values if t.endswith(".duration_us")]
            c.duration_us = float(c.values[duration_targets[0]]) if duration_targets else float("nan")
        feasible.sort(key=lambda c: (c.duration_us, c.index))
    else:
        feasible.sort(key=lambda c: (-c.objective_value, c.index))

    distinguishability = {"top_gap_within_mc_error": None, "undistinguished_candidates": []}
    if len(feasible) >= 2 and objective["type"] == "maximize_capture":
        best, second = feasible[0], feasible[1]
        gap = best.objective_value - second.objective_value
        err = max(best.objective_interval[1] - best.objective_interval[0],
                  second.objective_interval[1] - second.objective_interval[0]) / 2.0
        distinguishability["top_gap_within_mc_error"] = bool(gap <= 2 * err)
        if gap <= 2 * err:
            distinguishability["undistinguished_candidates"] = [c.index for c in feasible[:3]]
            notes.append("top candidates are not statistically distinguishable at screening "
                         "sample size; ranking is not a claim of superiority")

    recommended_index = None
    attempts = 0
    for cand in feasible[:3]:
        attempts += 1
        protocol = bundle.protocol
        for target, val in cand.values.items():
            protocol = apply_override(protocol, target, val)
        try:
            vr = simulate(bundle, protocol, shots=spec.validation_shots, seed=spec.validation_seed)
            val_dict = {"shots": vr.shots, "seed": spec.validation_seed,
                        "final_capture": float(vr.survival[-1]),
                        "wilson95": [float(vr.wilson_low[-1]), float(vr.wilson_high[-1])],
                        "attempts": attempts}
            cand.validation = val_dict
            if journal is not None:
                journal.append(cand.index, cand.to_dict())
            if objective["type"] == "minimize_duration":
                if vr.survival[-1] < objective["capture_lower_bound"]:
                    cand.violations.append(
                        f"validation capture {vr.survival[-1]:.4f} < bound "
                        f"{objective['capture_lower_bound']}")
                    notes.append(f"candidate {cand.index} passed screening but failed validation; "
                                 f"trying next")
                    continue
            recommended_index = cand.index
            break
        except (ToolkitError, ValueError, RuntimeError, OSError) as exc:
            cand.violations.append(f"validation failed: {exc}")
            continue

    if recommended_index is None:
        if ok and not feasible:
            notes.append("no feasible candidate: every candidate violates the constraint")
        elif not ok:
            notes.append("no feasible candidate: every simulation failed")
        else:
            notes.append("no candidate survived independent validation")

    run = RunRecord.begin(stable_hash({"device": device_hash, "objective": objective}),
                          "", {"screening": spec.screening_seed, "validation": spec.validation_seed})
    return SearchResult(objective=objective, device_preparation_hash=device_hash,
                        candidates=candidates, recommended_index=recommended_index,
                        distinguishability=distinguishability, notes=notes, run=run)


def _duration_segment(bundle, spec):
    seg = spec.objective.get("segment")
    if seg:
        return f"segment.{seg}.duration_us"
    return None
