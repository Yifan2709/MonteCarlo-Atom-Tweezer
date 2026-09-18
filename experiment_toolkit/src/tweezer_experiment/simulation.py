"""Forward prediction API (F3): simulate() and result assembly.

One core serves single transfers, repeated pickup/drop, long AOD transport
and combined sequences — anything expressible as compiled segments. Atom
IDs (array rows) are stable across the whole run; no atom is ever added,
resampled, cooled or renormalized. The denominator of every survival
fraction is the initial shot count (absorbing semantics) unless the caller
explicitly requests recapture-allowed diagnosis (absorb=False), in which
case the record says so.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .errors import ToolkitError
from .provenance import RunRecord, stable_hash, write_json, json_safe
from .schemas import ExperimentBundle, effective_config
from .protocol_compile import compile_protocol
from .physics.field import BeamContext
from .physics.initialization import draw_site_disorder, sample_site_bound, slm_site_energies
from .physics.propagate import propagate

MODEL_VERSION = "1.0.0"
Z95 = 1.959963984540054


def wilson_interval(successes: np.ndarray, n: int):
    p = np.asarray(successes, dtype=float) / n
    z2 = Z95 * Z95
    center = (p + z2 / (2 * n)) / (1 + z2 / n)
    half = Z95 * np.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / (1 + z2 / n)
    return np.maximum(0.0, center - half), np.minimum(1.0, center + half)


def _derive_seed(master: int, stream_key: int) -> int:
    """Documented child-seed rule: (master*1000003 + key) mod 2^63."""
    return (int(master) * 1_000_003 + int(stream_key)) % 2**63


@dataclass
class SimulationResult:
    run: RunRecord
    checkpoints: list                   # (label, trap, time_s)
    survival_counts: np.ndarray         # (n_ck+1,)
    shots: int
    absorb: bool
    survival: np.ndarray
    wilson_low: np.ndarray
    wilson_high: np.ndarray
    initial_unbound: int
    energy_stats: list                  # per checkpoint dict for survivors
    loss_histogram: dict                # segment name -> count
    noise_ledger: dict
    switch_work_mean: np.ndarray
    checkpoint_states: np.ndarray | None
    final_states: np.ndarray | None
    alive: np.ndarray | None = None
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "run": json_safe(self.run),
            "absorbing": self.absorb,
            "denominator": self.shots,
            "checkpoints": [{"label": label, "trap": trap, "time_s": t}
                            for label, trap, t in self.checkpoints],
            "survival_counts": self.survival_counts.tolist(),
            "survival": self.survival.tolist(),
            "wilson95_low": self.wilson_low.tolist(),
            "wilson95_high": self.wilson_high.tolist(),
            "initial_unbound_actual_site": self.initial_unbound,
            "energy_stats": self.energy_stats,
            "loss_histogram": self.loss_histogram,
            "noise_ledger": self.noise_ledger,
            "switch_work_mean_j": self.switch_work_mean.tolist(),
            "warnings": self.warnings,
            "interval_note": "Wilson 95% per checkpoint reflects Monte Carlo sampling "
                             "of independent atoms only (no common-path noise in phase 1)",
        }


def simulate(bundle: ExperimentBundle, protocol=None, *, shots: int | None = None,
             seed: int | None = None, absorb: bool = True,
             compiled=None) -> SimulationResult:
    shots = shots if shots is not None else bundle.numerics.shots
    master = seed if seed is not None else bundle.numerics.seed
    if compiled is None:
        spec = protocol if protocol is not None else bundle.protocol
        compiled = compile_protocol(bundle, spec)

    prep = bundle.preparation
    disorder = draw_site_disorder(prep, shots, _derive_seed(master, 1))
    ctx = BeamContext.from_device(bundle.device)
    pool = sample_site_bound(ctx, prep, shots, _derive_seed(master, 2),
                             disorder["slm_depth_factor"])
    pos, vel = pool["pos_m"], pool["vel_m_per_s"]
    jitter = (np.random.default_rng(_derive_seed(master, 3)).normal(
              0.0, prep.shot_jitter_sigma_m, (len(compiled.segment_bounds), shots))
              if prep.shot_jitter_sigma_m > 0 else
              np.zeros((len(compiled.segment_bounds), shots)))
    site = (disorder["slm_depth_factor"], disorder["aod_depth_factor"],
            disorder["aod_waist_factor"], disorder["alignment_offset_m"])
    raw = propagate(compiled, pos, vel, ctx, site, jitter, bundle.noise,
                    bundle.numerics.noise_dt_s,
                    _derive_seed(master, 4),
                    absorb=absorb, save_states=bundle.numerics.save_trajectories)

    counts = raw["alive"].sum(axis=0).astype(int)
    survival = counts / shots
    low, high = wilson_interval(counts, shots)
    initial_unbound = int(np.sum(pool["initial_energy_j"] >= 0))
    warnings = []
    if initial_unbound:
        warnings.append(f"{initial_unbound} atoms prepared unbound — sampler contract violated")
    energy_stats = []
    alive = raw["alive"]
    energies = raw["checkpoint_energy_j"]
    for k in range(alive.shape[1]):
        mask = alive[:, k]
        stats = {"checkpoint_index": k, "survivors": int(mask.sum())}
        if mask.any():
            e = energies[mask, k] / 1.380649e-29  # report in uK
            stats.update({"energy_uK_mean": float(np.mean(e)),
                          "energy_uK_median": float(np.median(e)),
                          "energy_uK_p90": float(np.quantile(e, 0.9))})
        energy_stats.append(stats)
    names = compiled.segment_names
    hist: dict[str, int] = {}
    for instance in raw["lost_instance"]:
        if instance >= 0:
            hist[names[instance]] = hist.get(names[instance], 0) + 1
    ck_meta = []
    for label, idx, trap in compiled.checkpoints:
        ck_meta.append((label, trap, float(compiled.controls[compiled.segment_bounds[idx][1], 0])))
    ledger = raw["noise_ledger_j"]
    inputs = effective_config(bundle, compiled)
    run = RunRecord.begin(stable_hash(inputs), stable_hash({} if bundle.calibration is None
                                                           else {"version": bundle.calibration.version}),
                          {"master": master, "streams": ["disorder", "initial", "jitter", "noise"]},
                          backend=bundle.numerics.backend)
    return SimulationResult(
        run=run, checkpoints=ck_meta, survival_counts=counts, shots=shots, absorb=absorb,
        survival=survival, wilson_low=low, wilson_high=high,
        initial_unbound=initial_unbound, energy_stats=energy_stats,
        loss_histogram=hist,
        noise_ledger={"mean_realized_j_per_atom": float(np.mean(ledger[:, 0])),
                      "mean_expected_j_per_atom": float(np.mean(ledger[:, 1])),
                      "model": bundle.noise.model},
        switch_work_mean=raw["switch_work_j"].mean(axis=0),
        checkpoint_states=raw["checkpoint_states"], final_states=raw["final_states"],
        alive=raw["alive"], warnings=warnings)


def save_result(directory: str | Path, result: SimulationResult, compiled=None) -> Path:
    """Persist result.json (+ per-atom alive labels / states npz when stored)."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory / "result.json", result.to_dict())
    if result.alive is not None:
        payload = {"alive": result.alive}
        if result.checkpoint_states is not None:
            payload["states"] = result.checkpoint_states
        if result.final_states is not None:
            payload["final_states"] = result.final_states
        np.savez_compressed(directory / "per_atom.npz", **payload)
    return directory / "result.json"
