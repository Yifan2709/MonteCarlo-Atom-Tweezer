# Level 1: 1D classical AOD -> SLM transfer

This directory was copied from the validated Level 0 project. The
`src/level0_static_trap` package and its baseline tests remain intact. The new
`src/level1_transfer_1d` package implements a finite-temperature classical,
one-dimensional AOD-to-SLM merge/drop-off and a 50-shot Monte Carlo baseline.

## Run

From this directory, using the same Python interpreter:

```bash
python -m pip install -e ".[dev]"
python -m pytest
python -m level0_static_trap.cli --config configs/default.yaml --trap both \
  --output-dir outputs/level0_static_trap/demo
python -m level1_transfer_1d.cli --config configs/transfer_1d.yaml
```

Level 1 writes the configured waveform, three potential snapshots, ideal and
thermal example trajectories, a per-shot Monte Carlo table, diagnostic plots,
`metrics.json`, and `summary.md` under
`outputs/level1_transfer_1d/demo`.

## Interpretation and scope

The reported `classical_capture_fraction` is a classical mechanical result for
this configured model. It is not experimental survival, AOD-SLM transfer
fidelity, quantum-state fidelity, or directly comparable with a literature
99.81% benchmark.

The model includes finite-temperature classical motion in a time-dependent 1D
Gaussian optical potential. It intentionally excludes hardware noise,
parameter scans, 2D/3D or axial motion, photon scattering, vacuum collisions,
multi-atom interactions, quantum internal states, coherence, and experimental
loss mechanisms.

The next appropriate extension is a transfer-duration scan; it is not included
in this Level 1 implementation.
