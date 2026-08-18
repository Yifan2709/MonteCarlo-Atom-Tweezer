# Level 1: 1D classical AOD -> SLM transfer

This is a finite-temperature classical mechanical model, not experimental survival, transfer fidelity, or quantum fidelity. It is not directly comparable to the literature 99.81% benchmark.

## Q1 - Time-dependent integrator static limit

PASS

## Q2 - Waveform

- Initial separation: 2.400 um
- Move duration: 312.000 us
- Ramp duration: 288.000 us
- AOD depth: 280.0 -> 0.0 uK
- SLM depth: 140.0 uK

## Q3 - Ideal atom

- Final energy / SLM depth: -0.99999966
- Captured: True

## Q4 - Timestep convergence

- dt / dt-half: 0.050 / 0.025 us
- Normalized final-energy difference: 1.045e-09
- Capture classification consistent: True

## Q5 - Thermal example

- Initial AOD energy: -279.245274 uK
- Final SLM energy: -139.477200 uK
- Excitation-energy change: -0.231926 uK
- Captured: True

## Q6 - Classical capture fraction

- Captured / total: 50 / 50
- classical_capture_fraction: 1.000000
- Wilson 95% interval: [0.928652, 1.000000]

## Q7 - Capture margin

- Mean / median capture margin: 0.9668631000656496 / 0.9818248249353418
- See `final_slm_energy_distribution.png`.

## Q8 - Classical motional excitation proxy

- Mean excitation-energy change: -1.924299 uK
- Median excitation-energy change: -1.061504 uK

## Validation

All passed: True
