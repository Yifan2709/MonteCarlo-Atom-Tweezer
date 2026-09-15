"""Frozen confirmation plan, chosen from the 256-shot scan only.

T=15 is the best shared scanned setting; T=5/15/25/25 are per-waveform minima.
No optimization or re-selection uses the confirmation pool. All four curves are
also rerun at each selected shared setting, preventing endpoint cherry-picking.
"""
from pathlib import Path
import numpy as np
import yaml
from numba import set_num_threads
from continuous_transfer.cli import ROOT, run_case, write_json

settings = yaml.safe_load((ROOT/"configs/continuous_level2c_345.yaml").read_text())
out = ROOT/"outputs/continuous_level2c_345/run_20260914"
set_num_threads(12)
plan = {"stage":"precision","temperatures_uK":[5,15,25],"arm":"matched",
        "shots":4096,"seed":23003,"dt_us":[0.025,0.0125],
        "selection":"256-shot seed 23001 scan only; all four curves at each selected T",
        "selection_rule":"minimize mean of four curve MSEs, equal curve weights",
        "per_waveform_best_scan_T":{"manual_200us":5,"manual_400us":15,
                                    "manual_600us":25,"ml_400us":25}}
write_json(out/"precision_plan.json",plan)
for temperature in plan["temperatures_uK"]:
    for key in settings["waveforms"]:
        for dt in plan["dt_us"]:
            run_case(settings,stage="precision",temperature=temperature,key=key,arm="matched",
                     shots=4096,seed=23003,dt_us=dt,out=out)
