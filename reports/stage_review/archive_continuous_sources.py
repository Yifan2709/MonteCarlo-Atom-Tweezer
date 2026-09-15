"""Archive exact source/config inputs alongside the new simulation outputs."""
from pathlib import Path
import hashlib
import json
from zipfile import ZipFile, ZIP_DEFLATED

repo = Path(__file__).resolve().parents[2]
project = repo/"simulation/level2_joint_transfer"
out = project/"outputs/continuous_level2c_345/run_20260914"
files = list((project/"src").rglob("*.py"))
files += list((project/"configs").rglob("*.yaml"))
files += list((project/"tests").glob("test*continuous*.py"))
files += [project/"pyproject.toml",project/"reference/fig6d_survival_digitized.csv",
          project/"outputs/level2c_pickup_survival/scan_stage1/fig6c_ml_trajectory_digitized.json",
          Path(__file__).parent/"continuous_precision.py"]
manifest = []
with ZipFile(out/"reproducibility_sources.zip","w",ZIP_DEFLATED) as archive:
    for file in sorted(set(files)):
        name = file.relative_to(repo).as_posix()
        archive.write(file,name)
        manifest.append({"path":name,"sha256":hashlib.sha256(file.read_bytes()).hexdigest()})
(out/"source_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
print(f"Archived {len(manifest)} source/configuration files")
