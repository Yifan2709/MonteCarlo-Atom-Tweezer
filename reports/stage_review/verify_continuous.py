"""Verify numerical archives, report hashes and native PPT evidence."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np

root = Path(__file__).resolve().parents[2]
project = root/"simulation/level2_joint_transfer"
report = root/"reports/stage_review"
stem = "AOD_SLM_continuous_level2c_345_2026-09-14_v10"
e = json.loads((report/"output"/(stem+".evidence.json")).read_text(encoding="utf-8"))
for item in e["result_manifest"]:
    assert hashlib.sha256((root/item["path"]).read_bytes()).hexdigest() == item["sha256"]
files = list((project/"outputs/continuous_level2c_345/run_20260914").glob("*_N*.json"))
for file in files:
    r = json.loads(file.read_text(encoding="utf-8"))
    with np.load(file.with_suffix(".npz")) as a:
        alive = a["alive"]
        assert alive.shape == (r["shots"],61)
        assert np.all(np.diff(alive.astype(int),axis=1)<=0)
        np.testing.assert_array_equal(alive.sum(0),r["survival_counts"])
        np.testing.assert_allclose(a["injected_j"],a["elapsed_s"]*r["heating_power_uK_per_s"]["total"]*1.380649e-29,
                                   rtol=1e-9,atol=1e-38)
        assert np.all(np.isfinite(a["states"][alive]))
        ab = a["spin_ab"][alive]
        assert np.max(abs(np.sum(abs(ab)**2,axis=1)-1)) < 1e-10
    for source,sha in r["sources"].items():
        assert hashlib.sha256((project/source).read_bytes()).hexdigest() == sha
v9 = report/"output/AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9_reviewed.pptx"
v10 = report/"output"/(stem+"_reviewed.pptx")
with ZipFile(v9) as z9, ZipFile(v10) as z10:
    old = [n for n in z9.namelist() if n.startswith("ppt/embeddings/") and n.endswith(".xlsx")]
    for n in old:
        assert z9.read(n) == z10.read(n)
    slides = [n for n in z10.namelist() if n.startswith("ppt/slides/slide") and n.endswith(".xml")]
    charts = [n for n in z10.namelist() if "/charts/" in n and n.endswith(".xml")
              and ET.fromstring(z10.read(n)).tag.endswith("}chartSpace")]
    workbooks = [n for n in z10.namelist() if n.startswith("ppt/embeddings/") and n.endswith(".xlsx")]
    assert (len(slides),len(charts),len(workbooks)) == (43,21,21)
    text = "".join(ET.fromstring(z10.read("ppt/slides/slide32.xml")).itertext())
    assert "最新波形连续接入、重算及数值复核已在后续页完成" in text
qa = {"cases_verified":len(files),"result_hashes_verified":len(e["result_manifest"]),
      "slides":len(slides),"charts":len(charts),"embedded_workbooks":len(workbooks),
      "original_workbooks_unchanged":len(old),"monotone_loss_and_heat_and_unitarity":"passed",
      "visual_review":"43-page overview plus all changed/new pages inspected; slide 32 corrected and re-rendered"}
(report/".build/v10/final_qa.json").write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(qa,ensure_ascii=False))
