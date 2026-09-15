"""Scientific comparison and source-overlay QA for the corrected PDF reference."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pypdfium2 as pdfium


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
e = json.loads((HERE / "output/Level5_vector_reference_v2.evidence.json").read_text("utf-8"))
vectors = e["reference_vectors"]
plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False, "font.size": 11})
names = ["手工 200 μs", "手工 400 μs", "手工 600 μs", "论文 ML 400 μs"]
fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.5), layout="constrained")
for ax, row, label in zip(axes.flat, e["selected"], names):
    s = vectors["series"][row["waveform"]]
    cp = row["comparison_paper_fit"]
    n = np.array([m["n"] for m in s["markers"]])
    y = np.array([m["survival"] for m in s["markers"]])
    low, high = (np.array([m[k] for m in s["markers"]]) for k in ["lower", "upper"])
    ax.plot(s["display_n"], s["display_fit"], color="#17324D", lw=2, label="论文拟合线（PDF 原路径）")
    ax.errorbar(n, y, yerr=np.vstack((y-low, high-y)), fmt="D", markersize=4, capsize=2,
                color="#17324D", elinewidth=1, label="论文散点与原图误差条")
    simn = np.arange(cp["n"][-1]+1)
    ax.plot(simn, np.asarray(row["survival"])[simn], color="#BB3E35", lw=2, label="连续 Monte Carlo")
    ax.set(xlim=(0, cp["n"][-1]), ylim=(0, 1.03), xlabel="单向转移数 n", ylabel="存活概率 S(n)",
           title=f"{label}，冻结 T={row['temperature_uK']} μK")
    ax.grid(alpha=.2)
    ax.legend(fontsize=9, loc="lower left")
fig.savefig(HERE / "output/Level5_修正参考曲线对照.png", dpi=180)
fig.savefig(HERE / "output/Level5_修正参考曲线对照.svg")
plt.close(fig)

# Overlay every extracted object onto its exact position in a paper crop.
doc = pdfium.PdfDocument(ROOT / vectors["source_pdf"])
im = doc[6].render(scale=5).to_pil()
crop = (410, 55, 520, 115)  # PDF screen coordinates (pt).
im = im.crop(tuple(int(v*5) for v in crop))
fig, ax = plt.subplots(figsize=(12, 6.6), layout="constrained")
ax.imshow(im, extent=[crop[0], crop[2], crop[3], crop[1]])
for s in vectors["series"].values():
    pts = np.array(s["fit_path_pdf_pt"])
    ax.plot(pts[:,0], pts[:,1], color="red", lw=.7, alpha=.8)
    centers = np.array([m["pdf_center_pt"] for m in s["markers"]])
    ax.plot(centers[:,0], centers[:,1], "+", color="red", markersize=6, markeredgewidth=.8)
ax.set(xlim=(crop[0],crop[2]), ylim=(crop[3],crop[1]), title="源图叠加核验：红线为拟合路径，红十字为独立散点中心")
ax.set_axis_off()
fig.savefig(HERE / ".build/vector_reference_overlay.png", dpi=180)
plt.close(fig)
print("Saved comparison PNG/SVG and source overlay")
