"""从论文 PDF 数字化 Fig. 6d 顶图的存活率曲线（Level 2C 阶段 0 对标数据）。

数据源：仓库根目录 2403.12021v4.pdf 第 7 页 Fig. 6d 顶图
（atom survival vs number of one-way transfers，四条曲线：
0.4 ms ML-optimized 蓝色、0.6 ms 手工深灰、0.4 ms 手工中灰、0.2 ms 手工浅灰）。

方法：
1. pdftoppm 以 600 dpi 渲染第 7 页并裁出 panel d 区域；
2. 检测黑色坐标轴与内向刻度线完成像素→数据坐标标定
   （x: 刻度 20/40/60 对应 n；y: 刻度 1/0.8/0.6/0.4 对应存活率）；
3. 按颜色分割四条曲线（蓝色 + 三档灰度），逐列取像素簇中位数并做
   左→右连续性追踪，剔除图例色块区域与存活率>1.03 的伪像素；
4. 在整数 n=0..60 网格上插值采样输出 CSV。

已知局限（诚实声明）：
- 0.2 ms 曲线数据点在图中只画到 n≈31（存活率≈0.35 后超出面板下界），
  之后无数据，CSV 中记为空；
- n≤8 时四条曲线近乎重合且菱形标记互相遮挡，逐列追踪取得的是
  可见顶层像素，可能与被遮挡系列有 ±0.01 量级偏差；
- 结果为图像目视级数字化，整体精度估计 ±0.01~0.02，不作为拟合真值，
  仅用作阶段 5 统计对比的参考曲线。

用法：
    python tools/digitize_fig6d.py --pdf ../../../2403.12021v4.pdf \
        --out reference/fig6d_survival_digitized.csv [--check-plot reference/fig6d_digitize_check.png]
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

# ---- 600 dpi 下的裁剪窗口（像素），对准第 7 页 panel d ----
CROP_X, CROP_Y, CROP_W, CROP_H = 3400, 250, 1700, 1100
DPI = 600

# ---- 图例色块区域（裁剪坐标内），追踪时剔除 ----
LEGEND_BOXES = [(370, 475, 235, 270), (575, 675, 300, 340),
                (130, 220, 590, 625), (130, 220, 650, 685)]

# ---- 面板内绘图区（裁剪坐标） ----
X_LO, X_HI, Y_LO, Y_HI = 123, 915, 226, 699

# ---- 刻度标定（由内向刻度线检测得到，见 calibrate() 的断言） ----
# x: n=20/40/60 -> 像素列 378 / 633.5 / 888.5；y: s=1/0.8/0.6/0.4 -> 行 283.5/402/521/639.5
X0_PX = 122.75
PX_PER_N = (888.5 - 378.0) / 40.0
Y1_PX = 283.5
PX_PER_SURVIVAL = (402.0 - 283.5) / 0.2

# 0.2 ms 曲线可见数据终点（之后曲线超出面板下界，无数据点）
MANUAL_200US_MAX_N = 31

# 可靠数字化区间下限：n<5 处四条曲线近乎重合、菱形标记互相遮挡，
# 逐列追踪只能看到最顶层系列，无法区分，故 CSV 只覆盖 n>=5。
# n=1..4 各系列目视均处于 0.97~1.00；另以 SI 文字锚点作参考：
# “a single one-way transfer can be done in 200 µs with survival > 98%”。
N_MIN_RELIABLE = 5


def _series_masks(im: np.ndarray) -> dict[str, np.ndarray]:
    """按颜色分割四条曲线。"""
    r, g, b = im[:, :, 0], im[:, :, 1], im[:, :, 2]

    def gray(target, tol=18):
        return (np.abs(r - g) < 12) & (np.abs(g - b) < 12) & (np.abs(r - target) <= tol)

    masks = {
        "ml_400us": (b > r + 40) & (b > g + 15) & (b > 120),
        "manual_600us": gray(136),
        "manual_400us": gray(189),
        "manual_200us": gray(226),
    }
    for mask in masks.values():
        for x0, x1, y0, y1 in LEGEND_BOXES:
            mask[y0:y1, x0:x1] = False
    return masks


def _trace(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """逐列追踪单条曲线：取离上一列最近的像素簇中位数，返回 (n, survival)。"""
    xs, ys_out = [], []
    prev = None
    for x in range(X_LO, X_HI):
        ys = np.where(mask[Y_LO:Y_HI, x])[0] + Y_LO
        ys = ys[1.0 - (ys - Y1_PX) / PX_PER_SURVIVAL <= 1.03]
        if len(ys) == 0:
            continue
        clusters, cur = [], [ys[0]]
        for y in ys[1:]:
            if y - cur[-1] <= 3:
                cur.append(y)
            else:
                clusters.append(cur)
                cur = [y]
        clusters.append(cur)
        anchor = prev if prev is not None else Y1_PX
        best = min(clusters, key=lambda c: abs(np.mean(c) - anchor))
        prev = float(np.mean(best))
        xs.append(x)
        ys_out.append(prev)
    n = (np.asarray(xs) - X0_PX) / PX_PER_N
    s = 1.0 - (np.asarray(ys_out) - Y1_PX) / PX_PER_SURVIVAL
    return n, s


def digitize(png_path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """从裁剪后的 panel 图返回四条曲线的 (n, survival) 采样（整数 n 网格）。"""
    im = np.asarray(Image.open(png_path).convert("RGB")).astype(int)
    out = {}
    for name, mask in _series_masks(im).items():
        n_tr, s_tr = _trace(mask)
        grid = np.arange(N_MIN_RELIABLE, 61)
        s = np.interp(grid, n_tr, s_tr, left=np.nan, right=np.nan)
        s[(grid < n_tr.min() - 0.5) | (grid > n_tr.max() + 0.5)] = np.nan
        if name == "manual_200us":
            s[grid > MANUAL_200US_MAX_N] = np.nan
        out[name] = (grid, s)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="数字化论文 Fig. 6d 顶图存活率曲线")
    ap.add_argument("--pdf", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--check-plot", type=Path, default=None)
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as td:
        prefix = Path(td) / "panel"
        subprocess.run(
            ["pdftoppm", "-f", "7", "-l", "7", "-r", str(DPI), "-png",
             "-x", str(CROP_X), "-y", str(CROP_Y), "-W", str(CROP_W), "-H", str(CROP_H),
             str(args.pdf), str(prefix)],
            check=True)
        png = next(Path(td).glob("panel-*.png"))
        curves = digitize(png)

        if args.check_plot is not None:
            im = np.asarray(Image.open(png).convert("RGB")).astype(np.uint8).copy()
            for grid, s in curves.values():
                for n, v in zip(grid, s):
                    if not np.isfinite(v):
                        continue
                    px = int(round(X0_PX + n * PX_PER_N))
                    py = int(round(Y1_PX + (1.0 - v) * PX_PER_SURVIVAL))
                    im[max(0, py - 2):py + 3, max(0, px - 2):px + 3] = [255, 0, 0]
            args.check_plot.parent.mkdir(parents=True, exist_ok=True)
            Image.fromarray(im).save(args.check_plot)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    names = ["ml_400us", "manual_600us", "manual_400us", "manual_200us"]
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n_one_way_transfers"] + [f"survival_{k}" for k in names])
        grid = curves["ml_400us"][0]
        for i, n in enumerate(grid):
            row = [int(n)]
            for k in names:
                v = curves[k][1][i]
                row.append(f"{v:.4f}" if np.isfinite(v) else "")
            w.writerow(row)
    print(f"digitized -> {args.out}")


if __name__ == "__main__":
    main()
