"""Generate the Level-1 AOD->SLM transfer Monte Carlo PPT (16:9, Chinese).

Embeds the figures produced by make_ppt_figures.py and the numbers from
outputs/level1_transfer_1d/demo/metrics.json. Run from this directory:

    python make_transfer_ppt.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"
CROSSCHECK = HERE / "scan_parameters"  # 独立副本 run_parameter_scans.py 的输出
METRICS = HERE.parent / "outputs" / "level1_transfer_1d" / "demo" / "metrics.json"
OUTPUT = HERE / "Level1_Transfer_MC.pptx"

NAVY = RGBColor(0x1F, 0x38, 0x64)
BLUE = RGBColor(0x2E, 0x74, 0xB5)
GRAY = RGBColor(0x40, 0x40, 0x40)
LIGHT = RGBColor(0xF2, 0xF6, 0xFB)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
STEP_COLORS = [
    RGBColor(0x44, 0x72, 0xC4), RGBColor(0x2E, 0xA0, 0x6E), RGBColor(0xE8, 0x7D, 0x2D),
    RGBColor(0xC0, 0x50, 0x4D), RGBColor(0x7B, 0x57, 0xA0),
]
FONT = "PingFang SC"

SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)


def _set_ea_font(run):
    """让中文字符也使用指定字体（python-pptx 默认只设置拉丁字体）。"""
    run.font.name = FONT
    rPr = run._r.get_or_add_rPr()
    for tag in ("a:ea", "a:cs"):
        element = rPr.find(qn(tag))
        if element is None:
            element = rPr.makeelement(qn(tag), {})
            rPr.append(element)
        element.set("typeface", FONT)


def add_text(slide, left, top, width, height, lines, size=16, color=GRAY,
             bold=False, align=PP_ALIGN.LEFT, line_spacing=1.15):
    box = slide.shapes.add_textbox(left, top, width, height)
    frame = box.text_frame
    frame.word_wrap = True
    if isinstance(lines, str):
        lines = [lines]
    for i, line in enumerate(lines):
        p = frame.paragraphs[0] if i == 0 else frame.add_paragraph()
        p.alignment = align
        p.line_spacing = line_spacing
        run = p.add_run()
        run.text = line
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = color
        _set_ea_font(run)
    return box


def add_picture_fit(slide, name, left, top, max_w, max_h, caption=None):
    # figures/ 下的成品图优先；否则按 pptplots/ 下相对路径（如扫描目录）解析
    path = FIGURES / name if (FIGURES / name).exists() else HERE / name
    with Image.open(path) as image:
        px_w, px_h = image.size
    scale = min(max_w / px_w, max_h / px_h)
    w, h = int(px_w * scale), int(px_h * scale)
    picture = slide.shapes.add_picture(str(path), int(left + (max_w - w) / 2), int(top + (max_h - h) / 2), w, h)
    if caption:
        add_text(slide, left, top + max_h + Emu(30000), max_w, Inches(0.35), caption,
                 size=11, color=BLUE, align=PP_ALIGN.CENTER)
    return picture


def blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def content_slide(prs, kicker, title):
    """带顶部标题条的标准内容页。"""
    slide = blank(prs)
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, SLIDE_W, Inches(1.0))
    bar.fill.solid()
    bar.fill.fore_color.rgb = NAVY
    bar.line.fill.background()
    bar.shadow.inherit = False
    add_text(slide, Inches(0.55), Inches(0.10), Inches(11), Inches(0.35), kicker, size=13, color=RGBColor(0x9D, 0xC3, 0xE6), bold=True)
    add_text(slide, Inches(0.55), Inches(0.38), Inches(12), Inches(0.55), title, size=26, color=WHITE, bold=True)
    return slide


def bullets(slide, items, left=Inches(0.55), top=Inches(1.35), width=Inches(4.9), size=15):
    lines = []
    for item in items:
        if isinstance(item, tuple):  # (text, indent-level)
            text, level = item
        else:
            text, level = item, 0
        prefix = "• " if level == 0 else "– "
        lines.append(f"{prefix}{text}")
    return add_text(slide, left, top, width, Inches(5.5), lines, size=size, line_spacing=1.3)


def bullet_runs(slide, items, left, top, width):
    """支持局部加粗/变色的项目符号列表：item = [(text, bold, color), ...]。"""
    box = slide.shapes.add_textbox(left, top, width, Inches(5.6))
    frame = box.text_frame
    frame.word_wrap = True
    for i, item in enumerate(items):
        p = frame.paragraphs[0] if i == 0 else frame.add_paragraph()
        p.line_spacing = 1.3
        lead = p.add_run()
        lead.text = "• "
        lead.font.size = Pt(15)
        lead.font.color.rgb = BLUE
        lead.font.bold = True
        _set_ea_font(lead)
        for text, bold, color in item:
            run = p.add_run()
            run.text = text
            run.font.size = Pt(15)
            run.font.bold = bold
            run.font.color.rgb = color
            _set_ea_font(run)
    return box


# ---------------------------------------------------------------- 幻灯片 1：封面
def slide_cover(prs):
    slide = blank(prs)
    background = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, SLIDE_W, SLIDE_H)
    background.fill.solid()
    background.fill.fore_color.rgb = NAVY
    background.line.fill.background()
    background.shadow.inherit = False
    add_text(slide, Inches(0.9), Inches(2.2), Inches(11.5), Inches(1.2),
             "一维经典 AOD → SLM 转移的蒙特卡洛模拟", size=36, color=WHITE, bold=True)
    add_text(slide, Inches(0.9), Inches(3.4), Inches(11.5), Inches(0.6),
             "四步流水线：预检查 → 热初态采样 → 轨迹推进 → 判定与统计", size=20, color=RGBColor(0x9D, 0xC3, 0xE6))
    add_text(slide, Inches(0.9), Inches(5.9), Inches(11.5), Inches(1.0),
             ["Cs-133 · 高斯光镊 · T = 5 μK · 600 μs 转移 · 50 shots",
              "数值框架：Level 0 已验证组件 + 时变 velocity-Verlet · 扩展 12 维扫描 + 独立交叉验证 / 42,000+ shots"],
             size=14, color=RGBColor(0xC9, 0xD6, 0xEA))


# ---------------------------------------------------------------- 幻灯片 3：流程图
def slide_pipeline(prs):
    slide = content_slide(prs, "方法总览", "四步流水线：每一步计算什么")
    steps = [
        ("Step 0", "预检查", "积分器静态极限\n波形端点\n理想转移 + dt 收敛", "不出模拟先自证"),
        ("Step 1", "热初态采样", "简谐热分布 σx, σv\n拒绝法保证 E<0", "T = 5 μK"),
        ("Step 2", "轨迹推进", "时变 Verlet 积分\n合力 = SLM + AOD(x,t)", "600 μs, dt=0.05 μs"),
        ("Step 3", "末态判定", "E_f = ½mv²+U_SLM\nE_f < 0 → 俘获", "记录 capture margin"),
        ("Step 4", "统计汇总", "俘获率 + Wilson CI\n激发能变化分布", "× 50 shots"),
    ]
    left, top = Inches(0.45), Inches(1.6)
    width, height, gap = Inches(2.42), Inches(3.3), Inches(0.12)
    for i, (tag, name, body, note) in enumerate(steps):
        x = left + i * (width + gap)
        shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, x, top, width, height)
        shape.fill.solid()
        shape.fill.fore_color.rgb = LIGHT
        shape.line.color.rgb = STEP_COLORS[i]
        shape.line.width = Pt(2)
        shape.shadow.inherit = False
        header = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, x, top, width, Inches(0.62))
        header.fill.solid()
        header.fill.fore_color.rgb = STEP_COLORS[i]
        header.line.fill.background()
        header.shadow.inherit = False
        add_text(slide, x, top + Inches(0.07), width, Inches(0.5), f"{tag} · {name}",
                 size=15, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
        add_text(slide, x + Inches(0.12), top + Inches(0.8), width - Inches(0.24), Inches(1.9),
                 body.split("\n"), size=13, color=GRAY, line_spacing=1.25)
        add_text(slide, x + Inches(0.12), top + height - Inches(0.55), width - Inches(0.24), Inches(0.45),
                 note, size=12, color=STEP_COLORS[i], bold=True)
    add_text(slide, Inches(0.45), Inches(5.35), Inches(12.4), Inches(0.5),
             "Step 1–3 循环执行 50 次（固定 seed，可复现）；Step 0 在一切之前拦截数值风险。",
             size=14, color=NAVY, bold=True)
    add_text(slide, Inches(0.45), Inches(5.9), Inches(12.4), Inches(1.0),
             ["输出：mc_shots.csv（逐 shot 初/末态、能量、俘获裕量）、metrics.json、summary.md 与 12 张诊断图。"],
             size=13, color=GRAY)


# ---------------------------------------------------------------- 参数扫描三页
def slide_scan_framework(prs):
    slide = content_slide(prs, "扩展 · 参数扫描（一）框架", "把四步流水线变成可批量重放的实验台")
    bullet_runs(slide, [
        [("扫什么：", True, NAVY),
         ("12 个单维参数（时长、移动配比、波形形状、ramp 指数、间距、双阱深度、束腰失配、对准偏移、初始激发能、dt、温度）+ 深度二维、MC 统计与相位检验", False, GRAY)],
        [("怎么跑：", True, NAVY),
         ("每点独立 MC（198 shots），共享一个 multiprocessing Pool；扫描点 = 基线配置 + 单字段覆盖", False, GRAY)],
        [("两种初态模式：", True, NAVY),
         ("thermal（同主流程）或 fixed_excitation（固定激发能，转折点间均匀采 x，能量守恒定 |v|）", False, GRAY)],
        [("波形可插拔：", True, NAVY),
         ("5 种移动形状（linear / smoothstep / smootherstep / sine / ease-in-quad）× 任意 ramp 指数", False, GRAY)],
        [("每点产出：", True, NAVY),
         ("逐 shot CSV + summary.csv + 六联诊断图；总量 34,000+ shots（温度扫描独占 16,000）", False, GRAY)],
    ], Inches(0.55), Inches(1.45), Inches(4.9))
    add_picture_fit(slide, "scan_duration/scan_duration_detail.png",
                    Inches(5.6), Inches(1.4), Inches(7.2), Inches(5.3),
                    caption="示例：转移时长扫描的六联诊断图（每个扫描点各出一张）")


def slide_scan_dashboard(prs):
    slide = content_slide(prs, "扩展 · 参数扫描（二）总览", "12 维单参数扫描：除温度外俘获率全程 100%")
    bullet_runs(slide, [
        [("平台极宽：", True, NAVY),
         ("11 个维度上俘获率恒为 1.0（每点 198 shots，Wilson 下界 98.1%）", False, GRAY)],
        [("基线是内点：", True, NAVY),
         ("绿线标注的默认配置落在平台的中心，不是精细调参的孤立解", False, GRAY)],
        [("数值已收敛：", True, NAVY),
         ("dt 从 0.05 → 0.0125 μs 结果不变 → 平台不是步长伪影", False, GRAY)],
        [("唯一敏感维度：", True, STEP_COLORS[3]),
         ("温度（右下图），交给下一页详查", False, GRAY)],
    ], Inches(0.55), Inches(1.5), Inches(4.0))
    add_picture_fit(slide, "scan_dashboard.png", Inches(4.75), Inches(1.3), Inches(8.2), Inches(5.5),
                    caption="12 个扫描维度的俘获率总览（3 × 4 面板，绿色竖线 = 基线点）")


def slide_scan_boundaries(prs):
    slide = content_slide(prs, "扩展 · 参数扫描（三）边界", "三个真正的边界：温度、能量悬崖、加热区")
    trio = [
        ("scan_sensitivity_ranking.png", "敏感度排序：各维度俘获率极差"),
        ("scan_temperature/capture_fraction_vs_temperature.png", "温度边界：20 μK 后开始损失"),
        ("scan_energy_boundary/energy_boundary_scatter.png", "能量悬崖：180 → 220 μK 崩塌"),
    ]
    for i, (name, caption) in enumerate(trio):
        add_picture_fit(slide, name, Inches(0.35) + i * Inches(4.28), Inches(1.35),
                        Inches(4.15), Inches(2.75), caption=caption)
    bullet_runs(slide, [
        [("温度边界：", True, NAVY),
         ("T ≤ 20 μK 全俘获；30 μK 起下降（99.7%）；250 μK → 76.9%（999 shots/点）", False, GRAY)],
        [("能量悬崖：", True, STEP_COLORS[3]),
         ("初始激发 180 μK 仍全俘获，220 μK 崩塌至 1.5% —— 判据本质是初激发 ≪ SLM 深度", False, GRAY)],
        [("加热区（俘获仍 100%，但 dE > 0）：", True, NAVY),
         ("SLM 过深 320 μK → +55 μK；AOD 过浅 140 μK → +12 μK；间距过近 0.8 μm → +7.8 μK", False, GRAY)],
        [("工程含义：", True, NAVY),
         ("安全工作区 = 温度 ≪ 30 μK 且初激发 ≪ SLM 深度，同时避开过深 SLM / 过浅 AOD / 过近间距", False, GRAY)],
    ], Inches(0.55), Inches(4.55), Inches(12.4))


def slide_scan_crosscheck(prs):
    """独立副本 run_parameter_scans.py 的 6 维扫描：交叉验证 + 快转移盲区。"""
    def read_rows(name):
        with open(CROSSCHECK / name, newline="") as handle:
            rows = []
            for raw in csv.DictReader(handle):
                row = {}
                for key, value in raw.items():
                    try:
                        row[key] = float(value)
                    except (TypeError, ValueError):
                        row[key] = value
                rows.append(row)
            return rows

    duration = sorted(read_rows("scan_duration_us.csv"), key=lambda r: r["value"])
    meta = json.loads((CROSSCHECK / "scan_meta.json").read_text())
    fastest = duration[0]
    first_full = next(r for r in duration if r["capture_fraction"] >= 0.999)
    plateau = [r["mean_excitation_change_uK"] for r in duration if r["value"] >= 200]
    plateau_mean = sum(plateau) / len(plateau)
    aod = {r["value"]: r for r in read_rows("scan_aod_depth_uK.csv")}
    slm = {r["value"]: r for r in read_rows("scan_slm_depth_uK.csv")}
    sep = {r["value"]: r for r in read_rows("scan_initial_center_um.csv")}

    slide = content_slide(prs, "扩展 · 参数扫描（四）交叉验证",
                          "独立副本复现同一结论，并补上快转移盲区")
    bullet_runs(slide, [
        [("是什么：", True, NAVY),
         ("另一份独立代码副本中的 run_parameter_scans.py：6 个单参数扫描"
          "（时长/配比/双阱深度/间距/波形形状）+ 温度合并面板，OAT 设计与本套件相同", False, GRAY)],
        [("怎么跑：", True, NAVY),
         (f"每点 {int(meta['shots_per_point'])} shots（8×25 并行），41 点共 {meta['total_shots']:,} shots、"
          f"{meta['elapsed_s']:.0f} s（{meta['n_procs']} 进程）；误差棒 = Wilson 95% CI，虚线 = 基线", False, GRAY)],
        [("新发现 · 快转移崩塌：", True, STEP_COLORS[3]),
         (f"{fastest['value']:.0f} μs → 俘获率 {fastest['capture_fraction']:.1%}"
          f"（{int(fastest['captured'])}/{int(fastest['shots'])}）、平均加热 "
          f"{fastest['mean_excitation_change_uK']:+.0f} μK；{first_full['value']:.0f} μs → 恢复 100%", False, GRAY)],
        [("   本套件时长扫描最快只到 100 μs，独立副本补上了 <100 μs 盲区", False, GRAY)],
        [("重叠区间互相印证：", True, NAVY),
         (f"≥{first_full['value']:.0f} μs 后两套结果一致（全俘获，平台 ΔE ≈ {plateau_mean:+.1f} μK）；"
          f"加热趋势也相同：浅 AOD {aod[100.0]['mean_excitation_change_uK']:+.0f} μK@100μK、"
          f"深 SLM {slm[280.0]['mean_excitation_change_uK']:+.0f} μK@280μK、"
          f"近间距 {sep[0.6]['mean_excitation_change_uK']:+.1f} μK@0.6μm", False, GRAY)],
        [("安全工作区新增下界：", True, NAVY), ("转移时长 ≳ 75 μs", True, STEP_COLORS[1])],
    ], Inches(0.55), Inches(1.45), Inches(6.0))
    add_picture_fit(slide, "scan_parameters/scan_duration_us.png",
                    Inches(6.8), Inches(1.4), Inches(6.1), Inches(5.2),
                    caption="独立副本时长扫描：50 μs 崩塌、75 μs 起全俘获（上：俘获率+CI，下：ΔE）")


def build():
    metrics = json.loads(METRICS.read_text())
    mc = metrics["monte_carlo"]
    energy = metrics["energy"]
    ideal = metrics["ideal"]

    prs = Presentation()
    prs.slide_width, prs.slide_height = SLIDE_W, SLIDE_H

    slide_cover(prs)

    # P2 物理场景
    slide = content_slide(prs, "物理场景", "模拟什么：一个原子、两张光镊、一次交接")
    bullet_runs(slide, [
        [("原子：", True, NAVY), ("Cs-133，经典质点，一维径向坐标 x", False, GRAY)],
        [("AOD 光镊（移动方）：", True, NAVY), ("初始位于 x = 2.4 μm，深度 280 μK", False, GRAY)],
        [("SLM 光镊（接收方）：", True, NAVY), ("静止于 x = 0，深度 140 μK", False, GRAY)],
        [("两阱同为高斯型：", True, NAVY), ("U(x) = −U₀ exp[−2(x−xc)²/w²]，w = 1.17 μm", False, GRAY)],
        [("任务：", True, NAVY), ("AOD 携带原子平移合并到 SLM 后释放，问原子最终是否留在 SLM", False, GRAY)],
    ], Inches(0.55), Inches(1.4), Inches(4.7))
    add_picture_fit(slide, "01_setup_two_traps.png", Inches(5.5), Inches(1.35), Inches(7.3), Inches(5.3),
                    caption="初始构型：间距 2.4 μm 的双高斯阱")

    # P3 流程图
    slide_pipeline(prs)

    # P4 Step 0 预检查
    slide = content_slide(prs, "Step 0 · 预检查", "先验证工具，再相信结果：三项前置检查全部通过")
    bullet_runs(slide, [
        [("① 静态极限回归：", True, NAVY),
         ("时变积分器退化为 F(x) 时，与 Level 0 已验证的 velocity-Verlet 轨迹逐点一致", False, GRAY)],
        [("② 波形端点检查：", True, NAVY),
         ("t=0 时 AOD 原位满深度；t=T 时到达 SLM 且深度为零", False, GRAY)],
        [("③ 理想转移 + 步长收敛：", True, NAVY),
         ("原子从阱底静止出发，dt 与 dt/2 两次运行", False, GRAY)],
        [("   末能差 / U₀ = ", False, GRAY), ("1.0 × 10⁻⁹", True, STEP_COLORS[3]),
         ("，俘获判定一致 → 数值误差远小于物理效应", False, GRAY)],
        [("任何一项失败：", True, NAVY), ("程序报错停止，不进入蒙特卡洛", False, GRAY)],
    ], Inches(0.55), Inches(1.5), Inches(6.1))
    add_picture_fit(slide, "12_validation_dashboard.png", Inches(6.9), Inches(1.4), Inches(6.0), Inches(5.2),
                    caption="七项自动验证面板（含 MC 可复现性）")

    # P5 Step 1 采样
    slide = content_slide(prs, "Step 1 · 热初态采样", "每个 shot 从有限温度的束缚态出发")
    bullet_runs(slide, [
        [("简谐近似热分布：", True, NAVY), ("x₀ ~ N(xc, σx²)，v₀ ~ N(0, σv²)", False, GRAY)],
        [("   σx = √(kBT/mω²) = ", False, GRAY), ("0.078 μm", True, STEP_COLORS[1]),
         ("，  σv = √(kBT/m) = ", False, GRAY), ("0.0177 m/s", True, STEP_COLORS[1]),
         ("（T = 5 μK，ω 取 AOD 阱频）", False, GRAY)],
        [("拒绝采样：", True, NAVY), ("仅当 AOD 静态势中总能量 E₀ < 0（束缚）才接受", False, GRAY)],
        [("   1000 样本抽查全部束缚；固定 seed = 20260811，逐 shot 可复现", False, GRAY)],
        [("物理含义：", True, NAVY), ("初态弥散 ≈ 0.08 μm ≪ 2.4 μm 转移距离，原子确定地在 AOD 阱内", False, GRAY)],
    ], Inches(0.55), Inches(1.5), Inches(6.0))
    add_picture_fit(slide, "06_thermal_distribution.png", Inches(6.8), Inches(1.4), Inches(6.1), Inches(5.2),
                    caption="初态采样分布：位置与速度（含束缚判定边界）")

    # P6 Step 2a 波形
    slide = content_slide(prs, "Step 2 · 轨迹推进（a）", "两段式波形：先移动合并，再绝热释放")
    bullet_runs(slide, [
        [("移动段（0 – 312 μs，52%）：", True, NAVY),
         ("AOD 中心 2.4 μm → 0", False, GRAY)],
        [("   smoothstep h(s) = 3s² − 2s³：起点/终点速度与加速度均为零，避免踢原子", False, GRAY)],
        [("释放段（312 – 600 μs，48%）：", True, NAVY),
         ("AOD 深度 280 μK → 0", False, GRAY)],
        [("   按 (1−s)² 二次下降：后段变化最慢，降低释放末端的激发", False, GRAY)],
        [("SLM 全程静止 140 μK", True, NAVY), ("，只通过势场叠加“接手”原子", False, GRAY)],
    ], Inches(0.55), Inches(1.5), Inches(6.1))
    add_picture_fit(slide, "02_waveform_schedule.png", Inches(6.9), Inches(1.4), Inches(6.0), Inches(5.2),
                    caption=" engineered 波形：AOD 中心与深度的时间表")

    # P7 Step 2b 推进
    slide = content_slide(prs, "Step 2 · 轨迹推进（b）", "时变 velocity-Verlet 在合力下积分整条轨迹")
    bullet_runs(slide, [
        [("合力构成：", True, NAVY), ("F(x,t) = F_SLM(x) + F_AOD(x,t)，两张高斯阱实时叠加", False, GRAY)],
        [("积分参数：", True, NAVY), ("dt = 0.05 μs，12001 步覆盖 600 μs", False, GRAY)],
        [("   步长合法性：ωmax·dt < 0.05（配置加载时强制）", False, GRAY)],
        [("注意：", True, STEP_COLORS[3]),
         ("外场显式随时间变化，机械能不守恒是物理而非数值误差", False, GRAY)],
        [("理想原子（阱底静止出发）：", True, NAVY),
         ("全程跟随 AOD 移动并入 SLM，末态 E/U₀ = −0.9999997", False, GRAY)],
    ], Inches(0.55), Inches(1.5), Inches(6.0))
    add_picture_fit(slide, "04_ideal_position.png", Inches(6.9), Inches(1.4), Inches(6.0), Inches(5.2),
                    caption="理想转移：原子位置 vs AOD 中心随时间演化")

    # P8 Step 3 判定
    slide = content_slide(prs, "Step 3 · 末态判定", "一个能量判据二分所有结局：俘获或逃逸")
    bullet_runs(slide, [
        [("末态能量：", True, NAVY), ("E_f = ½mv_f² + U_SLM(x_f)（以 SLM 静态势衡量）", False, GRAY)],
        [("判据：", True, NAVY), ("E_f < 0 → 俘获；E_f ≥ 0 → 逃逸", False, GRAY)],
        [("俘获裕量：", True, NAVY), ("margin = −E_f / U₀^SLM，越大越“稳”", False, GRAY)],
        [("   margin ≥ 1 表示动能已不足以爬出 SLM 阱口", False, GRAY)],
        [("同时记录激发能变化：", True, NAVY), ("ΔE_exc = E_f,exc − E₀,exc，衡量转移注入/带走的运动能量", False, GRAY)],
    ], Inches(0.55), Inches(1.5), Inches(6.0))
    add_picture_fit(slide, "09_capture_margin_per_shot.png", Inches(6.9), Inches(1.4), Inches(6.0), Inches(5.2),
                    caption="逐 shot 俘获裕量：50/50 全部深度俘获")

    # P9 Step 4a 统计
    slide = content_slide(prs, "Step 4 · 统计（a）", f"50 shots 结果：经典俘获率 {mc['classical_capture_fraction']:.0%}")
    low, high = mc["wilson_95_interval"]
    bullet_runs(slide, [
        [("俘获计数：", True, NAVY), (f"{mc['captured_shots']} / {mc['shots']} shots 全部俘获", False, GRAY)],
        [("经典俘获率：", True, NAVY), (f"{mc['classical_capture_fraction']:.1%}", True, STEP_COLORS[1])],
        [("Wilson 95% 置信区间：", True, NAVY),
         (f"[{low:.1%}, {high:.1%}]（小样本下界即由此而来）", False, GRAY)],
        [("理想原子对照：", True, NAVY),
         (f"末态 E/U₀ = {ideal['final_energy_over_slm_depth']:.8f}，激发变化仅 "
          f"{ideal['excitation_energy_change_uK']:.0e} μK", False, GRAY)],
        [("解读：", True, NAVY), ("该波形下转移对初温 5 μK 的经典热样品是“无损”的", False, GRAY)],
    ], Inches(0.55), Inches(1.5), Inches(6.0))
    add_picture_fit(slide, "11_capture_rate_summary.png", Inches(6.9), Inches(1.4), Inches(6.0), Inches(5.2),
                    caption="俘获率汇总与置信区间")

    # P10 Step 4b 激发能
    slide = content_slide(prs, "Step 4 · 统计（b）", "转移不是加热而是冷却：激发能平均下降 ~2 μK")
    bullet_runs(slide, [
        [("激发能变化：", True, NAVY),
         (f"均值 {energy['mean_excitation_change_uK']:+.2f} μK，中位数 {energy['median_excitation_change_uK']:+.2f} μK", False, GRAY)],
        [("末态能量：", True, NAVY),
         (f"平均 E_f/U₀ = {energy['mean_final_energy_over_slm_depth']:.3f}（接近阱底）", False, GRAY)],
        [("机制：", True, NAVY), ("合并减速 + 释放段 (1−s)² 绝热性 → 相空间收缩", False, GRAY)],
        [("对照初温：", True, NAVY), ("初态平均激发 6.6 μK → 末态 4.6 μK", False, GRAY)],
    ], Inches(0.55), Inches(1.45), Inches(4.6))
    add_picture_fit(slide, "08_phase_space_collapse.png", Inches(0.45), Inches(4.0), Inches(6.2), Inches(3.1))
    add_picture_fit(slide, "10_excitation_change_hist.png", Inches(6.8), Inches(1.4), Inches(6.1), Inches(5.6),
                    caption="逐 shot 激发能变化分布（负值 = 冷却）")

    # P11-P14 参数扫描扩展（含独立副本交叉验证）
    slide_scan_framework(prs)
    slide_scan_dashboard(prs)
    slide_scan_boundaries(prs)
    slide_scan_crosscheck(prs)

    # P14 结论
    slide = content_slide(prs, "结论", "四步流水线 + 参数扫描给出的四个答案")
    conclusions = [
        ("数值可信", "静态极限、波形端点、dt 收敛（10⁻⁹）、MC 可复现全部通过——结果不受数值伪影支配"),
        ("俘获完全", "50/50 经典俘获，Wilson 95% CI 下界 92.9%；理想情形激发变化 ~5×10⁻⁵ μK"),
        ("过程是冷却的", "热样品平均激发能 6.6 → 4.6 μK（−1.9 μK），波形设计对初始热误差有鲁棒性"),
        ("工作区很宽", "12+6 维扫描（两套独立实现）中除温度与快转移外俘获率恒 100%；边界清晰：时长 ≳ 75 μs、温度 ≪ 30 μK、初激发 ≪ SLM 深度"),
    ]
    top = Inches(1.55)
    for i, (title, body) in enumerate(conclusions):
        y = top + i * Inches(1.7)
        badge = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.7), y, Inches(0.85), Inches(0.85))
        badge.fill.solid()
        badge.fill.fore_color.rgb = STEP_COLORS[i]
        badge.line.fill.background()
        badge.shadow.inherit = False
        add_text(slide, Inches(0.7), y + Inches(0.14), Inches(0.85), Inches(0.6), str(i + 1),
                 size=24, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
        add_text(slide, Inches(1.85), y - Inches(0.02), Inches(2.6), Inches(0.5), title,
                 size=19, color=NAVY, bold=True)
        add_text(slide, Inches(1.85), y + Inches(0.45), Inches(10.9), Inches(1.0), body,
                 size=15, color=GRAY, line_spacing=1.25)

    # P12 边界
    slide = content_slide(prs, "边界与局限", "这个数字是什么、不是什么")
    bullet_runs(slide, [
        [("是：", True, STEP_COLORS[1]),
         ("一维经典力学模型在指定配置下的俘获率与激发能统计", False, GRAY)],
        [("不是：", True, STEP_COLORS[3]),
         ("实验存活率、AOD-SLM 转移保真度或量子态保真度", False, GRAY)],
        [("   不可直接对标文献 99.81% 基准——物理内容不同", False, GRAY)],
        [("未包含：", True, NAVY),
         ("硬件噪声、2D/3D 与轴向运动、光子散射、真空碰撞、多原子相互作用、量子内态与相干性", False, GRAY)],
        [("参数假设：", True, NAVY),
         ("转移时长 600 μs、温度 5 μK、AOD 束腰 1.17 μm 均为 assumed，已在配置中显式标记", False, GRAY)],
        [("下一步：", True, NAVY),
         ("单参数扫描已完成（两套实现互证）；可扩展方向：二维联合扫描（深度×时长）、噪声模型与多原子相互作用", False, GRAY)],
    ], Inches(0.55), Inches(1.5), Inches(12.2))
    add_picture_fit(slide, "03_potential_evolution.png", Inches(0.45), Inches(4.15), Inches(12.4), Inches(3.0))

    prs.save(OUTPUT)
    return OUTPUT, len(prs.slides._sldIdLst)


if __name__ == "__main__":
    path, count = build()
    print(f"Saved {count} slides -> {path}")
