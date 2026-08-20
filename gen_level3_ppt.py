# -*- coding: utf-8 -*-
"""生成 Level 3 总结 PPT：用 outputs/level3_3d_transport_lensing/demo 的 12 张图
+ metrics/summary 数据要点组装幻灯片。运行：
.venv/bin/python gen_level3_ppt.py
"""
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

REPO = Path(__file__).resolve().parent
FIG = REPO / "simulation" / "level2_joint_transfer" / "outputs" / \
    "level3_3d_transport_lensing" / "demo"
OUT = REPO / "Level3_三维AOD长距离运输与柱面透镜_总结报告.pptx"

BLUE = RGBColor(0x1F, 0x77, 0xB4)
DARK = RGBColor(0x22, 0x2A, 0x35)
RED = RGBColor(0xD6, 0x27, 0x28)
GREEN = RGBColor(0x2C, 0xA0, 0x2C)
ORANGE = RGBColor(0xE8, 0x7D, 0x0E)
GRAY = RGBColor(0x66, 0x66, 0x66)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT = RGBColor(0xDD, 0xE8, 0xF5)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]


def add_slide():
    return prs.slides.add_slide(BLANK)


def set_text(box, text, size=16, color=DARK, bold=False, align=PP_ALIGN.LEFT):
    lines = text.split("\n")
    box.text = lines[0]
    run = box.text_frame.paragraphs[0].runs[0]
    run.font.size = Pt(size)
    run.font.color.rgb = color
    run.font.bold = bold
    for line in lines[1:]:
        paragraph = box.text_frame.add_paragraph()
        run = paragraph.add_run()
        run.text = line
        run.font.size = Pt(size)
        run.font.color.rgb = color
        run.font.bold = bold
        paragraph.alignment = align
    box.text_frame.paragraphs[0].alignment = align
    return box.text_frame.paragraphs[0]


def add_textbox(slide, x, y, w, h, text, size=16, color=DARK, bold=False,
                align=PP_ALIGN.LEFT):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    box.text_frame.word_wrap = True
    return set_text(box, text, size, color, bold, align)


def add_bullets(slide, x, y, w, h, items, size=15, gap=Pt(6)):
    """items: (head, body) 列表，head 加粗、body 灰色小一号。"""
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.word_wrap = True
    for i, (head, body) in enumerate(items):
        paragraph = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        paragraph.space_after = gap
        run = paragraph.add_run()
        run.text = f"• {head}"
        run.font.size = Pt(size)
        run.font.bold = True
        run.font.color.rgb = DARK
        if body:
            run2 = paragraph.add_run()
            run2.text = f"　{body}"
            run2.font.size = Pt(size - 1)
            run2.font.color.rgb = GRAY
    return box


def title_bar(slide, title, subtitle=None):
    bar = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width,
                                 Inches(0.95))
    bar.fill.solid()
    bar.fill.fore_color.rgb = BLUE
    bar.line.fill.background()
    add_textbox(slide, 0.5, 0.12, 12.4, 0.6, title, size=26, color=WHITE,
                bold=True)
    if subtitle:
        add_textbox(slide, 0.55, 0.62, 12.3, 0.35, subtitle, size=13,
                    color=LIGHT)


def add_picture_fit(slide, path: Path, x, y, w, h):
    """等比缩放居中放入 (x,y,w,h) 框（英寸）。"""
    from PIL import Image
    with Image.open(path) as image:
        iw, ih = image.size
    scale = min(w / iw, h / ih)
    nw, nh = iw * scale, ih * scale
    left = Inches(x + (w - nw) / 2)
    top = Inches(y + (h - nh) / 2)
    return slide.shapes.add_picture(str(path), left, top, Inches(nw), Inches(nh))


def add_table(slide, x, y, w, rows, col_widths=None, size=12, header=True):
    """rows: list[list[str]]；第一行表头蓝底白字。返回 table shape。"""
    n_r, n_c = len(rows), len(rows[0])
    shape = slide.shapes.add_table(n_r, n_c, Inches(x), Inches(y), Inches(w),
                                   Inches(0.32 * n_r))
    table = shape.table
    if col_widths:
        for j, cw in enumerate(col_widths):
            table.columns[j].width = Inches(cw)
    for i, row in enumerate(rows):
        for j, text in enumerate(row):
            cell = table.cell(i, j)
            cell.text = str(text)
            for paragraph in cell.text_frame.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(size)
                    if header and i == 0:
                        run.font.bold = True
                        run.font.color.rgb = WHITE
                    else:
                        run.font.color.rgb = DARK
            if header and i == 0:
                cell.fill.solid()
                cell.fill.fore_color.rgb = BLUE
            elif i % 2 == 0:
                cell.fill.solid()
                cell.fill.fore_color.rgb = RGBColor(0xF0, 0xF4, 0xF9)
    return shape


def fig_caption(slide, text, y=6.95):
    add_textbox(slide, 0.55, y, 12.3, 0.4, text, size=11, color=GRAY)


# ================================================================ 1 封面
slide = add_slide()
bg = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width,
                            prs.slide_height)
bg.fill.solid()
bg.fill.fore_color.rgb = RGBColor(0x0E, 0x2A, 0x47)
bg.line.fill.background()
add_textbox(slide, 1.0, 1.9, 11.3, 0.7,
            "Monte Carlo for AOD→SLM 原子转移 · 论文 arXiv:2403.12021v4",
            size=22, color=RGBColor(0x9E, 0xC5, 0xE8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 2.8, 11.3, 1.6,
            "Level 3：三维 AOD 长距离运输\n与柱面透镜效应", size=40,
            color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 4.9, 11.3, 0.9,
            "三维经典模型 · crossed-AOD 柱面透镜势 · 270–610 μm 运输\n"
            "直线 vs 对角 · 三类轨迹 · 热初态蒙特卡洛 · 独立 2000-shot 验证",
            size=16, color=RGBColor(0xBB, 0xD3, 0xE8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 6.5, 11.3, 0.5,
            "preflight 16 项全通过 · 96 个粗扫描格点 · 22 项 Level 3 测试",
            size=13, color=RGBColor(0x88, 0xA8, 0xC8), align=PP_ALIGN.CENTER)

# ================================================================ 2 背景与定位
slide = add_slide()
title_bar(slide, "背景与 Level 3 的定位",
          "在 Level 0 静态势阱 / Level 1 顺序转移 / Level 2 同步波形基础上隔离研究长距离运输")
add_bullets(slide, 0.55, 1.25, 6.1, 5.5, [
    ("论文锚点", "280 μK AOD 光镊中 610 μm、1.6 ms 对角运输；直线运输受柱面透镜"
     "效应限制（Fig. 5、ED Fig. 10、Methods、SI §III）"),
    ("核心问题", "有限温度铯原子在 270–610 μm 运输中，方向、轨迹形状、时长和"
     "透镜强度如何影响经典留阱率、末态激发与轴向运动"),
    ("隔离研究", "暂不把 transfer 与 long move 组合，先回答『AOD 单阱已装原子后"
     "如何长距离移动』"),
    ("诚实边界", "经典留阱 = 末态仍被束缚，不含真空损失、光子散射、退相干；"
     "不得与论文 99.953(2)% 量子保真度等同"),
])
add_table(slide, 6.95, 1.4, 5.9, [
    ["物理量", "数值"],
    ["原子", "¹³³Cs，m = 132.905 u"],
    ["AOD 阱深 D", "k_B × 280 μK"],
    ["束腰 w₀", "1.17 μm（假设匹配 SLM 后孔径）"],
    ["波长 λ", "1055 nm"],
    ["Rayleigh 长度 z_R", "4.076 μm"],
    ["ω_r / 2π", "36.0 kHz"],
    ["ω_z / 2π", "7.31 kHz"],
    ["初态温度", "5 μK（D/kT = 56）"],
], col_widths=[2.3, 3.6], size=13)
add_textbox(slide, 6.95, 4.9, 5.9, 1.6,
            "留阱判据：运输终点静态势中 E_f < 0 → retained；\n"
            "100 μs 保持后二次判定并统计不一致数。", size=13, color=GRAY)

# ================================================================ 3 静态势
slide = add_slide()
title_bar(slide, "三维高斯光镊势：模型底座",
          "U₀ = −D/(1+(z/z_R)²) · exp(−2(ξ²+η²)/(w₀²(1+(z/z_R)²)))，z=0 切片严格回归 Level 0")
add_picture_fit(slide, FIG / "static_3d_trap_slices.png", 0.8, 1.15, 11.7, 5.3)
add_bullets(slide, 0.55, 6.3, 12.3, 0.9, [
    ("数值曲率复现解析频率", "ω_r、ω_z 相对误差 < 10⁻⁶（preflight 第 3 项）；"
     "三维解析力与势梯度有限差分一致到 10⁻⁵"),
], size=13)

# ================================================================ 4 透镜模型
slide = add_slide()
title_bar(slide, "Crossed-AOD 柱面透镜势（论文 SI 启发）",
          "快速扫频使柱面焦点沿光轴偏移：z_s/z_R = 2σ·ċ/v_s，U = −D·G₁(ξ,z;z_s₁)·G₂(η,z;z_s₂)")
add_bullets(slide, 0.55, 1.2, 6.2, 4.2, [
    ("三种校准模式", "lensing_off（理想）/ absolute_critical_velocity（用户给定 v_s）"
     "/ normalized_reference_scan（默认）"),
    ("无量纲 severity 扫描", "以 610 μm、1.6 ms 对角 sine 为参考，扫描"
     " v_max,axis/v_s ∈ {0.50, 0.75, 1.00, 1.25}，不伪造论文未提供的硬件校准"),
    ("同一硬件含义", "每个 severity 反解一个固定 v_s，用于所有距离/方向/时长/轨迹，"
     "不为每个候选重新缩放"),
    ("零偏移极限", "z_s₁ = z_s₂ = 0 时严格恢复标准三维高斯势（数值验证 rtol<10⁻¹²）"),
])
add_table(slide, 6.95, 1.4, 5.9, [
    ["severity v/v_s", "固定 v_s (m/s)"],
    ["0.50", "1.078"],
    ["0.75", "0.719"],
    ["1.00", "0.539"],
    ["1.25", "0.431"],
], col_widths=[3.0, 2.9], size=13)
add_textbox(slide, 6.95, 3.4, 5.9, 3.0,
            "对角同号（σ₁=σ₂=1）与论文『对角运动把柱面透镜\n变为球面透镜』一致；反号 (−1,+1)"
            " 仅作模型敏感性。\n\n默认轴符号写入配置 (axis_signs: [1,1])，\n"
            "severity 为无量纲场景参数，非论文测量值。", size=13, color=GRAY)

# ================================================================ 5 分岔
slide = add_slide()
title_bar(slide, "单阱 → 双阱分岔：直线运输的物理瓶颈",
          "直线：|z_s| > 2 z_R 时轴向势分裂为双阱（临界 ≈ 单轴速度达 v_s）；对角同号：保持单一球面阱")
add_picture_fit(slide, FIG / "lensing_potential_bifurcation.png", 0.8, 1.15,
                11.7, 5.0)
add_bullets(slide, 0.55, 6.25, 12.3, 0.9, [
    ("数值极值诊断", "axial_minima 在自适应 4001 点网格上判最低点数/位置/势垒；"
     "临界区 1.8z_R 单阱 → 2.2z_R 双阱，与解析关系互为回归测试"),
], size=13)

# ================================================================ 6 轨迹族
slide = add_slide()
title_bar(slide, "三类运输轨迹族",
          "adiabatic-sine（论文 Methods）· 四段 constant-jerk（J=32L/T³，分段解析积分）· minimum-jerk")
add_picture_fit(slide, FIG / "trajectory_profile_comparison.png", 0.8, 1.15,
                11.7, 4.6)
add_bullets(slide, 0.55, 5.9, 12.3, 1.3, [
    ("端点约束", "三类均满足起点/终点精确、端点速度与加速度为零、位置单调无过冲；"
     "constant-jerk 符号模式 +J/−J/−J/+J 与 J=32L/T³ 数值验证"),
    ("命名诚实", "minimum-jerk (10s³−15s⁴+6s⁵) 明确标注为通用平滑控制轨迹，"
     "不冒充论文 constant-jerk；相同 L/T 下用同一初态样本比较"),
], size=13)

# ================================================================ 7 straight vs diagonal
slide = add_slide()
title_bar(slide, "直线 vs 对角：几何决定透镜强度",
          "相同路径长度 L 与时长 T 下：对角把单轴速度降为 1/√2，两柱面焦点同移；直线单轴满速 → 焦点偏移最大")
add_picture_fit(slide, FIG / "straight_vs_diagonal.png", 0.8, 1.15, 11.7, 5.6)
fig_caption(slide, "左：轴 1 峰值速度 / v_s；中：峰值 |z_s₁|/z_R；右：瞬时最低点深度（μK）"
                   "——straight 在同一 severity 下焦点偏移更大、有效阱深更浅")

# ================================================================ 8 数值方法
slide = add_slide()
title_bar(slide, "数值方法与验证流水线",
          "复用 Level 0–2 的统计、I/O 与种子管理；16 项 preflight 失败即中止大规模 MC")
add_bullets(slide, 0.55, 1.2, 6.1, 5.6, [
    ("三维时变 velocity-Verlet", "向量化推进整批 shots；显式时变势下二阶精度，"
     "不声称辛能量守恒"),
    ("显式功-能核算", "∂U/∂t 解析链式导数分解为 移动项 + 焦点偏移项 + 深度项，"
     "梯形积分 W_ext，检查残差 R_W = E−E₀−W"),
    ("热初态采样", "简谐提议 + E_i<0 束缚拒绝（D/kT=56 时接受率≈1）；"
     "明确非完整有限深阱正则系综"),
    ("三池种子隔离", "scan 31001 / validation 31002 / convergence 31003 互不相同；"
     "条件间 common random numbers 配对比较"),
    ("分阶段预算", "A 确定性扫描 → B 粗 MC(128 shots, 0.10 μs) → C 冻结验证"
     "(2000 shots, 0.05 μs) → D 敏感性(500 shots)"),
])
add_table(slide, 6.95, 1.4, 5.9, [
    ["preflight 代表项", "结果"],
    ["Level 0–2 全部测试", "通过"],
    ["三维势回归 Level 0 切片", "rtol < 10⁻¹²"],
    ["力 vs 有限差分梯度", "相对误差 < 10⁻⁵"],
    ["功-能残差随步长收敛", "达标（0.05/0.025 μs）"],
    ["128 热初态步长配对", "标签翻转率 ≤ 2%"],
    ["100 μs 保持段能量误差", "< 10⁻³ 阱深"],
], col_widths=[3.6, 2.3], size=13)
add_textbox(slide, 6.95, 4.6, 5.9, 2.4,
            "冻结主条件（打开 validation 前锁定）：\n"
            "① paper_anchor：610 μm / 1.6 ms / 对角 / sine\n"
            "② 同条件 constant-jerk\n"
            "③ straight_control：510 μm / 1.6 ms / 直线 / sine\n"
            "④ 同条件 constant-jerk\n"
            "每个条件 × {lensing_off, sev1.0}", size=13, color=GRAY)

# ================================================================ 9 粗扫描相图
slide = add_slide()
title_bar(slide, "阶段 B：粗 MC 相图（96 格点 × 128 shots）",
          "留阱率 × 时长 × severity：直线在短时长 + 高 severity 下崩塌，对角全区间存活")
add_picture_fit(slide, FIG / "survival_phase_diagram.png", 0.55, 1.1, 7.4, 5.4)
add_textbox(slide, 8.2, 1.3, 4.7, 5.2,
            "关键数字（adiabatic-sine）：\n\n"
            "● 直线 510 μm / sev1.0\n"
            "   600–800 μs：留阱率 0.00\n"
            "   1000 μs：0.05 → 1200 μs 起 1.00\n\n"
            "● 对角 610 μm / sev1.25\n"
            "   600 μs 仍 0.93，≥800 μs 全 1.00\n\n"
            "● 无透镜（sev0）所有格点 1.00\n\n"
            "结论：柱面透镜是短时长直线运输的\n唯一损失通道；1.6 ms（论文锚点）\n"
            "下所有已测条件均安全。", size=14, color=DARK)

# ================================================================ 10 轴向动力学
slide = add_slide()
title_bar(slide, "轴向动力学与 lost 轨迹解剖",
          "左：验证主条件 retained 代表 z(t)；右：粗扫描失败格点（直线/510 μm/600 μs/sev1.25）的 lost 轨迹与瞬时最低点分支")
add_picture_fit(slide, FIG / "axial_dynamics.png", 0.8, 1.15, 11.7, 5.0)
add_bullets(slide, 0.55, 6.25, 12.3, 1.0, [
    ("分支诊断", "虚线为瞬时轴向最低点位置：单阱→双阱→单阱随速度演化；"
     "lost 轨迹用独立诊断种子（scan+9100）捕获，|mean_z|<0.1z_R 标 branch_ambiguous 不强行归类"),
], size=13)

# ================================================================ 11 冻结验证
slide = add_slide()
title_bar(slide, "阶段 C：冻结主条件独立验证（2000 shots × 8 条件）",
          "validation_pool 一次性使用；同初态配对比较 + 0.025 μs 子集复算")
add_table(slide, 0.55, 1.2, 7.6, [
    ["条件", "sev0.0", "sev1.0"],
    ["对角 610 μm · sine（论文锚点）", "2000/2000", "2000/2000"],
    ["对角 610 μm · constant-jerk", "2000/2000", "2000/2000"],
    ["直线 510 μm · sine", "2000/2000", "2000/2000"],
    ["直线 510 μm · constant-jerk", "2000/2000", "2000/2000"],
], col_widths=[4.0, 1.8, 1.8], size=13)
add_textbox(slide, 0.55, 3.3, 7.6, 1.6,
            "所有条件 Wilson 95% CI = [0.9981, 1.0000]；四组配对差异均为 0\n"
            "（sine vs cjerk、straight vs diagonal、lensing off vs sev1.0，均不显著）。\n"
            "1.6 ms 时长下模型处于深度安全区，差异需在更短时长才能分辨——与粗扫描一致。",
            size=14, color=DARK)
add_picture_fit(slide, FIG / "heating_comparison.png", 8.35, 1.2, 4.6, 3.0)
add_picture_fit(slide, FIG / "final_energy_distributions.png", 8.35, 4.3,
                4.6, 2.6)
fig_caption(slide, "右上图：retained shots 的 Δε 与末态激发分布；右下图：末态静态能量分布与 E_f=0 阈值")

# ================================================================ 12 加热标度律
slide = add_slide()
title_bar(slide, "加热标度律：一项与任务书不一致的诚实发现",
          "简谐、无透镜、小激发极限；干涉振荡取上包络拟合")
add_picture_fit(slide, FIG / "heating_scaling_trends.png", 0.55, 1.15, 8.0, 5.3)
add_bullets(slide, 8.8, 1.4, 4.3, 5.0, [
    ("实测斜率", "sine −6.64 / constant-jerk −6.09 / minimum-jerk −6.60；"
     "L 斜率均 +2.00"),
    ("T⁻⁴ 属于 m=1 族", "恒加速参照轨迹（三角速度）实测 −3.74，复现 T⁻⁴/ω³"),
    ("不一致如实报告", "四段 constant-jerk 加速度连续（m=2），自洽 T⁻⁶；"
     "与任务书 §5B 声称的 T⁻⁴ 不符，不改数据、不迎合"),
], size=14)

# ================================================================ 13 功-能 & 收敛
slide = add_slide()
title_bar(slide, "功-能核算与步长收敛",
          "ΔE(t)、W_ext(t)、R_W(t) 由真实轨迹时间序列记录（非首末连线）；三分解：移动 + 焦点偏移 + 深度")
add_picture_fit(slide, FIG / "work_energy_balance.png", 0.55, 1.15, 7.6, 2.9)
add_picture_fit(slide, FIG / "timestep_convergence.png", 0.55, 4.2, 7.6, 2.7)
add_textbox(slide, 8.5, 1.3, 4.6, 5.3,
            "● R_W 随步长二阶收敛；0.10/0.05/0.025 μs\n  配对复算末态能量差、残差、"
            "标签翻转率\n  均达标（preflight 第 11/12/14 项）\n\n"
            "● 粗扫描功-能残差 |R_W|/D < 5×10⁻⁴\n\n"
            "● 静态 100 μs 保持段峰峰值\n  能量误差 < 10⁻³ 阱深\n\n"
            "● validation 与 0.025 μs 子集\n  留阱标签 100% 一致", size=14,
            color=DARK)

# ================================================================ 14 敏感性
slide = add_slide()
title_bar(slide, "阶段 D：敏感性（one-factor-at-a-time，500 shots/点）",
          "基准：对角 610 μm / 1.6 ms / sine / sev1.0")
add_picture_fit(slide, FIG / "sensitivity_panels.png", 0.55, 1.15, 12.2, 4.3)
add_bullets(slide, 0.55, 5.6, 12.3, 1.6, [
    ("稳健区", "温度 3–10 μK、深度 160–920 μK、束腰 ±0.06 μm、severity 0.5–1.25 "
     "下留阱率均 1.0000（Wilson 下界 ≥0.9924）"),
    ("唯一退化通道", "AOD 反号 (−1,+1)：留阱率 0.944 [0.920, 0.961]——两柱面焦点"
     "反向偏移使势更浅；仅作模型敏感性，不冒充论文配置"),
], size=13)

# ================================================================ 15 深度补偿
slide = add_slide()
title_bar(slide, "探索性 lensing 深度补偿（非论文已实现功能）",
          "straight 510 μm / 1.6 ms / sine / sev1.0，500 shots 同初态对照；D_cmd(t) = D·clip(D/d_eff(t), 1, 上限)")
add_table(slide, 0.55, 1.25, 7.4, [
    ["模式", "留阱率", "retained 中位 Δε"],
    ["lensing_uncompensated", "1.0000", "15.83 μK"],
    ["exploratory_lensing_compensation", "1.0000", "17.65 μK"],
], col_widths=[4.2, 1.5, 1.7], size=13)
add_bullets(slide, 0.55, 2.9, 7.4, 3.6, [
    ("功率需求", "峰值所需倍率 2.36，实际上限 2.0，饱和时间比例 25.9%"),
    ("有效阱深", "未补偿最低 118.4 μK → 补偿后 236.8 μK（基准 280 μK）"),
    ("诚实声明", "达到功率上限，不得声称『保持恒定阱深』——报告明确标注 "
     "not_maintained_power_cap_reached"),
    ("结论", "该条件下补偿不改变留阱率、反而略增加热；论文式 RF 功率-位置校准"
     "留作 Level 4"),
], size=14)
add_textbox(slide, 8.2, 1.4, 4.9, 5.2,
            "实现要点\n\n● 补偿表在控制网格上预计算\n  （321 点）并线性插值，不逐\n  shot 逐步调用优化器——极值\n  诊断与轨迹推进分离\n\n● 补偿只依赖控制路径（轨迹、\n  焦点偏移），与原子状态无关\n\n● 端点倍率 = 1，运输终点与\n  未补偿条件末态定义一致\n\n● 输出 compensation.csv 与\n  compensation_schedule.csv\n  （倍率/命令深度/有效阱深时序）",
            size=13, color=GRAY)

# ================================================================ 16 软件工程
slide = add_slide()
title_bar(slide, "实现与质量保障",
          "simulation/level2_joint_transfer/src/level3_3d_transport_lensing/，共 2619 行")
add_table(slide, 0.55, 1.2, 12.2, [
    ["模块", "职责"],
    ["gaussian_3d.py / aod_lensing.py", "三维势、力、频率；G₁G₂ 透镜势、∂U/∂t 三分解、数值极值诊断、v_s 反解"],
    ["transport_profiles.py / integrators_3d.py", "三类轨迹与几何分解；向量化时变 velocity-Verlet + 功-能核算 + 静态保持"],
    ["thermal_sampling_3d.py / transport_statistics.py", "三维热采样；复用 Level 1 Wilson 区间与 Level 2 配对 bootstrap"],
    ["depth_compensation.py / heating_scaling.py", "探索性补偿调度；加热标度律检验（T⁻⁶ / T⁻⁴ 判别）"],
    ["level3_simulation.py / preflight3.py", "扫描-粗MC-验证-敏感性编排；16 项预检查（先跑 Level 0–2 全部测试）"],
    ["level3_cli.py / level3_visualization.py", "CLI 六个 stage；12 张诊断图 + 程序化验收记录（有限性/图例/坐标覆盖）"],
], col_widths=[4.2, 8.0], size=12)
add_bullets(slide, 0.55, 4.6, 12.2, 2.4, [
    ("测试", "test_level3.py 共 22 项测试函数，覆盖任务书 25 项要求：势/力/轨迹/采样/"
     "功-能/判据/种子隔离/CLI 冒烟/JSON 严格性/sources 保护"),
    ("产物", "26 个文件：CSV×8、PNG×12、JSON×4、Markdown×1、YAML×1，"
     "全部程序化复检 schema 与有限值；支持事后合并视觉审阅（--apply-visual-review）"),
], size=13)

# ================================================================ 17 结论与边界
slide = add_slide()
title_bar(slide, "结论、边界与 Level 4 建议")
add_bullets(slide, 0.55, 1.2, 6.1, 5.6, [
    ("透镜效应是直线短时长运输的\n唯一经典损失通道", "sev≥0.75 且 T≤800 μs 时"
     "直线留阱率归零；对角因单轴速度 1/√2 + 球面阱全区间存活"),
    ("论文锚点条件深度安全", "610 μm / 1.6 ms / 对角在全部 severity 下"
     "2000/2000 留阱，配对差异不显著"),
    ("轨迹族在安全区不可分辨", "sine 与 constant-jerk 差异为 0；"
     "差异只在粗扫描短时长区显现"),
    ("标度律修正", "实测 constant-jerk 加热 ∝ T⁻⁶（m=2 自洽），"
     "T⁻⁴ 属恒加速族——任务书 §5B 有误，已如实报告"),
], size=14)
add_textbox(slide, 6.95, 1.3, 5.9, 2.4,
            "边界声明\n\n● 经典留阱率 ≠ 论文 99.953(2)% 量子保真度\n"
            "● severity 为无量纲场景，非论文测量\n"
            "● 未含：真空损失、光子散射、退相干、\n  pick-up/drop-off 组合、重复往返",
            size=14, color=DARK)
add_textbox(slide, 6.95, 4.0, 5.9, 2.8,
            "Level 4 建议\n\n1. pick-up—长距离运输—drop-off 组合序列\n"
            "   （含深度 ramp 与 SLM 关断）\n"
            "2. 重复往返与损失累计/自洽存活模型\n"
            "3. 量子内态退相干（XY4、强度噪声）与成像误差\n"
            "4. 论文式 RF 功率-位置校准的 lensing 补偿",
            size=14, color=DARK)

prs.save(OUT)
print(f"saved: {OUT} ({len(prs.slides.__iter__.__self__._sldIdLst)} slides)"
      if False else f"saved: {OUT} ({len(prs.slides._sldIdLst)} slides)")
