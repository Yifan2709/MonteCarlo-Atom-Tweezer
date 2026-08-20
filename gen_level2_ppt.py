# -*- coding: utf-8 -*-
"""生成 Level 2 总结 PPT：用 level2_figures/ 的 8 张图 + 数据要点组装幻灯片。"""
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

REPO = Path(__file__).resolve().parent
FIG = REPO / "level2_figures"
DEMO = REPO / "simulation" / "level2_joint_transfer" / "outputs" / "level2_joint_transfer" / "demo"
OUT = REPO / "Level2_同步波形优化_总结报告.pptx"

BLUE = RGBColor(0x1F, 0x77, 0xB4)
DARK = RGBColor(0x22, 0x2A, 0x35)
RED = RGBColor(0xD6, 0x27, 0x28)
GREEN = RGBColor(0x2C, 0xA0, 0x2C)
GRAY = RGBColor(0x66, 0x66, 0x66)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]


def add_slide():
    return prs.slides.add_slide(BLANK)


def set_text(box, text, size=16, color=DARK, bold=False, align=PP_ALIGN.LEFT):
    """写入多行文本，返回首个 paragraph。"""
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


def add_bullets(slide, x, y, w, h, items, size=16, gap=Pt(6)):
    """带项目符号的多行文本（手动加 • 便于控制）。"""
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
    bar = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width, Inches(0.95))
    bar.fill.solid()
    bar.fill.fore_color.rgb = BLUE
    bar.line.fill.background()
    add_textbox(slide, 0.5, 0.12, 12.4, 0.6, title, size=26, color=WHITE, bold=True)
    if subtitle:
        add_textbox(slide, 0.55, 0.62, 12.3, 0.35, subtitle, size=13,
                    color=RGBColor(0xDD, 0xE8, 0xF5))


def add_picture_fit(slide, path: Path, x, y, w, h):
    """把图片等比缩放到 (x,y,w,h) 框内居中。"""
    from PIL import Image
    with Image.open(path) as image:
        iw, ih = image.size
    scale = min(w / iw, h / ih)
    nw, nh = iw * scale, ih * scale
    left = Inches(x + (w - nw) / 2)
    top = Inches(y + (h - nh) / 2)
    return slide.shapes.add_picture(str(path), left, top, Inches(nw), Inches(nh))


# ================================================================ 1 封面
slide = add_slide()
bg = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width, prs.slide_height)
bg.fill.solid()
bg.fill.fore_color.rgb = RGBColor(0x0E, 0x2A, 0x47)
bg.line.fill.background()
add_textbox(slide, 1.0, 2.1, 11.3, 1.2,
            "Monte Carlo for AOD→SLM 原子转移", size=30,
            color=RGBColor(0x9E, 0xC5, 0xE8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 3.0, 11.3, 1.4,
            "Level 2：同步移动-降深波形优化\n与鲁棒蒙特卡洛验证", size=40,
            color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 5.4, 11.3, 0.5,
            "600 μs 顺序基线 vs 400 μs 压缩基线 vs 400 μs 优化同步波形",
            size=16, color=RGBColor(0xBB, 0xD3, 0xE8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 6.6, 11.3, 0.4,
            "88 项测试全部通过 · 257 个优化候选 · 2000 shots 独立验证",
            size=13, color=RGBColor(0x88, 0xA8, 0xC8), align=PP_ALIGN.CENTER)

# ================================================================ 2 要解决什么问题
slide = add_slide()
title_bar(slide, "要回答的问题", "Level 2 的定位")
add_bullets(slide, 0.7, 1.4, 11.9, 3.4, [
    ("物理场景", "280 μK AOD 光镊抓住一颗铯原子，旁边是静止的 140 μK SLM 光镊；"
     "AOD 一边挪向 SLM 一边关灯，把原子留给 SLM。"),
    ("Level 1 基线", "顺序式两步走：前 52% 时间只移动（600 μs），后 48% 只降深。"),
    ("Level 2 问题", "允许位置与深度同时变化（重叠），能否把转移从 600 μs 压缩到 400 μs，"
     "同时保持或提高 SLM 俘获率、降低末态激发？"),
], size=18)
note = slide.shapes.add_shape(1, Inches(0.7), Inches(4.7), Inches(11.9), Inches(1.9))
note.fill.solid()
note.fill.fore_color.rgb = RGBColor(0xFD, 0xF3, 0xE7)
note.line.color.rgb = RGBColor(0xE8, 0xA8, 0x60)
box = add_textbox(slide, 1.0, 4.95, 11.3, 1.5,
                  "方法学底线（本项目的核心价值）\n"
                  "· 三批互不重叠的随机样本池，防止蒙特卡洛过拟合\n"
                  "· 优化后未优于基线 ≠ 失败：如实报告『未观察到显著提升』，不得用验证集回调参数",
                  size=14)
add_textbox(slide, 0.7, 6.8, 11.9, 0.4,
            "模型边界：一维经典动力学；不含二维/轴向、光子散射、技术噪声、量子内态与多原子相互作用。",
            size=12, color=GRAY)

# ================================================================ 3 物理模型
slide = add_slide()
title_bar(slide, "物理模型（沿用 Level 0/1 已验证实现）", "一维焦平面径向高斯势")
add_bullets(slide, 0.7, 1.35, 6.6, 4.4, [
    ("SLM 静止阱", "U_SLM(x) = −D_SLM·exp(−2x²/w²)，D_SLM = 140 μK，w = 1.17 μm"),
    ("AOD 时变阱", "U_AOD(x,t) = −D(t)·exp(−2[x−c(t)]²/w²)，初始 D₀ = 280 μK，c₀ = 2.4 μm"),
    ("动力学", "时变 velocity-Verlet 推进；转移后追加 100 μs 纯 SLM 保持段"),
    ("俘获判据", "末态 SLM 单阱能量 E_f < 0 ⇒ captured（严格符号判据）"),
], size=15)
panel = slide.shapes.add_shape(1, Inches(7.5), Inches(1.35), Inches(5.2), Inches(4.4))
panel.fill.solid()
panel.fill.fore_color.rgb = RGBColor(0xEA, 0xF3, 0xFB)
panel.line.color.rgb = BLUE
add_textbox(slide, 7.8, 1.55, 4.7, 4.0,
            "显式时变势中的功-能检查\n\n"
            "转移期间机械能不守恒——外部控制在做功：\n\n"
            "dE/dt = ∂U_AOD/∂t = 深度功 + 移动功\n\n"
            "数值上验证：\n"
            "R_W = E(t) − E(0) − W_ext(t) ≈ 0\n\n"
            "残差实测 ≪ 1e-3（相对阱深），\n"
            "保持段能量仅有界振荡、无漂移",
            size=14)
add_textbox(slide, 0.7, 6.15, 12.0, 0.8,
            "初态：5 μK 简谐近似热分布抽样 + 完整高斯阱束缚拒绝（E_AOD<0）；"
            "记录验收率、均值、标准差与种子。w_ADM=1.17 μm 为与 SLM 后孔径匹配的假设值。",
            size=13, color=GRAY)

# ================================================================ 4 方法总览
slide = add_slide()
title_bar(slide, "方法总览：三池隔离 + 公共随机数", "防止用验证集『作弊』的完整闭环")
add_picture_fit(slide, FIG / "01_流程总览.png", 0.25, 1.0, 12.8, 6.3)

# ================================================================ 5 波形对比
slide = add_slide()
title_bar(slide, "三类主波形", "优化器调整的 5 个参数：移动窗口×2、降深窗口×2、降深幂次")
add_picture_fit(slide, FIG / "02_波形对比.png", 0.25, 1.0, 12.8, 5.4)
add_textbox(slide, 0.7, 6.55, 12.0, 0.7,
            "『同步』的含义：优化波形让移动几乎贯穿全程（s∈[0.019, 1]），降深在 s≥0.607 才开始——"
            "移动与降深重叠约 39% 的时间；约束保证单调、无过冲、端点速度/加速度为零。",
            size=13, color=GRAY)

# ================================================================ 6 势能演化
slide = add_slide()
title_bar(slide, "势场演化：原子如何被『倒』进 SLM 阱")
add_picture_fit(slide, FIG / "03_势能剖面演化.png", 0.25, 1.0, 12.8, 6.3)

# ================================================================ 7 优化过程
slide = add_slide()
title_bar(slide, "优化过程", "拉丁超立方探索 → 精修选择 → PCHIP 样条 → 冻结")
add_picture_fit(slide, FIG / "04_优化过程.png", 0.25, 1.0, 8.9, 6.2)
panel = slide.shapes.add_shape(1, Inches(9.3), Inches(1.2), Inches(3.7), Inches(5.4))
panel.fill.solid()
panel.fill.fore_color.rgb = RGBColor(0xEE, 0xF8, 0xEE)
panel.line.color.rgb = GREEN
add_textbox(slide, 9.5, 1.4, 3.3, 5.0,
            "分阶段预算\n\n"
            "阶段1：128 shots\nLHS 探索 257 候选\n（有效229/无效28/失败0）\n\n"
            "阶段2：512 shots\n前 8 名入围精修\n+ 8 控制点单调样条\n\n"
            "冻结：唯一波形写入\nbest_waveform.yaml\n\n"
            "目标 = 俘获率(主)\n+ 裕量分位数\n+ 末态激发 + 光滑度\n（各组成项分开保存）",
            size=13)

# ================================================================ 8 独立验证
slide = add_slide()
title_bar(slide, "独立验证（2000 shots × 3 波形，同一初态数组）")
add_picture_fit(slide, FIG / "05_独立验证结果.png", 0.25, 1.0, 12.8, 5.3)
add_textbox(slide, 0.7, 6.45, 12.0, 0.8,
            "结果：三波形俘获率均为 2000/2000 = 1.0000（Wilson 95% CI [0.9981, 1.0000]），"
            "配对差全部为 0，末态激发均值≈0.025。功-能残差与保持段误差远低于 1e-3 验收阈值。",
            size=13, color=GRAY)

# ================================================================ 9 初末态
slide = add_slide()
title_bar(slide, "蒙特卡洛初态与末态分布")
add_picture_fit(slide, FIG / "06_初态与末态分布.png", 0.25, 1.1, 12.8, 5.6)
add_textbox(slide, 0.7, 6.8, 12.0, 0.5,
            "左：5 μK 热初态相空间（集中于 AOD 阱底附近）；右：末态位置集中于 SLM 阱（0 μm）内。",
            size=13, color=GRAY)

# ================================================================ 10 鲁棒性
slide = add_slide()
title_bar(slide, "鲁棒性敏感性（冻结参数后不再调优）")
add_picture_fit(slide, FIG / "07_鲁棒性.png", 0.25, 1.0, 12.8, 5.4)
add_textbox(slide, 0.7, 6.55, 12.0, 0.7,
            "温度 3/5/7/10 μK、末端对准偏差 ±0.10 μm（各 500 shots）下俘获率保持 100%；"
            "步长 0.10→0.025 μs 连续量收敛、俘获标签翻转≈0。",
            size=13, color=GRAY)

# ================================================================ 11 Preflight
slide = add_slide()
title_bar(slide, "Preflight：不过关不得启动优化")
add_picture_fit(slide, FIG / "08_预检查清单.png", 0.25, 1.0, 8.9, 6.2)
panel = slide.shapes.add_shape(1, Inches(9.3), Inches(1.2), Inches(3.7), Inches(5.4))
panel.fill.solid()
panel.fill.fore_color.rgb = RGBColor(0xFD, 0xEE, 0xEE)
panel.line.color.rgb = RED
add_textbox(slide, 9.5, 1.45, 3.3, 5.0,
            "为什么重要？\n\n"
            "· 时变势中 E(t) 变化\n  ≠ 积分器漂移，\n  必须用功-能关系区分\n\n"
            "· 步长减半收敛 +\n  MC 配对标签翻转率\n  实测 = 0\n\n"
            "· 任一项失败即终止，\n  不会生成『看似完整』\n  的最优结果\n\n"
            "本次：9/9 关键检查通过\n（284 s，含 88 项旧测试）",
            size=13)

# ================================================================ 12 测试与代码质量
slide = add_slide()
title_bar(slide, "测试与代码质量", "88 项测试全部通过")
add_bullets(slide, 0.7, 1.3, 7.3, 5.4, [
    ("波形", "端点严格正确、单调有界、非法参数/奇异端点拒绝、解析导数 vs 有限差分"),
    ("复现", "sequential_600us 与 Level 1 轨迹逐点一致；时变积分器静态极限与 Level 0 一致"),
    ("防泄漏", "三池种子互不相同；spy 测试验证优化全程只抽取 optimization/selection 池"),
    ("公平性", "CRN：三波形初态逐行完全相同；near-threshold 样本保留并按严格符号计数"),
    ("统计", "Wilson 区间、配对 bootstrap 在构造数据上验证正确"),
    ("工程", "CLI 各 stage 冒烟测试、检查点断点续跑、CSV/JSON schema 与图片程序化检查"),
], size=14)
panel = slide.shapes.add_shape(1, Inches(8.3), Inches(1.3), Inches(4.5), Inches(5.2))
panel.fill.solid()
panel.fill.fore_color.rgb = RGBColor(0xF2, 0xF2, 0xF2)
add_textbox(slide, 8.55, 1.5, 4.0, 4.8,
            "复现命令\n\n"
            "$ pip install -e \".[dev]\"\n"
            "$ pytest -q\n"
            "→ 88 passed\n\n"
            "$ level2-joint-transfer \\\n"
            "    --config configs/…yaml \\\n"
            "    --stage all\n\n"
            "阶段：preflight / optimize /\nvalidate / robustness / all",
            size=12)

# ================================================================ 12.5 功-能平衡
slide = add_slide()
title_bar(slide, "功-能核算：转移期间外部做功的显式验证", "demo 输出 · 多条代表轨迹")
add_picture_fit(slide, DEMO / "work_energy_balance.png", 0.25, 1.0, 12.8, 5.3)
add_textbox(slide, 0.7, 6.45, 12.0, 0.8,
            "黑线 ΔE(t) 与绿线外部总功 W_ext(t) 逐点重合，残差 R_W（红线）远小于 1e-3·D_SLM；"
            "深度功在降深段把能量取走，移动功在移动段注入能量，两者符号与量级符合物理直觉。",
            size=13, color=GRAY)

# ================================================================ 12.6 步长收敛
slide = add_slide()
title_bar(slide, "步长收敛与俘获标签稳定性", "demo 输出 · dt = 0.10 / 0.05 / 0.025 μs")
add_picture_fit(slide, DEMO / "timestep_convergence.png", 0.25, 1.0, 12.8, 5.3)
add_textbox(slide, 0.7, 6.45, 12.0, 0.8,
            "左：末态能量最大偏差随 dt 减小呈二阶下降（对数坐标）；右：三种波形俘获标签翻转率 ≈ 0，"
            "说明 dt=0.05 μs 的验证结果不是积分误差的产物。",
            size=13, color=GRAY)

# ================================================================ 13 结论
slide = add_slide()
title_bar(slide, "结论", "诚实报告是本 Level 的核心要求")
add_bullets(slide, 0.7, 1.35, 11.9, 3.6, [
    ("实现完成", "三类波形、功-能核算、分阶段优化、独立验证、鲁棒性检查全部实现并实际运行。"),
    ("物理结论", "5 μK 工况下问题处于饱和区：三波形俘获率均 100%，配对差为 0 —— "
     "未观察到统计上明确的提升，如实报告，不得据此调参。"),
    ("数值可信", "preflight 10/10 通过；功-能残差、保持段误差、步长收敛均远优于验收阈值。"),
    ("鉴别力限制", "如需区分波形优劣，应升温或进一步压缩时长，把系统推离 100% 饱和区。"),
], size=17)
bar = slide.shapes.add_shape(1, Inches(0.7), Inches(5.4), Inches(11.9), Inches(1.6))
bar.fill.solid()
bar.fill.fore_color.rgb = RGBColor(0xE8, 0xF0, 0xFE)
bar.line.color.rgb = BLUE
add_textbox(slide, 1.0, 5.6, 11.3, 1.2,
            "留给后续 Level 的物理过程：二维/轴向运动、AOD 柱面透镜效应、光子散射、技术噪声、"
            "量子内态与相干性、加热/损失随机过程、多原子相互作用、反向 pick-up 与 610 μm 长距离运输。",
            size=14)

prs.save(OUT)
print(f"已生成 {OUT}（{len(prs.slides.slides if hasattr(prs.slides, 'slides') else prs.slides._sldIdLst)} 页）"
      if False else f"已生成 {OUT}")
