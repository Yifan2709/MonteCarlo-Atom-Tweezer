# -*- coding: utf-8 -*-
"""生成 Level 2B 总结 PPT：用 level2b demo 输出的 8 张图 + 实际数据组装幻灯片。

图源：simulation/level2b_speed_modify/outputs/level2b_speed_limit/demo/
数据：metrics.json / speed_scan_results.csv / comparison_results.csv（demo 运行）
"""
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

REPO = Path(__file__).resolve().parent
DEMO = REPO / "simulation" / "level2b_speed_modify" / "outputs" / "level2b_speed_limit" / "demo"
OUT = REPO / "Level2B_移动速度边界与变速剖面验证_总结报告.pptx"

BLUE = RGBColor(0x1F, 0x77, 0xB4)
DARK = RGBColor(0x22, 0x2A, 0x35)
RED = RGBColor(0xD6, 0x27, 0x28)
GREEN = RGBColor(0x2C, 0xA0, 0x2C)
ORANGE = RGBColor(0xE8, 0xA8, 0x60)
GRAY = RGBColor(0x66, 0x66, 0x66)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

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


def add_bullets(slide, x, y, w, h, items, size=16, gap=Pt(6)):
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
    from PIL import Image
    with Image.open(path) as image:
        iw, ih = image.size
    scale = min(w / iw, h / ih)
    nw, nh = iw * scale, ih * scale
    left = Inches(x + (w - nw) / 2)
    top = Inches(y + (h - nh) / 2)
    return slide.shapes.add_picture(str(path), left, top, Inches(nw), Inches(nh))


def panel(slide, x, y, w, h, fill, line):
    shape = slide.shapes.add_shape(1, Inches(x), Inches(y), Inches(w), Inches(h))
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill
    shape.line.color.rgb = line
    return shape


# ================================================================ 1 封面
slide = add_slide()
bg = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width, prs.slide_height)
bg.fill.solid()
bg.fill.fore_color.rgb = RGBColor(0x0E, 0x2A, 0x47)
bg.line.fill.background()
add_textbox(slide, 1.0, 1.7, 11.3, 0.8,
            "Monte Carlo for AOD→SLM 原子转移 · 独立扩展工程 level2b_speed_modify",
            size=20, color=RGBColor(0x9E, 0xC5, 0xE8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 2.6, 11.3, 1.6,
            "Level 2B：AOD 移动速度边界\n与变速剖面验证", size=40,
            color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 4.4, 11.3, 0.6,
            "纯移动跟随上限 · 端到端 v95 临界速度 · 恒速 vs 变速双口径配对比较",
            size=17, color=RGBColor(0xBB, 0xD3, 0xE8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 5.2, 11.3, 0.5,
            "9 条曲线全部覆盖失效边界 · η_95 ≈ 0.29 坍缩 · 独立复核 30/30 CI 重叠",
            size=15, color=RGBColor(0x88, 0xA8, 0xC8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 5.8, 11.3, 0.5,
            "实验 4 噪声：默认幅度不移动边界（v95 完全不变）· ×100 参量通道使 v95 腰斩至 20.8 mm/s",
            size=15, color=RGBColor(0x88, 0xA8, 0xC8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 6.6, 11.3, 0.4,
            "preflight 8/8 · 测试 105 项全过 · 原 Level 2 工程零改动（校验和一致）",
            size=13, color=RGBColor(0x88, 0xA8, 0xC8), align=PP_ALIGN.CENTER)

# ================================================================ 2 问题定位
slide = add_slide()
title_bar(slide, "为什么做 Level 2B：Level 2 的饱和教训", "把『多快算太快』变成可测量的问题")
add_bullets(slide, 0.7, 1.35, 11.9, 2.6, [
    ("Level 2 结论回顾", "5 μK、400 μs 下三波形俘获率均 2000/2000 = 100%，配对差全为 0 —— "
     "问题处于饱和区，无法区分波形优劣。"),
    ("未测量过的量", "当前 AOD 移动速度：优化波形峰值 11.5 mm/s、压缩基线 17.3 mm/s；"
     "距离失效有多远、失效边界长什么样，完全未知。"),
    ("Level 2B 三个问题", "① 纯移动时原子多快会跟丢？② 完整 drop-off 的临界速度 v95/v50 是多少、"
     "随温度怎么变？③ 边界附近恒速巡航与变速剖面谁更优？"),
], size=17)
note = panel(slide, Inches(0.7), Inches(4.3), Inches(11.9), Inches(1.7),
             RGBColor(0xFD, 0xF3, 0xE7), ORANGE)
add_textbox(slide, 1.0, 4.5, 11.3, 1.3,
            "两条诚实底线（违反即失败）\n"
            "· 饱和是合法结论：若扫描区间未覆盖失效边界，必须写『v95 > 最大扫描速度』，禁止宣称已测得最大速度\n"
            "· 临界速度是估计量：必须带 bootstrap/logistic 双方法置信区间；边界点用从未用过的新种子独立复核",
            size=14)
add_textbox(slide, 0.7, 6.35, 11.9, 0.6,
            "工程定位：新文件夹 simulation/level2b_speed_modify/（tag: level2 + speed_modify），"
            "level0/1/2 作为逐字节副本基线，原工程零改动。", size=13, color=GRAY)

# ================================================================ 3 物理模型与派生尺度
slide = add_slide()
title_bar(slide, "物理模型（沿用 Level 2，逐参数一致）与派生尺度", "一维经典 + 显式速度变量")
add_bullets(slide, 0.7, 1.35, 6.6, 4.2, [
    ("模型不变", "Cs 一维焦平面双高斯阱：AOD 280 μK 时变（中心 c(t)、深度 D(t)），SLM 140 μK 静止；"
     "w=1.17 μm、d=2.4 μm、5 μK 热初态（Level 1 采样器）、velocity-Verlet、E_f<0 严格判据。"),
    ("绝热速度尺度", "v_ad = w_AOD·ω_AOD：阱频与束宽的乘积，量纲即速度。"),
    ("绝热参数", "η = v_peak / v_ad：所有速度结果用它无量纲化。"),
    ("每步位移数", "δ = dt·v_peak / w：步长充分性的显式指标，全程 ≤ 0.01（规则保证）。"),
    ("显式排除", "AOD 声学器件带宽/RF 扫描速率限制 —— 测的是动力学上限，不是器件规格。"),
], size=15)
panel(slide, Inches(7.5), Inches(1.35), Inches(5.2), Inches(4.2),
      RGBColor(0xEA, 0xF3, 0xFB), BLUE)
add_textbox(slide, 7.8, 1.55, 4.7, 3.8,
            "实测派生量（demo 运行）\n\n"
            "ν_AOD = 36.0 kHz（周期 27.8 μs）\n"
            "ν_SLM = 25.5 kHz（周期 39.2 μs）\n\n"
            "v_ad = 264.7 mm/s\n\n"
            "关键锚点：\n"
            "当前 Level 2 优化波形 v_peak = 11.5 mm/s\n"
            "→ η = 0.043，离失效边界约 7 倍余量\n"
            "（这正是 Level 2 100% 饱和的原因）",
            size=14)
add_textbox(slide, 0.7, 5.8, 12.0, 0.8,
            "速度三口径强制分开报告：v_avg = d/T_move、v_peak = max|ċ|、v_rms —— "
            "『移动速度』在任何场合都必须说明口径，这是恒速 vs 变速公平比较的前提。",
            size=13, color=GRAY)

# ================================================================ 4 方法总览
slide = add_slide()
title_bar(slide, "方法总览：三组递进实验", "恒速巡航剖面让『最大速度』成为直接参数")
cols = [
    ("实验 1 · 纯移动跟随上限", RGBColor(0xE8, 0xF0, 0xFE), BLUE,
     "深度恒定 280 μK\n中心梯形巡航移 2.4 μm\n\n扫描 v_cruise\n6 → 500 mm/s（11 个近似对数点）\n\n"
     "判定：|x−c| > 2w 或 ε/D₀ > 0.5\n输出：跟丢率 + 激发标度律"),
    ("实验 2 · 端到端速度上限", RGBColor(0xEE, 0xF8, 0xEE), GREEN,
     "三结构 × 三温度 × 12 时长\n(600 μs → 20 μs，等价 v_avg 4→120 mm/s)\n\n"
     "sequential / 冻结形状缩放 / 梯形+冻结降深\nCRN 同初态、曲线级早停(<90%)\n\n"
     "输出：v95/v50 + 双方法 CI\n失效三分类（跟丢/降深段/交接）"),
    ("实验 3 · 恒速 vs 变速", RGBColor(0xFD, 0xEE, 0xEE), RED,
     "边界附近 5 个速度点\n三剖面：梯形 f=0.30 / smootherstep / 梯形 f=0.15\n\n"
     "双口径配对比较：\n相同 v_avg（梯形峰值更低）\n相同 v_peak（梯形均值更低）\n\n"
     "配对 bootstrap + 独立种子复核"),
]
for i, (head, fill, line, body) in enumerate(cols):
    x = 0.6 + i * 4.15
    panel(slide, Inches(x), Inches(1.3), Inches(3.9), Inches(4.5), fill, line)
    add_textbox(slide, x + 0.2, 1.45, 3.5, 0.5, head, size=15, bold=True)
    add_textbox(slide, x + 0.2, 2.0, 3.5, 3.7, body, size=12)
add_textbox(slide, 0.7, 6.05, 11.9, 0.9,
            "防作弊设计：scan_seed=22001 与 boundary_seed=22002 互不相同（有测试守护）；"
            "逐点落盘断点续跑；每步位移规则 dt ≤ min(50 ns, T_move/2000, 0.01·w/v_peak)；"
            "势/积分器/采样器/统计/噪声模型全部 import 自 level0/1/2，不复制第二套。"
            "另有实验 4（噪声 × 速度）在冻结基线之上做敏感性评估（详见后文）。",
            size=13, color=GRAY)

# ================================================================ 5 剖面族
slide = add_slide()
title_bar(slide, "核心剖面：梯形恒速巡航（trapezoid_cruise）",
          "smootherstep 加/减速 + 恒速巡航 —— v_cruise 直接就是待测最大速度")
add_picture_fit(slide, DEMO / "trapezoid_vs_smoothstep.png", 0.25, 1.0, 12.8, 5.3)
add_textbox(slide, 0.7, 6.45, 12.0, 0.8,
            "左/中：两种剖面的中心 c(t) 与速度 |ċ(t)| —— 梯形有平坦巡航段，smootherstep 峰值居中。"
            "右：实验 3 的俘获率结果。数学要点：v_c = d/[(1−f)·T_move]，峰值速度严格等于 v_c（校准 <1e-6，有单元测试）；"
            "端点速度/加速度为零、位移严格闭合 2.4 μm。", size=13, color=GRAY)

# ================================================================ 6 预检查
slide = add_slide()
title_bar(slide, "Preflight：不过关不得启动扫描", "本级 6 项 + Level 2 十项复跑 + 全部既有测试")
add_bullets(slide, 0.7, 1.3, 7.6, 5.2, [
    ("既有测试", "本工程全部测试（level0/1/2 基线 + level2b 新增）一次通过。"),
    ("Level 2 复跑", "原 preflight 十项检查在副本工程全部通过 —— 证明基线未因迁移而退化。"),
    ("梯形校准", "v_peak == v_cruise（相对误差 <1e-6）、端点零速/零加速度、位移严格 d、单调无过冲（f=0.15/0.30/0.45 × 多时长）。"),
    ("导数一致性", "新剖面解析一/二阶导数 vs 中心差分一致。"),
    ("缩放不变性", "T→T/2 时 v_peak×2、a_max×4、归一化形状逐点一致。"),
    ("慢极限回归", "梯形 v_avg≈6.1 mm/s vs Level 2 冻结 400 μs 波形：同初态配对差 CI 含 0、激发差 <0.02 —— 新剖面在慢极限退化为已验证结果。"),
    ("最快点收敛", "v=500 mm/s 处 dt 与 dt/2：标签翻转率 = 0。"),
    ("静态极限", "静止极限下引擎积分路径与 Level 0 积分器逐点一致。"),
], size=14)
panel(slide, Inches(8.6), Inches(1.3), Inches(4.2), Inches(5.2),
      RGBColor(0xFD, 0xEE, 0xEE), RED)
add_textbox(slide, 8.85, 1.5, 3.7, 4.9,
            "结果：8/8 全部通过\n（456 s，含 Level 2\n十项复跑与全部测试）\n\n"
            "实现要点：\n\n· 预检查失败 → CLI 拒绝\n  进入 scan/compare 阶段\n\n"
            "· 冒烟旗标仅测试用：\n  --skip-level2-preflight\n  --skip-existing-tests\n  正式运行不得使用\n\n"
            "· 慢极限回归是『新代码\n  不改变已验证物理』的\n  直接证据",
            size=13)

# ================================================================ 7 实验 1 跟随
slide = add_slide()
title_bar(slide, "实验 1：纯移动跟随上限", "深度恒定 280 μK · 5 μK 初态 · 梯形巡航 6→500 mm/s")
add_picture_fit(slide, DEMO / "follow_loss_diagnostics.png", 0.25, 1.0, 12.8, 5.3)
add_textbox(slide, 0.7, 6.45, 12.0, 0.8,
            "左：失效模式构成随速度演变（纯移动下只有跟丢）；右：每步中心位移 δ 始终 ≤ 0.0015·w，"
            "远低于 0.01 规则上限 —— 快速点的结果不是积分误差的产物。", size=13, color=GRAY)

# ================================================================ 8 实验 1 结果+标度律
slide = add_slide()
title_bar(slide, "实验 1 结果：跟随边界与激发标度律", "v_peak ≤ 63 mm/s（η≤0.24）零跟丢；≥158 mm/s 全跟丢")
add_picture_fit(slide, DEMO / "excitation_vs_speed_loglog.png", 0.25, 1.0, 8.6, 5.3)
panel(slide, Inches(9.2), Inches(1.2), Inches(3.9), Inches(5.3),
      RGBColor(0xEE, 0xF8, 0xEE), GREEN)
add_textbox(slide, 9.45, 1.4, 3.4, 5.0,
            "关键数字\n\n"
            "v_peak=6–63 mm/s：\n跟丢率 = 0\nε_med 0.013 → 0.038\n\n"
            "v_peak=100 mm/s：\n跟丢率 75%\nε_med = 0.59\n\n"
            "v_peak ≥158 mm/s：\n跟丢率 100%\n\n"
            "标度律拟合：\nε_med ∝ v^α\nα = 1.33\n95% CI [1.08, 1.82]\nR² = 0.83\n\n"
            "换算 η：边界位于\nη ≈ 0.24 – 0.38",
            size=13)
add_textbox(slide, 0.7, 6.45, 8.4, 0.8,
            "左：激发中位数与 90 分位随 v_peak 双对数图 + 幂律拟合线；右：跟丢率陡升曲线。"
            "注：α 拟合覆盖全部点（含跟丢区），CI 较宽，如实报告。", size=13, color=GRAY)

# ================================================================ 9 实验 2 v95
slide = add_slide()
title_bar(slide, "实验 2：端到端俘获率 vs 速度（9 条曲线）", "三结构 × 三温度 · 每点 300 shots · Wilson 95% CI")
add_picture_fit(slide, DEMO / "capture_vs_speed.png", 0.25, 1.0, 12.8, 5.3)
add_textbox(slide, 0.7, 6.45, 12.0, 0.8,
            "所有曲线均真正穿越 95% 水平（boundary_covered = true）：v95(v_avg) = 41–56 mm/s，随温度升高而下降。"
            "俘获率从 100% 崩塌到 0 发生在约 2 倍速度区间内 —— 边界相当陡峭。", size=13, color=GRAY)

# ================================================================ 10 η 坍缩
slide = add_slide()
title_bar(slide, "关键发现：失效由峰值速度主导，三结构坍缩到同一 η_95", "排序与峰值/平均速度比完全一致")
add_picture_fit(slide, DEMO / "adiabatic_collapse.png", 0.25, 1.0, 8.6, 5.3)
panel(slide, Inches(9.2), Inches(1.2), Inches(3.9), Inches(5.3),
      RGBColor(0xEA, 0xF3, 0xFB), BLUE)
add_textbox(slide, 9.45, 1.4, 3.4, 5.0,
            "20 μK 的 v95 排序\n\n"
            "frozen_shape 41.2 mm/s\nsequential  53.2 mm/s\ntrapezoid   54.5 mm/s\n\n"
            "结构峰值/平均比：\nsmootherstep 1.875\nsequential   1.50\ntrapezoid    1.43\n\n"
            "换算成 v_peak 后：\n三结构 v95 峰值速度\n= 77 / 80 / 78 mm/s\n\n"
            "→ η_95 ≈ 0.29 – 0.30\n（三结构一致坍缩）\n\n"
            "失效机制 = 峰值速度\n驱动的绝热跟随崩塌\n（与实验 1 边界衔接）",
            size=13)
add_textbox(slide, 0.7, 6.45, 8.4, 0.8,
            "左：俘获率对 η = v_peak/v_ad 的坍缩情况（温度分标记、结构分色）。"
            "含义：控制波形只要管住峰值速度，平均速度可以放心压高 —— 这直接指导实验 3 的解释。", size=13, color=GRAY)

# ================================================================ 11 临界速度 vs 温度
slide = add_slide()
title_bar(slide, "临界速度 vs 温度（双方法估计）", "bootstrap 插值 + logistic 拟合，均带 95% CI")
add_picture_fit(slide, DEMO / "critical_speed_vs_temperature.png", 0.25, 1.0, 8.6, 5.3)
panel(slide, Inches(9.2), Inches(1.2), Inches(3.9), Inches(5.3),
      RGBColor(0xF2, 0xF2, 0xF2), GRAY)
add_textbox(slide, 9.45, 1.4, 3.4, 5.0,
            "读图要点\n\n"
            "· v95 随温度下降：\n54.9 → 47.9 → 41.2\n（frozen 形状）\n\n"
            "· sequential/trapezoid\n温度敏感性较弱\n\n"
            "· 部分曲线 v50 标记\nabove_scanned：早停\n（<90% 即停）截断了\n曲线，未及 50% 水平，\n如实记录不外推\n\n"
            "· 两种方法点估计\n一致（CSV 各自全列）",
            size=13)
add_textbox(slide, 0.7, 6.45, 8.4, 0.8,
            "工程含义：若要在更保守温区（5 μK）安全运行，按 v95≈55 mm/s 留 2 倍余量取 v_avg ≤ 27 mm/s（对应总时长 ≈ 89 μs），"
            "仍比当前 400 μs 快 4.5 倍。", size=13, color=GRAY)

# ================================================================ 12 实验 3
slide = add_slide()
title_bar(slide, "实验 3：恒速 vs 变速，双口径配对比较（20 μK）", "两种公平口径分别报告，缺一不可")
add_picture_fit(slide, DEMO / "trapezoid_vs_smoothstep.png", 0.25, 1.0, 8.6, 5.3)
panel(slide, Inches(9.2), Inches(1.2), Inches(3.9), Inches(5.3),
      RGBColor(0xFD, 0xF3, 0xE7), ORANGE)
add_textbox(slide, 9.45, 1.4, 3.4, 5.0,
            "口径 1 · 相同 v_avg\n（@61.8 mm/s，均显著）\n"
            "f=0.15 梯形 0.948\nf=0.30 梯形 0.738\nsmootherstep 0.302\n"
            "→ 低峰值者胜\n\n"
            "口径 2 · 相同 v_peak\n（@82.3 mm/s）\n"
            "smootherstep 0.998\nf=0.30 梯形 0.911\nf=0.15 梯形 0.730\n"
            "→ 低均值（更长时\n  长）者胜\n\n"
            "独立复核（新种子）：\n30/30 组合 CI 重叠\n\n"
            "对称结论与 η 主导\n失效的图像自洽",
            size=13)
add_textbox(slide, 0.7, 6.45, 8.4, 0.8,
            "不存在普适『最优剖面』：比什么口径就赢什么 —— 相同预算比平均时梯形占优，"
            "相同峰值限制时平滑变速占优。任何单口径宣称都有偏。", size=13, color=GRAY)

# ================================================================ 13 边界相空间
slide = add_slide()
title_bar(slide, "边界处发生了什么：末态相空间", "v_ref = 41.2 mm/s 梯形剖面 · 300 shots · 20 μK")
add_picture_fit(slide, DEMO / "boundary_phase_space.png", 0.9, 1.0, 8.0, 5.4)
panel(slide, Inches(9.2), Inches(1.2), Inches(3.9), Inches(5.3),
      RGBColor(0xEE, 0xF8, 0xEE), GREEN)
add_textbox(slide, 9.45, 1.4, 3.4, 5.0,
            "读图\n\n"
            "· 绿色 = captured\n  （E_f < 0）\n· 红色 = 未俘获\n· 黑色虚线 =\n  E_SLM = 0 分界线\n\n"
            "物理图像：\n\n"
            "未俘获原子集中在\n分界线附近 —— 它们\n跟到了终点，但携带\n的振荡能量刚好超过\nSLM 阱深（140 μK），\n在交接瞬间溢出。\n\n"
            "与失效三分类统计\n互相印证：边界区以\n『交接失败』为主，\n更快时才转为跟丢。",
            size=13)

# ================================================================ 13.5 实验 4 模型
slide = add_slide()
title_bar(slide, "实验 4：平均强度噪声加热 × 移动速度（模型与接线）",
          "复用 level2 noise_heating 三通道确定性模型 · 新增 CLI stage: noise")
add_bullets(slide, 0.7, 1.3, 7.0, 5.2, [
    ("模型性质", "系综平均加热功率的确定性注入：无随机性、能量随时间单调累积；"
     "每步速度踢 v ← sign(v)·√(v²+2ΔE/m)，逐通道累计注入能量。"),
    ("三通道公式", "光子反冲 P=2·E_r·Γ_sc（Cs D2 两能级）；参量 dε/dt=Γ_par·ε_osc，"
     "Γ_par=(π²/4)ν₀²·S_RIN(2ν₀)；指向 P=m·ω₀⁴·S_x(ν₀)/8。PSD 均为单边每 Hz。"),
    ("默认幅度", "反冲 0.355 μK/s（Γ_sc=2.77 s⁻¹, E_r=64 nK）；参量 Γ_par=16 s⁻¹"
     "（−80 dBc/Hz）；指向 1.31 μK/s（1 pm²/Hz）；ν₀=25.46 kHz（SLM 解析阱频）。"),
    ("接线方式", "speed_engine 两段积分（转移+保持）均切换到 velocity_verlet_time_dependent_"
     "heating；ε_osc 取双阱较深阱底（与 level2 逐字一致）；逐 shot 记录 4 列注入能量；"
     "功-能恒等式扩展为 ΔE = W_ext + W_noise + R_W。"),
    ("防作弊约定", "scan 阶段保持无噪声（冻结基线不被触碰）；噪声幅度全部 "
     "assumed_sensitivity_only；噪声结果只做敏感性评估、不回调波形参数。"),
], size=13.5)
panel(slide, Inches(7.9), Inches(1.3), Inches(4.9), Inches(5.2),
      RGBColor(0xEA, 0xF3, 0xFB), BLUE)
add_textbox(slide, 8.15, 1.5, 4.4, 4.9,
            "四个子实验\n\n"
            "4a 边界移动\nfrozen_shape @20 μK，12 时长\n× 档 {关, ×1, ×100}\n→ v95 是否移动\n\n"
            "4b 开/关配对\n4 峰值速度 × 2 温度\n同初态 CRN 配对 bootstrap\n+ 裕量退化/注入能量\n\n"
            "4c 幅度敏感性\n边界速度 55 mm/s\n×{1,10,100,1000}\n→ 退化起点\n\n"
            "4d 通道分解\n×100 下仅反冲/仅参量/仅指向\n→ 谁主导退化\n\n"
            "预算：300 shots/点\n61 行全部成功落盘",
            size=13)

# ================================================================ 13.6 实验 4 结果
slide = add_slide()
title_bar(slide, "实验 4 结果：默认幅度不移动边界；×100 参量通道使 v95 腰斩",
          "noise_v95.json / noise_speed_results.csv（300 shots/点）")
add_picture_fit(slide, DEMO / "noise_speed_validation.png", 0.25, 1.0, 8.4, 5.4)
panel(slide, Inches(9.0), Inches(1.2), Inches(4.1), Inches(5.4),
      RGBColor(0xFD, 0xF3, 0xE7), ORANGE)
add_textbox(slide, 9.25, 1.35, 3.6, 5.2,
            "4a v95（bootstrap）\n"
            "关   40.8 [40.2,41.7] mm/s\n"
            "×1   40.8 [40.2,41.7]（不变）\n"
            "×100 20.8 [6.1,21.4]（腰斩）\n\n"
            "4b 配对差（×1）\n"
            "8 组合全部 = 0，CI[0,0]\n"
            "裕量退化 ≤ 4.5e-4\n"
            "单次注入 ≈ 0.04–0.09 μK\n\n"
            "4c 敏感性（55 mm/s）\n"
            "×1/×10：100%（裕量 0.84）\n"
            "×100：50%（0.42）\n"
            "×1000：50%（−0.29）\n\n"
            "4d 通道（×100）\n"
            "仅参量 0.480 ≈ 三通道 0.477\n"
            "仅反冲 1.000 / 仅指向 1.000\n"
            "参量注入 11.6 μK（唯一显著）\n\n"
            "ΔE=W_ext+W_noise+R_W\n残差 ≤6.6e-6 ≪ 1e-3",
            size=11.5)
add_textbox(slide, 0.7, 6.5, 8.2, 0.8,
            "四联图：左上 v95 曲线（关/×1/×100）；右上开/关配对差；左下幅度敏感性（俘获率+裕量双轴）；"
            "右下通道分解。方法差异如实报告：×100 档 logistic 拟合给 15.6 mm/s（悬崖过陡尾部偏缓），以非参 bootstrap 为主口径。",
            size=12, color=GRAY)

# ================================================================ 14 数值可信度
slide = add_slide()
title_bar(slide, "数值可信度：高速工况的步长收敛", "自适应步长规则 + 最快点 dt 减半检查")
add_picture_fit(slide, DEMO / "timestep_convergence_fast.png", 0.25, 1.0, 8.6, 5.3)
panel(slide, Inches(9.2), Inches(1.2), Inches(3.9), Inches(5.3),
      RGBColor(0xEA, 0xF3, 0xFB), BLUE)
add_textbox(slide, 9.45, 1.4, 3.4, 5.0,
            "步长规则\n\n"
            "dt ≤ min(\n  50 ns 基准,\n  T_move / 2000,\n  0.01·w / v_peak )\n\n"
            "实测：\n\n"
            "· δ = dt·v_peak/w\n  全程 ≤ 0.0015\n  （上限 0.01）\n\n"
            "· 最快点 500 mm/s：\n  dt = 3.43 ns\n  dt/2 标签翻转 = 0\n\n"
            "· 末态连续量随 dt\n  二阶下降（左图 log-log）\n\n"
            "· 功-能残差中位数\n  远低于 1e-3 阈值\n  （逐点写入 CSV）",
            size=13)
add_textbox(slide, 0.7, 6.45, 8.4, 0.8,
            "左：末态位置/速度偏差相对最细步长随 dt 的收敛 + 俘获率稳定性；右下转写自 preflight 实测值。"
            "所有扫描点实际 dt 与 δ 都写入 speed_scan_results.csv，可逐点核查。", size=13, color=GRAY)

# ================================================================ 15 结论
slide = add_slide()
title_bar(slide, "结论", "三问三答 + 诚实边界")
add_bullets(slide, 0.7, 1.3, 11.9, 3.9, [
    ("跟随上限（实验 1）", "纯移动 v_peak ≤ 63 mm/s（η≤0.24）零跟丢；100 mm/s 跟丢 75%；≥158 mm/s 全丢。"
     "激发标度 ε_med ∝ v^1.33 [1.08, 1.82]。"),
    ("端到端临界速度（实验 2）", "9 条曲线全部覆盖边界：v95(v_avg) = 41–56 mm/s（300 shots/点，Wilson CI）；"
     "换算峰值速度后三结构坍缩到 η_95 ≈ 0.29–0.30 —— 失效由峰值速度主导。"),
    ("恒速 vs 变速（实验 3）", "口径决定胜负：同 v_avg 低峰值者胜、同 v_peak 低均值者胜（均显著）；"
     "独立复核 30/30 CI 重叠，结论可靠。"),
    ("噪声 × 速度（实验 4）", "默认幅度（×1）对单次转移完全不可分辨（v95 逐点不变、配对差全 0）；"
     "×100 时参量通道独自使 v95 从 40.8 腰斩到 20.8 mm/s —— 边界对强度噪声的"
     "放大系数约 100 倍裕量。"),
    ("对既有工作的解释", "Level 2 优化波形 v_peak=11.5 mm/s（η=0.043），距边界约 7 倍余量 —— "
     "定量解释了 Level 2 的 100% 饱和；同时指出提速空间：按 v95 留 2 倍余量可到 v_avg≈27 mm/s（≈89 μs）。"),
], size=14)
bar = panel(slide, Inches(0.7), Inches(5.4), Inches(11.9), Inches(1.5),
            RGBColor(0xE8, 0xF0, 0xFE), BLUE)
add_textbox(slide, 1.0, 5.55, 11.3, 1.3,
            "边界与诚实声明：一维经典模型，不含 AOD 器件带宽/RF 扫描速率（动力学上限 ≠ 器件规格）、"
            "二维/柱面透镜（Level 3）、量子内态/相干性（Level 4-5）；demo 预算 300 shots/点（正式运行建议 ≥1000）；"
            "部分曲线 v50 因早停为 above_scanned；噪声幅度全部 assumed_sensitivity_only（非实验测量值），"
            "×100 档 bootstrap 与 logistic 估计不一致（20.8 vs 15.6 mm/s）已如实双报告；"
            "图片 9/9 程序化检查通过，未做人工视觉审阅。",
            size=12)

# ================================================================ 16 工程与复现
slide = add_slide()
title_bar(slide, "工程与复现", "新文件夹隔离 · 复用不复制 · 全链路可核查")
add_bullets(slide, 0.7, 1.3, 7.3, 5.2, [
    ("目录隔离", "simulation/level2b_speed_modify/（tag: level2 + speed_modify）；"
     "level0/1/2 为逐字节副本，与原工程校验和一致（8d9ae671cfe5b53f），原工程 git 零改动。"),
    ("复用不复制", "势/常数/积分器/采样器/Wilson 与配对 bootstrap/功-能核算/三通道噪声模型"
     "（build_mean_heating/scaled_heating）全部 import 原包；smootherstep 模式与 "
     "OverlappedWaveform 逐点一致有测试守护。"),
    ("测试", "105 项全部通过（103 快 + 2 slow CLI）：剖面校准/缩放不变性/口径换算/步长规则/"
     "临界速度回收/慢极限回归/断点续跑/噪声通道隔离与注入守恒/带噪功-能恒等式/CLI 冒烟（含 noise 阶段）/schema。"),
    ("输出", "21 个文件：preflight.json、speed_scan_results.csv（119 行）、scan_flags.npz、"
     "critical_speeds.csv、comparison_results.csv（60 行）、noise_speed_results.csv（61 行，0 失败）、"
     "noise_v95.json、metrics.json、9 张 PNG 等。"),
], size=13.5)
panel(slide, Inches(8.3), Inches(1.3), Inches(4.5), Inches(5.2),
      RGBColor(0xF2, 0xF2, 0xF2), GRAY)
add_textbox(slide, 8.55, 1.5, 4.0, 4.8,
            "复现命令\n\n"
            "$ cd simulation/level2b_speed_modify\n"
            "$ .venv/bin/python -m pytest\n"
            "→ 105 passed\n\n"
            "$ .venv/bin/python -m \\\n"
            "    level2b_speed_limit.level2b_cli \\\n"
            "    --config configs/level2b_speed_limit.yaml \\\n"
            "    --stage all\n\n"
            "阶段：preflight / scan /\ncompare / noise / all\n\n"
            "断点续跑：中断后重跑同命令，\n已完成扫描点自动跳过",
            size=12)

prs.save(OUT)
print(f"已生成 {OUT}（{len(prs.slides._sldIdLst)} 页）")
