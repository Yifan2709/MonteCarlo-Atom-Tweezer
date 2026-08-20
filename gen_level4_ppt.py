# -*- coding: utf-8 -*-
"""生成 Level 4 总结 PPT：用 outputs/level4_end_to_end_roundtrip/demo 的
23 张图 + metrics/summary 数据要点组装幻灯片。运行：
.venv/bin/python gen_level4_ppt.py
"""
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

REPO = Path(__file__).resolve().parent
FIG = REPO / "simulation" / "level2_joint_transfer" / "outputs" / \
    "level4_end_to_end_roundtrip" / "demo"
OUT = REPO / "Level4_端到端往返协议与半经典相干性_总结报告.pptx"

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
    from PIL import Image
    with Image.open(path) as image:
        iw, ih = image.size
    scale = min(w / iw, h / ih)
    nw, nh = iw * scale, ih * scale
    left = Inches(x + (w - nw) / 2)
    top = Inches(y + (h - nh) / 2)
    return slide.shapes.add_picture(str(path), left, top, Inches(nw), Inches(nh))


def add_table(slide, x, y, w, rows, col_widths=None, size=12, header=True):
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
            "Level 4：端到端往返协议\n与半经典相干性", size=40,
            color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 4.9, 11.3, 0.9,
            "pick-up — split — 375 μm 长运输 — return — merge — drop-off\n"
            "Protocol A/B 组合 Level 2 冻结 transfer 波形与 Level 3 冻结三维运输",
            size=16, color=RGBColor(0xBB, 0xD3, 0xE8), align=PP_ALIGN.CENTER)
add_textbox(slide, 1.0, 6.5, 11.3, 0.5,
            "preflight 20 项全通过 · 2000-shot 单轮 · 500-shot × 10 轮 · 128-shot × 40 轮长尾",
            size=13, color=RGBColor(0x88, 0xA8, 0xC8), align=PP_ALIGN.CENTER)

# ================================================================ 2 背景与定位
slide = add_slide()
title_bar(slide, "背景与 Level 4 的定位",
          "前三级分别隔离验证：Level 0 静态阱 / Level 1 顺序转移 / Level 2 同步波形 / Level 3 长距离运输")
add_bullets(slide, 0.55, 1.25, 6.1, 5.5, [
    ("首次组合", "把 Level 2 冻结 transfer 波形与 Level 3 冻结三维 AOD 运输"
     "组合成完整 pick-up—split—长运输—返回—merge—drop-off 往返"),
    ("论文锚点", "Fig. 5、Fig. 6、ED Fig. 10e、Methods『Atom transport』"
     "『Atom transfer between SLM and AOD tweezers』『Combined atom transfer and move』"),
    ("半经典相干性", "追踪每条经典轨迹上两超精细钟态的相对相位，"
     "比较无 DD / spin echo / 理想瞬时 XY4 的相干对比度"),
    ("诚实边界", "不含 Clifford twirling、IRB 拟合、Raman 散射、真空碰撞、"
     "脉冲误差、成像误差；结果不得称为论文 transport/IRB fidelity"),
])
add_table(slide, 6.95, 1.4, 5.9, [
    ["物理量", "数值"],
    ["原子", "¹³³Cs，m = 132.905 u"],
    ["SLM 阱深", "k_B × 180 μK（运输期 60 μK）"],
    ["AOD 阱深", "k_B × 280 μK"],
    ["SLM / AOD 束腰", "1.17 μm / 1.17 μm"],
    ["SLM / AOD 波长", "1061 nm / 1055 nm"],
    ["split 距离", "2.4 μm"],
    ["长运输距离", "375 μm（对角）"],
    ["初态温度", "5 μK"],
], col_widths=[2.3, 3.6], size=13)
add_textbox(slide, 6.95, 4.7, 5.9, 1.8,
            "成功定义：所有必需 checkpoint 满足预期阱归属\n"
            "（bound_to_slm / bound_to_aod / shared_or_ambiguous / unbound），\n"
            "且 final hold 后 E_SLM < 0。", size=13, color=GRAY)

# ================================================================ 3 协议时序
slide = add_slide()
title_bar(slide, "Protocol A / B：论文时序启发的端到端往返",
          "总时长 5654 μs / 3954 μs；9 个分段、8 个 checkpoint、XY4 / spin-echo 脉冲")
add_picture_fit(slide, FIG / "protocol_timeline.png", 0.55, 1.1, 12.3, 3.9)
add_bullets(slide, 0.55, 5.15, 12.3, 1.9, [
    ("分段锚点", "initial_slm_hold(100μs) → pickup(850μs) → split(400μs) → "
     "outbound_long_move(1450μs) → remote_hold(54μs) → inbound_long_move(1450μs)"
     " → merge(400μs) → dropoff(850μs) → final_slm_hold(100μs)"),
    ("深度 ramp", "pickup/dropoff 用有零端点斜率的 monotone smootherstep（模型假设，"
     "非论文逐点波形）；split/merge 用 Level 3 已验证 piecewise constant-jerk；long move 用 adiabatic-sine"),
    ("Protocol B", "pickup/dropoff 用 Level 2 冻结最佳 dropoff 波形的时间反向，"
     "不重复 split/merge；SLM 恒 140 μK 与 Protocol A 的 180→60 μK 是两个不同 baseline"),
], size=13)

# ================================================================ 4 控制量与透镜
slide = add_slide()
title_bar(slide, "AOD 深度、中心、速度与 lensing 焦点偏移",
          "控制量由分段解析公式拼接，端点/导数/时间连续性经 preflight 第 3 项验证")
add_picture_fit(slide, FIG / "trap_path_and_depths.png", 0.55, 1.15, 12.2, 4.3)
add_bullets(slide, 0.55, 5.6, 12.3, 1.6, [
    ("复用 Level 3 模型", "AOD crossed-AOD lensing 势 U = −D·G₁·G₂，"
     "焦点偏移 z_s/z_R = 2σ·ċ/v_s 沿用冻结 Level 3 nominal v_s = 0.5392 m/s（severity 1.0）"),
    ("split→long move 衔接", "split 方向与 long-move 对角方向一致；outbound 终点为 AOD "
     "相对起始 SLM 的 2.4+375 μm 路径位置；inbound/merge/dropoff 为严格时间反向控制序列"),
], size=13)

# ================================================================ 5 初态与冻结
slide = add_slide()
title_bar(slide, "初态采样与上游冻结输入",
          "谐振子提议 + E<0 束缚拒绝；上游产物经 SHA-256 冻结并写入 upstream_manifest.json")
add_picture_fit(slide, FIG / "initial_ensemble.png", 0.55, 1.15, 8.0, 4.6)
add_picture_fit(slide, FIG / "upstream_provenance.png", 8.7, 1.3, 4.2, 4.4)
add_bullets(slide, 0.55, 5.9, 12.3, 1.2, [
    ("复现性", "validation pool（seed=41001，T=5 μK，D=180 μK）接受率 1.000；"
     "明确标注非完整有限深高斯阱的严格正则系综"),
    ("冻结产物", "Level 2 best_waveform.yaml + Level 2/3 metrics/config 均记录"
     "绝对路径、SHA-256、生成时间；reconstructed_from_config=false"),
], size=13)

# ================================================================ 6 preflight
slide = add_slide()
title_bar(slide, "preflight：20 项自检，失败即停止大样本",
          "先跑 Level 0–3 全部测试，再验证上游、连续性、回归、反向性、相位解析与步长收敛")
add_picture_fit(slide, FIG / "preflight_checks.png", 0.55, 1.15, 12.2, 4.6)
add_bullets(slide, 0.55, 5.9, 12.3, 1.2, [
    ("关键验证", "时间反向 + 速度反转误差 < 10⁻¹² m；Level 2 三维切片俘获率 1.0、"
     "Level 3 长运输留阱率 1.0；解析相位相对误差 < 10⁻⁹；分段功-能残差 / 阱深 < 10⁻⁴"),
    ("步长收敛", "0.10/0.05/0.025 μs 配对复算：标签翻转 0、max|Δpos| < 0.05 μm、"
     "max|ΔW| < 1 μK、max|Δφ| < 0.05 rad——repeated 用 0.10 μs 主结果合规"),
], size=13)

# ================================================================ 7 basin 分类
slide = add_slide()
title_bar(slide, "重叠阱 basin 归属分类",
          "两阱都显著时沿连线找总势最低点与分隔势垒；能量接近鞍点或深度接近零时标 ambiguous")
add_picture_fit(slide, FIG / "basin_classification.png", 0.55, 1.15, 12.2, 4.6)
add_bullets(slide, 0.55, 5.9, 12.3, 1.2, [
    ("四分类输出", "bound_to_slm / bound_to_aod / shared_or_ambiguous / unbound；"
     "无法稳定归类时标 ambiguous，不强行判成功；near-threshold 单独记录、不静默删除"),
    ("关键 checkpoint", "split_end（阱间距 2.40 μm）两阱均显著、atoms 主要在 AOD 区；"
     "merge_end（阱间距 0.00 μm）合并阱标 shared_or_ambiguous；ambiguous 比例 ≈ 25%"),
], size=13)

# ================================================================ 8 代表轨迹
slide = add_slide()
title_bar(slide, "代表轨迹：成功、失败阶段与重俘获",
          "压力池（50 μK，诊断种子）采集各类 first-failure 阶段代表；位置 + 能量 + 累积功")
add_picture_fit(slide, FIG / "representative_roundtrips.png", 0.55, 1.15, 12.2, 4.9)
add_bullets(slide, 0.55, 6.15, 12.3, 1.0, [
    ("成功轨迹", "AOD 携带原子沿对角往返 375 μm，能量在 −180 μK 左右振荡、末态回到 SLM"),
    ("失败/重俘获", "split 阶段失败的原子弹道发散、总机械能转正；recaptured 在后续 checkpoint "
     "恢复 SLM 束缚——三类协议结果（final_retained / roundtrip_success / absorbing_success）据此区分"),
], size=13)

# ================================================================ 9 单轮验证
slide = add_slide()
title_bar(slide, "Stage B：单轮 validation（2000 shots × 4 条件，dt=0.05 μs）",
          "validation_pool 一次性使用；同初态配对比较；冻结条件写入 frozen_validation_protocols.yaml")
add_table(slide, 0.55, 1.2, 7.6, [
    ["条件", "roundtrip success", "Wilson 95% CI"],
    ["protocol_a_nominal_lensing", "2000/2000 (1.0000)", "[0.9981, 1.0000]"],
    ["protocol_a_lensing_off", "2000/2000 (1.0000)", "[0.9981, 1.0000]"],
    ["protocol_b_frozen", "2000/2000 (1.0000)", "[0.9981, 1.0000]"],
    ["static_idle_same_duration", "2000/2000 (1.0000)", "[0.9981, 1.0000]"],
], col_widths=[3.8, 1.9, 1.9], size=13)
add_textbox(slide, 0.55, 3.4, 7.6, 1.5,
            "全部 4 条件 roundtrip success 1.0000，final retained 1.0000，recaptured 0。\n"
            "3 组配对差异（lensing_off vs nominal、protocol_b vs a、static_idle vs a）"
            "均为 0、CI 含 0——5 μK 初态下经典动力学深度安全，差异只在相位中显现。",
            size=14, color=DARK)
add_picture_fit(slide, FIG / "checkpoint_survival.png", 8.35, 1.2, 4.6, 3.4)
add_picture_fit(slide, FIG / "failure_stage_breakdown.png", 8.35, 4.7, 4.6, 2.1)
fig_caption(slide, "右上图：四个条件逐 checkpoint 累计/条件成功概率；右下图：单轮 first-failure 阶段堆叠（nominal 无失败）")

# ================================================================ 10 配对比较 + 能量演化
slide = add_slide()
title_bar(slide, "配对比较与 checkpoint 能量演化",
          "同初态配对四格计数与 bootstrap 95% CI；checkpoint 总能量分位数带随协议阶段演变")
add_picture_fit(slide, FIG / "paired_comparisons.png", 0.55, 1.15, 6.2, 4.9)
add_picture_fit(slide, FIG / "checkpoint_energy_evolution.png", 6.95, 1.15, 6.1, 4.9)
add_bullets(slide, 0.55, 6.2, 12.3, 1.0, [
    ("配对差异均为 0", "四格配对计数 2000/2000 全成功，无任何配对差异；"
     "checkpoint 能量带显示 pickup 后 SLM 变浅、transport 期 AOD 携带、merge 后回 SLM"),
    ("配对 bootstrap 机制", "同一初态重采样，报告差异点估计与 95% CI；"
     "paired_counts 四格（both/only_a/only_b/none）用于判别显著性"),
], size=13)

# ================================================================ 11 功-能账本
slide = add_slide()
title_bar(slide, "分段功-能账本：ΔE = W_ext + R_W",
          "外部功四通道分解（移动 / lensing 偏移 / AOD 深度 / SLM 深度）；静态 hold 段只有数值误差")
add_picture_fit(slide, FIG / "energy_ledger.png", 0.55, 1.15, 8.2, 5.3)
add_bullets(slide, 8.95, 1.4, 4.1, 5.0, [
    ("逐段一致", "pickup ΔE=+58.3 μK ↔ W_ext=+58.3 μK；merge ΔE=+154.9 μK ↔ "
     "W_ext=+154.9 μK；每段与全周期均检查，不靠正负抵消"),
    ("残差量级", "max|R_W| ≈ 2–4 × 10⁻³ μK，mean|R_W| ≈ 4 × 10⁻⁴ μK；"
     "相对 280 μK 阱深 < 1.5 × 10⁻⁵，验证 velocity-Verlet 时变精度"),
    ("静态 hold", "initial_slm_hold / final_slm_hold 的 W_ext≈0，残差即数值误差"),
], size=14)

# ================================================================ 12 半经典相位
slide = add_slide()
title_bar(slide, "半经典差分光移相位与 DD toggling",
          "φ_i = ∫ y(t)·[η_SLM·U_SLM + η_AOD·U_AOD]/ħ dt；A_SLM/A_AOD 分开累积，η 可线性重算")
add_picture_fit(slide, FIG / "dd_toggling_and_phase.png", 0.55, 1.15, 12.2, 4.7)
add_bullets(slide, 0.55, 6.0, 12.3, 1.1, [
    ("toggling 函数", "y(t) = ±1 表示理想瞬时 π 脉冲的符号翻转；"
     "脉冲严格落在时间样本内部（≥0.2·dt 远离边界），分段边界自动安全"),
    ("脉冲网格拆分", "步内脉冲把步拆成两段解析积分；_pulses_to_steps 验证"
     "无重复/漏积分；XY4 的 X/Y 轴顺序在纯退相位模型中等价"),
], size=13)

# ================================================================ 13 相位分布 + 对比度
slide = add_slide()
title_bar(slide, "相位分布与各 DD 模式对比度随轮次",
          "条件相干对比度 C_cond = |mean(e^{iφ})|；500 shots × 10 rounds")
add_picture_fit(slide, FIG / "phase_and_contrast.png", 0.55, 1.15, 7.4, 4.3)
add_picture_fit(slide, FIG / "usable_coherent_fraction.png", 8.2, 1.2, 4.8, 4.1)
add_table(slide, 0.55, 5.6, 7.4, [
    ["DD 模式", "C_cond（第 1 轮）", "C_cond（第 10 轮）"],
    ["none", "0.8957", "0.0568"],
    ["spin_echo", "1.0000", "0.0351"],
    ["xy4_per_long_move", "0.9735", "0.0496"],
], col_widths=[2.4, 2.5, 2.5], size=12)
add_textbox(slide, 8.2, 5.6, 4.8, 1.5,
            "单轮 nominal：none C=0.8797，echo C=1.0000，XY4 C=0.9695。\n"
            "usable_coherent_fraction = P(success)×C_cond 仅为工程代理，非量子过程保真度。",
            size=13, color=GRAY)

# ================================================================ 14 分段相位分解
slide = add_slide()
title_bar(slide, "分段相位贡献分解",
          "各 DD 模式按分段的 mean|φ_seg| 与 std（条件于 roundtrip success）")
add_picture_fit(slide, FIG / "segment_phase_breakdown.png", 0.55, 1.15, 12.2, 4.6)
add_bullets(slide, 0.55, 5.9, 12.3, 1.2, [
    ("相位来源", "长运输段贡献最大（mean|φ| ≈ 6.7 rad）；pickup/dropoff/merge 次之；"
     "remote_hold 54 μs 贡献最小；spin_echo 对整体对称 quasistatic shift 消除最好"),
    ("分段分解的价值", "相位按段累积，可定位退相干主要来自长运输还是 transfer；"
     "为 Level 5 的有限脉冲与真实噪声谱提供分段诊断基准"),
], size=13)

# ================================================================ 15 多轮生存 + 长尾
slide = add_slide()
title_bar(slide, "Stage C：重复往返 survival 与 40 轮长尾",
          "500 shots × 10 rounds（Protocol A/B）+ 128 shots × 40 rounds 长尾（Protocol A）")
add_picture_fit(slide, FIG / "survival_vs_round.png", 0.55, 1.15, 7.6, 4.6)
add_picture_fit(slide, FIG / "long_tail_survival.png", 8.35, 1.15, 4.6, 4.6)
add_bullets(slide, 0.55, 6.0, 12.3, 1.1, [
    ("absorbing 语义", "一旦某 checkpoint 失败，该 shot 后续轮次永久记为失败；"
     "10 轮 absorbing survival 1.0000（两协议），非吸收 final retention 同为 1.0000"),
    ("长尾稳健", "128 shots × 40 rounds 的 absorbing survival 1.0000；"
     "经验生存拟合（clipped-Boltzmann / exponential）因数据几乎全 1、参数不可辨识而安全跳过"),
], size=13)

# ================================================================ 16 加热
slide = add_slide()
title_bar(slide, "多轮加热：留阱 shots 的 SLM 激发随轮次",
          "仅成功/留阱 shots；中位数、10%/90%/95% 分位数；在线性拟合残差合理时才报告斜率")
add_picture_fit(slide, FIG / "heating_vs_round.png", 0.55, 1.15, 8.2, 4.6)
add_bullets(slide, 9.0, 1.4, 4.0, 4.6, [
    ("Protocol A", "第 1 轮中位激发 12.9 μK，第 10 轮 13.0 μK；"
     "90% 分位 24.2 μK，95% 分位 28.1 μK；斜率 ≈ +0.01 μK/轮"),
    ("Protocol B", "第 1 轮中位激发 51.5 μK，第 10 轮 51.4 μK；"
     "略高于 A（SLM 恒 140 μK，运输期深度未先降后升）"),
    ("混用警告", "未留阱粒子的正能量不与留阱态激发混入同一均值"),
], size=14)

# ================================================================ 17 消融
slide = add_slide()
title_bar(slide, "Stage D：消融归因（不改变预注册主协议）",
          "transfer-only / transport-only / no-lensing / no-DLS / no-noise / static-idle / no-remote-hold / constant-jerk-long-move")
add_picture_fit(slide, FIG / "ablation_comparison.png", 0.55, 1.15, 12.2, 4.4)
add_bullets(slide, 0.55, 5.7, 12.3, 1.5, [
    ("经典成功率全 1.0", "8 个消融条件 roundtrip success 均 1.0（Wilson 下界 ≥0.9924）；"
     "transfer-only 中位激发 13.1 μK，transport-only −83.5 μK（阱深定义差异）"),
    ("对比度差异", "transfer-only C(none)=C(xy4)=0.988（无长运输相位）；"
     "no-DLS C=1.0（相位恒为 0）；删除 remote_hold 对经典成功率无影响"),
    ("归因纪律", "消融只用于归因，不从消融中挑一个更好协议后重新命名为预注册主协议"),
], size=13)

# ================================================================ 18 敏感性
slide = add_slide()
title_bar(slide, "Stage E：one-factor-at-a-time 敏感性（500 shots/点）",
          "温度、深度、对准、lensing severity、η、强度噪声、DD 脉冲时刻")
add_picture_fit(slide, FIG / "sensitivity_panels.png", 0.55, 1.15, 12.2, 4.3)
add_bullets(slide, 0.55, 5.6, 12.3, 1.6, [
    ("经典稳健区", "温度 3–10 μK、SLM/AOD 深度 ±10%、横向对准 ±0.10 μm、lensing "
     "severity 0.75/1.25 下 roundtrip success 均 1.0000"),
    ("相位敏感项", "η 从 1.2→1.6 × 10⁻⁴ 时 C(none) 从 0.896→0.829；"
     "assumed 噪声 rms=0.05 quasistatic 时 C(none)=0.393（仅敏感性演示，非 nominal）"),
    ("DD 脉冲时刻", "±0.5 μs 整体平移对 success 与 C_cond 无影响（理想瞬时脉冲模型）"),
], size=13)

# ================================================================ 19 步长收敛
slide = add_slide()
title_bar(slide, "步长收敛与数值验证",
          "状态、外部功、相位、success 标签随 dt 收敛；0.10 μs 主结果与 0.05/0.025 μs 一致")
add_picture_fit(slide, FIG / "timestep_convergence.png", 0.55, 1.15, 8.2, 4.6)
add_bullets(slide, 9.0, 1.4, 4.0, 4.6, [
    ("二阶收敛", "max|Δpos|、max|Δvel|、max|ΔW|、max|Δφ| 随 dt 二阶下降；"
     "0.05 μs 相对 0.025 μs 标签翻转 0"),
    ("主结果合规", "若用 0.10 μs 做多轮主结果，必须在固定子集上证明 survival、"
     "heating、phase、work residual 与 0.05/0.025 μs 一致到预设容差——preflight 第 19 项已通过"),
    ("网格处理", "DD 脉冲时刻经 _snap_off_grid 处理，保证与最近时间样本距离 ≥ 0.2·dt"),
], size=14)

# ================================================================ 20 噪声模型
slide = add_slide()
title_bar(slide, "可选技术噪声模型（assumed，仅敏感性模式）",
          "quasistatic 每 shot 一次抽样 / Ornstein-Uhlenbeck 时变；SLM/AOD 可独立或相关")
add_picture_fit(slide, FIG / "noise_model_diagnostic.png", 0.55, 1.15, 12.2, 4.3)
add_bullets(slide, 0.55, 5.6, 12.3, 1.5, [
    ("参数来源", "RMS、相关时间、相关系数均来自配置并标记 assumed，"
     "不得声称由论文唯一确定；噪声种子与初态种子分离"),
    ("OU 更新", "ε(t+dt) = ε·exp(−dt/τ) + √(1−exp(−2dt/τ))·σ·n，"
     "边界处从 quasistatic 抽样热启动；测试验证均值/方差/相关时间与 seed 复现"),
    ("nominal 不含噪声", "nominal 运行默认不加入随机技术噪声；"
     "敏感性网格只用于检验对比度对假设噪声的依赖"),
], size=13)

# ================================================================ 21 软件工程
slide = add_slide()
title_bar(slide, "实现与质量保障",
          "simulation/level2_joint_transfer/src/level4_roundtrip_coherence/，23 张诊断图全程序化验收")
add_table(slide, 0.55, 1.2, 12.2, [
    ["模块", "职责"],
    ["protocol_segments.py / upstream_artifacts.py", "分段时序与连续拼接；上游产物发现、hash、冻结与 manifest"],
    ["roundtrip_simulation.py / trap_assignment.py", "单轮/多轮向量化推进；重叠阱 basin 与 checkpoint 归属"],
    ["energy_ledger.py / dephasing.py / dd_sequences.py", "分段功-能账本；差分光移相位、toggling、contrast；none/echo/XY4/custom 脉冲"],
    ["noise_models.py / roundtrip_statistics.py", "quasistatic/OU 强度噪声；Wilson、配对 bootstrap、生存拟合、加热斜率"],
    ["preflight4.py / sampling4.py", "20 项预检查；三维谐振提议 + 束缚拒绝初态采样"],
    ["level4_cli.py / level4_visualization.py", "CLI 五个 stage；23 张诊断图 + 程序化验收记录（有限性/图例/坐标覆盖/视觉审阅）"],
], col_widths=[4.2, 8.0], size=12)
add_bullets(slide, 0.55, 4.6, 12.2, 2.4, [
    ("测试", "test_level4.py 共 32 项测试函数，覆盖 prompt §17 的 30 项要求："
     "分段连续/控制 roundtrip/Protocol B 不重复/L2-L3 回归/时间反向/总势与力一致/basin 分类/"
     "checkpoint 判据/first-failure 状态机/recaptured 区分/多轮状态连续/功-能/静态 hold/"
     "解析相位/toggling/脉冲网格/contrast 案例/生存条件/噪声统计/种子隔离/共同随机数/"
     "冻结无泄漏/统计/拟合安全失败/步长收敛/near-threshold/CLI 冒烟/schema/sources 保护"),
    ("产物", "23 张 PNG + CSV×10 + JSON×6 + Markdown×1 + YAML×2，"
     "全部程序化复检 schema、有限值、键唯一性与图例一致性；视觉审阅 23/23 通过"),
], size=13)

# ================================================================ 22 结论与边界
slide = add_slide()
title_bar(slide, "结论、边界与 Level 5 建议")
add_bullets(slide, 0.55, 1.2, 6.1, 5.6, [
    ("经典动力学深度安全", "5 μK 初态下 Protocol A/B 单轮 2000/2000、10 轮 500/500、"
     "40 轮长尾 128/128 的 roundtrip success 均 1.0000，配对差异为 0"),
    ("差异只在相位中显现", "spin_echo 对整体对称 shift 消除最好（C≈1.0000）；"
     "XY4-per-long-move 次之（C≈0.9695→0.0496）；none 退相干最快"),
    ("长运输是主要退相干源", "分段相位分解显示长运输段贡献 mean|φ| ≈ 6.7 rad，"
     "远高于 transfer 段；为 Level 5 的有限脉冲与真实噪声谱提供分段诊断基准"),
    ("消融与敏感性稳健", "8 个消融条件与温度/深度/对准/lensing 敏感性下 success 均 1.0；"
     "η 与 assumed 强度噪声是主要相位敏感项"),
], size=14)
add_textbox(slide, 6.95, 1.3, 5.9, 2.6,
            "边界声明\n\n● 经典成功概率 ≠ 论文 IRB/transport fidelity\n"
            "● C_cond 与 usable_coherent_fraction ≠ 量子过程保真度\n"
            "● η、噪声参数为假设值；lensing severity 为无量纲场景\n"
            "● remote 54 μs 仅为相位积累 hold，不输出 gate fidelity",
            size=14, color=DARK)
add_textbox(slide, 6.95, 4.1, 5.9, 2.8,
            "Level 5 建议\n\n1. 有限时长/含脉冲误差的 Rabi 驱动与完整自旋哈密顿量\n"
            "2. 真实噪声谱（磁场、强度 OU 谱）与 pulse-area error\n"
            "3. Clifford/IRB 随机基准协议拟合\n"
            "4. 多原子并行操作的相互作用与成像误差",
            size=14, color=DARK)

prs.save(OUT)
print(f"saved: {OUT} ({len(prs.slides._sldIdLst)} slides)")
