# -*- coding: utf-8 -*-
"""生成 Level 2 总结 PPT（简明版）：一页一张图 + 一行说明。

图片由 visualize_level2_slides.py 从 demo 输出数据生成（level2_slide_figures/）。
"""
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

REPO = Path(__file__).resolve().parent
FIG = REPO / "level2_slide_figures"
OUT = REPO / "Level2_同步波形优化_总结报告.pptx"

BLUE = RGBColor(0x1F, 0x77, 0xB4)
DARK = RGBColor(0x22, 0x2A, 0x35)
GRAY = RGBColor(0x66, 0x66, 0x66)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]


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


def figure_slide(figure, title, caption):
    """一页一张图：标题栏 + 大图 + 底部一行说明。"""
    slide = prs.slides.add_slide(BLANK)
    bar = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width, Inches(0.86))
    bar.fill.solid()
    bar.fill.fore_color.rgb = BLUE
    bar.line.fill.background()
    box = slide.shapes.add_textbox(Inches(0.5), Inches(0.10), Inches(12.4), Inches(0.66))
    box.text_frame.word_wrap = True
    paragraph = box.text_frame.paragraphs[0]
    run = paragraph.add_run()
    run.text = title
    run.font.size = Pt(25)
    run.font.bold = True
    run.font.color.rgb = WHITE
    add_picture_fit(slide, FIG / figure, 0.30, 1.00, 12.73, 5.35)
    note = slide.shapes.add_textbox(Inches(0.7), Inches(6.55), Inches(11.9), Inches(0.6))
    note.text_frame.word_wrap = True
    note.text_frame.paragraphs[0].alignment = PP_ALIGN.CENTER
    run = note.text_frame.paragraphs[0].add_run()
    run.text = caption
    run.font.size = Pt(14)
    run.font.color.rgb = GRAY


# ================================================================ 封面
slide = prs.slides.add_slide(BLANK)
bg = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width, prs.slide_height)
bg.fill.solid()
bg.fill.fore_color.rgb = RGBColor(0x0E, 0x2A, 0x47)
bg.line.fill.background()


def cover_text(x, y, w, h, text, size, color, bold=False):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    box.text_frame.word_wrap = True
    box.text_frame.paragraphs[0].alignment = PP_ALIGN.CENTER
    run = box.text_frame.paragraphs[0].add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = color


cover_text(1.0, 1.55, 11.3, 0.5, "Monte Carlo for AOD→SLM 原子转移", 18,
           RGBColor(0x9E, 0xC5, 0xE8))
cover_text(1.0, 2.25, 11.3, 1.0, "Level 2：同步波形优化", 40, WHITE, bold=True)
cover_text(1.0, 3.35, 11.3, 0.5, "600 μs 顺序基线 → 400 μs 同步重叠波形的鲁棒蒙特卡洛验证", 17,
           RGBColor(0xBB, 0xD3, 0xE8))

chips = [
    (1.55, "2000 / 2000", "独立验证俘获率"),
    (5.30, "257", "优化候选（192 LHS + 65 精修）"),
    (9.05, "10 / 10", "preflight 预检查通过"),
]
for x, big, small in chips:
    chip = slide.shapes.add_shape(5, Inches(x), Inches(4.35), Inches(2.75), Inches(1.35))
    chip.fill.solid()
    chip.fill.fore_color.rgb = RGBColor(0x1B, 0x3A, 0x5C)
    chip.line.color.rgb = RGBColor(0x2C, 0x5A, 0x8C)
    cover_text(x, 4.50, 2.75, 0.6, big, 26, WHITE, bold=True)
    cover_text(x, 5.18, 2.75, 0.4, small, 12, RGBColor(0x9E, 0xC5, 0xE8))

cover_text(1.0, 6.45, 11.3, 0.4, "结论预览：5 μK 工况处于饱和区，三波形俘获率均 100%，如实报告", 13,
           RGBColor(0x88, 0xA8, 0xC8))

# ================================================================ 一页一图
figure_slide(
    "s01_问题.png", "问题：600 μs 能压到 400 μs 吗？",
    "Level 1 先移动后降深（600 μs）；Level 2 让移动与降深同步重叠，压缩到 400 μs（提速 33%）。")

figure_slide(
    "s02_物理模型.png", "物理模型：一维双高斯阱",
    "280 μK AOD 阱载着 5 μK 铯原子移向 140 μK 静止 SLM 阱并逐渐变浅；末态 E_f < 0 判定俘获。")

figure_slide(
    "s03_方法.png", "方法：三池隔离 + 公共随机数",
    "三批互不重叠的样本池各司其职，validation 只报告一次；所有波形用同一组初态（CRN）比较。")

figure_slide(
    "s04_波形对比.png", "三类波形：『同步』= 重叠 39%",
    "优化波形移动几乎贯穿全程（s ≥ 0.019），降深 s ≥ 0.607 才开始——移动与降深重叠约 39% 的时间。")

figure_slide(
    "s05_势场演化.png", "势场演化：原子被『倒』进 SLM",
    "AOD 阱边移动边变浅，阱底（原子跟随）一路滑向 SLM；交接时刻 AOD 深度恰好为 0。")

figure_slide(
    "s11_噪声模型.png", "噪声模型：三通道常数功率加热",
    "系综平均功率的确定性注入（无随机性）：反冲 / 强度（线性化）/ 指向三通道均为常数功率，注入能量计入功-能账本 ΔE = W_ext + W_noise + R_W。")

figure_slide(
    "s12_噪声结果.png", "噪声 × 速度：×1000 才有影响，参量主导",
    "默认 ×1 注入 ≤ 0.2 μK：v95 逐点不变（40.8 mm/s）；×100 无显著移动（42.6，CI 重叠）；×1000 全程 < 0.95，@55 mm/s 俘获率 0.50，其中仅参量通道（0.477）≈ 三通道（0.477），仅反冲/仅指向仍 1.00。")

figure_slide(
    "s13_最高速度测算.png", "最高速度测算：硬件条件下的 v95（噪声 ×1）",
    "v95 随 AOD 深度按 √D 缩放（140/280/560 μK → 29.6/41.8/58.1 mm/s），η95 ≈ 0.29 与深度无关（v_ad ∝ √D）；束腰收窄到 0.90 μm 提速到 54.9 mm/s（阱更硬），放宽到 1.50 μm 降到 38.4 mm/s。")

figure_slide(
    "s14_同时长方案对比.png", "同时长移动方案：峰值速度低者胜",
    "AOD 起动→停止时长与位移完全相同（44/39/35 μs × 2.4 μm，噪声 ×1）：快加速梯形（v_peak=1.18·v_avg）0.99/0.96/0.78 全面最优，三角（2.0×）0.48/0.09/0.00 最差——俘获由瞬时峰值速度决定，加速平缓与否是次要因素。")

prs.save(OUT)
print(f"已生成 {OUT}（{len(prs.slides._sldIdLst)} 页）")
