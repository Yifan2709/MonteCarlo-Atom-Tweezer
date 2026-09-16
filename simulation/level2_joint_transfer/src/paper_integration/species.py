"""三个模型配置体系：cs_existing / rb87_2022_native / yb171_2026_native。

物种与实验参数严格分开（任务书 §3）。每个参数带溯源（§5）。
数值来源：R = Bluvstein et al., Nature 604, 451 (2022)（arXiv:2112.03923）；
Y = Zhang et al., Nat. Phys. 22, 910 (2026)（arXiv:2506.13724v2）；
C = Manetsch et al.（arXiv:2403.12021v4，仓库现有对标）。
注：本环境网络受限，无法下载 Y 的 Zenodo 源数据与 R 的期刊版补充材料，
论文数值取自 arXiv 公开版本正文/Methods/ED；未能获取的参数显式标 assumed。
"""
from __future__ import annotations

import math

from .provenance import ParamSpec, Provenance, register

KB = 1.380649e-23
H_PLANCK = 6.62607015e-34
U_TO_KG = 1.66053906660e-27


def _mk(config, name, value, unit, prov, source, unc=None, note=""):
    return config.add(ParamSpec(name, float(value), unit, Provenance(prov),
                                source, unc, note))


# ----------------------------------------------------------------------------
# cs_existing：保留原 Cs 配置（Manetsch / 仓库现有），用于 B0 回归
# ----------------------------------------------------------------------------
cs = register("cs_existing")
_mk(cs, "mass_u", 132.90545196, "u", "external_measured", "Cs-133 同位素质量")
_mk(cs, "aod_wavelength_m", 1055e-9, "m", "paper_measured", "C Ext. Data Fig.2")
_mk(cs, "aod_waist_m", 1.17e-6, "m", "paper_measured", "C AOD 阱参数")
_mk(cs, "aod_depth_uK", 280.0, "uK", "paper_measured", "C Methods/Fig.6")
_mk(cs, "slm_depth_uK", 140.0, "uK", "paper_measured", "C Methods（转移时）")
_mk(cs, "temperature_uK", 25.0, "uK", "assumed_sensitivity",
    "诊断假设，非实验标定；初温-加热耦合见 v10 报告")
_mk(cs, "move_distance_m", 2.4e-6, "m", "paper_measured", "C Fig.6 caption")

# ----------------------------------------------------------------------------
# rb87_2022_native：R 论文 87Rb 移动光镊运输
# ----------------------------------------------------------------------------
rb = register("rb87_2022_native")
_mk(rb, "mass_u", 86.909180531, "u", "external_measured", "Rb-87 同位素质量")
_mk(rb, "wavelength_m", 828e-9, "m", "paper_measured", "R Methods：AOD 光镊 828 nm")
_mk(rb, "depth_loading_hz", 16e6, "Hz(U/h)", "paper_measured",
    "R Methods：装载阱深 2π×16 MHz")
_mk(rb, "depth_circuit_hz", 4e6, "Hz(U/h)", "paper_measured",
    "R Methods：线路运行阱深 2π×4 MHz（~100 能级）")
_mk(rb, "omega_r_circuit_rad_s", 2 * math.pi * 40e3, "rad/s", "paper_measured",
    "R Methods：线路径向频率 2π×40 kHz")
_mk(rb, "omega_axial_rad_s", 2 * math.pi * 6e3, "rad/s", "paper_measured",
    "R Methods：轴向 2π×6 kHz")
# paper_derived：由 U=h×4 MHz 与 ω_r 反推束腰 w=sqrt(4U/(m ω²))
_m = 86.909180531 * U_TO_KG
_U = 4e6 * H_PLANCK
_w = math.sqrt(4 * _U / (_m * (2 * math.pi * 40e3) ** 2))
_mk(rb, "waist_m", _w, "m", "paper_derived",
    "w=sqrt(4U/(m·ω_r²)), U=h·4MHz, ω_r=2π·40kHz", note="≈1.08 μm")
_mk(rb, "zr_m", math.pi * _w**2 / 828e-9, "m", "paper_derived",
    "zR=π w²/λ")
_mk(rb, "bell_move_distance_m", 110e-6, "m", "paper_measured",
    "R 图1：Bell 对分离 110 μm")
_mk(rb, "bell_move_time_s", 300e-6, "s", "paper_measured", "R 图1：300 μs")
_mk(rb, "bell_flat_speed_m_per_s", 0.55e-6 / 1e-6, "m/s", "paper_measured",
    "R 图1d：b≲0.55 μm/μs 保真度不受影响")
_mk(rb, "bell_raw_fidelity", 0.948, "", "paper_measured", "R：94.8(2)%")
_mk(rb, "t2_star_s", 4e-3, "s", "paper_measured", "R Methods")
_mk(rb, "t2_dd_s", 1.49, "s", "paper_measured", "R Methods：XY16，unc 0.08 s", 0.08)
_mk(rb, "t1_s", 4.0, "s", "paper_measured", "R Methods")
_mk(rb, "cz_layer_error_y", 0.002, "", "paper_measured",
    "R ED 线路 MC：CZ 层 Y 误差 0.2%/qubit")
_mk(rb, "cz_layer_error_x", 0.002, "", "paper_measured", "同上 X 0.2%")
_mk(rb, "cz_layer_error_z", 0.005, "", "paper_measured", "同上 Z 0.5%")
_mk(rb, "cz_layer_error_loss", 0.005, "", "paper_measured", "同上 loss 0.5%")
_mk(rb, "oneq_ambient_x", 0.001, "", "paper_measured", "R ED：1Q 环境误差 X 0.1%")
_mk(rb, "oneq_ambient_y", 0.001, "", "paper_measured", "R ED：1Q 环境误差 Y 0.1%")
_mk(rb, "oneq_ambient_z", 0.004, "", "paper_measured", "R ED：1Q 环境误差 Z 0.4%")
_mk(rb, "oneq_ambient_loss", 0.002, "", "paper_measured", "R ED：1Q loss 0.2%")
_mk(rb, "initial_loss_prob", 0.01, "", "paper_measured", "R ED：初态 loss 1%")
_mk(rb, "rydberg_omega_hz", 3.6e6, "Hz", "paper_measured", "R Methods：Ω/2π=3.6 MHz")
_mk(rb, "rydberg_delta_hz", -1.36e6, "Hz", "paper_derived",
    "Δ=−0.377371×2π MHz（R Methods 数值）")
_mk(rb, "rydberg_xi_hz", -0.621089e6, "Hz", "paper_measured", "R Methods ξ")
_mk(rb, "cz_duration_s", 0.683201 / 3.6e6, "s", "paper_derived",
    "τ=0.683201/Ω，≈190 ns")
_mk(rb, "steane_xl_raw", 0.71, "", "paper_measured", "R 图3：0.71(2)", 0.02)
_mk(rb, "steane_xl_corrected", 0.991, "", "paper_measured", "R 图3")
_mk(rb, "steane_no_error_frac", 0.66, "", "paper_measured", "R：66(1)%", 0.01)
_mk(rb, "nmax_analytic_low", 26.0, "", "paper_measured",
    "R Methods：erf 留存模型 Nmax≈26")
_mk(rb, "nmax_analytic_high", 33.0, "", "paper_measured", "R Methods：Nmax 上界 33")
_mk(rb, "trap_drop_s", 500e-9, "s", "paper_measured", "R：每次 CZ 关阱 500 ns")

# ----------------------------------------------------------------------------
# yb171_2026_native：Y 论文 171Yb 亚稳态分区架构
# ----------------------------------------------------------------------------
yb = register("yb171_2026_native")
_mk(yb, "mass_u", 170.9363316, "u", "external_measured", "Yb-171 同位素质量")
_mk(yb, "wavelength_m", 487e-9, "m", "paper_measured", "Y：487 nm 光阱")
_mk(yb, "waist_m", 0.77e-6, "m", "paper_measured", "Y Methods：w=0.77 μm")
_mk(yb, "zr_m", math.pi * 0.77e-6**2 / 487e-9, "m", "paper_derived", "zR=π w²/λ")
# 存储阱深由 fr=30 kHz 反推（paper_derived）；门区调制平均深度为实测
_ybm = 170.9363316 * U_TO_KG
_ybw = 0.77e-6
_ybzr = math.pi * _ybw**2 / 487e-9
_U_yb = _ybm * (2 * math.pi * 30e3) ** 2 * _ybw**2 / 4.0
_mk(yb, "storage_depth_uK", _U_yb / KB / 1e-6, "uK", "paper_derived",
    "U=m·ω_r²·w²/4，ω_r=2π·30 kHz（fr=30 kHz）", note="≈109 μK")
_mk(yb, "gate_zone_avg_depth_uK", 37.0, "uK", "paper_measured",
    "Y ED Fig.5：50% 占空比 V_avg=37 μK")
_mk(yb, "omega_r_rad_s", 2 * math.pi * 30e3, "rad/s", "paper_measured",
    "Y ED Fig.5：fr=30 kHz")
_mk(yb, "move_time_s", 0.89e-3, "s", "paper_measured", "Y：单程移动 0.89 ms")
_mk(yb, "handoff_time_s", 0.78e-3, "s", "paper_measured",
    "Y：每次交接 ramp 0.78 ms（单程两次）")
_mk(yb, "zone_distance_m", 50e-6, "m", "assumed_sensitivity",
    "存储区-门区间距：论文文本未给出精确值，默认 50 μm，Y1 中作灵敏度扫描")
# 有效透镜增益：zs(t)=tau_eff·zR·ẋ(t)。映射 level3 的 zs=2σ·zR·ċ/v_s（τ_eff=2σ/v_s）；
# 仓库 Level 3 的有效场景 v_s≈0.54 m/s、σ=1 → τ_eff≈3.7 s/m，量级一致。
# 论文的 τ_AOD=w_AOD/v_s 具体数值在可获取文本中未给出；按“失焦在快速运输区间
# 起主导”取默认 4.0 s/m，Y1 中显式扫描 {0, 2, 4, 8}。此参数是 assumed，
# 不得用于绝对曲线验收，只用于结构判据。
_tau_eff = 4.0  # s/m：δz/zR = τ_eff·ẋ；ẋ=0.1 m/s 时 δz/zR=0.4
_mk(yb, "lensing_tau_eff_s_per_m", _tau_eff, "s/m", "assumed_sensitivity",
    "zs=τ_eff·zR·ẋ；对应 level3 的 2σ/v_s≈3.7；扫描 {0, 2, 4, 8}",
    note="默认值仅用于结构对照，不得标为论文标定")
_mk(yb, "eta_vls", 3e-4, "", "paper_measured", "Y：487 nm 矢量光移比")
_mk(yb, "synthetic_B_G", 0.05, "G", "paper_measured", "Y：合成磁场 50 mG")
_mk(yb, "splitting_parallel_hz", 57.0, "Hz", "paper_measured",
    "Y：偏振平行时 57 Hz 劈裂（50 mG）")
_mk(yb, "optimal_config_splitting_hz", 0.3, "Hz", "paper_measured",
    "Y：最优配置 B_SLM∥B_AOD⊥B₀")
_mk(yb, "trip_flip_dephase_prob", 6e-4, "", "paper_measured",
    "Y：最优配置每单程失相干+翻转 <6×10⁻⁴（上界）")
_mk(yb, "vls_phase_per_nm_rad", 7e-3, "rad/nm", "paper_measured",
    "Y：每单程相位灵敏度 7 mrad/nm")
_mk(yb, "alignment_sigma_m", 80e-9, "m", "paper_measured", "Y：SLM-AOD 对准 80 nm")
_mk(yb, "transverse_offset_m", 30e-12, "m", "paper_measured", "Y：δy≈30 pm")
_mk(yb, "transverse_sigma_m", 31e-9, "m", "paper_measured", "Y：σy=31 nm")
_mk(yb, "scatter_prob_per_handoff", 0.001, "", "paper_measured",
    "Y：交接散射 0.1%/次")
_mk(yb, "scatter_prob_per_move", 0.001, "", "paper_measured",
    "Y：移动散射 0.1%/次（存储深度）")
_mk(yb, "scattering_depth_exponent", 2.0, "", "paper_measured",
    "Y ED Fig.3：散射率对阱深二次拟合")
_mk(yb, "leakage_erasure_fraction", 0.72, "", "paper_measured",
    "Y：泄漏中可擦除占比 r_e,l≈0.72")
_mk(yb, "erasure_check_fp", 0.014, "", "paper_measured", "Y：擦除检查 FP 1.4%")
_mk(yb, "erasure_check_fn", 0.014, "", "paper_measured", "Y：擦除检查 FN 1.4%")
_mk(yb, "terminal_imaging_fp", 0.001, "", "paper_measured", "Y：终端成像 FP 0.1%")
_mk(yb, "terminal_imaging_fn", 0.005, "", "paper_measured", "Y：终端成像 FN 0.5%")
_mk(yb, "roundtrip1_loss_prob", 0.011, "", "paper_measured",
    "Y：首次往返 ³P₀ 损失 1.1(4)%", 0.004)
_mk(yb, "roundtrip1_erasure_frac", 0.5, "", "paper_measured",
    "Y：首次往返损失中约一半可擦除检测")
_mk(yb, "roundtrip5_erasure_frac", 0.10, "", "paper_measured",
    "Y：第 5 次往返擦除占比降至 ~10%")
_mk(yb, "gate_ar_error", 0.016, "", "paper_measured", "Y：ϵ_AR=0.016(1)", 0.001)
_mk(yb, "gate_ar_conditional_error", 0.010, "", "paper_measured",
    "Y：无擦除条件 ϵ_c=0.010(1)", 0.001)
_mk(yb, "gate_to_error", 0.011, "", "paper_measured", "Y：ϵ_TO=0.011(1)", 0.001)
_mk(yb, "gate_ar_erasure_frac", 0.38, "", "paper_measured",
    "Y：AR 误差 38(6)% 可擦除", 0.06)
_mk(yb, "gate_to_loss_frac", 0.75, "", "paper_measured", "Y：TO 误差 75% 为丢失")
_mk(yb, "gate_ar_scaling_exp", 4.0, "", "paper_measured",
    "Y ED Fig.6c：AR ϵ∝(ΔI/I)⁴")
_mk(yb, "gate_to_scaling_exp", 2.0, "", "paper_measured",
    "Y ED Fig.6c：TO ϵ∝(ΔI/I)²")
_mk(yb, "oneq_rf_fidelity", 0.9990, "", "paper_measured", "Y：F1Q=0.9990(1)")
_mk(yb, "rydberg_omega_hz", 2.5e6, "Hz", "paper_measured",
    "Y：302 nm 单光子 Ω_r/2π=2.5 MHz")
_mk(yb, "rydberg_lifetime_s", 88e-6, "s", "paper_measured",
    "Y ED Fig.6c：τ_Ryd=88 μs")
_mk(yb, "code_prep_fid_raw", 0.981, "", "paper_measured", "Y：0.981(2)", 0.002)
_mk(yb, "code_prep_fid_flag", 0.990, "", "paper_measured", "Y：flag 0.990(1)", 0.001)
_mk(yb, "code_prep_fid_erasure", 0.995, "", "paper_measured",
    "Y：擦除后选 0.995(1)", 0.001)
_mk(yb, "decode_unconditional", 0.874, "", "paper_measured", "Y：0.874(4)", 0.004)
_mk(yb, "decode_with_erasure", 0.902, "", "paper_measured", "Y：0.902(3)", 0.003)
_mk(yb, "teleport_fid_raw", 0.771, "", "paper_measured", "Y：0.771(9)", 0.009)
_mk(yb, "teleport_fid_erasure_ps", 0.87, "", "paper_measured", "Y：0.87(2)", 0.02)
_mk(yb, "teleport_fid_adaptive", 0.802, "", "paper_measured", "Y：0.802(8)", 0.008)
_mk(yb, "transport_error_share_logical", 0.66, "", "paper_measured",
    "Y：t=0 逻辑误差中运输占 66%（数值模型）")
_mk(yb, "moves_per_atom", 4.0, "", "paper_measured", "Y：每原子 4 次移动")
_mk(yb, "sequential_move_overhead_s", 29e-3, "s", "paper_measured",
    "Y：串行移动 29 ms/原子")
_mk(yb, "detection_latency_s", 25e-3, "s", "paper_measured",
    "Y：EMCCD→自适应控制 25 ms")
_mk(yb, "gamma_zero_jerk", 1.5625, "", "paper_measured",
    "Y：zero-jerk 轨迹 γ=1.5625")
_mk(yb, "gamma_min_jerk", 1.875, "", "paper_measured",
    "Y：minimum-jerk γ=1.875")
_mk(yb, "trap_mod_freq_hz", 400e3, "Hz", "paper_measured",
    "Y ED Fig.5：f_m=400 kHz，50% 占空")
_mk(yb, "trap_mod_cycles", 1e6, "", "paper_measured", "Y ED Fig.5：>10⁶ 周期")
_mk(yb, "trap_mod_adab_ramp_s", 100e-6, "s", "paper_measured",
    "Y：占空比绝热渐变 100 μs")


def species_summary() -> dict:
    from .provenance import REGISTRIES
    return {k: v.as_dict() for k, v in REGISTRIES.items()}
