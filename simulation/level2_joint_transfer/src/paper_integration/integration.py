"""I1：动力学 → 量子线路的实际接通与参数扰动传播。

链路（真实数据依赖，非硬编码）：
- Y 侧：yb_transport 扫描 → 平均加热率/末态能量 → yb_leakage 机械/泄漏
  通道 → yb_code 存储/传送逻辑观测量。扰动 τ_eff ×2 → 重算链路 →
  逻辑观测量必须移动。
- R 侧：rb_transport 存活率–时长曲线 → 运输层 loss 率 → rb_circuit
  Steane 指标。扰动运输速度（T 减半）→ 指标必须移动。
"""
from __future__ import annotations

import numpy as np

from .provenance import REGISTRIES
from . import yb_transport, yb_leakage, yb_code, rb_transport, rb_circuit


def yb_chain(tau_eff: float, move_time_s: float, distance_m: float,
             shots: int = 256, seed: int = 4242, dt_s: float = 0.1e-6,
             trials_circuit: int = 20000) -> dict:
    """τ_eff → 运输 → 通道 → 逻辑。输出各级数值供扰动对比。

    耦合声明（assumed，结构模型）：首往返瞬态机械份额按
    1+4×(1−存活率) 随运输损失放大；加热率随末态能量与损失增长。
    """
    move = yb_transport.run_single_move(distance_m, move_time_s,
                                        "sixth_zero_jerk", tau_eff, shots,
                                        seed, dt_s)
    mean_energy_uK = abs(move["checkpoints"][-1]["mean_energy_bottom_uK"])
    heating_per_trip = 0.5 * mean_energy_uK * (1.0 - move["survival"] + 0.1)
    rt = yb_leakage.simulate_roundtrips(
        4096, 5, seed, heating_per_trip_uK=heating_per_trip,
        rt1_extra_mech_scale=1.0 + 4.0 * (1.0 - move["survival"]))
    per_atom_loss = 1 - (1 - rt["first_roundtrip_loss_frac"]) ** 0.25
    storage = yb_code.storage_experiment(
        trials_circuit, seed, 0.5,
        transport_loss=4 * per_atom_loss * 0.25 + 0.01)
    tele = yb_code.teleportation(trials_circuit, seed + 1,
                                 transport_loss_per_move=per_atom_loss)
    return {"tau_eff": tau_eff, "move_survival": move["survival"],
            "mean_energy_uK": mean_energy_uK,
            "heating_per_trip_uK": heating_per_trip,
            "rt1_loss_frac": rt["first_roundtrip_loss_frac"],
            "storage_unconditional": storage["unconditional"],
            "storage_erasure_informed": storage["erasure_informed"],
            "teleport_success": tele["success"]}


def yb_perturbation(move_time_s: float = 0.6e-3, distance_m: float = 100e-6,
                    shots: int = 256, seed: int = 4242) -> dict:
    """τ_eff = 4 vs 8：全链路重算，验证扰动传播到最终观测量。"""
    base = yb_chain(4.0, move_time_s, distance_m, shots, seed)
    pert = yb_chain(8.0, move_time_s, distance_m, shots, seed)
    moved = {
        k: bool(abs(pert[k] - base[k]) > 1e-9)
        for k in ("move_survival", "rt1_loss_frac", "teleport_success")}
    return {"baseline": base, "perturbed": pert, "propagated": moved,
            "all_propagated": all(moved.values())}


def rb_chain(move_duration_s: float, distance_m: float = 110e-6,
             shots: int = 256, seed: int = 777,
             trials_circuit: int = 20000) -> dict:
    res = rb_transport.run_single_transport(distance_m, move_duration_s,
                                            shots, seed, 0.05e-6)
    transport_layer_loss = max(0.0, 1.0 - res["survival"])
    steane = rb_circuit.steane_circuit(trials_circuit, seed + 5,
                                       transport_layer_loss=transport_layer_loss)
    return {"move_duration_s": move_duration_s,
            "move_survival": res["survival"],
            "transport_layer_loss": transport_layer_loss,
            "steane_raw_xbar": steane["raw_xbar"],
            "steane_no_error_frac": steane["no_detected_error_frac"]}


def rb_perturbation(distance_m: float = 110e-6, shots: int = 256,
                    seed: int = 777) -> dict:
    """T=300 μs vs 150 μm/μs 的快速档：验证运输损失传播进 Steane 指标。"""
    base = rb_chain(300e-6, distance_m, shots, seed)
    pert = rb_chain(150e-6, distance_m, shots, seed)
    moved = {
        k: bool(abs(pert[k] - base[k]) > 1e-9)
        for k in ("move_survival", "transport_layer_loss", "steane_raw_xbar")}
    return {"baseline": base, "perturbed": pert, "propagated": moved,
            "all_propagated": all(moved.values())}
