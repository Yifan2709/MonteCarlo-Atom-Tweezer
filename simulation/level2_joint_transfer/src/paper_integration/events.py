"""M1 事件账本与随机流组织。

- 真实事件（true event）与检测记录（detection record）分开保存；
  策略/解码器只能读取 detection 记录（任务书 §7）。
- 随机流按 (trial_id, atom_id, channel) 键控，删除死原子或数组重排
  不改变其他原子的噪声实现（§6.4）。
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

import numpy as np

TRUE_KINDS = ("mechanical_loss", "leakage", "photoion", "pauli_error",
              "gate_error", "idle_error")
DETECTION_KINDS = ("erasure_check", "terminal_image", "stabilizer_readout")


def stream_seed(base_seed: int, *keys) -> int:
    """稳定键控种子：同一键序列始终得到同一流。"""
    payload = json.dumps([int(base_seed), *[str(k) for k in keys]],
                         sort_keys=True).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def channel_rng(base_seed: int, trial_id: int, atom_id: int, channel: str):
    return np.random.default_rng(stream_seed(base_seed, trial_id, atom_id,
                                             channel))


@dataclass
class TrueEvent:
    trial_id: int
    atom_id: int
    t_s: float
    kind: str
    channel: str
    detail: dict = field(default_factory=dict)


@dataclass
class DetectionRecord:
    trial_id: int
    atom_id: int
    t_s: float
    kind: str          # erasure_check / terminal_image / ...
    detected: bool | None
    meta: dict = field(default_factory=dict)


class Ledger:
    """一次运行的事件账本。"""

    def __init__(self):
        self.true_events: list[TrueEvent] = []
        self.detections: list[DetectionRecord] = []

    def add_true(self, ev: TrueEvent) -> None:
        if ev.kind not in TRUE_KINDS:
            raise ValueError(f"true event kind {ev.kind!r} 不在登记表")
        self.true_events.append(ev)

    def add_detection(self, rec: DetectionRecord) -> None:
        if rec.kind not in DETECTION_KINDS:
            raise ValueError(f"detection kind {rec.kind!r} 不在登记表")
        self.detections.append(rec)

    def detections_for(self, trial_id: int, upto_t_s: float | None = None):
        out = [d for d in self.detections if d.trial_id == trial_id
               and (upto_t_s is None or d.t_s <= upto_t_s)]
        return out

    def summary(self) -> dict:
        by_kind: dict[str, int] = {}
        for e in self.true_events:
            by_kind[e.kind] = by_kind.get(e.kind, 0) + 1
        det = {}
        for d in self.detections:
            det[d.kind] = det.get(d.kind, 0) + 1
        return {"true_events": by_kind, "detections": det}


def erasure_check_detection(trial_id, atom_id, t_s, has_leakage: bool,
                            fp: float, fn: float, rng) -> DetectionRecord:
    """擦除检查：泄漏→以 1−FN 检出；无泄漏→以 FP 误报。仅输出观测。"""
    if has_leakage:
        detected = bool(rng.random() >= fn)
    else:
        detected = bool(rng.random() < fp)
    return DetectionRecord(trial_id, atom_id, t_s, "erasure_check", detected)
