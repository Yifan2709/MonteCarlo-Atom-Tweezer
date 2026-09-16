"""参数溯源与登记：每个影响结论的参数必须声明来源类别。

来源类别（任务书 §5）：
- paper_measured      论文直接测量/引用的数值（含正文、Methods、ED 表格）
- paper_derived       由论文数值推导（保存推导式与误差传播）
- external_measured   论文之外的公开测量
- calibrated_on_training  在声明训练集上拟合（冻结后不得再改）
- assumed_sensitivity 假设值/灵敏度参数（不得标为论文实测）

单位约定：内部一律 SI（m, s, J, kg, rad, Hz 以 fields 后缀标注）；
对外便捷单位仅做展示换算，不进入计算核心。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class Provenance(str, Enum):
    PAPER_MEASURED = "paper_measured"
    PAPER_DERIVED = "paper_derived"
    EXTERNAL_MEASURED = "external_measured"
    CALIBRATED = "calibrated_on_training"
    ASSUMED = "assumed_sensitivity"


@dataclass(frozen=True)
class ParamSpec:
    name: str
    value: float
    unit: str
    provenance: Provenance
    source: str          # 论文位置或推导式，或假设理由
    uncertainty: float | None = None
    note: str = ""

    def as_dict(self) -> dict:
        return {
            "name": self.name, "value": self.value, "unit": self.unit,
            "provenance": self.provenance.value, "source": self.source,
            "uncertainty": self.uncertainty, "note": self.note,
        }


@dataclass
class ParamRegistry:
    """按 (config, name) 登记参数；合同与报告直接序列化本表。"""
    config: str
    params: dict[str, ParamSpec] = field(default_factory=dict)

    def add(self, spec: ParamSpec) -> ParamSpec:
        if spec.name in self.params:
            raise ValueError(f"duplicate param {spec.name} in {self.config}")
        self.params[spec.name] = spec
        return spec

    def get(self, name: str) -> ParamSpec:
        if name not in self.params:
            raise KeyError(
                f"param {name!r} not defined for config {self.config!r}; "
                "禁止静默回退到其他物种默认值")
        return self.params[name]

    def require_provenance(self, *names: str) -> None:
        for n in names:
            self.get(n)

    def as_dict(self) -> dict:
        return {"config": self.config,
                "params": {k: v.as_dict() for k, v in self.params.items()}}


# 全局登记表：验证脚本按合同核对这里的内容
REGISTRIES: dict[str, ParamRegistry] = {}


def register(config: str) -> ParamRegistry:
    reg = ParamRegistry(config)
    REGISTRIES[config] = reg
    return reg
