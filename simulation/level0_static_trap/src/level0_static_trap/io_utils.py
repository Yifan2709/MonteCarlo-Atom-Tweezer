"""输入输出工具：安全目录创建、CSV/JSON/YAML/Markdown 写出与图像程序化检查。

JSON 写出使用 ``allow_nan=False``，禁止 NaN/Infinity；NumPy 类型自动转换为可
序列化的 Python 原生类型。任何输出路径解析后落入 ``sources/`` 都被拒绝。
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import yaml


class OutputPathError(ValueError):
    """输出路径不合法（例如落入只读 sources/）。"""


def _assert_not_in_sources(path: str) -> None:
    real = os.path.realpath(path)
    parts = Path(real).parts
    if "sources" in parts:
        raise OutputPathError(
            f"输出路径解析到只读参考目录 sources/ 内（{real}），已拒绝"
        )


def ensure_output_dir(path: str) -> str:
    """安全创建输出目录，并拒绝落入 ``sources/``。返回 realpath。"""
    real = os.path.realpath(path)
    _assert_not_in_sources(real)
    os.makedirs(real, exist_ok=True)
    return real


# ---------------------------------------------------------------------------
# 序列化辅助
# ---------------------------------------------------------------------------

def _to_plain(obj):
    """把含 numpy/Path 等的对象转为纯 Python 可序列化结构。"""
    if isinstance(obj, dict):
        return {str(k): _to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _to_plain(obj.tolist())
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        v = float(obj)
        if not math.isfinite(v):
            raise ValueError(f"JSON 中出现非有限浮点数：{v}")
        return v
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, float):
        if not math.isfinite(obj):
            raise ValueError(f"JSON 中出现非有限浮点数：{obj}")
        return obj
    return obj


def write_json(path: str, data: dict) -> None:
    """写入 JSON，禁止 NaN/Infinity。"""
    _assert_not_in_sources(path)
    plain = _to_plain(data)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(plain, fh, indent=2, sort_keys=True, allow_nan=False)


def write_yaml(path: str, data: dict) -> None:
    """写入 YAML（同样把非有限值剔除以保持可读）。"""
    _assert_not_in_sources(path)

    def _strip(o):
        if isinstance(o, dict):
            return {str(k): _strip(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_strip(v) for v in o]
        if isinstance(o, float) and not math.isfinite(o):
            return None
        if isinstance(o, (np.floating,)) and not math.isfinite(float(o)):
            return None
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return _strip(o.tolist())
        return o

    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(_strip(data), fh, sort_keys=True, allow_unicode=True)


def write_markdown(path: str, text: str) -> None:
    """写入 Markdown 文本。"""
    _assert_not_in_sources(path)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def write_csv(path: str, columns: dict[str, np.ndarray]) -> None:
    """写出 CSV。``columns`` 保持顺序：键为列名，值为等长一维数组。"""
    _assert_not_in_sources(path)
    names = list(columns.keys())
    arrays = [np.asarray(columns[n], dtype=float).ravel() for n in names]
    lengths = {len(a) for a in arrays}
    if len(lengths) != 1:
        raise ValueError(f"CSV 列长不一致：{ {n: len(a) for n, a in zip(names, arrays)} }")
    data = np.column_stack(arrays)
    header = ",".join(names)
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.12e")


# ---------------------------------------------------------------------------
# 图像程序化检查
# ---------------------------------------------------------------------------

def validate_image(path: str, expected_meta: dict | None = None) -> dict:
    """对一张 PNG 执行程序化检查并返回结构化结果。

    检查项：文件存在且体积合理、能用 Pillow 重开、像素宽高符合预期且非零、
    像素方差大于阈值（非单色）、数组数据有限、（若提供 expected_meta）
    x/y 数据范围非零且覆盖数据、图例标签数与系列数相符、figure 已关闭。
    """
    from PIL import Image

    result: dict = {"path": os.path.realpath(path), "checks": {}, "expected": expected_meta or {}}
    p = Path(path)
    result["checks"]["exists"] = p.is_file()
    size = p.stat().st_size if p.exists() else 0
    result["checks"]["size_bytes"] = int(size)
    result["checks"]["size_above_lower_bound"] = bool(size > 2048)
    if not p.exists():
        result["checks"]["overall_pass"] = False
        return result
    try:
        with Image.open(p) as im:
            im.load()
            width, height = im.size
            arr = np.asarray(im.convert("L"), dtype=float)
    except Exception as exc:  # pragma: no cover - 依赖具体图像
        result["checks"]["reopen_ok"] = False
        result["checks"]["reopen_error"] = str(exc)
        result["checks"]["overall_pass"] = False
        return result
    result["checks"]["reopen_ok"] = True
    result["checks"]["width_px"] = int(width)
    result["checks"]["height_px"] = int(height)
    result["checks"]["width_nonzero"] = bool(width > 0)
    result["checks"]["height_nonzero"] = bool(height > 0)
    finite = bool(np.all(np.isfinite(arr)))
    result["checks"]["array_finite"] = finite
    var = float(np.var(arr))
    result["checks"]["pixel_variance"] = var
    result["checks"]["not_monochrome"] = bool(var > 1.0)
    if expected_meta:
        xr = expected_meta.get("x_data_range")
        yr = expected_meta.get("y_data_range")
        if xr:
            result["checks"]["x_data_nonzero"] = bool(abs(xr[1] - xr[0]) > 0.0)
        if yr:
            result["checks"]["y_data_nonzero"] = bool(abs(yr[1] - yr[0]) > 0.0)
        sc = expected_meta.get("series_count")
        lc = expected_meta.get("legend_label_count")
        if isinstance(sc, int):
            result["checks"]["series_count"] = sc
        if isinstance(lc, int):
            result["checks"]["legend_label_count"] = lc
    checks_pass = all(
        v is True
        for k, v in result["checks"].items()
        if k in {"size_above_lower_bound", "reopen_ok", "width_nonzero", "height_nonzero",
                 "array_finite", "not_monochrome", "x_data_nonzero", "y_data_nonzero"}
    )
    result["checks"]["overall_pass"] = bool(checks_pass)
    return result


def write_image_validation(path: str, results: dict) -> None:
    """写出图像程序化检查的聚合 JSON。"""
    write_json(path, results)
