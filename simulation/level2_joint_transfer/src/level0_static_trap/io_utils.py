"""安全输出 CSV/JSON/YAML/Markdown 及 PNG 程序化验收。"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from PIL import Image
import yaml


def ensure_output_directory(path: str | Path, project_root: Path) -> Path:
    """创建输出目录，并拒绝解析后落入项目 sources/ 的路径。"""

    target = Path(path).expanduser().resolve()
    sources = (project_root / "sources").resolve()
    if target == sources or sources in target.parents:
        raise ValueError(f"拒绝向只读参考目录 sources/ 写入输出: {target}")
    target.mkdir(parents=True, exist_ok=True)
    return target


def to_builtin(value: Any) -> Any:
    """递归转换 NumPy/Path 类型，并拒绝非有限浮点数。"""

    if isinstance(value, Mapping):
        return {str(key): to_builtin(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_builtin(item) for item in value]
    if isinstance(value, np.ndarray):
        return [to_builtin(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return to_builtin(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("输出包含 NaN 或 Infinity")
    return value


def save_json(path: Path, data: Any) -> None:
    """以严格 JSON 保存，不允许 NaN 或 Infinity。"""

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(to_builtin(data), handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


def save_yaml(path: Path, data: Any) -> None:
    """以 UTF-8 YAML 保存配置。"""

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        yaml.safe_dump(to_builtin(data), handle, allow_unicode=True, sort_keys=False)


def save_markdown(path: Path, text: str) -> None:
    """保存 UTF-8 Markdown。"""

    path.write_text(text, encoding="utf-8", newline="\n")


def save_csv(path: Path, columns: Mapping[str, Any]) -> None:
    """按列保存 CSV，并检查长度一致及数据有限。"""

    names = list(columns)
    arrays = [np.asarray(columns[name]) for name in names]
    if not arrays or len({array.size for array in arrays}) != 1:
        raise ValueError("CSV 各列长度必须一致且至少包含一列")
    for name, array in zip(names, arrays):
        if np.issubdtype(array.dtype, np.number) and not np.all(np.isfinite(array)):
            raise ValueError(f"CSV 列 {name} 含 NaN 或 Infinity")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        writer.writerows(zip(*(array.tolist() for array in arrays)))


def validate_pngs(plot_records: list[dict]) -> dict[str, Any]:
    """重新打开每张 PNG，检查尺寸、方差、数据、坐标、图例和关闭状态。"""

    images: dict[str, Any] = {}
    for record in plot_records:
        path = Path(record["path"])
        checks: dict[str, Any] = {
            "exists": path.is_file(),
            "file_size_over_5k": path.is_file() and path.stat().st_size > 5000,
            "reopen_success": False,
            "dimensions_valid": False,
            "pixel_variance": None,
            "not_nearly_monochrome": False,
            "data_finite": record["data_finite"],
            "data_x_range_nonzero": record["data_x_range_nonzero"],
            "data_y_range_nonzero": record["data_y_range_nonzero"],
            "axes_cover_data": record["axes_cover_data"],
            "legend_matches_series": record["legend_matches_series"],
            "figure_closed": record["figure_closed"],
        }
        try:
            with Image.open(path) as image:
                image.load()
                pixels = np.asarray(image.convert("RGB"), dtype=float)
                variance = float(np.var(pixels))
                checks["reopen_success"] = True
                checks["width_px"], checks["height_px"] = image.size
                checks["dimensions_valid"] = image.size[0] >= 1000 and image.size[1] >= 600
                checks["pixel_variance"] = variance
                checks["not_nearly_monochrome"] = variance > 1.0
        except Exception as exc:  # pragma: no cover - failure is serialized for diagnosis
            checks["reopen_error"] = str(exc)
        boolean_checks = [value for value in checks.values() if isinstance(value, bool)]
        checks["passed"] = all(boolean_checks)
        images[path.name] = checks
    return {
        "all_passed": all(item["passed"] for item in images.values()),
        "visual_review_performed": False,
        "visual_review_findings": [],
        "programmatic_checks_scope": (
            "文件有效性、像素非空、数据有限性、非零范围、坐标覆盖、图例数量与 figure 关闭；"
            "不证明文字无重叠、曲线易读或布局美观。"
        ),
        "images": images,
    }

