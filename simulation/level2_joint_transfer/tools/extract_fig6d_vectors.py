"""Read independent fit paths, marker centers and error bars from paper PDF.

No raster thresholding, smoothing, monotonic regression, or fit to simulation.
The archived raster reference remains unchanged. This creates versioned outputs.
"""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
import pdfplumber


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parents[1]
PDF = ROOT / "2403.12021v4.pdf"
OUT = PROJECT / "reference" / "fig6d_vector_v2"
COLORS = {
    "manual_200us": (0.886719,) * 3,
    "manual_400us": (0.740234,) * 3,
    "manual_600us": (0.533203,) * 3,
    "ml_400us": (0.364746, 0.662109, 0.781250),
}


def color_equal(a, b):
    return isinstance(a, (tuple, list)) and len(a) == 3 and np.max(np.abs(np.array(a) - b)) < 1e-5


def write_csv(file, header, rows):
    with file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def extract():
    with pdfplumber.open(PDF) as doc:
        page = doc.pages[6]
        # These screen-coordinate bounds identify only Fig. 6d's survival panel.
        # Exact calibration is recovered from vector tick marks, not these bounds.
        black = [l for l in page.lines if color_equal(l["stroking_color"], (0, 0, 0))]
        xticks = sorted(l["x0"] for l in black if 420 < l["x0"] < 520
                        and abs(l["x1"] - l["x0"]) < 1e-8
                        and 111 < l["top"] < 113 and 113 < l["bottom"] < 115)
        yticks = sorted(l["top"] for l in black if 421 < l["x0"] < 423
                        and 1 < l["x1"] - l["x0"] < 3
                        and abs(l["bottom"] - l["top"]) < 1e-8
                        and 60 < l["top"] < 110)
        assert len(xticks) == 4 and len(yticks) == 4
        xcoef = np.polyfit(xticks, [0, 20, 40, 60], 1)
        ycoef = np.polyfit(yticks, [1, .8, .6, .4], 1)
        nx = lambda x: np.polyval(xcoef, x)
        sy = lambda y: np.polyval(ycoef, y)
        series = {}
        for name, color in COLORS.items():
            curves = [(i, c) for i, c in enumerate(page.curves)
                      if color_equal(c["stroking_color"], color)
                      and c["x0"] > 420 and c["x1"] < 522 and c["top"] < 70]
            fits = [(i, c) for i, c in curves if c["stroke"] and not c["fill"]
                    and c["x1"] - c["x0"] > 60]
            assert len(fits) == 1, (name, len(fits))
            fit_index, fit = fits[0]
            assert all(op[0] in ("m", "l") for op in fit["path"])
            points = np.array(fit["pts"])
            nfit, sfit = nx(points[:, 0]), sy(points[:, 1])
            assert np.all(np.diff(nfit) > 0) and np.all(np.diff(sfit) <= 0)
            markers = []
            for idx, c in enumerate(page.curves):
                if not (c["fill"] and c["stroke"] and color_equal(c["stroking_color"], color)
                        and 3.5 < c["x1"] - c["x0"] < 4.8
                        and 420 < c["x0"] < 517 and 60 < c["top"] < 114):
                    continue
                cx, cy = (c["x0"] + c["x1"]) / 2, (c["top"] + c["bottom"]) / 2
                n = float(nx(cx))
                if not (0 <= n <= 60.01 and abs(n - round(n)) < .005):
                    continue  # Reject legend diamonds by their off-grid x.
                error_lines = [(j, l) for j, l in enumerate(page.lines)
                               if color_equal(l["stroking_color"], color)
                               and abs(l["x1"] - l["x0"]) < 1e-8
                               and abs(l["x0"] - cx) < .005
                               and abs((l["top"] + l["bottom"]) / 2 - cy) < .006]
                assert len(error_lines) == 1, (name, n, len(error_lines))
                error_idx, err = error_lines[0]
                markers.append({"n": int(round(n)), "n_unrounded": n, "survival": float(sy(cy)),
                                "lower": float(sy(err["bottom"])), "upper": float(sy(err["top"])),
                                "pdf_center_pt": [cx, cy], "curve_index": idx, "error_line_index": error_idx})
            markers.sort(key=lambda m: m["n"])
            expected = [2, 4, 8, 12, 16, 20, 30] + ([] if name == "manual_200us" else [40, 50, 60])
            assert [m["n"] for m in markers] == expected, (name, markers)
            max_n = 31 if name == "manual_200us" else 60
            # Preserve the prior comparison interval to isolate reference changes.
            grid = np.arange(5, max_n + 1)
            reference = np.interp(grid, nfit, sfit)
            # Use original PDF polyline vertices (plus cropped endpoints) for display.
            keep = (nfit >= 0) & (nfit <= max_n)
            display_n = np.r_[nfit[keep], max_n]
            display_s = np.interp(display_n, nfit, sfit)
            series[name] = {"fit_path_curve_index": fit_index, "color": color,
                            "fit_path_pdf_pt": points.tolist(), "fit_path_n": nfit.tolist(),
                            "fit_path_survival": sfit.tolist(), "markers": markers,
                            "comparison_n": grid.tolist(), "comparison_fit": reference.tolist(),
                            "display_n": display_n.tolist(), "display_fit": display_s.tolist()}
        result = {"schema": "fig6d_vector_v2", "source_pdf": str(PDF.relative_to(ROOT)),
                  "source_sha256": hashlib.sha256(PDF.read_bytes()).hexdigest(), "page": 7,
                  "method": "Separate stroked fit polylines, filled diamond centers, and vertical error-bar paths; affine vector tick calibration",
                  "not_performed": ["raster tracing", "smoothing", "refitting", "monotonic correction", "fitting to Monte Carlo"],
                  "limitations": "Coordinates describe the published figure, not authors' raw experimental files. Error-bar endpoints are recovered geometrically; do not assign a confidence convention not established by the paper. Hidden vector objects remain recoverable. Dense fit samples are not independent measurements.",
                  "calibration": {"x_ticks_pdf_pt": xticks, "x_values": [0, 20, 40, 60],
                                  "y_ticks_pdf_pt": yticks, "y_values": [1, .8, .6, .4],
                                  "n_from_pdf_x": xcoef.tolist(), "survival_from_pdf_top": ycoef.tolist(),
                                  "max_x_tick_residual_n": float(max(abs(nx(xticks) - [0, 20, 40, 60]))),
                                  "max_y_tick_residual": float(max(abs(sy(yticks) - [1, .8, .6, .4])))},
                  "series": series}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "fig6d_vectors.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csv(OUT / "experimental_marker_centers.csv", ["waveform", "n", "survival", "error_lower", "error_upper", "pdf_curve_index"],
              [(k, m["n"], m["survival"], m["lower"], m["upper"], m["curve_index"]) for k, s in series.items() for m in s["markers"]])
    write_csv(OUT / "published_fit_comparison_grid.csv", ["waveform", "n", "survival_fit"],
              [(k, n, v) for k, s in series.items() for n, v in zip(s["comparison_n"], s["comparison_fit"])])
    names = list(COLORS)
    lookups = {k: dict(zip(s["comparison_n"], s["comparison_fit"])) for k, s in series.items()}
    write_csv(OUT / "published_fit_reference.csv", ["n_one_way_transfers"] + ["survival_" + k for k in names],
              [(n, *[lookups[k].get(n, "") for k in names]) for n in range(5, 61)])
    print(json.dumps({"calibration": result["calibration"], "series": {k: {"markers": len(s["markers"]), "fit_endpoint": s["comparison_fit"][-1]} for k, s in series.items()}}, ensure_ascii=False))
    return result


if __name__ == "__main__":
    extract()
