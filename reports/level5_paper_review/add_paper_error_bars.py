"""Add recovered asymmetric paper error bars to native scatter series.

Artifact Tool authors all slide/chart content. This OOXML step supplies custom
per-point errors unavailable in its public series error-bar configuration.
"""
from pathlib import Path
import json
import sys
import zipfile
import xml.etree.ElementTree as ET


C = "http://schemas.openxmlformats.org/drawingml/2006/chart"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
ET.register_namespace("c", C)
ET.register_namespace("a", A)
NS = {"c": C, "a": A}
keys = ["manual_200us", "manual_400us", "manual_600us", "ml_400us"]


def child(parent, tag, **attrs):
    return ET.SubElement(parent, f"{{{C}}}{tag}", attrs)


def main(source, dest, vector_file):
    data = json.loads(Path(vector_file).read_text("utf-8"))["series"]
    count = 0
    with zipfile.ZipFile(source) as zin, zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            raw = zin.read(item.filename)
            if "/charts/" in item.filename and item.filename.endswith(".xml"):
                root = ET.fromstring(raw)
                for series in root.findall(".//c:scatterChart/c:ser", NS):
                    values = series.findall("c:tx//c:v", NS)
                    if not any(v.text == "论文散点" for v in values):
                        continue
                    assert count < 4
                    markers = data[keys[count]]["markers"]
                    errors = ET.Element(f"{{{C}}}errBars")
                    child(errors, "errDir", val="y")
                    child(errors, "errBarType", val="both")
                    child(errors, "errValType", val="cust")
                    child(errors, "noEndCap", val="0")
                    for sign in ("plus", "minus"):
                        lit = child(child(errors, sign), "numLit")
                        child(lit, "formatCode").text = "0.00000000"
                        child(lit, "ptCount", val=str(len(markers)))
                        for i, m in enumerate(markers):
                            v = m["upper"] - m["survival"] if sign == "plus" else m["survival"] - m["lower"]
                            child(child(lit, "pt", idx=str(i)), "v").text = f"{v:.12g}"
                    sp = child(errors, "spPr")
                    ln = ET.SubElement(sp, f"{{{A}}}ln", {"w": "9525"})
                    fill = ET.SubElement(ln, f"{{{A}}}solidFill")
                    ET.SubElement(fill, f"{{{A}}}srgbClr", {"val": "17324D"})
                    series.insert(list(series).index(series.find("c:xVal", NS)), errors)
                    count += 1
                raw = ET.tostring(root, encoding="utf-8", xml_declaration=True)
            zout.writestr(item, raw)
    assert count == 4, count
    print(f"Added custom paper error bars to {count} scatter series")


if __name__ == "__main__":
    main(*sys.argv[1:])
