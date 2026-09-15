"""修复导出工具丢失的原始图表工作簿；逐 slide/chart 对齐并原样复制，禁止重造数据。"""
import argparse
from pathlib import PurePosixPath
import posixpath
from zipfile import ZipFile, ZIP_DEFLATED
import xml.etree.ElementTree as ET

REL = "http://schemas.openxmlformats.org/package/2006/relationships"
OFFICE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
CHART = "http://schemas.openxmlformats.org/drawingml/2006/chart"
CT = "http://schemas.openxmlformats.org/package/2006/content-types"


def relpart(part):
    p = PurePosixPath(part)
    return str(p.parent / "_rels" / (p.name + ".rels"))


def target(part, rel):
    return posixpath.normpath(posixpath.join(posixpath.dirname(part), rel)).lstrip("/")


def charts(zipfile, slide):
    part = f"ppt/slides/slide{slide}.xml"
    relationships = ET.fromstring(zipfile.read(relpart(part)))
    refs = {x.get("Id"): x.get("Target") for x in relationships}
    return [target(part, refs[x.get(f"{{{OFFICE}}}id")])
            for x in ET.fromstring(zipfile.read(part)).iter(f"{{{CHART}}}chart")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source")
    parser.add_argument("candidate")
    parser.add_argument("output")
    args = parser.parse_args()
    with ZipFile(args.source) as original, ZipFile(args.candidate) as draft:
        payload = {n: draft.read(n) for n in draft.namelist()}
        copied = set()

        def restore(part):
            if part in copied:
                return
            copied.add(part)
            payload[part] = original.read(part)
            rp = relpart(part)
            if rp in original.namelist():
                payload[rp] = original.read(rp)
                for rel in ET.fromstring(payload[rp]):
                    if rel.get("TargetMode") != "External":
                        restore(target(part, rel.get("Target")))

        source_slides = [n for n in original.namelist()
                         if n.startswith("ppt/slides/slide") and n.endswith(".xml")]
        for i in range(1, len(source_slides) + 1):
            before, after = charts(original, i), charts(draft, i)
            if before != after:
                raise ValueError(f"Source/draft chart mapping changed on slide {i}")
            for part in before:
                restore(part)
        types = ET.fromstring(payload["[Content_Types].xml"])
        existing = {(x.tag, x.get("PartName"), x.get("Extension")) for x in types}
        for item in ET.fromstring(original.read("[Content_Types].xml")):
            needed = (item.get("PartName", "").lstrip("/") in copied
                      or item.get("Extension") == "xlsx")
            key = (item.tag, item.get("PartName"), item.get("Extension"))
            if needed and key not in existing:
                types.append(item)
                existing.add(key)
        ET.register_namespace("", CT)
        payload["[Content_Types].xml"] = ET.tostring(types, encoding="utf-8", xml_declaration=True)
        with ZipFile(args.output, "w", ZIP_DEFLATED) as result:
            for name, content in payload.items():
                result.writestr(name, content)
        print(f"Restored {sum(p.endswith('.xlsx') for p in copied)} original workbooks byte-for-byte")


if __name__ == "__main__":
    main()
