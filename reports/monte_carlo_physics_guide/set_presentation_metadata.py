"""Set the report title and leave author fields empty before final validation."""
from pathlib import Path
from zipfile import ZipFile
import sys
import xml.etree.ElementTree as ET

source, destination = map(Path, sys.argv[1:3])
namespaces = {
    '': 'http://schemas.openxmlformats.org/package/2006/metadata/core-properties',
    'dc': 'http://purl.org/dc/elements/1.1/',
    'dcterms': 'http://purl.org/dc/terms/',
    'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
}
for prefix, uri in namespaces.items():
    ET.register_namespace(prefix, uri)
with ZipFile(source) as zin, ZipFile(destination, 'w') as zout:
    for info in zin.infolist():
        payload = zin.read(info.filename)
        if info.filename == 'docProps/core.xml':
            root = ET.fromstring(payload)
            for uri, tag, value in [
                (namespaces['dc'], 'title', '光镊原子转移的 Monte Carlo 模拟：方法与结果汇总'),
                (namespaces['dc'], 'creator', ''),
                (namespaces[''], 'lastModifiedBy', ''),
            ]:
                element = root.find(f'{{{uri}}}{tag}')
                if element is None:
                    element = ET.SubElement(root, f'{{{uri}}}{tag}')
                element.text = value
            payload = ET.tostring(root, encoding='utf-8', xml_declaration=True)
        zout.writestr(info, payload)
