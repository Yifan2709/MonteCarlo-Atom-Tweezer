"""Freeze source provenance and hash all deliverables without deleting anything."""
from pathlib import Path
import csv,json,hashlib,subprocess,zipfile,datetime,importlib.metadata
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parent.parent

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda:f.read(2**20),b''):h.update(data)
    return h.hexdigest()

def main():
    code_files=list((REPO/'simulation/level2_joint_transfer/src').rglob('*.py'))
    code_files+=list((REPO/'simulation/level2_joint_transfer/tests').rglob('*.py'))
    code_files+=list(ROOT.glob('*.py'))+list(ROOT.glob('*.ps1'))
    code_files+=[ROOT/'requirements_exact.txt',REPO/'simulation/level2_joint_transfer/pyproject.toml']
    hashes={str(p.relative_to(REPO)).replace('\\','/'):sha(p) for p in sorted(code_files)}
    with zipfile.ZipFile(ROOT/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(code_files):
            info=zipfile.ZipInfo(str(p.relative_to(REPO)).replace('\\','/'),(1980,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED;z.writestr(info,p.read_bytes())
    git=lambda *args:subprocess.check_output(['git',*args],cwd=REPO,text=True).strip()
    paper_dir=REPO/'tmp/pdfs/paper_relation_20260915'
    meta=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_commit=git('rev-parse','HEAD'),branch=git('branch','--show-current'),
        physics_commit='63dc14e29c41e229d88cc97d8a33549fe8ada4ef',
        baseline_remote_commit='1fed6329e95c4439820aff31d644376b817b5d6f',
        source_files_sha256=hashes,source_snapshot_sha256=sha(ROOT/'source_snapshot.zip'),
        papers_sha256={p.name:sha(p) for p in [paper_dir/'bluvstein_2022.pdf',paper_dir/'zhang_2026.pdf',REPO/'2403.12021v4.pdf']},
        dependency_versions={n:importlib.metadata.version(n) for n in ['numpy','scipy','numba','llvmlite','matplotlib','pandas','pdfplumber','pypdf','pypdfium2','pytest','PyYAML','Pillow','setuptools','wheel']},
        raw_data_policy='All runs/*.npz kept locally and hashed below, excluded only from Git; copy them with the study for full archival transfer.',
        experimental_validation='Incomplete: quantitative paper residuals retained; Cs apparatus parameters uncalibrated.',
        author_data=dict(doi='10.5281/zenodo.19491381',license='CC-BY-4.0',sha256=sha(ROOT/'zhang_dataset.zip')))
    (ROOT/'provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    rows=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or p.name=='artifact_manifest.csv' or '__pycache__' in p.parts:continue
        rows.append(dict(path=str(p.relative_to(ROOT)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)))
    with (ROOT/'artifact_manifest.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=['path','bytes','sha256']);writer.writeheader();writer.writerows(rows)
    print('Hashed',len(rows),'files;',sum(r['bytes'] for r in rows),'bytes. Source:',meta['source_commit'])

if __name__=='__main__':main()
