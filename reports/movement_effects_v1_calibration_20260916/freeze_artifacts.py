"""Freeze source, references, environment, and an inventory of all computed atoms."""
from run_study import *
from datetime import datetime,timezone
import subprocess,zipfile,importlib.metadata

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','numba','matplotlib','PyYAML','pdfplumber','pytest']}
    (HERE/'requirements_exact.txt').write_text('\n'.join(f'{k}=={v}' for k,v in versions.items())+'\n',encoding='utf-8')
    files=[]
    for directory in ['src','configs','reference','tests']:
        files.extend(p for p in (ROOT/directory).rglob('*') if p.is_file()
                     and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.nbc','.nbi'))
    files.extend([ROOT/'pyproject.toml',ROOT/'README.md',ROOT/SETTINGS['ml_waveform']])
    files.extend(p for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py','.ps1','.yaml','.md') )
    files.extend(p for p in HERE.glob('*queue.json'))
    files.extend([HERE/'ml_vector_nodes.json',HERE/'final_lock.json',HERE/'base_commit.txt',
                  HERE/'requirements_exact.txt',HERE/'extra_timestep_config.json'])
    files=sorted(set(files))
    with zipfile.ZipFile(HERE/'source_final.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,str(p.relative_to(WORK)))
    refs=[WORK/'2403.12021v4.pdf',ROOT/'reference/fig6d_vector_v2/experimental_marker_centers.csv',
          ROOT/'reference/fig6d_vector_v2/published_fit_reference.csv',ROOT/SETTINGS['ml_waveform'],
          HERE/'ml_vector_nodes.json']
    inventory=[]
    for p in sorted([*(HERE/'runs').rglob('*'),*(HERE/'timestep_extra').rglob('*')]):
        if p.is_file():inventory.append(dict(path=str(p.relative_to(HERE)),bytes=p.stat().st_size,sha256=sha(p)))
    with (HERE/'raw_inventory.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256']);w.writeheader();w.writerows(inventory)
    manifest=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=WORK,text=True).strip(),
        branch=subprocess.check_output(['git','branch','--show-current'],cwd=WORK,text=True).strip(),
        audited_remote_commit=(HERE/'base_commit.txt').read_text().strip(),python=platform.python_version(),
        versions=versions,source_archive=dict(path='source_final.zip',sha256=sha(HERE/'source_final.zip')),
        source_files={str(p.relative_to(WORK)):sha(p) for p in files},
        references={str(p.relative_to(WORK)):sha(p) for p in refs},
        raw_inventory=dict(path='raw_inventory.csv',sha256=sha(HERE/'raw_inventory.csv'),
            files=len(inventory),bytes=sum(r['bytes'] for r in inventory)),
        note='Large per-atom NPZ files are retained locally and listed by hash; excluded from git to avoid bloating history. No remote push.')
    write(HERE/'source_manifest.json',manifest)
    print(json.dumps(manifest['raw_inventory']),flush=True)

if __name__=='__main__':main()
