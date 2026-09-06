"""Fetch pinned upstream source and full official case inputs, into ignored _work."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import tarfile
import urllib.request
import time

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'_work'
OF='2895884d2be01862173c88d70f86b358d2f1a50a'
RT='dd5feaaaa500ba7283140107806300d551cff0a7'
CASE='glue-codes/openfast/IEA_LB_RWT-AeroAcoustics/'
WORK.mkdir(exist_ok=True)

def fetch(url,dest):
    if dest.exists():return
    for attempt in range(3):
        try:
            data=urllib.request.urlopen(url,timeout=90).read()
            dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(data)
            return
        except Exception:
            if attempt==2:raise
            time.sleep(1)

archive=WORK/'openfast.tar.gz'
fetch(f'https://codeload.github.com/OpenFAST/openfast/tar.gz/{OF}',archive)
source=WORK/'openfast'
if not (source/'.source_complete').exists():
    with tarfile.open(archive) as stream:
        for member in stream.getmembers():
            target=(WORK/member.name).resolve()
            if not target.is_relative_to(WORK.resolve()) or member.issym() or member.islnk():
                raise ValueError('Unsafe archive member')
        for member in stream.getmembers():
            member.name='openfast/'+member.name.split('/',1)[-1]
            if member.name.startswith('openfast/docs/'):continue
            stream.extract(member,WORK,filter='data')
    (source/'.source_complete').write_text(OF)
print('OpenFAST source ready:',source,flush=True)
tree_file=WORK/'rtest_tree.json'
fetch(f'https://api.github.com/repos/OpenFAST/r-test/git/trees/{RT}?recursive=1',tree_file)
tree=json.loads(tree_file.read_text())['tree']
entries=[x for x in tree if x['type']=='blob' and x['path'].startswith(CASE)]
inputs=[x for x in entries if not x['path'].endswith(('.out','.outb','.log'))]
def download(entry):
    name=entry['path']
    fetch(f'https://raw.githubusercontent.com/OpenFAST/r-test/{RT}/{name}',WORK/'r-test'/name)
    data=(WORK/'r-test'/name).read_bytes()
    actual=hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()
    if actual!=entry['sha']:raise RuntimeError(f'Git blob hash mismatch: {name}')
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
    list(pool.map(download,inputs))
print('Official case inputs ready:',len(inputs),flush=True)
provenance={'openfast_commit':OF,'rtest_commit':RT,
    'openfast_archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
    'input_git_blobs':{x['path']:x['sha'] for x in inputs},'all_input_blobs_verified':True}
(ROOT/'docs').mkdir(exist_ok=True)
(ROOT/'docs/source-provenance.json').write_text(json.dumps(provenance,indent=2))
print('Reference output sizes:',[(x['path'].split('/')[-1],x.get('size')) for x in entries if x not in inputs],flush=True)
