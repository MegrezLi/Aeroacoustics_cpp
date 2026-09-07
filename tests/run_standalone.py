"""Run the bundled case and compare it with saved outputs from original Fortran."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from compare_full_case import run, CASE

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('executable',type=Path)
p.add_argument('output',type=Path)
p.add_argument('--report',type=Path,required=True)
args=p.parse_args()
baseline=ROOT/'reference/results'/CASE
case=ROOT/'examples'/CASE
meta=json.loads((baseline/'reference.json').read_text(encoding='utf-8'))
for folder,entries in [(case,meta['inputs_sha256']),(baseline,meta['output_sha256'])]:
    for name,digest in entries.items():
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest,f'Changed reference input/output: {name}'
subprocess.run([str(args.executable.resolve()),str(case/(CASE+'.fst')),str(args.output.resolve())],check=True)
runtime=json.loads((args.output/'run.json').read_text(encoding='utf-8'))
assert runtime['solver']=='standalone C++' and runtime['steps']==3200
assert runtime['dt']==0.00625 and runtime['duration']==20 and runtime['acoustic_samples']==201
report=run(baseline,args.output)
report['identical_input_files']=len(meta['inputs_sha256'])
report['runtime']=runtime
report['cpp_executable_sha256']=hashlib.sha256(args.executable.read_bytes()).hexdigest()
report['reference_openfast_commit']=meta['openfast_commit']
args.report.parent.mkdir(parents=True,exist_ok=True)
args.report.write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
