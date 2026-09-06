"""Run the unshortened official 20 s case in its own input/output directory."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
ROOT=Path(__file__).resolve().parents[1]
CASE='IEA_LB_RWT-AeroAcoustics'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--label',required=True)
    parser.add_argument('--case',type=Path,default=ROOT/'_work/r-test/glue-codes/openfast'/CASE)
    args=parser.parse_args()
    if not args.label.replace('-','').replace('_','').isalnum():raise ValueError('Invalid run label')
    destination=ROOT/'_work/runs'/args.label
    if destination.exists():raise FileExistsError(f'Preserve earlier run: {destination}')
    shutil.copytree(args.case,destination)
    inputs={str(p.relative_to(destination)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in destination.rglob('*') if p.is_file()}
    executable=args.executable.resolve();started=time.perf_counter()
    print(f'Running {executable} in {destination}',flush=True)
    with (destination/'run.log').open('w') as log:
        result=subprocess.run([str(executable),CASE+'.fst'],cwd=destination,stdout=log,stderr=subprocess.STDOUT)
    metadata={'executable_sha256':hashlib.sha256(executable.read_bytes()).hexdigest(),
              'returncode':result.returncode,'elapsed_seconds':time.perf_counter()-started,'inputs_sha256':inputs}
    (destination/'run.json').write_text(json.dumps(metadata,indent=2))
    tail=(destination/'run.log').read_text(errors='replace')[-6000:];print(tail)
    if result.returncode!=0 or 'OpenFAST terminated normally' not in tail:
        raise RuntimeError('OpenFAST did not terminate normally; inspect run.log')
    print(f'Completed: {destination}',flush=True)

if __name__=='__main__':main()
