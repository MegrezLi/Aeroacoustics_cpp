"""Compare all four acoustic output files and OpenFAST's binary dynamics output."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
CASE='IEA_LB_RWT-AeroAcoustics'

def acoustic(path):
    lines=path.read_text().splitlines()
    start=next(i for i,line in enumerate(lines) if re.match(r'^\s*0\.0+\s',line))
    data=np.loadtxt(lines[start:])
    if not np.all(np.isfinite(data)):raise AssertionError(f'Nonfinite acoustic output: {path}')
    return data,lines[start-2].split(),lines[start-1].split()

def run(baseline,candidate):
    inputs=[]
    for folder in (baseline,candidate):
        meta=json.loads((folder/'run.json').read_text())
        assert meta['returncode']==0
        inputs.append(meta['inputs_sha256'])
        assert 'OpenFAST terminated normally' in (folder/'run.log').read_text()
    assert inputs[0]==inputs[1],'Inputs differ between baseline and candidate'
    match=re.search(r'C\+\+ acoustic kernel calls:\s*(\d+)',(candidate/'run.log').read_text())
    assert match and int(match[1])>0,'No proof the full case called C++'
    report={'case':CASE,'input_files':len(inputs[0]),'identical_inputs':True,
            'cpp_kernel_calls':int(match[1]),'files':{}}
    for k in range(1,5):
        name=f'{CASE}_{k}.out'
        a,an,au=acoustic(baseline/name);b,bn,bu=acoustic(candidate/name)
        assert an==bn and au==bu,'Channel headers differ'
        np.testing.assert_array_equal(a[:,0],b[:,0])
        assert a.shape[0]==201 and a[0,0]==0 and a[-1,0]==20,'Case duration was changed'
        np.testing.assert_allclose(a[:,1:],b[:,1:],rtol=0,atol=2e-5,err_msg=name)
        report['files'][name]={'rows':a.shape[0],'channels':a.shape[1]-1,
            'compared_values':a[:,1:].size,'max_abs_error_dB':float(np.max(abs(a[:,1:]-b[:,1:])))}
    # OpenFAST embeds wall-clock generation time in its ASCII description.
    # Normalize only that exact field; all scaling, channel metadata and binary
    # time-series bytes must remain identical.
    name=CASE+'.outb';a=(baseline/name).read_bytes();b=(candidate/name).read_bytes()
    pattern=rb'Predictions were generated on \d{2}-[A-Za-z]{3}-\d{4} at \d{2}:\d{2}:\d{2}'
    aa,count_a=re.subn(pattern,b'Predictions were generated at <wall-clock time>',a)
    bb,count_b=re.subn(pattern,b'Predictions were generated at <wall-clock time>',b)
    assert count_a==count_b==1,'Unexpected binary description format'
    assert aa==bb,'Dynamics output differs beyond its generation timestamp'
    report['dynamics_output']={'byte_identical_except_generation_timestamp':True,'bytes':len(a),
        'normalized_sha256':hashlib.sha256(aa).hexdigest(),
        'baseline_sha256':hashlib.sha256(a).hexdigest(),'candidate_sha256':hashlib.sha256(b).hexdigest()}
    report['acoustic_tolerance_dB']=2e-5
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('candidate',type=Path)
    p.add_argument('--report',type=Path,default=ROOT/'docs/validation-full-case.json')
    args=p.parse_args();result=run(args.baseline,args.candidate)
    args.report.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
