"""Compare the four acoustic outputs from original Fortran and standalone C++."""
from pathlib import Path
import argparse
import json
import re
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
CASE='IEA_LB_RWT-AeroAcoustics'

def acoustic(path):
    lines=path.read_text().splitlines()
    start=next(i for i,line in enumerate(lines) if re.match(r'^\s*0(?:\.0*)?\s',line))
    data=np.loadtxt(lines[start:])
    if not np.all(np.isfinite(data)):raise AssertionError(f'Nonfinite acoustic output: {path}')
    return data,lines[start-2].split(),lines[start-1].split()

def run(baseline,candidate):
    report={'case':CASE,'comparison':'original Fortran versus standalone C++','files':{}}
    # Full-case bounds cover solver convergence and the reference text precision.
    tolerances=[0.001,0.01,0.01,0.05]
    for k,tolerance in enumerate(tolerances,1):
        name=f'{CASE}_{k}.out'
        a,_,_=acoustic(baseline/name);b,_,_=acoustic(candidate/name)
        assert a.shape==b.shape,'Output dimensions differ'
        np.testing.assert_allclose(a[:,0],b[:,0],rtol=0,atol=1e-12)
        assert a.shape[0]==201 and a[0,0]==0 and a[-1,0]==20,'Expected the unshortened 20 s case'
        difference=np.abs(a[:,1:]-b[:,1:]);row,column=np.unravel_index(np.argmax(difference),difference.shape)
        report['files'][name]={'rows':a.shape[0],'channels':a.shape[1]-1,'compared_values':difference.size,
            'max_abs_error_dB':float(difference[row,column]),'mean_abs_error_dB':float(np.mean(difference)),
            'worst_time_s':float(a[row,0]),'worst_channel':int(column+1),'tolerance_dB':tolerance}
        np.testing.assert_allclose(a[:,1:],b[:,1:],rtol=0,atol=tolerance,err_msg=name)
    report['compared_values']=sum(x['compared_values'] for x in report['files'].values())
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('candidate',type=Path)
    p.add_argument('--report',type=Path,default=ROOT/'docs/validation-full-case.json')
    args=p.parse_args();result=run(args.baseline,args.candidate)
    args.report.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
