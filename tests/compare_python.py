"""Independent executable-library comparison against the supplied Python port."""
from pathlib import Path
import argparse
import ctypes as ct
import json
import os
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference/python'))
import aeroacoustics as py

def run(library,fortran_compiler=None):
    handles=[]
    runtime_dirs = os.environ.get('AEROACOUSTICS_RUNTIME_DIRS', '').split(os.pathsep)
    for folder in filter(None, runtime_dirs):
        if os.name=='nt' and Path(folder).exists():handles.append(os.add_dll_directory(folder))
    os.environ.setdefault('MKL_THREADING_LAYER','SEQUENTIAL')
    lib=ct.CDLL(str(Path(library).resolve()))
    fp=ct.POINTER(ct.c_double);ip=ct.POINTER(ct.c_int)
    lib.aeroacoustics_kernel.argtypes=[ct.c_int,ct.c_int,fp,fp,ip,fp]
    lib.aeroacoustics_last_error.restype=ct.c_char_p
    lib.aeroacoustics_backend.restype=ct.c_char_p
    lib.aeroacoustics_call_count.restype=ct.c_ulonglong
    maximum={};counts={};rng=np.random.default_rng(391432)
    fortran_maximum={}
    if fortran_compiler:
        from fortran_reference import build
        reference_library,command=build(fortran_compiler,'double')
        reference_fn=reference_library.__acoustic_reference_MOD_evaluate
        reference_fn.restype=None
    for index in range(204):
        x=np.array([rng.uniform(-3,22),rng.uniform(.1,4),rng.uniform(15,100),rng.uniform(10,170),
            rng.uniform(15,165),rng.uniform(.1,5),rng.uniform(1,500),12.5,.00746583,.0110586,.0012,
            rng.uniform(.01,.2),rng.uniform(.0002,.008),rng.uniform(0,14),.02,.12,340.,1.48e-5,1.225,40.,1.,
            .000378576,.00198438,.0110586,1.,1.],dtype=float)
        if index>=180:x[0]=[-3,0,1.333,3,5,7.5,12.5,20][(index-180)//3];x[11]=0.
        flags=np.array([1+index%2,index%3,index%2],dtype=np.int32)
        p=py.Parameters(x_blmethod=int(flags[0]),itrip=int(flags[1]),round=bool(flags[2]))
        a,c,u,theta,phi,span,r,stall,delta,ds,dp,ti,h,psi,t1,t10=x[:16]
        expected={1:[py.lblvs(a,c,u,theta,phi,span,r,p,delta,ds,dp,stall)],
            2:py.tblte(a,c,u,theta,phi,span,r,p,delta,ds,dp,stall),
            3:[py.tipnois(a,1.,c,u,theta,phi,r,p)],
            4:[py.inflownoise(np.radians(a),c,u,theta,phi,span,r,ti,p)],
            5:[py.blunt(a,c,u,theta,phi,span,r,h,psi,p,delta,ds,dp,stall)],
            6:[py.simple_guidati(u,c,t10,t1,p)]}
        if index<16:expected[7]=py.tblte_tno(u,theta,phi,span,r,x[21:23],(x[23],delta),x[24:26],p)
        if fortran_compiler:
            reference_input=x[:22].copy()
            reference_flags=np.array([*flags,int(index<16)],dtype=np.int32)
            reference_result=np.zeros((len(p.freqlist),10),order='F')
            n=ct.c_int(len(p.freqlist))
            reference_fn(reference_input.ctypes.data_as(fp),reference_flags.ctypes.data_as(ip),
                         p.freqlist.ctypes.data_as(fp),ct.byref(n),reference_result.ctypes.data_as(fp))
            mapping={1:slice(0,1),2:slice(1,4),3:slice(4,5),4:slice(5,6),5:slice(6,7),6:slice(7,8),7:slice(8,10)}
        for op,arrays in expected.items():
            args=x.copy()
            if op==4:args[0]=np.radians(args[0])
            out=np.zeros((3,len(p.freqlist)))
            status=lib.aeroacoustics_kernel(op,len(p.freqlist),p.freqlist.ctypes.data_as(fp),args.ctypes.data_as(fp),flags.ctypes.data_as(ip),out.ctypes.data_as(fp))
            if status:raise RuntimeError(lib.aeroacoustics_last_error().decode())
            actual=out[:len(arrays)];reference=np.array(arrays)
            if not np.all(np.isfinite(actual)):raise AssertionError(f'Nonfinite result: case {index}, op {op}')
            error=float(np.max(abs(actual-reference)))
            maximum[op]=max(maximum.get(op,0.),error);counts[op]=counts.get(op,0)+actual.size
            np.testing.assert_allclose(actual,reference,rtol=0,atol=2e-8,err_msg=f'case={index}, op={op}')
            if fortran_compiler:
                original=reference_result[:,mapping[op]].T
                np.testing.assert_allclose(actual,original,rtol=0,atol=2e-8,err_msg=f'Fortran case={index}, op={op}')
                fortran_maximum[op]=max(fortran_maximum.get(op,0.),float(np.max(abs(actual-original))))
    return {'backend':lib.aeroacoustics_backend().decode(),'cases':204,'TNO_cases':16,
            'evaluated_values':sum(counts.values()),'values_by_operation':counts,
            'max_abs_error_dB':maximum,'fortran_max_abs_error_dB':fortran_maximum,
            'fortran_literal_precision':'promoted to double' if fortran_compiler else None,
            'tolerance_dB':2e-8,'kernel_calls':lib.aeroacoustics_call_count()}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('library');parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--fortran-compiler')
    args=parser.parse_args();result=run(args.library,args.fortran_compiler)
    args.report.parent.mkdir(exist_ok=True,parents=True);args.report.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
