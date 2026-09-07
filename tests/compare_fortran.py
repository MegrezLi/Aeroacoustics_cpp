"""Compare the C++ shared library directly with the original Fortran routines."""
from pathlib import Path
import argparse
import ctypes as ct
import json
import os
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def run(library,fortran_compiler):
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
    from fortran_reference import build
    reference_library,command=build(fortran_compiler,'double')
    reference_fn=reference_library.__acoustic_reference_MOD_evaluate
    reference_fn.restype=None
    frequencies=np.array([10,12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000],dtype=float)
    mapping={1:slice(0,1),2:slice(1,4),3:slice(4,5),4:slice(5,6),5:slice(6,7),6:slice(7,8),7:slice(8,10)}
    for index in range(204):
        x=np.array([rng.uniform(-3,22),rng.uniform(.1,4),rng.uniform(15,100),rng.uniform(10,170),
            rng.uniform(15,165),rng.uniform(.1,5),rng.uniform(1,500),12.5,.00746583,.0110586,.0012,
            rng.uniform(.01,.2),rng.uniform(.0002,.008),rng.uniform(0,14),.02,.12,340.,1.48e-5,1.225,40.,1.,
            .000378576,.00198438,.0110586,1.,1.],dtype=float)
        if index>=180:x[0]=[-3,0,1.333,3,5,7.5,12.5,20][(index-180)//3];x[11]=0.
        flags=np.array([1+index%2,index%3,index%2],dtype=np.int32)
        reference_input=x[:22].copy()
        reference_flags=np.array([*flags,int(index<16)],dtype=np.int32)
        reference_result=np.zeros((len(frequencies),10),order='F')
        n=ct.c_int(len(frequencies))
        reference_fn(reference_input.ctypes.data_as(fp),reference_flags.ctypes.data_as(ip),
                     frequencies.ctypes.data_as(fp),ct.byref(n),reference_result.ctypes.data_as(fp))
        operations=range(1,8 if index<16 else 7)
        for op in operations:
            args=x.copy()
            if op==4:args[0]=np.radians(args[0])
            out=np.zeros((3,len(frequencies)))
            status=lib.aeroacoustics_kernel(op,len(frequencies),frequencies.ctypes.data_as(fp),args.ctypes.data_as(fp),flags.ctypes.data_as(ip),out.ctypes.data_as(fp))
            if status:raise RuntimeError(lib.aeroacoustics_last_error().decode())
            reference=reference_result[:,mapping[op]].T;actual=out[:reference.shape[0]]
            if not np.all(np.isfinite(actual)):raise AssertionError(f'Nonfinite result: case {index}, op {op}')
            error=float(np.max(abs(actual-reference)))
            maximum[op]=max(maximum.get(op,0.),error);counts[op]=counts.get(op,0)+actual.size
            np.testing.assert_allclose(actual,reference,rtol=0,atol=2e-8,err_msg=f'case={index}, op={op}')
    return {'backend':lib.aeroacoustics_backend().decode(),'cases':204,'TNO_cases':16,
            'evaluated_values':sum(counts.values()),'values_by_operation':counts,
            'fortran_max_abs_error_dB':maximum,
            'fortran_literal_precision':'promoted to double',
            'tolerance_dB':2e-8,'kernel_calls':lib.aeroacoustics_call_count()}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('library');parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--fortran-compiler',required=True)
    args=parser.parse_args();result=run(args.library,args.fortran_compiler)
    args.report.parent.mkdir(exist_ok=True,parents=True);args.report.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
