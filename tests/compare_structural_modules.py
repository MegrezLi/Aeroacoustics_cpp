"""Compile original Fortran module probes and compare their C++ counterparts."""
import argparse
import json
from pathlib import Path
import subprocess
import os
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--openfast-build',type=Path,required=True)
    p.add_argument('--cpp-build',type=Path,required=True)
    p.add_argument('--case',type=Path,default=ROOT/'examples/IEA_LB_RWT-AeroAcoustics')
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--link-library',action='append',default=[])
    p.add_argument('--report',type=Path,required=True)
    args=p.parse_args();args.work.mkdir(parents=True,exist_ok=True)
    suffix='.exe' if os.name=='nt' else ''
    libs={'structure':['elastodyn/libelastodynlib.a'],'mesh':[],
          'bem':['aerodyn/libaerodynlib.a','aerodyn/libbasicaerolib.a']}
    report={'reference':'original pinned OpenFAST modules, double precision','modules':{}}
    for name,dependencies in libs.items():
        exe=(args.work/(name+'-reference'+suffix)).resolve()
        libraries=[args.openfast_build/'modules'/d for d in dependencies+['nwtc-library/libnwtclibs.a']]
        subprocess.run([args.compiler,'-ffree-line-length-none','-fdefault-real-8','-fdefault-double-8',
                        f'-I{args.openfast_build / "ftnmods"}',str(ROOT/f'tests/fortran/turbine_{name}_reference.f90'),
                        *map(str,libraries),*(args.link_library or ['-llapack','-lblas']),'-o',str(exe)],check=True)
        cpp=args.work/(name+'-cpp.csv');ftn=args.work/(name+'-fortran.csv');fixture=args.work/(name+'-fixture.txt')
        probe=(args.cpp_build/f'turbine_{name}_probe{suffix}').resolve()
        case=args.case/'IEA_LB_RWT-AeroAcoustics.fst'
        if name=='structure':
            subprocess.run([str(probe),str(case),str(cpp)],check=True)
            subprocess.run([str(exe),str(args.case/'RotorSE_FAST_IEA_landBased_RWT_ElastoDyn.dat'),str(ftn)],check=True)
        else:
            inputs=[str(case)] if name=='bem' else []
            subprocess.run([str(probe),*inputs,str(fixture),str(cpp)],check=True)
            subprocess.run([str(exe),str(fixture),str(ftn)],check=True)
        a=np.genfromtxt(ftn,delimiter=',',names=True);b=np.genfromtxt(cpp,delimiter=',',names=True)
        assert a.shape==b.shape and a.dtype.names==b.dtype.names
        errors={}
        for column in a.dtype.names:
            assert np.all(np.isfinite(a[column])) and np.all(np.isfinite(b[column]))
            np.testing.assert_allclose(a[column],b[column],rtol=1e-11,atol=1e-9,err_msg=f'{name}: {column}')
            errors[column]=float(np.max(abs(a[column]-b[column])))
        report['modules'][name]={'rows':len(a),'max_absolute_error':errors,'absolute_tolerance':1e-9,'relative_tolerance':1e-11}
    args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
