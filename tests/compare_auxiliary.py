"""Compare actual C++ driver, geometry, BL interpolation and quadrature to Python."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference/python'))
import aeroacoustics as py

def run(executable):
    with tempfile.TemporaryDirectory() as folder:
        table_path=Path(folder)/'bl.dat'
        text=['BL fixture','normalized thickness', '2 nRe','3 nAoA']
        for r in (1.,4.):
            text += [f'{r} Reynolds','header','units']
            for a in (-10.,0.,20.):
                v=np.arange(1,9)*.001*(1+r+a*.01)
                text.append(' '.join(map(str,[a,*v])))
        table_path.write_text('\n'.join(text))
        output=subprocess.check_output([str(Path(executable).resolve()),str(table_path)],text=True)
        table=py.BLTable.read(table_path)
    expected={};actual={}
    for line in output.splitlines():
        key,values=line.split(' ',1);actual[key]=np.fromstring(values,sep=' ')
    expected['bl']=[]
    for a in (-30.,-2.,15.,40.):
        for r in (.1e6,2e6,9e6):
            bl=table.interpolate(a,r,2.7)
            expected['bl'].extend((*bl.dstar,*bl.d99,*bl.cf,*bl.edge_velocity_ratio))
    expected['elements']=[]
    for percent in (20.,70.,100.):
        first,lengths=py.blade_elements([0,1,5,10,16],percent)
        expected['elements'].extend((first,*lengths))
    expected['quadrature']=np.array([py.qk61(lambda x:np.exp(x)*np.cos(3*x),-.3,u,True) for u in (-1.,0.,.5,3.)]).ravel()
    expected['thickness']=py.guidati_thickness([[.25,0],[1,0],[.1,-.05],[.01,-.02],[0,0],[.01,.03],[.1,.07],[1,0]])
    p=py.Parameters(itrip=0,lammod=1,tipmod=1,bluntmod=1,timod=2,aweighting=True)
    for method in (1,2):
        key=f'driver{method}';expected[key]=[]
        driver=py.AcousticDriver(p,[1,5,10],2,[[175,0,2],[0,175,2]],dt=1.,hub_height=10.,ti_method=method)
        for step in range(9):
            blades=[]
            for b in range(2):
                nodes=[]
                for j in range(3):
                    a=.1*(step+b);c,s=np.cos(a),np.sin(a)
                    nodes.append(dict(speed=40+3*j+b,chord=.5+.2*j,alpha_deg=2+j,
                        aero_center=[0.,b*2.,10.+(j+1)*3],inflow=[8.+np.sin(step+b+j),.2*step,0],
                        global_to_local=[[c,-s,0],[s,c,0],[0,0,1]]))
                blades.append(nodes)
            expected[key].extend(driver.step(step,blades).ravel())
            expected[key].extend(driver.state.values.ravel())
    report={}
    for key,wanted in expected.items():
        wanted=np.array(wanted);got=actual[key]
        np.testing.assert_allclose(got,wanted,rtol=0,atol=2e-8,err_msg=key)
        finite=np.isfinite(wanted)
        report[key]={'values':len(wanted),'max_abs_error':float(np.max(abs(got[finite]-wanted[finite])))}
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('executable');parser.add_argument('--report',required=True)
    args=parser.parse_args();report=run(args.executable)
    Path(args.report).write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
