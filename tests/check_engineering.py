"""Run C++ engineering configurations and compare outputs; no solver implemented in Python."""
import argparse
import csv
import json
import re
import subprocess
from pathlib import Path
import numpy as np
from check_performance import variant, ROOT, CASE

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    exe,out=args.executable.resolve(),args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'examples'/CASE/(CASE+'.fst')
    examples=ROOT/'examples/engineering'
    normal=out/'normal-controller.dat'
    text=(examples/'controller.dat').read_text()
    normal.write_text(re.sub(r'(?m)^10 NoiseStart','100 NoiseStart',text))
    def execute(name,case=examples/'closed_loop.fst',extra=(),duration=20):
        destination=out/name
        proc=subprocess.run([str(exe),str(case),str(destination),str(duration),*extra],capture_output=True,text=True)
        assert proc.returncode==0,proc.stderr
        meta=json.loads((destination/'run.json').read_text())
        levels=np.loadtxt(destination/(CASE+'_1.out'),skiprows=3)
        operation=None
        if(destination/'operation.csv').exists():
            operation=np.genfromtxt(destination/'operation.csv',delimiter=',',names=True)
            for key in operation.dtype.names:
                assert np.isfinite(operation[key]).all(),key
            assert operation['rotor_speed_rad_s'].min()>0
            assert np.abs(operation['pitch_rate_rad_s']).max()<=np.deg2rad(8)+1e-12
            assert np.abs(operation['yaw_rate_rad_s']).max()<=np.deg2rad(.3)+1e-12
            assert operation['generator_torque_Nm'].min()>=0
            assert operation['generator_torque_Nm'].max()<=45000
            np.testing.assert_allclose(operation['electrical_power_W'],.95*operation['generator_torque_Nm']*operation['generator_speed_rad_s'],rtol=1e-14)
        return meta,levels,operation
    wind='--wind-grid='+str(examples/'gust_veer.wind')
    propagation='--propagation='+str(examples/'propagation.dat')
    control='--controller='+str(examples/'controller.dat')
    ref=execute('fixed',case=base)
    normal_run=execute('normal',extra=[wind,'--controller='+str(normal)])
    noise=execute('noise',extra=[wind,control])
    propagated=execute('propagated',extra=[wind,control,propagation])
    assert not normal_run[2]['noise_mode'].any()
    assert noise[2]['noise_mode'][-1]==1
    assert noise[2]['pitch_rad'][-1]>normal_run[2]['pitch_rad'][-1]
    assert noise[2]['rotor_speed_rad_s'][-1]<normal_run[2]['rotor_speed_rad_s'][-1]
    assert np.max(np.abs(noise[1][:,1:]-normal_run[1][:,1:]))>.01
    assert (propagated[0]['engineering']['propagation_parameters']['ground']=='rigid')
    # On the original ground-level receivers, a rigid image duplicates the source
    # path, with absorption superposed. Check the net change stays finite.
    assert np.isfinite(propagated[1]).all()
    assert np.max(np.abs(propagated[1][:,1:]-noise[1][:,1:]))>1
    # Resolve explicit aero/controller exchange by halving the main step twice.
    levels=[]; speeds=[]
    for div in (1,2,4,8):
        case=variant(out/f'inputs-dt{div}',{})
        text=case.read_text(encoding='utf-8-sig')
        text,count=re.subn(r'(?m)^(\s*)\S+(\s+DT\b)',lambda m:m[1]+str(.00625/div)+m[2],text)
        assert count==1
        case.write_text(text,encoding='utf-8')
        run=execute(f'dt{div}',case,[wind,control],duration=2)
        levels.append(run[1][:,1:]);speeds.append(run[2]['rotor_speed_rad_s'][::div])
    coarse=float(np.max(np.abs(levels[0]-levels[3])))
    fine=float(np.max(np.abs(levels[1]-levels[3])))
    quarter=float(np.max(np.abs(levels[2]-levels[3])))
    assert coarse<1 and quarter<fine<.08,(coarse,fine,quarter)
    # Compare the whole common time grid; one final point can have cancellation.
    speed_errors=[float(np.max(np.abs(s-speeds[3]))) for s in speeds[:3]]
    assert speed_errors[2]<speed_errors[1]<speed_errors[0]<.001
    # Input validation: bounds are errors, never silent wind extrapolation.
    proc=subprocess.run([str(exe),str(base),str(out/'outside'),'21',wind],capture_output=True,text=True)
    assert proc.returncode!=0 and 'outside grid' in proc.stderr
    assert (out/'outside/run.json').stat().st_size==0
    report={'closed_loop_normal_and_noise':'passed','coupled_steps':12800,
            'noise_final_rotor_rad_s':float(noise[2]['rotor_speed_rad_s'][-1]),
            'normal_final_rotor_rad_s':float(normal_run[2]['rotor_speed_rad_s'][-1]),
            'noise_vs_normal_max_abs_OASPL_dB':float(np.max(np.abs(noise[1][:,1:]-normal_run[1][:,1:]))),
            'time_step_OASPL_max_abs_to_eighth_dB':{'dt':coarse,'half_dt':fine,'quarter_dt':quarter},
            'time_step_rotor_speed_max_abs_to_eighth_rad_s':speed_errors,
            'controller_parameters_calibrated':False,'new_models_compared_to_field_measurements':False,
            'propagation':propagated[0]['engineering'],'wind_extrapolation':'rejected'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
