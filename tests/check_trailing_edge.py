"""E4/E9 CLI and coupled regressions. Python supplies inputs and checks C++ outputs only."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from check_performance import ROOT, CASE, variant


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    build, out = args.build.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    suffix = '.exe' if (build/'aeroacoustics_turbine.exe').exists() else ''
    exe = build/('aeroacoustics_turbine'+suffix)
    case = ROOT/'examples'/CASE/(CASE+'.fst')
    fixtures = ROOT/'examples/trailing-edge'
    surfaces = '--surfaces='+str(fixtures/'surfaces.dat')
    def run(name, config=None, duration=2, input_case=case, extras=(), surface=True, fails=False):
        command = [str(exe), str(input_case), str(out/name), str(duration)]
        if surface: command.append(surfaces)
        if config: command.append('--trailing-edge='+str(config))
        p = subprocess.run(command+list(extras),capture_output=True,text=True)
        assert (p.returncode != 0) == fails, (name,p.stdout,p.stderr)
        if fails:
            meta=out/name/'run.json'
            assert not meta.exists() or not meta.read_text().strip(), name
            return p.stderr
        assert json.loads((out/name/'run.json').read_text())['solver']=='standalone C++'
        return out/name
    def data(folder,number):return np.loadtxt(folder/f'{CASE}_{number}.out',skiprows=3)
    bpm=run('bpm')
    straight=run('straight',fixtures/'howe-straight.dat')
    serrated=run('serrated',fixtures/'howe-serrated.dat')
    assert (bpm/'dynamics.csv').read_bytes()==(straight/'dynamics.csv').read_bytes()==(serrated/'dynamics.csv').read_bytes()
    assert not np.array_equal(data(bpm,1),data(straight,1))
    assert not np.array_equal(data(straight,1),data(serrated,1))
    # Only the two trailing-edge sides and their replacement of BPM separation may change.
    a,b,c=(data(folder,3)[:,1:].reshape(-1,2,34,7) for folder in (bpm,straight,serrated))
    np.testing.assert_array_equal(a[:,:,:,[0,4,5,6]],b[:,:,:,[0,4,5,6]])
    np.testing.assert_array_equal(b[:,:,:,[0,4,5,6]],c[:,:,:,[0,4,5,6]])
    mask=np.loadtxt(serrated/f'{CASE}_3.out.mask',skiprows=3) if (serrated/f'{CASE}_3.out.mask').exists() else np.loadtxt(serrated/f'{CASE}_3.mask',skiprows=3)
    assert np.all(mask[:,1:].reshape(-1,2,34,7)[:,:,:,3]==0)
    blocked=run('blocked',fixtures/'howe-serrated.dat',extras=['--observer-block-size=3'])
    for path in serrated.glob('*.out'):assert path.read_bytes()==(blocked/path.name).read_bytes()
    reference=run('tno-reference',fixtures/'tno-reference.dat',duration=.2)
    supplied=run('tno-input',fixtures/'tno-input.dat',duration=.2)
    assert (reference/'dynamics.csv').read_bytes()==(supplied/'dynamics.csv').read_bytes()
    assert not np.array_equal(data(reference,2),data(supplied,2))
    for folder,mode in ((reference,'reference'),(supplied,'input')):
        assert json.loads((folder/'run.json').read_text())['trailing_edge_model']['tno_edge_velocity']==mode
    legacy_case=variant(out/'tno-case',{'TBLTEMod':2})
    default_tno=run('tno-default',duration=.2,input_case=legacy_case)
    for path in reference.glob('*.out'):assert path.read_bytes()==(default_tno/path.name).read_bytes()
    # Also exercise the optional AA-file setting, independently of CLI overrides.
    aa=out/'tno-case/AeroAcousticsInput.dat'
    with aa.open('a') as f:f.write('\n"input" TNOEdgeVelocity\n')
    aa_input=run('tno-aa-input',duration=.2,input_case=legacy_case)
    for path in supplied.glob('*.out'):assert path.read_bytes()==(aa_input/path.name).read_bytes()
    run('missing-bl',fixtures/'howe-straight.dat',surface=False,fails=True)
    bad=out/'bad.dat';bad.write_text((fixtures/'howe-serrated.dat').read_text().replace('0.04 WavelengthM','0 WavelengthM'))
    run('bad-geometry',bad,fails=True)
    run('duplicate',fixtures/'howe-serrated.dat',extras=['--trailing-edge='+str(fixtures/'howe-straight.dat')],fails=True)
    bad.write_text('"test" Provenance\n"tno" Model\n"unknown" TNOEdgeVelocity\n')
    run('bad-mode',bad,fails=True)
    bad.write_text((fixtures/'howe-serrated.dat').read_text()+'\n0.02 HalfHeightM\n')
    run('duplicate-field',bad,fails=True)
    # E1/E2/E7 combination: retained receiver-time histories consume the new source normally.
    metrics=(ROOT/'examples/engineering/metrics.dat').read_text()
    import re
    for key,value in {'Start':.2,'End':2,'RetardedTime':'false','NumTones':0,'ObserverFile':'none'}.items():
        metrics,n=re.subn(r'(?m)^\S+(\s+'+key+r'\b)',lambda m:str(value)+m[1],metrics);assert n==1
    (out/'metrics.dat').write_text(metrics)
    air=(ROOT/'examples/engineering/propagation.dat').read_text().replace('rigid GroundModel','none GroundModel')
    (out/'air.dat').write_text(air)
    run('combined',fixtures/'howe-serrated.dat',extras=['--controller='+str(ROOT/'examples/engineering/controller.dat'),
        '--tower='+str(ROOT/'examples/farm/tower-test.dat'),'--propagation='+str(out/'air.dat'),'--metrics='+str(out/'metrics.dat')])
    # Farm-unit configuration must reach its independent Simulation.
    unit=out/'unit.dat';unit.write_text(f'"one" Name\n"{case.as_posix()}" Case\n0 X\n0 Y\n"none" Controller\n"none" Tower\n"{(fixtures/"surfaces.dat").as_posix()}" Surfaces\n"{(fixtures/"howe-serrated.dat").as_posix()}" TrailingEdge\n')
    farm=out/'farm.dat';farm.write_text(f'"Synthetic E4 farm integration" Provenance\n0.2 Duration\n0 StatisticsStart\nfalse Wakes\n0.05 WakeExpansion\n0.95 MaxCt\n32000000 MaxValues\n"{(ROOT/"examples/farm/receivers.dat").as_posix()}" ObserverFile\n"none" Propagation\n1 NumTurbines\n"unit.dat" TurbineFiles\n0 NumSources\n')
    p=subprocess.run([str(build/('aeroacoustics_farm'+suffix)),str(farm),str(out/'farm')],capture_output=True,text=True)
    assert p.returncode==0,p.stderr
    assert json.loads((out/'farm/turbine_1/run.json').read_text())['trailing_edge_model']['selection']==3
    report={'howe_dynamics':'byte-identical to the same supplied-BL BPM case',
            'unrelated_sources':'unchanged; no BPM separation in Howe channels',
            'observer_blocks':'byte-identical',
            'serration_total_change_dB_range':[float(np.min(data(serrated,1)[:,1:]-data(straight,1)[:,1:])),float(np.max(data(serrated,1)[:,1:]-data(straight,1)[:,1:]))],
            'tno':'default/reference identical; input changes spectra; AA field and CLI agree; dynamics unchanged',
            'integration':'controller, tower, surface, absorption, metrics and farm entry passed',
            'guards':'missing attached BL, invalid geometry/mode and duplicate CLI rejected',
            'field_validated':False,'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
