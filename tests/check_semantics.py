"""Check S4-S5 output metadata and unsupported module rejection; no acoustic model in Python."""
import argparse
import json
import subprocess
from pathlib import Path
from check_performance import variant, ROOT, CASE

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    exe, out = args.executable.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {}
    for weighting in ('False', 'True'):
        case = variant(out / ('inputs-' + weighting), {'AWeighting': weighting})
        destination = out / weighting
        subprocess.run([str(exe), str(case), str(destination), '.1'], check=True, capture_output=True)
        metadata = json.loads((destination / 'run.json').read_text())
        semantics = metadata['acoustic_metadata']
        assert semantics['linear_quantity'] == 'relative_mean_square_pressure'
        assert semantics['level_quantity'] == 'sound_pressure_level'
        assert semantics['reference_pressure_Pa'] == 2e-5
        assert semantics['weighting'] == ('A' if weighting == 'True' else 'unweighted')
        bands = semantics['bands_center_lower_upper_Hz']
        assert len(bands) == 34
        for i, (center, low, high) in enumerate(bands):
            assert 0 < low < center < high
            if i:
                assert abs(low - bands[i-1][2]) < 1e-10
        for number in range(1,5):
            lines = (destination / f'{CASE}_{number}.out').read_text().splitlines()
            unit = '(dBA)' if weighting == 'True' else '(dB)'
            assert set(lines[2].split()[1:]) == {unit}
        assert metadata['module_configuration']['dof_count'] == 9
        dofs = metadata['module_configuration']['dofs']
        assert len({d['name'] for d in dofs}) == 9
        assert all(d['unit'] == 'm' and d['acceleration_scale'] == 1 for d in dofs)
        blocks = metadata['module_configuration']['coupling_blocks']
        assert [(b['offset'], b['size']) for b in blocks] == [(0, 3), (3, 3), (6, 3)]
        assert len((destination / 'dynamics.csv').read_text().splitlines()[0].split(',')) == 19
        report[weighting] = {'acoustic_metadata': semantics, 'module_configuration': metadata['module_configuration']}
    for key, value, error in [('CompServo','1','Unsupported enabled module'),
                              ('TwFADOF1','True','Unsupported active DOF'),
                              ('DBEMT_Mod','1','Unsupported wake model')]:
        case = variant(out / key, {key:value})
        result = subprocess.run([str(exe),str(case),str(out / (key+'-out')),'0'],capture_output=True,text=True)
        assert result.returncode != 0 and error in result.stderr, result.stderr
    report['unsupported_modules'] = 'servo, tower DOF, dynamic wake rejected before solving'
    (out / 'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('S4-S5 CLI metadata and module guards passed')

if __name__ == '__main__':
    main()
