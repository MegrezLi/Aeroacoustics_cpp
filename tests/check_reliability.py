"""Exercise failure handling, zero-energy masks and lookup policies through the CLI."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASE = 'IEA_LB_RWT-AeroAcoustics'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    exe, output = args.executable.resolve(), args.output.resolve()
    # Isolated new fixtures avoid overwriting prior run evidence.
    output.mkdir(parents=True, exist_ok=False)
    fixture = output / 'inputs'
    shutil.copytree(ROOT / 'examples' / CASE, fixture)
    case = fixture / (CASE + '.fst')

    def run(folder, *options, success=True):
        result = subprocess.run([str(exe), str(case), str(folder), '0', *options], capture_output=True, text=True)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        if not success:
            assert 'completed' not in result.stdout
            meta = folder / 'run.json'
            if meta.is_file():
                try:
                    json.loads(meta.read_text())
                except json.JSONDecodeError:
                    pass
                else:
                    raise AssertionError('Failed run left a valid completion record')
        return result

    regular = output / 'regular'
    run(regular)
    masks = {}
    for data in regular.glob('*.out'):
        values = np.loadtxt(data, skiprows=3, ndmin=2)
        mask = np.loadtxt(data.with_suffix('.mask'), skiprows=3, ndmin=2)
        assert mask.shape == values.shape
        np.testing.assert_array_equal(mask[:, 0], values[:, 0])
        flags = mask[:, 1:]
        assert np.all((flags == 0) | (flags == 1))
        assert np.all(values[:, 1:][flags == 0] == 0)
        masks[data.name] = {'zero': int((flags == 0).sum()), 'positive': int((flags == 1).sum())}
    assert len(masks) == 4 and any(x['zero'] for x in masks.values())
    meta = json.loads((regular / 'run.json').read_text())
    assert 'matching .mask' in meta['zero_energy_encoding']

    failed_outputs = []
    prefix = CASE + '_'
    for name in ['run.json', 'dynamics.csv', 'lookup_diagnostics.csv', prefix+'1.out', prefix+'1.mask']:
        folder = output / ('blocked-' + name)
        (folder / name).mkdir(parents=True)
        result = run(folder, success=False)
        assert name in result.stderr, result.stderr
        failed_outputs.append(name)

    # Force real Airfoil::at queries outside a small but still valid sorted polar.
    polar = fixture / 'Airfoils/RotorSE_FAST_IEA_landBased_RWT_AeroDyn_Polar_10.dat'
    original = polar.read_text(encoding='utf-8')
    lines = original.splitlines()
    changes = 0
    for i, line in enumerate(lines):
        tokens = line.split()
        if len(tokens) == 4:
            try:
                values = [float(x) for x in tokens]
            except ValueError:
                continue
            values[0] *= 1e-4
            lines[i] = ' '.join(format(x, '.17g') for x in values)
            changes += 1
    assert changes > 10
    polar.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    clamped = output / 'clamped'
    result = run(clamped, '--lookup-policy=clamp')
    assert 'out-of-range lookup' in result.stderr
    entries = list(csv.DictReader((clamped / 'lookup_diagnostics.csv').open()))
    assert entries and all(int(x['blade']) > 0 and int(x['node']) > 0 for x in entries)
    assert all(x['stage'] in {'BEM_trial', 'UA_evaluate', 'UA_advance'} for x in entries)
    total = sum(int(x['calls']) for x in entries)
    assert total == json.loads((clamped / 'run.json').read_text())['lookup_out_of_range_calls']
    result = run(output / 'strict', '--lookup-policy=error', success=False)
    for token in ['Lookup out of range', 'blade=', 'node=', 'stage=', 'time=', 'alpha_deg']:
        assert token in result.stderr, result.stderr
    polar.write_text(original, encoding='utf-8')

    # Overflow an active tip model with a finite, positive parameter.
    aa = fixture / 'AeroAcousticsInput.dat'
    original = aa.read_text(encoding='utf-8')
    lines = original.splitlines()
    for i, line in enumerate(lines):
        tokens = line.split()
        if len(tokens) > 1 and tokens[1].lower() in {'tipmod', 'alprat'}:
            tokens[0] = '1' if tokens[1].lower() == 'tipmod' else '1e308'
            lines[i] = ' '.join(tokens)
    aa.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    result = run(output / 'invalid-spectrum', success=False)
    for token in ['Acoustic time=', 'observer=', 'blade=', 'node=', 'mechanism=tip', 'frequency_Hz=']:
        assert token in result.stderr, result.stderr
    report = dict(masks=masks, failed_output_paths=failed_outputs, lookup_report_groups=len(entries),
                  lookup_out_of_range_calls=total, strict_lookup='passed', invalid_spectrum='passed')
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
