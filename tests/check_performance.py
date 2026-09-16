"""P4-P7 CLI regressions; Python only schedules C++ and compares saved Fortran output."""
import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path
from compare_full_case import CASE, run

ROOT = Path(__file__).resolve().parents[1]


def variant(destination, changes):
    shutil.copytree(ROOT / 'examples' / CASE, destination)
    for key, value in changes.items():
        count = 0
        for path in destination.rglob('*'):
            if not path.is_file():
                continue
            text, n = re.subn(r'(?m)^(\s*)\S+(\s+' + key + r'\b)',
                              lambda m: m[1] + str(value) + m[2], path.read_text(encoding='utf-8-sig'))
            if n:
                path.write_text(text, encoding='utf-8')
                count += n
        assert count == 1, key
    return destination / (CASE + '.fst')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('build', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    build, output = args.build.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    suffix = '.exe' if (build / 'aeroacoustics_turbine.exe').exists() else ''
    exe = build / ('aeroacoustics_turbine' + suffix)
    report = {}

    def execute(case, destination, *options):
        result = subprocess.run([str(exe), str(case), str(destination), *map(str, options)],
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        return json.loads((destination / 'run.json').read_text())

    for wind in (8, 9):
        case = variant(output / f'wind{wind}-inputs', {'HWindSpeed': wind})
        destination = output / f'wind{wind}-scaled'
        metadata = execute(case, destination, 20, '--solver=scaled', '--observer-block-size=3')
        assert metadata['structural_mode'] == 'scaled'
        assert metadata['structural_last_scaled_residual'] <= 1
        reference = ROOT / 'reference/results' / (CASE + ('-wind9' if wind == 9 else ''))
        report[f'scaled_wind{wind}'] = {'fortran': run(reference, destination), 'runtime': metadata}

    outputs = {}
    for count in (1, 2, 3, 4):
        case = variant(output / f'count{count}-inputs', {'NrOutFile': count})
        destination = output / f'count{count}'
        execute(case, destination, .1, '--observer-block-size=1')
        outputs[count] = destination
        for missing in range(count + 1, 5):
            assert not (destination / f'{CASE}_{missing}.out').exists()
    for count, destination in outputs.items():
        for number in range(1, count + 1):
            for extension in ('out', 'mask'):
                name = f'{CASE}_{number}.{extension}'
                assert (destination / name).read_bytes() == (outputs[4] / name).read_bytes()
    report['output_counts_1_to_4'] = 'requested files byte-identical; unused files absent'

    tno = variant(output / 'tno-inputs', {'TBLTEMod': 2})
    result = subprocess.run([str(build / ('performance_api_probe' + suffix)), str(tno),
                             str(output / 'tno-api')], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    report['tno_batch_api'] = result.stdout
    for option in ('--solver=wrong', '--observer-block-size=0', '--observer-block-size=-1',
                   '--observer-block-size=1 --observer-block-size=2'):
        result = subprocess.run([str(exe), str(tno), str(output / 'invalid'), '0', *option.split()],
                                capture_output=True, text=True)
        assert result.returncode != 0, option
    report['invalid_cli_options'] = 'passed'
    (output / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
