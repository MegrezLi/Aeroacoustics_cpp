"""Paired P4-P7 timings. Run without concurrent builds/tests; no physics in this script."""
import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import time
from pathlib import Path
from check_performance import variant, ROOT, CASE


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('before', type=Path, help='fd193c7 build directory')
    p.add_argument('after', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--rounds', type=int, default=5)
    args = p.parse_args()
    assert args.rounds > 0
    before, after, output = args.before.resolve(), args.after.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    suffix = '.exe' if (after / 'aeroacoustics_turbine.exe').exists() else ''
    executable = lambda build, name: str(build / (name + suffix))
    report = {'baseline_commit': 'fd193c7f9fb6b467e5f9f06c25e1b187f327691c',
              'platform': platform.platform(), 'logical_cpus': os.cpu_count(),
              'rounds': args.rounds, 'results': {}}
    report['executable_sha256'] = {name: hashlib.sha256(Path(executable(build, 'aeroacoustics_turbine')).read_bytes()).hexdigest()
                                   for name, build in [('before', before), ('after', after)]}

    def run(command):
        start = time.perf_counter()
        result = subprocess.run(list(map(str, command)), check=True, capture_output=True, text=True)
        return time.perf_counter() - start, result.stdout

    def paired(name, commands, kernel=False):
        times = {key: [] for key in commands}
        for command in commands.values():
            run(command)  # Warm-up outside timed rounds.
        checksums = {}
        for trial in range(args.rounds):
            keys = list(commands)
            if trial % 2:
                keys.reverse()
            for key in keys:
                wall, stdout = run(commands[key])
                if kernel:
                    payload = json.loads(stdout)
                    times[key].append(payload['seconds'])
                    checksums[key] = payload['checksum']
                else:
                    times[key].append(wall)
        if kernel:
            assert len(set(checksums.values())) == 1
        medians = {key: statistics.median(values) for key, values in times.items()}
        keys = list(medians)
        report['results'][name] = {'seconds': times, 'median_seconds': medians,
                                   'speedup': medians[keys[0]] / medians[keys[1]],
                                   'clock': 'C++ warmed loop' if kernel else 'external wall clock'}
        (output / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(name, medians, flush=True)

    def identical(a, b):
        for file in a.iterdir():
            if file.name != 'run.json':
                assert file.read_bytes() == (b / file.name).read_bytes(), file.name

    for model in ('bpm', 'tno'):
        for observers in (2, 16):
            paired(f'{model}_{observers}_observers', {
                name: [executable(build, 'acoustic_workspace_probe'), 'benchmark', model, observers,
                       3 if model == 'tno' else 30]
                for name, build in [('before', before), ('after', after)]}, kernel=True)
    for count in (1, 4):
        case = variant(output / f'count{count}-inputs', {'NrOutFile': count})
        paired(f'full_case_output{count}', {name: [executable(build, 'aeroacoustics_turbine'), case,
                                                  output / f'count{count}-{name}', 20]
                                           for name, build in [('before', before), ('after', after)]})
        identical(output / f'count{count}-before', output / f'count{count}-after')
    case = ROOT / 'examples' / CASE / (CASE + '.fst')
    paired('structural_mode_full_case', {mode: [executable(after, 'aeroacoustics_turbine'), case,
                                               output / mode, 20, '--solver=' + mode]
                                        for mode in ('reference', 'scaled')})
    report['structural_diagnostics'] = {
        mode: {k: v for k, v in json.loads((output / mode / 'run.json').read_text()).items()
               if k.startswith('structural_')} for mode in ('reference', 'scaled')}
    # Four complete independent cases; never parallelize dependent solver time steps.
    commands = {}
    for workers in (1, 4):
        commands[str(workers)] = [executable(after, 'aeroacoustics_batch'), workers]
        for job in range(4):
            commands[str(workers)] += [case, output / f'batch{workers}-{job}']
    paired('four_case_batch_workers', commands)
    for job in range(4):
        identical(output / f'batch1-{job}', output / f'batch4-{job}')
    report['numerical_files'] = 'byte-identical before/after and serial/parallel; scaled mode tested separately'
    (output / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
