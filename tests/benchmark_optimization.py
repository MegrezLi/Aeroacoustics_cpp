"""Paired before/after benchmarks and numerical checks; no Python physics model."""
import argparse
import hashlib
import json
import platform
import statistics
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASE = 'IEA_LB_RWT-AeroAcoustics'


def run(command):
    return subprocess.run([str(x) for x in command], check=True, capture_output=True, text=True).stdout


def executable(directory, name):
    return directory / (name + ('.exe' if platform.system() == 'Windows' else ''))


def compare(old, new, skip=0, delimiter=None):
    a = np.loadtxt(old, skiprows=skip, delimiter=delimiter)
    b = np.loadtxt(new, skiprows=skip, delimiter=delimiter)
    np.testing.assert_allclose(a, b, atol=2e-8, rtol=0)
    finite = np.isfinite(a) & np.isfinite(b)
    return dict(values=int(a.size), finite_values=int(finite.sum()),
                max_abs_error=float(np.max(np.abs(a[finite] - b[finite]))),
                byte_identical=old.read_bytes() == new.read_bytes())


def paired(measure, repeats):
    samples = {'before': [], 'after': []}
    for name in samples:
        measure(name)  # Untimed warm-up, including executable/library initialization.
    for i in range(repeats):
        order = ['before', 'after'] if i % 2 == 0 else ['after', 'before']
        for name in order:
            samples[name].append(measure(name))
    before = statistics.median(samples['before'])
    after = statistics.median(samples['after'])
    return dict(seconds=samples, median_before_s=before, median_after_s=after,
                speedup=before / after, time_reduction_percent=100 * (1 - after / before))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline_build', type=Path)
    parser.add_argument('candidate_build', type=Path)
    parser.add_argument('--wind9-case', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--baseline-commit', required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    assert args.repeats >= 3
    builds = dict(before=args.baseline_build.resolve(), after=args.candidate_build.resolve())
    result = dict(baseline_commit=args.baseline_commit, platform=platform.platform(),
                  processor=platform.processor(), paired_repeats=args.repeats,
                  method='Alternating paired runs, medians after warm-up; portable Release build',
                  executable_sha256={}, full_case={}, acoustics={})
    for name, directory in builds.items():
        result['executable_sha256'][name] = hashlib.sha256(executable(directory, 'aeroacoustics_turbine').read_bytes()).hexdigest()
        run([executable(directory, 'acoustic_workspace_probe'), 'verify', directory / 'workspace.csv'])
    result['workspace_comparison'] = compare(builds['before'] / 'workspace.csv', builds['after'] / 'workspace.csv')
    cases = {'official8': ROOT / 'examples' / CASE / (CASE + '.fst'), 'wind9': args.wind9_case.resolve()}
    for case, path in cases.items():
        outputs = {name: directory / ('benchmark-' + case) for name, directory in builds.items()}

        def measure(name):
            run([executable(builds[name], 'aeroacoustics_turbine'), path, outputs[name]])
            metadata = json.loads((outputs[name] / 'run.json').read_text())
            assert metadata['steps'] == 3200 and metadata['acoustic_samples'] == 201
            return metadata['elapsed_seconds']

        entry = paired(measure, args.repeats)
        entry['outputs'] = {}
        for name in ['dynamics.csv', *[f'{CASE}_{k}.out' for k in range(1, 5)]]:
            entry['outputs'][name] = compare(outputs['before'] / name, outputs['after'] / name,
                                            skip=1 if name.endswith('.csv') else 3,
                                            delimiter=',' if name.endswith('.csv') else None)
        result['full_case'][case] = entry
        print(case, json.dumps({k: v for k, v in entry.items() if k != 'outputs'}), flush=True)
    for model in ['bpm', 'tno']:
        for observers in [2, 16]:
            iterations = 20 if model == 'bpm' else 2
            checksums = {}

            def measure(name):
                data = json.loads(run([executable(builds[name], 'acoustic_workspace_probe'),
                                       'benchmark', model, observers, iterations]))
                checksums[name] = data['checksum']
                return data['seconds'] / iterations

            entry = paired(measure, args.repeats)
            np.testing.assert_allclose(checksums['before'], checksums['after'], rtol=1e-12, atol=1e-8)
            entry.update(nodes=30, observers=observers, frequencies=34, seconds_basis='one snapshot',
                         checksum_difference=abs(checksums['before'] - checksums['after']))
            result['acoustics'][f'{model}-{observers}'] = entry
            print(model, observers, entry['speedup'], flush=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
