"""Compare pre-refactor and current C++ files; no physics implementation in Python."""
import argparse
import hashlib
import json
import platform
import statistics
import subprocess
from pathlib import Path

from support import CASE, variant


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('before', type=Path)
    p.add_argument('after', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--baseline-commit', required=True)
    p.add_argument('--rounds', type=int, default=1, help='Use 3+ rounds for alternating timings after warm-up')
    args = p.parse_args()
    if args.rounds < 1:
        p.error('--rounds must be positive')
    args.output.mkdir(parents=True, exist_ok=False)
    executables = {'before': args.before.resolve(), 'after': args.after.resolve()}
    result = {'baseline_commit': args.baseline_commit,
              'platform': platform.platform(),
              'rounds': args.rounds,
              'executables_sha256': {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in executables.items()},
              'cases': {}}
    for name, key, value, duration in [('official8', None, None, 20), ('wind9', 'HWindSpeed', '9', 20),
                                       ('table_bl2', 'BLMod', '2', 20), ('tno', 'TBLTEMod', '2', .1)]:
        folder = args.output / name
        inputs = folder / 'inputs'
        variant(inputs, {key: value} if key else {})
        def execute(version):
            subprocess.run([str(executables[version]), str((inputs / (CASE + '.fst')).resolve()),
                            str((folder / version).resolve()), str(duration)], check=True,
                           capture_output=True, text=True)
            meta = json.loads((folder / version / 'run.json').read_text())
            assert meta['steps'] == (16 if name == 'tno' else 3200)
            assert meta['acoustic_samples'] == (2 if name == 'tno' else 201)
            return meta['elapsed_seconds']

        timings = {version: [] for version in executables}
        if args.rounds > 1:
            for version in executables:
                execute(version)
        for trial in range(args.rounds):
            for version in (['before', 'after'] if trial % 2 == 0 else ['after', 'before']):
                timings[version].append(execute(version))
        hashes = {}
        assert {f.name for f in (folder / 'before').iterdir()} == {
            f.name for f in (folder / 'after').iterdir()}, (name, 'output file set')
        for f in (folder / 'before').iterdir():
            if f.name == 'run.json':
                a, b = [json.loads((folder / version / f.name).read_text()) for version in executables]
                a.pop('elapsed_seconds'); b.pop('elapsed_seconds')
                a.pop('acoustic_metadata', None); b.pop('acoustic_metadata', None)
                a.pop('module_configuration', None); b.pop('module_configuration', None)
                # Newly added diagnostic counters do not change numerical output.
                a = {k: v for k, v in a.items() if not k.startswith('structural_')}
                b = {k: v for k, v in b.items() if not k.startswith('structural_')}
                assert a == b, (name, 'runtime metadata')
                continue
            before, after = f.read_bytes(), (folder / 'after' / f.name).read_bytes()
            assert before == after, (name, f.name)
            hashes[f.name] = hashlib.sha256(after).hexdigest()
        result['cases'][name] = {'duration_s': duration, 'byte_identical_sha256': hashes,
                                'seconds': timings}
        if args.rounds > 1:
            medians = {k: statistics.median(v) for k, v in timings.items()}
            result['cases'][name].update(median_seconds=medians,
                                        speedup=medians['before'] / medians['after'])
        print(name, 'all numerical files and lookup report byte-identical', flush=True)
        (args.output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
