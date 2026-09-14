"""Paired P1-P3 whole-turbine runs; Python only schedules and compares C++ output."""
import argparse
import hashlib
import json
import platform
import re
import shutil
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASE = 'IEA_LB_RWT-AeroAcoustics'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('before', type=Path)
    p.add_argument('after', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--repeats', type=int, default=5)
    p.add_argument('--baseline-commit', required=True)
    args = p.parse_args()
    assert args.repeats >= 3
    args.output.mkdir(parents=True, exist_ok=False)
    executables = {'before': args.before.resolve(), 'after': args.after.resolve()}
    result = dict(baseline_commit=args.baseline_commit, platform=platform.platform(),
                  method='Alternating paired runs, one warm-up per executable, median internal wall time; includes output I/O',
                  executable_sha256={k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in executables.items()},
                  cases={})
    for name in ['official8', 'wind9', 'table_bl2', 'tno_short']:
        folder = args.output / name
        inputs = folder / 'inputs'
        shutil.copytree(ROOT / 'examples' / CASE, inputs)
        if name == 'wind9':
            key, value = 'HWindSpeed', '9'
        elif name == 'table_bl2':
            key, value = 'BLMod', '2'
        elif name == 'tno_short':
            key, value = 'TBLTEMod', '2'
        else:
            key = None
        if key:
            replacements = 0
            for file in inputs.rglob('*'):
                if file.is_file():
                    text = file.read_text(encoding='utf-8-sig')
                    text, n = re.subn(r'(?m)^(\s*)\S+(\s+' + key + r'\b)',
                                      lambda m: m[1] + value + m[2], text)
                    if n:
                        file.write_text(text, encoding='utf-8')
                        replacements += n
            assert replacements == 1, (key, replacements)
        outputs = {k: folder / k for k in executables}
        timings = {k: [] for k in executables}

        def run(k):
            command = [str(executables[k]), str((inputs / (CASE + '.fst')).resolve()),
                       str(outputs[k].resolve())]
            if name == 'tno_short':
                command += ['0.1']
            subprocess.run(command, check=True, capture_output=True, text=True)
            meta = json.loads((outputs[k] / 'run.json').read_text())
            assert meta['steps'] == (16 if name == 'tno_short' else 3200)
            assert meta['acoustic_samples'] == (2 if name == 'tno_short' else 201)
            return meta['elapsed_seconds']

        for k in executables:
            run(k)
        repeats = args.repeats if name != 'tno_short' else 0
        for i in range(repeats):
            for k in (['before', 'after'] if i % 2 == 0 else ['after', 'before']):
                timings[k].append(run(k))
        files = ['dynamics.csv', *[CASE + '_' + str(i) + ext for i in range(1, 5) for ext in ['.out', '.mask']]]
        hashes = {}
        for file in files:
            old, new = outputs['before'] / file, outputs['after'] / file
            assert old.read_bytes() == new.read_bytes(), (name, file, 'Numerical output changed')
            hashes[file] = hashlib.sha256(new.read_bytes()).hexdigest()
        entry = dict(byte_identical_outputs_sha256=hashes, seconds=timings)
        if repeats:
            before, after = [statistics.median(timings[k]) for k in ['before', 'after']]
            entry.update(median_before_s=before, median_after_s=after, speedup=before / after,
                         time_reduction_percent=100 * (1 - after / before))
        result['cases'][name] = entry
        print(name, json.dumps({k: v for k, v in entry.items() if k != 'byte_identical_outputs_sha256'}), flush=True)
        (args.output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
