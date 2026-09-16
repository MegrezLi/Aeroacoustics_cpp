"""Compare pre-refactor and current C++ files; no physics implementation in Python."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASE = 'IEA_LB_RWT-AeroAcoustics'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('before', type=Path)
    p.add_argument('after', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--baseline-commit', default='2d7fda532b389371a854c49daad9e96282c938b1')
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    executables = {'before': args.before.resolve(), 'after': args.after.resolve()}
    result = {'baseline_commit': args.baseline_commit,
              'executables_sha256': {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in executables.items()},
              'cases': {}}
    for name, key, value, duration in [('official8', None, None, 20), ('wind9', 'HWindSpeed', '9', 20),
                                       ('table_bl2', 'BLMod', '2', 20), ('tno', 'TBLTEMod', '2', .1)]:
        folder = args.output / name
        inputs = folder / 'inputs'
        shutil.copytree(ROOT / 'examples' / CASE, inputs)
        if key:
            count = 0
            for f in inputs.rglob('*'):
                if f.is_file():
                    text, n = re.subn(r'(?m)^(\s*)\S+(\s+' + key + r'\b)',
                                      lambda m: m[1] + value + m[2], f.read_text(encoding='utf-8-sig'))
                    if n:
                        f.write_text(text, encoding='utf-8')
                        count += n
            assert count == 1
        for version, exe in executables.items():
            subprocess.run([str(exe), str((inputs / (CASE + '.fst')).resolve()),
                            str((folder / version).resolve()), str(duration)], check=True,
                           capture_output=True, text=True)
        hashes = {}
        for f in (folder / 'before').iterdir():
            if f.name == 'run.json':
                a, b = [json.loads((folder / version / f.name).read_text()) for version in executables]
                a.pop('elapsed_seconds'); b.pop('elapsed_seconds')
                # Newly added diagnostic counters do not change numerical output.
                a = {k: v for k, v in a.items() if not k.startswith('structural_')}
                b = {k: v for k, v in b.items() if not k.startswith('structural_')}
                assert a == b, (name, 'runtime metadata')
                continue
            before, after = f.read_bytes(), (folder / 'after' / f.name).read_bytes()
            assert before == after, (name, f.name)
            hashes[f.name] = hashlib.sha256(after).hexdigest()
        result['cases'][name] = {'duration_s': duration, 'byte_identical_sha256': hashes}
        print(name, 'all numerical files and lookup report byte-identical', flush=True)
        (args.output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
