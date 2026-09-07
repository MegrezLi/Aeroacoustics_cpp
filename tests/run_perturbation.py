"""Run a 9 m/s variant and compare only original Fortran with C++."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np
from compare_full_case import CASE, acoustic, run

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('work', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    case = args.work / 'inputs'
    shutil.copytree(ROOT / 'examples' / CASE, case, dirs_exist_ok=True)
    wind = case / 'RotorSE_FAST_IEA_landBased_RWT_InflowFile.dat'
    # Preserve line endings, so the changed input is identical to the Fortran run.
    content, count = re.subn(rb'(?m)^[ \t]*8[ \t]+(HWindSpeed\b)', rb'9            \1', wind.read_bytes())
    assert count == 1
    wind.write_bytes(content.replace(b'\r\n', b'\n').replace(b'\n', b'\r\n'))
    baseline = ROOT / 'reference/results' / (CASE + '-wind9')
    metadata = json.loads((baseline / 'reference.json').read_text(encoding='utf-8'))
    for folder, hashes in [(case, metadata['inputs_sha256']), (baseline, metadata['output_sha256'])]:
        for name, digest in hashes.items():
            assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == digest, name
    output = args.work / 'outputs'
    subprocess.run([str(args.executable.resolve()), str(case.resolve() / (CASE + '.fst')),
                    str(output.resolve())], check=True)
    result = run(baseline, output)
    original = acoustic(ROOT / 'reference/results' / CASE / (CASE + '_1.out'))[0]
    changed = acoustic(output / (CASE + '_1.out'))[0]
    response = float(np.max(np.abs(original[:, 1:] - changed[:, 1:])))
    assert response > 0.1, 'Wind-speed change did not affect the computed sound levels'
    result.update(wind_speed_m_s=9, changed_input='HWindSpeed: 8 -> 9',
                  max_oaspl_change_from_8m_s_dB=response, identical_input_files=len(metadata['inputs_sha256']),
                  cpp_executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest())
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
