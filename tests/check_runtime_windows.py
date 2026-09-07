"""Audit PE imports and run the turbine with only C++ runtime DLLs on Windows."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from compare_full_case import CASE, run

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('--gcc-bin', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    args = p.parse_args()
    assert os.name == 'nt', 'This PE/runtime audit is Windows-specific'
    # A new directory prevents a DLL left by a previous run from hiding a dependency.
    args.work.mkdir(parents=True, exist_ok=False)
    binary = args.work / args.executable.name
    shutil.copyfile(args.executable, binary)
    runtime_names = {'libgcc_s_seh-1.dll', 'libstdc++-6.dll', 'libwinpthread-1.dll'}
    imports = {}
    pending = [binary]
    while pending:
        path = pending.pop()
        data = subprocess.check_output([str(args.gcc_bin / 'objdump.exe'), '-p', str(path)], text=True)
        names = re.findall(r'DLL Name:\s*(\S+)', data)
        imports[path.name] = names
        for name in names:
            normalized = name.lower()
            assert not any(x in normalized for x in ['fortran', 'openfast', 'python', 'mkl']), name
            if normalized in runtime_names:
                target = args.work / name
                if not target.exists():
                    shutil.copyfile(args.gcc_bin / name, target)
                    pending.append(target)
            else:
                assert (Path(os.environ['SystemRoot']) / 'System32' / name).exists() or normalized.startswith('api-ms-win-'), name
    case = args.work / 'inputs'
    shutil.copytree(ROOT / 'examples' / CASE, case)
    environment = os.environ.copy()
    environment['PATH'] = str(Path(os.environ['SystemRoot']) / 'System32')
    for name in ['PYTHONPATH', 'PYTHONHOME', 'AEROACOUSTICS_RUNTIME_DIRS', 'MKLROOT', 'MKL_DIR']:
        environment.pop(name, None)
    output = args.work / 'outputs'
    subprocess.run([str(binary.resolve()), str(case.resolve() / (CASE + '.fst')), str(output.resolve())],
                   cwd=args.work.resolve(), env=environment, check=True)
    comparison = run(ROOT / 'reference/results' / CASE, output)
    result = dict(platform='Windows x86_64', executable_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  imports=imports, child_path=environment['PATH'], duration_s=20,
                  steps=3200, acoustic_comparison=comparison,
                  description='Only executable, C++ runtime DLLs and input files staged; reference outputs used after execution.')
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
