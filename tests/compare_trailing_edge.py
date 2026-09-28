"""Compare C++ to original Fortran TNO and an independent Fortran Howe complex-step oracle."""
import argparse
import csv
import ctypes as ct
import json
import os
import subprocess
from pathlib import Path
import numpy as np
from fortran_reference import ROOT, build


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('build', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--fortran-compiler', required=True)
    args = p.parse_args()
    handles = []
    if os.name == 'nt':
        for directory in filter(None, os.environ.get('AEROACOUSTICS_RUNTIME_DIRS', '').split(os.pathsep)):
            handles.append(os.add_dll_directory(directory))
    folder = args.output.resolve()
    folder.mkdir(parents=True, exist_ok=False)
    build_dir = args.build.resolve()
    suffix = '.exe' if os.name == 'nt' else ''
    subprocess.run([str(build_dir / ('trailing_edge_probe' + suffix)), str(folder / 'rows'),
                    str(ROOT / 'examples/IEA_LB_RWT-AeroAcoustics/IEA_LB_RWT-AeroAcoustics.fst'),
                    str(ROOT / 'examples/trailing-edge')], check=True)
    lib, command = build(args.fortran_compiler, 'double')
    f = lib.__acoustic_reference_MOD_evaluate_tno
    pointer = ct.POINTER(ct.c_double)
    f.restype = None
    f.argtypes = [pointer] * 9 + [ct.POINTER(ct.c_int), pointer]
    groups = {}
    with (folder/'rows/tno.csv').open() as stream:
        for row in csv.DictReader(stream):
            key = (row['mode'], float(row['suction_ratio']), float(row['pressure_ratio']))
            groups.setdefault(key, []).append(row)
    maximum = 0
    values = 0
    for (mode, rs, rp), rows in groups.items():
        scalar = [ct.c_double(x) for x in (63.92, 90, 90, .509, 1.22)]
        cf = np.array([.000378576, .00198438])
        delta = np.array([.0110586, .00746583])
        ratios = np.array([1., 1.] if mode == 'reference' else [rs, rp])
        freq = np.array([float(r['frequency']) for r in rows])
        result = np.zeros((len(rows), 2), order='F')
        n = ct.c_int(len(rows))
        f(*[ct.byref(x) for x in scalar], cf.ctypes.data_as(pointer), delta.ctypes.data_as(pointer),
          ratios.ctypes.data_as(pointer), freq.ctypes.data_as(pointer), ct.byref(n), result.ctypes.data_as(pointer))
        actual = np.array([[float(r['pressure']), float(r['suction'])] for r in rows])
        np.testing.assert_allclose(actual, result, rtol=0, atol=2e-8)
        maximum = max(maximum, float(np.max(np.abs(actual-result))))
        values += actual.size
    path = folder / ('howe-oracle.dll' if os.name == 'nt' else 'howe-oracle.so')
    cmd = [args.fortran_compiler, '-shared', '-O0', '-fcheck=all']
    if os.name != 'nt': cmd += ['-fPIC']
    subprocess.run(cmd + [str(ROOT/'tests/howe_reference.f90'), '-o', str(path)], cwd=folder, check=True)
    oracle = ct.CDLL(str(path)).howe_reference_shape
    oracle.restype = ct.c_double
    oracle.argtypes = [ct.c_double] * 4
    error = 0
    count = 0
    with (folder/'rows/howe.csv').open() as stream:
        for row in csv.DictReader(stream):
            expected = oracle(*[float(row[k]) for k in ('q', 'h', 'lambda', 'delta')])
            actual = float(row['shape'])
            difference = abs(actual-expected)/max(actual, expected)
            assert difference < 1e-9, (row, expected, difference)
            error = max(error, difference)
            count += 1
    report = {'TNO': {'values': values, 'max_abs_error_dB': maximum,
                       'reference': 'untouched OpenFAST TBLTE_TNO body; new shim supplies both side velocity ratios',
                       'modes': 'reference/input; unit, unequal and signed input ratios'},
              'Howe': {'cases': count, 'max_relative_shape_error': error,
                       'reference': 'new independent Fortran complex-step differentiation of published f; not OpenFAST'},
              'physical_validation': 'not field or wind-tunnel validation'}
    (folder/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
