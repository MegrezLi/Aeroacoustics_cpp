"""Compare full coupled trajectories from original Fortran and C++ probes."""
import argparse
import json
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('fortran_csv', type=Path)
    p.add_argument('cpp_csv', type=Path)
    p.add_argument('--report', type=Path, required=True)
    args = p.parse_args()
    with args.fortran_csv.open(encoding='utf-8') as stream:
        names = stream.readline().strip().split(',')
    a = np.loadtxt(args.fortran_csv, delimiter=',', skiprows=1)
    b = np.loadtxt(args.cpp_csv, delimiter=',', skiprows=1)
    b = b[b[:, 0] > 0]  # Fortran's diagnostic begins after the first step.
    assert a.shape == b.shape and len(a) == 3200 * 90
    assert np.isfinite(a).all() and np.isfinite(b).all()
    np.testing.assert_allclose(a[:, :3], b[:, :3], rtol=0, atol=1e-12)
    # Full trajectories cross polar/UA piecewise boundaries. These bounds
    # are distinct from the much tighter equal-input module tests.
    bounds = dict(q1=2e-5, q2=2e-5, q3=1e-5, qd1=2e-4, qd2=3e-4, qd3=1e-4,
                  phi=5e-6, alpha=5e-6, speed=2e-4, a=3e-5, ap=1e-5,
                  cl=0.002, cd=0.001, cm=0.001, x=3e-5, y=3e-5, z=3e-5,
                  fx=5, fy=5, fz=5, mx=5, my=5, mz=5)
    result = dict(reference='Original Fortran solver with read-only trace', rows=len(a),
                  compared_values=int(a[:, 3:].size), duration_s=20, channels={})
    for j, name in enumerate(names[3:], 3):
        d = np.abs(a[:, j] - b[:, j])
        i = int(d.argmax())
        result['channels'][name] = dict(max_absolute_error=float(d[i]), tolerance=bounds[name],
                                        worst_time_s=float(a[i, 0]), blade=int(a[i, 1]) + 1,
                                        node=int(a[i, 2]) + 1)
        assert d[i] <= bounds[name], f'{name}: {d[i]} > {bounds[name]}'
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
