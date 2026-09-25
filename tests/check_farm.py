"""E6 orchestration: all wake, source, turbine and statistical calculations run in C++."""
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from check_performance import ROOT, CASE, variant


def rows(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    build, out = args.build.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    suffix = '.exe' if (build / 'aeroacoustics_farm.exe').exists() else ''
    exe, single = build / ('aeroacoustics_farm' + suffix), build / ('aeroacoustics_turbine' + suffix)
    case = ROOT / 'examples' / CASE / (CASE + '.fst')
    fixtures = ROOT / 'examples/farm'
    default_units = [('upstream', 0, 0), ('downstream', 800, 0)]

    def run(name, units=default_units, wakes=True, tower=False, sources=(), shift=(0, 0),
            input_case=case, memory=32000000, propagation='none', controller='none', surfaces='none', fail=False):
        folder = out / (name + '-inputs')
        folder.mkdir()
        (folder / 'observers.dat').write_text(f'2\n! global receivers\n{300+shift[0]} {400+shift[1]} 2\n{1100+shift[0]} {400+shift[1]} 2\n')
        files = []
        for i, (label, x, y) in enumerate(units):
            path = folder / f'unit{i}.dat'
            tower_path = (fixtures / 'tower-test.dat').as_posix() if tower else 'none'
            path.write_text(f'"{label}" Name\n"{input_case.as_posix()}" Case\n{x+shift[0]} X\n{y+shift[1]} Y\n"{controller}" Controller\n"{tower_path}" Tower\n"{surfaces}" Surfaces\n')
            files.append(path.name)
        config = f'"Synthetic E6 regression" Provenance\n2 Duration\n0.2 StatisticsStart\n{str(wakes).lower()} Wakes\n0.05 WakeExpansion\n0.95 MaxCt\n{memory} MaxValues\n"observers.dat" ObserverFile\n"{propagation}" Propagation\n{len(files)} NumTurbines\n'
        config += '\n'.join(f'"{f}"' + (' TurbineFiles' if i == 0 else '') for i, f in enumerate(files))
        config += f'\n{len(sources)} NumSources\n'
        config += '\n'.join('"' + str(p).replace('\\', '/') + '"' + (' SourceFiles' if i == 0 else '') for i, p in enumerate(sources)) + '\n'
        (folder / 'farm.dat').write_text(config)
        p = subprocess.run([str(exe), str(folder / 'farm.dat'), str(out / name)], capture_output=True, text=True)
        assert (p.returncode != 0) == fail, (name, p.stdout, p.stderr)
        if fail:
            assert not (out / name / 'farm.json').exists()
            return p.stderr
        meta = json.loads((out / name / 'farm.json').read_text())
        assert not meta['field_validated']
        return out / name

    base = run('single', units=default_units[:1], wakes=False)
    # Single-turbine degeneration uses exactly the same receiver map and time window.
    template = (ROOT / 'examples/engineering/metrics.dat').read_text()
    observer_path = (out / 'single-inputs/observers.dat').as_posix()
    for key, value in {'Start': .2, 'End': 2, 'RetardedTime': 'false', 'NumTones': 0,
                       'ObserverFile': '"' + observer_path + '"'}.items():
        template, count = re.subn(r'(?m)^\S+(\s+' + key + r'\b)', lambda m: str(value) + m[1], template)
        assert count == 1, key
    (out / 'metrics.dat').write_text(template)
    p = subprocess.run([str(single), str(case), str(out / 'standalone'), '2', '--metrics=' + str(out / 'metrics.dat')], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    for file in (base / 'turbine_1').glob('*.out'):
        assert file.read_bytes() == (out / 'standalone' / file.name).read_bytes(), file.name
    for actual, expected in zip(rows(base / 'farm_receivers.csv'), rows(out / 'standalone/receiver_map.csv')):
        assert abs(float(actual['LAeq_dBA']) - float(expected['LAeq_turbine_dB'])) < 1e-10
    free = run('free', wakes=False)
    wake = run('wake')
    diagnostics = rows(wake / 'wake_diagnostics.csv')
    upstream = [float(r['hub_inflow_m_s']) for r in diagnostics if r['turbine'] == 'upstream']
    downstream = [float(r['hub_inflow_m_s']) for r in diagnostics if r['turbine'] == 'downstream']
    assert all(b < a for a, b in zip(upstream, downstream))
    assert (wake / 'turbine_1/dynamics.csv').read_bytes() == (free / 'turbine_1/dynamics.csv').read_bytes()
    assert (wake / 'turbine_2/dynamics.csv').read_bytes() != (free / 'turbine_2/dynamics.csv').read_bytes()
    reverse = run('reverse', units=list(reversed(default_units)))
    assert (wake / 'farm_history.csv').read_bytes() == (reverse / 'farm_history.csv').read_bytes()
    translated = run('translated', shift=(10000, -2000))
    translation_error = max(abs(float(a['LA_dBA']) - float(b['LA_dBA']))
                            for a, b in zip(rows(wake / 'farm_history.csv'), rows(translated / 'farm_history.csv')))
    assert translation_error < 1e-6
    tower = run('tower', units=default_units[:1], tower=True, wakes=False)
    assert (tower / 'turbine_1/dynamics.csv').read_bytes() != (base / 'turbine_1/dynamics.csv').read_bytes()
    assert json.loads((tower / 'turbine_1/run.json').read_text())['tower_influence']['potential']
    both = run('sources', sources=[fixtures / 'cooling-test.dat', fixtures / 'tone-test.dat'])
    # Incoherent energy totals must match the separately averaged source contributions.
    sums = {}
    for r in rows(both / 'source_contributions.csv'):
        sums[r['observer']] = sums.get(r['observer'], 0) + 10 ** (float(r['LAeq_dBA']) / 10)
    for r in rows(both / 'farm_receivers.csv'):
        assert abs(10 * math.log10(sums[r['observer']]) - float(r['LAeq_dBA'])) < 1e-10
    assert (both / 'turbine_2/dynamics.csv').read_bytes() == (wake / 'turbine_2/dynamics.csv').read_bytes()
    weighted_case = variant(out / 'weighted-case', {'AWeighting': 'True'})
    weighted = run('weighted', input_case=weighted_case, sources=[fixtures / 'cooling-test.dat', fixtures / 'tone-test.dat'])
    weighting_error = max(abs(float(a['LAeq_dBA'])-float(b['LAeq_dBA'])) for a,b in zip(rows(both/'farm_receivers.csv'),rows(weighted/'farm_receivers.csv')))
    assert weighting_error < 1e-8
    air_text = (ROOT / 'examples/engineering/propagation.dat').read_text().replace('rigid GroundModel', 'none GroundModel')
    (out / 'air.dat').write_text(air_text)
    air = run('air', propagation=(out/'air.dat').as_posix(), sources=[fixtures / 'cooling-test.dat'])
    assert all(math.isfinite(float(r['LAeq_dBA'])) for r in rows(air/'farm_receivers.csv'))
    combined = run('combined', tower=True, controller=(ROOT/'examples/engineering/controller.dat').as_posix(),
                   surfaces=(ROOT/'examples/engineering/surfaces.dat').as_posix(),
                   propagation=(out/'air.dat').as_posix(), sources=[fixtures/'cooling-test.dat'])
    assert (combined/'turbine_2/operation.csv').exists()
    run('memory', memory=10, fail=True)
    run('close', units=[('a', 0, 0), ('b', 150, 0)], fail=True)
    run('duplicate', units=[('a', 0, 0), ('a', 800, 0)], fail=True)
    run('ground', propagation=(ROOT / 'examples/engineering/propagation.dat').as_posix(), fail=True)
    report = {'single_turbine': 'original per-turbine outputs byte-identical; LAeq agrees with receiver metrics',
              'wake_feedback': 'downstream wind and dynamics change; upstream dynamics unchanged',
              'downstream_hub_wind_range_m_s': [min(downstream), max(downstream)],
              'ct_limited_samples': json.loads((wake / 'farm.json').read_text())['ct_limited_samples'],
              'order_independence': 'farm history byte-identical after reversing input order',
              'translation_max_abs_dB': translation_error,
              'a_weighting_max_abs_dB': weighting_error,
              'combined': 'two controlled turbines with surface data, tower influence, absorption and supplied cooling spectrum',
              'tower': 'leaf inflow/dynamics respond; source metadata preserved',
              'source_energy': 'separate LAeq contributions reproduce aggregate LAeq in energy',
              'guards': 'memory, near-wake spacing, duplicate names and unsupported ground rejected',
              'field_validated': False, 'executable_sha256': hashlib.sha256(exe.read_bytes()).hexdigest()}
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
