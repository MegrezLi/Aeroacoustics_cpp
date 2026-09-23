"""E3/E5/E7 CLI regression orchestration. All simulation and metrics are computed in C++."""
import argparse
import csv
import json
import re
import subprocess
from pathlib import Path
import numpy as np
from check_performance import ROOT, CASE, variant


def rows(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    exe, out = args.executable.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    fixtures = ROOT / 'examples/engineering'
    case = ROOT / 'examples' / CASE / (CASE + '.fst')
    report = {}

    def config(name, **changes):
        text = (fixtures / 'metrics.dat').read_text()
        text = text.replace('"receiver-map.dat"', '"' + (fixtures / 'receiver-map.dat').as_posix() + '"')
        text = text.replace('"tone-test.dat"', '"' + (fixtures / 'tone-test.dat').as_posix() + '"')
        for key, value in changes.items():
            text, count = re.subn(r'(?m)^\S+(\s+' + key + r'\b)', lambda m: str(value) + m[1], text)
            assert count == 1, key
        path = out / (name + '.dat')
        path.write_text(text)
        return '--metrics=' + str(path)

    def execute(name, input_case=case, duration=20, options=(), fails=False):
        destination = out / name
        process = subprocess.run([str(exe), str(input_case), str(destination), str(duration), *options],
                                 capture_output=True, text=True)
        if fails:
            assert process.returncode != 0, name
            if (destination / 'run.json').exists():
                assert (destination / 'run.json').read_text() == '', name
            return process.stderr
        assert process.returncode == 0, process.stderr
        json.loads((destination / 'run.json').read_text())
        return destination

    surface = '--surfaces=' + str(fixtures / 'surfaces.dat')
    base = execute('baseline', duration=2)
    changed = execute('surface', duration=2, options=[surface])
    assert (base / 'dynamics.csv').read_bytes() == (changed / 'dynamics.csv').read_bytes()
    assert (base / (CASE + '_1.out')).read_bytes() != (changed / (CASE + '_1.out')).read_bytes()
    metadata = rows(changed / 'surface_datasets.csv')[0]
    assert metadata['polar_replaced'] == '0' and 'Synthetic' in metadata['provenance']
    report['surface_sensitivity'] = 'BL-only input changes noise and preserves dynamics byte-for-byte'
    polar_file = ROOT / 'examples' / CASE / 'Airfoils/RotorSE_FAST_IEA_landBased_RWT_AeroDyn_Polar_29.dat'
    polar_state = (fixtures / 'surface-test.dat').read_text().replace('"surface-test-bl.dat"',
                  '"' + (fixtures / 'surface-test-bl.dat').as_posix() + '"').replace('"none" Polar', '"' + polar_file.as_posix() + '" Polar')
    (out / 'polar-state.dat').write_text(polar_state)
    (out / 'polar-states.dat').write_text('1 NumStates\n"polar-state.dat" StateFiles\n')
    supplied = execute('supplied-polar', duration=2, options=['--surfaces=' + str(out / 'polar-states.dat')])
    assert rows(supplied / 'surface_datasets.csv')[0]['polar_replaced'] == '1'
    assert (supplied / 'dynamics.csv').read_bytes() == (changed / 'dynamics.csv').read_bytes()
    assert (supplied / (CASE + '_1.out')).read_bytes() == (changed / (CASE + '_1.out')).read_bytes()
    report['full_polar_input'] = 'supplied original polar/UA/coordinates reproduces BL-only state exactly'

    standard = config('standard', RetardedTime='false', NumTones=0, ObserverFile='none', Start=0, End=2)
    compared = execute('metrics-default', duration=2, options=[standard])
    for file in base.glob('*.out'):
        assert file.read_bytes() == (compared / file.name).read_bytes()
    # Statistics work even if the ordinary files request only total sound level.
    minimal_case = variant(out / 'minimal-inputs', {'NrOutFile': 1})
    minimal = execute('minimal', input_case=minimal_case, duration=2, options=[standard])
    assert (minimal / 'receiver_history.csv').read_bytes() == (compared / 'receiver_history.csv').read_bytes()

    unweighted = execute('unweighted', options=[config('unweighted', NumTones=0, RetardedTime='false', BackgroundLAeq=30)])
    weighted_case = variant(out / 'weighted-inputs', {'AWeighting': 'True'})
    weighted = execute('weighted', input_case=weighted_case, options=[config('weighted', NumTones=0, RetardedTime='false', BackgroundLAeq=30)])
    for a, b in zip(rows(unweighted / 'receiver_map.csv'), rows(weighted / 'receiver_map.csv')):
        for key in ('LAeq_turbine_dB', 'L5_dBA', 'L50_dBA', 'L95_dBA'):
            assert abs(float(a[key]) - float(b[key])) < 1e-9, key
        expected = 10 * np.log10(10 ** (float(a['LAeq_turbine_dB']) / 10) + 1000)
        assert abs(float(a['LAeq_with_background_dB']) - expected) < 1e-10
    report['weighting_background'] = 'A-weighted and unweighted source runs yield the same LAeq; background energy addition verified'

    result = execute('receiver-time', options=['--metrics=' + str(fixtures / 'metrics.dat'), surface])
    summary = rows(result / 'receiver_map.csv')
    assert len(summary) == 4
    assert len(rows(result / 'am_windows.csv')) == 4
    tones = rows(result / 'receiver_tones.csv')
    frequencies = [float(row['frequency_received_Hz']) for row in tones]
    assert min(frequencies) < 990 and max(frequencies) > 1010
    bins = rows(result / 'wind_bins.csv')
    for row in summary:
        subset = [b for b in bins if b['observer'] == row['observer']]
        duration = sum(float(b['duration_s']) for b in subset)
        energy = sum(float(b['duration_s']) * 10 ** (float(b['LAeq_turbine_dB']) / 10) for b in subset)
        assert abs(duration - float(row['duration_s'])) < 1e-9
        assert abs(10 * np.log10(energy / duration) - float(row['LAeq_turbine_dB'])) < 1e-9
    report['moving_tone_frequency_range_Hz'] = [min(frequencies), max(frequencies)]
    report['receiver_map'] = summary

    refined_case = variant(out / 'refined-inputs', {'DT_AA': .05})
    refined = execute('refined', input_case=refined_case, options=[config('refined', ReceiverDT=.05), surface])
    errors = [abs(float(a['LAeq_turbine_dB']) - float(b['LAeq_turbine_dB']))
              for a, b in zip(summary, rows(refined / 'receiver_map.csv'))]
    assert max(errors) < .15, errors
    report['DT_AA_0p1_vs_0p05_LAeq_max_abs_dB'] = max(errors)

    # Integrate with E1 closed-loop operation and E2 direct-path air absorption.
    propagation_file = out / 'air-only.dat'
    propagation_file.write_text((fixtures / 'propagation.dat').read_text().replace('rigid Ground', 'none Ground'))
    combined = execute('closed-loop', input_case=fixtures / 'closed_loop.fst',
                       options=[surface, config('closed-loop'), '--controller=' + str(fixtures / 'controller.dat'),
                                '--wind-grid=' + str(fixtures / 'gust_veer.wind'), '--propagation=' + str(propagation_file)])
    combined_bins = rows(combined / 'wind_bins.csv')
    assert len({row['wind_lower_m_s'] for row in combined_bins}) >= 2
    assert (combined / 'operation.csv').exists()
    report['closed_loop_wind_absorption'] = 'surface datasets, variable-speed tones, retarded metrics and wind bins completed together'

    apparent = execute('apparent', duration=2, options=[config('apparent', NumTones=0, RetardedTime='false', ApparentPower='true', Start=0, End=2)])
    assert all(row['apparent_LWA_freefield_dB'] for row in rows(apparent / 'receiver_map.csv'))
    failed = execute('duplicate-rejected', duration=2, options=[standard, '--metrics=duplicate'], fails=True)
    assert 'duplicate' in failed
    error = execute('ground-delay', duration=2, options=['--metrics=' + str(fixtures / 'metrics.dat'),
                    '--propagation=' + str(fixtures / 'propagation.dat')], fails=True)
    assert 'RetardedTime' in error
    execute('memory-limit', duration=2, options=[config('memory', MaxValues=1)], fails=True)
    execute('no-window', duration=2, options=[config('no-window', Start=10, End=12)], fails=True)

    bad_surface = (fixtures / 'surface-test.dat').read_text().replace('"surface-test-bl.dat"',
                   '"' + (fixtures / 'surface-test-bl.dat').as_posix() + '"').replace('1000 ReMin', '90000000 ReMin')
    (out / 'bad-surface.dat').write_text(bad_surface)
    (out / 'bad-surfaces.dat').write_text('1 NumStates\n"bad-surface.dat" StateFiles\n')
    error = execute('surface-range', duration=2, options=['--surfaces=' + str(out / 'bad-surfaces.dat')], fails=True)
    assert 'outside validated' in error
    report['guards'] = 'duplicate flags, incompatible ground delay, memory cap, empty window, surface validity interval rejected'
    report['field_validated'] = False
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
