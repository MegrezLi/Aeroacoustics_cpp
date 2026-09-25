"""E8 test orchestration only: C++ computes comparisons, sensitivities and uncertainty."""
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


def write(path, records):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    build, out = args.build.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    suffix = '.exe' if (build / 'aeroacoustics_validate.exe').exists() else ''
    exe = build / ('aeroacoustics_validate' + suffix)
    turbine = build / ('aeroacoustics_turbine' + suffix)
    fixtures = ROOT / 'examples/validation'

    def run(*arguments, fail=False):
        p = subprocess.run([str(exe), *map(str, arguments)], capture_output=True, text=True)
        assert (p.returncode != 0) == fail, (arguments, p.stdout, p.stderr)
        return p

    obs, pred = fixtures / 'observations.csv', fixtures / 'predictions.csv'
    run('compare', obs, pred, out / 'analytic', 2)
    results = rows(out / 'analytic/groups.csv')
    held = next(r for r in results if r['split'] == 'validation')
    assert float(held['bias_db']) == 0 and float(held['rmse_db']) == 2 and held['count'] == '2'
    report = json.loads((out / 'analytic/report.json').read_text())
    assert report['declared_external_validation_rows'] == 0 and not report['field_accuracy_established']
    # Calibration residual is 10 dB; it must not contaminate the held-out score or fit a bias.
    assert next(r for r in results if r['split'] == 'calibration')['bias_db'] == '10'
    for name, changes in [('leakage', {'group': 'training_fixture'}), ('wrong_window', {'end_s': '19'}),
                          ('negative_u', {'standard_u_db': '-1'}), ('nan', {'level_db': 'nan'}),
                          ('wrong_weighting', {'weighting': 'Z'}), ('unknown_signal', {'signal': 'unknown'}),
                          ('duplicate', {'id': 'cal1'})]:
        records = rows(obs)
        records[1].update(changes)
        path = out / (name + '.csv')
        write(path, records)
        run('compare', path, pred, out / name, 2, fail=True)
        assert not (out / name / 'report.json').exists()
    run('compare', obs, pred, out / 'analytic', 2, fail=True)
    run('compare', obs, pred, out / 'bad_k', '2junk', fail=True)
    # Quoted commas, quotes and UTF-8 provenance survive the copied input audit trail.
    records = rows(obs)
    records[0]['provenance'] = 'Synthetic "quoted", 合成测试'
    write(out / 'quoted.csv', records)
    run('compare', out / 'quoted.csv', pred, out / 'quoted-result', 2)
    assert rows(out / 'quoted-result/observations.csv')[0]['provenance'] == records[0]['provenance']
    records[1]['standard_u_db'] = ''
    write(out / 'unknown-u.csv', records)
    run('compare', out / 'unknown-u.csv', pred, out / 'unknown-u', 2)
    unknown = next(r for r in rows(out / 'unknown-u/residuals.csv') if r['id'] == 'val1')
    assert unknown['combined_standard_u_db'] == '' and unknown['normalized_error'] == ''
    band_obs, band_pred = rows(obs), rows(pred)
    for records in (band_obs, band_pred):
        for r in records:
            r.update(quantity='band_Leq', weighting='Z', frequency_hz='1000', bandwidth_hz='230')
    write(out / 'band-obs.csv', band_obs)
    write(out / 'band-pred.csv', band_pred)
    run('compare', out / 'band-obs.csv', out / 'band-pred.csv', out / 'bands', 2)
    band_pred[1]['bandwidth_hz'] = '200'
    write(out / 'band-pred-mismatch.csv', band_pred)
    run('compare', out / 'band-obs.csv', out / 'band-pred-mismatch.csv', out / 'wrong-band', 2, fail=True)
    run('budget', fixtures / 'perturbations.csv', fixtures / 'correlations.csv', out / 'budget', 2)
    assert math.isclose(float(rows(out / 'budget/budgets.csv')[0]['standard_u_db']), math.sqrt(13), rel_tol=1e-12)
    corr = rows(fixtures / 'correlations.csv')
    corr[0]['parameter_b'] = 'missing'
    write(out / 'bad-corr.csv', corr)
    run('budget', fixtures / 'perturbations.csv', out / 'bad-corr.csv', out / 'bad-budget', 2, fail=True)
    run('ensemble', fixtures / 'samples.csv', out / 'ensemble')
    ensemble = rows(out / 'ensemble/ensemble.csv')[0]
    assert float(ensemble['p025_db']) == .5 and float(ensemble['p975_db']) == 19.5
    samples = rows(fixtures / 'samples.csv')
    samples.append(dict(samples[0], output_id='incomplete'))
    write(out / 'incomplete.csv', samples)
    run('ensemble', out / 'incomplete.csv', out / 'bad-ensemble', fail=True)

    # Five actual C++ turbine runs: two symmetric wind perturbation sizes, with identical
    # receiver/window/surface setup. This is a sensitivity test, never a field measurement.
    template = (ROOT / 'examples/engineering/metrics.dat').read_text()
    for key, value in {'Start': 0, 'End': 2, 'RetardedTime': 'false', 'NumTones': 0, 'ObserverFile': 'none'}.items():
        template, count = re.subn(r'(?m)^\S+(\s+' + key + r'\b)', lambda m: str(value) + m[1], template)
        assert count == 1, key
    (out / 'metrics.dat').write_text(template)
    maps = {}
    for wind in (7.8, 7.9, 8., 8.1, 8.2):
        name = f'wind{wind:g}'
        case = variant(out / (name + '-inputs'), {'HWindSpeed': wind})
        dest = out / name
        p = subprocess.run([str(turbine), str(case), str(dest), '2', '--metrics=' + str(out / 'metrics.dat')],
                           capture_output=True, text=True)
        assert p.returncode == 0, p.stderr
        maps[wind] = rows(dest / 'receiver_map.csv')
    context = []
    for r in maps[8.]:
        item = dict(rows(pred)[1])
        item.update(id='receiver' + r['observer'], group='turbine_synthetic', origin='simulation',
                    receiver=r['observer'], start_s=r['start_receiver_s'], end_s=r['end_receiver_s'],
                    x_m=r['x_m'], y_m=r['y_m'], z_m=r['z_m'], level_db='0', surface='original_case')
        context.append(item)
    write(out / 'context.csv', context)
    run('import-map', out / 'wind8/receiver_map.csv', out / 'context.csv', out / 'imported.csv',
        'C++ actual run wind8/run.json; original_case; fixed 8 m/s; 0 deg')
    imported = rows(out / 'imported.csv')
    for r, m in zip(imported, maps[8.]):
        assert float(r['level_db']) == float(m['LAeq_turbine_dB'])
    # Exact synthetic measurement fixture checks the adapter and zero residual limit.
    write(out / 'synthetic-observations.csv', [dict(r, origin='synthetic', provenance='C++ output copied for interface test only') for r in imported])
    run('compare', out / 'synthetic-observations.csv', out / 'imported.csv', out / 'turbine-comparison', 2)
    assert all(float(r['rmse_db']) == 0 for r in rows(out / 'turbine-comparison/groups.csv'))
    context[0]['x_m'] = str(float(context[0]['x_m']) + 1)
    write(out / 'bad-context.csv', context)
    run('import-map', out / 'wind8/receiver_map.csv', out / 'bad-context.csv', out / 'bad-import.csv', 'test', fail=True)
    (out / 'independent.csv').write_text('output_id,parameter_a,parameter_b,rho\n')
    perturbations = []
    for step in (.2, .1):
        for i, r in enumerate(maps[8.]):
            perturbations.append(dict(output_id=f'receiver{r["observer"]}_step{step}', parameter='HWindSpeed', unit='m/s',
                x_minus=8-step, x0=8, x_plus=8+step, u_x=.1,
                y_minus_db=maps[round(8-step, 1)][i]['LAeq_turbine_dB'], y0_db=r['LAeq_turbine_dB'],
                y_plus_db=maps[round(8+step, 1)][i]['LAeq_turbine_dB'],
                provenance='Actual C++ turbine runs, 2s; assumed u_wind=0.1 m/s for test only'))
    write(out / 'turbine-perturbations.csv', perturbations)
    run('budget', out / 'turbine-perturbations.csv', out / 'independent.csv', out / 'turbine-budget', 2)
    sensitivities = rows(out / 'turbine-budget/sensitivities.csv')
    assert all(math.isfinite(float(r['derivative_db_per_unit'])) and abs(float(r['derivative_db_per_unit'])) > .01 for r in sensitivities)
    report = {'analytic_comparison': 'bias=0, RMSE=2 dB; calibration excluded',
              'correlated_budget': 'u=sqrt(13), singular and invalid matrices checked by C++ probe',
              'ensemble': 'equal-weight 2.5/50/97.5 percentiles; incomplete outputs rejected',
              'turbine_wind_sensitivities': sensitivities,
              'guards': 'split leakage, context, uncertainty, IDs, CSV, output reuse and missing draws',
              'field_accuracy_established': False,
              'executable_sha256': hashlib.sha256(exe.read_bytes()).hexdigest()}
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
