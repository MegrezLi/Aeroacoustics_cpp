"""Paired timings and allocation counts for the five October optimizations."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('before', type=Path)
    parser.add_argument('after', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--rounds', type=int, default=5)
    args = parser.parse_args()
    if args.rounds < 1:
        parser.error('--rounds must be positive')
    executables = {'before': args.before.resolve(), 'after': args.after.resolve()}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {
        'baseline': 'ffc85ed (four correctness fixes included)',
        'rounds': args.rounds,
        'executables_sha256': {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in executables.items()},
        'raw': {'before': [], 'after': []},
        'scope': 'Microbenchmarks; cumulative requested bytes, not peak RSS. Timings exclude compilation.',
    }
    for trial in range(args.rounds + 1):
        rows = {}
        for version in (['before', 'after'] if trial % 2 == 0 else ['after', 'before']):
            directory = output / f'{version}-{trial}'
            result = subprocess.run([str(executables[version]), str(directory)], check=True,
                                    capture_output=True, text=True)
            (output / f'{version}-{trial}.csv').write_text(result.stdout)
            rows[version] = list(csv.DictReader(io.StringIO(result.stdout)))
            if trial:  # Whole-process warm-up excluded from reported rounds.
                report['raw'][version].append(rows[version])
        for before, after in zip(rows['before'], rows['after']):
            if before['case'] != after['case']:
                raise AssertionError('Benchmark case order differs')
            if abs(float(before['checksum']) - float(after['checksum'])) > 1e-7:
                raise AssertionError('Benchmark numerical checksum differs: ' + before['case'])
        # Receiver history, maps, bins, metadata and AM must retain exact output.
        for reference in (output / f'before-{trial}').iterdir():
            candidate = output / f'after-{trial}' / reference.name
            if reference.read_bytes() != candidate.read_bytes():
                raise AssertionError(f'Metrics output differs: {reference.name}')
    report['metrics_files_byte_identical'] = True
    summary = {}
    for index, row in enumerate(report['raw']['before'][0]):
        name = row['case']
        values = {}
        for version in executables:
            samples = [trial[index] for trial in report['raw'][version]]
            values[version] = {
                'median_seconds': statistics.median(float(s['seconds']) for s in samples),
                'allocations': int(samples[0]['allocations']),
                'requested_bytes': int(samples[0]['requested_bytes']),
            }
        values['elapsed_reduction_percent'] = 100 * (1 - values['after']['median_seconds'] / values['before']['median_seconds'])
        summary[name] = values
    report['summary'] = summary
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
