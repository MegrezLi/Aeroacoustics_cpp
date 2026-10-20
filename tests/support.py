"""Shared input fixtures for CLI regressions; no acoustic formulas live here."""
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
CASE = 'IEA_LB_RWT-AeroAcoustics'


def variant(destination, changes):
    destination = Path(destination)
    shutil.copytree(ROOT / 'examples' / CASE, destination)
    if changes:
        counts = dict.fromkeys(changes, 0)
        labels = '|'.join(re.escape(key) for key in changes)
        pattern = re.compile(r'(?m)^([ \t]*)\S+([ \t]+)(' + labels + r')\b')

        def replace(match):
            key = match[3]
            counts[key] += 1
            return match[1] + str(changes[key]) + match[2] + key

        for path in destination.rglob('*'):
            if path.is_file():
                text, count = pattern.subn(replace, path.read_text(encoding='utf-8-sig'))
                if count:
                    path.write_text(text, encoding='utf-8')
        if any(count != 1 for count in counts.values()):
            raise ValueError(f'Expected exactly one input value per label: {counts}')
    return destination / (CASE + '.fst')
