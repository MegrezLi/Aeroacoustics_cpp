"""Exercise the literature example and its analysis script together."""
from pathlib import Path
import subprocess
import sys

executable, output = sys.argv[1:]
output = Path(output)
csv = output / 'spectrum.csv'
subprocess.run([executable, str(csv), '--check'], check=True)
subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / 'examples/bpm1989/analyze.py'),
                str(csv), '--output-dir', str(output / 'analysis')], check=True)
