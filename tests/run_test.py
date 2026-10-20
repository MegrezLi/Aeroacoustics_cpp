"""Give each CTest invocation its own directory; preserve reports and failed fixtures."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    root, name, *command = sys.argv[1:]
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=name + '-', dir=root))
    command = [arg.replace('@WORK@', str(work)) for arg in command]
    env = dict(os.environ, AEROACOUSTICS_TEST_WORKDIR=str(work))
    print(f'Test artifacts: {work}', flush=True)
    status = subprocess.run(command, cwd=work, env=env).returncode
    if status:
        print(f'Command exited with {status} (0x{status & 0xffffffff:08x})', flush=True)
    return status


if __name__ == '__main__':
    sys.exit(main())
