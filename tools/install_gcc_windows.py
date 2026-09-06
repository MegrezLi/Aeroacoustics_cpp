"""Install a pinned WinLibs toolchain alongside existing compilers (no PATH edits)."""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

TAG = '16.2.0posix-14.0.0-ucrt-r1'
NAME = 'winlibs-x86_64-posix-seh-gcc-16.2.0-mingw-w64ucrt-14.0.0-r1.zip'
URL = f'https://github.com/brechtsanders/winlibs_mingw/releases/download/{TAG}/{NAME}'
SHA256 = 'c1f52294597c0b73786b2a78eb5d176d89226d2f21875eab75e783a8b1cefcc4'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    cache = Path(__file__).resolve().parents[1] / '_work/downloads'
    cache.mkdir(parents=True, exist_ok=True)
    archive = cache / NAME
    if not archive.exists():
        partial = archive.with_suffix('.partial')
        print(f'Downloading {URL}', flush=True)
        urllib.request.urlretrieve(URL, partial)
        partial.rename(archive)
    actual = hashlib.file_digest(archive.open('rb'), 'sha256').hexdigest()
    if actual != SHA256:
        raise RuntimeError(f'SHA256 mismatch: {actual}')
    print(f'SHA256 verified: {actual}', flush=True)
    destination = args.destination.resolve()
    if destination.exists():
        raise RuntimeError(f'Refusing to overwrite existing installation: {destination}')
    with zipfile.ZipFile(archive) as package:
        for item in package.infolist():
            target = (destination / item.filename).resolve()
            if not target.is_relative_to(destination):
                raise RuntimeError(f'Unsafe archive entry: {item.filename}')
        print(f'Extracting to {destination}', flush=True)
        package.extractall(destination)
    for name in ['gcc.exe', 'g++.exe', 'gfortran.exe', 'mingw32-make.exe']:
        assert (destination / 'mingw64/bin' / name).is_file(), name
    (destination / 'installation.json').write_text(json.dumps(
        {'source': URL, 'sha256': actual, 'release': TAG}, indent=2))
    print('Installation complete.', flush=True)

if __name__ == '__main__':
    main()
