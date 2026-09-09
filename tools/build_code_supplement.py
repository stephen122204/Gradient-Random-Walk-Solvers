"""Build a self-contained source/data supplement without generated solver output.

Usage: python tools/build_code_supplement.py --output output/code-supplement.zip
Optional --check-outputs DIR includes existing verification logs as supporting
records. The archive manifest hashes every included source and log file.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRS = ('checks', 'configs', 'docs', 'figure_data', 'figure_scripts',
               'pinned_ensembles', 'provenance', 'studies', 'tests', 'tools')


def build(output, check_outputs=None):
    files = {}
    for path in ROOT.iterdir():
        if path.is_file() and (path.suffix == '.py' or path.name in {
            'README.md', 'LICENSE', 'requirements.txt', 'CITATION.cff',
            'expected_values.json', 'config_template.jsonc', '.gitignore'
        }):
            files['source/' + path.name] = path
    for directory in SOURCE_DIRS:
        for path in (ROOT / directory).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts and path.name != '.DS_Store':
                files['source/' + path.relative_to(ROOT).as_posix()] = path
    if check_outputs:
        for path in check_outputs.rglob('*'):
            if path.is_file() and path.suffix in {'.txt', '.json', '.log'}:
                files['check-outputs/' + path.relative_to(check_outputs).as_posix()] = path
    try:
        if not (ROOT / '.git').exists():
            raise FileNotFoundError('Extracted source has no repository metadata')
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                         stderr=subprocess.DEVNULL, text=True).strip()
        status = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT,
                                         stderr=subprocess.DEVNULL, text=True).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit, status = None, 'Git metadata unavailable in extracted source'
    manifest = {
        'built_at_utc': datetime.now(timezone.utc).isoformat(),
        'git_head': commit, 'working_tree_status': status,
        'source_of_truth': 'Included file hashes identify this package, including uncommitted changes.',
        'files': {name: {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                         'bytes': path.stat().st_size} for name, path in sorted(files.items())},
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name, path in sorted(files.items()):
            archive.write(path, name)
        archive.writestr('PROVENANCE.json', json.dumps(manifest, indent=2) + '\n')
        archive.writestr('README.txt',
            'Revised GRW paper code and data supplement\n\n'
            'Extract this ZIP and work inside source/. Create a Python 3.11 environment\n'
            'and install requirements.txt. Follow source/README.md and\n'
            'source/docs/METHODS_AND_REPRODUCIBILITY.md for all reproduction commands.\n'
            'The package includes the two-step extension and historical design record.\n'
            'PROVENANCE.json hashes every included source and verification record.\n'
            'Check logs document the run that built this package, when supplied.\n'
            'No Zenodo version is created or published by this local packaging command.\n')
    print(f'Wrote {output} with {len(files)} source/log files and a SHA-256 manifest')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'output/code-supplement.zip')
    parser.add_argument('--check-outputs', type=Path)
    args = parser.parse_args()
    build(args.output.resolve(), args.check_outputs.resolve() if args.check_outputs else None)
