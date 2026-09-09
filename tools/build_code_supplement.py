"""Package the paper source and committed data."""
import argparse
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRS = ('checks', 'configs', 'figure_data', 'figure_scripts',
               'pinned_ensembles', 'provenance', 'studies', 'tools')
ROOT_FILES = {'README.md', 'LICENSE', 'requirements.txt', 'CITATION.cff',
              'expected_values.json', 'config_template.jsonc', '.gitignore'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'output/code-supplement.zip')
    output = parser.parse_args().output.resolve()
    files = [p for p in ROOT.iterdir()
             if p.is_file() and (p.suffix == '.py' or p.name in ROOT_FILES)]
    for directory in SOURCE_DIRS:
        files.extend(p for p in (ROOT / directory).rglob('*')
                     if p.is_file() and '__pycache__' not in p.parts and p.name != '.DS_Store')
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(files):
            archive.write(path, 'source/' + path.relative_to(ROOT).as_posix())
    print(f'Wrote {output} with {len(files)} source and data files')


if __name__ == '__main__':
    main()
