"""Reproduce and verify the two-step heat study against committed data.

Usage from the repository root:
  python checks/verify_two_step.py                 # predict, run, validate, compare
  python checks/verify_two_step.py --no-rerun      # compare existing reproduction
  python checks/verify_two_step.py --pinned-only   # check archived design and moments

Fresh files go to output/heat_two_step_reproduction/. The committed data and
historical design lock are read-only inputs. Rerunning a known design reproduces
the original experiment; it does not create new evidence of prospective design.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from statistics import median
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from studies import study_t9_heat_two_step as study

PIN = ROOT / 'pinned_ensembles/heat_two_step'
OUT = ROOT / 'output/heat_two_step_reproduction'
FILES = ('predictions.json', 'validation_spec.json', 'ensembles.json', 'validation.json')
REL_TOL, ABS_TOL = 1e-9, 1e-12
# These are the only machine-dependent fields in the four study files.
TIMING_FIELDS = {'runtime', 'mean_run_s'}


def compare(expected, actual, path='root'):
    """Strict shapes/types; numeric tolerance only for finite floating values."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            raise AssertionError(f'{path}: expected a dictionary')
        keys = set(expected) - TIMING_FIELDS
        if keys != set(actual) - TIMING_FIELDS:
            raise AssertionError(f'{path}: keys differ')
        for key in sorted(keys):
            compare(expected[key], actual[key], f'{path}/{key}')
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f'{path}: list length/type differs')
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f'{path}[{i}]')
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise AssertionError(f'{path}: expected a number')
        if not (math.isfinite(expected) and math.isfinite(actual)
                and math.isclose(expected, actual, rel_tol=REL_TOL, abs_tol=ABS_TOL)):
            raise AssertionError(f'{path}: {expected!r} != {actual!r}')
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f'{path}: {expected!r} != {actual!r}')


def read(directory, name):
    return json.loads((directory / name).read_text())


def check_rows(rows, counts, require_realizations=True):
    expected = {(n, m, ref, conv) for n in counts for m in study.M_LIST
                for ref in ('finite', 'infinite') for conv in ('edge', 'center')}
    keys = [(r['N'], r['M'], r['reference'], r['convention']) for r in rows]
    if len(keys) != len(expected) or set(keys) != expected:
        raise AssertionError('Missing, duplicate, or unexpected comparison configurations')
    for row in rows:
        values = row.get('per_realization_sq')
        if row['S'] != study.S:
            raise AssertionError('Realization count differs from the study design')
        if values is None:
            if require_realizations:
                raise AssertionError('Missing per-realization squared errors')
        else:
            if len(values) != study.S:
                raise AssertionError('Realization vector length differs from the study design')
            if not all(math.isfinite(v) and v >= 0 for v in values):
                raise AssertionError('Invalid squared realization error')
            compare(sum(values) / len(values), row['E_total']**2, 'RMS from realizations')
        compare(row['E_bias']**2 + row['E_spread']**2, row['E_total']**2,
                'bias/spread/total identity')


def check_archived_design():
    lock = ROOT / 'provenance/heat_two_step/DESIGN_LOCKED_AT.txt'
    lines = lock.read_text().splitlines()
    for name in ('predictions.json', 'validation_spec.json'):
        matches = [line.split()[0] for line in lines[1:] if line.endswith('/' + name)]
        if len(matches) != 1 or hashlib.sha256((PIN / name).read_bytes()).hexdigest() != matches[0]:
            raise AssertionError(f'Historical design hash mismatch: {name}')
    pred = read(PIN, 'predictions.json')
    compare(pred['per_M'], [study.predictions(m) for m in study.M_LIST], 'predictions')
    compare(read(PIN, 'validation_spec.json'), study.validation_spec(), 'count rule')
    for m in study.M_LIST:
        # The covariance and scalar variance are implemented separately.
        covariance = study.covariance(m, 400)
        compare(float(covariance.trace()) * 400, study.predictions(m)['V_h'], 'covariance trace')
    print('PASS historical prediction/spec hashes, recomputed predictions, count rule, covariance trace')


def verify(rerun=True, pinned_only=False):
    try:
        check_archived_design()
        data = PIN if pinned_only else OUT
        if rerun and not pinned_only:
            for mode in ('predict', 'run', 'validate'):
                subprocess.run([sys.executable, str(ROOT / 'studies/study_t9_heat_two_step.py'),
                                mode, '--output-dir', str(OUT)], cwd=ROOT, check=True)
        if not pinned_only:
            for name in FILES:
                expected, actual = read(PIN, name), read(data, name)
                if name == 'validation.json':
                    # The historical validation stores summary moments only.
                    # Fresh runs also retain realization errors, checked below.
                    for expected_row, actual_row in zip(expected['rows'], actual['rows']):
                        if 'per_realization_sq' not in expected_row:
                            actual_row.pop('per_realization_sq', None)
                compare(expected, actual, name)
                print(f'PASS {name} matches committed data (timings excluded)')
        ensemble = read(data, 'ensembles.json')
        validation = read(data, 'validation.json')
        compare(read(data, 'validation_spec.json'), validation['spec'], 'validation design')
        check_rows(ensemble['rows'], study.N_LIST)
        check_rows(validation['rows'], [validation['spec']['N']], require_realizations=not pinned_only)
        pred_by_m = {p['M']: p for p in read(data, 'predictions.json')['per_M']}
        ratios = []
        for row in ensemble['rows']:
            p = pred_by_m[row['M']]
            b2 = p[f"{row['reference']}_{row['convention']}"]['B2']
            ratios.append(row['E_total'] / math.sqrt(b2 + p['V_h']/row['N']))
        selected = next(row for row in validation['rows'] if row['M'] == validation['spec']['M']
                        and row['reference'] == 'infinite' and row['convention'] == 'edge')
        predicted = validation['spec']['edge_convention']['predicted_E_total_at_N']
        ratios.sort()
        print(f'PASS {len(ensemble["rows"])} production configurations and '
              f'{len(validation["rows"])} validation configurations with 30 realizations each')
        print(f'Measured/predicted total range {ratios[0]:.6f} to {ratios[-1]:.6f}, '
              f'median {median(ratios):.6f}')
        print(f'Validation N={selected["N"]}: measured {selected["E_total"]:.9f}, '
              f'predicted {predicted:.9f}, target {validation["spec"]["target_rms_error"]:.9f}')
        print('The count rule concerns expected squared error; the finite ensemble is not a target guarantee.')
        print('TWO-STEP VERIFY: PASS' + (' (archived data checks, no solver rerun)' if pinned_only else ''))
        return 0
    except (AssertionError, OSError, ValueError, KeyError, StopIteration, subprocess.CalledProcessError) as exc:
        print(f'TWO-STEP VERIFY: FAIL: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--no-rerun', action='store_true')
    group.add_argument('--pinned-only', action='store_true')
    args = parser.parse_args()
    raise SystemExit(verify(rerun=not args.no_rerun, pinned_only=args.pinned_only))
