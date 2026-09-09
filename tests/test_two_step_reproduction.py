"""Regression checks for omissions or corrupted data in the reproduction gate."""
import copy
import subprocess
import sys
import tempfile
from pathlib import Path
import unittest

from checks.verify_two_step import PIN, check_rows, compare, read
from studies.study_t9_heat_two_step import N_LIST, covariance, run_one


class TwoStepReproductionTests(unittest.TestCase):
    def test_rejects_modified_realization_error(self):
        rows = copy.deepcopy(read(PIN, 'ensembles.json')['rows'])
        rows[0]['per_realization_sq'][0] += 0.001
        with self.assertRaises(AssertionError):
            check_rows(rows, N_LIST)

    def test_rejects_missing_configuration(self):
        rows = read(PIN, 'ensembles.json')['rows'][:-1]
        with self.assertRaises(AssertionError):
            check_rows(rows, N_LIST)

    def test_rejects_nonfinite_numeric_match(self):
        with self.assertRaises(AssertionError):
            compare(1.0, float('nan'))

    def test_rejects_integer_as_boolean(self):
        with self.assertRaises(AssertionError):
            compare(1, True)

    def test_historical_validation_requires_explicit_summary_only_mode(self):
        data = read(PIN, 'validation.json')
        check_rows(data['rows'], [data['spec']['N']], require_realizations=False)
        with self.assertRaises(AssertionError):
            check_rows(data['rows'], [data['spec']['N']])


class StudyInputTests(unittest.TestCase):
    def test_equal_allocation_rejects_invalid_counts_before_running(self):
        for count in (0, 1, 3, -2, 4.5, True):
            with self.subTest(count=count):
                with self.assertRaises(ValueError):
                    run_one(count, 5000)
                with self.assertRaises(ValueError):
                    covariance(50, count)

    def test_analysis_import_does_not_parse_cli_or_write_files(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            code = (
                f"import sys; sys.path.insert(0, {str(root)!r}); "
                "sys.argv = ['import-check', '--not-a-study-option']; "
                "import studies.analyze_t9_heat_two_step"
            )
            result = subprocess.run(
                [sys.executable, '-c', code],
                cwd=directory,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, '')
            self.assertEqual(list(Path(directory).iterdir()), [])


if __name__ == '__main__':
    unittest.main()
