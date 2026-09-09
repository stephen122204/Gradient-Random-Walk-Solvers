"""Regression checks for omissions or corrupted data in the reproduction gate."""
import copy
import unittest

from checks.verify_two_step import PIN, check_rows, compare, read
from studies.study_t9_heat_two_step import N_LIST


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


if __name__ == '__main__':
    unittest.main()
