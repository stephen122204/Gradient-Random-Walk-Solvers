"""Fast release tests for solver and experiment invariants used in Paper 1."""

import unittest

import numpy as np

from config import validate_config_dict
from simulation import _reflect_arrays, _validated_step_count
from studies.study_t7_heat_grid_paired import _decompose, _operator_control
from studies.study_t8_burgers_controls import (
    _fd_reference,
    _phi0_from_grid,
    validate_pipeline,
)
from verify_solver import exact_heat_step


class TimeAndBoundaryTests(unittest.TestCase):
    def test_fixed_step_count_rejects_silent_truncation(self):
        self.assertEqual(_validated_step_count(0.5, 0.005), 100)
        with self.assertRaises(ValueError):
            _validated_step_count(1.0, 0.3)

    def test_arbitrarily_large_reflection_and_neumann_sign(self):
        positions = np.array([25.0, -14.0, 1.0])
        weights = np.ones(3)
        x, w = _reflect_arrays(
            positions, weights, 4.0, "dirichlet", "neumann"
        )
        self.assertTrue(np.all((0.0 <= x) & (x <= 4.0)))
        np.testing.assert_allclose(x, [1.0, 2.0, 1.0])
        np.testing.assert_allclose(w, [-1.0, 1.0, 1.0])


class FHNFormulaTests(unittest.TestCase):
    def test_reaction_statistic_matches_derivative_and_has_zero_integral(self):
        diffusion = 0.5
        theta = np.sqrt(2.0) * 0.25
        c2 = -1.5 * diffusion
        c1 = 1.5 * diffusion - theta
        c0 = 0.5 * theta - 0.25 * diffusion

        def reaction(u):
            return u * (1.0 - u) * (
                theta / 2.0 - diffusion * (1.0 - 2.0 * u) / 4.0
            )

        for u in np.linspace(0.1, 0.9, 9):
            eps = 1e-6
            numerical = (reaction(u + eps) - reaction(u - eps)) / (2.0 * eps)
            analytic = c2 * u**2 + c1 * u + c0
            self.assertAlmostEqual(numerical, analytic, places=9)
        self.assertAlmostEqual(c2 / 3.0 + c1 / 2.0 + c0, 0.0, places=14)


class HeatAttributionTests(unittest.TestCase):
    def test_bias_spread_total_identity(self):
        runs = np.array([[0.0, 0.8, 1.0], [0.0, 1.0, 1.0],
                         [0.0, 1.2, 1.0]])
        reference = np.array([0.0, 1.0, 1.0])
        bias, spread, total, residual = _decompose(runs, reference, 0.5)
        self.assertAlmostEqual(total**2, bias**2 + spread**2, places=14)
        self.assertLess(residual, 1e-14)

    def test_operator_control_localizes_center_alignment_error(self):
        length, center, alpha, final_time = 4.0, 2.0, 0.5, 0.5
        fixed = {}
        for bins in (300, 400):
            edges = np.linspace(0.0, length, bins + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])
            fixed[bins] = {
                'edges': edges,
                'centers': centers,
                'dx': float(edges[1] - edges[0]),
                'u_exact_center': exact_heat_step(
                    centers, final_time, center, 0.0, 1.0, alpha
                ),
                'u_exact_edge': exact_heat_step(
                    edges[1:], final_time, center, 0.0, 1.0, alpha
                ),
            }
        control = _operator_control(length, center, alpha, final_time, fixed)
        for bins in (300, 400):
            item = control[f'M{bins}']
            self.assertGreater(item['center_compare'], item['edge_compare'])
            self.assertLess(item['edge_compare'], 0.0012)
        self.assertLess(control['finite_domain_reference_gap'], 0.0012)


class BurgersControlTests(unittest.TestCase):
    def test_parameterized_pipeline_matches_packaged_solver(self):
        self.assertLess(validate_pipeline(), 1e-12)

    def test_exact_transformed_boundary_data_removes_boundary_model_error(self):
        pinned = _fd_reference(400, 4.0, boundary='pinned')
        exact = _fd_reference(400, 4.0, boundary='exact')
        self.assertGreater(pinned['rmse_u'], 0.1)
        self.assertLess(exact['rmse_u'], 5e-4)
        self.assertGreater(pinned['rmse_u'] / exact['rmse_u'], 400.0)

    def test_transformed_gradient_weights_telescope(self):
        x = np.linspace(0.0, 4.0, 400)
        u = -np.tanh(x - 2.0)
        phi = _phi0_from_grid(x, u, 0.5)
        self.assertAlmostEqual(float(np.diff(phi).sum()),
                               float(phi[-1] - phi[0]), places=14)

    def test_config_rejects_nonpositive_cole_hopf_viscosity(self):
        data = {
            'equation_type': 'burgers', 'domain_type': 'Finite',
            'domain_size': 4.0,
            'boundary_conditions': {
                'LEFT': {'type': 'Dirichlet', 'value': 0.0},
                'RIGHT': {'type': 'Dirichlet', 'value': 0.0},
            },
            'diff_constant': 0.0, 'time_step': 0.005,
            'total_time': 0.5, 'num_points': 400,
            'burgers_initial_condition': {'type': 'stationary_shock'},
        }
        with self.assertRaises(ValueError):
            validate_config_dict(data)


if __name__ == '__main__':
    unittest.main()
