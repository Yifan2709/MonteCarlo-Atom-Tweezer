import numpy as np

from level0_static_trap.potentials import gaussian_force, gaussian_potential


DEPTH = 2.5e-27
WAIST = 1.17e-6
CENTER = 0.31e-6


def test_center_potential_is_negative_depth():
    assert gaussian_potential(CENTER, DEPTH, WAIST, CENTER) == -DEPTH


def test_center_force_is_zero():
    assert gaussian_force(CENTER, DEPTH, WAIST, CENTER) == 0.0


def test_force_on_both_sides_points_to_center():
    assert gaussian_force(CENTER - 0.2 * WAIST, DEPTH, WAIST, CENTER) > 0.0
    assert gaussian_force(CENTER + 0.2 * WAIST, DEPTH, WAIST, CENTER) < 0.0


def test_force_matches_negative_finite_difference_gradient():
    x = np.linspace(CENTER - WAIST, CENTER + WAIST, 2001)
    potential = gaussian_potential(x, DEPTH, WAIST, CENTER)
    numerical_force = -np.gradient(potential, x, edge_order=2)
    exact_force = gaussian_force(x, DEPTH, WAIST, CENTER)
    scale = np.max(np.abs(exact_force))
    assert np.max(np.abs(numerical_force[2:-2] - exact_force[2:-2])) / scale < 3e-6


def test_scalar_and_array_inputs_are_supported():
    scalar = gaussian_potential(CENTER, DEPTH, WAIST, CENTER)
    array = gaussian_potential(np.array([CENTER]), DEPTH, WAIST, CENTER)
    assert np.asarray(scalar).shape == ()
    assert array.shape == (1,)

