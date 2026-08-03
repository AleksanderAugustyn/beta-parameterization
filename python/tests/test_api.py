"""Contract tests for the 3.0.0 Python surface.

Shape-validation failures are results, not exceptions; only usage errors
(closed handle, non-1-D input, failed create) raise BetaParamError.
"""
from __future__ import annotations

import numpy as np
import pytest

import beta_parameterization as bp

SPHERE_4 = np.zeros(4)
ASYMMETRIC_4 = [0.10, 0.30, 0.15, 0.05]
INVALID_4 = [0.0, -1.8, 0.0, 0.0]


# --- theta_grid ------------------------------------------------------------

def test_theta_grid_is_open_and_uniform() -> None:
    thetas = bp.theta_grid(721)
    assert thetas.size == 721
    assert thetas[0] > 0.0
    assert thetas[-1] < np.pi
    spacing = np.diff(thetas)
    assert np.allclose(spacing, np.pi / 722.0, rtol=0.0, atol=1.0e-12)


def test_theta_grid_excludes_both_poles_for_small_n() -> None:
    thetas = bp.theta_grid(3)
    assert np.allclose(thetas, np.array([1.0, 2.0, 3.0]) * np.pi / 4.0)


# --- module-level (standalone) tier ---------------------------------------

def test_sphere_radius_grid_is_unit() -> None:
    res = bp.radius_grid(SPHERE_4, bp.theta_grid(181))
    assert res.ok
    assert res.status == bp.Status.valid
    assert np.allclose(res.radii, 1.0, rtol=0.0, atol=1.0e-14)


def test_sphere_derivative_is_zero() -> None:
    res = bp.radius_and_derivative(SPHERE_4, bp.theta_grid(64))
    assert res.ok
    assert np.allclose(res.radii, 1.0, rtol=0.0, atol=1.0e-14)
    assert np.allclose(res.dr_dtheta, 0.0, rtol=0.0, atol=1.0e-14)


def test_standalone_matches_cache() -> None:
    thetas = bp.theta_grid(64)
    standalone = bp.radius_grid(ASYMMETRIC_4, thetas,
                                conserve_volume=True, apply_com=True)
    with bp.Cache(4, thetas, conserve_volume=True, apply_com=True) as cache:
        cached = cache.radius_grid(ASYMMETRIC_4)
    assert standalone.ok and cached.ok
    assert np.array_equal(standalone.radii, cached.radii)


# --- cached tier -----------------------------------------------------------

def test_cache_resolves_asymmetric_shape() -> None:
    with bp.Cache(4, bp.theta_grid(64),
                  conserve_volume=True, apply_com=True) as cache:
        resolved = cache.resolve_shape(ASYMMETRIC_4)
    assert resolved.ok
    assert resolved.status == bp.Status.valid
    assert resolved.volume_factor > 0.0
    assert resolved.r_north > 0.0
    assert resolved.r_south > 0.0
    assert resolved.corrected_beta10 != ASYMMETRIC_4[0]


def test_invalid_shape_returns_status_not_exception() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        res = cache.radius_grid(INVALID_4)
    assert res.status != bp.Status.valid
    assert not res.ok
    assert np.all(res.radii == 0.0)
    assert res.message


def test_unchecked_returns_negative_radii_for_invalid_shape() -> None:
    with bp.Cache(4, bp.theta_grid(64), apply_com=False) as cache:
        res = cache.radius_grid_unchecked(INVALID_4)
    assert res.ok
    assert np.any(res.radii < 0.0)


def test_wrong_param_count_is_a_status_not_an_exception() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        res = cache.radius_grid([0.0, 0.2, 0.0])
    assert res.status == bp.Status.wrong_param_count
    assert not res.ok
    assert np.all(res.radii == 0.0)


def test_too_many_params_raises() -> None:
    with pytest.raises(bp.BetaParamError, match="too many"):
        bp.Cache(9, bp.theta_grid(64))


def test_cache_derivative_matches_finite_difference() -> None:
    thetas = bp.theta_grid(2001)
    with bp.Cache(4, thetas) as cache:
        res = cache.radius_and_derivative(ASYMMETRIC_4)
    assert res.ok
    fd = np.gradient(res.radii, thetas)
    assert np.allclose(res.dr_dtheta[5:-5], fd[5:-5], rtol=1.0e-4, atol=1.0e-6)


# --- usage errors ----------------------------------------------------------

def test_use_after_close_raises() -> None:
    cache = bp.Cache(4, bp.theta_grid(16))
    cache.close()
    cache.close()  # idempotent
    with pytest.raises(bp.BetaParamError, match="closed"):
        cache.radius_grid(SPHERE_4)


def test_non_1d_params_raise() -> None:
    with bp.Cache(4, bp.theta_grid(16)) as cache:
        with pytest.raises(bp.BetaParamError):
            cache.radius_grid(np.zeros((2, 2)))


def test_pole_theta_rejected_at_create() -> None:
    with pytest.raises(bp.BetaParamError):
        bp.Cache(4, np.linspace(0.0, np.pi, 32))


# --- diagnostics -----------------------------------------------------------

def test_status_message_covers_every_code() -> None:
    for status in bp.Status:
        message = bp.status_message(status)
        assert message and "unknown" not in message
    assert "north" in bp.status_message(bp.Status.north_pole)
    assert bp.Status.node_set_mismatch == 106


# --- removed 2.x surface ---------------------------------------------------

@pytest.mark.parametrize("name", [
    "NodeSet", "MESSAGE_BUFFER_SIZE", "radius_grid_with_com_shift",
    "build_node_set", "radius_grid_standalone",
    "radius_grid_standalone_with_com_shift",
])
def test_removed_names_absent(name: str) -> None:
    assert not hasattr(bp, name)
    assert name not in bp.__all__
    assert not hasattr(bp.Cache, name)
