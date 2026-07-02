"""Python == Fortran: bindings must reproduce the library bit-for-bit."""
from __future__ import annotations

import numpy as np
import pytest

import beta_parameterization as bp

# Golden shapes shared with tests/beta_param_golden_test.f08 (Task 8).
# The _EXPECTED constants are pasted from golden_capture output after Task 8.
G1_PARAMS = [0.0, 0.215, 0.0, 0.095]
G2_PARAMS = [0.0, 0.85, 0.35, 0.18, 0.05, 0.02]
GOLDEN_IDX = [0, 30, 60, 90, 120, 150, 180]


def test_sphere_exact() -> None:
    with bp.Cache(max_beta_params=8, n_grid=181) as cache:
        res = cache.radius_grid(np.zeros(4))
    assert res.status == bp.Status.VALID
    assert res.ok
    np.testing.assert_array_equal(res.radii, np.ones(181))


def test_cache_equals_standalone() -> None:
    params = [0.0, 0.25, 0.10, 0.05]
    with bp.Cache(max_beta_params=4, n_grid=181) as cache:
        res_cache = cache.radius_grid(params)
    res_sa = bp.radius_grid_standalone(params, n_grid=181)
    np.testing.assert_array_equal(res_cache.radii, res_sa.radii)


def test_error_codes_surface_as_status() -> None:
    with bp.Cache(max_beta_params=4, n_grid=181) as cache:
        res = cache.radius_grid([0.0, 4.0])          # interior negative
        assert res.status == bp.Status.ERROR_INTERIOR_NEGATIVE
        assert not res.ok
        assert res.message                            # non-empty, content not asserted
        res = cache.radius_grid([0.1] * 5)            # too many params
        assert res.status == bp.Status.ERROR_TOO_MANY_PARAMS


def test_com_shift_returns_corrected_beta10() -> None:
    with bp.Cache(max_beta_params=4, n_grid=181) as cache:
        res = cache.radius_grid_with_com_shift([0.30, 0.60, 0.40, 0.10])
    assert res.status == bp.Status.VALID
    assert res.corrected_beta10 is not None
    assert abs(res.corrected_beta10 - 0.30) > 1e-6


def test_constructor_raises_on_bad_init() -> None:
    with pytest.raises(bp.BetaParamError):
        bp.Cache(max_beta_params=0, n_grid=181)
    with pytest.raises(bp.BetaParamError):
        bp.Cache(max_beta_params=bp.MAX_BETA_PARAMS_LIMIT + 1, n_grid=181)


def test_theta_grid() -> None:
    t = bp.theta_grid(181)
    assert t.shape == (181,)
    assert t[0] == 0.0
    assert t[-1] == pytest.approx(np.pi, rel=1e-15)


def test_use_after_close_raises() -> None:
    cache = bp.Cache(max_beta_params=4, n_grid=41)
    cache.close()
    with pytest.raises(bp.BetaParamError):
        cache.radius_grid([0.0, 0.2])
