"""Contract tests for the 4.0.0 Python surface.

Shape-validation failures are results, not exceptions; only usage errors
(closed handle, non-1-D input, failed create) raise BetaParamError.
"""
from __future__ import annotations

import itertools
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import beta_parameterization as bp

SPHERE_4 = np.zeros(4)
ASYMMETRIC_4 = [0.10, 0.30, 0.15, 0.05]
INVALID_4 = [0.0, -1.8, 0.0, 0.0]
OPTION_COMBINATIONS = list(itertools.product([False, True], repeat=2))


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


# --- module-level (one-shot) tier ------------------------------------------

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


@pytest.mark.parametrize("conserve_volume,apply_com", OPTION_COMBINATIONS)
@pytest.mark.parametrize("max_params", [4, 8, 64])
def test_one_shot_matches_cache_bitwise(max_params: int, conserve_volume: bool,
                                        apply_com: bool) -> None:
    thetas = bp.theta_grid(64)
    one_shot = bp.radius_and_derivative(
        ASYMMETRIC_4, thetas, conserve_volume=conserve_volume, apply_com=apply_com)
    with bp.Cache(max_params, thetas) as cache:
        cached = cache.radius_and_derivative(
            ASYMMETRIC_4, conserve_volume=conserve_volume, apply_com=apply_com)
    assert one_shot.ok and cached.ok
    assert one_shot.radii.tobytes() == cached.radii.tobytes()
    assert one_shot.dr_dtheta.tobytes() == cached.dr_dtheta.tobytes()


def test_one_shot_empty_params_is_wrong_param_count() -> None:
    res = bp.radius_grid([], bp.theta_grid(16))
    assert res.status == bp.Status.wrong_param_count
    assert np.all(res.radii == 0.0)


def test_one_shot_too_many_params() -> None:
    res = bp.radius_grid(np.full(65, 0.001), bp.theta_grid(16))
    assert res.status == bp.Status.too_many_params
    assert np.all(res.radii == 0.0)


# --- cached tier -----------------------------------------------------------

def test_cache_attributes() -> None:
    with bp.Cache(8, bp.theta_grid(32)) as cache:
        assert cache.max_params == 8
        assert cache.n_thetas == 32
        assert not hasattr(cache, "n_params")
        assert not hasattr(cache, "conserve_volume")
        assert not hasattr(cache, "apply_com")


def test_cache_resolves_asymmetric_shape() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        resolved = cache.resolve_shape(ASYMMETRIC_4, conserve_volume=True,
                                       apply_com=True)
        raw = cache.resolve_shape(ASYMMETRIC_4)
    assert resolved.ok
    assert resolved.status == bp.Status.valid
    assert resolved.volume_factor > 0.0
    assert resolved.r_north > 0.0
    assert resolved.r_south > 0.0
    assert resolved.corrected_beta10 != ASYMMETRIC_4[0]
    # options default to off: beta10 untouched, factor exactly 1
    assert raw.ok
    assert raw.corrected_beta10 == ASYMMETRIC_4[0]
    assert raw.volume_factor == 1.0


def test_one_cache_serves_every_option_combination() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        results = {
            (cv, com): cache.radius_grid(ASYMMETRIC_4, conserve_volume=cv,
                                         apply_com=com)
            for cv, com in OPTION_COMBINATIONS
        }
    assert all(res.ok for res in results.values())
    blobs = {res.radii.tobytes() for res in results.values()}
    assert len(blobs) == 4


def test_short_vector_equals_zero_padded_bitwise() -> None:
    with bp.Cache(8, bp.theta_grid(64)) as cache:
        short = cache.radius_and_derivative([0.05, 0.2, 0.1], conserve_volume=True,
                                            apply_com=True)
        for padded_params in ([0.05, 0.2, 0.1, 0.0, 0.0],
                              [0.05, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
                              [0.05, 0.2, 0.1, -0.0, -0.0]):
            padded = cache.radius_and_derivative(padded_params, conserve_volume=True,
                                                 apply_com=True)
            assert short.ok and padded.ok
            assert short.radii.tobytes() == padded.radii.tobytes()
            assert short.dr_dtheta.tobytes() == padded.dr_dtheta.tobytes()


def test_invalid_shape_returns_status_not_exception() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        res = cache.radius_grid(INVALID_4)
    assert res.status == bp.Status.north_pole
    assert not res.ok
    assert np.all(res.radii == 0.0)
    assert res.message


def test_rejected_call_leaves_no_trace() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        before = cache.radius_grid(ASYMMETRIC_4, conserve_volume=True, apply_com=True)
        assert not cache.radius_grid(INVALID_4).ok
        assert not cache.radius_grid(np.zeros(5)).ok
        after = cache.radius_grid(ASYMMETRIC_4, conserve_volume=True, apply_com=True)
    assert before.ok and after.ok
    assert before.radii.tobytes() == after.radii.tobytes()


def test_unchecked_returns_negative_radii_for_invalid_shape() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        res = cache.radius_grid_unchecked(INVALID_4)
    assert res.ok
    assert np.any(res.radii < 0.0)


def test_unchecked_takes_no_conserve_volume() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        with pytest.raises(TypeError):
            cache.radius_grid_unchecked(ASYMMETRIC_4, conserve_volume=True)  # type: ignore[call-arg]


def test_com_failure_on_negative_volume_integral() -> None:
    # beta20 = -20 makes the volume integral negative; with COM that is a 103,
    # never a crash, and without COM the north pole fails first.
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        with_com = cache.radius_grid([0.0, -20.0], apply_com=True)
        without_com = cache.radius_grid([0.0, -20.0])
        unchecked = cache.radius_grid_unchecked([0.0, -20.0], apply_com=True)
    assert with_com.status == bp.Status.com_not_converged
    assert without_com.status == bp.Status.north_pole
    assert unchecked.status == bp.Status.com_not_converged
    assert np.all(unchecked.radii == 0.0)


def test_param_count_outside_range_is_a_status_not_an_exception() -> None:
    with bp.Cache(4, bp.theta_grid(64)) as cache:
        too_long = cache.radius_grid(np.zeros(5))
        empty = cache.radius_grid([])
    for res in (too_long, empty):
        assert res.status == bp.Status.wrong_param_count
        assert not res.ok
        assert np.all(res.radii == 0.0)


def test_max_params_limits() -> None:
    thetas = bp.theta_grid(16)
    with bp.Cache(bp.MAX_BETA_PARAMS_LIMIT, thetas) as cache:
        assert cache.max_params == 64
    with pytest.raises(bp.BetaParamError, match=r"too many.*status 1"):
        bp.Cache(bp.MAX_BETA_PARAMS_LIMIT + 1, thetas)
    with pytest.raises(bp.BetaParamError, match=r"status 5"):
        bp.Cache(0, thetas)


def test_cache_derivative_matches_finite_difference() -> None:
    thetas = bp.theta_grid(2001)
    with bp.Cache(4, thetas) as cache:
        res = cache.radius_and_derivative(ASYMMETRIC_4)
    assert res.ok
    fd = np.gradient(res.radii, thetas)
    assert np.allclose(res.dr_dtheta[5:-5], fd[5:-5], rtol=1.0e-4, atol=1.0e-6)


def test_inputs_need_not_be_contiguous_float64() -> None:
    """Lists of ints, float32 arrays and strided views are converted, not
    reinterpreted: they give the result of the equivalent float64 array."""
    thetas = bp.theta_grid(64)
    reference_params = np.array([0.0, 0.25, 0.125, 0.0625])  # exact in float32
    with bp.Cache(4, thetas) as cache:
        reference = cache.radius_grid(reference_params)
        as_float32 = cache.radius_grid(reference_params.astype(np.float32))
        strided = np.zeros(8)
        strided[::2] = reference_params
        as_view = cache.radius_grid(strided[::2])
        as_list = cache.radius_grid(list(reference_params))
        int_sphere = cache.radius_grid([0, 0, 0, 0])
    with bp.Cache(4, list(thetas)) as cache_from_list:
        from_list_thetas = cache_from_list.radius_grid(reference_params)
    assert reference.ok
    for other in (as_float32, as_view, as_list, from_list_thetas):
        assert other.ok
        assert other.radii.tobytes() == reference.radii.tobytes()
    assert int_sphere.ok
    assert np.allclose(int_sphere.radii, 1.0, rtol=0.0, atol=1.0e-14)


# --- concurrency -----------------------------------------------------------

def test_threads_share_one_cache() -> None:
    """Contract family 3 from Python: ctypes releases the GIL during a call,
    so the workers really do run inside the library at the same time."""
    thetas = bp.theta_grid(256)
    rng = np.random.default_rng(20261001)
    shapes = [rng.uniform(-0.05, 0.05, size=int(n))
              for n in rng.integers(1, 9, size=200)]
    jobs = [(shape, cv, com)
            for shape, (cv, com) in zip(shapes, itertools.cycle(OPTION_COMBINATIONS))]

    with bp.Cache(8, thetas) as cache:
        def run(job: tuple[np.ndarray, bool, bool]) -> tuple[bp.Status, bytes, bytes]:
            shape, cv, com = job
            res = cache.radius_and_derivative(shape, conserve_volume=cv, apply_com=com)
            return res.status, res.radii.tobytes(), res.dr_dtheta.tobytes()

        serial = [run(job) for job in jobs]
        with ThreadPoolExecutor(max_workers=8) as pool:
            # every job submitted eight times, interleaved across the workers
            threaded = list(pool.map(run, jobs * 8))

    assert all(status == bp.Status.valid for status, _, _ in serial)
    assert threaded == serial * 8


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
    with pytest.raises(bp.BetaParamError, match=r"status 105"):
        bp.Cache(4, np.linspace(0.0, np.pi, 32))


def test_single_theta_rejected_at_create() -> None:
    with pytest.raises(bp.BetaParamError, match=r"status 3"):
        bp.Cache(4, [1.0])


# --- diagnostics -----------------------------------------------------------

def test_status_message_covers_every_code() -> None:
    for status in bp.Status:
        message = bp.status_message(status)
        assert message and "unknown" not in message
    assert "north" in bp.status_message(bp.Status.north_pole)
    assert bp.Status.node_set_mismatch == 106
    # code 6 is retired: no enum member, and the library no longer names it
    assert 6 not in {int(status) for status in bp.Status}
    assert "unknown" in bp.status_message(6)


# --- removed surface -------------------------------------------------------

@pytest.mark.parametrize("name", [
    "NodeSet", "MESSAGE_BUFFER_SIZE", "radius_grid_with_com_shift",
    "build_node_set", "radius_grid_standalone",
    "radius_grid_standalone_with_com_shift",
    "CACHE_MAX_PARAMS", "Tables",
])
def test_removed_names_absent(name: str) -> None:
    assert not hasattr(bp, name)
    assert name not in bp.__all__
    assert not hasattr(bp.Cache, name)


def test_removed_status_member_absent() -> None:
    assert not hasattr(bp.Status, "tables_not_initialized")


def test_cache_constructor_takes_no_options() -> None:
    with pytest.raises(TypeError):
        bp.Cache(4, bp.theta_grid(16), conserve_volume=True)  # type: ignore[call-arg]
