"""High-level API: Status, result objects, Cache, one-shot functions.

Two tiers, matching the C API. The module-level functions build, use and
discard their own cache per call; :class:`Cache` holds the tables that depend
only on ``max_params`` and the theta set, built once and reused for any number
of shapes.

A :class:`Cache` is immutable after construction and may be shared between
threads. ``conserve_volume`` and ``apply_com`` are per-call options.

Shape-validation failures come back as result objects carrying a
:class:`Status`; :class:`BetaParamError` is raised only for usage errors — a
closed handle, non-1-D input, or a failed create.
"""
from __future__ import annotations

import ctypes
from dataclasses import dataclass
from enum import IntEnum
from typing import Optional

import numpy as np
import numpy.typing as npt

from ._cdefs import MAX_BETA_PARAMS_LIMIT, c_dbl_p, configure
from ._libloader import load_library

_lib: Optional[ctypes.CDLL] = None


def _get_lib() -> ctypes.CDLL:
    """Load and configure the shared library once per process."""
    global _lib
    if _lib is None:
        _lib = configure(load_library())
    return _lib


class BetaParamError(RuntimeError):
    """Raised for usage errors: failed create, closed handle, non-1-D input."""


class Status(IntEnum):
    """Status codes of the C API (``BETA_PARAM_*`` in the header).

    0-5 are the shared shape-parameterization contract codes (6 is retired);
    100+ are this library's own codes, append-only.
    """

    valid = 0
    too_many_params = 1
    cache_not_initialized = 2
    invalid_grid = 3
    wrong_param_count = 4
    invalid_init = 5
    north_pole = 100
    south_pole = 101
    interior_negative = 102
    com_not_converged = 103
    invalid_buffer_size = 104
    pole_node = 105
    node_set_mismatch = 106


def status_message(status: int) -> str:
    """Fixed description of a status code.

    Parameters
    ----------
    status : int
        A :class:`Status` value or its integer code.

    Returns
    -------
    str
        The library's static message; ``"unknown status code"`` for codes it
        does not know.
    """
    return _get_lib().beta_param_status_message(int(status)).decode()


@dataclass(frozen=True)
class RadiusGridResult:
    """One radius evaluation. Radii are zero-filled when ``status`` is nonzero."""

    radii: npt.NDArray[np.float64]
    status: Status

    @property
    def ok(self) -> bool:
        """True when the evaluation succeeded."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class RadiusDerivativeResult:
    """R(theta) and dR/dtheta. Both buffers are zero-filled on failure."""

    radii: npt.NDArray[np.float64]
    dr_dtheta: npt.NDArray[np.float64]
    status: Status

    @property
    def ok(self) -> bool:
        """True when the evaluation succeeded."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class ResolvedShape:
    """A shape resolved without evaluating a grid.

    ``corrected_beta10`` is a beta-space value and is never volume-scaled;
    ``r_north`` and ``r_south`` are. ``volume_factor`` is exactly 1.0 when the
    shape was resolved with ``conserve_volume=False``.
    """

    corrected_beta10: float
    r_north: float
    r_south: float
    volume_factor: float
    status: Status

    @property
    def ok(self) -> bool:
        """True when the shape resolved successfully."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


def theta_grid(n: int) -> npt.NDArray[np.float64]:
    """The open uniform theta grid ``theta_i = i*pi/(n+1)``, ``i = 1..n``.

    Open: neither pole is included, because the pole guard rejects a pole
    theta as cache input. Pole radii come analytically from
    :meth:`Cache.resolve_shape`.

    Parameters
    ----------
    n : int
        Number of interior points.

    Returns
    -------
    numpy.ndarray
        ``n`` thetas in radians, strictly inside ``(0, pi)``.
    """
    return np.arange(1, int(n) + 1, dtype=np.float64) * np.pi / (int(n) + 1.0)


def _as_1d(values: npt.ArrayLike, name: str) -> npt.NDArray[np.float64]:
    """Contiguous 1-D float64 view of ``values``; BetaParamError otherwise."""
    arr = np.ascontiguousarray(values, dtype=np.float64)
    if arr.ndim != 1:
        raise BetaParamError(f"{name} must be 1-D, got shape {arr.shape}")
    return arr


def _ptr(arr: npt.NDArray[np.float64]) -> ctypes.POINTER(ctypes.c_double):  # type: ignore[valid-type]
    """Raw double pointer to a contiguous float64 array."""
    return arr.ctypes.data_as(c_dbl_p)


def radius_grid(params: npt.ArrayLike, thetas: npt.ArrayLike,
                conserve_volume: bool = False,
                apply_com: bool = False) -> RadiusGridResult:
    """R(theta) for one shape, a cache built and discarded per call.

    Parameters
    ----------
    params : array_like
        Beta parameters, 1 .. 64 entries.
    thetas : array_like
        Evaluation angles in radians; at least 2, none at a pole.
    conserve_volume : bool, optional
        Rescale radii to fixed volume.
    apply_com : bool, optional
        Apply the centre-of-mass correction.

    Returns
    -------
    RadiusGridResult
        Radii and status; radii are zero-filled on failure.

    Notes
    -----
    ``params`` and ``thetas`` must be finite. Non-finite input is undefined
    behavior: the library cannot detect NaN under fast-math, so the call
    returns ``Status.valid`` with NaN radii instead of an error. Screen inputs
    before calling.

    The result is bitwise identical to :meth:`Cache.radius_grid` on a cache
    with ``max_params >= len(params)`` over the same thetas.
    """
    p = _as_1d(params, "params")
    t = _as_1d(thetas, "thetas")
    radii = np.zeros(t.size, dtype=np.float64)
    status = _get_lib().beta_param_radius_grid_standalone(
        _ptr(p), p.size, _ptr(t), t.size,
        int(bool(conserve_volume)), int(bool(apply_com)), _ptr(radii))
    return RadiusGridResult(radii=radii, status=Status(status))


def radius_and_derivative(params: npt.ArrayLike, thetas: npt.ArrayLike,
                          conserve_volume: bool = False,
                          apply_com: bool = False) -> RadiusDerivativeResult:
    """R(theta) and dR/dtheta for one shape, a cache discarded per call.

    Parameters
    ----------
    params : array_like
        Beta parameters, 1 .. 64 entries.
    thetas : array_like
        Evaluation angles in radians; at least 2, none at a pole.
    conserve_volume : bool, optional
        Rescale radii to fixed volume.
    apply_com : bool, optional
        Apply the centre-of-mass correction.

    Returns
    -------
    RadiusDerivativeResult
        Radii, derivatives and status; both buffers zero-filled on failure.

    Notes
    -----
    ``params`` and ``thetas`` must be finite; see :func:`radius_grid`.
    """
    p = _as_1d(params, "params")
    t = _as_1d(thetas, "thetas")
    radii = np.zeros(t.size, dtype=np.float64)
    dr_dtheta = np.zeros(t.size, dtype=np.float64)
    status = _get_lib().beta_param_radius_and_derivative_standalone(
        _ptr(p), p.size, _ptr(t), t.size,
        int(bool(conserve_volume)), int(bool(apply_com)),
        _ptr(radii), _ptr(dr_dtheta))
    return RadiusDerivativeResult(
        radii=radii, dr_dtheta=dr_dtheta, status=Status(status))


class Cache:
    """Read-only tables over a fixed theta set, reused for any number of shapes.

    Create once, evaluate many shapes. Nothing parameter-dependent is stored:
    every method computes its result from ``params`` alone, so calls are
    independent and a ``Cache`` may be shared between threads. ``close()`` must
    not race with a compute.

    Parameters
    ----------
    max_params : int
        Longest ``params`` array the cache accepts, 1 .. ``MAX_BETA_PARAMS_LIMIT``.
        A shorter array is accepted; its missing trailing parameters are zero.
    thetas : array_like
        Evaluation angles in radians; at least 2, none at a pole. See
        :func:`theta_grid`.

    Raises
    ------
    BetaParamError
        If the library rejects the arguments; the message carries the status.

    Notes
    -----
    Every ``params`` array (and ``thetas``) must be finite. Non-finite input is
    undefined behavior: the library cannot detect NaN under fast-math, so
    compute calls return ``Status.valid`` with NaN outputs instead of an error.
    Screen inputs before calling.
    """

    def __init__(self, max_params: int, thetas: npt.ArrayLike) -> None:
        lib = _get_lib()
        t = _as_1d(thetas, "thetas")
        create_status = ctypes.c_int(0)
        handle = lib.beta_param_cache_create(
            int(max_params), _ptr(t), t.size, ctypes.byref(create_status))
        if not handle:
            raise BetaParamError(
                f"Cache creation failed: {status_message(create_status.value)} "
                f"(status {create_status.value})")
        self._handle: Optional[ctypes.c_void_p] = ctypes.c_void_p(handle)
        # Bound here so __del__ never needs module globals during shutdown.
        self._destroy = lib.beta_param_cache_destroy
        self.max_params = int(max_params)
        self.n_thetas = int(t.size)

    def close(self) -> None:
        """Destroy the underlying handle. Idempotent."""
        handle = getattr(self, "_handle", None)
        if handle is not None:
            self._handle = None
            self._destroy(handle)

    def __enter__(self) -> "Cache":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()

    def _require_handle(self) -> ctypes.c_void_p:
        handle = getattr(self, "_handle", None)
        if handle is None:
            raise BetaParamError("Cache is closed")
        return handle

    def radius_grid(self, params: npt.ArrayLike,
                    conserve_volume: bool = False,
                    apply_com: bool = False) -> RadiusGridResult:
        """R(theta) at the cache's thetas, validated.

        Parameters
        ----------
        params : array_like
            Beta parameters, 1 .. ``max_params`` entries.
        conserve_volume : bool, optional
            Rescale radii to fixed volume.
        apply_com : bool, optional
            Apply the centre-of-mass correction.

        Returns
        -------
        RadiusGridResult
            Radii and status; radii are zero-filled on failure.
        """
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = np.zeros(self.n_thetas, dtype=np.float64)
        status = _get_lib().beta_param_cache_radius_grid(
            handle, _ptr(p), p.size,
            int(bool(conserve_volume)), int(bool(apply_com)),
            _ptr(radii), radii.size)
        return RadiusGridResult(radii=radii, status=Status(status))

    def radius_and_derivative(self, params: npt.ArrayLike,
                              conserve_volume: bool = False,
                              apply_com: bool = False) -> RadiusDerivativeResult:
        """R(theta) and dR/dtheta at the cache's thetas.

        Parameters
        ----------
        params : array_like
            Beta parameters, 1 .. ``max_params`` entries.
        conserve_volume : bool, optional
            Rescale radii and derivatives to fixed volume.
        apply_com : bool, optional
            Apply the centre-of-mass correction.

        Returns
        -------
        RadiusDerivativeResult
            Radii, derivatives and status; both zero-filled on failure.
        """
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = np.zeros(self.n_thetas, dtype=np.float64)
        dr_dtheta = np.zeros(self.n_thetas, dtype=np.float64)
        status = _get_lib().beta_param_cache_radius_and_derivative(
            handle, _ptr(p), p.size,
            int(bool(conserve_volume)), int(bool(apply_com)),
            _ptr(radii), _ptr(dr_dtheta), radii.size)
        return RadiusDerivativeResult(
            radii=radii, dr_dtheta=dr_dtheta, status=Status(status))

    def radius_grid_unchecked(self, params: npt.ArrayLike,
                              apply_com: bool = False) -> RadiusGridResult:
        """R(theta) with no validation gates and no volume scaling.

        The rendering path: a shape the checked path rejects still yields its
        (partly negative) outline instead of a zero-filled buffer. Only usage
        errors and a failed COM correction are reported. The radii are never
        volume-scaled, which is why there is no ``conserve_volume`` option;
        scale by :attr:`ResolvedShape.volume_factor` if they must match
        :meth:`radius_grid`.

        Parameters
        ----------
        params : array_like
            Beta parameters, 1 .. ``max_params`` entries.
        apply_com : bool, optional
            Apply the centre-of-mass correction.

        Returns
        -------
        RadiusGridResult
            Unscaled radii and status.
        """
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = np.zeros(self.n_thetas, dtype=np.float64)
        status = _get_lib().beta_param_cache_radius_grid_unchecked(
            handle, _ptr(p), p.size, int(bool(apply_com)),
            _ptr(radii), radii.size)
        return RadiusGridResult(radii=radii, status=Status(status))

    def resolve_shape(self, params: npt.ArrayLike,
                      conserve_volume: bool = False,
                      apply_com: bool = False) -> ResolvedShape:
        """Resolve a shape without evaluating a grid.

        Parameters
        ----------
        params : array_like
            Beta parameters, 1 .. ``max_params`` entries.
        conserve_volume : bool, optional
            Compute the volume factor and scale the polar radii by it.
        apply_com : bool, optional
            Apply the centre-of-mass correction.

        Returns
        -------
        ResolvedShape
            COM-corrected beta10, the analytic polar radii, and the applied
            volume factor.
        """
        handle = self._require_handle()
        p = _as_1d(params, "params")
        corrected_beta10 = ctypes.c_double(0.0)
        r_north = ctypes.c_double(0.0)
        r_south = ctypes.c_double(0.0)
        volume_factor = ctypes.c_double(0.0)
        status = _get_lib().beta_param_cache_resolve_shape(
            handle, _ptr(p), p.size,
            int(bool(conserve_volume)), int(bool(apply_com)),
            ctypes.byref(corrected_beta10), ctypes.byref(r_north),
            ctypes.byref(r_south), ctypes.byref(volume_factor))
        return ResolvedShape(
            corrected_beta10=corrected_beta10.value,
            r_north=r_north.value, r_south=r_south.value,
            volume_factor=volume_factor.value, status=Status(status))


__all__ = [
    "BetaParamError", "Cache", "RadiusDerivativeResult", "RadiusGridResult",
    "ResolvedShape", "Status", "radius_and_derivative", "radius_grid",
    "status_message", "theta_grid",
    "MAX_BETA_PARAMS_LIMIT",
]
