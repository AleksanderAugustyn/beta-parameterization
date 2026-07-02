"""High-level API: Status, RadiusGridResult, Cache, standalone functions."""
from __future__ import annotations

import ctypes
from dataclasses import dataclass
from enum import IntEnum
from typing import Optional

import numpy as np
import numpy.typing as npt

from ._cdefs import MAX_BETA_PARAMS_LIMIT, MESSAGE_BUFFER_SIZE, c_dbl_p, configure
from ._libloader import load_library

_lib: Optional[ctypes.CDLL] = None


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is None:
        _lib = configure(load_library())
    return _lib


class BetaParamError(RuntimeError):
    """Raised for lifecycle errors (init failure, use-after-close)."""


class Status(IntEnum):
    """Mirrors BETA_PARAM_* codes in include/beta_parameterization.h."""
    VALID = 0
    ERROR_NORTH_POLE = 1
    ERROR_SOUTH_POLE = 2
    ERROR_EMPTY_PARAMS = 3
    ERROR_INTERIOR_NEGATIVE = 4
    ERROR_INVALID_MAX_PARAMS = 5
    ERROR_TOO_MANY_PARAMS = 6
    ERROR_COM_NOT_CONVERGED = 7
    ERROR_INVALID_BUFFER_SIZE = 8


@dataclass(frozen=True)
class RadiusGridResult:
    """One radius-grid evaluation. Invalid shapes are normal flow, not exceptions."""
    radii: npt.NDArray[np.float64]
    status: Status
    message: str
    corrected_beta10: Optional[float] = None

    @property
    def ok(self) -> bool:
        return self.status == Status.VALID


def theta_grid(n_grid: int) -> npt.NDArray[np.float64]:
    """The θ grid the library evaluates on: n_grid uniform points over [0, π]."""
    return np.linspace(0.0, np.pi, n_grid)


def _as_params(params: npt.ArrayLike) -> npt.NDArray[np.float64]:
    arr = np.ascontiguousarray(params, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"params must be 1-D, got shape {arr.shape}")
    return arr


class Cache:
    """Reusable precomputed tables; mirrors beta_param::Cache in the C++ wrapper.

    Construction raises BetaParamError on invalid arguments; the hot path
    returns RadiusGridResult with a Status instead of raising, so parameter
    sweeps can render invalid shapes. Safe to share across threads.
    """

    def __init__(self, max_beta_params: int, n_grid: int) -> None:
        lib = _get_lib()
        buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
        handle = lib.beta_param_cache_create(
            int(max_beta_params), int(n_grid), MESSAGE_BUFFER_SIZE, buf)
        if not handle:
            raise BetaParamError(
                f"Cache init failed: {buf.value.decode(errors='replace')}")
        self._handle: Optional[ctypes.c_void_p] = ctypes.c_void_p(handle)
        self.max_beta_params = int(max_beta_params)
        self.n_grid = int(n_grid)

    def close(self) -> None:
        if getattr(self, "_handle", None) is not None:
            _get_lib().beta_param_cache_destroy(self._handle)
            self._handle = None

    def __enter__(self) -> "Cache":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()

    def _require_handle(self) -> ctypes.c_void_p:
        if self._handle is None:
            raise BetaParamError("Cache is closed")
        return self._handle

    def radius_grid(self, params: npt.ArrayLike) -> RadiusGridResult:
        """R(θ) on the cache grid, no COM shift."""
        handle = self._require_handle()
        arr = _as_params(params)
        radii = np.zeros(self.n_grid, dtype=np.float64)
        buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
        status = _get_lib().beta_param_cache_compute_radius_grid(
            handle,
            arr.ctypes.data_as(c_dbl_p), arr.size,
            radii.ctypes.data_as(c_dbl_p),
            MESSAGE_BUFFER_SIZE, buf)
        return RadiusGridResult(
            radii=radii, status=Status(status),
            message=buf.value.decode(errors="replace"))

    def radius_grid_with_com_shift(self, params: npt.ArrayLike) -> RadiusGridResult:
        """R(θ) with the COM shift applied to β₁₀ first."""
        handle = self._require_handle()
        arr = _as_params(params)
        radii = np.zeros(self.n_grid, dtype=np.float64)
        corrected = ctypes.c_double(0.0)
        buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
        status = _get_lib().beta_param_cache_compute_radius_grid_with_com_shift(
            handle,
            arr.ctypes.data_as(c_dbl_p), arr.size,
            radii.ctypes.data_as(c_dbl_p), ctypes.byref(corrected),
            MESSAGE_BUFFER_SIZE, buf)
        return RadiusGridResult(
            radii=radii, status=Status(status),
            message=buf.value.decode(errors="replace"),
            corrected_beta10=corrected.value)


def radius_grid_standalone(params: npt.ArrayLike, n_grid: int) -> RadiusGridResult:
    """One-off evaluation (builds and discards a cache internally)."""
    arr = _as_params(params)
    radii = np.zeros(int(n_grid), dtype=np.float64)
    buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
    status = _get_lib().beta_param_compute_radius_grid_standalone(
        arr.ctypes.data_as(c_dbl_p), arr.size, int(n_grid),
        radii.ctypes.data_as(c_dbl_p),
        MESSAGE_BUFFER_SIZE, buf)
    return RadiusGridResult(
        radii=radii, status=Status(status),
        message=buf.value.decode(errors="replace"))


def radius_grid_standalone_with_com_shift(
        params: npt.ArrayLike, n_grid: int) -> RadiusGridResult:
    """One-off evaluation with COM shift."""
    arr = _as_params(params)
    radii = np.zeros(int(n_grid), dtype=np.float64)
    corrected = ctypes.c_double(0.0)
    buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
    status = _get_lib().beta_param_compute_radius_grid_standalone_with_com_shift(
        arr.ctypes.data_as(c_dbl_p), arr.size, int(n_grid),
        radii.ctypes.data_as(c_dbl_p), ctypes.byref(corrected),
        MESSAGE_BUFFER_SIZE, buf)
    return RadiusGridResult(
        radii=radii, status=Status(status),
        message=buf.value.decode(errors="replace"),
        corrected_beta10=corrected.value)
