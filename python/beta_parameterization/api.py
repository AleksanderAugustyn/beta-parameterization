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
    POLE_NODE = 9
    ERROR_NO_UNIFORM_GRID = 10


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


@dataclass(frozen=True)
class ResolvedShape:
    """One resolve_shape call: COM-corrected coefficients plus analytic pole radii."""
    beta_con: npt.NDArray[np.float64]
    corrected_beta10: float
    r_north: float
    r_south: float
    status: Status
    message: str

    @property
    def ok(self) -> bool:
        return self.status == Status.VALID


@dataclass(frozen=True)
class RadiusDerivativeResult:
    """One radius_and_derivative call: R and dR/dθ at a node set's thetas."""
    radii: npt.NDArray[np.float64]
    dr_dtheta: npt.NDArray[np.float64]
    status: Status
    message: str

    @property
    def ok(self) -> bool:
        return self.status == Status.VALID


class NodeSet:
    """Precomputed Legendre tables at fixed thetas. Create via Cache.build_node_set."""

    def __init__(self, handle: ctypes.c_void_p, n_nodes: int) -> None:
        self._handle: Optional[ctypes.c_void_p] = handle
        self.n_nodes = n_nodes

    def close(self) -> None:
        if getattr(self, "_handle", None) is not None:
            _get_lib().beta_param_node_set_destroy(self._handle)
            self._handle = None

    def __enter__(self) -> "NodeSet":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()

    def _require_handle(self) -> ctypes.c_void_p:
        if self._handle is None:
            raise BetaParamError("NodeSet is closed")
        return self._handle


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

    n_grid=None builds a node-set-only cache; the uniform-grid entry points
    then return Status.ERROR_NO_UNIFORM_GRID.
    """

    def __init__(self, max_beta_params: int, n_grid: int | None = None) -> None:
        lib = _get_lib()
        buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
        c_n_grid = 0 if n_grid is None else int(n_grid)
        handle = lib.beta_param_cache_create(
            int(max_beta_params), c_n_grid, MESSAGE_BUFFER_SIZE, buf)
        if not handle:
            raise BetaParamError(
                f"Cache init failed: {buf.value.decode(errors='replace')}")
        self._handle: Optional[ctypes.c_void_p] = ctypes.c_void_p(handle)
        self.max_beta_params = int(max_beta_params)
        self.n_grid = c_n_grid

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

    def build_node_set(self, thetas: npt.ArrayLike) -> NodeSet:
        """Precompute P_k and P_k' tables at the given thetas (radians).

        Raises BetaParamError on failure — including any pole node
        (Status.POLE_NODE); pole radii come from resolve_shape instead.
        """
        handle = self._require_handle()
        arr = np.ascontiguousarray(thetas, dtype=np.float64)
        if arr.ndim != 1:
            raise ValueError(f"thetas must be 1-D, got shape {arr.shape}")
        buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
        ns_handle = _get_lib().beta_param_node_set_create(
            handle,
            arr.ctypes.data_as(c_dbl_p), arr.size,
            MESSAGE_BUFFER_SIZE, buf)
        if not ns_handle:
            raise BetaParamError(
                f"NodeSet build failed: {buf.value.decode(errors='replace')}")
        return NodeSet(ctypes.c_void_p(ns_handle), int(arr.size))

    def resolve_shape(self, params: npt.ArrayLike,
                      apply_com_correction: bool = True) -> ResolvedShape:
        """Resolve a shape once (COM iteration + polar pre-check).

        Node-set-independent: feed the returned beta_con to
        radius_and_derivative for any number of node sets.

        apply_com_correction=False skips the COM iteration
        (corrected_beta10 = input beta10), matching radius_grid's
        no-COM semantics.
        """
        handle = self._require_handle()
        arr = _as_params(params)
        beta_con = np.zeros(self.max_beta_params, dtype=np.float64)
        corrected = ctypes.c_double(0.0)
        r_north = ctypes.c_double(0.0)
        r_south = ctypes.c_double(0.0)
        buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
        status = _get_lib().beta_param_cache_resolve_shape(
            handle,
            arr.ctypes.data_as(c_dbl_p), arr.size,
            beta_con.ctypes.data_as(c_dbl_p),
            ctypes.byref(corrected), ctypes.byref(r_north), ctypes.byref(r_south),
            1 if apply_com_correction else 0,
            MESSAGE_BUFFER_SIZE, buf)
        return ResolvedShape(
            beta_con=beta_con, corrected_beta10=corrected.value,
            r_north=r_north.value, r_south=r_south.value,
            status=Status(status), message=buf.value.decode(errors="replace"))

    def radius_and_derivative(
            self, beta_con: npt.ArrayLike, node_set: NodeSet) -> RadiusDerivativeResult:
        """R and dR/dθ at a node set, for a shape resolved by resolve_shape."""
        handle = self._require_handle()
        ns_handle = node_set._require_handle()
        arr = np.ascontiguousarray(beta_con, dtype=np.float64)
        radii = np.zeros(node_set.n_nodes, dtype=np.float64)
        dr_dtheta = np.zeros(node_set.n_nodes, dtype=np.float64)
        buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
        status = _get_lib().beta_param_cache_compute_radius_and_derivative(
            handle, ns_handle,
            arr.ctypes.data_as(c_dbl_p), arr.size,
            radii.ctypes.data_as(c_dbl_p), dr_dtheta.ctypes.data_as(c_dbl_p),
            node_set.n_nodes,
            MESSAGE_BUFFER_SIZE, buf)
        return RadiusDerivativeResult(
            radii=radii, dr_dtheta=dr_dtheta,
            status=Status(status), message=buf.value.decode(errors="replace"))


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
