"""ctypes signatures for the 3.0.0 beta-parameterization C API.

Mirrors ``include/beta_parameterization.h``: opaque handles, ``_create``
functions returning NULL plus a status out-parameter, every other function
returning its status directly.
"""
from __future__ import annotations

import ctypes

c_dbl_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int)

#: Highest ``n_params`` the standalone tier accepts (BETA_PARAM_MAX_PARAMS_LIMIT).
MAX_BETA_PARAMS_LIMIT = 64
#: Highest ``n_params`` a cache accepts (BETA_PARAM_CACHE_MAX_PARAMS).
CACHE_MAX_PARAMS = 8


def configure(lib: ctypes.CDLL) -> ctypes.CDLL:
    """Set argtypes/restypes on the loaded library (idempotent).

    Parameters
    ----------
    lib : ctypes.CDLL
        Freshly loaded ``libbeta_parameterization.so``.

    Returns
    -------
    ctypes.CDLL
        The same handle, with every exported symbol declared.
    """
    # --- diagnostics ---
    lib.beta_param_status_message.argtypes = [ctypes.c_int]
    lib.beta_param_status_message.restype = ctypes.c_char_p

    # --- tables lifecycle ---
    lib.beta_param_tables_create.argtypes = [
        ctypes.c_int, c_dbl_p, ctypes.c_int, c_int_p]
    lib.beta_param_tables_create.restype = ctypes.c_void_p

    lib.beta_param_tables_destroy.argtypes = [ctypes.c_void_p]
    lib.beta_param_tables_destroy.restype = None

    # --- cache lifecycle ---
    lib.beta_param_cache_create.argtypes = [
        ctypes.c_int, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, c_int_p]
    lib.beta_param_cache_create.restype = ctypes.c_void_p

    lib.beta_param_cache_create_shared.argtypes = [
        ctypes.c_void_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, c_int_p]
    lib.beta_param_cache_create_shared.restype = ctypes.c_void_p

    lib.beta_param_cache_destroy.argtypes = [ctypes.c_void_p]
    lib.beta_param_cache_destroy.restype = None

    # --- node-set lifecycle (declared for completeness; not exposed in Python) ---
    lib.beta_param_node_set_create.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int, c_int_p]
    lib.beta_param_node_set_create.restype = ctypes.c_void_p

    lib.beta_param_node_set_destroy.argtypes = [ctypes.c_void_p]
    lib.beta_param_node_set_destroy.restype = None

    # --- cached computes ---
    lib.beta_param_cache_radius_grid.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int, c_dbl_p, ctypes.c_int]
    lib.beta_param_cache_radius_grid.restype = ctypes.c_int

    lib.beta_param_cache_radius_and_derivative.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        c_dbl_p, c_dbl_p, ctypes.c_int]
    lib.beta_param_cache_radius_and_derivative.restype = ctypes.c_int

    lib.beta_param_cache_radius_grid_unchecked.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int, c_dbl_p, ctypes.c_int]
    lib.beta_param_cache_radius_grid_unchecked.restype = ctypes.c_int

    lib.beta_param_cache_resolve_shape.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p]
    lib.beta_param_cache_resolve_shape.restype = ctypes.c_int

    lib.beta_param_cache_node_radius_and_derivative.argtypes = [
        ctypes.c_void_p, ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        c_dbl_p, c_dbl_p, ctypes.c_int]
    lib.beta_param_cache_node_radius_and_derivative.restype = ctypes.c_int

    # --- standalone computes ---
    lib.beta_param_radius_grid_standalone.argtypes = [
        c_dbl_p, ctypes.c_int, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, c_dbl_p]
    lib.beta_param_radius_grid_standalone.restype = ctypes.c_int

    lib.beta_param_radius_and_derivative_standalone.argtypes = [
        c_dbl_p, ctypes.c_int, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, c_dbl_p, c_dbl_p]
    lib.beta_param_radius_and_derivative_standalone.restype = ctypes.c_int
    return lib
