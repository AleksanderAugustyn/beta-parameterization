"""ctypes signatures for the 4.0.0 beta-parameterization C API.

Mirrors ``include/beta_parameterization.h``: opaque handles, ``_create``
functions returning NULL plus a status out-parameter, every other function
returning its status directly.
"""
from __future__ import annotations

import ctypes

c_dbl_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int)

#: Longest parameter vector either tier accepts, and the highest ``max_params``
#: (BETA_PARAM_MAX_PARAMS_LIMIT).
MAX_BETA_PARAMS_LIMIT = 64


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

    # --- cache lifecycle ---
    lib.beta_param_cache_create.argtypes = [
        ctypes.c_int, c_dbl_p, ctypes.c_int, c_int_p]
    lib.beta_param_cache_create.restype = ctypes.c_void_p

    lib.beta_param_cache_destroy.argtypes = [ctypes.c_void_p]
    lib.beta_param_cache_destroy.restype = None

    # --- node-set lifecycle (declared for completeness; not exposed in Python) ---
    lib.beta_param_node_set_create.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int, c_int_p]
    lib.beta_param_node_set_create.restype = ctypes.c_void_p

    lib.beta_param_node_set_destroy.argtypes = [ctypes.c_void_p]
    lib.beta_param_node_set_destroy.restype = None

    # --- cached computes: (cache, params, n_params, options..., outputs...) ---
    lib.beta_param_cache_radius_grid.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, c_dbl_p, ctypes.c_int]
    lib.beta_param_cache_radius_grid.restype = ctypes.c_int

    lib.beta_param_cache_radius_and_derivative.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, c_dbl_p, c_dbl_p, ctypes.c_int]
    lib.beta_param_cache_radius_and_derivative.restype = ctypes.c_int

    lib.beta_param_cache_radius_grid_unchecked.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        ctypes.c_int, c_dbl_p, ctypes.c_int]
    lib.beta_param_cache_radius_grid_unchecked.restype = ctypes.c_int

    lib.beta_param_cache_resolve_shape.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int,
        c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p]
    lib.beta_param_cache_resolve_shape.restype = ctypes.c_int

    lib.beta_param_cache_node_radius_and_derivative.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int, ctypes.c_void_p,
        ctypes.c_int, ctypes.c_int, c_dbl_p, c_dbl_p, ctypes.c_int]
    lib.beta_param_cache_node_radius_and_derivative.restype = ctypes.c_int

    # --- one-shot computes ---
    lib.beta_param_radius_grid_standalone.argtypes = [
        c_dbl_p, ctypes.c_int, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, c_dbl_p]
    lib.beta_param_radius_grid_standalone.restype = ctypes.c_int

    lib.beta_param_radius_and_derivative_standalone.argtypes = [
        c_dbl_p, ctypes.c_int, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, c_dbl_p, c_dbl_p]
    lib.beta_param_radius_and_derivative_standalone.restype = ctypes.c_int
    return lib
