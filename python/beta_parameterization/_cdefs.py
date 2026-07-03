"""ctypes signatures for the beta-parameterization C API."""
from __future__ import annotations

import ctypes

c_dbl_p = ctypes.POINTER(ctypes.c_double)

MAX_BETA_PARAMS_LIMIT = 64
MESSAGE_BUFFER_SIZE = 256


def configure(lib: ctypes.CDLL) -> ctypes.CDLL:
    """Set argtypes/restypes on the loaded library (idempotent)."""
    lib.beta_param_cache_create.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_char_p]
    lib.beta_param_cache_create.restype = ctypes.c_void_p

    lib.beta_param_cache_destroy.argtypes = [ctypes.c_void_p]
    lib.beta_param_cache_destroy.restype = None

    lib.beta_param_cache_compute_radius_grid.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int, c_dbl_p,
        ctypes.c_int, ctypes.c_char_p]
    lib.beta_param_cache_compute_radius_grid.restype = ctypes.c_int

    lib.beta_param_cache_compute_radius_grid_with_com_shift.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int, c_dbl_p, c_dbl_p,
        ctypes.c_int, ctypes.c_char_p]
    lib.beta_param_cache_compute_radius_grid_with_com_shift.restype = ctypes.c_int

    lib.beta_param_compute_radius_grid_standalone.argtypes = [
        c_dbl_p, ctypes.c_int, ctypes.c_int, c_dbl_p,
        ctypes.c_int, ctypes.c_char_p]
    lib.beta_param_compute_radius_grid_standalone.restype = ctypes.c_int

    lib.beta_param_compute_radius_grid_standalone_with_com_shift.argtypes = [
        c_dbl_p, ctypes.c_int, ctypes.c_int, c_dbl_p, c_dbl_p,
        ctypes.c_int, ctypes.c_char_p]
    lib.beta_param_compute_radius_grid_standalone_with_com_shift.restype = ctypes.c_int

    lib.beta_param_node_set_create.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int, ctypes.c_int, ctypes.c_char_p]
    lib.beta_param_node_set_create.restype = ctypes.c_void_p

    lib.beta_param_node_set_destroy.argtypes = [ctypes.c_void_p]
    lib.beta_param_node_set_destroy.restype = None

    lib.beta_param_cache_resolve_shape.argtypes = [
        ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p,
        ctypes.c_int, ctypes.c_char_p]
    lib.beta_param_cache_resolve_shape.restype = ctypes.c_int

    lib.beta_param_cache_compute_radius_and_derivative.argtypes = [
        ctypes.c_void_p, ctypes.c_void_p, c_dbl_p, ctypes.c_int,
        c_dbl_p, c_dbl_p, ctypes.c_int,
        ctypes.c_int, ctypes.c_char_p]
    lib.beta_param_cache_compute_radius_and_derivative.restype = ctypes.c_int
    return lib
