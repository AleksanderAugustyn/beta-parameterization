/**
 * @file beta_parameterization.h
 * @brief C API for the Fortran beta parameterization library (v2.0.0).
 *
 * Two-tier API:
 *   - Cache-backed (recommended for batch use): create one cache up front,
 *     reuse it across many shape evaluations. The cache is immutable after
 *     creation and safe to share across threads.
 *   - Standalone single-shape (one-off, inefficient for batch): builds and
 *     discards a cache internally per call.
 *
 * All radius-grid functions return a status code; 0 (BETA_PARAM_VALID) means
 * success. The message buffer is always null-terminated (truncated to fit
 * if needed). Recommended buffer size: 256.
 *
 * Thread safety:
 *   A `beta_param_cache_t*` is immutable after `beta_param_cache_create`.
 *   Multiple threads may concurrently call any combination of the
 *   `beta_param_cache_compute_*` functions on the same cache. Only
 *   `beta_param_cache_create` and `beta_param_cache_destroy` require
 *   exclusive access.
 */

#ifndef BETA_PARAMETERIZATION_H
#define BETA_PARAMETERIZATION_H

#ifdef __cplusplus
extern "C" {
#endif

/* --- Limits --- */
#define BETA_PARAM_MAX_PARAMS_LIMIT 64

/* --- Status codes (mirror Fortran LEGENDRE_* parameters) --- */
#define BETA_PARAM_VALID                     0
#define BETA_PARAM_ERROR_NORTH_POLE          1
#define BETA_PARAM_ERROR_SOUTH_POLE          2
#define BETA_PARAM_ERROR_EMPTY_PARAMS        3
#define BETA_PARAM_ERROR_INTERIOR_NEGATIVE   4
#define BETA_PARAM_ERROR_INVALID_MAX_PARAMS  5
#define BETA_PARAM_ERROR_TOO_MANY_PARAMS     6
#define BETA_PARAM_ERROR_COM_NOT_CONVERGED   7
#define BETA_PARAM_ERROR_INVALID_BUFFER_SIZE 8

/* --- Opaque cache handle --- */
typedef struct beta_param_cache beta_param_cache_t;

/* --- Cache lifecycle --- */

/**
 * Create a cache. Returns NULL on failure (reason in message_buf).
 *
 * @param max_beta_params  1 .. BETA_PARAM_MAX_PARAMS_LIMIT
 * @param n_grid           Number of θ grid points (>= 2)
 * @param message_buf_len  Size of message_buf including null terminator
 * @param message_buf      Buffer to receive failure reason (empty on success)
 */
beta_param_cache_t* beta_param_cache_create(
        int max_beta_params, int n_grid,
        int message_buf_len, char* message_buf);

/** Destroy a cache. Null-safe. */
void beta_param_cache_destroy(beta_param_cache_t* cache);

/* --- Cache hot path --- */

/**
 * Compute R(θ) on the cache's θ grid.
 *
 * @param cache            Cache (must be non-NULL)
 * @param params           Deformation parameters β₁..βₙ
 * @param n_params         Number of deformation parameters
 * @param radii            Output buffer; must hold at least cache's n_grid doubles
 * @param message_buf_len  Size of message_buf
 * @param message_buf      Buffer to receive validation message (empty on success)
 * @return                 BETA_PARAM_VALID (0) on success, error code otherwise
 */
int beta_param_cache_compute_radius_grid(
        const beta_param_cache_t* cache,
        const double* params, int n_params,
        double* radii,
        int message_buf_len, char* message_buf);

/**
 * Compute R(θ) with COM shift applied to β₁₀ first.
 * `corrected_beta10` receives the post-shift β₁₀.
 */
int beta_param_cache_compute_radius_grid_with_com_shift(
        const beta_param_cache_t* cache,
        const double* params, int n_params,
        double* radii, double* corrected_beta10,
        int message_buf_len, char* message_buf);

/* --- Standalone single-shape (one-off, inefficient for batch) --- */

int beta_param_compute_radius_grid_standalone(
        const double* params, int n_params, int n_grid,
        double* radii,
        int message_buf_len, char* message_buf);

int beta_param_compute_radius_grid_standalone_with_com_shift(
        const double* params, int n_params, int n_grid,
        double* radii, double* corrected_beta10,
        int message_buf_len, char* message_buf);

#ifdef __cplusplus
}
#endif

#endif /* BETA_PARAMETERIZATION_H */
