/**
 * @file beta_parameterization.h
 * @brief C API for the Fortran beta parameterization library (v2.1.0).
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
 * Status-code notes:
 *   BETA_PARAM_ERROR_INVALID_MAX_PARAMS (5) covers every invalid cache-init
 *   parameter (max_beta_params outside [1, 64] or n_grid < 2), a NULL cache
 *   handle passed to a compute function, and - through the standalone entry
 *   points, where max_beta_params is taken as n_params - an empty params
 *   array. On any failure, output buffers (radii, corrected_beta10) are
 *   zero-filled.
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
#define BETA_PARAM_ERROR_POLE_NODE           9

/* --- Opaque cache handle --- */
typedef struct beta_param_cache beta_param_cache_t;

/* --- Opaque node-set handle (precomputed Legendre tables at fixed thetas) --- */
typedef struct beta_param_node_set beta_param_node_set_t;

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

/* --- Node-set API (arbitrary thetas, R + analytic dR/dtheta) --- */

/**
 * Build a node set: precomputed Legendre P_k and P_k' tables at the given
 * thetas, sized to the cache's max_beta_params. Returns NULL on failure
 * (reason in message_buf) — including any pole node (theta where
 * cos(theta)^2 == 1 in double precision), which is rejected with
 * BETA_PARAM_ERROR_POLE_NODE; use the resolve_shape polar radii instead.
 *
 * Thread safety: a node set is immutable after create; multiple threads may
 * concurrently pass it to beta_param_cache_compute_radius_and_derivative.
 *
 * @param cache            Cache the tables are sized to (must be non-NULL)
 * @param thetas           Node angles in radians (need not be uniform)
 * @param n_thetas         Number of nodes (>= 1)
 * @param message_buf_len  Size of message_buf including null terminator
 * @param message_buf      Buffer to receive failure reason (empty on success)
 */
beta_param_node_set_t* beta_param_node_set_create(
        const beta_param_cache_t* cache, const double* thetas, int n_thetas,
        int message_buf_len, char* message_buf);

/** Destroy a node set. Null-safe. */
void beta_param_node_set_destroy(beta_param_node_set_t* node_set);

/**
 * Resolve a shape once: pad + normalize params, run the COM iteration,
 * polar pre-check. Node-set-independent — feed the resulting beta_con to
 * beta_param_cache_compute_radius_and_derivative for any number of node sets.
 * On failure all outputs are zero-filled.
 *
 * @param cache             Cache (must be non-NULL)
 * @param params            Deformation parameters (beta10 is the COM dipole)
 * @param n_params          Number of deformation parameters
 * @param beta_con          Output; must hold at least the cache's max_beta_params doubles
 * @param corrected_beta10  Receives the post-shift beta10
 * @param r_north           Receives analytic R(0)
 * @param r_south           Receives analytic R(pi)
 * @param message_buf_len   Size of message_buf
 * @param message_buf       Buffer to receive validation message (empty on success)
 * @return                  BETA_PARAM_VALID (0) on success, error code otherwise
 */
int beta_param_cache_resolve_shape(
        const beta_param_cache_t* cache, const double* params, int n_params,
        double* beta_con, double* corrected_beta10, double* r_north, double* r_south,
        int message_buf_len, char* message_buf);

/**
 * Evaluate R and dR/dtheta at a node set for a shape already resolved by
 * beta_param_cache_resolve_shape. Dot products only — no iteration.
 * On failure radii and dr_dtheta are zero-filled.
 *
 * @param cache            Cache (must be non-NULL)
 * @param node_set         Node set built by this cache (must be non-NULL)
 * @param beta_con         Resolved coefficients; n_beta_con must equal the
 *                         cache's max_beta_params
 * @param n_beta_con       Number of beta_con entries
 * @param radii            Output R(theta_i); n_nodes doubles
 * @param dr_dtheta        Output dR/dtheta(theta_i); n_nodes doubles
 * @param n_nodes          Must equal the node set's node count
 * @param message_buf_len  Size of message_buf
 * @param message_buf      Buffer to receive validation message (empty on success)
 * @return                 BETA_PARAM_VALID (0) on success, error code otherwise
 */
int beta_param_cache_compute_radius_and_derivative(
        const beta_param_cache_t* cache, const beta_param_node_set_t* node_set,
        const double* beta_con, int n_beta_con,
        double* radii, double* dr_dtheta, int n_nodes,
        int message_buf_len, char* message_buf);

/* --- Standalone single-shape (one-off, inefficient for batch) ---
 *
 * Return BETA_PARAM_VALID (0) on success, an error code otherwise.
 * max_beta_params is taken as n_params, so n_params must be in [1, 64]
 * and n_grid >= 2; violations return BETA_PARAM_ERROR_INVALID_MAX_PARAMS.
 * On failure, radii (and corrected_beta10) are zero-filled.
 */

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
