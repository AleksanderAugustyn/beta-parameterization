/**
 * @file beta_parameterization.h
 * @brief C API for the Fortran beta parameterization library (v4.0.0).
 *
 * Two tiers of computation:
 *   - One-shot: the `*_standalone` functions compute one shape per call. All
 *     workspace is internal and discarded on return.
 *   - Read-only cache: create a `beta_param_cache_t` once, then call the
 *     `beta_param_cache_*` functions for any number of shapes.
 *
 * Two handle types:
 *   - `beta_param_cache_t`    — everything determined by `max_params` and the
 *     primary theta set (normalization constants, Gauss-Legendre tables,
 *     Legendre tables at the primary thetas).
 *   - `beta_param_node_set_t` — an extra set of evaluation thetas with its own
 *     Legendre tables, built from a cache.
 *
 * Thread safety:
 *   Both handles are immutable after creation. Every compute function takes
 *   them `const` and may be called concurrently, from any number of threads,
 *   on the same handles. Create and destroy must not race with any other call
 *   on the same handle. The library holds no mutable global state.
 *
 * Lifetime:
 *   A cache must outlive every node set built from it. Destroy order: node
 *   sets first, their cache last. A node set serves a cache only if it was
 *   built from a cache with at least that cache's `max_params`.
 *
 * Parameters:
 *   A cached call accepts 1 .. max_params parameters; a one-shot call accepts
 *   1 .. BETA_PARAM_MAX_PARAMS_LIMIT. Missing trailing parameters are zero: a
 *   short vector and its zero-padded form give bitwise-identical results.
 *   `conserve_volume` and `apply_com` are per-call options (nonzero = on); one
 *   cache serves every combination.
 *
 * Diagnostics:
 *   There are no message buffers. `_create` functions return NULL on failure
 *   and write the reason to the nullable `int* status` out-parameter (pass
 *   NULL to ignore it). Every other function returns the status code directly;
 *   BETA_PARAM_VALID (0) means success. beta_param_status_message() maps a code
 *   to a fixed, static, null-terminated string; the returned pointer is owned
 *   by the library, never freed by the caller, and is safe to read from any
 *   thread.
 *
 * Failure behavior:
 *   On any nonzero status from a compute function, every output buffer is
 *   zero-filled. No state exists, so the next call is unaffected. A NULL
 *   handle passed into a compute function returns
 *   BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED (2).
 *
 * Size arguments:
 *   Every size argument must be the ACTUAL extent of the caller's buffer. A
 *   wrong size that is honest about the caller's own memory is reported with a
 *   status code, however large; a stated size larger than the real buffer is
 *   undefined behavior. Negative counts are treated as zero.
 *
 * Precondition — finite input:
 *   `params` and `thetas` must be finite. Non-finite input is undefined
 *   behavior: the library cannot detect NaN under fast-math, so no check
 *   rejects it, and no particular result is promised. A call may return
 *   BETA_PARAM_VALID (0) with NaN outputs; a NaN trailing parameter is
 *   trimmed like a zero, giving the finite outputs of the shorter vector; a
 *   Debug build may trap. Screen inputs before calling. Magnitudes must be
 *   physical as well: with `apply_com` the centre-of-mass quadrature evaluates
 *   R^4 before any validity gate and overflows for |beta| beyond about 1e70.
 */

#ifndef BETA_PARAMETERIZATION_H
#define BETA_PARAMETERIZATION_H

#ifdef __cplusplus
extern "C" {
#endif

/* --- Limits --- */
/** Longest parameter vector either tier accepts; the highest `max_params`. */
#define BETA_PARAM_MAX_PARAMS_LIMIT 64

/* --- Shared contract status codes (0-99, identical numbers in every
 *     shape-parameterization library; 6 is retired and never reused) --- */
#define BETA_PARAM_VALID                        0
#define BETA_PARAM_ERROR_TOO_MANY_PARAMS        1
#define BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED  2
#define BETA_PARAM_ERROR_INVALID_GRID           3
#define BETA_PARAM_ERROR_WRONG_PARAM_COUNT      4
#define BETA_PARAM_ERROR_INVALID_INIT           5

/* --- Library status codes (>= 100, append-only) --- */
#define BETA_PARAM_ERROR_NORTH_POLE          100
#define BETA_PARAM_ERROR_SOUTH_POLE          101
#define BETA_PARAM_ERROR_INTERIOR_NEGATIVE   102
#define BETA_PARAM_ERROR_COM_NOT_CONVERGED   103
#define BETA_PARAM_ERROR_INVALID_BUFFER_SIZE 104
#define BETA_PARAM_ERROR_POLE_NODE           105
#define BETA_PARAM_ERROR_NODE_SET_MISMATCH   106

/* --- Opaque handles --- */
typedef struct beta_param_cache beta_param_cache_t;
typedef struct beta_param_node_set beta_param_node_set_t;

/* --- Diagnostics --- */

/**
 * Fixed description of a status code. Never NULL; unknown codes map to an
 * "unknown status code" string. The pointer is to static storage: do not free
 * it, and it stays valid for the life of the process.
 */
const char* beta_param_status_message(int status);

/* --- Cache lifecycle --- */

/**
 * Build the read-only cache. Returns NULL on failure.
 *
 * @param max_params  Longest parameter vector the cache accepts,
 *                    1 .. BETA_PARAM_MAX_PARAMS_LIMIT; also the Legendre order
 *                    of its tables. Below 1: BETA_PARAM_ERROR_INVALID_INIT;
 *                    above the limit: BETA_PARAM_ERROR_TOO_MANY_PARAMS.
 * @param thetas      Primary theta set in radians; at least 2 entries
 *                    (BETA_PARAM_ERROR_INVALID_GRID), none at or beyond a pole
 *                    (cos(theta)^2 == 1 in double precision;
 *                    BETA_PARAM_ERROR_POLE_NODE)
 * @param n_thetas    Number of entries in thetas
 * @param status      Nullable; receives BETA_PARAM_VALID or the rejecting code
 */
beta_param_cache_t* beta_param_cache_create(
        int max_params, const double* thetas, int n_thetas, int* status);

/** Destroy a cache. NULL-safe. Every node set built from it must already be
 *  destroyed. */
void beta_param_cache_destroy(beta_param_cache_t* cache);

/* --- Node-set lifecycle --- */

/**
 * Build an extra evaluation set (thetas plus Legendre P_k and P_k' tables)
 * sized to `cache`. Returns NULL on failure — including a pole node, which is
 * rejected with BETA_PARAM_ERROR_POLE_NODE; use the resolve_shape polar radii
 * for the poles instead. A NULL cache gives
 * BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED.
 *
 * @param cache     Cache handle; supplies max_params. Must outlive the node
 *                  set.
 * @param thetas    Node angles in radians; any order, need not be uniform
 * @param n_thetas  Number of nodes (at least 2)
 * @param status    Nullable; receives the rejecting code on failure
 */
beta_param_node_set_t* beta_param_node_set_create(
        const beta_param_cache_t* cache, const double* thetas, int n_thetas,
        int* status);

/** Destroy a node set. NULL-safe. */
void beta_param_node_set_destroy(beta_param_node_set_t* node_set);

/* --- Cached computes ---
 *
 * None of these modifies the cache. `params` holds 1 .. max_params entries,
 * else BETA_PARAM_ERROR_WRONG_PARAM_COUNT. Output buffer lengths must equal
 * the cache's theta count (or the node set's node count), else
 * BETA_PARAM_ERROR_INVALID_BUFFER_SIZE. `conserve_volume` and `apply_com` are
 * nonzero to switch the option on.
 */

/** R(theta) at the cache's primary thetas. `n_radii` must equal that count. */
int beta_param_cache_radius_grid(
        const beta_param_cache_t* cache, const double* params, int n_params,
        int conserve_volume, int apply_com, double* radii, int n_radii);

/** R(theta) and dR/dtheta at the cache's primary thetas. */
int beta_param_cache_radius_and_derivative(
        const beta_param_cache_t* cache, const double* params, int n_params,
        int conserve_volume, int apply_com,
        double* radii, double* dr_dthetas, int n_radii);

/**
 * R(theta) at the primary thetas with NO validation gates and NO volume
 * scaling — the rendering/diagnostic path. It reports usage errors and a
 * failed COM correction only, so a shape rejected by the checked path still
 * yields its (partly negative) outline instead of a zero-filled buffer.
 *
 * CAUTION: the radii are UNSCALED — this path never computes the volume
 * factor, which is why it takes no `conserve_volume`. Mixing it with
 * beta_param_cache_radius_grid(..., conserve_volume = 1, ...) draws two
 * outlines of different size for one shape; scale by
 * beta_param_cache_resolve_shape()'s volume_factor if the sizes must agree.
 */
int beta_param_cache_radius_grid_unchecked(
        const beta_param_cache_t* cache, const double* params, int n_params,
        int apply_com, double* radii, int n_radii);

/**
 * Resolve a shape without evaluating a grid: the COM-corrected beta10, the
 * analytic polar radii, and the applied volume factor.
 *
 * `corrected_beta10` is a beta-space value and is never volume-scaled;
 * `r_north` and `r_south` are scaled. `volume_factor` is exactly 1.0 when
 * conserve_volume = 0.
 */
int beta_param_cache_resolve_shape(
        const beta_param_cache_t* cache, const double* params, int n_params,
        int conserve_volume, int apply_com,
        double* corrected_beta10, double* r_north, double* r_south,
        double* volume_factor);

/**
 * R(theta) and dR/dtheta at a node set's thetas. The node set must come from
 * a cache with at least this cache's max_params, else
 * BETA_PARAM_ERROR_NODE_SET_MISMATCH; a node set built from this cache always
 * qualifies. `n_nodes` must equal the node set's node count. To fold several
 * theta grids into one call, build one node set over their concatenation.
 */
int beta_param_cache_node_radius_and_derivative(
        const beta_param_cache_t* cache, const double* params, int n_params,
        const beta_param_node_set_t* node_set,
        int conserve_volume, int apply_com,
        double* radii, double* dr_dthetas, int n_nodes);

/* --- One-shot computes (tier 1) ---
 *
 * A cache is built, used and discarded per call — no handle, nothing to free.
 * `n_params` may be 1 .. BETA_PARAM_MAX_PARAMS_LIMIT: above that,
 * BETA_PARAM_ERROR_TOO_MANY_PARAMS, never silent truncation; below 1,
 * BETA_PARAM_ERROR_WRONG_PARAM_COUNT. Output buffers hold n_thetas doubles.
 * Outputs are zero-filled on failure. Results are bitwise identical to the
 * cached functions on a cache with max_params >= n_params.
 */

int beta_param_radius_grid_standalone(
        const double* params, int n_params,
        const double* thetas, int n_thetas,
        int conserve_volume, int apply_com, double* radii);

int beta_param_radius_and_derivative_standalone(
        const double* params, int n_params,
        const double* thetas, int n_thetas,
        int conserve_volume, int apply_com,
        double* radii, double* dr_dthetas);

#ifdef __cplusplus
}
#endif

#endif /* BETA_PARAMETERIZATION_H */
