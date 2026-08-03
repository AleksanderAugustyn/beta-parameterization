/**
 * @file beta_parameterization.h
 * @brief C API for the Fortran beta parameterization library (v3.0.0).
 *
 * Three handle types, from shared to per-shape:
 *   - `beta_param_tables_t`   — everything determined by `max_l` and the
 *     primary theta set (normalization constants, Gauss-Legendre tables,
 *     Legendre tables at the primary thetas). Built once, shared read-only.
 *   - `beta_param_node_set_t` — an extra set of evaluation thetas with its own
 *     Legendre tables, sized to one tables handle. Immutable after creation.
 *   - `beta_param_cache_t`    — per-shape working state: the recompute engine
 *     and every per-shape buffer. Mutated by every compute call.
 *
 * Two tiers of computation:
 *   - Cached: create a cache once, call the `beta_param_cache_*` functions for
 *     many shapes. Only the intermediates invalidated by the changed
 *     parameters are recomputed.
 *   - Standalone: one-off calls that build, use and discard their own tables.
 *
 * Thread safety (BREAKING CHANGE at 3.0.0 — the 2.x promise is withdrawn):
 *   A `beta_param_cache_t*` is THREAD-CONFINED. It is mutated by every compute
 *   call; concurrent use of one cache from more than one thread is undefined.
 *   Give every thread its own cache. `beta_param_tables_t*` and
 *   `beta_param_node_set_t*` are immutable after creation and may be shared
 *   across threads for concurrent reads, including as the backing tables of
 *   per-thread caches created with beta_param_cache_create_shared().
 *
 * Lifetime:
 *   A tables handle passed to beta_param_cache_create_shared() or
 *   beta_param_node_set_create() MUST outlive every cache and node set built
 *   from it — the cache/node set holds a reference, not a copy. Destroy order:
 *   node sets and caches first, their tables last.
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
 *   zero-filled and the cache is returned to a cold state — the next call
 *   recomputes from scratch. No partially updated results are ever visible.
 *   A NULL handle passed into a compute function returns
 *   BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED (2).
 *
 * Precondition — finite input:
 *   `params` and `thetas` must be finite. Non-finite input is undefined
 *   behavior: the library cannot detect NaN under fast-math, so validation
 *   comparisons silently pass and the call returns BETA_PARAM_VALID (0) with
 *   NaN outputs. Screen inputs before calling.
 */

#ifndef BETA_PARAMETERIZATION_H
#define BETA_PARAMETERIZATION_H

#ifdef __cplusplus
extern "C" {
#endif

/* --- Limits --- */
/** Highest Legendre order a tables handle may be built for. */
#define BETA_PARAM_MAX_PARAMS_LIMIT 64
/** Highest n_params a cache (cached tier) accepts. */
#define BETA_PARAM_CACHE_MAX_PARAMS 8

/* --- Shared contract status codes (0-99, identical numbers in every
 *     shape-parameterization library) --- */
#define BETA_PARAM_VALID                        0
#define BETA_PARAM_ERROR_TOO_MANY_PARAMS        1
#define BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED  2
#define BETA_PARAM_ERROR_INVALID_GRID           3
#define BETA_PARAM_ERROR_WRONG_PARAM_COUNT      4
#define BETA_PARAM_ERROR_INVALID_INIT           5
#define BETA_PARAM_ERROR_TABLES_NOT_INITIALIZED 6

/* --- Library status codes (>= 100, append-only after 3.0.0) --- */
#define BETA_PARAM_ERROR_NORTH_POLE          100
#define BETA_PARAM_ERROR_SOUTH_POLE          101
#define BETA_PARAM_ERROR_INTERIOR_NEGATIVE   102
#define BETA_PARAM_ERROR_COM_NOT_CONVERGED   103
#define BETA_PARAM_ERROR_INVALID_BUFFER_SIZE 104
#define BETA_PARAM_ERROR_POLE_NODE           105
#define BETA_PARAM_ERROR_NODE_SET_MISMATCH   106

/* --- Opaque handles --- */
typedef struct beta_param_tables beta_param_tables_t;
typedef struct beta_param_cache beta_param_cache_t;
typedef struct beta_param_node_set beta_param_node_set_t;

/* --- Diagnostics --- */

/**
 * Fixed description of a status code. Never NULL; unknown codes map to an
 * "unknown status code" string. The pointer is to static storage: do not free
 * it, and it stays valid for the life of the process.
 */
const char* beta_param_status_message(int status);

/* --- Tables lifecycle --- */

/**
 * Build the shared immutable level. Returns NULL on failure.
 *
 * @param max_l     1 .. BETA_PARAM_MAX_PARAMS_LIMIT
 * @param thetas    Primary theta set in radians; at least 2 entries, none at
 *                  or beyond a pole (cos(theta)^2 == 1 in double precision)
 * @param n_thetas  Number of entries in thetas
 * @param status    Nullable; receives BETA_PARAM_VALID or the rejecting code
 */
beta_param_tables_t* beta_param_tables_create(
        int max_l, const double* thetas, int n_thetas, int* status);

/** Destroy a tables handle. NULL-safe. Every cache and node set built from it
 *  must already be destroyed. */
void beta_param_tables_destroy(beta_param_tables_t* tables);

/* --- Cache lifecycle --- */

/**
 * Create a cache that owns private tables built with max_l = n_params.
 * Returns NULL on failure.
 *
 * @param n_params         1 .. BETA_PARAM_CACHE_MAX_PARAMS; the exact length
 *                         every later params array must have
 * @param thetas           Primary theta set in radians (see tables_create)
 * @param n_thetas         Number of entries in thetas
 * @param conserve_volume  Nonzero: rescale radii to fixed volume
 * @param apply_com        Nonzero: apply the centre-of-mass correction
 * @param status           Nullable; receives the rejecting code on failure
 */
beta_param_cache_t* beta_param_cache_create(
        int n_params, const double* thetas, int n_thetas,
        int conserve_volume, int apply_com, int* status);

/**
 * Create a cache over caller-owned shared tables. Returns NULL on failure.
 * The tables handle MUST outlive the cache; it is referenced, not copied, and
 * beta_param_cache_destroy() never frees it.
 *
 * @param tables           Tables handle; n_params must not exceed its max_l
 * @param n_params         1 .. BETA_PARAM_CACHE_MAX_PARAMS
 * @param conserve_volume  Nonzero: rescale radii to fixed volume
 * @param apply_com        Nonzero: apply the centre-of-mass correction
 * @param status           Nullable; receives the rejecting code on failure
 */
beta_param_cache_t* beta_param_cache_create_shared(
        const beta_param_tables_t* tables, int n_params,
        int conserve_volume, int apply_com, int* status);

/** Destroy a cache. NULL-safe. Shared tables are left untouched. */
void beta_param_cache_destroy(beta_param_cache_t* cache);

/* --- Node-set lifecycle --- */

/**
 * Build an extra evaluation set (thetas plus Legendre P_k and P_k' tables)
 * sized to `tables`. Returns NULL on failure — including a pole node, which is
 * rejected with BETA_PARAM_ERROR_POLE_NODE; use the resolve_shape polar radii
 * for the poles instead.
 *
 * @param tables    Tables handle; must outlive the node set
 * @param thetas    Node angles in radians; any order, need not be uniform
 * @param n_thetas  Number of nodes (at least 2)
 * @param status    Nullable; receives the rejecting code on failure
 */
beta_param_node_set_t* beta_param_node_set_create(
        const beta_param_tables_t* tables, const double* thetas, int n_thetas,
        int* status);

/** Destroy a node set. NULL-safe. */
void beta_param_node_set_destroy(beta_param_node_set_t* node_set);

/* --- Cached computes ---
 *
 * Every one of these mutates the cache (thread-confined, see above). `params`
 * must hold exactly the cache's n_params entries, else
 * BETA_PARAM_ERROR_WRONG_PARAM_COUNT. Output buffer lengths must equal the
 * cache's theta count (or the node set's node count), else
 * BETA_PARAM_ERROR_INVALID_BUFFER_SIZE. Radii and derivatives are
 * COM-corrected and volume-scaled per the cache's flags.
 */

/** R(theta) at the cache's primary thetas. `n_radii` must equal that count. */
int beta_param_cache_radius_grid(
        beta_param_cache_t* cache, const double* params, int n_params,
        double* radii, int n_radii);

/** R(theta) and dR/dtheta at the cache's primary thetas. */
int beta_param_cache_radius_and_derivative(
        beta_param_cache_t* cache, const double* params, int n_params,
        double* radii, double* dr_dthetas, int n_radii);

/**
 * R(theta) at the primary thetas with NO validation gates and NO volume
 * scaling — the rendering/diagnostic path. It reports usage errors only (not
 * initialized, wrong parameter count, buffer size, COM non-convergence), so a
 * shape rejected by the checked path still yields its (partly negative)
 * outline instead of a zero-filled buffer.
 *
 * CAUTION: the radii are UNSCALED even on a cache created with
 * conserve_volume = 1 — this path never computes the volume factor. Mixing it
 * with beta_param_cache_radius_grid() on such a cache draws two outlines of
 * different size for one shape; scale by beta_param_cache_resolve_shape()'s
 * volume_factor if the sizes must agree.
 */
int beta_param_cache_radius_grid_unchecked(
        beta_param_cache_t* cache, const double* params, int n_params,
        double* radii, int n_radii);

/**
 * Resolve a shape without evaluating a grid: the COM-corrected beta10, the
 * analytic polar radii, and the applied volume factor.
 *
 * `corrected_beta10` is a beta-space value and is never volume-scaled;
 * `r_north` and `r_south` are scaled. `volume_factor` is exactly 1.0 when the
 * cache was created with conserve_volume = 0.
 */
int beta_param_cache_resolve_shape(
        beta_param_cache_t* cache, const double* params, int n_params,
        double* corrected_beta10, double* r_north, double* r_south,
        double* volume_factor);

/**
 * R(theta) and dR/dtheta at a node set's thetas. The node set must be built
 * from tables whose max_l is at least the cache's n_params, else
 * BETA_PARAM_ERROR_NODE_SET_MISMATCH. `n_nodes` must equal the node set's node
 * count. This evaluation is not cached; the resolve/validation/volume
 * intermediates behind it are.
 */
int beta_param_cache_node_radius_and_derivative(
        beta_param_cache_t* cache, const beta_param_node_set_t* node_set,
        const double* params, int n_params,
        double* radii, double* dr_dthetas, int n_nodes);

/* --- Standalone computes (tier 1) ---
 *
 * One-off: tables are built, used and discarded per call — no handle, nothing
 * to free, no engine. `n_params` may be 1 .. BETA_PARAM_MAX_PARAMS_LIMIT (the
 * cached-tier cap does not apply); above that,
 * BETA_PARAM_ERROR_TOO_MANY_PARAMS, never silent truncation. Output buffers
 * hold n_thetas doubles. Outputs are zero-filled on failure.
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
