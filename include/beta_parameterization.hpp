/**
 * @file beta_parameterization.hpp
 * @brief Header-only C++20 RAII wrapper around the 3.0.0 C API.
 *
 * Requires C++20 (std::span). C++17 callers use the C API directly.
 *
 * Thread model (BREAKING CHANGE at 3.0.0 — the 2.x promise is withdrawn):
 *   The 2.x header promised that one Cache could serve concurrent calls from
 *   many threads. It cannot, and no longer claims to.
 *     - `Tables` and `NodeSet` are immutable after construction. Share them
 *       freely across threads for concurrent reads.
 *     - `Cache` is THREAD-CONFINED. Every compute call mutates it, so its
 *       compute methods are non-const and concurrent use from more than one
 *       thread is undefined. Give every thread its own Cache.
 *   The intended pattern is one shared `Tables` plus one `Cache` per thread,
 *   built with the shared-tables constructor.
 *
 * Lifetime:
 *   A `Tables` passed to `Cache` (shared constructor) or to `NodeSet` MUST
 *   outlive them — the handle is referenced, not copied. Destroy caches and
 *   node sets before their tables. Moving a `Tables` keeps the underlying
 *   handle address, so a move does not invalidate dependents; destroying or
 *   assigning over the last owner does.
 *
 * Diagnostics:
 *   Constructors throw `std::runtime_error` (the only failure mode with no
 *   return value). Everything else returns a `Status`; `Status::valid` is
 *   success. Buffer-size and parameter-count mismatches are Status codes, not
 *   exceptions: the wrapper forwards the spans' sizes as the C length
 *   arguments and lets the library reject them.
 *
 * Failure behavior (from the C layer): on any nonzero status, output buffers
 * are zero-filled and the cache returns to a cold state.
 *
 * Precondition — finite input: `params` (and every theta) must be finite.
 * Non-finite input is undefined behavior; the library cannot detect NaN under
 * fast-math, so a NaN parameter yields Status::valid with NaN outputs instead
 * of an error. Screen inputs before calling.
 */

#ifndef BETA_PARAMETERIZATION_HPP
#define BETA_PARAMETERIZATION_HPP

#include "beta_parameterization.h"

#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace beta_param {

/** Highest Legendre order a Tables may be built for. */
inline constexpr int max_params_limit = BETA_PARAM_MAX_PARAMS_LIMIT;
/** Highest n_params a Cache accepts. */
inline constexpr int cache_max_params = BETA_PARAM_CACHE_MAX_PARAMS;

enum class Status : int {
    valid                  = BETA_PARAM_VALID,
    too_many_params        = BETA_PARAM_ERROR_TOO_MANY_PARAMS,
    cache_not_initialized  = BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED,
    invalid_grid           = BETA_PARAM_ERROR_INVALID_GRID,
    wrong_param_count      = BETA_PARAM_ERROR_WRONG_PARAM_COUNT,
    invalid_init           = BETA_PARAM_ERROR_INVALID_INIT,
    tables_not_initialized = BETA_PARAM_ERROR_TABLES_NOT_INITIALIZED,
    north_pole             = BETA_PARAM_ERROR_NORTH_POLE,
    south_pole             = BETA_PARAM_ERROR_SOUTH_POLE,
    interior_negative      = BETA_PARAM_ERROR_INTERIOR_NEGATIVE,
    com_not_converged      = BETA_PARAM_ERROR_COM_NOT_CONVERGED,
    invalid_buffer_size    = BETA_PARAM_ERROR_INVALID_BUFFER_SIZE,
    pole_node              = BETA_PARAM_ERROR_POLE_NODE,
    node_set_mismatch      = BETA_PARAM_ERROR_NODE_SET_MISMATCH,
};

/**
 * Fixed description of a status code. The view spans a static, null-terminated
 * string owned by the library, so it outlives any caller and is safe to read
 * from any thread.
 */
[[nodiscard]] inline std::string_view status_message(const Status status) noexcept {
    return beta_param_status_message(static_cast<int>(status));
}

namespace detail {

[[noreturn]] inline void throw_status(const char* what, const int status) {
    throw std::runtime_error(std::string{what} + ": "
                             + beta_param_status_message(status));
}

constexpr int as_int(const std::size_t n) noexcept { return static_cast<int>(n); }

}  // namespace detail

/**
 * The shared immutable level: normalization constants, Gauss-Legendre tables
 * and Legendre tables at the primary thetas. Build once, share read-only
 * across threads. Move-only.
 */
class Tables {
public:
    /**
     * @param max_l  1 .. max_params_limit
     * @param thetas Primary theta set in radians; at least 2 entries, none at
     *               a pole
     * @throws std::runtime_error if the library rejects the arguments
     */
    Tables(const int max_l, const std::span<const double> thetas) {
        int status = BETA_PARAM_VALID;
        handle_ = beta_param_tables_create(max_l, thetas.data(),
                                           detail::as_int(thetas.size()), &status);
        if (handle_ == nullptr) {  // create failed: nothing was allocated
            detail::throw_status("beta_param::Tables", status);
        }
        n_thetas_ = detail::as_int(thetas.size());
        max_l_    = max_l;
    }

    Tables(const Tables&)            = delete;
    Tables& operator=(const Tables&) = delete;

    Tables(Tables&& other) noexcept
        : handle_{std::exchange(other.handle_, nullptr)},
          max_l_{std::exchange(other.max_l_, 0)},
          n_thetas_{std::exchange(other.n_thetas_, 0)} {}

    Tables& operator=(Tables&& other) noexcept {
        if (this != &other) {
            beta_param_tables_destroy(handle_);
            handle_   = std::exchange(other.handle_, nullptr);
            max_l_    = std::exchange(other.max_l_, 0);
            n_thetas_ = std::exchange(other.n_thetas_, 0);
        }
        return *this;
    }

    ~Tables() { beta_param_tables_destroy(handle_); }

    [[nodiscard]] int max_l()    const noexcept { return max_l_; }
    [[nodiscard]] int n_thetas() const noexcept { return n_thetas_; }

    /** Raw handle for C-API interop; NULL only in a moved-from object. */
    [[nodiscard]] const beta_param_tables_t* native_handle() const noexcept {
        return handle_;
    }

private:
    beta_param_tables_t* handle_ = nullptr;
    int                  max_l_ = 0;
    int                  n_thetas_ = 0;
};

/**
 * An extra evaluation set (thetas plus their Legendre tables) sized to one
 * Tables. Immutable after construction, shareable across threads. Move-only.
 * The Tables it was built from must outlive it.
 */
class NodeSet {
public:
    /**
     * @param tables Must outlive this node set
     * @param thetas Node angles in radians, any order, at least 2, no poles
     * @throws std::runtime_error if the library rejects the arguments
     *         (a pole node gives Status::pole_node)
     */
    NodeSet(const Tables& tables, const std::span<const double> thetas) {
        int status = BETA_PARAM_VALID;
        handle_ = beta_param_node_set_create(tables.native_handle(), thetas.data(),
                                             detail::as_int(thetas.size()), &status);
        if (handle_ == nullptr) {
            detail::throw_status("beta_param::NodeSet", status);
        }
        n_nodes_ = detail::as_int(thetas.size());
    }

    NodeSet(const NodeSet&)            = delete;
    NodeSet& operator=(const NodeSet&) = delete;

    NodeSet(NodeSet&& other) noexcept
        : handle_{std::exchange(other.handle_, nullptr)},
          n_nodes_{std::exchange(other.n_nodes_, 0)} {}

    NodeSet& operator=(NodeSet&& other) noexcept {
        if (this != &other) {
            beta_param_node_set_destroy(handle_);
            handle_  = std::exchange(other.handle_, nullptr);
            n_nodes_ = std::exchange(other.n_nodes_, 0);
        }
        return *this;
    }

    ~NodeSet() { beta_param_node_set_destroy(handle_); }

    [[nodiscard]] int n_nodes() const noexcept { return n_nodes_; }

    /** Raw handle for C-API interop; NULL only in a moved-from object. */
    [[nodiscard]] const beta_param_node_set_t* native_handle() const noexcept {
        return handle_;
    }

private:
    beta_param_node_set_t* handle_ = nullptr;
    int                    n_nodes_ = 0;
};

/**
 * Per-shape working state: the recompute engine and every per-shape buffer.
 *
 * THREAD-CONFINED. Every compute call mutates the cache, which is why they are
 * non-const. One Cache per thread; never share one.
 */
class Cache {
public:
    /**
     * Private-tables cache: owns tables built with max_l = n_params.
     *
     * @param n_params        1 .. cache_max_params; the exact length every
     *                        later params span must have
     * @param thetas          Primary theta set in radians
     * @param conserve_volume Rescale radii to fixed volume
     * @param apply_com       Apply the centre-of-mass correction
     * @throws std::runtime_error if the library rejects the arguments
     */
    Cache(const int n_params, const std::span<const double> thetas,
          const bool conserve_volume, const bool apply_com) {
        int status = BETA_PARAM_VALID;
        handle_ = beta_param_cache_create(n_params, thetas.data(),
                                          detail::as_int(thetas.size()),
                                          conserve_volume ? 1 : 0, apply_com ? 1 : 0,
                                          &status);
        if (handle_ == nullptr) {
            detail::throw_status("beta_param::Cache", status);
        }
        n_params_ = n_params;
        n_thetas_ = detail::as_int(thetas.size());
    }

    /**
     * Shared-tables cache: the per-thread constructor. `tables` is referenced,
     * not copied, and MUST outlive this cache.
     *
     * @param tables   Backing tables; n_params must not exceed its max_l
     * @throws std::runtime_error if the library rejects the arguments
     */
    Cache(const Tables& tables, const int n_params,
          const bool conserve_volume, const bool apply_com) {
        int status = BETA_PARAM_VALID;
        handle_ = beta_param_cache_create_shared(tables.native_handle(), n_params,
                                                 conserve_volume ? 1 : 0,
                                                 apply_com ? 1 : 0, &status);
        if (handle_ == nullptr) {
            detail::throw_status("beta_param::Cache", status);
        }
        n_params_ = n_params;
        n_thetas_ = tables.n_thetas();
    }

    Cache(const Cache&)            = delete;
    Cache& operator=(const Cache&) = delete;

    Cache(Cache&& other) noexcept
        : handle_{std::exchange(other.handle_, nullptr)},
          n_params_{std::exchange(other.n_params_, 0)},
          n_thetas_{std::exchange(other.n_thetas_, 0)} {}

    Cache& operator=(Cache&& other) noexcept {
        if (this != &other) {
            beta_param_cache_destroy(handle_);
            handle_   = std::exchange(other.handle_, nullptr);
            n_params_ = std::exchange(other.n_params_, 0);
            n_thetas_ = std::exchange(other.n_thetas_, 0);
        }
        return *this;
    }

    ~Cache() { beta_param_cache_destroy(handle_); }

    [[nodiscard]] int n_params() const noexcept { return n_params_; }
    [[nodiscard]] int n_thetas() const noexcept { return n_thetas_; }

    /** Raw handle for C-API interop; NULL only in a moved-from object. */
    [[nodiscard]] beta_param_cache_t* native_handle() noexcept { return handle_; }

    /** R(theta) at the primary thetas. `radii` must hold n_thetas() doubles. */
    [[nodiscard]] Status radius_grid(const std::span<const double> params,
                                     const std::span<double> radii) {
        return static_cast<Status>(beta_param_cache_radius_grid(
                handle_, params.data(), detail::as_int(params.size()),
                radii.data(), detail::as_int(radii.size())));
    }

    /** R(theta) and dR/dtheta at the primary thetas. Both buffers must hold
     *  exactly n_thetas() doubles; unequal span sizes are rejected here with
     *  Status::invalid_buffer_size rather than silently truncated. */
    [[nodiscard]] Status radius_and_derivative(const std::span<const double> params,
                                               const std::span<double> radii,
                                               const std::span<double> dr_dthetas) {
        if (radii.size() != dr_dthetas.size()) return Status::invalid_buffer_size;
        return static_cast<Status>(beta_param_cache_radius_and_derivative(
                handle_, params.data(), detail::as_int(params.size()),
                radii.data(), dr_dthetas.data(), detail::as_int(radii.size())));
    }

    /** R(theta) with no validation gates and no volume scaling — the
     *  rendering/diagnostic path. Reports usage errors only.
     *  CAUTION: on a conserve_volume Cache these radii are still UNSCALED, so
     *  mixing this call with radius_grid() draws two outlines of different size
     *  for one shape. Scale by resolve_shape()'s volume_factor to match. */
    [[nodiscard]] Status radius_grid_unchecked(const std::span<const double> params,
                                               const std::span<double> radii) {
        return static_cast<Status>(beta_param_cache_radius_grid_unchecked(
                handle_, params.data(), detail::as_int(params.size()),
                radii.data(), detail::as_int(radii.size())));
    }

    /** Outcome of resolve_shape(). `corrected_beta10` is a beta-space value
     *  and is never volume-scaled; the polar radii are. `volume_factor` is
     *  exactly 1.0 when the cache was built with conserve_volume = false. */
    struct Resolved {
        double corrected_beta10;
        double r_north;
        double r_south;
        double volume_factor;
        Status status;
    };

    /** Resolve a shape without evaluating a grid. */
    [[nodiscard]] Resolved resolve_shape(const std::span<const double> params) {
        Resolved out{};
        const int s = beta_param_cache_resolve_shape(
                handle_, params.data(), detail::as_int(params.size()),
                &out.corrected_beta10, &out.r_north, &out.r_south,
                &out.volume_factor);
        out.status = static_cast<Status>(s);
        return out;
    }

    /** R(theta) and dR/dtheta at a node set's thetas. The node set's tables
     *  must have max_l >= n_params(), else Status::node_set_mismatch. Both
     *  buffers must hold exactly node_set.n_nodes() doubles; unequal span sizes
     *  are rejected here with Status::invalid_buffer_size. */
    [[nodiscard]] Status node_radius_and_derivative(
            const NodeSet& node_set, const std::span<const double> params,
            const std::span<double> radii, const std::span<double> dr_dthetas) {
        if (radii.size() != dr_dthetas.size()) return Status::invalid_buffer_size;
        return static_cast<Status>(beta_param_cache_node_radius_and_derivative(
                handle_, node_set.native_handle(),
                params.data(), detail::as_int(params.size()),
                radii.data(), dr_dthetas.data(), detail::as_int(radii.size())));
    }

private:
    beta_param_cache_t* handle_ = nullptr;
    int                 n_params_ = 0;
    int                 n_thetas_ = 0;
};

/* --- Standalone computes: build, use and discard their own tables. --- */

/** One-off R(theta); `radii` holds exactly thetas.size() doubles. n_params may
 *  go up to max_params_limit. */
[[nodiscard]] inline Status radius_grid_standalone(
        const std::span<const double> params, const std::span<const double> thetas,
        const bool conserve_volume, const bool apply_com,
        const std::span<double> radii) {
    if (radii.size() != thetas.size()) return Status::invalid_buffer_size;
    return static_cast<Status>(beta_param_radius_grid_standalone(
            params.data(), detail::as_int(params.size()),
            thetas.data(), detail::as_int(thetas.size()),
            conserve_volume ? 1 : 0, apply_com ? 1 : 0, radii.data()));
}

/** One-off R(theta) and dR/dtheta; both buffers hold exactly thetas.size()
 *  doubles. */
[[nodiscard]] inline Status radius_and_derivative_standalone(
        const std::span<const double> params, const std::span<const double> thetas,
        const bool conserve_volume, const bool apply_com,
        const std::span<double> radii, const std::span<double> dr_dthetas) {
    if (radii.size() != thetas.size() || dr_dthetas.size() != thetas.size()) {
        return Status::invalid_buffer_size;
    }
    return static_cast<Status>(beta_param_radius_and_derivative_standalone(
            params.data(), detail::as_int(params.size()),
            thetas.data(), detail::as_int(thetas.size()),
            conserve_volume ? 1 : 0, apply_com ? 1 : 0,
            radii.data(), dr_dthetas.data()));
}

}  // namespace beta_param

#endif  // BETA_PARAMETERIZATION_HPP
