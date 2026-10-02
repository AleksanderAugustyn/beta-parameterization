/**
 * @file beta_parameterization.hpp
 * @brief Header-only C++20 RAII wrapper around the 4.0.0 C API.
 *
 * Requires C++20 (std::span). C++17 callers use the C API directly.
 *
 * Two tiers:
 *   - One-shot: the free functions `radius_grid_standalone` and
 *     `radius_and_derivative_standalone`.
 *   - Read-only cache: `Cache`, built once, plus optional `NodeSet`s for
 *     extra evaluation thetas.
 *
 * Thread model:
 *   `Cache` and `NodeSet` are immutable after construction. Every compute
 *   method is `const` and may be called concurrently from any number of
 *   threads on the same object. Construction, move and destruction must not
 *   race with any other use of the same object.
 *
 * Lifetime:
 *   A `Cache` must outlive every `NodeSet` built from it. Destroy node sets
 *   first. Moving a `Cache` keeps the underlying handle, so a move does not
 *   affect its node sets.
 *
 * Diagnostics:
 *   Constructors throw `std::runtime_error` (the only failure mode with no
 *   return value). Everything else returns a `Status`; `Status::valid` is
 *   success. Buffer-size and parameter-count mismatches are Status codes, not
 *   exceptions: the wrapper forwards the spans' sizes as the C length
 *   arguments and lets the library reject them.
 *
 * Failure behavior (from the C layer): on any nonzero status, output buffers
 * are zero-filled. No state exists, so the next call is unaffected.
 *
 * Precondition — finite input: `params` (and every theta) must be finite.
 * Non-finite input is undefined behavior; the library cannot detect NaN under
 * fast-math, so no check rejects it, and no particular result is promised. A
 * NaN parameter may yield Status::valid with NaN outputs; a NaN trailing
 * parameter is trimmed like a zero, giving the finite outputs of the shorter
 * vector. Screen inputs before calling.
 */

#ifndef BETA_PARAMETERIZATION_HPP
#define BETA_PARAMETERIZATION_HPP

#include "beta_parameterization.h"

#include <algorithm>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace beta_param {

/** Longest parameter vector either tier accepts; the highest max_params. */
inline constexpr int max_params_limit = BETA_PARAM_MAX_PARAMS_LIMIT;

enum class Status : int {
    valid                  = BETA_PARAM_VALID,
    too_many_params        = BETA_PARAM_ERROR_TOO_MANY_PARAMS,
    cache_not_initialized  = BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED,
    invalid_grid           = BETA_PARAM_ERROR_INVALID_GRID,
    wrong_param_count      = BETA_PARAM_ERROR_WRONG_PARAM_COUNT,
    invalid_init           = BETA_PARAM_ERROR_INVALID_INIT,
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

constexpr int as_flag(const bool on) noexcept { return on ? 1 : 0; }

}  // namespace detail

class Cache;

/**
 * An extra evaluation set (thetas plus their Legendre tables) built from one
 * Cache. Immutable after construction, shareable across threads. Move-only.
 * It serves any Cache whose max_params() does not exceed that of the Cache it
 * was built from.
 */
class NodeSet {
public:
    /**
     * @param cache  Supplies max_params; must outlive this node set
     * @param thetas Node angles in radians, any order, at least 2, no poles
     * @throws std::runtime_error if the library rejects the arguments
     *         (a pole node gives Status::pole_node)
     */
    NodeSet(const Cache& cache, std::span<const double> thetas);

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
 * The read-only cache: normalization constants, Gauss-Legendre tables and
 * Legendre tables at the primary thetas. Build once, share across threads.
 * Every compute method is `const`. Move-only.
 */
class Cache {
public:
    /**
     * @param max_params Longest parameter vector the cache accepts,
     *                   1 .. max_params_limit
     * @param thetas     Primary theta set in radians; at least 2 entries, none
     *                   at a pole
     * @throws std::runtime_error if the library rejects the arguments
     */
    Cache(const int max_params, const std::span<const double> thetas) {
        int status = BETA_PARAM_VALID;
        handle_ = beta_param_cache_create(max_params, thetas.data(),
                                          detail::as_int(thetas.size()), &status);
        if (handle_ == nullptr) {  // create failed: nothing was allocated
            detail::throw_status("beta_param::Cache", status);
        }
        max_params_ = max_params;
        n_thetas_   = detail::as_int(thetas.size());
    }

    Cache(const Cache&)            = delete;
    Cache& operator=(const Cache&) = delete;

    Cache(Cache&& other) noexcept
        : handle_{std::exchange(other.handle_, nullptr)},
          max_params_{std::exchange(other.max_params_, 0)},
          n_thetas_{std::exchange(other.n_thetas_, 0)} {}

    Cache& operator=(Cache&& other) noexcept {
        if (this != &other) {
            beta_param_cache_destroy(handle_);
            handle_     = std::exchange(other.handle_, nullptr);
            max_params_ = std::exchange(other.max_params_, 0);
            n_thetas_   = std::exchange(other.n_thetas_, 0);
        }
        return *this;
    }

    ~Cache() { beta_param_cache_destroy(handle_); }

    [[nodiscard]] int max_params() const noexcept { return max_params_; }
    [[nodiscard]] int n_thetas()   const noexcept { return n_thetas_; }

    /** Raw handle for C-API interop; NULL only in a moved-from object. */
    [[nodiscard]] const beta_param_cache_t* native_handle() const noexcept {
        return handle_;
    }

    /** R(theta) at the primary thetas. `radii` must hold n_thetas() doubles. */
    [[nodiscard]] Status radius_grid(const std::span<const double> params,
                                     const bool conserve_volume, const bool apply_com,
                                     const std::span<double> radii) const {
        return static_cast<Status>(beta_param_cache_radius_grid(
                handle_, params.data(), detail::as_int(params.size()),
                detail::as_flag(conserve_volume), detail::as_flag(apply_com),
                radii.data(), detail::as_int(radii.size())));
    }

    /** R(theta) and dR/dtheta at the primary thetas. Both buffers must hold
     *  exactly n_thetas() doubles; unequal span sizes are rejected here with
     *  Status::invalid_buffer_size rather than silently truncated, and both
     *  output spans are zero-filled as on any other rejection. */
    [[nodiscard]] Status radius_and_derivative(const std::span<const double> params,
                                               const bool conserve_volume,
                                               const bool apply_com,
                                               const std::span<double> radii,
                                               const std::span<double> dr_dthetas) const {
        if (radii.size() != dr_dthetas.size()) {
            std::ranges::fill(radii, 0.0);
            std::ranges::fill(dr_dthetas, 0.0);
            return Status::invalid_buffer_size;
        }
        return static_cast<Status>(beta_param_cache_radius_and_derivative(
                handle_, params.data(), detail::as_int(params.size()),
                detail::as_flag(conserve_volume), detail::as_flag(apply_com),
                radii.data(), dr_dthetas.data(), detail::as_int(radii.size())));
    }

    /** R(theta) with no validation gates and no volume scaling — the
     *  rendering/diagnostic path. Reports usage errors and a failed COM
     *  correction only.
     *  CAUTION: these radii are UNSCALED, so mixing this call with
     *  radius_grid(..., conserve_volume = true, ...) draws two outlines of
     *  different size for one shape. Scale by resolve_shape()'s volume_factor
     *  to match. */
    [[nodiscard]] Status radius_grid_unchecked(const std::span<const double> params,
                                               const bool apply_com,
                                               const std::span<double> radii) const {
        return static_cast<Status>(beta_param_cache_radius_grid_unchecked(
                handle_, params.data(), detail::as_int(params.size()),
                detail::as_flag(apply_com),
                radii.data(), detail::as_int(radii.size())));
    }

    /** Outcome of resolve_shape(). `corrected_beta10` is a beta-space value
     *  and is never volume-scaled; the polar radii are. `volume_factor` is
     *  exactly 1.0 when conserve_volume = false. */
    struct Resolved {
        double corrected_beta10;
        double r_north;
        double r_south;
        double volume_factor;
        Status status;
    };

    /** Resolve a shape without evaluating a grid. */
    [[nodiscard]] Resolved resolve_shape(const std::span<const double> params,
                                         const bool conserve_volume,
                                         const bool apply_com) const {
        Resolved out{};
        const int s = beta_param_cache_resolve_shape(
                handle_, params.data(), detail::as_int(params.size()),
                detail::as_flag(conserve_volume), detail::as_flag(apply_com),
                &out.corrected_beta10, &out.r_north, &out.r_south,
                &out.volume_factor);
        out.status = static_cast<Status>(s);
        return out;
    }

    /** R(theta) and dR/dtheta at a node set's thetas. The node set must come
     *  from a Cache with at least this cache's max_params(), else
     *  Status::node_set_mismatch. Both buffers must hold exactly
     *  node_set.n_nodes() doubles; unequal span sizes are rejected here with
     *  Status::invalid_buffer_size, both output spans zero-filled as on any
     *  other rejection. */
    [[nodiscard]] Status node_radius_and_derivative(
            const std::span<const double> params, const NodeSet& node_set,
            const bool conserve_volume, const bool apply_com,
            const std::span<double> radii, const std::span<double> dr_dthetas) const {
        if (radii.size() != dr_dthetas.size()) {
            std::ranges::fill(radii, 0.0);
            std::ranges::fill(dr_dthetas, 0.0);
            return Status::invalid_buffer_size;
        }
        return static_cast<Status>(beta_param_cache_node_radius_and_derivative(
                handle_, params.data(), detail::as_int(params.size()),
                node_set.native_handle(),
                detail::as_flag(conserve_volume), detail::as_flag(apply_com),
                radii.data(), dr_dthetas.data(), detail::as_int(radii.size())));
    }

private:
    beta_param_cache_t* handle_ = nullptr;
    int                 max_params_ = 0;
    int                 n_thetas_ = 0;
};

inline NodeSet::NodeSet(const Cache& cache, const std::span<const double> thetas) {
    int status = BETA_PARAM_VALID;
    handle_ = beta_param_node_set_create(cache.native_handle(), thetas.data(),
                                         detail::as_int(thetas.size()), &status);
    if (handle_ == nullptr) {
        detail::throw_status("beta_param::NodeSet", status);
    }
    n_nodes_ = detail::as_int(thetas.size());
}

/* --- One-shot computes: build, use and discard their own cache. --- */

/** One-off R(theta); `radii` holds exactly thetas.size() doubles. params may
 *  hold 1 .. max_params_limit entries. A size mismatch is rejected here with
 *  Status::invalid_buffer_size and `radii` zero-filled. */
[[nodiscard]] inline Status radius_grid_standalone(
        const std::span<const double> params, const std::span<const double> thetas,
        const bool conserve_volume, const bool apply_com,
        const std::span<double> radii) {
    if (radii.size() != thetas.size()) {
        std::ranges::fill(radii, 0.0);
        return Status::invalid_buffer_size;
    }
    return static_cast<Status>(beta_param_radius_grid_standalone(
            params.data(), detail::as_int(params.size()),
            thetas.data(), detail::as_int(thetas.size()),
            detail::as_flag(conserve_volume), detail::as_flag(apply_com),
            radii.data()));
}

/** One-off R(theta) and dR/dtheta; both buffers hold exactly thetas.size()
 *  doubles. A size mismatch is rejected here with Status::invalid_buffer_size
 *  and both output spans zero-filled. */
[[nodiscard]] inline Status radius_and_derivative_standalone(
        const std::span<const double> params, const std::span<const double> thetas,
        const bool conserve_volume, const bool apply_com,
        const std::span<double> radii, const std::span<double> dr_dthetas) {
    if (radii.size() != thetas.size() || dr_dthetas.size() != thetas.size()) {
        std::ranges::fill(radii, 0.0);
        std::ranges::fill(dr_dthetas, 0.0);
        return Status::invalid_buffer_size;
    }
    return static_cast<Status>(beta_param_radius_and_derivative_standalone(
            params.data(), detail::as_int(params.size()),
            thetas.data(), detail::as_int(thetas.size()),
            detail::as_flag(conserve_volume), detail::as_flag(apply_com),
            radii.data(), dr_dthetas.data()));
}

}  // namespace beta_param

#endif  // BETA_PARAMETERIZATION_HPP
