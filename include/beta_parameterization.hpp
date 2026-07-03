/**
 * @file beta_parameterization.hpp
 * @brief Header-only C++20 RAII wrapper around the C API.
 *
 * Requires C++20 (uses std::span). For C++17 callers, use the C API directly.
 *
 * Thread safety mirrors the C API: a Cache is immutable after construction;
 * multiple threads may concurrently call compute_radius_grid on the same Cache.
 */

#ifndef BETA_PARAMETERIZATION_HPP
#define BETA_PARAMETERIZATION_HPP

#include "beta_parameterization.h"

#include <array>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

namespace beta_param {

inline constexpr int max_params_limit = BETA_PARAM_MAX_PARAMS_LIMIT;

enum class Status : int {
    Valid                  = BETA_PARAM_VALID,
    ErrorNorthPole         = BETA_PARAM_ERROR_NORTH_POLE,
    ErrorSouthPole         = BETA_PARAM_ERROR_SOUTH_POLE,
    ErrorEmptyParams       = BETA_PARAM_ERROR_EMPTY_PARAMS,
    ErrorInteriorNegative  = BETA_PARAM_ERROR_INTERIOR_NEGATIVE,
    ErrorInvalidMaxParams  = BETA_PARAM_ERROR_INVALID_MAX_PARAMS,
    ErrorTooManyParams     = BETA_PARAM_ERROR_TOO_MANY_PARAMS,
    ErrorComNotConverged   = BETA_PARAM_ERROR_COM_NOT_CONVERGED,
    ErrorInvalidBufferSize = BETA_PARAM_ERROR_INVALID_BUFFER_SIZE,
    ErrorPoleNode          = BETA_PARAM_ERROR_POLE_NODE,
};

namespace detail {

inline constexpr int message_buffer_size = 256;

}  // namespace detail

class NodeSet;

class Cache {
public:
    /// Throws std::runtime_error if init fails.
    Cache(int max_beta_params, int n_grid)
        : max_beta_params_{max_beta_params}, n_grid_{n_grid}
    {
        std::array<char, detail::message_buffer_size> buf{};
        handle_ = beta_param_cache_create(
                max_beta_params, n_grid,
                static_cast<int>(buf.size()), buf.data());
        if (handle_ == nullptr) {
            throw std::runtime_error(
                    std::string{"beta_param::Cache init failed: "} + buf.data());
        }
    }

    Cache(const Cache&)            = delete;
    Cache& operator=(const Cache&) = delete;

    Cache(Cache&& other) noexcept
        : handle_{std::exchange(other.handle_, nullptr)},
          max_beta_params_{other.max_beta_params_},
          n_grid_{other.n_grid_}
    {}

    Cache& operator=(Cache&& other) noexcept {
        if (this != &other) {
            beta_param_cache_destroy(handle_);
            handle_          = std::exchange(other.handle_, nullptr);
            max_beta_params_ = other.max_beta_params_;
            n_grid_          = other.n_grid_;
        }
        return *this;
    }

    ~Cache() { beta_param_cache_destroy(handle_); }

    [[nodiscard]] int max_beta_params() const noexcept { return max_beta_params_; }
    [[nodiscard]] int n_grid()          const noexcept { return n_grid_; }

    [[nodiscard]] Status compute_radius_grid(
            std::span<const double> params,
            std::span<double>       radii,
            std::string&            message) const
    {
        if (static_cast<int>(radii.size()) != n_grid_) {
            message = "radii buffer size " + std::to_string(radii.size())
                    + " does not match cache n_grid " + std::to_string(n_grid_);
            return Status::ErrorInvalidBufferSize;
        }
        std::array<char, detail::message_buffer_size> buf{};
        const int s = beta_param_cache_compute_radius_grid(
                handle_,
                params.data(), static_cast<int>(params.size()),
                radii.data(),
                static_cast<int>(buf.size()), buf.data());
        message.assign(buf.data());
        return static_cast<Status>(s);
    }

    [[nodiscard]] Status compute_radius_grid_with_com_shift(
            std::span<const double> params,
            std::span<double>       radii,
            double&                 corrected_beta10,
            std::string&            message) const
    {
        if (static_cast<int>(radii.size()) != n_grid_) {
            message = "radii buffer size " + std::to_string(radii.size())
                    + " does not match cache n_grid " + std::to_string(n_grid_);
            corrected_beta10 = 0.0;
            return Status::ErrorInvalidBufferSize;
        }
        std::array<char, detail::message_buffer_size> buf{};
        const int s = beta_param_cache_compute_radius_grid_with_com_shift(
                handle_,
                params.data(), static_cast<int>(params.size()),
                radii.data(), &corrected_beta10,
                static_cast<int>(buf.size()), buf.data());
        message.assign(buf.data());
        return static_cast<Status>(s);
    }

    /// Resolve a shape once (COM iteration + polar pre-check); beta_con must
    /// hold max_beta_params() doubles. Feed the result to
    /// compute_radius_and_derivative for any number of node sets.
    [[nodiscard]] Status resolve_shape(
            std::span<const double> params,
            std::span<double>       beta_con,
            double&                 corrected_beta10,
            double&                 r_north,
            double&                 r_south,
            std::string&            message) const
    {
        if (static_cast<int>(beta_con.size()) != max_beta_params_) {
            message = "beta_con buffer size " + std::to_string(beta_con.size())
                    + " does not match cache max_beta_params " + std::to_string(max_beta_params_);
            corrected_beta10 = 0.0;
            r_north          = 0.0;
            r_south          = 0.0;
            return Status::ErrorInvalidBufferSize;
        }
        std::array<char, detail::message_buffer_size> buf{};
        const int s = beta_param_cache_resolve_shape(
                handle_,
                params.data(), static_cast<int>(params.size()),
                beta_con.data(), &corrected_beta10, &r_north, &r_south,
                static_cast<int>(buf.size()), buf.data());
        message.assign(buf.data());
        return static_cast<Status>(s);
    }

    /// Evaluate R and dR/dtheta at a node set (defined after NodeSet below).
    [[nodiscard]] Status compute_radius_and_derivative(
            std::span<const double> beta_con,
            const NodeSet&          node_set,
            std::span<double>       radii,
            std::span<double>       dr_dtheta,
            std::string&            message) const;

private:
    friend class NodeSet;

    beta_param_cache_t* handle_           = nullptr;
    int                 max_beta_params_  = 0;
    int                 n_grid_           = 0;
};

/// RAII wrapper over a node-set handle: precomputed Legendre tables at fixed
/// thetas. Immutable after construction; safe to share across threads.
class NodeSet {
public:
    /// Throws std::runtime_error if the build fails (e.g. a pole node).
    NodeSet(const Cache& cache, std::span<const double> thetas)
        : n_nodes_{static_cast<int>(thetas.size())}
    {
        std::array<char, detail::message_buffer_size> buf{};
        handle_ = beta_param_node_set_create(
                cache.handle_, thetas.data(), static_cast<int>(thetas.size()),
                static_cast<int>(buf.size()), buf.data());
        if (handle_ == nullptr) {
            throw std::runtime_error(
                    std::string{"beta_param::NodeSet build failed: "} + buf.data());
        }
    }

    NodeSet(const NodeSet&)            = delete;
    NodeSet& operator=(const NodeSet&) = delete;

    NodeSet(NodeSet&& other) noexcept
        : handle_{std::exchange(other.handle_, nullptr)},
          n_nodes_{other.n_nodes_}
    {}

    NodeSet& operator=(NodeSet&& other) noexcept {
        if (this != &other) {
            beta_param_node_set_destroy(handle_);
            handle_  = std::exchange(other.handle_, nullptr);
            n_nodes_ = other.n_nodes_;
        }
        return *this;
    }

    ~NodeSet() { beta_param_node_set_destroy(handle_); }

    [[nodiscard]] int n_nodes() const noexcept { return n_nodes_; }

private:
    friend class Cache;

    beta_param_node_set_t* handle_  = nullptr;
    int                    n_nodes_ = 0;
};

[[nodiscard]] inline Status Cache::compute_radius_and_derivative(
        std::span<const double> beta_con,
        const NodeSet&          node_set,
        std::span<double>       radii,
        std::span<double>       dr_dtheta,
        std::string&            message) const
{
    if (static_cast<int>(radii.size()) != node_set.n_nodes()
            || static_cast<int>(dr_dtheta.size()) != node_set.n_nodes()) {
        message = "radii/dr_dtheta buffers must match node set n_nodes = "
                + std::to_string(node_set.n_nodes());
        return Status::ErrorInvalidBufferSize;
    }
    std::array<char, detail::message_buffer_size> buf{};
    const int s = beta_param_cache_compute_radius_and_derivative(
            handle_, node_set.handle_,
            beta_con.data(), static_cast<int>(beta_con.size()),
            radii.data(), dr_dtheta.data(), static_cast<int>(radii.size()),
            static_cast<int>(buf.size()), buf.data());
    message.assign(buf.data());
    return static_cast<Status>(s);
}

[[nodiscard]] inline Status compute_radius_grid_standalone(
        std::span<const double> params, int n_grid,
        std::span<double>       radii,
        std::string&            message)
{
    if (static_cast<int>(radii.size()) != n_grid) {
        message = "radii buffer size " + std::to_string(radii.size())
                + " does not match n_grid " + std::to_string(n_grid);
        return Status::ErrorInvalidBufferSize;
    }
    std::array<char, detail::message_buffer_size> buf{};
    const int s = beta_param_compute_radius_grid_standalone(
            params.data(), static_cast<int>(params.size()), n_grid,
            radii.data(),
            static_cast<int>(buf.size()), buf.data());
    message.assign(buf.data());
    return static_cast<Status>(s);
}

[[nodiscard]] inline Status compute_radius_grid_standalone_with_com_shift(
        std::span<const double> params, int n_grid,
        std::span<double>       radii,
        double&                 corrected_beta10,
        std::string&            message)
{
    if (static_cast<int>(radii.size()) != n_grid) {
        message = "radii buffer size " + std::to_string(radii.size())
                + " does not match n_grid " + std::to_string(n_grid);
        corrected_beta10 = 0.0;
        return Status::ErrorInvalidBufferSize;
    }
    std::array<char, detail::message_buffer_size> buf{};
    const int s = beta_param_compute_radius_grid_standalone_with_com_shift(
            params.data(), static_cast<int>(params.size()), n_grid,
            radii.data(), &corrected_beta10,
            static_cast<int>(buf.size()), buf.data());
    message.assign(buf.data());
    return static_cast<Status>(s);
}

}  // namespace beta_param

#endif  // BETA_PARAMETERIZATION_HPP
