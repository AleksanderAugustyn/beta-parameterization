// Smoke test for the 3.0.0 C API (raw C surface only — the C++ RAII wrapper
// has its own suite). Exercises the SHARED library: the same binary the Python
// bindings load.
//
// Scope: every prototype in the header is called at least once on its happy
// path, plus the three failure modes the C layer itself owns — create failure
// with a status out-parameter, NULL handles into computes, and the static
// status-message strings.
#include "beta_parameterization.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <numbers>
#include <string>
#include <vector>

namespace {
int failures = 0;

void check(const bool ok, const char* label) {
    if (!ok) {
        std::printf("FAIL: %s\n", label);
        ++failures;
    }
}

// Open uniform grid theta_i = i*pi/(n+1), i = 1..n — no pole nodes.
std::vector<double> open_theta_grid(const int n) {
    std::vector<double> thetas(static_cast<std::size_t>(n));
    for (int i = 1; i <= n; ++i) {
        thetas[static_cast<std::size_t>(i - 1)] =
                static_cast<double>(i) * std::numbers::pi / static_cast<double>(n + 1);
    }
    return thetas;
}

bool all_positive(const std::vector<double>& v) {
    for (const double x : v) {
        if (!(x > 0.0)) return false;
    }
    return true;
}
}  // namespace

int main() {
    constexpr int n_thetas = 32;
    constexpr int max_l = 8;
    const std::vector<double> thetas = open_theta_grid(n_thetas);
    const std::vector<double> params{0.0, 0.25, 0.10, 0.05};
    const int n_params = static_cast<int>(params.size());

    // --- Tables ---
    int status = -1;
    beta_param_tables_t* tables =
            beta_param_tables_create(max_l, thetas.data(), n_thetas, &status);
    check(tables != nullptr, "tables_create succeeds");
    check(status == BETA_PARAM_VALID, "tables_create reports VALID");

    // NULL status pointer is accepted (nullable out-parameter).
    beta_param_tables_t* tables_no_status =
            beta_param_tables_create(max_l, thetas.data(), n_thetas, nullptr);
    check(tables_no_status != nullptr, "tables_create accepts a NULL status pointer");
    beta_param_tables_destroy(tables_no_status);

    // Pole node in the theta set is rejected; the handle is NULL.
    const double pole_thetas[2] = {0.0, 1.0};
    status = -1;
    beta_param_tables_t* bad_tables =
            beta_param_tables_create(max_l, pole_thetas, 2, &status);
    check(bad_tables == nullptr, "tables_create rejects a pole node with a NULL handle");
    check(status == BETA_PARAM_ERROR_POLE_NODE, "tables_create pole node -> code 105");

    // --- Caches: shared (tables-backed) and private, both flag combinations ---
    status = -1;
    beta_param_cache_t* shared_cache =
            beta_param_cache_create_shared(tables, n_params, 1, 1, &status);
    check(shared_cache != nullptr, "cache_create_shared succeeds");
    check(status == BETA_PARAM_VALID, "cache_create_shared reports VALID");

    status = -1;
    beta_param_cache_t* plain_cache = beta_param_cache_create(
            n_params, thetas.data(), n_thetas, 0, 0, &status);
    check(plain_cache != nullptr, "cache_create succeeds");
    check(status == BETA_PARAM_VALID, "cache_create reports VALID");

    beta_param_cache_t* cache_no_status = beta_param_cache_create(
            n_params, thetas.data(), n_thetas, 1, 0, nullptr);
    check(cache_no_status != nullptr, "cache_create accepts a NULL status pointer");
    beta_param_cache_destroy(cache_no_status);

    // Failed create: n_params above the cached-tier cap.
    status = -1;
    beta_param_cache_t* too_many = beta_param_cache_create(
            9, thetas.data(), n_thetas, 0, 0, &status);
    check(too_many == nullptr, "cache_create with n_params=9 returns NULL");
    check(status == BETA_PARAM_ERROR_TOO_MANY_PARAMS,
          "cache_create with n_params=9 reports TOO_MANY_PARAMS");

    // --- Node set ---
    const std::vector<double> node_thetas{0.4, 1.5707963267948966, 2.7};
    status = -1;
    beta_param_node_set_t* nodes = beta_param_node_set_create(
            tables, node_thetas.data(), static_cast<int>(node_thetas.size()), &status);
    check(nodes != nullptr, "node_set_create succeeds");
    check(status == BETA_PARAM_VALID, "node_set_create reports VALID");

    // --- Cached computes (happy path) ---
    std::vector<double> radii(n_thetas, -1.0);
    int s = beta_param_cache_radius_grid(
            shared_cache, params.data(), n_params, radii.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "cache_radius_grid VALID");
    check(all_positive(radii), "cache_radius_grid radii positive");

    std::vector<double> radii2(n_thetas, -1.0);
    std::vector<double> dr(n_thetas, 0.0);
    s = beta_param_cache_radius_and_derivative(
            shared_cache, params.data(), n_params, radii2.data(), dr.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "cache_radius_and_derivative VALID");
    check(radii2 == radii, "cache_radius_and_derivative radii match radius_grid");

    // Private cache, both flags off — the standalone comparison partner below
    // (same max_l = n_params, so the normalization tables are bit-identical).
    std::vector<double> radii_plain(n_thetas, -1.0);
    s = beta_param_cache_radius_grid(
            plain_cache, params.data(), n_params, radii_plain.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "cache_radius_grid on a private cache VALID");

    std::vector<double> radii_unchecked(n_thetas, -1.0);
    s = beta_param_cache_radius_grid_unchecked(
            plain_cache, params.data(), n_params, radii_unchecked.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "cache_radius_grid_unchecked VALID");
    check(all_positive(radii_unchecked), "cache_radius_grid_unchecked radii positive");

    double corrected_beta10 = -1.0, r_north = -1.0, r_south = -1.0, volume_factor = -1.0;
    s = beta_param_cache_resolve_shape(shared_cache, params.data(), n_params,
                                       &corrected_beta10, &r_north, &r_south,
                                       &volume_factor);
    check(s == BETA_PARAM_VALID, "cache_resolve_shape VALID");
    check(r_north > 0.0 && r_south > 0.0, "cache_resolve_shape polar radii positive");
    check(volume_factor > 0.0, "cache_resolve_shape volume factor positive");

    const int n_nodes = static_cast<int>(node_thetas.size());
    std::vector<double> node_radii(static_cast<std::size_t>(n_nodes), -1.0);
    std::vector<double> node_dr(static_cast<std::size_t>(n_nodes), 0.0);
    s = beta_param_cache_node_radius_and_derivative(
            shared_cache, nodes, params.data(), n_params,
            node_radii.data(), node_dr.data(), n_nodes);
    check(s == BETA_PARAM_VALID, "cache_node_radius_and_derivative VALID");
    check(all_positive(node_radii), "node radii positive");

    // --- Buffer-size mismatch is caught by the Fortran layer ---
    std::vector<double> short_radii(n_thetas - 1, 0.0);
    s = beta_param_cache_radius_grid(
            shared_cache, params.data(), n_params, short_radii.data(), n_thetas - 1);
    check(s == BETA_PARAM_ERROR_INVALID_BUFFER_SIZE,
          "wrong radii length -> INVALID_BUFFER_SIZE");

    // --- Wrong parameter count is caught by the engine ---
    const std::vector<double> five{0.0, 0.1, 0.1, 0.1, 0.1};
    s = beta_param_cache_radius_grid(shared_cache, five.data(), 5, radii.data(), n_thetas);
    check(s == BETA_PARAM_ERROR_WRONG_PARAM_COUNT,
          "wrong params length -> WRONG_PARAM_COUNT");

    // --- NULL handles into computes ---
    s = beta_param_cache_radius_grid(nullptr, params.data(), n_params, radii.data(), n_thetas);
    check(s == BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED,
          "NULL cache into radius_grid -> CACHE_NOT_INITIALIZED");
    s = beta_param_cache_node_radius_and_derivative(
            shared_cache, nullptr, params.data(), n_params,
            node_radii.data(), node_dr.data(), n_nodes);
    check(s == BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED,
          "NULL node set into node compute -> CACHE_NOT_INITIALIZED");

    // --- Standalone entry points ---
    std::vector<double> sa_radii(n_thetas, -1.0);
    s = beta_param_radius_grid_standalone(params.data(), n_params, thetas.data(),
                                          n_thetas, 0, 0, sa_radii.data());
    check(s == BETA_PARAM_VALID, "radius_grid_standalone VALID");
    check(sa_radii == radii_plain,
          "standalone matches the same-max_l no-flags cache path bit for bit");

    std::vector<double> sa_radii2(n_thetas, -1.0);
    std::vector<double> sa_dr(n_thetas, 0.0);
    s = beta_param_radius_and_derivative_standalone(
            params.data(), n_params, thetas.data(), n_thetas, 1, 1,
            sa_radii2.data(), sa_dr.data());
    check(s == BETA_PARAM_VALID, "radius_and_derivative_standalone (both flags) VALID");
    check(all_positive(sa_radii2), "standalone radii positive");

    // --- Status messages are static, non-empty, and code-specific ---
    const char* msg_north = beta_param_status_message(BETA_PARAM_ERROR_NORTH_POLE);
    check(msg_north != nullptr && std::string(msg_north).find("north") != std::string::npos,
          "status_message(100) mentions the north pole");
    check(msg_north == beta_param_status_message(BETA_PARAM_ERROR_NORTH_POLE),
          "status_message returns the same static pointer every call");
    check(std::strcmp(beta_param_status_message(BETA_PARAM_VALID), "valid") == 0,
          "status_message(0) is 'valid'");
    check(std::string(beta_param_status_message(-12345)).find("unknown") != std::string::npos,
          "status_message on an unknown code falls back");

    // --- Teardown (destroys are NULL-safe by contract) ---
    beta_param_node_set_destroy(nodes);
    beta_param_node_set_destroy(nullptr);
    beta_param_cache_destroy(shared_cache);
    beta_param_cache_destroy(plain_cache);
    beta_param_cache_destroy(nullptr);
    beta_param_tables_destroy(tables);
    beta_param_tables_destroy(nullptr);

    std::printf("c_api_smoke_test: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
