// Smoke test for the 4.0.0 C API (raw C surface only — the C++ RAII wrapper
// has its own suite). Exercises the SHARED library: the same binary the Python
// bindings load.
//
// Scope: every prototype in the header is called at least once on its happy
// path, plus the failure modes the C layer itself owns — create failure with
// a status out-parameter, NULL handles, negative counts, a large wrong size
// argument, and the static status-message strings.
#include "beta_parameterization.h"

#include <algorithm>
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
    return std::ranges::all_of(v, [](const double x) { return x > 0.0; });
}

bool all_zero(const std::vector<double>& v) {
    return std::ranges::all_of(
            v, [](const double x) { return std::fpclassify(x) == FP_ZERO; });
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size()
           && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}
}  // namespace

int main() {
    constexpr int n_thetas = 32;
    constexpr int max_params = 8;
    const std::vector<double> thetas = open_theta_grid(n_thetas);
    const std::vector<double> params{0.0, 0.25, 0.10, 0.05};
    const int n_params = static_cast<int>(params.size());

    // --- Cache lifecycle ---
    int status = -1;
    beta_param_cache_t* cache =
            beta_param_cache_create(max_params, thetas.data(), n_thetas, &status);
    check(cache != nullptr, "cache_create succeeds");
    check(status == BETA_PARAM_VALID, "cache_create reports VALID");

    // NULL status pointer is accepted (nullable out-parameter).
    beta_param_cache_t* cache_no_status =
            beta_param_cache_create(max_params, thetas.data(), n_thetas, nullptr);
    check(cache_no_status != nullptr, "cache_create accepts a NULL status pointer");
    beta_param_cache_destroy(cache_no_status);

    // Failed creates: each rejection keeps its own code, the handle is NULL.
    const double pole_thetas[2] = {0.0, 1.0};
    status = -1;
    check(beta_param_cache_create(max_params, pole_thetas, 2, &status) == nullptr
                  && status == BETA_PARAM_ERROR_POLE_NODE,
          "cache_create pole node -> NULL, code 105");
    status = -1;
    check(beta_param_cache_create(BETA_PARAM_MAX_PARAMS_LIMIT + 1, thetas.data(),
                                  n_thetas, &status) == nullptr
                  && status == BETA_PARAM_ERROR_TOO_MANY_PARAMS,
          "cache_create max_params = 65 -> NULL, code 1");
    status = -1;
    check(beta_param_cache_create(0, thetas.data(), n_thetas, &status) == nullptr
                  && status == BETA_PARAM_ERROR_INVALID_INIT,
          "cache_create max_params = 0 -> NULL, code 5");
    status = -1;
    check(beta_param_cache_create(max_params, thetas.data(), 1, &status) == nullptr
                  && status == BETA_PARAM_ERROR_INVALID_GRID,
          "cache_create one theta -> NULL, code 3");
    status = -1;
    check(beta_param_cache_create(max_params, thetas.data(), -5, &status) == nullptr
                  && status == BETA_PARAM_ERROR_INVALID_GRID,
          "cache_create negative n_thetas -> NULL, code 3");

    // The limit itself is accepted.
    status = -1;
    beta_param_cache_t* cache_64 = beta_param_cache_create(
            BETA_PARAM_MAX_PARAMS_LIMIT, thetas.data(), n_thetas, &status);
    check(cache_64 != nullptr && status == BETA_PARAM_VALID,
          "cache_create max_params = 64 succeeds");
    beta_param_cache_destroy(cache_64);

    // --- Node set ---
    const std::vector<double> node_thetas{0.4, 1.5707963267948966, 2.7};
    const int n_nodes = static_cast<int>(node_thetas.size());
    status = -1;
    beta_param_node_set_t* nodes =
            beta_param_node_set_create(cache, node_thetas.data(), n_nodes, &status);
    check(nodes != nullptr, "node_set_create succeeds");
    check(status == BETA_PARAM_VALID, "node_set_create reports VALID");

    status = -1;
    check(beta_param_node_set_create(nullptr, node_thetas.data(), n_nodes, &status)
                  == nullptr
                  && status == BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED,
          "node_set_create NULL cache -> NULL, code 2");

    // --- Cached computes (happy path); the handle is const throughout ---
    const beta_param_cache_t* const_cache = cache;

    std::vector<double> radii(n_thetas, -1.0);
    int s = beta_param_cache_radius_grid(
            const_cache, params.data(), n_params, 1, 1, radii.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "cache_radius_grid VALID");
    check(all_positive(radii), "cache_radius_grid radii positive");

    std::vector<double> radii2(n_thetas, -1.0);
    std::vector<double> dr(n_thetas, 0.0);
    s = beta_param_cache_radius_and_derivative(
            const_cache, params.data(), n_params, 1, 1,
            radii2.data(), dr.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "cache_radius_and_derivative VALID");
    check(same_bits(radii2, radii), "cache_radius_and_derivative radii match radius_grid");

    // Options are per call: the same cache serves the raw shape too.
    std::vector<double> radii_raw(n_thetas, -1.0);
    s = beta_param_cache_radius_grid(
            const_cache, params.data(), n_params, 0, 0, radii_raw.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "cache_radius_grid without options VALID");
    check(!same_bits(radii_raw, radii), "options change the output on one cache");

    std::vector<double> radii_unchecked(n_thetas, -1.0);
    s = beta_param_cache_radius_grid_unchecked(
            const_cache, params.data(), n_params, 0, radii_unchecked.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "cache_radius_grid_unchecked VALID");
    check(same_bits(radii_unchecked, radii_raw),
          "unchecked without COM equals the raw checked radii for a valid shape");

    double corrected_beta10 = -1.0, r_north = -1.0, r_south = -1.0, volume_factor = -1.0;
    s = beta_param_cache_resolve_shape(const_cache, params.data(), n_params, 1, 1,
                                       &corrected_beta10, &r_north, &r_south,
                                       &volume_factor);
    check(s == BETA_PARAM_VALID, "cache_resolve_shape VALID");
    check(r_north > 0.0 && r_south > 0.0, "cache_resolve_shape polar radii positive");
    check(volume_factor > 0.0, "cache_resolve_shape volume factor positive");

    std::vector<double> node_radii(static_cast<std::size_t>(n_nodes), -1.0);
    std::vector<double> node_dr(static_cast<std::size_t>(n_nodes), 0.0);
    s = beta_param_cache_node_radius_and_derivative(
            const_cache, params.data(), n_params, nodes, 1, 1,
            node_radii.data(), node_dr.data(), n_nodes);
    check(s == BETA_PARAM_VALID, "cache_node_radius_and_derivative VALID");
    check(all_positive(node_radii), "node radii positive");

    // --- Short vector: accepted, and identical to its zero-padded form ---
    const std::vector<double> short_params{0.0, 0.25};
    const std::vector<double> padded_params{0.0, 0.25, 0.0, 0.0, 0.0};
    std::vector<double> radii_short(n_thetas, -1.0);
    std::vector<double> radii_padded(n_thetas, -2.0);
    s = beta_param_cache_radius_grid(const_cache, short_params.data(), 2, 1, 1,
                                     radii_short.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "short vector VALID");
    s = beta_param_cache_radius_grid(const_cache, padded_params.data(), 5, 1, 1,
                                     radii_padded.data(), n_thetas);
    check(s == BETA_PARAM_VALID, "zero-padded vector VALID");
    check(same_bits(radii_short, radii_padded), "short == zero-padded, bitwise");

    // --- Buffer-size mismatch is caught by the Fortran layer ---
    std::vector<double> short_radii(n_thetas - 1, 1.0);
    s = beta_param_cache_radius_grid(
            const_cache, params.data(), n_params, 1, 1, short_radii.data(), n_thetas - 1);
    check(s == BETA_PARAM_ERROR_INVALID_BUFFER_SIZE,
          "wrong radii length -> INVALID_BUFFER_SIZE");
    check(all_zero(short_radii), "wrong radii length zero-fills the buffer");

    // --- A LARGE wrong size is a status code, not a crash ---
    // 10M doubles = 80 MB: far above the default 8 MB stack, so with automatic
    // marshalling buffers (3.0.0, Release -fstack-arrays) this call segfaulted
    // before the size ever reached the check. The buffer is honest: it really
    // holds n_big doubles.
    constexpr int n_big = 10'000'000;
    std::vector<double> big(static_cast<std::size_t>(n_big), 1.0);
    s = beta_param_cache_radius_grid(
            const_cache, params.data(), n_params, 1, 1, big.data(), n_big);
    check(s == BETA_PARAM_ERROR_INVALID_BUFFER_SIZE,
          "10M-element wrong buffer -> INVALID_BUFFER_SIZE, no crash");
    check(all_zero(big), "10M-element wrong buffer is zero-filled");

    std::vector<double> big2(static_cast<std::size_t>(n_big), 1.0);
    std::ranges::fill(big, 1.0);
    s = beta_param_cache_radius_and_derivative(
            const_cache, params.data(), n_params, 1, 1, big.data(), big2.data(), n_big);
    check(s == BETA_PARAM_ERROR_INVALID_BUFFER_SIZE,
          "10M-element wrong derivative buffers -> INVALID_BUFFER_SIZE, no crash");
    check(all_zero(big) && all_zero(big2), "10M-element derivative buffers zero-filled");

    // --- Parameter count outside 1 .. max_params ---
    const std::vector<double> nine(9, 0.01);
    s = beta_param_cache_radius_grid(const_cache, nine.data(), 9, 1, 1,
                                     radii.data(), n_thetas);
    check(s == BETA_PARAM_ERROR_WRONG_PARAM_COUNT, "nine params on an 8-cache -> code 4");
    s = beta_param_cache_radius_grid(const_cache, params.data(), 0, 1, 1,
                                     radii.data(), n_thetas);
    check(s == BETA_PARAM_ERROR_WRONG_PARAM_COUNT, "zero params -> code 4");
    s = beta_param_cache_radius_grid(const_cache, params.data(), -3, 1, 1,
                                     radii.data(), n_thetas);
    check(s == BETA_PARAM_ERROR_WRONG_PARAM_COUNT, "negative n_params -> code 4");

    // A NULL params pointer with a zero count is an empty vector, not a crash.
    s = beta_param_cache_radius_grid(const_cache, nullptr, 0, 1, 1,
                                     radii.data(), n_thetas);
    check(s == BETA_PARAM_ERROR_WRONG_PARAM_COUNT, "NULL params, zero count -> code 4");
    s = beta_param_radius_grid_standalone(nullptr, 0, thetas.data(), n_thetas,
                                          0, 0, radii.data());
    check(s == BETA_PARAM_ERROR_WRONG_PARAM_COUNT,
          "one-shot NULL params, zero count -> code 4");

    // --- NULL handles into computes ---
    s = beta_param_cache_radius_grid(nullptr, params.data(), n_params, 1, 1,
                                     radii.data(), n_thetas);
    check(s == BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED,
          "NULL cache into radius_grid -> CACHE_NOT_INITIALIZED");
    s = beta_param_cache_node_radius_and_derivative(
            const_cache, params.data(), n_params, nullptr, 1, 1,
            node_radii.data(), node_dr.data(), n_nodes);
    check(s == BETA_PARAM_ERROR_CACHE_NOT_INITIALIZED,
          "NULL node set into node compute -> CACHE_NOT_INITIALIZED");

    // --- A node set built for a smaller cache is rejected ---
    status = -1;
    beta_param_cache_t* small_cache =
            beta_param_cache_create(2, thetas.data(), n_thetas, &status);
    check(small_cache != nullptr, "small cache_create succeeds");
    beta_param_node_set_t* small_nodes = beta_param_node_set_create(
            small_cache, node_thetas.data(), n_nodes, &status);
    check(small_nodes != nullptr, "small node_set_create succeeds");
    s = beta_param_cache_node_radius_and_derivative(
            const_cache, params.data(), n_params, small_nodes, 1, 1,
            node_radii.data(), node_dr.data(), n_nodes);
    check(s == BETA_PARAM_ERROR_NODE_SET_MISMATCH, "small node set -> NODE_SET_MISMATCH");
    // Destroy order: node sets first, their cache last.
    beta_param_node_set_destroy(small_nodes);
    beta_param_cache_destroy(small_cache);

    // --- One-shot entry points ---
    std::vector<double> sa_radii(n_thetas, -1.0);
    s = beta_param_radius_grid_standalone(params.data(), n_params, thetas.data(),
                                          n_thetas, 1, 1, sa_radii.data());
    check(s == BETA_PARAM_VALID, "radius_grid_standalone VALID");
    check(same_bits(sa_radii, radii2),
          "one-shot matches the cached path bit for bit (cache max_params = 8, n = 4)");

    std::vector<double> sa_radii2(n_thetas, -1.0);
    std::vector<double> sa_dr(n_thetas, 0.0);
    s = beta_param_radius_and_derivative_standalone(
            params.data(), n_params, thetas.data(), n_thetas, 1, 1,
            sa_radii2.data(), sa_dr.data());
    check(s == BETA_PARAM_VALID, "radius_and_derivative_standalone VALID");
    check(same_bits(sa_radii2, radii2) && same_bits(sa_dr, dr),
          "one-shot derivative path matches the cached path bit for bit");

    s = beta_param_radius_grid_standalone(params.data(), 0, thetas.data(), n_thetas,
                                          0, 0, sa_radii.data());
    check(s == BETA_PARAM_ERROR_WRONG_PARAM_COUNT, "one-shot with zero params -> code 4");
    s = beta_param_radius_grid_standalone(params.data(), -1, thetas.data(), n_thetas,
                                          0, 0, sa_radii.data());
    check(s == BETA_PARAM_ERROR_WRONG_PARAM_COUNT, "one-shot with negative n_params -> code 4");
    const std::vector<double> too_many(BETA_PARAM_MAX_PARAMS_LIMIT + 1, 0.001);
    s = beta_param_radius_grid_standalone(too_many.data(),
                                          BETA_PARAM_MAX_PARAMS_LIMIT + 1,
                                          thetas.data(), n_thetas, 0, 0, sa_radii.data());
    check(s == BETA_PARAM_ERROR_TOO_MANY_PARAMS, "one-shot with 65 params -> code 1");
    check(all_zero(sa_radii), "one-shot rejection zero-fills the buffer");

    // --- Status messages are static, non-empty, and code-specific ---
    const char* msg_north = beta_param_status_message(BETA_PARAM_ERROR_NORTH_POLE);
    check(msg_north != nullptr && std::string(msg_north).find("north") != std::string::npos,
          "status_message(100) mentions the north pole");
    check(msg_north == beta_param_status_message(BETA_PARAM_ERROR_NORTH_POLE),
          "status_message returns the same static pointer every call");
    check(std::strcmp(beta_param_status_message(BETA_PARAM_VALID), "valid") == 0,
          "status_message(0) is 'valid'");
    check(std::string(beta_param_status_message(6)).find("unknown") != std::string::npos,
          "status_message on the retired code 6 falls back");
    check(std::string(beta_param_status_message(-12345)).find("unknown") != std::string::npos,
          "status_message on an unknown code falls back");

    // --- Teardown (destroys are NULL-safe by contract) ---
    beta_param_node_set_destroy(nodes);
    beta_param_node_set_destroy(nullptr);
    beta_param_cache_destroy(cache);
    beta_param_cache_destroy(nullptr);

    std::printf("c_api_smoke_test: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
