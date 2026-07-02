// Smoke test for the C API (raw) and the C++20 RAII wrapper.
// Exercises the SHARED library — the same binary the Python bindings load.
#include "beta_parameterization.h"
#include "beta_parameterization.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {
int failures = 0;

void check(bool ok, const char* label) {
    if (!ok) {
        std::printf("FAIL: %s\n", label);
        ++failures;
    }
}
}  // namespace

int main() {
    std::array<char, 256> buf{};

    // --- Raw C API ---
    check(beta_param_cache_create(0, 181, static_cast<int>(buf.size()), buf.data()) == nullptr,
          "C: create with max_beta_params=0 returns NULL");
    check(buf[0] != '\0', "C: failure message is non-empty");

    beta_param_cache_t* cache = beta_param_cache_create(8, 181, static_cast<int>(buf.size()), buf.data());
    check(cache != nullptr, "C: create(8, 181) succeeds");

    std::vector<double> params{0.0, 0.25, 0.10, 0.05};
    std::vector<double> radii(181, -1.0);
    int s = beta_param_cache_compute_radius_grid(
            cache, params.data(), static_cast<int>(params.size()),
            radii.data(), static_cast<int>(buf.size()), buf.data());
    check(s == BETA_PARAM_VALID, "C: compute returns VALID");
    // R(0) = 1 + sum(beta_l * N_l), N_l = sqrt((2l+1)/(4pi)), 4pi = 16*atan(1).
    const double r0_expected = 1.0
            + 0.25 * std::sqrt(5.0 / (16.0 * std::atan(1.0)))
            + 0.10 * std::sqrt(7.0 / (16.0 * std::atan(1.0)))
            + 0.05 * std::sqrt(9.0 / (16.0 * std::atan(1.0)));
    check(std::fabs(radii[0] - r0_expected) < 1e-12, "C: R(0) matches analytic pole sum");

    std::vector<double> too_many(9, 0.1);
    s = beta_param_cache_compute_radius_grid(
            cache, too_many.data(), 9, radii.data(), static_cast<int>(buf.size()), buf.data());
    check(s == BETA_PARAM_ERROR_TOO_MANY_PARAMS, "C: 9 params on max-8 cache -> code 6");

    // Message truncation: tiny buffer must still be null-terminated within bounds.
    std::array<char, 8> tiny{};
    tiny.fill('X');
    s = beta_param_cache_compute_radius_grid(
            cache, too_many.data(), 9, radii.data(), static_cast<int>(tiny.size()), tiny.data());
    check(std::memchr(tiny.data(), '\0', tiny.size()) != nullptr,
          "C: truncated message is null-terminated within the buffer");

    // Standalone bit-parity contract: standalone(params) equals a cache built
    // with max_beta_params = n_params. (A max=8 cache differs by <= 1 ulp at
    // some grid points — the FF norm-constant tables are max-dependent in the
    // last bit — so parity is defined against the same-max cache.)
    beta_param_cache_t* cache4 = beta_param_cache_create(4, 181, static_cast<int>(buf.size()), buf.data());
    check(cache4 != nullptr, "C: create(4, 181) succeeds");
    beta_param_cache_compute_radius_grid(cache4, params.data(), static_cast<int>(params.size()),
                                         radii.data(), static_cast<int>(buf.size()), buf.data());
    std::vector<double> radii_sa(181, -1.0);
    s = beta_param_compute_radius_grid_standalone(
            params.data(), static_cast<int>(params.size()), 181,
            radii_sa.data(), static_cast<int>(buf.size()), buf.data());
    check(s == BETA_PARAM_VALID, "C: standalone VALID");
    bool identical = true;
    for (int i = 0; i < 181; ++i) identical = identical && (radii_sa[i] == radii[i]);
    check(identical, "C: standalone bit-identical to same-max cache path");
    beta_param_cache_destroy(cache4);

    // Regression for review finding F1 (seed S-A): failed standalone must
    // write zeros, not garbage.
    std::vector<double> radii_fail(41, -7.0);
    s = beta_param_compute_radius_grid_standalone(
            params.data(), 0, 41, radii_fail.data(), static_cast<int>(buf.size()), buf.data());
    check(s != BETA_PARAM_VALID, "C: standalone with n_params=0 fails");
    bool all_defined = true;
    for (double r : radii_fail) all_defined = all_defined && (r == 0.0);
    check(all_defined, "C: failed standalone writes zeros, not garbage (F1 regression)");

    beta_param_cache_destroy(cache);
    beta_param_cache_destroy(nullptr);   // null-safe by contract

    // --- C++ wrapper ---
    bool threw = false;
    try {
        beta_param::Cache bad(0, 181);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    check(threw, "hpp: constructor throws on invalid max_beta_params");

    // Same max as the cache4 result held in `radii` — parity is per-max (see above).
    beta_param::Cache cxx(4, 181);
    std::vector<double> radii_cxx(181);
    std::string message;
    auto status = cxx.compute_radius_grid(params, radii_cxx, message);
    check(status == beta_param::Status::Valid, "hpp: compute Valid");
    identical = true;
    for (int i = 0; i < 181; ++i) identical = identical && (radii_cxx[i] == radii[i]);
    check(identical, "hpp: bit-identical to raw C path");

    std::vector<double> wrong_size(180);
    status = cxx.compute_radius_grid(params, wrong_size, message);
    check(status == beta_param::Status::ErrorInvalidBufferSize, "hpp: wrong radii size pre-checked");

    double corrected = -1.0;
    std::vector<double> asym{0.30, 0.60, 0.40, 0.10};
    status = cxx.compute_radius_grid_with_com_shift(asym, radii_cxx, corrected, message);
    check(status == beta_param::Status::Valid, "hpp: com-shift Valid");
    check(std::fabs(corrected - 0.30) > 1e-6, "hpp: corrected beta10 moved");

    std::printf("c_api_smoke_test: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
