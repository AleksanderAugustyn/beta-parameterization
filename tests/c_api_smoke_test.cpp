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

    // --- Node-set-only cache: n_grid <= 0 sentinel ---
    beta_param_cache_t* lean = beta_param_cache_create(
            8, 0, static_cast<int>(buf.size()), buf.data());
    check(lean != nullptr, "C: create(8, 0) node-set-only cache succeeds");
    s = beta_param_cache_compute_radius_grid(
            lean, params.data(), static_cast<int>(params.size()),
            radii.data(), static_cast<int>(buf.size()), buf.data());
    check(s == BETA_PARAM_ERROR_NO_UNIFORM_GRID,
          "C: uniform entry point on node-set-only cache -> code 10");
    check(buf[0] != '\0', "C: NO_UNIFORM_GRID message is non-empty");
    beta_param_cache_destroy(lean);

    // --- C++ wrapper ---
    bool threw = false;
    try {
        beta_param::Cache bad(0, 181);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    check(threw, "hpp: constructor throws on invalid max_beta_params");

    beta_param::Cache lean_cxx{8};
    check(lean_cxx.n_grid() == 0, "C++: single-arg Cache reports n_grid 0");

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

    // --- Node-set API (C) ---
    {
        std::array<char, 256> nbuf{};
        const double thetas[3] = {0.4, 1.5707963267948966, 2.7};
        beta_param_cache_t* ns_cache =
                beta_param_cache_create(8, 181, static_cast<int>(nbuf.size()), nbuf.data());
        check(ns_cache != nullptr, "C: node-set cache create");

        beta_param_node_set_t* nodes = beta_param_node_set_create(
                ns_cache, thetas, 3, static_cast<int>(nbuf.size()), nbuf.data());
        check(nodes != nullptr, "C: node_set_create succeeds");

        const double ns_params[4] = {0.1, 0.2, 0.05, 0.1};
        double beta_con[8];
        double corrected_beta10 = -1.0, r_north = -1.0, r_south = -1.0;
        int st = beta_param_cache_resolve_shape(
                ns_cache, ns_params, 4, beta_con, &corrected_beta10, &r_north, &r_south,
                1, static_cast<int>(nbuf.size()), nbuf.data());
        check(st == BETA_PARAM_VALID, "C: resolve_shape VALID");
        check(r_north > 0.0 && r_south > 0.0, "C: resolve_shape pole radii positive");

        double corrected_no_com = -1.0;
        st = beta_param_cache_resolve_shape(
                ns_cache, ns_params, 4, beta_con, &corrected_no_com, &r_north, &r_south,
                0, static_cast<int>(nbuf.size()), nbuf.data());
        check(st == BETA_PARAM_VALID, "C: resolve_shape no-COM VALID");
        check(corrected_no_com == 0.1, "C: no-COM corrected_beta10 == input beta10");

        double ns_radii[3], ns_dr[3];
        st = beta_param_cache_compute_radius_and_derivative(
                ns_cache, nodes, beta_con, 8, ns_radii, ns_dr, 3,
                static_cast<int>(nbuf.size()), nbuf.data());
        check(st == BETA_PARAM_VALID, "C: radius_and_derivative VALID");
        for (int i = 0; i < 3; ++i) check(ns_radii[i] > 0.0, "C: node-set radius positive");

        // Pole node must be rejected with the new code (NULL handle, message set)
        const double pole_theta[1] = {0.0};
        beta_param_node_set_t* bad = beta_param_node_set_create(
                ns_cache, pole_theta, 1, static_cast<int>(nbuf.size()), nbuf.data());
        check(bad == nullptr, "C: pole node rejected with NULL handle");
        check(nbuf[0] != '\0', "C: pole rejection message non-empty");

        beta_param_node_set_destroy(nodes);
        beta_param_node_set_destroy(nullptr);  // null-safe by contract
        beta_param_cache_destroy(ns_cache);
    }

    // --- Node-set API (C++ wrapper) ---
    {
        beta_param::Cache c8(8, 181);
        const std::array<double, 3> thetas{0.4, 1.5707963267948966, 2.7};
        beta_param::NodeSet nodes(c8, thetas);
        check(nodes.n_nodes() == 3, "hpp: NodeSet n_nodes");

        const std::vector<double> ns_params{0.1, 0.2, 0.05, 0.1};
        std::vector<double> beta_con(8);
        double corrected_beta10 = 0.0, r_north = 0.0, r_south = 0.0;
        std::string msg;
        auto st = c8.resolve_shape(ns_params, beta_con, corrected_beta10, r_north, r_south, msg);
        check(st == beta_param::Status::Valid, "hpp: resolve_shape Valid");

        std::vector<double> r(3), dr(3);
        st = c8.compute_radius_and_derivative(beta_con, nodes, r, dr, msg);
        check(st == beta_param::Status::Valid, "hpp: radius_and_derivative Valid");
        check(r[0] > 0.0 && r[1] > 0.0 && r[2] > 0.0, "hpp: node radii positive");

        bool pole_threw = false;
        try {
            const std::array<double, 1> pole{0.0};
            beta_param::NodeSet bad(c8, pole);
        } catch (const std::runtime_error&) {
            pole_threw = true;
        }
        check(pole_threw, "hpp: NodeSet constructor throws on pole node");
    }

    std::printf("c_api_smoke_test: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
