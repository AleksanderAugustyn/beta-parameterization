// Thread model of 3.0.0, positive test: SHARE TABLES, CONFINE CACHES.
//
// One immutable beta_param::Tables is built on the main thread and read
// concurrently by 8 worker threads. Every thread owns its Cache; no cache
// handle ever crosses a thread boundary. Each thread walks 1000 shapes through
// its incremental cache, and every 100th iteration re-computes the same shape
// on a freshly built COLD cache (same tables, same flags, single call) and
// requires the two results to be bit-identical — memcmp on the raw doubles, no
// tolerance. That is the incremental-vs-cold invariant the Fortran bitwise
// suite proves single-threaded, re-checked here under concurrency.
//
// A data race on the shared tables, or cache state leaking between threads,
// shows up as a status mismatch or a byte mismatch. Any failure in any thread
// makes the process exit nonzero.
#include "beta_parameterization.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <numbers>
#include <span>
#include <thread>
#include <vector>

namespace {

constexpr int n_thetas = 64;
constexpr int max_l = 8;
constexpr int n_params = 6;
constexpr int n_threads = 8;
constexpr int n_iters = 1000;
constexpr int verify_every = 100;

// Open uniform grid theta_i = i*pi/(n+1), i = 1..n — no pole nodes.
std::vector<double> open_theta_grid(const int n) {
    std::vector<double> thetas(static_cast<std::size_t>(n));
    for (int i = 1; i <= n; ++i) {
        thetas[static_cast<std::size_t>(i - 1)] =
                static_cast<double>(i) * std::numbers::pi / static_cast<double>(n + 1);
    }
    return thetas;
}

// Deterministic per-(thread, iteration) shape. Amplitudes stay small enough
// that the shapes resolve, but the exact values are irrelevant: the test
// compares two evaluations of the SAME parameters, whatever their status.
std::vector<double> shape_for(const int thread_id, const int iter) {
    const double t = static_cast<double>(thread_id + 1) * 0.011;
    const double u = static_cast<double>(iter % 37) * 0.004;
    return {0.02 * t, 0.30 - u, 0.10 + t, 0.05 - 0.5 * u, 0.02 + 0.5 * t, 0.01};
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size()
           && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

}  // namespace

int main() {
    const std::vector<double> thetas = open_theta_grid(n_thetas);

    beta_param::Tables tables{max_l, thetas};  // shared, immutable, read-only

    std::vector<int> failures(static_cast<std::size_t>(n_threads), 0);
    // char, not bool: std::vector<bool> packs bits, so concurrent writes to
    // distinct "elements" would be a genuine data race.
    std::vector<char> finished(static_cast<std::size_t>(n_threads), 0);

    // Counted so the bitwise comparison cannot pass vacuously: a rejected
    // shape zero-fills both buffers, and two zero buffers always match.
    std::vector<int> valid_computes(static_cast<std::size_t>(n_threads), 0);

    const auto worker = [&tables, &failures, &finished, &valid_computes](const int id) {
        int& fails = failures[static_cast<std::size_t>(id)];
        try {
            // This thread's own cache — never touched by any other thread.
            beta_param::Cache cache{tables, n_params, true, true};

            std::vector<double> radii(n_thetas);
            std::vector<double> dr(n_thetas);
            std::vector<double> cold_radii(n_thetas);
            std::vector<double> cold_dr(n_thetas);

            for (int it = 0; it < n_iters; ++it) {
                const std::vector<double> params = shape_for(id, it);

                const beta_param::Status s =
                        cache.radius_and_derivative(params, radii, dr);
                if (s == beta_param::Status::valid) {
                    ++valid_computes[static_cast<std::size_t>(id)];
                }

                if (it % verify_every != 0) continue;

                // Cold reference: fresh cache over the same shared tables,
                // one call, no incremental history.
                beta_param::Cache cold{tables, n_params, true, true};
                const beta_param::Status cold_s =
                        cold.radius_and_derivative(params, cold_radii, cold_dr);

                if (s != cold_s) {
                    ++fails;
                    std::printf("thread %d iter %d: status %d vs cold %d (%s)\n", id, it,
                                static_cast<int>(s), static_cast<int>(cold_s),
                                beta_param::status_message(s).data());
                    continue;
                }
                if (!same_bits(radii, cold_radii) || !same_bits(dr, cold_dr)) {
                    ++fails;
                    std::printf("thread %d iter %d: incremental result differs from cold\n",
                                id, it);
                }
            }
            finished[static_cast<std::size_t>(id)] = 1;
        } catch (const std::exception& e) {
            ++fails;
            std::printf("thread %d threw: %s\n", id, e.what());
        }
    };

    const auto t0 = std::chrono::steady_clock::now();
    {
        std::vector<std::jthread> pool;
        pool.reserve(n_threads);
        for (int t = 0; t < n_threads; ++t) pool.emplace_back(worker, t);
    }  // all joined here
    const auto t1 = std::chrono::steady_clock::now();

    int total = 0;
    int valid_total = 0;
    for (int t = 0; t < n_threads; ++t) {
        total += failures[static_cast<std::size_t>(t)];
        valid_total += valid_computes[static_cast<std::size_t>(t)];
        if (!finished[static_cast<std::size_t>(t)]) {
            std::printf("thread %d did not report success\n", t);
            ++total;
        }
    }
    if (valid_total != n_threads * n_iters) {
        std::printf("only %d of %d computes were VALID — comparison would be "
                    "partly vacuous\n", valid_total, n_threads * n_iters);
        ++total;
    }

    // Wrapper surface outside the thread model: the standalone free functions
    // must agree bit for bit with a same-max_l, no-flags cache.
    {
        const std::vector<double> params{0.0, 0.25, 0.10, 0.05};
        beta_param::Cache plain{static_cast<int>(params.size()), thetas, false, false};
        std::vector<double> cached(n_thetas);
        std::vector<double> standalone(n_thetas);
        const beta_param::Status a = plain.radius_grid(params, cached);
        const beta_param::Status b = beta_param::radius_grid_standalone(
                params, thetas, false, false, standalone);
        if (a != beta_param::Status::valid || b != beta_param::Status::valid
            || !same_bits(cached, standalone)) {
            std::printf("standalone wrapper disagrees with the cached path\n");
            ++total;
        }
    }

    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    std::printf("thread_stress_test: %d thread(s) x %d iters, %d failure(s), %.2f s\n",
                n_threads, n_iters, total, seconds);
    return total == 0 ? 0 : 1;
}
