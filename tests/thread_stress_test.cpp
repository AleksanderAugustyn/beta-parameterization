// Contract family 3: concurrency. MANY THREADS, ONE CACHE.
//
// One beta_param::Cache and one beta_param::NodeSet are built on the main
// thread and never touched again except through const compute methods. A
// serial pass first computes a reference for every job; then 8 worker threads
// each recompute EVERY job on the same shared objects, concurrently, and
// require their results to equal the serial reference bit for bit — memcmp on
// the raw doubles, no tolerance.
//
// The jobs cycle through all five outputs, all four (conserve_volume,
// apply_com) combinations and three vector lengths, so threads are always
// running different outputs and options against the same tables.
//
// A data race on the shared cache, or any state leaking between calls, shows
// up as a status mismatch or a byte mismatch. Any failure in any thread makes
// the process exit nonzero.
#include "beta_parameterization.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <numbers>
#include <span>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace {

constexpr int n_thetas = 64;
constexpr int n_nodes = 24;
constexpr int max_params = 8;
constexpr int n_threads = 8;
constexpr int n_jobs = 6000;
constexpr int n_routines = 5;

// Open uniform grid theta_i = i*pi/(n+1), i = 1..n — no pole nodes.
std::vector<double> open_theta_grid(const int n) {
    std::vector<double> thetas(static_cast<std::size_t>(n));
    for (int i = 1; i <= n; ++i) {
        thetas[static_cast<std::size_t>(i - 1)] =
                static_cast<double>(i) * std::numbers::pi / static_cast<double>(n + 1);
    }
    return thetas;
}

// Deterministic per-job shape of 6, 4 or 2 parameters. Amplitudes stay small
// enough that every shape resolves in every option combination.
std::vector<double> shape_for(const int job) {
    const double t = static_cast<double>(job % 23 + 1) * 0.004;
    const double u = static_cast<double>(job % 37) * 0.004;
    std::vector<double> params{0.02 * t, 0.30 - u, 0.10 + t,
                               0.05 - 0.5 * u, 0.02 + 0.5 * t, 0.01};
    params.resize(static_cast<std::size_t>(6 - 2 * (job % 3)));
    return params;
}

// Everything one call returns, flattened, plus its status.
struct Result {
    beta_param::Status status = beta_param::Status::valid;
    std::vector<double> values;
};

bool same_result(const Result& a, const Result& b) {
    return a.status == b.status && a.values.size() == b.values.size()
           && std::memcmp(a.values.data(), b.values.data(),
                          a.values.size() * sizeof(double)) == 0;
}

// One job: routine = job % 5, option combination = job % 4. 3, 4 and 5 are
// pairwise coprime, so 60 consecutive jobs cover every (length, options,
// routine) combination.
Result run_job(const beta_param::Cache& cache, const beta_param::NodeSet& nodes,
               const int job) {
    const std::vector<double> params = shape_for(job);
    const bool conserve_volume = (job % 4 == 1) || (job % 4 == 3);
    const bool apply_com = (job % 4 == 2) || (job % 4 == 3);

    Result out;
    switch (job % n_routines) {
        case 0: {
            out.values.assign(n_thetas, -1.0);
            out.status = cache.radius_grid(params, conserve_volume, apply_com, out.values);
            break;
        }
        case 1: {
            out.values.assign(2 * n_thetas, -1.0);
            const std::span<double> all{out.values};
            out.status = cache.radius_and_derivative(
                    params, conserve_volume, apply_com,
                    all.first(n_thetas), all.last(n_thetas));
            break;
        }
        case 2: {
            const beta_param::Cache::Resolved r =
                    cache.resolve_shape(params, conserve_volume, apply_com);
            out.status = r.status;
            out.values = {r.corrected_beta10, r.r_north, r.r_south, r.volume_factor};
            break;
        }
        case 3: {
            out.values.assign(2 * n_nodes, -1.0);
            const std::span<double> all{out.values};
            out.status = cache.node_radius_and_derivative(
                    params, nodes, conserve_volume, apply_com,
                    all.first(n_nodes), all.last(n_nodes));
            break;
        }
        default: {
            out.values.assign(n_thetas, -1.0);
            out.status = cache.radius_grid_unchecked(params, apply_com, out.values);
            break;
        }
    }
    return out;
}

}  // namespace

int main() {
    const std::vector<double> thetas = open_theta_grid(n_thetas);
    const std::vector<double> node_thetas = open_theta_grid(n_nodes);

    // Shared, immutable, read-only: const from here on.
    const beta_param::Cache cache{max_params, thetas};
    const beta_param::NodeSet nodes{cache, node_thetas};

    // Serial reference, computed before any thread exists.
    std::vector<Result> reference;
    reference.reserve(n_jobs);
    int valid_reference = 0;
    for (int job = 0; job < n_jobs; ++job) {
        reference.push_back(run_job(cache, nodes, job));
        if (reference.back().status == beta_param::Status::valid) ++valid_reference;
    }

    std::vector<int> failures(static_cast<std::size_t>(n_threads), 0);
    // char, not bool: std::vector<bool> packs bits, so concurrent writes to
    // distinct "elements" would be a genuine data race.
    std::vector<char> finished(static_cast<std::size_t>(n_threads), 0);

    const auto worker = [&cache, &nodes, &reference, &failures, &finished](const int id) {
        int& fails = failures[static_cast<std::size_t>(id)];
        try {
            // Each thread starts at a different job so that, at any moment,
            // the threads are in different routines and option combinations.
            for (int k = 0; k < n_jobs; ++k) {
                const int job = (k + id * 7) % n_jobs;
                const Result got = run_job(cache, nodes, job);
                if (!same_result(got, reference[static_cast<std::size_t>(job)])) {
                    ++fails;
                    std::printf("thread %d job %d: differs from the serial reference "
                                "(status %d vs %d)\n", id, job,
                                static_cast<int>(got.status),
                                static_cast<int>(
                                        reference[static_cast<std::size_t>(job)].status));
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
    for (int t = 0; t < n_threads; ++t) {
        total += failures[static_cast<std::size_t>(t)];
        if (!finished[static_cast<std::size_t>(t)]) {
            std::printf("thread %d did not report success\n", t);
            ++total;
        }
    }
    // The comparison must not pass vacuously: a rejected shape zero-fills its
    // buffers, and two zero buffers always match.
    if (valid_reference != n_jobs) {
        std::printf("only %d of %d reference jobs were VALID — comparison would be "
                    "partly vacuous\n", valid_reference, n_jobs);
        ++total;
    }

    // Wrapper surface outside the thread model: the one-shot free functions
    // agree bit for bit with the shared cache (max_params 8, vector length 4).
    {
        const std::vector<double> params{0.0, 0.25, 0.10, 0.05};
        std::vector<double> cached(n_thetas);
        std::vector<double> one_shot(n_thetas);
        const beta_param::Status a = cache.radius_grid(params, true, true, cached);
        const beta_param::Status b = beta_param::radius_grid_standalone(
                params, thetas, true, true, one_shot);
        if (a != beta_param::Status::valid || b != beta_param::Status::valid
            || std::memcmp(cached.data(), one_shot.data(),
                           cached.size() * sizeof(double)) != 0) {
            std::printf("one-shot wrapper disagrees with the cached path\n");
            ++total;
        }
    }

    // A moved-from Cache holds no handle: computes report it, they do not crash.
    {
        beta_param::Cache source{4, thetas};
        const beta_param::Cache target{std::move(source)};
        std::vector<double> radii(n_thetas, 1.0);
        const std::vector<double> params{0.0, 0.25};
        // NOLINTNEXTLINE(bugprone-use-after-move): the moved-from state is the point.
        const beta_param::Status moved = source.radius_grid(params, false, false, radii);
        const beta_param::Status live = target.radius_grid(params, false, false, radii);
        if (moved != beta_param::Status::cache_not_initialized
            || live != beta_param::Status::valid || target.max_params() != 4) {
            std::printf("moved-from Cache: status %d, moved-to status %d\n",
                        static_cast<int>(moved), static_cast<int>(live));
            ++total;
        }
    }

    // Constructor failures throw with the library's message.
    try {
        const beta_param::Cache too_many{beta_param::max_params_limit + 1, thetas};
        std::printf("Cache with max_params = 65 did not throw\n");
        ++total;
    } catch (const std::runtime_error& e) {
        if (std::string_view{e.what()}.find("too many") == std::string_view::npos) {
            std::printf("unexpected constructor message: %s\n", e.what());
            ++total;
        }
    }

    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    std::printf("thread_stress_test: %d thread(s) x %d jobs, %d failure(s), %.2f s\n",
                n_threads, n_jobs, total, seconds);
    return total == 0 ? 0 : 1;
}
