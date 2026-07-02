// Thread-safety stress: one shared cache, 8 threads x 200 iterations x 16 shapes.
// Every concurrent result must be bit-identical to the serial reference.
#include "beta_parameterization.h"

#include <array>
#include <atomic>
#include <cstdio>
#include <thread>
#include <vector>

int main() {
    constexpr int n_grid = 181;
    constexpr int n_threads = 8;
    constexpr int n_iters = 200;

    std::array<char, 256> buf{};
    beta_param_cache_t* cache = beta_param_cache_create(8, n_grid, static_cast<int>(buf.size()), buf.data());
    if (cache == nullptr) { std::printf("FAIL: cache create\n"); return 1; }

    std::vector<std::vector<double>> shapes;
    for (int k = 0; k < 16; ++k) {
        shapes.push_back({0.0, 0.05 * k, 0.02 * k, 0.01 * k});
    }

    // Serial reference.
    std::vector<std::vector<double>> reference(shapes.size(), std::vector<double>(n_grid));
    std::vector<int> ref_status(shapes.size());
    for (size_t k = 0; k < shapes.size(); ++k) {
        std::array<char, 256> mb{};
        ref_status[k] = beta_param_cache_compute_radius_grid(
                cache, shapes[k].data(), static_cast<int>(shapes[k].size()),
                reference[k].data(), static_cast<int>(mb.size()), mb.data());
    }

    std::atomic<int> failures{0};
    auto worker = [&]() {
        std::vector<double> radii(n_grid);
        std::array<char, 256> mb{};
        for (int it = 0; it < n_iters; ++it) {
            for (size_t k = 0; k < shapes.size(); ++k) {
                const int s = beta_param_cache_compute_radius_grid(
                        cache, shapes[k].data(), static_cast<int>(shapes[k].size()),
                        radii.data(), static_cast<int>(mb.size()), mb.data());
                if (s != ref_status[k]) { ++failures; continue; }
                if (s != BETA_PARAM_VALID) continue;
                for (int i = 0; i < n_grid; ++i) {
                    if (radii[i] != reference[k][i]) { ++failures; break; }
                }
            }
        }
    };

    std::vector<std::thread> pool;
    for (int t = 0; t < n_threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    beta_param_cache_destroy(cache);

    std::printf("thread_stress_test: %d mismatch(es)\n", failures.load());
    return failures.load() == 0 ? 0 : 1;
}
