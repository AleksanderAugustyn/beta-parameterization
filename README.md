# beta-parameterization

Axially symmetric nuclear shapes in the spherical-harmonic (beta) expansion: a Fortran 2018 library with a C API, a C++20 wrapper and a Python wheel. It implements the two-tier shape parameterization contract shared with its sibling libraries.

## The shape

The surface radius, in units of the spherical radius R₀, is

    R(θ) = c · [ 1 + Σ_{λ=1}^{n} β_λ · C_λ · P_λ(cos θ) ],    C_λ = √((2λ+1) / 4π)

so that `C_λ · P_λ(cos θ)` is the spherical harmonic Y_λ0. `params(λ)` is β_λ0; position λ always means order λ. Up to 64 orders are supported.

Two per-call options change the shape:

- `apply_com` — β₁ is replaced by the value that puts the centre of mass at the origin (Newton iteration, tolerance 10⁻⁵ R₀).
- `conserve_volume` — the scale `c = (2 / ∫ R³ d cos θ)^(1/3)` restores the volume of the unit sphere. Without it, `c = 1`.

A shape is valid when R > 10⁻⁶ at both poles and on an internal 512-node Gauss-Legendre grid. An invalid shape returns a status code and zero-filled outputs.

## Two tiers

1. **One-shot.** One call computes one shape. All workspace is internal and discarded on return. For scripts and one-off plots.
2. **Read-only cache.** Build a cache once for a theta grid, then request any number of shapes against it. The cache holds only what does not depend on the parameters (Legendre tables, quadrature). It is immutable after creation and may be shared across threads. For loops and hot paths.

Both tiers return bitwise-identical results. Nothing parameter-dependent is stored: every call is independent of the calls before it.

Rules common to both tiers:

- A cache accepts `1 .. max_params` parameters per call; a one-shot call accepts `1 .. 64`.
- Missing trailing parameters are zero. A short vector and its zero-padded form give identical bits.
- Thetas are in radians, at least two, none at a pole (the open grid `θ_i = i·π/(n+1)` is the usual choice). Pole radii come from `resolve_shape`.
- Inputs must be finite and of physical magnitude. Release builds use fast-math and cannot detect NaN, and with `apply_com` the quadrature overflows for |β| beyond about 10⁷⁰.
- Every output is zero-filled on a nonzero status.

### Fortran

```fortran
program beta_example
    use precision_utilities_mod, only: ik, rk
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_resolve_shape_s, &
            compute_radius_grid_standalone_s, SHAPE_VALID
    implicit none

    integer(kind = ik), parameter :: N = 180_ik
    real(kind = rk), parameter :: PI = 3.141592653589793_rk
    type(cache_t) :: cache
    real(kind = rk) :: thetas(N), radii(N)
    real(kind = rk) :: beta10, r_north, r_south, volume_factor
    integer(kind = ik) :: i, status

    do i = 1_ik, N
        thetas(i) = real(i, rk) * PI / real(N + 1_ik, rk)   ! open grid: no pole nodes
    end do

    ! Tier 1: one-shot. params = (beta1, beta2), both options on.
    call compute_radius_grid_standalone_s([0.0_rk, 0.25_rk], thetas, .true., .true., &
            radii, status)
    if (status /= SHAPE_VALID) error stop 'one-shot failed'

    ! Tier 2: build once, share, compute many.
    call cache_init_s(cache, 8_ik, thetas, status)          ! max_params = 8
    call cache_radius_grid_s(cache, [0.0_rk, 0.25_rk], .true., .true., radii, status)
    call cache_resolve_shape_s(cache, [0.0_rk, 0.25_rk], .true., .true., &
            beta10, r_north, r_south, volume_factor, status)
    call cache_free_s(cache)
end program beta_example
```

Cached outputs: `cache_radius_grid_s`, `cache_radius_and_derivative_s`, `cache_resolve_shape_s` (corrected β₁, pole radii, volume factor), `cache_node_radius_and_derivative_s` (radius and derivative at a `node_set_t`: extra thetas built from the cache with `node_set_build_s`), and `cache_radius_grid_unchecked_s` (no validity gate and no volume scaling: the outline of a shape the checked path rejects). To evaluate several theta grids in one call, build one node set over their concatenation.

### C

```c
#include "beta_parameterization.h"
#include <stdio.h>

int main(void) {
    enum { N = 180 };
    const double pi = 3.141592653589793;
    double thetas[N], radii[N];
    const double params[2] = {0.0, 0.25};
    for (int i = 0; i < N; ++i) thetas[i] = (i + 1) * pi / (N + 1);

    /* Tier 1: one-shot. */
    int status = beta_param_radius_grid_standalone(params, 2, thetas, N, 1, 1, radii);
    if (status != BETA_PARAM_VALID) {
        printf("%s\n", beta_param_status_message(status));
        return 1;
    }

    /* Tier 2: build once, share, compute many. */
    beta_param_cache_t* cache = beta_param_cache_create(8, thetas, N, &status);
    if (cache == NULL) {
        printf("%s\n", beta_param_status_message(status));
        return 1;
    }
    status = beta_param_cache_radius_grid(cache, params, 2, 1, 1, radii, N);
    beta_param_cache_destroy(cache);
    return status;
}
```

C++20 callers can use the RAII wrapper in `beta_parameterization.hpp` (`beta_param::Cache`, `beta_param::NodeSet`, `const` compute methods, `std::span` arguments).

### Python

```python
import beta_parameterization as bp

thetas = bp.theta_grid(180)                     # open grid, no pole nodes

# Tier 1: one-shot.
res = bp.radius_grid([0.0, 0.25], thetas, conserve_volume=True, apply_com=True)
assert res.ok, res.message

# Tier 2: build once, share, compute many.
with bp.Cache(8, thetas) as cache:              # max_params = 8
    res = cache.radius_and_derivative([0.0, 0.25], conserve_volume=True, apply_com=True)
    shape = cache.resolve_shape([0.0, 0.25], conserve_volume=True, apply_com=True)
    print(res.radii[:3], shape.volume_factor)
```

Shape-validation failures come back as results carrying a `Status`; `BetaParamError` is raised only for usage errors (failed create, closed handle, non-1-D input).

## Status codes

| Code | Name | Meaning |
|---|---|---|
| 0 | `SHAPE_VALID` | success |
| 1 | `SHAPE_ERROR_TOO_MANY_PARAMS` | `max_params` or a one-shot vector exceeds 64 |
| 2 | `SHAPE_ERROR_CACHE_NOT_INITIALIZED` | uninitialized cache (or NULL handle) |
| 3 | `SHAPE_ERROR_INVALID_GRID` | fewer than 2 thetas, or a table that cannot be allocated |
| 4 | `SHAPE_ERROR_WRONG_PARAM_COUNT` | vector length outside `1..max_params`, or an empty one-shot vector |
| 5 | `SHAPE_ERROR_INVALID_INIT` | `max_params < 1` |
| 100 | `BETA_PARAM_ERROR_NORTH_POLE` | radius not positive at θ = 0 |
| 101 | `BETA_PARAM_ERROR_SOUTH_POLE` | radius not positive at θ = π |
| 102 | `BETA_PARAM_ERROR_INTERIOR_NEGATIVE` | radius not positive in the interior |
| 103 | `BETA_PARAM_ERROR_COM_NOT_CONVERGED` | centre-of-mass correction did not converge |
| 104 | `BETA_PARAM_ERROR_INVALID_BUFFER_SIZE` | output buffer size does not match the theta or node count |
| 105 | `BETA_PARAM_ERROR_POLE_NODE` | a theta at or beyond a pole |
| 106 | `BETA_PARAM_ERROR_NODE_SET_MISMATCH` | node set unbuilt, or built for a smaller cache |

Codes 0–5 are shared by every shape parameterization library; the C header and the Python `Status` enum mirror all of them under `BETA_PARAM_*` names.

## Building

Requires GCC (gfortran, g++) and CMake ≥ 3.20. Dependencies are fetched by CMake.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

As a CMake dependency:

```cmake
FetchContent_Declare(
        beta-parameterization
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/beta-parameterization.git
        GIT_TAG 4.0.0
)
FetchContent_MakeAvailable(beta-parameterization)
target_link_libraries(my_target PRIVATE BetaParameterization::beta_parameterization)
```

Targets: `BetaParameterization::beta_parameterization` (static, Fortran), `::beta_parameterization_shared` (shared, C API), `::beta_parameterization_cxx` (header-only C++ wrapper over the shared library).

Python:

```bash
pip install beta-parameterization==4.0.0
```

The wheel is self-contained (manylinux2014, x86-64). To run the Python tests against a local build instead: `BETA_PARAM_LIB=$PWD/build/libbeta_parameterization.so PYTHONPATH=python python -m pytest python/tests`.

## Dependencies

| Dependency | Version |
|---|---|
| [fortran-foundations](https://github.com/AleksanderAugustyn/fortran-foundations) | 3.0.0 |
| [gcc-compiler-options](https://github.com/AleksanderAugustyn/gcc-compiler-options) | 2.0.0 |
| numpy (wheel only) | ≥ 1.21 |

## Design

- **Stateless computes.** A cache holds only parameter-independent tables. Every compute takes it read-only and is `pure`; per-call scratch lives on the stack. One cache serves every thread and every option combination.
- **One pipeline.** A one-shot call builds a local cache and runs the cached routine, so the two tiers cannot drift apart.
- **Reproducible bits.** The library's own objects are compiled without link-time optimization. Under `-flto -ffast-math` a kernel inlined into two call sites may round differently in each; without LTO every caller executes the same machine code, and equal inputs give equal bits.
- **No stops.** Every failure is a status code; nothing in the library calls `error stop`.
