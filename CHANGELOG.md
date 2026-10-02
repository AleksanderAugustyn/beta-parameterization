# Changelog

All notable changes to this project are documented here. Versions follow
semantic versioning; the format follows [Keep a Changelog](https://keepachangelog.com).

## 4.0.0

Adoption of the two-tier shape parameterization contract on
`shape_core_mod` (fortran-foundations 3.0.0). The incremental tier is gone:
one-shot functions plus one read-only cache, built once and shared across
threads. Every public surface — Fortran, C, C++, Python — changes; nothing in
3.x compiles against 4.0.0 unchanged.

### Changed (breaking)

- **One cache type.** `tables_t` is renamed `cache_t`; the mutable 3.x
  `cache_t`, `cache_init_shared_s`, the shared-tables mode and the `target`
  requirement are removed. `cache_init_s(cache, max_params, thetas, status)`.
  C: `beta_param_cache_create(max_params, thetas, n_thetas, status)`; the
  `tables` handle and `beta_param_cache_create_shared` are removed. C++:
  `Tables` is removed, `Cache(max_params, thetas)`. Python:
  `Cache(max_params, thetas)`.
- **Caches are immutable and shareable.** Every compute takes the cache
  read-only (`intent(in)`, `const` handle, `const` method) and is `pure` in
  Fortran. Any number of threads may compute on one cache. The 3.0.0
  thread-confinement rule is withdrawn.
- **Options are per call.** `conserve_volume` and `apply_com` move from cache
  creation to each compute, after `params` and before the outputs. One cache
  serves every combination. `cache_radius_grid_unchecked_s` takes `apply_com`
  only: it never volume-scales.
- **`max_params` replaces `n_params`.** A cache accepts any vector of
  `1 .. max_params` entries, up to 64 (the 8-parameter cap is gone). Missing
  trailing parameters are zero.
- **Node sets are built from a cache**:
  `node_set_build_s(node_set, cache, thetas, status)`; argument order of the
  node compute is `(cache, params, node_set, conserve_volume, apply_com, radii,
  dr_dthetas, status)`. A node set serves a cache when it was built from a
  cache with at least that cache's `max_params`
  (`BETA_PARAM_ERROR_NODE_SET_MISMATCH` otherwise) — judged on the two
  objects, no longer on the vector length.
- **Status codes.** Code 6 (`SHAPE_ERROR_TABLES_NOT_INITIALIZED`) is retired
  and never reused; an uninitialized cache passed to `node_set_build_s` now
  returns 2. An empty one-shot vector returns 4 (was 5). `max_params > 64` at
  init returns 1. C: a negative `n_thetas` returns 3 (was 5).
- **Getters.** `cache_n_params_f` becomes `cache_max_params_f`; C++
  `Cache::n_params()` becomes `max_params()`; Python `Cache.n_params` becomes
  `max_params`, and the `conserve_volume` / `apply_com` attributes are gone.

### Removed

- The recompute engine and everything built on it: `cache_recompute_count_f`,
  the `BETA_PARAM_I_*` indices, the dependency map, the minimality test suite.
- `SHAPE_CACHE_MAX_PARAMS`, `SHAPE_STANDALONE_MAX_PARAMS`
  (use `SHAPE_MAX_PARAMS`), `BETA_PARAM_CACHE_MAX_PARAMS`, Python
  `CACHE_MAX_PARAMS`, `Status.tables_not_initialized`.

### Added

- **Trailing zeros are trimmed.** Both tiers reduce a vector to its last
  nonzero entry before any arithmetic, so a short vector and its zero-padded
  form return identical bits. Interior zeros are kept.
- **README.md** with the shape definition, both tiers in Fortran, C and
  Python, build instructions and version pins.
- **CI on push and pull request** (`.github/workflows/tests.yml`): every suite
  in Debug and Release, pytest against the built library, and the
  manylinux2014 wheel build.
- Test suites for the contract's families: `equivalence` (one-shot ≡ cached,
  short ≡ zero-padded), `statelessness`, `boundary`, and a concurrency suite
  running eight threads on one shared cache.

### Fixed

- **COM correction on a non-positive volume integral.** With `apply_com`, a
  shape far outside the valid domain (for example `[0, -20]`) made the Newton
  step take the cube root of a negative number: a floating-point trap in Debug
  builds, and in Release a result that depended on a NaN comparison. The
  iteration now reports non-convergence (103). Without `apply_com` such a
  shape is still rejected by the validity gate (100).
- **C API stack overflow on a large wrong size.** Marshalling buffers were
  caller-sized automatic arrays; under Release a large wrong size argument
  overflowed the stack before the size was checked. They are heap-allocated
  now: the call returns `BETA_PARAM_ERROR_INVALID_BUFFER_SIZE`.
- **Allocation failure in cache or node-set build** returns
  `SHAPE_ERROR_INVALID_GRID` (3) instead of terminating the process.

### Build

- The library's own objects are compiled with `-fno-lto`. Under
  `-flto -ffast-math` the summation kernels were inlined into each caller and
  optimized per call site, so the same call could return different last bits
  from different call sites. Without LTO every caller executes one machine-code
  body per kernel, which is what makes the bitwise guarantees hold. Cost: about
  30% on shape evaluation time. Consumers keep LTO for their own code.
- fortran-foundations pin 2.4.0 → 3.0.0.

## 3.0.0

Adoption of the three-tier shape parameterization contract against
`shape_core_mod` (fortran-foundations 2.4.0). Every public surface — Fortran,
C, C++, Python — changes. Read the breaking-change list before upgrading;
nothing in 2.x compiles against 3.0.0 unchanged.

### BREAKING: thread-safety promise withdrawn

**The 2.x guarantee that one cache serves concurrent computes from many threads
is gone.** It was never true once the cache started tracking recompute state.

| object | 3.0.0 rule |
|---|---|
| `tables_t` / `beta_param_tables_t*` / `Tables` | immutable after creation; share across threads for concurrent reads |
| `node_set_t` / `beta_param_node_set_t*` / `NodeSet` | immutable after build; share across threads for concurrent reads |
| `cache_t` / `beta_param_cache_t*` / `Cache` / Python `Cache` | **THREAD-CONFINED**; every compute mutates it; concurrent use from more than one thread is undefined |

The supported pattern is one shared `tables_t` plus one cache per thread, built
with `cache_init_shared_s` / `beta_param_cache_create_shared()` /
`Tables`-taking `Cache` constructor. Callers that fanned out one 2.x cache
across an OpenMP region **must** hoist a per-thread cache; the old code will
race silently. C++ `Cache` compute methods are no longer `const` — the compiler
catches most of these at the call site.

**Fortran callers of `cache_init_shared_s` MUST declare their `tables_t` with
the `target` attribute.** The cache stores a pointer to it, and a pointer
associated with a non-target dummy becomes undefined when that procedure
returns (F2018 8.5.17). Without `target` the code compiles, links and usually
appears to work — then reads freed memory. Declare
`type(tables_t), target :: tables`.

### Added

- **In-library volume conservation.** `conserve_volume` is a cache/standalone
  creation flag. The factor `c = (2 / Σᵢ wᵢ Rᵢ³)^(1/3)` is computed on the
  internal GL-512 set over the COM-corrected, unscaled shape, after validation,
  and applied to radii, dR/dθ and the polar radii — never to β-space values.
  Consumers that rescaled shapes themselves (wmmm's
  `compute_volume_factor_from_dense_f`) should drop their copy.
- `apply_com` as a creation flag, replacing the separate
  `*_with_com_shift` entry points.
- `tables_t` — a shared, immutable tables level (normalization constants,
  Gauss-Legendre set, Legendre tables at the primary thetas) with
  `tables_init_s` / `tables_free_s` / `tables_max_l_f` / `tables_n_thetas_f`.
- `status_message_f(status)` (C: `beta_param_status_message`) — one pure lookup
  returning a fixed string per code.
- `cache_radius_grid_unchecked_s` — the rendering/diagnostic path: no validity
  gate, no volume scaling, usage errors only, so a rejected shape still yields
  its outline. Its radii are UNSCALED even on a `conserve_volume` cache; mixing
  it with the checked path draws two outlines of different size for one shape.
- `cache_recompute_count_f` plus the public intermediate indices
  `BETA_PARAM_I_RESOLVED`, `_I_MIN_RADIUS`, `_I_VOLUME`, `_I_RADII`, `_I_DERIV`
  — always-on recompute counters for cache-minimality testing.
- `BETA_PARAM_ERROR_NODE_SET_MISMATCH` (106) for an unbuilt node set or one
  whose `max_l` is below the cache's `n_params`.
- Public constants for the two tier caps, so callers can check before creating
  anything: `SHAPE_CACHE_MAX_PARAMS` (8) and `SHAPE_STANDALONE_MAX_PARAMS` (64)
  re-exported from Fortran, `BETA_PARAM_CACHE_MAX_PARAMS` (8) and
  `BETA_PARAM_MAX_PARAMS_LIMIT` (64) in the C header,
  `beta_param::cache_max_params` / `beta_param::max_params_limit` in the C++
  header, and `CACHE_MAX_PARAMS` (8) / `MAX_BETA_PARAMS_LIMIT` (64) in Python.
- The documented dependency map lives in the module header of
  `src/beta_parameterization_mod.f08`: five intermediates, every mask = all
  `n_params` bits.

### Changed

- **BREAKING: the cached tier now caps `n_params` at 8, down from 64.** A cache
  carries the recompute engine, whose per-intermediate dependency masks are
  fixed-width, so `shape_core`'s `SHAPE_CACHE_MAX_PARAMS` = 8 is the hard limit.
  `cache_init_s` / `cache_init_shared_s` reject `n_params > 8` with
  `SHAPE_ERROR_TOO_MANY_PARAMS` (1) and `n_params < 1` with
  `SHAPE_ERROR_INVALID_INIT` (5) — the engine check runs before any table is
  built, so nothing is allocated. **The standalone tier is unaffected and still
  accepts up to 64** (`SHAPE_STANDALONE_MAX_PARAMS`): it constructs no engine.
  2.x accepted `max_beta_params` up to 64 for the cache, so any consumer that
  created a cache with more than 8 parameters — wmmm's
  `number_of_deformation_parameters = 20` is the known case — fails at cache
  creation and must either reduce its parameter count or move to the standalone
  tier. Beta PES calculations realistically stay within 8 dimensions; wmmm's
  20 existed for fos→beta conversion, which is obsolete now that fos works on
  the radius grid directly.
- **BREAKING: the internal uniform grid is gone.** 2.x caches took `n_grid` and
  generated their own uniform theta grid. 3.0.0 takes an explicit **primary
  theta set** — `tables_init_s(tables, max_l, thetas, status)`,
  `cache_init_s(cache, n_params, thetas, ...)`. The caller owns the grid.
  Pole values (θ = 0, π) are rejected as cache input; get them analytically
  from `cache_resolve_shape_s`'s `r_north` / `r_south`.
- **BREAKING: all `message` arguments removed** from every entry point in
  Fortran, C, C++ and Python. Diagnostics are the integer status plus
  `status_message_f`. No per-call message formatting exists anywhere.
- **BREAKING (bitwise): the COM correction is now Newton's method** on the
  numerator of the COM integral (`N(β10) = Σᵢ wᵢ xᵢ Rᵢ⁴`,
  `N'(β10) = 4·C₁·Σᵢ wᵢ xᵢ² Rᵢ³`), replacing 2.x fixed-point iteration. The
  tolerance (`|z_cm| < 1e-5`) and iteration cap (50) are unchanged, so the same
  shapes converge — but **corrected β₁₀ and every downstream radius differ in
  the last bits from 2.x**. Goldens were re-baselined; consumers pinning exact
  values must re-baseline too.
- **BREAKING: `beta_con` is off the public surface.** 2.x
  `cache_resolve_shape` returned the resolved coefficient array and
  `compute_radius_and_derivative` took it back as an argument. 3.0.0 keeps
  resolved coefficients inside the cache as intermediate 1;
  `cache_resolve_shape_s` returns only `corrected_beta10`, `r_north`, `r_south`
  and `volume_factor`, and node evaluation takes `(cache, node_set, params)`.
- Type-bound methods (`cache%init`, `cache%compute_radius_grid`, …) replaced by
  free subroutines (`cache_init_s`, `cache_radius_grid_s`, …), matching the
  contract's naming.
- Standalone tier: `compute_radius_grid_standalone_s(params, thetas,
  conserve_volume, apply_com, radii, status)` and
  `compute_radius_and_derivative_standalone_s(...)`. `n_params = size(params)`,
  accepted up to the **standalone** cap of 64
  (`SHAPE_STANDALONE_MAX_PARAMS` = `MAX_BETA_PARAMS_LIMIT`); above →
  `SHAPE_ERROR_TOO_MANY_PARAMS`, never silent truncation. The
  `max_beta_params` padding argument is removed. Tier 1 never constructs an
  engine, which is why it keeps the higher cap.
- Failure semantics are now uniform: on ANY nonzero status inside a checked
  cached compute the library zero-fills every output and invalidates the whole
  engine, so the next call runs cold. Usage errors follow the contract's
  normative order — a buffer-size mismatch (104) is reported before a node-set
  mismatch (106).
- Minimum fortran-foundations is **2.4.0** (`shape_core_mod`); gcc-opts stays at
  **2.0.0**.

### Changed: status codes

Shared contract codes (0–99) come from `shape_core_mod`: `SHAPE_VALID` 0,
`SHAPE_ERROR_TOO_MANY_PARAMS` 1, `SHAPE_ERROR_CACHE_NOT_INITIALIZED` 2,
`SHAPE_ERROR_INVALID_GRID` 3 (theta floor: `size(thetas) >= 2`),
`SHAPE_ERROR_WRONG_PARAM_COUNT` 4, `SHAPE_ERROR_INVALID_INIT` 5,
`SHAPE_ERROR_TABLES_NOT_INITIALIZED` 6.

Library codes are renamed `LEGENDRE_*` → `BETA_PARAM_ERROR_*` and renumbered
≥ 100. This range is **append-only** from 3.0.0 on:

| old (2.x) | new (3.0.0) |
|---|---|
| `LEGENDRE_ERROR_NORTH_POLE` 1 | `BETA_PARAM_ERROR_NORTH_POLE` 100 |
| `LEGENDRE_ERROR_SOUTH_POLE` 2 | `BETA_PARAM_ERROR_SOUTH_POLE` 101 |
| `LEGENDRE_ERROR_INTERIOR_NEGATIVE` 4 | `BETA_PARAM_ERROR_INTERIOR_NEGATIVE` 102 |
| `LEGENDRE_ERROR_COM_NOT_CONVERGED` 7 | `BETA_PARAM_ERROR_COM_NOT_CONVERGED` 103 |
| `LEGENDRE_ERROR_INVALID_BUFFER_SIZE` 8 | `BETA_PARAM_ERROR_INVALID_BUFFER_SIZE` 104 |
| `LEGENDRE_ERROR_POLE_NODE` 9 | `BETA_PARAM_ERROR_POLE_NODE` 105 |
| — (new) | `BETA_PARAM_ERROR_NODE_SET_MISMATCH` 106 |
| `LEGENDRE_ERROR_EMPTY_PARAMS` 3 | → shared 4 (compute) / 5 (init) |
| `LEGENDRE_ERROR_INVALID_MAX_PARAMS` 5 | → shared 5 / 3 / 2 by cause |
| `LEGENDRE_ERROR_TOO_MANY_PARAMS` 6 | → shared 1 |
| `LEGENDRE_ERROR_NO_UNIFORM_GRID` 10 | removed (no uniform grid) |

Note that the old and new numbering overlap: 2.x code 1 meant "north pole",
3.0.0 code 1 means "too many params". Callers comparing raw integers must be
updated, not merely recompiled.

### Changed: C API

Handle-based, with statuses reported through out-parameters instead of
in-band returns:

- `beta_param_tables_create(max_l, thetas, n_thetas, int* status)` /
  `beta_param_tables_destroy` — new handle type `beta_param_tables_t*`.
- `beta_param_cache_create(n_params, thetas, n_thetas, conserve_volume,
  apply_com, int* status)` and
  `beta_param_cache_create_shared(tables, n_params, conserve_volume, apply_com,
  int* status)` — both return `NULL` on failure and write the rejecting code to
  the nullable `status` out-param. The 2.x `n_grid` + `message` arguments are
  gone.
- Computes drop the `compute_` infix and the `_with_com_shift` variants:
  `beta_param_cache_radius_grid`, `_radius_and_derivative`,
  `_radius_grid_unchecked`, `_resolve_shape`, `_node_radius_and_derivative`,
  `beta_param_radius_grid_standalone`, `_radius_and_derivative_standalone`.
  Each returns the status code directly.
- `beta_param_status_message(int)` returns a static string; no caller-supplied
  message buffers anywhere.
- Lifetime rule: a tables handle passed to `beta_param_cache_create_shared()`
  or `beta_param_node_set_create()` must outlive everything built from it;
  `beta_param_cache_destroy()` never frees shared tables.

### Changed: C++ header

`beta_parameterization.hpp` (C++20, header-only) exposes `Tables`, `Cache` and
`NodeSet` RAII wrappers. `Cache` is move-only and thread-confined; its compute
methods are non-`const`. Status strings come from `beta_param_status_message`.
The 2.x concurrent-compute documentation is replaced, not amended.

### Changed: Python wheel

- **BREAKING: `NodeSet` is removed.** BetaRender's single-set usage maps onto
  the primary theta set. Node sets return to Python only when a consumer needs
  a genuine second set — the same principle that keeps `tables_t` out of
  Python.
- **BREAKING: `theta_grid(n)` is now the OPEN uniform grid**
  `theta_i = i·π/(n+1), i = 1..n`. The 2.x closed `linspace(0, π, n)` includes
  both poles, which the pole guard now rejects as cache input. Close plotted
  curves with the analytic `r_north` / `r_south` from `resolve_shape`.
- Module-level tier 1: `radius_grid(params, thetas, conserve_volume=False,
  apply_com=False)`, `radius_and_derivative(...)`.
- `Cache(n_params, thetas, conserve_volume=False, apply_com=False)` with
  `.radius_grid`, `.radius_and_derivative`, `.radius_grid_unchecked`,
  `.resolve_shape`. Thread-confined.
- **Returns are result objects, not exceptions**: `RadiusGridResult`,
  `RadiusDerivativeResult`, `ResolvedShape`, each with `.ok`, `.status` and a
  `.message` property. Shape-validation failures come back as a status-carrying
  result; `BetaParamError` is raised only for usage errors (closed handle,
  non-1-D params, create failure).
- The `Status` enum is renumbered per the table above.

### Removed

- Concurrent-compute-on-one-cache support (see the breaking section above).
- The internally generated uniform theta grid and
  `LEGENDRE_ERROR_NO_UNIFORM_GRID`.
- All `message` output arguments and the message-formatting code behind them.
- `*_with_com_shift` entry points (Fortran, C), superseded by the `apply_com`
  creation flag.
- `beta_con` from every public signature.
- The `max_beta_params` padding argument on the standalone entry points.
- Python `NodeSet`.

### Migration sketch

```fortran
! 2.x
call cache%init(max_beta_params = 4_ik, n_grid = 80_ik, error_code = ec, message = msg)
call cache%compute_radius_grid_with_com_shift(params, radii, ec, msg)
call cache%destroy()

! 3.0.0 — the caller owns the theta set; flags are set once, at creation
thetas = [(i * PI / 81.0_rk, i = 1, 80)]          ! open grid: no poles
call cache_init_s(cache, n_params = 4_ik, thetas = thetas, &
        conserve_volume = .false., apply_com = .true., status = status)
call cache_radius_grid_s(cache, params, radii, status)
call cache_free_s(cache)
```

## 2.3.4 and earlier

No changelog was kept before 3.0.0. See the git history.
