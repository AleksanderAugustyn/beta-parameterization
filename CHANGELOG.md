# Changelog

All notable changes to this project are documented here. Versions follow
semantic versioning; the format follows [Keep a Changelog](https://keepachangelog.com).

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
- The documented dependency map lives in the module header of
  `src/beta_parameterization_mod.f08`: five intermediates, every mask = all
  `n_params` bits.

### Changed

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
  accepted up to 64; above → `SHAPE_ERROR_TOO_MANY_PARAMS`, never silent
  truncation. The `max_beta_params` padding argument is removed. Tier 1 never
  constructs an engine.
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
