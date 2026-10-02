!> Public Fortran API for the beta parameterization library.
!!
!! Two tiers (shape parameterization contract, two-tier revision):
!!
!!   - One-shot: `compute_*_standalone_s` computes one shape per call. All
!!     workspace is internal and discarded on return.
!!   - Read-only cache: `cache_t` holds everything determined by `max_params`
!!     and the primary theta set (normalization constants, Gauss-Legendre
!!     quadrature, Legendre tables). Built once, then shared.
!!
!! `node_set_t` is a caller-owned set of extra theta nodes with its own
!! Legendre tables, built from a cache.
!!
!! Every entry point reports through the shared status contract
!! (`SHAPE_*` codes, library codes >= 100); none of them stop.
!!
!! ## Usage pattern
!!
!! ```fortran
!! call cache_init_s(cache, max_params, thetas, status)        ! once
!! do i = 1, n_shapes                                          ! any thread
!!     call cache_radius_grid_s(cache, params(:, i), conserve_volume, &
!!             apply_com, radii, status)
!! end do
!! call cache_free_s(cache)
!! ```
!!
!! ## Thread model
!!
!! A cache and a node set are immutable after `cache_init_s` /
!! `node_set_build_s`. Every compute takes them `intent(in)` and is `pure`, so
!! any number of threads may compute on one cache concurrently. `cache_init_s`
!! and `cache_free_s` on a cache must not race with any other call on it. The
!! library holds no mutable module state; parameter-dependent intermediates are
!! per-call stack scratch.
!!
!! ## Parameters
!!
!! A cached call accepts `1 <= size(params) <= max_params`; a one-shot call
!! accepts `1 <= size(params) <= 64`. Missing trailing parameters are zero.
!! Both tiers trim trailing zero parameters before any arithmetic, so a short
!! vector and its zero-padded form run the same loops on the same data and
!! return identical bits. Options (`conserve_volume`, `apply_com`) are
!! per-call arguments: one cache serves every combination.
!!
!! ## Precondition: finite input
!!
!! `params` (and every theta) must be finite. Non-finite input is UNDEFINED
!! BEHAVIOR — the library cannot detect NaN under fast-math, so no check
!! rejects it, and no particular result is promised: a call may return
!! `SHAPE_VALID` (0) with NaN outputs; a NaN trailing parameter is
!! trimmed like a zero, giving the finite outputs of the shorter vector; a
!! Debug build may trap. No runtime NaN check exists on any path; screen
!! inputs before calling.
!!
!! Magnitudes must be physical as well. With `apply_com` the COM quadrature
!! evaluates R**4 before any validity gate, which overflows for |beta| beyond
!! about 1e70 — a floating-point trap in Debug builds, Inf arithmetic in
!! Release. Shapes anywhere near the valid domain have |beta| of order 1.
!!
!! ## Failure semantics
!!
!! Every output argument is zero-filled on entry and stays zero on any nonzero
!! status. Usage errors are checked in a fixed order and depend on sizes only,
!! never on parameter values: uninitialized cache (2), parameter count (4),
!! output buffer size (104), node set (106). The buffer check runs before the
!! node-set check: an unbuilt node set has zero nodes, so any non-empty buffer
!! trips `BETA_PARAM_ERROR_INVALID_BUFFER_SIZE` (104) and only a zero-length
!! buffer reaches `BETA_PARAM_ERROR_NODE_SET_MISMATCH` (106). Shape codes
!! follow: 103 (COM correction, only with `apply_com`), then 100/101/102 from
!! the validity gate.
module beta_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_utilities_mod, only: &
            compute_spherical_harmonics_normalization_constants_s, &
            compute_gauss_legendre_quadrature_s
    use shape_core_mod, only: &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_MAX_PARAMS
    use beta_parameterization_workers_mod, only: &
            precompute_legendre_table_s, &
            precompute_legendre_derivative_table_s, &
            eval_polar_radii_s, eval_radius_grid_s, eval_radius_derivative_s, &
            find_min_radius_s, &
            compute_com_integrals_s, newton_com_correction_s

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public types
    !---------------------------------------------------------------------------
    public :: cache_t
    public :: node_set_t

    !---------------------------------------------------------------------------
    ! Cache lifecycle
    !---------------------------------------------------------------------------
    public :: cache_init_s, cache_free_s
    public :: cache_max_params_f, cache_n_thetas_f, cache_is_initialized_f

    !---------------------------------------------------------------------------
    ! Node-set lifecycle
    !---------------------------------------------------------------------------
    public :: node_set_build_s, node_set_free_s, node_set_n_nodes_f

    !---------------------------------------------------------------------------
    ! Cached computes (tier 2)
    !---------------------------------------------------------------------------
    public :: cache_resolve_shape_s
    public :: cache_radius_grid_s, cache_radius_and_derivative_s
    public :: cache_node_radius_and_derivative_s
    public :: cache_radius_grid_unchecked_s

    !---------------------------------------------------------------------------
    ! One-shot computes (tier 1)
    !---------------------------------------------------------------------------
    public :: compute_radius_grid_standalone_s
    public :: compute_radius_and_derivative_standalone_s

    !---------------------------------------------------------------------------
    ! Public limits
    !---------------------------------------------------------------------------
    !> Highest Legendre order the library supports (the contract's N_max).
    integer(kind = ik), parameter, public :: MAX_BETA_PARAMS_LIMIT = 64_ik

    !> Longest parameter vector either tier accepts.
    integer(kind = ik), parameter :: PARAM_LIMIT = min(SHAPE_MAX_PARAMS, MAX_BETA_PARAMS_LIMIT)

    !> Fixed Gauss-Legendre quadrature order for volume/COM integrals.
    integer(kind = ik), parameter :: N_QUAD = 512_ik

    ! Shared contract codes re-exported for consumers
    public :: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS
    public :: SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID
    public :: SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT
    public :: SHAPE_MAX_PARAMS

    ! Library codes (contract range >= 100, append-only)
    integer(kind = ik), parameter, public :: BETA_PARAM_ERROR_NORTH_POLE          = 100_ik
    integer(kind = ik), parameter, public :: BETA_PARAM_ERROR_SOUTH_POLE          = 101_ik
    integer(kind = ik), parameter, public :: BETA_PARAM_ERROR_INTERIOR_NEGATIVE   = 102_ik
    integer(kind = ik), parameter, public :: BETA_PARAM_ERROR_COM_NOT_CONVERGED   = 103_ik
    integer(kind = ik), parameter, public :: BETA_PARAM_ERROR_INVALID_BUFFER_SIZE = 104_ik
    integer(kind = ik), parameter, public :: BETA_PARAM_ERROR_POLE_NODE           = 105_ik
    integer(kind = ik), parameter, public :: BETA_PARAM_ERROR_NODE_SET_MISMATCH   = 106_ik

    integer(kind = ik), parameter, public :: STATUS_MESSAGE_LEN = 64_ik
    public :: status_message_f

    !---------------------------------------------------------------------------
    ! Validation thresholds (policy — workers compute, API decides)
    !---------------------------------------------------------------------------
    real(kind = rk), parameter :: R_MIN_THRESHOLD = 1.0e-6_rk

    !> Read-only cache: everything that depends only on `max_params` and the
    !! primary theta set. Built once, then shared read-only by every compute and
    !! every thread — no shape data inside.
    !!
    !! Every component is allocatable or scalar, so intrinsic assignment is a
    !! deep copy and re-initializing a live cache releases the old contents
    !! (`intent(out)`); nothing can leak.
    !!
    !! Invariants:
    !!   - After a successful `cache_init_s`: is_initialized == .true. and every
    !!     allocatable component is allocated to its declared shape.
    !!   - A failed `cache_init_s` and `cache_free_s` leave the default
    !!     (uninitialized) state.
    type :: cache_t
        private
        logical            :: is_initialized = .false.
        integer(kind = ik) :: max_params = 0_ik
        integer(kind = ik) :: n_thetas   = 0_ik

        real(kind = rk), allocatable :: norm_constants(:)             ! (max_params)
        real(kind = rk), allocatable :: gl_nodes(:)                   ! (N_QUAD)
        real(kind = rk), allocatable :: gl_weights(:)                 ! (N_QUAD)
        real(kind = rk), allocatable :: legendre_gl(:, :)             ! (N_QUAD, max_params + 1)
        real(kind = rk), allocatable :: thetas(:)                     ! (n_thetas)
        real(kind = rk), allocatable :: sin_thetas(:)                 ! (n_thetas)
        real(kind = rk), allocatable :: legendre_primary(:, :)        ! (n_thetas, max_params + 1)
        real(kind = rk), allocatable :: legendre_primary_deriv(:, :)  ! (n_thetas, max_params + 1)
    end type cache_t

    !> A caller-owned set of theta nodes with precomputed Legendre P_k and P_k'
    !! tables sized to one cache's max_params. Built once (startup), then reused
    !! across shapes. Carries no shape data; immutable after build and safe to
    !! share across threads. Deliberately generic: no dense/folding/coulomb
    !! vocabulary — the caller owns what a set means.
    !!
    !! Lifetime (contract rule): the cache a node set was built from must
    !! outlive it. Free node sets first.
    !!
    !! Invariants:
    !!   - After a successful `node_set_build_s`: is_built == .true., every
    !!     allocatable component allocated, max_l == the source cache's
    !!     max_params.
    !!   - A failed build leaves the default (unbuilt) state.
    type :: node_set_t
        private
        logical            :: is_built = .false.
        integer(kind = ik) :: n_nodes  = 0_ik
        integer(kind = ik) :: max_l    = 0_ik
        real(kind = rk), allocatable :: thetas(:)
        real(kind = rk), allocatable :: sin_thetas(:)
        real(kind = rk), allocatable :: legendre_table(:, :)        ! P_k(cos theta_i)
        real(kind = rk), allocatable :: legendre_deriv_table(:, :)  ! P_k'(cos theta_i)
    end type node_set_t

contains

    !> Fixed diagnostic string for a status code.
    pure function status_message_f(status) result(msg)
        integer(kind = ik), intent(in) :: status
        character(len = STATUS_MESSAGE_LEN) :: msg
        select case (status)
        case (SHAPE_VALID);                          msg = 'valid'
        case (SHAPE_ERROR_TOO_MANY_PARAMS);          msg = 'too many parameters'
        case (SHAPE_ERROR_CACHE_NOT_INITIALIZED);    msg = 'cache not initialized'
        case (SHAPE_ERROR_INVALID_GRID);             msg = 'theta grid below minimum size (2)'
        case (SHAPE_ERROR_WRONG_PARAM_COUNT);        msg = 'params length outside 1..max_params'
        case (SHAPE_ERROR_INVALID_INIT);             msg = 'invalid init arguments'
        case (BETA_PARAM_ERROR_NORTH_POLE);          msg = 'north pole radius not positive'
        case (BETA_PARAM_ERROR_SOUTH_POLE);          msg = 'south pole radius not positive'
        case (BETA_PARAM_ERROR_INTERIOR_NEGATIVE);   msg = 'interior radius not positive'
        case (BETA_PARAM_ERROR_COM_NOT_CONVERGED);   msg = 'COM correction did not converge'
        case (BETA_PARAM_ERROR_INVALID_BUFFER_SIZE); msg = 'output buffer size mismatch'
        case (BETA_PARAM_ERROR_POLE_NODE);           msg = 'theta at or beyond a pole'
        case (BETA_PARAM_ERROR_NODE_SET_MISMATCH);   msg = 'node set unbuilt or max_l too small'
        case default;                                msg = 'unknown status code'
        end select
    end function status_message_f

    !===========================================================================
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Reject theta sets that cannot carry Legendre derivative tables.
    !!
    !! @param[in]  thetas  Candidate theta nodes (radians)
    !! @param[out] status  SHAPE_VALID, SHAPE_ERROR_INVALID_GRID or
    !!                     BETA_PARAM_ERROR_POLE_NODE
    pure subroutine validate_theta_set_s(thetas, status)
        real(kind = rk),    intent(in)  :: thetas(:)
        integer(kind = ik), intent(out) :: status
        integer(kind = ik) :: i
        status = SHAPE_VALID
        if (size(thetas, kind = ik) < 2_ik) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        do i = 1_ik, size(thetas, kind = ik)
            ! Guard on the rounded cosine, not theta: the derivative recurrence
            ! divides by 1 - x^2, and a tiny theta can round to cos(theta) == 1.
            if (1.0_rk - cos(thetas(i))**2 <= 0.0_rk) then
                status = BETA_PARAM_ERROR_POLE_NODE
                return
            end if
        end do
    end subroutine validate_theta_set_s

    !> Build the read-only cache for one `max_params` and one primary theta set.
    !!
    !! The Legendre order of every table is `max_params`. Re-initializing a
    !! live cache is safe: `intent(out)` releases the old contents first.
    !!
    !! @param[out] cache       Fully populated on success; default
    !!                         (uninitialized) state otherwise
    !! @param[in]  max_params  Longest parameter vector the cache accepts,
    !!                         1 <= max_params <= 64
    !! @param[in]  thetas      Primary theta nodes (radians), at least 2, none polar
    !! @param[out] status      SHAPE_VALID on success, else the rejecting code:
    !!                         5 (max_params < 1), 1 (max_params > 64),
    !!                         3 (fewer than 2 thetas, or a table that cannot
    !!                         be allocated), 105 (pole node)
    pure subroutine cache_init_s(cache, max_params, thetas, status)
        type(cache_t),      intent(out) :: cache
        integer(kind = ik), intent(in)  :: max_params
        real(kind = rk),    intent(in)  :: thetas(:)
        integer(kind = ik), intent(out) :: status
        real(kind = rk), allocatable :: x_primary(:)
        integer(kind = ik) :: i, n
        integer :: alloc_stat

        status = SHAPE_VALID
        if (max_params < 1_ik) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if
        if (max_params > PARAM_LIMIT) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if
        call validate_theta_set_s(thetas, status)
        if (status /= SHAPE_VALID) return

        n = size(thetas, kind = ik)
        allocate(cache%norm_constants(max_params), &
                cache%gl_nodes(N_QUAD), cache%gl_weights(N_QUAD), &
                cache%legendre_gl(N_QUAD, max_params + 1_ik), &
                cache%thetas(n), cache%sin_thetas(n), x_primary(n), &
                cache%legendre_primary(n, max_params + 1_ik), &
                cache%legendre_primary_deriv(n, max_params + 1_ik), &
                stat = alloc_stat)
        if (alloc_stat /= 0) then
            ! A multi-object allocate may have succeeded partway: release it.
            call cache_free_s(cache)
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        cache%max_params = max_params
        cache%n_thetas   = n
        call compute_spherical_harmonics_normalization_constants_s( &
                cache%norm_constants, max_params)
        ! ff signature is (n, nodes, weights); nodes come back DESCENDING in x:
        ! nodes(1) ~ +1 (theta ~ 0, north), nodes(N_QUAD) ~ -1 (theta ~ pi, south)
        call compute_gauss_legendre_quadrature_s(N_QUAD, cache%gl_nodes, cache%gl_weights)
        call precompute_legendre_table_s(cache%gl_nodes, max_params, cache%legendre_gl)

        cache%thetas = thetas
        do i = 1_ik, n
            x_primary(i) = cos(thetas(i))
            cache%sin_thetas(i) = sin(thetas(i))
        end do
        call precompute_legendre_table_s(x_primary, max_params, cache%legendre_primary)
        call precompute_legendre_derivative_table_s(x_primary, max_params, &
                cache%legendre_primary, cache%legendre_primary_deriv)
        cache%is_initialized = .true.
    end subroutine cache_init_s

    !> Release the cache. Infallible, and safe on an uninitialized cache:
    !! intent(out) deallocates every component and restores the defaults.
    pure subroutine cache_free_s(cache)
        type(cache_t), intent(out) :: cache
        cache%is_initialized = .false.   ! intent(out) already did this; explicit for clarity
    end subroutine cache_free_s

    !> Longest parameter vector the cache accepts; 0 when uninitialized.
    pure function cache_max_params_f(cache) result(max_params)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: max_params
        max_params = 0_ik
        if (cache%is_initialized) max_params = cache%max_params
    end function cache_max_params_f

    !> Number of primary theta nodes; 0 when uninitialized.
    pure function cache_n_thetas_f(cache) result(n)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: n
        n = 0_ik
        if (cache%is_initialized) n = cache%n_thetas
    end function cache_n_thetas_f

    !> .true. only after a successful init and before `cache_free_s`.
    pure function cache_is_initialized_f(cache) result(ok)
        type(cache_t), intent(in) :: cache
        logical :: ok
        ok = cache%is_initialized
    end function cache_is_initialized_f

    !===========================================================================
    ! NODE-SET LIFECYCLE
    !===========================================================================

    !> Precompute P_k and P_k' tables at caller-supplied theta nodes, sized to
    !! the source cache's max_params.
    !!
    !! Pole nodes are rejected (`BETA_PARAM_ERROR_POLE_NODE`): the derivative
    !! recurrence divides by 1 - x**2. Pole radii come analytically from
    !! `cache_resolve_shape_s` instead.
    !!
    !! @param[out] node_set  Filled on success; unbuilt otherwise
    !! @param[in]  cache     Initialized cache; supplies max_params. Must
    !!                       outlive the node set (contract lifetime rule).
    !! @param[in]  thetas    Node angles (radians); any order, need not be uniform
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code:
    !!                       2 (uninitialized cache), 3 (fewer than 2 thetas, or
    !!                       a table that cannot be allocated), 105 (pole node)
    pure subroutine node_set_build_s(node_set, cache, thetas, status)
        type(node_set_t),   intent(out) :: node_set
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: thetas(:)
        integer(kind = ik), intent(out) :: status
        real(kind = rk), allocatable :: x(:)
        integer(kind = ik) :: i, n
        integer :: alloc_stat

        status = SHAPE_VALID
        if (.not. cache%is_initialized) then
            status = SHAPE_ERROR_CACHE_NOT_INITIALIZED
            return
        end if
        call validate_theta_set_s(thetas, status)
        if (status /= SHAPE_VALID) return

        n = size(thetas, kind = ik)
        allocate(node_set%thetas(n), node_set%sin_thetas(n), x(n), &
                node_set%legendre_table(n, cache%max_params + 1_ik), &
                node_set%legendre_deriv_table(n, cache%max_params + 1_ik), &
                stat = alloc_stat)
        if (alloc_stat /= 0) then
            call node_set_free_s(node_set)
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        node_set%n_nodes = n
        node_set%max_l   = cache%max_params   ! recorded for the consumer-side mismatch check
        node_set%thetas  = thetas
        do i = 1_ik, n
            x(i) = cos(thetas(i))
            node_set%sin_thetas(i) = sin(thetas(i))
        end do
        call precompute_legendre_table_s(x, cache%max_params, node_set%legendre_table)
        call precompute_legendre_derivative_table_s(x, cache%max_params, &
                node_set%legendre_table, node_set%legendre_deriv_table)
        node_set%is_built = .true.
    end subroutine node_set_build_s

    !> Release the node set. Infallible: intent(out) deallocates every component
    !! and restores the default component values.
    pure subroutine node_set_free_s(node_set)
        type(node_set_t), intent(out) :: node_set
        node_set%is_built = .false.   ! intent(out) already did this; explicit for clarity
    end subroutine node_set_free_s

    !> Number of nodes in the set; 0 when unbuilt.
    pure function node_set_n_nodes_f(node_set) result(n)
        type(node_set_t), intent(in) :: node_set
        integer(kind = ik) :: n
        n = 0_ik
        if (node_set%is_built) n = node_set%n_nodes
    end function node_set_n_nodes_f

    !===========================================================================
    ! SHARED COMPUTE CORES (cache + plain arrays, nothing stored)
    !===========================================================================

    !> Normalize the parameters, optionally COM-correct them, take the poles.
    !!
    !! @param[in]    cache             Initialized cache (norms + GL data)
    !! @param[in]    apply_com         Apply the centre-of-mass correction
    !! @param[inout] beta_local        Parameters in, COM-corrected values out
    !! @param[out]   beta_con          beta_local(k) x norm_constants(k)
    !! @param[out]   corrected_beta10  beta_local(1) after the correction
    !! @param[out]   r_north           R(theta = 0), UNSCALED
    !! @param[out]   r_south           R(theta = pi), UNSCALED
    !! @param[out]   status            SHAPE_VALID or BETA_PARAM_ERROR_COM_NOT_CONVERGED
    pure subroutine resolve_core_s(cache, apply_com, beta_local, beta_con, &
            corrected_beta10, r_north, r_south, status)
        type(cache_t),      intent(in)    :: cache
        logical,            intent(in)    :: apply_com
        real(kind = rk),    intent(inout) :: beta_local(:)
        real(kind = rk),    intent(out)   :: beta_con(:)
        real(kind = rk),    intent(out)   :: corrected_beta10, r_north, r_south
        integer(kind = ik), intent(out)   :: status

        integer(kind = ik) :: n, n_iter
        logical            :: converged

        status           = SHAPE_VALID
        corrected_beta10 = 0.0_rk
        r_north          = 0.0_rk
        r_south          = 0.0_rk
        n = size(beta_local, kind = ik)

        beta_con(:) = beta_local(:) * cache%norm_constants(1:n)

        if (apply_com) then
            call newton_com_correction_s(beta_local, beta_con, &
                    cache%norm_constants(1:n), cache%gl_nodes, &
                    cache%gl_weights, cache%legendre_gl, converged, n_iter)
            if (.not. converged) then
                status = BETA_PARAM_ERROR_COM_NOT_CONVERGED
                return
            end if
        end if

        corrected_beta10 = beta_local(1)
        call eval_polar_radii_s(beta_con, r_north, r_south)
    end subroutine resolve_core_s

    !> Reject shapes whose radius is not positive everywhere.
    !!
    !! Poles come from the analytic values `resolve_core_s` produced; the
    !! interior is scanned on the Gauss-Legendre grid, which is DESCENDING in
    !! x — index 1 is x ~ +1 (theta ~ 0, north) and index N_QUAD is x ~ -1
    !! (theta ~ pi, south), so a minimum at either end is attributed to that pole.
    !!
    !! Every value here is unscaled: the volume factor is positive, so it cannot
    !! change any sign the check looks at.
    !!
    !! @param[in]  cache    Initialized cache (GL nodes + Legendre table)
    !! @param[in]  beta_con Resolved beta x norm products
    !! @param[in]  r_north  Unscaled north pole radius
    !! @param[in]  r_south  Unscaled south pole radius
    !! @param[out] status   SHAPE_VALID or the rejecting BETA_PARAM_ERROR_* code
    pure subroutine validate_core_s(cache, beta_con, r_north, r_south, status)
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: beta_con(:)
        real(kind = rk),    intent(in)  :: r_north, r_south
        integer(kind = ik), intent(out) :: status

        real(kind = rk)    :: r_gl(N_QUAD)
        real(kind = rk)    :: r_min
        integer(kind = ik) :: i_min

        status = SHAPE_VALID
        if (r_north <= R_MIN_THRESHOLD) then
            status = BETA_PARAM_ERROR_NORTH_POLE
            return
        end if
        if (r_south <= R_MIN_THRESHOLD) then
            status = BETA_PARAM_ERROR_SOUTH_POLE
            return
        end if

        call eval_radius_grid_s(beta_con, cache%legendre_gl, r_gl)
        call find_min_radius_s(r_gl, r_min, i_min)

        if (r_min <= R_MIN_THRESHOLD) then
            if (i_min == 1_ik) then
                status = BETA_PARAM_ERROR_NORTH_POLE
            else if (i_min == N_QUAD) then
                status = BETA_PARAM_ERROR_SOUTH_POLE
            else
                status = BETA_PARAM_ERROR_INTERIOR_NEGATIVE
            end if
        end if
    end subroutine validate_core_s

    !> The radial scale that restores the unit-sphere volume, or 1.
    !!
    !! Infallible — `validate_core_s` has already guaranteed R > 0 at every
    !! quadrature node, so the volume integral is positive and the cube root is
    !! real.
    !!
    !! @param[in]  cache            Initialized cache (GL data)
    !! @param[in]  beta_con         Resolved beta x norm products
    !! @param[in]  conserve_volume  .false. yields a factor of exactly 1
    !! @param[out] volume_factor    (2 / volume_integral)^(1/3), or 1
    pure subroutine volume_core_s(cache, beta_con, conserve_volume, volume_factor)
        type(cache_t),   intent(in)  :: cache
        real(kind = rk), intent(in)  :: beta_con(:)
        logical,         intent(in)  :: conserve_volume
        real(kind = rk), intent(out) :: volume_factor

        real(kind = rk) :: volume_integral, z_mean_integral

        if (conserve_volume) then
            call compute_com_integrals_s(beta_con, cache%gl_nodes, &
                    cache%gl_weights, cache%legendre_gl, &
                    volume_integral, z_mean_integral)
            volume_factor = (2.0_rk / volume_integral)**(1.0_rk / 3.0_rk)
        else
            volume_factor = 1.0_rk
        end if
    end subroutine volume_core_s

    !===========================================================================
    ! PER-CALL PIPELINE
    !===========================================================================

    !> The usage checks every cached compute shares, in contract order.
    !!
    !! @param[in]  cache   Cache the caller passed
    !! @param[in]  params  Parameter vector the caller passed
    !! @param[out] status  SHAPE_VALID, 2 (uninitialized cache) or 4 (length
    !!                     outside 1..max_params)
    pure subroutine check_call_s(cache, params, status)
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: params(:)
        integer(kind = ik), intent(out) :: status
        integer(kind = ik) :: n

        status = SHAPE_VALID
        if (.not. cache%is_initialized) then
            status = SHAPE_ERROR_CACHE_NOT_INITIALIZED
            return
        end if
        n = size(params, kind = ik)
        if (n < 1_ik .or. n > cache%max_params) status = SHAPE_ERROR_WRONG_PARAM_COUNT
    end subroutine check_call_s

    !> Index of the last nonzero parameter, at least 1.
    !!
    !! Trailing zeros are trimmed so that a short vector and its zero-padded
    !! form reach the kernels as the same data with the same trip count: their
    !! outputs are then bitwise identical whatever the compiler does with the
    !! summation loops. `abs(x) > 0` treats -0.0 as zero. Interior zeros stay.
    !!
    !! @param[in] params  Parameter vector, size >= 1
    pure function active_length_f(params) result(n_active)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik) :: n_active
        n_active = size(params, kind = ik)
        do while (n_active > 1_ik)
            if (abs(params(n_active)) > 0.0_rk) exit
            n_active = n_active - 1_ik
        end do
    end function active_length_f

    !> The front half of every compute: trim, resolve, and (checked paths only)
    !! validate and take the volume factor. Everything it produces is per-call
    !! scratch owned by the caller; nothing is stored.
    !!
    !! @param[in]  cache             Initialized cache, already through check_call_s
    !! @param[in]  params            Parameter vector, 1 <= size <= cache max_params
    !! @param[in]  conserve_volume   Compute the volume factor (checked paths)
    !! @param[in]  apply_com         Apply the centre-of-mass correction
    !! @param[in]  checked           .false. stops after the resolve stage
    !! @param[out] n_active          Active length after trimming
    !! @param[out] beta_con          Resolved beta x norm products in 1..n_active
    !! @param[out] corrected_beta10  beta10 after the COM correction
    !! @param[out] r_north           R(theta = 0), UNSCALED
    !! @param[out] r_south           R(theta = pi), UNSCALED
    !! @param[out] volume_factor     Radial scale; 1 when not conserved or unchecked
    !! @param[out] status            SHAPE_VALID or the first failing stage's code
    pure subroutine prepare_shape_s(cache, params, conserve_volume, apply_com, checked, &
            n_active, beta_con, corrected_beta10, r_north, r_south, volume_factor, status)
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: params(:)
        logical,            intent(in)  :: conserve_volume, apply_com, checked
        integer(kind = ik), intent(out) :: n_active
        real(kind = rk),    intent(out) :: beta_con(SHAPE_MAX_PARAMS)
        real(kind = rk),    intent(out) :: corrected_beta10, r_north, r_south, volume_factor
        integer(kind = ik), intent(out) :: status

        real(kind = rk) :: beta_local(SHAPE_MAX_PARAMS)

        volume_factor = 1.0_rk
        beta_con(:)   = 0.0_rk
        n_active = active_length_f(params)
        beta_local(1:n_active) = params(1:n_active)

        call resolve_core_s(cache, apply_com, beta_local(1:n_active), &
                beta_con(1:n_active), corrected_beta10, r_north, r_south, status)
        if (status /= SHAPE_VALID) return
        if (.not. checked) return

        call validate_core_s(cache, beta_con(1:n_active), r_north, r_south, status)
        if (status /= SHAPE_VALID) return
        call volume_core_s(cache, beta_con(1:n_active), conserve_volume, volume_factor)
    end subroutine prepare_shape_s

    !===========================================================================
    ! CACHED COMPUTES (TIER 2)
    !===========================================================================

    !> Resolve one shape: COM-corrected beta10, both pole radii, volume factor.
    !!
    !! The pole radii come out scaled by the volume factor (they are lengths);
    !! `corrected_beta10` is a deformation parameter and is NOT scaled. Every
    !! failure zero-fills all four outputs.
    !!
    !! @param[in]  cache             Initialized cache
    !! @param[in]  params            Parameter vector, 1 <= size <= cache max_params
    !! @param[in]  conserve_volume   Renormalize radii to fixed volume
    !! @param[in]  apply_com         Apply the centre-of-mass correction
    !! @param[out] corrected_beta10  beta10 after the COM correction
    !! @param[out] r_north           R(theta = 0) x volume_factor
    !! @param[out] r_south           R(theta = pi) x volume_factor
    !! @param[out] volume_factor     Radial scale (exactly 1 when not conserved)
    !! @param[out] status            SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_resolve_shape_s(cache, params, conserve_volume, apply_com, &
            corrected_beta10, r_north, r_south, volume_factor, status)
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: params(:)
        logical,            intent(in)  :: conserve_volume, apply_com
        real(kind = rk),    intent(out) :: corrected_beta10, r_north, r_south, volume_factor
        integer(kind = ik), intent(out) :: status

        real(kind = rk)    :: beta_con(SHAPE_MAX_PARAMS)
        real(kind = rk)    :: beta10, north, south, factor
        integer(kind = ik) :: n_active

        corrected_beta10 = 0.0_rk
        r_north          = 0.0_rk
        r_south          = 0.0_rk
        volume_factor    = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        call prepare_shape_s(cache, params, conserve_volume, apply_com, .true., &
                n_active, beta_con, beta10, north, south, factor, status)
        if (status /= SHAPE_VALID) return

        corrected_beta10 = beta10
        volume_factor    = factor
        r_north          = north * factor
        r_south          = south * factor
    end subroutine cache_resolve_shape_s

    !> R(theta) on the cache's primary theta grid, scaled by the volume factor.
    !!
    !! Every failure zero-fills `radii`.
    !!
    !! @param[in]  cache            Initialized cache
    !! @param[in]  params           Parameter vector, 1 <= size <= cache max_params
    !! @param[in]  conserve_volume  Renormalize radii to fixed volume
    !! @param[in]  apply_com        Apply the centre-of-mass correction
    !! @param[out] radii            R(theta_i) x volume_factor; size == cache n_thetas
    !! @param[out] status           SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_radius_grid_s(cache, params, conserve_volume, apply_com, &
            radii, status)
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: params(:)
        logical,            intent(in)  :: conserve_volume, apply_com
        real(kind = rk),    intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk)    :: beta_con(SHAPE_MAX_PARAMS)
        real(kind = rk)    :: beta10, north, south, factor
        integer(kind = ik) :: n_active

        radii(:) = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        if (size(radii, kind = ik) /= cache%n_thetas) then
            status = BETA_PARAM_ERROR_INVALID_BUFFER_SIZE
            return
        end if

        call prepare_shape_s(cache, params, conserve_volume, apply_com, .true., &
                n_active, beta_con, beta10, north, south, factor, status)
        if (status /= SHAPE_VALID) return

        call eval_radius_grid_s(beta_con(1:n_active), cache%legendre_primary, radii)
        radii(:) = radii(:) * factor
    end subroutine cache_radius_grid_s

    !> R(theta) and dR/dtheta on the primary theta grid, both volume-scaled.
    !!
    !! Both buffers are checked before any compute, and every failure zero-fills
    !! both.
    !!
    !! @param[in]  cache            Initialized cache
    !! @param[in]  params           Parameter vector, 1 <= size <= cache max_params
    !! @param[in]  conserve_volume  Renormalize radii to fixed volume
    !! @param[in]  apply_com        Apply the centre-of-mass correction
    !! @param[out] radii            R(theta_i) x volume_factor; size == cache n_thetas
    !! @param[out] dr_dthetas       dR/dtheta at theta_i x volume_factor; same size
    !! @param[out] status           SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_radius_and_derivative_s(cache, params, conserve_volume, &
            apply_com, radii, dr_dthetas, status)
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: params(:)
        logical,            intent(in)  :: conserve_volume, apply_com
        real(kind = rk),    intent(out) :: radii(:), dr_dthetas(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk)    :: beta_con(SHAPE_MAX_PARAMS)
        real(kind = rk)    :: beta10, north, south, factor
        integer(kind = ik) :: n_active

        radii(:)      = 0.0_rk
        dr_dthetas(:) = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        if (size(radii, kind = ik) /= cache%n_thetas .or. &
                size(dr_dthetas, kind = ik) /= cache%n_thetas) then
            status = BETA_PARAM_ERROR_INVALID_BUFFER_SIZE
            return
        end if

        call prepare_shape_s(cache, params, conserve_volume, apply_com, .true., &
                n_active, beta_con, beta10, north, south, factor, status)
        if (status /= SHAPE_VALID) return

        call eval_radius_grid_s(beta_con(1:n_active), cache%legendre_primary, radii)
        call eval_radius_derivative_s(beta_con(1:n_active), &
                cache%legendre_primary_deriv, cache%sin_thetas, dr_dthetas)
        radii(:)      = radii(:) * factor
        dr_dthetas(:) = dr_dthetas(:) * factor
    end subroutine cache_radius_and_derivative_s

    !> R(theta) and dR/dtheta at a caller-owned node set, both volume-scaled.
    !!
    !! The node set must be built and cover every vector the cache accepts
    !! (`node_set max_l >= cache max_params`, else
    !! `BETA_PARAM_ERROR_NODE_SET_MISMATCH`): the check depends on the two
    !! objects only, never on the vector, so a short vector and its zero-padded
    !! form get the same status. A node set built from this cache always
    !! passes. Buffers are checked against the node count first. Every failure
    !! zero-fills both buffers.
    !!
    !! @param[in]  cache            Initialized cache
    !! @param[in]  params           Parameter vector, 1 <= size <= cache max_params
    !! @param[in]  node_set         Built node set with max_l >= cache max_params
    !! @param[in]  conserve_volume  Renormalize radii to fixed volume
    !! @param[in]  apply_com        Apply the centre-of-mass correction
    !! @param[out] radii            R(theta_i) x volume_factor; size == node count
    !! @param[out] dr_dthetas       dR/dtheta at theta_i x volume_factor; same size
    !! @param[out] status           SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_node_radius_and_derivative_s(cache, params, node_set, &
            conserve_volume, apply_com, radii, dr_dthetas, status)
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: params(:)
        type(node_set_t),   intent(in)  :: node_set
        logical,            intent(in)  :: conserve_volume, apply_com
        real(kind = rk),    intent(out) :: radii(:), dr_dthetas(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk)    :: beta_con(SHAPE_MAX_PARAMS)
        real(kind = rk)    :: beta10, north, south, factor
        integer(kind = ik) :: n_active, n

        radii(:)      = 0.0_rk
        dr_dthetas(:) = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        n = node_set_n_nodes_f(node_set)
        if (size(radii, kind = ik) /= n .or. size(dr_dthetas, kind = ik) /= n) then
            status = BETA_PARAM_ERROR_INVALID_BUFFER_SIZE
            return
        end if
        if (.not. node_set%is_built .or. node_set%max_l < cache%max_params) then
            status = BETA_PARAM_ERROR_NODE_SET_MISMATCH
            return
        end if

        call prepare_shape_s(cache, params, conserve_volume, apply_com, .true., &
                n_active, beta_con, beta10, north, south, factor, status)
        if (status /= SHAPE_VALID) return

        call eval_radius_grid_s(beta_con(1:n_active), node_set%legendre_table, radii)
        call eval_radius_derivative_s(beta_con(1:n_active), &
                node_set%legendre_deriv_table, node_set%sin_thetas, dr_dthetas)
        radii(:)      = radii(:) * factor
        dr_dthetas(:) = dr_dthetas(:) * factor
    end subroutine cache_node_radius_and_derivative_s

    !> R(theta) on the primary grid with validation and volume scaling skipped.
    !!
    !! The rendering path: it resolves the parameters (COM correction included
    !! when asked for) and evaluates the Legendre sum, nothing more. A shape the
    !! validity check would reject comes back as a broken outline with
    !! SHAPE_VALID instead of an error code — which is the point, since a plot
    !! of the rejected shape is what explains the rejection.
    !!
    !! CAUTION: the radii are UNSCALED — this path never computes the volume
    !! factor, which is why it takes no `conserve_volume`. A caller that mixes
    !! this route with `cache_radius_grid_s(..., conserve_volume = .true., ...)`
    !! draws two outlines of different size for the same shape. Scale by the
    !! `volume_factor` from `cache_resolve_shape_s` if the sizes must agree.
    !!
    !! @param[in]  cache      Initialized cache
    !! @param[in]  params     Parameter vector, 1 <= size <= cache max_params
    !! @param[in]  apply_com  Apply the centre-of-mass correction
    !! @param[out] radii      R(theta_i), unscaled; size == cache n_thetas
    !! @param[out] status     SHAPE_VALID, or a usage code, or 103 (COM)
    pure subroutine cache_radius_grid_unchecked_s(cache, params, apply_com, radii, status)
        type(cache_t),      intent(in)  :: cache
        real(kind = rk),    intent(in)  :: params(:)
        logical,            intent(in)  :: apply_com
        real(kind = rk),    intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk)    :: beta_con(SHAPE_MAX_PARAMS)
        real(kind = rk)    :: beta10, north, south, factor
        integer(kind = ik) :: n_active

        radii(:) = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        if (size(radii, kind = ik) /= cache%n_thetas) then
            status = BETA_PARAM_ERROR_INVALID_BUFFER_SIZE
            return
        end if

        ! Resolve only: no validity scan, no volume integral.
        call prepare_shape_s(cache, params, .false., apply_com, .false., &
                n_active, beta_con, beta10, north, south, factor, status)
        if (status /= SHAPE_VALID) return

        call eval_radius_grid_s(beta_con(1:n_active), cache%legendre_primary, radii)
    end subroutine cache_radius_grid_unchecked_s

    !===========================================================================
    ! ONE-SHOT COMPUTES (TIER 1)
    !===========================================================================

    !> The usage checks both one-shot entries share, in contract order. Nothing
    !! is built until they all pass.
    !!
    !! @param[in]  n_params  size(params)
    !! @param[in]  n_thetas  size(thetas)
    !! @param[in]  n_radii   size(radii)
    !! @param[in]  n_derivs  size(dr_dthetas); pass n_radii again when the
    !!                       entry has no derivative buffer
    !! @param[out] status    SHAPE_VALID, 4 (empty), 1 (more than 64) or 104
    pure subroutine check_standalone_s(n_params, n_thetas, n_radii, n_derivs, status)
        integer(kind = ik), intent(in)  :: n_params, n_thetas, n_radii, n_derivs
        integer(kind = ik), intent(out) :: status

        status = SHAPE_VALID
        if (n_params < 1_ik) then
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
        else if (n_params > PARAM_LIMIT) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
        else if (n_radii /= n_thetas .or. n_derivs /= n_thetas) then
            status = BETA_PARAM_ERROR_INVALID_BUFFER_SIZE
        end if
    end subroutine check_standalone_s

    !> R(theta) on a caller-supplied theta set, with nothing kept between calls.
    !!
    !! Builds a local cache with `max_params = size(params)`, calls
    !! `cache_radius_grid_s` on it and discards it: one pipeline serves both
    !! tiers, so the results are identical to the cached ones bit for bit. The
    !! cost is one table build per call — use a `cache_t` for repeated
    !! evaluations. Every failure zero-fills `radii`.
    !!
    !! @param[in]  params           Parameter vector, 1 <= size <= 64
    !! @param[in]  thetas           Theta nodes (radians), at least 2, none polar
    !! @param[in]  conserve_volume  Renormalize radii to fixed volume
    !! @param[in]  apply_com        Apply the centre-of-mass correction
    !! @param[out] radii            R(theta_i) x volume_factor; size == size(thetas)
    !! @param[out] status           SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_radius_grid_standalone_s(params, thetas, conserve_volume, &
            apply_com, radii, status)
        real(kind = rk),    intent(in)  :: params(:)
        real(kind = rk),    intent(in)  :: thetas(:)
        logical,            intent(in)  :: conserve_volume
        logical,            intent(in)  :: apply_com
        real(kind = rk),    intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache
        integer(kind = ik) :: n_radii

        radii(:) = 0.0_rk
        n_radii  = size(radii, kind = ik)

        call check_standalone_s(size(params, kind = ik), size(thetas, kind = ik), &
                n_radii, n_radii, status)
        if (status /= SHAPE_VALID) return

        ! cache_init_s owns the theta-set contract (3 / 105).
        call cache_init_s(cache, size(params, kind = ik), thetas, status)
        if (status /= SHAPE_VALID) return

        call cache_radius_grid_s(cache, params, conserve_volume, apply_com, radii, status)
    end subroutine compute_radius_grid_standalone_s

    !> R(theta) and dR/dtheta on a caller-supplied theta set, nothing kept.
    !!
    !! The derivative twin of `compute_radius_grid_standalone_s`: same contract,
    !! same wrapper structure, both buffers checked before anything is built and
    !! both zero-filled on every failure.
    !!
    !! @param[in]  params           Parameter vector, 1 <= size <= 64
    !! @param[in]  thetas           Theta nodes (radians), at least 2, none polar
    !! @param[in]  conserve_volume  Renormalize radii to fixed volume
    !! @param[in]  apply_com        Apply the centre-of-mass correction
    !! @param[out] radii            R(theta_i) x volume_factor; size == size(thetas)
    !! @param[out] dr_dthetas       dR/dtheta at theta_i x volume_factor; same size
    !! @param[out] status           SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_radius_and_derivative_standalone_s(params, thetas, &
            conserve_volume, apply_com, radii, dr_dthetas, status)
        real(kind = rk),    intent(in)  :: params(:)
        real(kind = rk),    intent(in)  :: thetas(:)
        logical,            intent(in)  :: conserve_volume
        logical,            intent(in)  :: apply_com
        real(kind = rk),    intent(out) :: radii(:)
        real(kind = rk),    intent(out) :: dr_dthetas(:)
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        radii(:)      = 0.0_rk
        dr_dthetas(:) = 0.0_rk

        call check_standalone_s(size(params, kind = ik), size(thetas, kind = ik), &
                size(radii, kind = ik), size(dr_dthetas, kind = ik), status)
        if (status /= SHAPE_VALID) return

        call cache_init_s(cache, size(params, kind = ik), thetas, status)
        if (status /= SHAPE_VALID) return

        call cache_radius_and_derivative_s(cache, params, conserve_volume, apply_com, &
                radii, dr_dthetas, status)
    end subroutine compute_radius_and_derivative_standalone_s

end module beta_parameterization_mod
