!> Public Fortran API for the beta parameterization library.
!!
!! Three levels, from shared to per-shape:
!!
!!   - `tables_t`   — shared immutable level: everything determined by `max_l`
!!                    and the primary theta set. Built once, shared read-only.
!!   - `node_set_t` — caller-owned extra theta nodes with their own Legendre
!!                    tables, sized to one `tables_t`.
!!   - `cache_t`    — per-shape working level: one owner, one thread. Holds the
!!                    `shape_engine_t` recompute tracker plus every per-shape
!!                    buffer.
!!
!! Every entry point reports through the shared status contract
!! (`SHAPE_*` codes, library codes >= 100); none of them stop.
module beta_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_utilities_mod, only: &
            compute_spherical_harmonics_normalization_constants_s, &
            compute_gauss_legendre_quadrature_s
    use shape_core_mod, only: &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
            SHAPE_CACHE_MAX_PARAMS, SHAPE_STANDALONE_MAX_PARAMS, &
            shape_engine_t, shape_engine_init_s
    use beta_parameterization_workers_mod, only: &
            precompute_legendre_table_s, &
            precompute_legendre_derivative_table_s

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public types
    !---------------------------------------------------------------------------
    public :: cache_t
    public :: tables_t
    public :: node_set_t

    !---------------------------------------------------------------------------
    ! Tables lifecycle
    !---------------------------------------------------------------------------
    public :: tables_init_s, tables_free_s, tables_max_l_f, tables_n_thetas_f

    !---------------------------------------------------------------------------
    ! Node-set lifecycle
    !---------------------------------------------------------------------------
    public :: node_set_build_s, node_set_free_s, node_set_n_nodes_f

    !---------------------------------------------------------------------------
    ! Cache lifecycle
    !---------------------------------------------------------------------------
    public :: cache_init_s, cache_init_shared_s, cache_free_s
    public :: cache_n_params_f, cache_n_thetas_f, cache_is_initialized_f

    !---------------------------------------------------------------------------
    ! Public limits
    !---------------------------------------------------------------------------
    integer(kind = ik), parameter, public :: MAX_BETA_PARAMS_LIMIT = 64_ik

    !> Fixed Gauss-Legendre quadrature order for volume/COM integrals.
    integer(kind = ik), parameter :: N_QUAD = 512_ik

    !---------------------------------------------------------------------------
    ! Cached intermediates tracked by the embedded shape_engine_t
    !---------------------------------------------------------------------------
    integer(kind = ik), parameter :: I_RESOLVED = 1_ik, I_MIN_RADIUS = 2_ik, &
            I_VOLUME = 3_ik, I_RADII = 4_ik, I_DERIV = 5_ik, N_INTERMEDIATES = 5_ik

    ! Shared contract codes re-exported for consumers
    public :: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS
    public :: SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID
    public :: SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT
    public :: SHAPE_ERROR_TABLES_NOT_INITIALIZED
    public :: SHAPE_CACHE_MAX_PARAMS, SHAPE_STANDALONE_MAX_PARAMS

    ! Library codes (contract range >= 100, append-only after 3.0.0)
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

    !> Shared immutable level: everything that depends only on `max_l` and the
    !! primary theta set. Built once, then shared read-only by every consumer
    !! (node sets, caches, standalone entry points) — no shape data inside.
    !!
    !! A `tables_t` that will be shared with `cache_init_shared_s` MUST be
    !! declared with the `target` attribute by the caller: the cache keeps a
    !! pointer to it, and without `target` the association is undefined.
    !!
    !! Invariants:
    !!   - After a successful `tables_init_s`: is_initialized == .true. and every
    !!     allocatable component is allocated to its declared shape.
    !!   - `tables_free_s` restores the default (uninitialized) state.
    type :: tables_t
        private
        logical            :: is_initialized = .false.
        integer(kind = ik) :: max_l    = 0_ik
        integer(kind = ik) :: n_thetas = 0_ik

        real(kind = rk), allocatable :: norm_constants(:)             ! (max_l)
        real(kind = rk), allocatable :: gl_nodes(:)                   ! (N_QUAD)
        real(kind = rk), allocatable :: gl_weights(:)                 ! (N_QUAD)
        real(kind = rk), allocatable :: legendre_gl(:, :)             ! (N_QUAD, max_l + 1)
        real(kind = rk), allocatable :: thetas(:)                     ! (n_thetas)
        real(kind = rk), allocatable :: sin_thetas(:)                 ! (n_thetas)
        real(kind = rk), allocatable :: legendre_primary(:, :)        ! (n_thetas, max_l + 1)
        real(kind = rk), allocatable :: legendre_primary_deriv(:, :)  ! (n_thetas, max_l + 1)
    end type tables_t

    !> A caller-owned set of theta nodes with precomputed Legendre P_k and P_k'
    !! tables sized to one `tables_t`'s max_l. Built once (startup), then reused
    !! across shapes. Carries no shape data; immutable after build and safe to
    !! share across threads read-only. Deliberately generic: no dense/folding/
    !! coulomb vocabulary — the caller owns what a set means.
    !!
    !! Invariants:
    !!   - After a successful `node_set_build_s`: is_built == .true., every
    !!     allocatable component allocated, max_l == the source tables' max_l.
    !!   - A failed build leaves the default (unbuilt) state — `intent(out)`.
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

    !> Per-shape working level: the embedded `shape_engine_t` recompute tracker,
    !! the resolved coefficients, and every per-shape buffer.
    !!
    !! ## Ownership
    !!
    !! `tp` points at the shared tables. Two modes:
    !!   - `cache_init_s` (private): the cache heap-allocates its own `tables_t`
    !!     through `tp` and sets `owns_tables = .true.`; `cache_free_s` frees it.
    !!   - `cache_init_shared_s`: `tp` points at a caller-owned `tables_t` and
    !!     `owns_tables = .false.`. **The caller MUST declare that `tables_t`
    !!     with the `target` attribute** and must keep it alive (and not free it)
    !!     for the whole lifetime of the cache.
    !!
    !! ## `cache_free_s` is mandatory — there is no finalizer
    !!
    !! `cache_t` has no `final` binding, and both init routines take
    !! `intent(out) :: cache`, which default-initializes the cache on entry and
    !! so nulls `tp` before the previous target can be released. A private-mode
    !! cache therefore leaks its heap `tables_t` (and every allocatable inside
    !! it) if you either
    !!   - re-initialize a live cache (`cache_init_s(c, 4, ...)` then
    !!     `cache_init_s(c, 6, ...)` with no intervening free), or
    !!   - let the cache go out of scope without freeing it.
    !!
    !! Call `cache_free_s` before every re-initialization and before the cache
    !! goes out of scope. (A `final` binding is not the fix: `cache_free_s`
    !! resets through `cache = cache_t()`, whose LHS finalization would recurse.)
    !!
    !! ## Never copy-assign a cache_t
    !!
    !! Intrinsic assignment (`b = a`, passing by value, storing in an array that
    !! gets reallocated) copies `tp` and `owns_tables` shallowly, producing a
    !! dangling pointer or a double free. One cache has one owner and belongs to
    !! one thread; share the `tables_t` instead, not the cache.
    !!
    !! Invariants:
    !!   - After a successful init: is_initialized == .true., `tp` associated,
    !!     `radii`/`dr_dthetas` allocated to the primary theta count.
    !!   - Any failed init leaves the default state (both init routines take
    !!     `intent(out) :: cache`): is_initialized == .false., `tp` null,
    !!     owns_tables == .false.
    type :: cache_t
        private
        logical            :: is_initialized  = .false.
        integer(kind = ik) :: n_params        = 0_ik
        logical            :: conserve_volume = .false.
        logical            :: apply_com       = .false.

        type(shape_engine_t) :: engine                  !! Recompute tracker

        type(tables_t), pointer :: tp => null()         !! Shared tables (see Ownership)
        logical :: owns_tables = .false.                !! .true. => cache_free_s frees tp

        real(kind = rk) :: beta_local(SHAPE_CACHE_MAX_PARAMS) = 0.0_rk
        real(kind = rk) :: beta_con(SHAPE_CACHE_MAX_PARAMS)   = 0.0_rk
        real(kind = rk) :: corrected_beta10 = 0.0_rk
        real(kind = rk) :: r_north          = 0.0_rk
        real(kind = rk) :: r_south          = 0.0_rk
        real(kind = rk) :: r_min            = 0.0_rk
        integer(kind = ik) :: i_min         = 0_ik
        real(kind = rk) :: volume_factor    = 1.0_rk

        real(kind = rk), allocatable :: radii(:)        ! (tables n_thetas)
        real(kind = rk), allocatable :: dr_dthetas(:)   ! (tables n_thetas)
    end type cache_t

contains

    !> Fixed diagnostic string for a status code (spec 3.5).
    pure function status_message_f(status) result(msg)
        integer(kind = ik), intent(in) :: status
        character(len = STATUS_MESSAGE_LEN) :: msg
        select case (status)
        case (SHAPE_VALID);                          msg = 'valid'
        case (SHAPE_ERROR_TOO_MANY_PARAMS);          msg = 'too many parameters for this tier'
        case (SHAPE_ERROR_CACHE_NOT_INITIALIZED);    msg = 'cache not initialized'
        case (SHAPE_ERROR_INVALID_GRID);             msg = 'theta grid below minimum size (2)'
        case (SHAPE_ERROR_WRONG_PARAM_COUNT);        msg = 'params length differs from n_params'
        case (SHAPE_ERROR_INVALID_INIT);             msg = 'invalid init arguments'
        case (SHAPE_ERROR_TABLES_NOT_INITIALIZED);   msg = 'tables not initialized'
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
    ! TABLES LIFECYCLE
    !===========================================================================

    !> Reject theta sets that cannot carry Legendre derivative tables.
    !!
    !! @param[in]  thetas  Candidate theta nodes (radians)
    !! @param[out] status  SHAPE_VALID, SHAPE_ERROR_INVALID_GRID or
    !!                     BETA_PARAM_ERROR_POLE_NODE
    subroutine validate_theta_set_s(thetas, status)
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

    !> Build the shared immutable level for one `max_l` and one primary theta set.
    !!
    !! @param[out] tables  Fully populated on success; untouched-by-default state
    !!                     otherwise (intent(out) resets it on entry)
    !! @param[in]  max_l   1 <= max_l <= MAX_BETA_PARAMS_LIMIT
    !! @param[in]  thetas  Primary theta nodes (radians), at least 2, none polar
    !! @param[out] status  SHAPE_VALID on success, else the rejecting code
    subroutine tables_init_s(tables, max_l, thetas, status)
        type(tables_t),     intent(out) :: tables
        integer(kind = ik), intent(in)  :: max_l
        real(kind = rk),    intent(in)  :: thetas(:)
        integer(kind = ik), intent(out) :: status
        real(kind = rk), allocatable :: x_primary(:)
        integer(kind = ik) :: i, n

        status = SHAPE_VALID
        if (max_l < 1_ik .or. max_l > MAX_BETA_PARAMS_LIMIT) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if
        call validate_theta_set_s(thetas, status)
        if (status /= SHAPE_VALID) return

        n = size(thetas, kind = ik)
        tables%max_l = max_l
        tables%n_thetas = n
        allocate(tables%norm_constants(max_l))
        call compute_spherical_harmonics_normalization_constants_s( &
                tables%norm_constants, max_l)
        allocate(tables%gl_nodes(N_QUAD), tables%gl_weights(N_QUAD))
        ! ff signature is (n, nodes, weights); nodes come back DESCENDING in x:
        ! nodes(1) ~ +1 (theta ~ 0, north), nodes(N_QUAD) ~ -1 (theta ~ pi, south)
        call compute_gauss_legendre_quadrature_s(N_QUAD, tables%gl_nodes, tables%gl_weights)
        allocate(tables%legendre_gl(N_QUAD, max_l + 1_ik))
        call precompute_legendre_table_s(tables%gl_nodes, max_l, tables%legendre_gl)

        allocate(tables%thetas(n), tables%sin_thetas(n), x_primary(n))
        tables%thetas = thetas
        do i = 1_ik, n
            x_primary(i) = cos(thetas(i))
            tables%sin_thetas(i) = sin(thetas(i))
        end do
        allocate(tables%legendre_primary(n, max_l + 1_ik))
        allocate(tables%legendre_primary_deriv(n, max_l + 1_ik))
        call precompute_legendre_table_s(x_primary, max_l, tables%legendre_primary)
        call precompute_legendre_derivative_table_s(x_primary, max_l, &
                tables%legendre_primary, tables%legendre_primary_deriv)
        tables%is_initialized = .true.
    end subroutine tables_init_s

    !> Release the tables. Infallible: intent(out) deallocates every component
    !! and restores the default component values.
    pure subroutine tables_free_s(tables)
        type(tables_t), intent(out) :: tables
        tables%is_initialized = .false.   ! intent(out) already did this; explicit for clarity
    end subroutine tables_free_s

    !> Highest Legendre order the tables were built for; 0 when uninitialized.
    pure function tables_max_l_f(tables) result(max_l)
        type(tables_t), intent(in) :: tables
        integer(kind = ik) :: max_l
        max_l = 0_ik
        if (tables%is_initialized) max_l = tables%max_l
    end function tables_max_l_f

    !> Number of primary theta nodes; 0 when uninitialized.
    pure function tables_n_thetas_f(tables) result(n)
        type(tables_t), intent(in) :: tables
        integer(kind = ik) :: n
        n = 0_ik
        if (tables%is_initialized) n = tables%n_thetas
    end function tables_n_thetas_f

    !===========================================================================
    ! NODE-SET LIFECYCLE
    !===========================================================================

    !> Precompute P_k and P_k' tables at caller-supplied theta nodes, sized to
    !! the source tables' max_l.
    !!
    !! Pole nodes are rejected (`BETA_PARAM_ERROR_POLE_NODE`): the derivative
    !! recurrence divides by 1 - x**2. Pole radii come analytically from the
    !! resolve step instead.
    !!
    !! @param[out] node_set  Filled on success; unbuilt otherwise (intent(out)
    !!                       resets it on entry)
    !! @param[in]  tables    Initialized shared tables; supplies max_l
    !! @param[in]  thetas    Node angles (radians); any order, need not be uniform
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    subroutine node_set_build_s(node_set, tables, thetas, status)
        type(node_set_t),   intent(out) :: node_set
        type(tables_t),     intent(in)  :: tables
        real(kind = rk),    intent(in)  :: thetas(:)
        integer(kind = ik), intent(out) :: status
        real(kind = rk), allocatable :: x(:)
        integer(kind = ik) :: i, n

        status = SHAPE_VALID
        if (.not. tables%is_initialized) then
            status = SHAPE_ERROR_TABLES_NOT_INITIALIZED
            return
        end if
        call validate_theta_set_s(thetas, status)
        if (status /= SHAPE_VALID) return

        n = size(thetas, kind = ik)
        node_set%n_nodes = n
        node_set%max_l = tables%max_l   ! recorded for the consumer-side mismatch check
        allocate(node_set%thetas(n), node_set%sin_thetas(n), x(n))
        node_set%thetas = thetas
        do i = 1_ik, n
            x(i) = cos(thetas(i))
            node_set%sin_thetas(i) = sin(thetas(i))
        end do
        allocate(node_set%legendre_table(n, tables%max_l + 1_ik))
        allocate(node_set%legendre_deriv_table(n, tables%max_l + 1_ik))
        call precompute_legendre_table_s(x, tables%max_l, node_set%legendre_table)
        call precompute_legendre_derivative_table_s(x, tables%max_l, &
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
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Dependency masks for the embedded engine: every intermediate depends on
    !! every parameter, so each mask is the low `n_params` bits set.
    !!
    !! Out-of-range `n_params` yields an all-zero mask: `ishft` past the integer
    !! width is not defined, and the engine rejects the count before the mask
    !! value can matter.
    !!
    !! @param[in]  n_params  Requested parameter count (may be out of range)
    !! @param[out] masks     One mask per tracked intermediate
    pure subroutine engine_masks_s(n_params, masks)
        integer(kind = ik), intent(in)  :: n_params
        integer(kind = ik), intent(out) :: masks(N_INTERMEDIATES)
        masks(:) = 0_ik
        if (n_params >= 1_ik .and. n_params <= SHAPE_CACHE_MAX_PARAMS) then
            masks(:) = int(ishft(1_ik, n_params) - 1_ik, ik)
        end if
    end subroutine engine_masks_s

    !> Initialize a cache that owns its tables (private mode).
    !!
    !! The tables are heap-allocated through `cache%tp` and released by
    !! `cache_free_s`. Use `cache_init_shared_s` when many caches should share
    !! one `tables_t`.
    !!
    !! **Call `cache_free_s` first when re-initializing a live cache, and again
    !! before the cache goes out of scope.** There is no finalizer, and
    !! `intent(out) :: cache` nulls `tp` on entry, so re-initializing without a
    !! free leaks the previous heap `tables_t` and everything inside it.
    !!
    !! @param[out] cache            Ready on success; default (uninitialized)
    !!                              state on any failure
    !! @param[in]  n_params         1 <= n_params <= SHAPE_CACHE_MAX_PARAMS
    !! @param[in]  thetas           Primary theta nodes (radians), at least 2,
    !!                              none polar
    !! @param[in]  conserve_volume  Renormalize radii to fixed volume
    !! @param[in]  apply_com        Apply the centre-of-mass correction
    !! @param[out] status           SHAPE_VALID on success, else the rejecting code
    subroutine cache_init_s(cache, n_params, thetas, conserve_volume, apply_com, status)
        type(cache_t),      intent(out) :: cache
        integer(kind = ik), intent(in)  :: n_params
        real(kind = rk),    intent(in)  :: thetas(:)
        logical,            intent(in)  :: conserve_volume
        logical,            intent(in)  :: apply_com
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: masks(N_INTERMEDIATES)
        integer(kind = ik) :: n

        ! Engine first: it owns the n_params contract (< 1 or > 8 rejected here).
        call engine_masks_s(n_params, masks)
        call shape_engine_init_s(cache%engine, n_params, masks, status)
        if (status /= SHAPE_VALID) return

        ! Tables next; tables_init_s also validates the theta set.
        allocate(cache%tp)
        call tables_init_s(cache%tp, n_params, thetas, status)
        if (status /= SHAPE_VALID) then
            deallocate(cache%tp)
            nullify(cache%tp)
            return
        end if
        cache%owns_tables = .true.

        n = tables_n_thetas_f(cache%tp)
        allocate(cache%radii(n), cache%dr_dthetas(n))
        cache%radii(:)      = 0.0_rk
        cache%dr_dthetas(:) = 0.0_rk

        cache%n_params        = n_params
        cache%conserve_volume = conserve_volume
        cache%apply_com       = apply_com
        cache%is_initialized  = .true.
    end subroutine cache_init_s

    !> Initialize a cache over caller-owned shared tables.
    !!
    !! The cache stores a pointer to `tables` and never frees it. **The caller
    !! MUST declare `tables` with the `target` attribute**; without it the
    !! association is undefined once this routine returns. `tables` must also
    !! outlive the cache and must not be freed while the cache is in use.
    !!
    !! @param[in]  tables           Initialized shared tables; declared `target`
    !!                              by the caller
    !! @param[out] cache            Ready on success; default (uninitialized)
    !!                              state on any failure
    !! @param[in]  n_params         1 <= n_params <= min(SHAPE_CACHE_MAX_PARAMS,
    !!                              tables max_l)
    !! @param[in]  conserve_volume  Renormalize radii to fixed volume
    !! @param[in]  apply_com        Apply the centre-of-mass correction
    !! @param[out] status           SHAPE_VALID on success, else the rejecting code
    subroutine cache_init_shared_s(cache, tables, n_params, conserve_volume, apply_com, status)
        type(cache_t),      intent(out)        :: cache
        type(tables_t),     intent(in), target :: tables
        integer(kind = ik), intent(in)         :: n_params
        logical,            intent(in)         :: conserve_volume
        logical,            intent(in)         :: apply_com
        integer(kind = ik), intent(out)        :: status

        integer(kind = ik) :: masks(N_INTERMEDIATES)
        integer(kind = ik) :: n

        status = SHAPE_VALID
        if (.not. tables%is_initialized) then
            status = SHAPE_ERROR_TABLES_NOT_INITIALIZED
            return
        end if

        call engine_masks_s(n_params, masks)
        call shape_engine_init_s(cache%engine, n_params, masks, status)
        if (status /= SHAPE_VALID) return

        ! The shared tables cap n_params independently of the engine cap.
        if (n_params > tables_max_l_f(tables)) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if

        cache%tp => tables
        cache%owns_tables = .false.

        n = tables_n_thetas_f(tables)
        allocate(cache%radii(n), cache%dr_dthetas(n))
        cache%radii(:)      = 0.0_rk
        cache%dr_dthetas(:) = 0.0_rk

        cache%n_params        = n_params
        cache%conserve_volume = conserve_volume
        cache%apply_com       = apply_com
        cache%is_initialized  = .true.
    end subroutine cache_init_shared_s

    !> Release the cache. Infallible, and safe on an uninitialized cache.
    !!
    !! Frees the tables only when the cache owns them (private mode); shared
    !! tables are left untouched for their owner. `intent(inout)`, not
    !! `intent(out)`: the ownership flag and the pointer must still be readable
    !! on entry.
    !!
    !! @param[inout] cache  Reset to the default (uninitialized) state
    subroutine cache_free_s(cache)
        type(cache_t), intent(inout) :: cache
        if (cache%owns_tables .and. associated(cache%tp)) then
            call tables_free_s(cache%tp)
            deallocate(cache%tp)
        end if
        nullify(cache%tp)
        cache = cache_t()   ! deallocates radii/dr_dthetas, restores defaults
    end subroutine cache_free_s

    !> Parameter count fixed at init; 0 when uninitialized.
    pure function cache_n_params_f(cache) result(n)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: n
        n = 0_ik
        if (cache%is_initialized) n = cache%n_params
    end function cache_n_params_f

    !> Number of primary theta nodes behind this cache; 0 when uninitialized.
    pure function cache_n_thetas_f(cache) result(n)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: n
        n = 0_ik
        if (cache%is_initialized) then
            if (associated(cache%tp)) n = tables_n_thetas_f(cache%tp)
        end if
    end function cache_n_thetas_f

    !> .true. only after a successful init and before `cache_free_s`.
    pure function cache_is_initialized_f(cache) result(ok)
        type(cache_t), intent(in) :: cache
        logical :: ok
        ok = cache%is_initialized
    end function cache_is_initialized_f

end module beta_parameterization_mod
