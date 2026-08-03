!> Public Fortran API for the beta parameterization library.
!!
!! Defines `cache_t` — an immutable-after-init container for all precomputed
!! data (normalization constants, Legendre polynomial tables, Gauss–Legendre
!! nodes/weights). Each `cache_t` is safe to share across threads after init;
!! all per-call workspace is stack-allocated automatic-array storage.
!!
!! ## Typical use
!!
!! ```fortran
!! type(cache_t) :: cache
!! integer(ik)   :: error_code
!! character(len = 256) :: message
!! real(rk)      :: radii(80)
!!
!! call cache%init(max_beta_params = 4_ik, n_grid = 80_ik, &
!!         error_code = error_code, message = message)
!! ! ... loop over many shapes:
!! call cache%compute_radius_grid(params, radii, error_code, message)
!! ```
!!
!! ## Validation order (hot path)
!!
!! 1. Param-count check (size(params) ≤ max_beta_params).
!! 2. Buffer-size check (size(radii) == n_grid).
!! 3. Pad params; compute beta_con.
!! 4. (with-shift only) COM iteration; non-convergence is a validation failure.
!! 5. Polar pre-check on the **final** beta_con (post-COM if applicable).
!! 6. Grid evaluation.
!! 7. Defense-in-depth check via min radius — distinguishes north/south/interior.
module beta_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: &
            compute_spherical_harmonics_normalization_constants_s, &
            compute_gauss_legendre_quadrature_s
    use shape_core_mod, only: &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
            SHAPE_CACHE_MAX_PARAMS, SHAPE_STANDALONE_MAX_PARAMS
    use beta_parameterization_workers_mod, only: &
            precompute_legendre_table_s, &
            precompute_legendre_derivative_table_s, &
            eval_polar_radii_s, &
            eval_radius_grid_s, &
            eval_radius_derivative_s, &
            find_min_radius_s, &
            iterate_com_correction_s

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public types
    !---------------------------------------------------------------------------
    public :: cache_t
    public :: tables_t

    !---------------------------------------------------------------------------
    ! Public standalone procedures (one-off, build a temporary cache internally)
    !---------------------------------------------------------------------------
    public :: compute_radius_grid_standalone_s
    public :: compute_radius_grid_standalone_with_com_shift_s

    !---------------------------------------------------------------------------
    ! Tables lifecycle
    !---------------------------------------------------------------------------
    public :: tables_init_s, tables_free_s, tables_max_l_f, tables_n_thetas_f

    !---------------------------------------------------------------------------
    ! Public limits
    !---------------------------------------------------------------------------
    integer(kind = ik), parameter, public :: MAX_BETA_PARAMS_LIMIT = 64_ik

    !> Fixed Gauss-Legendre quadrature order for volume/COM integrals.
    integer(kind = ik), parameter :: N_QUAD = 512_ik

    !---------------------------------------------------------------------------
    ! Validation status codes
    !---------------------------------------------------------------------------
    integer(kind = ik), parameter, public :: LEGENDRE_VALID                     = 0_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_NORTH_POLE          = 1_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_SOUTH_POLE          = 2_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_EMPTY_PARAMS        = 3_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_INTERIOR_NEGATIVE   = 4_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_INVALID_MAX_PARAMS  = 5_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_TOO_MANY_PARAMS     = 6_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_COM_NOT_CONVERGED   = 7_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_INVALID_BUFFER_SIZE = 8_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_POLE_NODE          = 9_ik
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_NO_UNIFORM_GRID    = 10_ik

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

    !---------------------------------------------------------------------------
    ! Cache type
    !---------------------------------------------------------------------------
    !! Invariants:
    !!   - After successful init: is_initialized == .true. and all allocatable
    !!     components are allocated to their declared sizes.
    !!   - After destroy: is_initialized == .false. and all components deallocated.
    !!   - n_quad is fixed at 512 (not configurable).
    !!   - All hot-path methods declare `intent(in) :: self` — the cache is
    !!     immutable after init and safe to share across threads.
    type, public :: cache_t
        private
        integer(kind = ik) :: max_beta_params = 0_ik
        integer(kind = ik) :: n_grid          = 0_ik
        integer(kind = ik) :: n_quad          = 512_ik

        real(kind = rk), allocatable :: norm_constants(:)         ! (max_beta_params)
        real(kind = rk), allocatable :: legendre_theta_grid(:, :) ! (n_grid, max_beta_params + 1)
        real(kind = rk), allocatable :: gl_nodes(:)               ! (n_quad)
        real(kind = rk), allocatable :: gl_weights(:)             ! (n_quad)
        real(kind = rk), allocatable :: legendre_gl(:, :)         ! (n_quad, max_beta_params + 1)

        logical :: is_initialized = .false.
    contains
        procedure, pass(self) :: init    => cache_init_s
        procedure, pass(self) :: destroy => cache_destroy_s
        procedure, pass(self) :: compute_radius_grid &
                              => cache_compute_radius_grid_s
        procedure, pass(self) :: compute_radius_grid_with_com_shift &
                              => cache_compute_radius_grid_with_com_shift_s
        procedure, pass(self) :: build_node_set => cache_build_node_set_s
        procedure, pass(self) :: resolve_shape  => cache_resolve_shape_s
        procedure, pass(self) :: compute_radius_and_derivative &
                              => cache_compute_radius_and_derivative_s
        procedure, pass(self) :: max_beta_params_get
        procedure, pass(self) :: n_grid_get
        procedure, pass(self) :: is_initialized_get
        final :: cache_final
    end type cache_t

    !> A caller-owned set of theta nodes with precomputed Legendre P_k and P_k'
    !! tables sized to one cache's max_beta_params. Built once (startup), then
    !! reused across shapes. Carries no shape data; immutable after build and
    !! safe to share across threads read-only. Deliberately generic: no
    !! dense/folding/coulomb vocabulary — the caller owns what a set means.
    type, public :: node_set_t
        private
        logical            :: is_built = .false.
        integer(kind = ik) :: n_nodes  = 0_ik
        real(kind = rk), allocatable :: thetas(:)
        real(kind = rk), allocatable :: sin_thetas(:)
        real(kind = rk), allocatable :: legendre_table(:, :)        ! P_k(cos theta_i)
        real(kind = rk), allocatable :: legendre_deriv_table(:, :)  ! P_k'(cos theta_i)
    contains
        procedure, pass(self) :: is_built_get => node_set_is_built_get
        procedure, pass(self) :: n_nodes_get  => node_set_n_nodes_get
        procedure, pass(self) :: destroy      => node_set_destroy_s
    end type node_set_t

    !> Shared immutable level: everything that depends only on `max_l` and the
    !! primary theta set. Built once, then shared read-only by every consumer
    !! (node sets, caches, standalone entry points) — no shape data inside.
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
    ! LIFECYCLE
    !===========================================================================

    !> Initialize the cache: allocate components, fill precomputed tables.
    !!
    !! Idempotent — calling on an already-initialized cache destroys + re-inits
    !! (the `intent(out) :: self` triggers finalization of the old contents).
    !!
    !! @param[out] self            The cache (intent(out) auto-finalizes any prior state)
    !! @param[in]  max_beta_params 1 ≤ N ≤ MAX_BETA_PARAMS_LIMIT
    !! @param[in]  n_grid          Optional number of θ grid points (≥ 2) for
    !!                             the legacy uniform-grid entry points. Absent
    !!                             ⇒ node-set-only cache: the uniform table is
    !!                             not allocated and compute_radius_grid[_with_
    !!                             com_shift] return LEGENDRE_ERROR_NO_UNIFORM_GRID.
    !! @param[out] error_code      LEGENDRE_VALID on success, else error code
    !! @param[out] message         Human-readable error message (empty on success)
    subroutine cache_init_s(self, max_beta_params, n_grid, error_code, message)

        class(cache_t),     intent(out) :: self
        integer(kind = ik), intent(in)  :: max_beta_params
        integer(kind = ik), intent(in), optional :: n_grid
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message

        integer(kind = ik) :: i
        real(kind = rk)    :: h
        real(kind = rk), allocatable :: theta_x(:)

        error_code = LEGENDRE_VALID
        message    = ''

        if (max_beta_params < 1_ik .or. max_beta_params > MAX_BETA_PARAMS_LIMIT) then
            error_code = LEGENDRE_ERROR_INVALID_MAX_PARAMS
            write(message, '(A,I0,A,I0,A)') &
                    'max_beta_params must be in [1, ', MAX_BETA_PARAMS_LIMIT, &
                    '], got ', max_beta_params, '.'
            return
        end if

        if (present(n_grid)) then
            if (n_grid < 2_ik) then
                error_code = LEGENDRE_ERROR_INVALID_MAX_PARAMS  ! reuse — bad init parameter
                write(message, '(A,I0)') 'n_grid must be >= 2, got ', n_grid
                return
            end if
            self%n_grid = n_grid
        else
            self%n_grid = 0_ik    ! node-set-only cache: no uniform table
        end if

        self%max_beta_params = max_beta_params
        self%n_quad          = 512_ik

        ! Normalization constants C_λ for λ = 1..max_beta_params
        allocate(self%norm_constants(max_beta_params))
        call compute_spherical_harmonics_normalization_constants_s( &
                self%norm_constants, max_beta_params)

        ! Gauss–Legendre nodes and weights (fixed n_quad = 512)
        allocate(self%gl_nodes(self%n_quad), self%gl_weights(self%n_quad))
        call compute_gauss_legendre_quadrature_s( &
                self%n_quad, self%gl_nodes, self%gl_weights)

        ! Legendre P_k at the θ grid — only when a uniform grid was requested.
        if (present(n_grid)) then
            allocate(self%legendre_theta_grid(n_grid, max_beta_params + 1_ik))
            allocate(theta_x(n_grid))
            h = PI_C / real(n_grid - 1_ik, rk)
            do i = 1_ik, n_grid
                theta_x(i) = cos(real(i - 1_ik, rk) * h)
            end do
            call precompute_legendre_table_s(theta_x, max_beta_params, self%legendre_theta_grid)
            deallocate(theta_x)
        end if

        ! Legendre P_k at the GL nodes
        allocate(self%legendre_gl(self%n_quad, max_beta_params + 1_ik))
        call precompute_legendre_table_s(self%gl_nodes, max_beta_params, self%legendre_gl)

        self%is_initialized = .true.

    end subroutine cache_init_s

    !> Deallocate all components and reset flags. Safe on uninitialized cache.
    pure subroutine cache_destroy_s(self)
        class(cache_t), intent(inout) :: self
        if (allocated(self%norm_constants))      deallocate(self%norm_constants)
        if (allocated(self%legendre_theta_grid)) deallocate(self%legendre_theta_grid)
        if (allocated(self%gl_nodes))            deallocate(self%gl_nodes)
        if (allocated(self%gl_weights))          deallocate(self%gl_weights)
        if (allocated(self%legendre_gl))         deallocate(self%legendre_gl)
        self%max_beta_params = 0_ik
        self%n_grid          = 0_ik
        self%is_initialized  = .false.
    end subroutine cache_destroy_s

    !> Final binding — automatic cleanup on out-of-scope.
    impure elemental subroutine cache_final(self)
        type(cache_t), intent(inout) :: self
        call cache_destroy_s(self)
    end subroutine cache_final

    !===========================================================================
    ! ACCESSORS
    !===========================================================================

    pure function max_beta_params_get(self) result(n)
        class(cache_t), intent(in) :: self
        integer(kind = ik) :: n
        n = self%max_beta_params
    end function

    pure function n_grid_get(self) result(n)
        class(cache_t), intent(in) :: self
        integer(kind = ik) :: n
        n = self%n_grid
    end function

    pure function is_initialized_get(self) result(b)
        class(cache_t), intent(in) :: self
        logical :: b
        b = self%is_initialized
    end function

    !===========================================================================
    ! NODE SETS
    !===========================================================================

    pure function node_set_is_built_get(self) result(b)
        class(node_set_t), intent(in) :: self
        logical :: b
        b = self%is_built
    end function node_set_is_built_get

    pure function node_set_n_nodes_get(self) result(n)
        class(node_set_t), intent(in) :: self
        integer(kind = ik) :: n
        n = self%n_nodes
    end function node_set_n_nodes_get

    !> Deallocate all components and reset flags. Safe on an unbuilt set.
    pure subroutine node_set_destroy_s(self)
        class(node_set_t), intent(inout) :: self
        if (allocated(self%thetas))               deallocate(self%thetas)
        if (allocated(self%sin_thetas))           deallocate(self%sin_thetas)
        if (allocated(self%legendre_table))       deallocate(self%legendre_table)
        if (allocated(self%legendre_deriv_table)) deallocate(self%legendre_deriv_table)
        self%is_built = .false.
        self%n_nodes  = 0_ik
    end subroutine node_set_destroy_s

    !> Precompute P_k and P_k' tables at caller-supplied theta nodes.
    !!
    !! Pole nodes (1 - cos(theta)**2 == 0 in double precision) are rejected with
    !! LEGENDRE_ERROR_POLE_NODE: the derivative recursion divides by 1 - x**2.
    !! Gauss-Legendre nodes never land there; pole radii come analytically from
    !! resolve_shape instead.
    !!
    !! @param[in]  thetas      Node angles (radians); any order, need not be uniform
    !! @param[out] node_set    Filled tables (intent(out) resets any prior build)
    !! @param[out] error_code  LEGENDRE_VALID on success
    !! @param[out] message     Empty on success
    subroutine cache_build_node_set_s(self, thetas, node_set, error_code, message)

        class(cache_t),     intent(in)  :: self
        real(kind = rk),    intent(in)  :: thetas(:)
        type(node_set_t),   intent(out) :: node_set
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message

        integer(kind = ik) :: n, i
        real(kind = rk), allocatable :: x(:)

        error_code = LEGENDRE_VALID
        message    = ''

        if (.not. self%is_initialized) then
            error_code = LEGENDRE_ERROR_INVALID_MAX_PARAMS  ! reuse — API misuse, like cache_init_s
            message    = 'build_node_set: cache not initialized'
            return
        end if

        n = size(thetas, kind = ik)
        if (n < 1_ik) then
            error_code = LEGENDRE_ERROR_INVALID_BUFFER_SIZE
            message    = 'build_node_set: thetas is empty'
            return
        end if

        allocate(x(n))
        x = cos(thetas)
        do i = 1_ik, n
            if (1.0_rk - x(i)**2 <= 0.0_rk) then
                error_code = LEGENDRE_ERROR_POLE_NODE
                write(message, '(A,ES12.4,A)') &
                        'build_node_set: node at theta = ', thetas(i), &
                        ' is a pole; use resolve_shape polar radii instead'
                return
            end if
        end do

        node_set%n_nodes = n
        allocate(node_set%thetas(n), node_set%sin_thetas(n))
        node_set%thetas     = thetas
        node_set%sin_thetas = sin(thetas)
        allocate(node_set%legendre_table(n, self%max_beta_params + 1_ik))
        allocate(node_set%legendre_deriv_table(n, self%max_beta_params + 1_ik))
        call precompute_legendre_table_s(x, self%max_beta_params, node_set%legendre_table)
        call precompute_legendre_derivative_table_s(x, self%max_beta_params, &
                node_set%legendre_table, node_set%legendre_deriv_table)
        node_set%is_built = .true.

    end subroutine cache_build_node_set_s

    !> Resolve a shape once: pad + normalize params, run the COM iteration,
    !! polar pre-check. Node-set-independent — feed the resulting beta_con to
    !! compute_radius_and_derivative for any number of node sets.
    !!
    !! On any failure all outputs are zero-filled.
    !!
    !! @param[in]  params            Deformation parameters (beta10 is the COM dipole)
    !! @param[out] beta_con          Normalized, COM-corrected coefficients;
    !!                               size must equal cache max_beta_params
    !! @param[out] corrected_beta10  Post-shift beta10
    !! @param[out] r_north           Analytic R(0) from eval_polar_radii_s
    !! @param[out] r_south           Analytic R(pi) from eval_polar_radii_s
    !! @param[out] error_code        LEGENDRE_VALID on success
    !! @param[out] message           Empty on success
    !! @param[in]  apply_com_correction  Optional (default .true.); .false.
    !!                               skips the COM iteration (corrected_beta10
    !!                               = input beta10), matching
    !!                               compute_radius_grid's no-COM semantics.
    subroutine cache_resolve_shape_s(self, params, beta_con, corrected_beta10, &
            r_north, r_south, error_code, message, apply_com_correction)

        class(cache_t),     intent(in)  :: self
        real(kind = rk),    intent(in)  :: params(:)
        real(kind = rk),    intent(out) :: beta_con(:)
        real(kind = rk),    intent(out) :: corrected_beta10
        real(kind = rk),    intent(out) :: r_north
        real(kind = rk),    intent(out) :: r_south
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message
        logical,            intent(in), optional :: apply_com_correction

        real(kind = rk)    :: beta_local(self%max_beta_params)
        integer(kind = ik) :: n_params, n_iter
        logical            :: converged, do_com

        error_code       = LEGENDRE_VALID
        message          = ''
        beta_con(:)      = 0.0_rk
        corrected_beta10 = 0.0_rk
        r_north          = 0.0_rk
        r_south          = 0.0_rk

        n_params = size(params, kind = ik)
        if (n_params < 1_ik) then
            error_code = LEGENDRE_ERROR_EMPTY_PARAMS
            message    = 'params is empty'
            return
        end if
        if (n_params > self%max_beta_params) then
            error_code = LEGENDRE_ERROR_TOO_MANY_PARAMS
            write(message, '(A,I0,A,I0)') &
                    'params has ', n_params, &
                    ' entries but cache built for max_beta_params = ', self%max_beta_params
            return
        end if
        if (size(beta_con, kind = ik) /= self%max_beta_params) then
            error_code = LEGENDRE_ERROR_INVALID_BUFFER_SIZE
            write(message, '(A,I0,A,I0)') &
                    'beta_con buffer size ', size(beta_con), &
                    ' does not match cache max_beta_params ', self%max_beta_params
            return
        end if

        beta_local(:)          = 0.0_rk
        beta_local(1:n_params) = params(1:n_params)
        beta_con(:)            = beta_local(:) * self%norm_constants(:)

        do_com = .true.
        if (present(apply_com_correction)) do_com = apply_com_correction

        if (do_com) then
            call iterate_com_correction_s( &
                    beta_local, beta_con, self%norm_constants, &
                    self%gl_nodes, self%gl_weights, self%legendre_gl, &
                    converged, n_iter)

            corrected_beta10 = beta_local(1)

            if (.not. converged) then
                error_code = LEGENDRE_ERROR_COM_NOT_CONVERGED
                write(message, '(A,I0,A)') &
                        'COM correction failed to converge after ', n_iter, ' iterations'
                beta_con(:)      = 0.0_rk
                corrected_beta10 = 0.0_rk
                return
            end if
        else
            corrected_beta10 = beta_local(1)
        end if

        call eval_polar_radii_s(beta_con, r_north, r_south)
        if (r_north <= R_MIN_THRESHOLD) then
            error_code = LEGENDRE_ERROR_NORTH_POLE
            write(message, '(A,ES12.4)') 'North pole radius not positive: R(0) = ', r_north
            beta_con(:)      = 0.0_rk
            corrected_beta10 = 0.0_rk
            r_north          = 0.0_rk
            r_south          = 0.0_rk
            return
        end if
        if (r_south <= R_MIN_THRESHOLD) then
            error_code = LEGENDRE_ERROR_SOUTH_POLE
            write(message, '(A,ES12.4)') 'South pole radius not positive: R(pi) = ', r_south
            beta_con(:)      = 0.0_rk
            corrected_beta10 = 0.0_rk
            r_north          = 0.0_rk
            r_south          = 0.0_rk
            return
        end if

    end subroutine cache_resolve_shape_s

    !> Evaluate R and dR/dtheta at a node set, for a shape already resolved by
    !! resolve_shape. Dot products only — no iteration, no allocation.
    !!
    !! Defense-in-depth: the interior-positivity scan always reports
    !! LEGENDRE_ERROR_INTERIOR_NEGATIVE (node sets exclude the poles, so the
    !! uniform-grid pole-index mapping does not apply). Zero-fills on failure.
    !!
    !! @param[in]  beta_con    Resolved coefficients from resolve_shape;
    !!                         size must equal cache max_beta_params
    !! @param[in]  node_set    Built by this cache's build_node_set
    !! @param[out] radii       R(theta_i); size must equal node set n_nodes
    !! @param[out] dr_dthetas  dR/dtheta(theta_i); size must equal node set n_nodes
    !! @param[out] error_code  LEGENDRE_VALID on success
    !! @param[out] message     Empty on success
    subroutine cache_compute_radius_and_derivative_s(self, beta_con, node_set, &
            radii, dr_dthetas, error_code, message)

        class(cache_t),     intent(in)  :: self
        real(kind = rk),    intent(in)  :: beta_con(:)
        type(node_set_t),   intent(in)  :: node_set
        real(kind = rk),    intent(out) :: radii(:)
        real(kind = rk),    intent(out) :: dr_dthetas(:)
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message

        real(kind = rk)    :: r_min
        integer(kind = ik) :: i_min

        error_code    = LEGENDRE_VALID
        message       = ''
        radii(:)      = 0.0_rk
        dr_dthetas(:) = 0.0_rk

        if (.not. node_set%is_built) then
            error_code = LEGENDRE_ERROR_INVALID_BUFFER_SIZE
            message    = 'compute_radius_and_derivative: node set not built'
            return
        end if
        if (size(node_set%legendre_table, 2, kind = ik) /= self%max_beta_params + 1_ik) then
            error_code = LEGENDRE_ERROR_INVALID_BUFFER_SIZE
            message    = 'compute_radius_and_derivative: node set built by a different cache'
            return
        end if
        if (size(beta_con, kind = ik) /= self%max_beta_params) then
            error_code = LEGENDRE_ERROR_INVALID_BUFFER_SIZE
            write(message, '(A,I0,A,I0)') &
                    'beta_con size ', size(beta_con), &
                    ' does not match cache max_beta_params ', self%max_beta_params
            return
        end if
        if (size(radii, kind = ik) /= node_set%n_nodes .or. &
                size(dr_dthetas, kind = ik) /= node_set%n_nodes) then
            error_code = LEGENDRE_ERROR_INVALID_BUFFER_SIZE
            write(message, '(A,I0)') &
                    'radii/dr_dthetas buffers must match node set n_nodes = ', node_set%n_nodes
            return
        end if

        call eval_radius_grid_s(beta_con, node_set%legendre_table, radii)
        call eval_radius_derivative_s(beta_con, node_set%legendre_deriv_table, &
                node_set%sin_thetas, dr_dthetas)

        call find_min_radius_s(radii, r_min, i_min)
        if (r_min <= R_MIN_THRESHOLD) then
            error_code = LEGENDRE_ERROR_INTERIOR_NEGATIVE
            write(message, '(A,F8.4,A,ES12.4)') &
                    'Radius not positive at theta = ', node_set%thetas(i_min), ', R = ', r_min
            radii(:)      = 0.0_rk
            dr_dthetas(:) = 0.0_rk
            return
        end if

    end subroutine cache_compute_radius_and_derivative_s

    !===========================================================================
    ! HOT PATH
    !===========================================================================

    !> Compute R(θ) on the cache's θ grid, no COM shift.
    !!
    !! params shorter than max_beta_params is silently zero-padded (treat as
    !! a lower-order shape). Longer than max_beta_params returns
    !! LEGENDRE_ERROR_TOO_MANY_PARAMS.
    !!
    !! @param[in]  params      Deformation parameters β
    !! @param[out] radii       R(θ) values; size must equal cache n_grid
    !! @param[out] error_code  LEGENDRE_VALID on success
    !! @param[out] message     Empty on success
    subroutine cache_compute_radius_grid_s(self, params, radii, error_code, message)

        class(cache_t),     intent(in)  :: self
        real(kind = rk),    intent(in)  :: params(:)
        real(kind = rk),    intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message

        real(kind = rk) :: beta_local(self%max_beta_params)
        real(kind = rk) :: beta_con(self%max_beta_params)
        real(kind = rk) :: r_north, r_south, r_min, h, theta
        integer(kind = ik) :: n_params, i_min

        error_code = LEGENDRE_VALID
        message    = ''
        radii(:)   = 0.0_rk

        if (self%n_grid == 0_ik) then
            error_code = LEGENDRE_ERROR_NO_UNIFORM_GRID
            message    = 'cache built without a uniform grid (init n_grid absent); use the node-set API'
            return
        end if

        ! Step 1: param-count check
        n_params = size(params, kind = ik)
        if (n_params < 1_ik) then
            error_code = LEGENDRE_ERROR_EMPTY_PARAMS
            message    = 'params is empty'
            return
        end if
        if (n_params > self%max_beta_params) then
            error_code = LEGENDRE_ERROR_TOO_MANY_PARAMS
            write(message, '(A,I0,A,I0)') &
                    'params has ', n_params, &
                    ' entries but cache built for max_beta_params = ', self%max_beta_params
            return
        end if

        ! Step 2: buffer-size check
        if (size(radii, kind = ik) /= self%n_grid) then
            error_code = LEGENDRE_ERROR_INVALID_BUFFER_SIZE
            write(message, '(A,I0,A,I0)') &
                    'radii buffer size ', size(radii), &
                    ' does not match cache n_grid ', self%n_grid
            return
        end if

        ! Step 3: pad and form beta_con
        beta_local(:)            = 0.0_rk
        beta_local(1:n_params)   = params(1:n_params)
        beta_con(:)              = beta_local(:) * self%norm_constants(:)

        ! Step 4: skipped (no COM shift in this variant)

        ! Step 5: polar pre-check
        call eval_polar_radii_s(beta_con, r_north, r_south)
        if (r_north <= R_MIN_THRESHOLD) then
            error_code = LEGENDRE_ERROR_NORTH_POLE
            write(message, '(A,ES12.4)') 'North pole radius not positive: R(0) = ', r_north
            return
        end if
        if (r_south <= R_MIN_THRESHOLD) then
            error_code = LEGENDRE_ERROR_SOUTH_POLE
            write(message, '(A,ES12.4)') 'South pole radius not positive: R(pi) = ', r_south
            return
        end if

        ! Step 6: grid evaluation
        call eval_radius_grid_s(beta_con, self%legendre_theta_grid, radii)

        ! Step 7: defense-in-depth
        call find_min_radius_s(radii, r_min, i_min)
        if (r_min <= R_MIN_THRESHOLD) then
            if (i_min == 1_ik) then
                error_code = LEGENDRE_ERROR_NORTH_POLE
                write(message, '(A,ES12.4)') &
                        'North pole radius not positive: R(0) = ', radii(1)
            else if (i_min == self%n_grid) then
                error_code = LEGENDRE_ERROR_SOUTH_POLE
                write(message, '(A,ES12.4)') &
                        'South pole radius not positive: R(pi) = ', radii(self%n_grid)
            else
                error_code = LEGENDRE_ERROR_INTERIOR_NEGATIVE
                h     = PI_C / real(self%n_grid - 1_ik, rk)
                theta = real(i_min - 1_ik, rk) * h
                write(message, '(A,F8.4,A,ES12.4)') &
                        'Interior radius not positive at theta = ', theta, ', R = ', r_min
            end if
            return
        end if

    end subroutine cache_compute_radius_grid_s

    !> Compute R(θ) with COM shift applied to β₁₀ first.
    !!
    !! On non-convergence of the COM iteration, returns
    !! LEGENDRE_ERROR_COM_NOT_CONVERGED. v1 silently accepted non-convergence;
    !! v2 reports it.
    !!
    !! @param[in]  params            Deformation parameters β (β₁₀ is the COM dipole)
    !! @param[out] radii             R(θ) values; size must equal cache n_grid
    !! @param[out] corrected_beta10  Post-shift β₁₀
    !! @param[out] error_code        LEGENDRE_VALID on success
    !! @param[out] message           Empty on success
    subroutine cache_compute_radius_grid_with_com_shift_s( &
            self, params, radii, corrected_beta10, error_code, message)

        class(cache_t),     intent(in)  :: self
        real(kind = rk),    intent(in)  :: params(:)
        real(kind = rk),    intent(out) :: radii(:)
        real(kind = rk),    intent(out) :: corrected_beta10
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message

        real(kind = rk) :: beta_con(self%max_beta_params)
        real(kind = rk) :: r_north, r_south, r_min, h, theta
        integer(kind = ik) :: n_params, i_min

        error_code       = LEGENDRE_VALID
        message          = ''
        radii(:)         = 0.0_rk
        corrected_beta10 = 0.0_rk

        if (self%n_grid == 0_ik) then
            error_code = LEGENDRE_ERROR_NO_UNIFORM_GRID
            message    = 'cache built without a uniform grid (init n_grid absent); use the node-set API'
            return
        end if

        ! Step 1: param-count check
        n_params = size(params, kind = ik)
        if (n_params < 1_ik) then
            error_code = LEGENDRE_ERROR_EMPTY_PARAMS
            message    = 'params is empty'
            return
        end if
        if (n_params > self%max_beta_params) then
            error_code = LEGENDRE_ERROR_TOO_MANY_PARAMS
            write(message, '(A,I0,A,I0)') &
                    'params has ', n_params, &
                    ' entries but cache built for max_beta_params = ', self%max_beta_params
            return
        end if

        ! Step 2: buffer-size check
        if (size(radii, kind = ik) /= self%n_grid) then
            error_code = LEGENDRE_ERROR_INVALID_BUFFER_SIZE
            write(message, '(A,I0,A,I0)') &
                    'radii buffer size ', size(radii), &
                    ' does not match cache n_grid ', self%n_grid
            return
        end if

        ! Steps 3-5 (normalization, COM iteration, polar pre-check) live in
        ! resolve_shape; error codes and messages unchanged. Steps 1-2 above
        ! stay inline so error precedence is untouched.
        call self%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, &
                error_code, message)
        if (error_code /= LEGENDRE_VALID) return

        ! Step 6: grid evaluation
        call eval_radius_grid_s(beta_con, self%legendre_theta_grid, radii)

        ! Step 7: defense-in-depth
        call find_min_radius_s(radii, r_min, i_min)
        if (r_min <= R_MIN_THRESHOLD) then
            if (i_min == 1_ik) then
                error_code = LEGENDRE_ERROR_NORTH_POLE
                write(message, '(A,ES12.4)') &
                        'North pole radius not positive: R(0) = ', radii(1)
            else if (i_min == self%n_grid) then
                error_code = LEGENDRE_ERROR_SOUTH_POLE
                write(message, '(A,ES12.4)') &
                        'South pole radius not positive: R(pi) = ', radii(self%n_grid)
            else
                error_code = LEGENDRE_ERROR_INTERIOR_NEGATIVE
                h     = PI_C / real(self%n_grid - 1_ik, rk)
                theta = real(i_min - 1_ik, rk) * h
                write(message, '(A,F8.4,A,ES12.4)') &
                        'Interior radius not positive at theta = ', theta, ', R = ', r_min
            end if
            return
        end if

    end subroutine cache_compute_radius_grid_with_com_shift_s

    !===========================================================================
    ! STANDALONE
    !===========================================================================

    !> Standalone — builds a temporary cache, runs, lets `final` clean up.
    !!
    !! Inefficient if called repeatedly; use `cache_t` directly for batch runs.
    !!
    !! @param[in]  params           Deformation parameters β
    !! @param[in]  n_grid           Number of θ grid points
    !! @param[out] radii            R(θ) values; size must equal n_grid
    !! @param[out] error_code       LEGENDRE_VALID on success
    !! @param[out] message          Empty on success
    !! @param[in]  max_beta_params  Optional; defaults to size(params)
    subroutine compute_radius_grid_standalone_s( &
            params, n_grid, radii, error_code, message, max_beta_params)

        real(kind = rk),    intent(in)  :: params(:)
        integer(kind = ik), intent(in)  :: n_grid
        real(kind = rk),    intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message
        integer(kind = ik), intent(in), optional :: max_beta_params

        type(cache_t)      :: cache
        integer(kind = ik) :: mbp

        ! Outputs must be defined on every return path (v1 contract); the C
        ! wrappers copy them out unconditionally.
        radii(:) = 0.0_rk

        if (present(max_beta_params)) then
            mbp = max_beta_params
        else
            mbp = size(params, kind = ik)
        end if

        call cache%init(mbp, n_grid, error_code, message)
        if (error_code /= LEGENDRE_VALID) return

        call cache%compute_radius_grid(params, radii, error_code, message)
        ! cache%final fires here on scope exit.

    end subroutine compute_radius_grid_standalone_s

    subroutine compute_radius_grid_standalone_with_com_shift_s( &
            params, n_grid, radii, corrected_beta10, error_code, message, max_beta_params)

        real(kind = rk),    intent(in)  :: params(:)
        integer(kind = ik), intent(in)  :: n_grid
        real(kind = rk),    intent(out) :: radii(:)
        real(kind = rk),    intent(out) :: corrected_beta10
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message
        integer(kind = ik), intent(in), optional :: max_beta_params

        type(cache_t)      :: cache
        integer(kind = ik) :: mbp

        ! Outputs must be defined on every return path (v1 contract); the C
        ! wrappers copy them out unconditionally.
        radii(:)         = 0.0_rk
        corrected_beta10 = 0.0_rk

        if (present(max_beta_params)) then
            mbp = max_beta_params
        else
            mbp = size(params, kind = ik)
        end if

        call cache%init(mbp, n_grid, error_code, message)
        if (error_code /= LEGENDRE_VALID) return

        call cache%compute_radius_grid_with_com_shift( &
                params, radii, corrected_beta10, error_code, message)

    end subroutine compute_radius_grid_standalone_with_com_shift_s

end module beta_parameterization_mod
