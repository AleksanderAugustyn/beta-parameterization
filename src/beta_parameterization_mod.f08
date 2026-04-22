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
    use beta_parameterization_workers_mod, only: &
            precompute_legendre_table_s, &
            eval_polar_radii_s, &
            eval_radius_grid_s, &
            find_min_radius_s, &
            iterate_com_correction_s

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public types
    !---------------------------------------------------------------------------
    public :: cache_t

    !---------------------------------------------------------------------------
    ! Public standalone procedures (one-off, build a temporary cache internally)
    !---------------------------------------------------------------------------
    public :: compute_radius_grid_standalone_s
    public :: compute_radius_grid_standalone_with_com_shift_s

    !---------------------------------------------------------------------------
    ! Public limits
    !---------------------------------------------------------------------------
    integer(kind = ik), parameter, public :: MAX_BETA_PARAMS_LIMIT = 64_ik

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
        procedure, pass(self) :: max_beta_params_get
        procedure, pass(self) :: n_grid_get
        procedure, pass(self) :: is_initialized_get
        final :: cache_final
    end type cache_t

contains

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
    !! @param[in]  n_grid          Number of θ grid points (≥ 2)
    !! @param[out] error_code      LEGENDRE_VALID on success, else error code
    !! @param[out] message         Human-readable error message (empty on success)
    subroutine cache_init_s(self, max_beta_params, n_grid, error_code, message)

        class(cache_t),     intent(out) :: self
        integer(kind = ik), intent(in)  :: max_beta_params
        integer(kind = ik), intent(in)  :: n_grid
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

        if (n_grid < 2_ik) then
            error_code = LEGENDRE_ERROR_INVALID_MAX_PARAMS  ! reuse — bad init parameter
            write(message, '(A,I0)') 'n_grid must be >= 2, got ', n_grid
            return
        end if

        self%max_beta_params = max_beta_params
        self%n_grid          = n_grid
        self%n_quad          = 512_ik

        ! Normalization constants C_λ for λ = 1..max_beta_params
        allocate(self%norm_constants(max_beta_params))
        call compute_spherical_harmonics_normalization_constants_s( &
                self%norm_constants, max_beta_params)

        ! Gauss–Legendre nodes and weights (fixed n_quad = 512)
        allocate(self%gl_nodes(self%n_quad), self%gl_weights(self%n_quad))
        call compute_gauss_legendre_quadrature_s( &
                self%n_quad, self%gl_nodes, self%gl_weights)

        ! Legendre P_k at the θ grid: x_i = cos(θ_i), θ_i = (i-1) * π/(n_grid-1)
        allocate(self%legendre_theta_grid(n_grid, max_beta_params + 1_ik))
        allocate(theta_x(n_grid))
        h = PI_C / real(n_grid - 1_ik, rk)
        do i = 1_ik, n_grid
            theta_x(i) = cos(real(i - 1_ik, rk) * h)
        end do
        call precompute_legendre_table_s(theta_x, max_beta_params, self%legendre_theta_grid)
        deallocate(theta_x)

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

        real(kind = rk) :: beta_local(self%max_beta_params)
        real(kind = rk) :: beta_con(self%max_beta_params)
        real(kind = rk) :: r_north, r_south, r_min, h, theta
        integer(kind = ik) :: n_params, i_min, n_iter
        logical :: converged

        error_code       = LEGENDRE_VALID
        message          = ''
        radii(:)         = 0.0_rk
        corrected_beta10 = 0.0_rk

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
        beta_local(:)          = 0.0_rk
        beta_local(1:n_params) = params(1:n_params)
        beta_con(:)            = beta_local(:) * self%norm_constants(:)

        ! Step 4: COM iteration (modifies beta_local(1) and beta_con(1))
        call iterate_com_correction_s( &
                beta_local, beta_con, self%norm_constants, &
                self%gl_nodes, self%gl_weights, self%legendre_gl, &
                converged, n_iter)

        corrected_beta10 = beta_local(1)

        if (.not. converged) then
            error_code = LEGENDRE_ERROR_COM_NOT_CONVERGED
            write(message, '(A,I0,A)') &
                    'COM correction failed to converge after ', n_iter, ' iterations'
            return
        end if

        ! Step 5: polar pre-check (post-COM)
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
