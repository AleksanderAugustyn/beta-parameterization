!> Module for Legendre polynomial nuclear shape parameterization.
!!
!! ## Architecture
!!
!! This module is responsible for:
!! 1. Computing the Legendre/β-parameter shape at grid points
!! 2. Validating shape parameters (polar radii positive, interior radii positive)
!! 3. Optionally applying center-of-mass correction to β₁₀
!! 4. Providing R(θ) values for the radius grid
!!
!! ## Defense in Depth: Interior Point Validation
!!
!! Beyond checking polar radii, this module now validates that ALL interior
!! points have positive radii. This catches pathological shapes where large
!! higher-order deformations could create self-intersecting or "pinched" shapes
!! even when the poles appear valid.
!!
!! ## Physical Model
!!
!! The nuclear radius is expressed as:
!! \[ R(\theta) = R_0 \left[ 1 + \sum_{\lambda=1}^{\lambda_{max}} \beta_\lambda C_\lambda P_\lambda(\cos\theta) \right] \]
!!
!! where:
!! - \( \beta_\lambda \) are the deformation parameters
!! - \( C_\lambda = \sqrt{(2\lambda+1)/(4\pi)} \) are normalization constants
!! - \( P_\lambda \) are Legendre polynomials
!!
!! ## Common Deformations
!!
!! | Parameter | Name | Physical Meaning |
!! |-----------|------|------------------|
!! | β₁ | Dipole | Center-of-mass shift |
!! | β₂ | Quadrupole | Prolate (β₂>0) or oblate (β₂<0) deformation |
!! | β₃ | Octupole | Pear-shaped (reflection asymmetric) nuclei |
!! | β₄ | Hexadecapole | Higher-order shape refinement |
!!
module beta_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public Interface - Main Entry Point
    !---------------------------------------------------------------------------

    !> Main entry point: validates shape and computes radius grid
    public :: compute_legendre_radius_grid_s

    !> Standalone wrapper: internalizes all precomputation
    public :: compute_legendre_radius_grid_standalone_s

    !---------------------------------------------------------------------------
    ! Public Interface - Components
    !---------------------------------------------------------------------------

    !> Direct radius computation (legacy interface)
    public :: compute_radius_legendre_s

    !> Center-of-mass correction
    public :: correct_beta10_for_center_of_mass_s

    !> Shape validation
    public :: validate_legendre_params_s

    !---------------------------------------------------------------------------
    ! Validation Error Codes
    !---------------------------------------------------------------------------

    !> Shape is valid
    integer(kind = ik), parameter, public :: LEGENDRE_VALID = 0_ik

    !> Polar radius at θ=0 is not positive
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_NORTH_POLE = 1_ik

    !> Polar radius at θ=π is not positive
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_SOUTH_POLE = 2_ik

    !> Empty parameter array
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_EMPTY_PARAMS = 3_ik

    !> Interior radius not positive (shape self-intersects or pinches)
    integer(kind = ik), parameter, public :: LEGENDRE_ERROR_INTERIOR_NEGATIVE = 4_ik

    !---------------------------------------------------------------------------
    ! Defense in Depth: Minimum Radius Threshold
    !---------------------------------------------------------------------------
    !! Minimum allowed radius at any point. A small positive value prevents
    !! near-singular shapes that could cause numerical issues downstream.
    real(kind = rk), parameter :: R_MIN_THRESHOLD = 1.0e-6_rk

contains

    !===========================================================================
    ! MAIN ENTRY POINT
    !===========================================================================

    !> Validates Legendre shape and computes radius grid for radius_grid_mod.
    !!
    !! This is the main entry point for Legendre/β-parameter shape processing.
    !! It:
    !! 1. Optionally applies center-of-mass correction to β₁₀
    !! 2. Validates that polar radii are positive
    !! 3. Computes R(θ) values at grid points
    !! 4. Validates that ALL interior radii are positive (defense in depth)
    !!
    !! @param[in]  params              Shape parameters (β values)
    !! @param[in]  n_grid              Number of grid points for R(θ) output
    !! @param[out] radii               R(θ) values at grid points (size n_grid)
    !! @param[out] corrected_beta10    The (possibly corrected) β₁₀ value
    !! @param[out] is_valid            True if shape is valid
    !! @param[out] message             Error message if invalid
    !! @param[in]  apply_com_correction Optional: apply COM correction (default .false.)
    subroutine compute_legendre_radius_grid_s(&
            params, n_grid, norm_constants, legendre_theta_grid, &
            gl_nodes, gl_weights, legendre_gl, &
            radii, corrected_beta10, is_valid, message, apply_com_correction)

        implicit none

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_grid
        real(kind = rk), intent(in) :: norm_constants(:)
        real(kind = rk), intent(in) :: legendre_theta_grid(:, :)
        real(kind = rk), intent(in) :: gl_nodes(:)
        real(kind = rk), intent(in) :: gl_weights(:)
        real(kind = rk), intent(in) :: legendre_gl(:, :)
        real(kind = rk), intent(out) :: radii(n_grid)
        real(kind = rk), intent(out) :: corrected_beta10
        logical, intent(out) :: is_valid
        character(len = *), intent(out) :: message
        logical, intent(in), optional :: apply_com_correction

        ! Local variables
        integer(kind = ik) :: n_def_params, n_grid_precomp
        real(kind = rk), allocatable :: beta_local(:), beta_con(:)
        real(kind = rk) :: h, theta, r_raw, r_min
        logical :: do_com_correction
        integer(kind = ik) :: i, k, n_params, i_min
        integer(kind = ik) :: error_code

        n_def_params = size(norm_constants)
        n_grid_precomp = size(legendre_theta_grid, 1)
        allocate(beta_local(n_def_params), beta_con(n_def_params))

        ! Initialize outputs
        radii = 0.0_rk
        corrected_beta10 = 0.0_rk
        is_valid = .false.
        message = ''

        ! Determine if COM correction is requested
        do_com_correction = .false.
        if (present(apply_com_correction)) do_com_correction = apply_com_correction

        ! Copy parameters to local array (with zero-padding if needed)
        n_params = size(params, kind = ik)
        beta_local(:) = 0.0_rk
        beta_local(1:min(n_params, n_def_params)) = &
                params(1:min(n_params, n_def_params))

        ! Step 1: Apply center-of-mass correction if requested
        if (do_com_correction) then
            call correct_beta10_for_center_of_mass_s(beta_local, norm_constants, gl_nodes, gl_weights, legendre_gl)
        end if

        ! Store the (possibly corrected) β₁₀
        corrected_beta10 = beta_local(1)

        ! Step 2: Validate parameters before computing grid (poles check)
        call validate_legendre_params_s(beta_local, norm_constants, error_code, message)
        if (error_code /= LEGENDRE_VALID) then
            return
        end if

        ! Step 3: Compute R(θ) at grid points using pre-computed Legendre polynomials
        ! Precompute beta × normalization constant product
        beta_con(:) = beta_local(:) * norm_constants(:)

        h = PI_C / real(n_grid - 1_ik, rk)

        do i = 1_ik, n_grid
            ! Use pre-computed Legendre polynomials if available and grid matches
            if (n_grid == n_grid_precomp) then
                r_raw = 1.0_rk
                do k = 1, n_def_params
                    r_raw = r_raw + beta_con(k) * legendre_theta_grid(i, k + 1)
                end do
            else
                ! Compute directly for non-standard grid sizes
                theta = real(i - 1_ik, rk) * h
                call compute_radius_at_theta_s(beta_local, norm_constants, theta, r_raw)
            end if
            radii(i) = r_raw
        end do

        !-----------------------------------------------------------------------
        ! Step 4: Defense in Depth - Validate ALL radii are positive
        !-----------------------------------------------------------------------
        ! This catches pathological shapes that might pass polar checks but
        ! have negative or very small radii at interior points due to large
        ! higher-order deformations.

        r_min = huge(1.0_rk)
        i_min = 1_ik

        do i = 1_ik, n_grid
            if (radii(i) < r_min) then
                r_min = radii(i)
                i_min = i
            end if
        end do

        if (r_min <= R_MIN_THRESHOLD) then
            is_valid = .false.
            if (i_min == 1_ik) then
                write(message, '(A,ES12.4)') 'North pole radius not positive: R(0) = ', radii(1)
            else if (i_min == n_grid) then
                write(message, '(A,ES12.4)') 'South pole radius not positive: R(pi) = ', radii(n_grid)
            else
                theta = real(i_min - 1_ik, rk) * h
                write(message, '(A,F8.4,A,ES12.4)') &
                        'Interior radius not positive at theta = ', theta, ', R = ', r_min
            end if
            return
        end if

        is_valid = .true.
        message = ''

    end subroutine compute_legendre_radius_grid_s

    !===========================================================================
    ! STANDALONE WRAPPER
    !===========================================================================

    !> Standalone wrapper that internalizes all precomputation.
    !!
    !! Unlike compute_legendre_radius_grid_s, this routine does not require
    !! the caller to precompute normalization constants, Legendre polynomial
    !! grids, or Gauss-Legendre quadrature nodes/weights.
    !!
    !! @param[in]  params              Shape parameters (β values), λ=1..size(params)
    !! @param[in]  n_grid              Number of grid points for R(θ) output
    !! @param[out] radii               R(θ) values at grid points (size n_grid)
    !! @param[out] corrected_beta10    The (possibly corrected) β₁₀ value
    !! @param[out] is_valid            True if shape is valid
    !! @param[out] message             Error message if invalid
    !! @param[in]  apply_com_correction Optional: apply COM correction (default .false.)
    !! @param[in]  n_quad_in           Optional: number of GL quadrature points
    subroutine compute_legendre_radius_grid_standalone_s( &
            params, n_grid, &
            radii, corrected_beta10, is_valid, message, &
            apply_com_correction, n_quad_in)

        use mathematical_utilities_mod, only: &
                compute_spherical_harmonics_normalization_constants_s, &
                compute_gauss_legendre_quadrature_s, &
                compute_legendre_polynomials_s

        implicit none

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_grid
        real(kind = rk), intent(out) :: radii(n_grid)
        real(kind = rk), intent(out) :: corrected_beta10
        logical, intent(out) :: is_valid
        character(len = *), intent(out) :: message
        logical, intent(in), optional :: apply_com_correction
        integer(kind = ik), intent(in), optional :: n_quad_in

        ! Local variables
        integer(kind = ik) :: n_def_params, n_quad, i
        real(kind = rk), allocatable :: norm_constants(:)
        real(kind = rk), allocatable :: gl_nodes(:), gl_weights(:)
        real(kind = rk), allocatable :: legendre_theta_grid(:, :)
        real(kind = rk), allocatable :: legendre_gl(:, :)
        real(kind = rk) :: h, theta, x

        n_def_params = size(params, kind = ik)

        ! Guard against empty parameter array
        if (n_def_params < 1_ik) then
            radii = 0.0_rk
            corrected_beta10 = 0.0_rk
            is_valid = .false.
            message = 'Empty parameter array'
            return
        end if

        ! Compute normalization constants C_λ
        allocate(norm_constants(n_def_params))
        call compute_spherical_harmonics_normalization_constants_s( &
                norm_constants, n_def_params)

        ! Determine number of GL quadrature points
        if (present(n_quad_in)) then
            n_quad = n_quad_in
        else
            n_quad = max(2_ik * n_def_params + 1_ik, 20_ik)
        end if

        ! Compute GL nodes and weights
        allocate(gl_nodes(n_quad), gl_weights(n_quad))
        call compute_gauss_legendre_quadrature_s(n_quad, gl_nodes, gl_weights)

        ! Precompute Legendre polynomials on the θ grid
        h = PI_C / real(n_grid - 1_ik, rk)
        allocate(legendre_theta_grid(n_grid, n_def_params + 1))
        do i = 1_ik, n_grid
            theta = real(i - 1_ik, rk) * h
            x = cos(theta)
            call compute_legendre_polynomials_s( &
                    n_def_params + 1_ik, x, legendre_theta_grid(i, :))
        end do

        ! Precompute Legendre polynomials at GL nodes
        allocate(legendre_gl(n_quad, n_def_params + 1))
        do i = 1_ik, n_quad
            call compute_legendre_polynomials_s( &
                    n_def_params + 1_ik, gl_nodes(i), legendre_gl(i, :))
        end do

        ! Delegate to the existing routine
        call compute_legendre_radius_grid_s( &
                params, n_grid, norm_constants, legendre_theta_grid, &
                gl_nodes, gl_weights, legendre_gl, &
                radii, corrected_beta10, is_valid, message, &
                apply_com_correction)

    end subroutine compute_legendre_radius_grid_standalone_s

    !===========================================================================
    ! SHAPE VALIDATION
    !===========================================================================

    !> Validates Legendre/β-parameter shape parameters.
    !!
    !! Checks that the shape described by the parameters is physically valid:
    !! 1. Polar radii must be positive
    !! 2. No severe shape pathologies
    !!
    !! @param[in]  params      Shape parameters (β values)
    !! @param[out] error_code  Validation error code
    !! @param[out] message     Human-readable error message
    subroutine validate_legendre_params_s(params, norm_constants, error_code, message)

        use mathematical_utilities_mod, only: compute_legendre_polynomials_s

        implicit none

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: norm_constants(:)
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message

        integer(kind = ik) :: n_def_params, k
        real(kind = rk), allocatable :: legendre_polys(:)
        real(kind = rk) :: r_north, r_south, r_raw

        n_def_params = size(norm_constants)
        allocate(legendre_polys(n_def_params + 1))

        error_code = LEGENDRE_VALID
        message = ''

        ! Compute radius at north pole (x = cos(0) = 1)
        call compute_legendre_polynomials_s(n_def_params + 1_ik, 1.0_rk, legendre_polys)

        r_raw = 1.0_rk
        do k = 1, n_def_params
            r_raw = r_raw + params(k) * norm_constants(k) * legendre_polys(k + 1)
        end do
        r_north = r_raw

        if (r_north <= 0.0_rk) then
            error_code = LEGENDRE_ERROR_NORTH_POLE
            write(message, '(A,ES12.4)') 'North pole radius not positive: R(0) = ', r_north
            return
        end if

        ! Compute radius at south pole (x = cos(π) = -1)
        call compute_legendre_polynomials_s(n_def_params + 1_ik, -1.0_rk, legendre_polys)

        r_raw = 1.0_rk
        do k = 1, n_def_params
            r_raw = r_raw + params(k) * norm_constants(k) * legendre_polys(k + 1)
        end do
        r_south = r_raw

        if (r_south <= 0.0_rk) then
            error_code = LEGENDRE_ERROR_SOUTH_POLE
            write(message, '(A,ES12.4)') 'South pole radius not positive: R(pi) = ', r_south
            return
        end if

    end subroutine validate_legendre_params_s

    !===========================================================================
    ! RADIUS COMPUTATION
    !===========================================================================

    !> Computes radius at a single theta angle.
    !!
    !! @param[in]  params   Shape parameters (β values)
    !! @param[in]  theta    Polar angle θ ∈ [0, π]
    !! @param[out] r        Radius R(θ)
    subroutine compute_radius_at_theta_s(params, norm_constants, theta, r)

        use mathematical_utilities_mod, only: compute_legendre_polynomials_s

        implicit none

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: norm_constants(:)
        real(kind = rk), intent(in) :: theta
        real(kind = rk), intent(out) :: r

        integer(kind = ik) :: n_def_params, k
        real(kind = rk) :: x
        real(kind = rk), allocatable :: legendre_polys(:)

        n_def_params = size(norm_constants)
        allocate(legendre_polys(n_def_params + 1))

        x = cos(theta)
        call compute_legendre_polynomials_s(n_def_params + 1_ik, x, legendre_polys)

        r = 1.0_rk
        do k = 1, n_def_params
            r = r + params(k) * norm_constants(k) * legendre_polys(k + 1)
        end do

    end subroutine compute_radius_at_theta_s

    !> Computes nuclear radius using Legendre polynomial parameterization.
    !!
    !! Legacy interface conforming to the `radius_interface` abstract interface.
    !!
    !! @param[in]  params   Shape parameters (uses first number_of_deformation_parameters elements)
    !! @param[in]  x        Evaluation point cos(θ), -1 ≤ x ≤ 1
    !! @param[out] r        Nuclear radius R(x)
    pure subroutine compute_radius_legendre_s(params, norm_constants, x, r)

        use mathematical_utilities_mod, only: compute_legendre_polynomials_s

        implicit none

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: norm_constants(:)
        real(kind = rk), intent(in), value :: x
        real(kind = rk), intent(out) :: r

        integer(kind = ik) :: n_def_params, k
        real(kind = rk), allocatable :: legendre_polys(:), beta_con(:)

        n_def_params = size(norm_constants)
        allocate(legendre_polys(n_def_params + 1), beta_con(n_def_params))

        ! Compute Legendre polynomials
        call compute_legendre_polynomials_s(n_def_params + 1_ik, x, legendre_polys)

        ! Precompute deformation parameter times normalization constant
        beta_con(:) = params(1:n_def_params) * norm_constants(:)

        ! Compute radius: R(x) = 1 + Σ β_λ C_λ P_λ(x)
        r = 1.0_rk
        do k = 1, n_def_params
            r = r + beta_con(k) * legendre_polys(k + 1)
        end do

    end subroutine compute_radius_legendre_s

    !===========================================================================
    ! CENTER-OF-MASS CORRECTION
    !===========================================================================

    !> Corrects β₁₀ (dipole deformation parameter) for center-of-mass motion.
    !!
    !! This subroutine iteratively adjusts the dipole deformation parameter β₁₀
    !! to ensure that the center-of-mass of the nucleus is at the origin (z = 0).
    !!
    !! ## Algorithm
    !!
    !! 1. Compute z_cm = 3⟨z⟩β₀³/(8π)
    !! 2. If |z_cm| > tolerance: β₁ ← β₁ - z_cm/1.5
    !! 3. Repeat until converged or max iterations reached
    !!
    !! @param[inout] beta_params  Shape parameters; beta_params(1) is modified
    pure subroutine correct_beta10_for_center_of_mass_s(&
            beta_params, norm_constants, gl_nodes, gl_weights, legendre_gl)

        implicit none

        real(kind = rk), intent(inout) :: beta_params(:)
        real(kind = rk), intent(in) :: norm_constants(:)
        real(kind = rk), intent(in) :: gl_nodes(:)
        real(kind = rk), intent(in) :: gl_weights(:)
        real(kind = rk), intent(in) :: legendre_gl(:, :)

        ! Local constants
        real(kind = rk), parameter :: CM_TOLERANCE = 1.0e-5_rk
        real(kind = rk), parameter :: BETA10_ADJUSTMENT_FACTOR = 1.5_rk
        integer(kind = ik), parameter :: MAX_ITERATIONS = 50_ik

        ! Local variables
        integer(kind = ik) :: n_def_params, n_quad
        real(kind = rk) :: z_mean, z_center_of_mass, volume_conservation_factor
        real(kind = rk) :: radius, volume_integral, z_mean_integral
        real(kind = rk), allocatable :: beta_con(:)
        integer(kind = ik) :: iteration_count, i, k

        n_def_params = size(norm_constants)
        n_quad = size(gl_nodes)
        allocate(beta_con(n_def_params))

        ! Precompute beta × normalization constant product
        beta_con(:) = beta_params(:) * norm_constants(:)

        ! Initial calculation
        volume_integral = 0.0_rk
        z_mean_integral = 0.0_rk
        do i = 1, n_quad
            radius = 1.0_rk
            do k = 1, n_def_params
                radius = radius + beta_con(k) * legendre_gl(i, k + 1)
            end do
            volume_integral = volume_integral + (radius**3) * gl_weights(i)
            z_mean_integral = z_mean_integral + gl_nodes(i) * (radius**4) * gl_weights(i)
        end do
        volume_conservation_factor = (2.0_rk / volume_integral)**(1.0_rk / 3.0_rk)
        z_mean = z_mean_integral * PI_C
        z_center_of_mass = 3.0_rk * z_mean * volume_conservation_factor**3 / (8.0_rk * PI_C)

        ! Iteratively adjust beta_10
        iteration_count = 0_ik
        do while (abs(z_center_of_mass) >= CM_TOLERANCE .and. iteration_count < MAX_ITERATIONS)
            iteration_count = iteration_count + 1_ik
            beta_params(1) = beta_params(1) - z_center_of_mass / BETA10_ADJUSTMENT_FACTOR

            ! Update beta_con for the modified beta_params(1)
            beta_con(1) = beta_params(1) * norm_constants(1)

            ! Recompute
            volume_integral = 0.0_rk
            z_mean_integral = 0.0_rk
            do i = 1, n_quad
                radius = 1.0_rk
                do k = 1, n_def_params
                    radius = radius + beta_con(k) * legendre_gl(i, k + 1)
                end do
                volume_integral = volume_integral + (radius**3) * gl_weights(i)
                z_mean_integral = z_mean_integral + gl_nodes(i) * (radius**4) * gl_weights(i)
            end do
            volume_conservation_factor = (2.0_rk / volume_integral)**(1.0_rk / 3.0_rk)
            z_mean = z_mean_integral * PI_C
            z_center_of_mass = 3.0_rk * z_mean * volume_conservation_factor**3 / (8.0_rk * PI_C)
        end do

    end subroutine correct_beta10_for_center_of_mass_s

end module beta_parameterization_mod