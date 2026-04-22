!> Pure numerical kernels for the beta parameterization library.
!!
!! This module contains the inner-loop math: radius evaluation on a grid,
!! polar radius from analytic Legendre values, COM quadrature integrals, and
!! the COM iteration. No derived types, no I/O, no `iso_c_binding`,
!! no allocation in the hot path.
module beta_parameterization_workers_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_legendre_polynomials_s

    implicit none

    private

    ! Procedures (uncommented as each is added in Tasks 2–4)
    public :: precompute_legendre_table_s
    public :: eval_polar_radii_s
    public :: eval_radius_grid_s
    public :: find_min_radius_s
    public :: compute_com_integrals_s
    public :: iterate_com_correction_s

    !---------------------------------------------------------------------------
    ! COM-iteration algorithmic parameters
    !---------------------------------------------------------------------------
    real(kind = rk),    parameter :: CM_TOLERANCE             = 1.0e-5_rk
    real(kind = rk),    parameter :: BETA10_ADJUSTMENT_FACTOR = 1.5_rk
    integer(kind = ik), parameter :: MAX_ITERATIONS           = 50_ik

contains

    !> Fill `table(i, k+1) = P_k(x_values(i))` for k = 0..max_lambda.
    !!
    !! Used twice from cache_init: once for the theta grid, once for GL nodes.
    !!
    !! @param[in]  x_values    Evaluation points (any range; typically cos(theta) ∈ [-1, 1])
    !! @param[in]  max_lambda  Highest Legendre order to compute
    !! @param[out] table       Output table; column k+1 holds P_k. Shape (size(x_values), max_lambda + 1)
    pure subroutine precompute_legendre_table_s(x_values, max_lambda, table)

        real(kind = rk),    intent(in)  :: x_values(:)
        integer(kind = ik), intent(in)  :: max_lambda
        real(kind = rk),    intent(out) :: table(:, :)

        integer(kind = ik) :: i

        do i = 1_ik, size(x_values, kind = ik)
            call compute_legendre_polynomials_s(max_lambda + 1_ik, x_values(i), table(i, :))
        end do

    end subroutine precompute_legendre_table_s

    !> Compute R(0) and R(π) analytically from beta×normalization products.
    !!
    !! Uses P_λ(1) = 1 and P_λ(-1) = (-1)^λ — no Legendre evaluation needed.
    !! O(size(beta_con)). Called for early rejection before the grid evaluation.
    !!
    !! @param[in]  beta_con  beta_local(k) * norm_constants(k), for k = 1..max_beta_params
    !! @param[out] r_north   R(theta = 0) = 1 + Σ beta_con(k)
    !! @param[out] r_south   R(theta = π) = 1 + Σ beta_con(k) * (-1)^k
    pure subroutine eval_polar_radii_s(beta_con, r_north, r_south)

        real(kind = rk), intent(in)  :: beta_con(:)
        real(kind = rk), intent(out) :: r_north
        real(kind = rk), intent(out) :: r_south

        integer(kind = ik) :: k

        r_north = 1.0_rk
        r_south = 1.0_rk
        do k = 1_ik, size(beta_con, kind = ik)
            r_north = r_north + beta_con(k)
            if (mod(k, 2_ik) == 0_ik) then
                r_south = r_south + beta_con(k)
            else
                r_south = r_south - beta_con(k)
            end if
        end do

    end subroutine eval_polar_radii_s

    !> Evaluate R(theta) on the precomputed theta grid.
    !!
    !! radii(i) = 1 + Σ_k beta_con(k) * legendre_theta_grid(i, k+1)
    !!
    !! The inner loop. Pure, no allocation, no branch. The compiler is free to
    !! vectorize over `i` and unroll over `k` (small, ≤ 64).
    !!
    !! @param[in]  beta_con              beta×norm products (size = max_beta_params)
    !! @param[in]  legendre_theta_grid   Precomputed P_k(cos(theta_i)). Shape (n_grid, max_beta_params + 1)
    !! @param[out] radii                 R(theta) values (size = n_grid)
    pure subroutine eval_radius_grid_s(beta_con, legendre_theta_grid, radii)

        real(kind = rk), intent(in)  :: beta_con(:)
        real(kind = rk), intent(in)  :: legendre_theta_grid(:, :)
        real(kind = rk), intent(out) :: radii(:)

        integer(kind = ik) :: i, k, n_grid, n_lambda

        n_grid   = size(radii,    kind = ik)
        n_lambda = size(beta_con, kind = ik)

        do i = 1_ik, n_grid
            radii(i) = 1.0_rk
            do k = 1_ik, n_lambda
                radii(i) = radii(i) + beta_con(k) * legendre_theta_grid(i, k + 1_ik)
            end do
        end do

    end subroutine eval_radius_grid_s

    !> Return the minimum radius and its grid index.
    !!
    !! Caller decides whether `r_min` crosses the validity threshold and which
    !! error code to use based on `i_min` (1 → north pole, n_grid → south pole,
    !! interior otherwise).
    !!
    !! @param[in]  radii  Radius values
    !! @param[out] r_min  Minimum value
    !! @param[out] i_min  1-based index of the minimum (first occurrence)
    pure subroutine find_min_radius_s(radii, r_min, i_min)

        real(kind = rk),    intent(in)  :: radii(:)
        real(kind = rk),    intent(out) :: r_min
        integer(kind = ik), intent(out) :: i_min

        integer(kind = ik) :: i

        r_min = radii(1)
        i_min = 1_ik
        do i = 2_ik, size(radii, kind = ik)
            if (radii(i) < r_min) then
                r_min = radii(i)
                i_min = i
            end if
        end do

    end subroutine find_min_radius_s

    !> One pass of Gauss–Legendre quadrature for ∫R(θ)³ and ∫z·R(θ)⁴.
    !!
    !! The two integrals feed the volume-conservation factor and the
    !! center-of-mass z coordinate used by `iterate_com_correction_s`.
    !!
    !! @param[in]  beta_con          beta×norm products (size = max_beta_params)
    !! @param[in]  gl_nodes          Gauss–Legendre nodes (size = n_quad)
    !! @param[in]  gl_weights        Gauss–Legendre weights (size = n_quad)
    !! @param[in]  legendre_gl       P_k at the GL nodes. Shape (n_quad, max_beta_params + 1)
    !! @param[out] volume_integral   Σ R(x_i)³ * w_i
    !! @param[out] z_mean_integral   Σ x_i * R(x_i)⁴ * w_i
    pure subroutine compute_com_integrals_s( &
            beta_con, gl_nodes, gl_weights, legendre_gl, &
            volume_integral, z_mean_integral)

        real(kind = rk), intent(in)  :: beta_con(:)
        real(kind = rk), intent(in)  :: gl_nodes(:)
        real(kind = rk), intent(in)  :: gl_weights(:)
        real(kind = rk), intent(in)  :: legendre_gl(:, :)
        real(kind = rk), intent(out) :: volume_integral
        real(kind = rk), intent(out) :: z_mean_integral

        integer(kind = ik) :: i, k, n_quad, n_lambda
        real(kind = rk)    :: radius

        n_quad   = size(gl_nodes, kind = ik)
        n_lambda = size(beta_con, kind = ik)

        volume_integral = 0.0_rk
        z_mean_integral = 0.0_rk
        do i = 1_ik, n_quad
            radius = 1.0_rk
            do k = 1_ik, n_lambda
                radius = radius + beta_con(k) * legendre_gl(i, k + 1_ik)
            end do
            volume_integral = volume_integral + radius**3 * gl_weights(i)
            z_mean_integral = z_mean_integral + gl_nodes(i) * radius**4 * gl_weights(i)
        end do

    end subroutine compute_com_integrals_s

    !> Iteratively adjust β₁₀ until the center-of-mass z-coordinate is < CM_TOLERANCE.
    !!
    !! Modifies `beta_local(1)` and `beta_con(1)` in place. Reports convergence
    !! status and iteration count to the caller; **the caller decides** whether
    !! non-convergence is acceptable. The Fortran API in v2 treats it as a
    !! validation failure (`LEGENDRE_ERROR_COM_NOT_CONVERGED`).
    !!
    !! Algorithm: fixed-point iteration.
    !!   z_cm = 3 * z_mean_integral * vol_factor³ / 8     [PI_C cancels in original form]
    !!   β₁₀ ← β₁₀ - z_cm / BETA10_ADJUSTMENT_FACTOR
    !!
    !! @param[inout] beta_local     β values; only beta_local(1) is modified
    !! @param[inout] beta_con       beta×norm products; only beta_con(1) is modified
    !! @param[in]    norm_constants C_λ values (size = max_beta_params)
    !! @param[in]    gl_nodes       Gauss–Legendre nodes
    !! @param[in]    gl_weights     Gauss–Legendre weights
    !! @param[in]    legendre_gl    P_k at the GL nodes
    !! @param[out]   converged      .true. if |z_cm| < CM_TOLERANCE on exit
    !! @param[out]   n_iter         Number of iterations performed
    pure subroutine iterate_com_correction_s( &
            beta_local, beta_con, norm_constants, &
            gl_nodes, gl_weights, legendre_gl, &
            converged, n_iter)

        real(kind = rk),    intent(inout) :: beta_local(:)
        real(kind = rk),    intent(inout) :: beta_con(:)
        real(kind = rk),    intent(in)    :: norm_constants(:)
        real(kind = rk),    intent(in)    :: gl_nodes(:)
        real(kind = rk),    intent(in)    :: gl_weights(:)
        real(kind = rk),    intent(in)    :: legendre_gl(:, :)
        logical,            intent(out)   :: converged
        integer(kind = ik), intent(out)   :: n_iter

        real(kind = rk) :: volume_integral, z_mean_integral, z_cm, vol_factor

        n_iter = 0_ik

        call compute_com_integrals_s(beta_con, gl_nodes, gl_weights, legendre_gl, &
                volume_integral, z_mean_integral)
        vol_factor = (2.0_rk / volume_integral)**(1.0_rk / 3.0_rk)
        z_cm = 3.0_rk * z_mean_integral * vol_factor**3 / 8.0_rk

        do while (abs(z_cm) >= CM_TOLERANCE .and. n_iter < MAX_ITERATIONS)
            n_iter = n_iter + 1_ik
            beta_local(1) = beta_local(1) - z_cm / BETA10_ADJUSTMENT_FACTOR
            beta_con(1)   = beta_local(1) * norm_constants(1)

            call compute_com_integrals_s(beta_con, gl_nodes, gl_weights, legendre_gl, &
                    volume_integral, z_mean_integral)
            vol_factor = (2.0_rk / volume_integral)**(1.0_rk / 3.0_rk)
            z_cm = 3.0_rk * z_mean_integral * vol_factor**3 / 8.0_rk
        end do

        converged = abs(z_cm) < CM_TOLERANCE

    end subroutine iterate_com_correction_s

end module beta_parameterization_workers_mod
