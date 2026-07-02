!> COM-shift round-trip: corrected β₁₀ must re-yield |z_cm| < 1e-5,
!! recomputed independently of the library's own quadrature tables.
program beta_param_com_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_gauss_legendre_quadrature_s
    use beta_parameterization_mod, only: cache_t, LEGENDRE_VALID
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary

    implicit none

    integer(kind = ik), parameter :: N_QUAD = 512_ik
    real(kind = rk) :: gl_x(N_QUAD), gl_w(N_QUAD)

    call compute_gauss_legendre_quadrature_s(N_QUAD, gl_x, gl_w)

    call test_symmetric_shape_unchanged()
    call test_asymmetric_round_trip()
    call test_summary()

contains

    pure function norm_c(lambda) result(n)
        integer(kind = ik), intent(in) :: lambda
        real(kind = rk) :: n
        n = sqrt(real(2_ik * lambda + 1_ik, rk) / (4.0_rk * PI_C))
    end function norm_c

    pure function legendre_p(lambda, x) result(p)
        integer(kind = ik), intent(in) :: lambda
        real(kind = rk),    intent(in) :: x
        real(kind = rk) :: p
        select case (lambda)
        case (1_ik); p = x
        case (2_ik); p = (3.0_rk * x**2 - 1.0_rk) / 2.0_rk
        case (3_ik); p = (5.0_rk * x**3 - 3.0_rk * x) / 2.0_rk
        case (4_ik); p = (35.0_rk * x**4 - 30.0_rk * x**2 + 3.0_rk) / 8.0_rk
        case default; p = 0.0_rk
        end select
    end function legendre_p

    !> z_cm of the shape (β₁..β₄), same formula the workers use, own quadrature.
    function z_cm_of(betas) result(z_cm)
        real(kind = rk), intent(in) :: betas(4)
        real(kind = rk) :: z_cm
        real(kind = rk) :: r, vol_int, z_int, vol_factor
        integer(kind = ik) :: i, lam

        vol_int = 0.0_rk
        z_int   = 0.0_rk
        do i = 1_ik, N_QUAD
            r = 1.0_rk
            do lam = 1_ik, 4_ik
                r = r + betas(lam) * norm_c(lam) * legendre_p(lam, gl_x(i))
            end do
            vol_int = vol_int + r**3 * gl_w(i)
            z_int   = z_int + gl_x(i) * r**4 * gl_w(i)
        end do
        vol_factor = (2.0_rk / vol_int)**(1.0_rk / 3.0_rk)
        z_cm = 3.0_rk * z_int * vol_factor**3 / 8.0_rk
    end function z_cm_of

    !> Even-λ-only shape: z_cm = 0 by parity → zero iterations, β₁₀ unchanged (= 0).
    subroutine test_symmetric_shape_unchanged()
        type(cache_t)        :: cache
        real(kind = rk)      :: params(4), radii(181), corrected
        integer(kind = ik)   :: code
        character(len = 256) :: message

        call cache%init(4_ik, 181_ik, code, message)
        params = [0.0_rk, 0.50_rk, 0.0_rk, 0.10_rk]
        call cache%compute_radius_grid_with_com_shift(params, radii, corrected, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'com symmetric: valid')
        call assert_close(corrected, 0.0_rk, 1.0e-12_rk, 'com symmetric: beta10 stays 0')
    end subroutine test_symmetric_shape_unchanged

    !> Asymmetric shape: library's corrected β₁₀ must satisfy |z_cm| < CM_TOLERANCE.
    subroutine test_asymmetric_round_trip()
        type(cache_t)        :: cache
        real(kind = rk)      :: params(4), betas_corr(4), radii(181), corrected, z_cm
        integer(kind = ik)   :: code
        character(len = 256) :: message

        call cache%init(4_ik, 181_ik, code, message)
        params = [0.30_rk, 0.60_rk, 0.40_rk, 0.10_rk]
        call cache%compute_radius_grid_with_com_shift(params, radii, corrected, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'com asymmetric: valid')
        call assert_true(abs(corrected - params(1)) > 1.0e-6_rk, 'com asymmetric: beta10 actually moved')

        betas_corr    = params
        betas_corr(1) = corrected
        z_cm = z_cm_of(betas_corr)
        call assert_true(abs(z_cm) < 1.0e-5_rk, 'com asymmetric: |z_cm(corrected)| < 1e-5')
    end subroutine test_asymmetric_round_trip

end program beta_param_com_test
