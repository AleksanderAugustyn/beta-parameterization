!> Core correctness: sphere, analytic low-λ shapes, padding equivalence.
program beta_param_core_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use beta_parameterization_mod, only: cache_t, LEGENDRE_VALID, &
            compute_radius_grid_standalone_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary

    implicit none

    call test_sphere()
    call test_analytic_low_lambda()
    call test_padding_equivalence()
    call test_summary()

contains

    !> All-zero betas must give R(θ) ≡ 1 exactly (no rounding enters the sum).
    subroutine test_sphere()
        type(cache_t)        :: cache
        real(kind = rk)      :: params(4), radii(181)
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call cache%init(8_ik, 181_ik, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'sphere: cache init')

        params(:) = 0.0_rk
        call cache%compute_radius_grid(params, radii, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'sphere: compute valid')
        do i = 1_ik, 181_ik
            call assert_close(radii(i), 1.0_rk, 1.0e-15_rk, 'sphere: R == 1')
        end do

        ! Standalone path must agree with the cache path.
        call compute_radius_grid_standalone_s(params, 181_ik, radii, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'sphere: standalone valid')
        call assert_close(radii(91), 1.0_rk, 1.0e-15_rk, 'sphere: standalone equator')
    end subroutine test_sphere

    !> Analytic N_λ = sqrt((2λ+1)/(4π)) — independent of the library and of FF.
    pure function norm_c(lambda) result(n)
        integer(kind = ik), intent(in) :: lambda
        real(kind = rk) :: n
        n = sqrt(real(2_ik * lambda + 1_ik, rk) / (4.0_rk * PI_C))
    end function norm_c

    !> Analytic P_λ(x) for λ = 1..4 — closed forms, no recurrence.
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

    !> R(θ) for a (β₂, β₃, β₄) shape must match the closed-form Legendre sum.
    subroutine test_analytic_low_lambda()
        type(cache_t)        :: cache
        real(kind = rk)      :: params(4), radii(181), theta, x, r_ref
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call cache%init(4_ik, 181_ik, code, message)
        params = [0.0_rk, 0.25_rk, 0.10_rk, 0.05_rk]
        call cache%compute_radius_grid(params, radii, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'analytic: compute valid')

        do i = 1_ik, 181_ik
            theta = real(i - 1_ik, rk) * (PI_C / 180.0_rk)
            x     = cos(theta)
            r_ref = 1.0_rk + params(2) * norm_c(2_ik) * legendre_p(2_ik, x) &
                           + params(3) * norm_c(3_ik) * legendre_p(3_ik, x) &
                           + params(4) * norm_c(4_ik) * legendre_p(4_ik, x)
            call assert_close(radii(i), r_ref, 1.0e-13_rk, 'analytic: R(theta) matches closed form')
        end do
    end subroutine test_analytic_low_lambda

    !> Short params must be zero-padded: (β₂) alone ≡ (β₂, 0, 0, 0) on the same cache.
    subroutine test_padding_equivalence()
        type(cache_t)        :: cache
        real(kind = rk)      :: short_params(2), full_params(4)
        real(kind = rk)      :: radii_short(91), radii_full(91)
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call cache%init(4_ik, 91_ik, code, message)
        short_params = [0.0_rk, 0.30_rk]
        full_params  = [0.0_rk, 0.30_rk, 0.0_rk, 0.0_rk]
        call cache%compute_radius_grid(short_params, radii_short, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'padding: short valid')
        call cache%compute_radius_grid(full_params, radii_full, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'padding: full valid')
        do i = 1_ik, 91_ik
            call assert_close(radii_short(i), radii_full(i), 0.0_rk, 'padding: bit-identical')
        end do
    end subroutine test_padding_equivalence

end program beta_param_core_test
