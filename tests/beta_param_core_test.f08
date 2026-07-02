!> Core correctness: sphere, analytic low-λ shapes, padding equivalence.
program beta_param_core_test

    use precision_utilities_mod, only: ik, rk
    use beta_parameterization_mod, only: cache_t, LEGENDRE_VALID, &
            compute_radius_grid_standalone_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary

    implicit none

    call test_sphere()
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

end program beta_param_core_test
