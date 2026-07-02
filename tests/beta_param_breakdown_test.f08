!> Breakdown-boundary tests: the validity boundary the library reports must
!! coincide with the mathematical breakdown (R(θ) crossing zero), and inside
!! the boundary volume/surface must remain geometrically meaningful.
program beta_param_breakdown_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use beta_parameterization_mod, only: cache_t, LEGENDRE_VALID, &
            LEGENDRE_ERROR_INTERIOR_NEGATIVE, LEGENDRE_ERROR_NORTH_POLE, &
            LEGENDRE_ERROR_SOUTH_POLE
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_GRID = 2001_ik   ! odd: θ=π/2 on-grid
    real(kind = rk) :: n2

    n2 = sqrt(5.0_rk / (4.0_rk * PI_C))

    call test_interior_crossing_beta2()
    call test_pole_crossing_beta2_negative()
    call test_volume_surface_positive_inside()
    call test_summary()

contains

    !> Trapezoid rule on the uniform θ grid.
    pure function trapz(values, h) result(s)
        real(kind = rk), intent(in) :: values(:), h
        real(kind = rk) :: s
        integer(kind = ik) :: i
        s = 0.5_rk * (values(1) + values(size(values, kind = ik)))
        do i = 2_ik, size(values, kind = ik) - 1_ik
            s = s + values(i)
        end do
        s = s * h
    end function trapz

    !> β₂ > 0: R(π/2) = 1 - β₂N₂/2 crosses zero at β*₂ = 2/N₂ ≈ 3.1707.
    subroutine test_interior_crossing_beta2()
        type(cache_t)        :: cache
        real(kind = rk)      :: params(2), radii(N_GRID), beta_star
        integer(kind = ik)   :: code
        character(len = 256) :: message

        call cache%init(2_ik, N_GRID, code, message)
        beta_star = 2.0_rk / n2

        params = [0.0_rk, beta_star * (1.0_rk - 1.0e-3_rk)]
        call cache%compute_radius_grid(params, radii, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'breakdown: 0.999 beta* valid')

        params = [0.0_rk, beta_star * (1.0_rk + 1.0e-3_rk)]
        call cache%compute_radius_grid(params, radii, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_INTERIOR_NEGATIVE, &
                'breakdown: 1.001 beta* rejected as interior-negative')
    end subroutine test_interior_crossing_beta2

    !> β₂ < 0: both poles hit zero first at β₂ = -1/N₂ (P₂(±1) = 1).
    subroutine test_pole_crossing_beta2_negative()
        type(cache_t)        :: cache
        real(kind = rk)      :: params(2), radii(N_GRID), beta_pole
        integer(kind = ik)   :: code
        character(len = 256) :: message

        call cache%init(2_ik, N_GRID, code, message)
        beta_pole = -1.0_rk / n2

        params = [0.0_rk, beta_pole * (1.0_rk - 1.0e-3_rk)]
        call cache%compute_radius_grid(params, radii, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'breakdown: 0.999 pole crossing valid')

        params = [0.0_rk, beta_pole * (1.0_rk + 1.0e-3_rk)]
        call cache%compute_radius_grid(params, radii, code, message)
        call assert_true(code == LEGENDRE_ERROR_NORTH_POLE .or. code == LEGENDRE_ERROR_SOUTH_POLE, &
                'breakdown: 1.001 pole crossing rejected as pole error')
    end subroutine test_pole_crossing_beta2_negative

    !> Scan β₂ across (pole crossing, interior crossing): every shape the library
    !! accepts must yield V > 0 and S > 0, finite, from the returned radii.
    !! V = (2π/3)∫R³sinθ dθ,  S = 2π∫R sinθ √(R² + R'²) dθ  (R' via central diff).
    subroutine test_volume_surface_positive_inside()
        type(cache_t)        :: cache
        real(kind = rk)      :: params(2), radii(N_GRID)
        real(kind = rk)      :: h, theta, vol, surf, sphere_vol
        real(kind = rk)      :: integrand_v(N_GRID), integrand_s(N_GRID), dr(N_GRID)
        integer(kind = ik)   :: code, i, k, n_valid
        character(len = 256) :: message

        call cache%init(2_ik, N_GRID, code, message)
        h = PI_C / real(N_GRID - 1_ik, rk)
        sphere_vol = 4.0_rk * PI_C / 3.0_rk
        n_valid = 0_ik

        do k = -30_ik, 30_ik
            params = [0.0_rk, real(k, rk) * 0.1_rk]
            call cache%compute_radius_grid(params, radii, code, message)
            if (code /= LEGENDRE_VALID) cycle
            n_valid = n_valid + 1_ik

            dr(1)      = (radii(2) - radii(1)) / h
            dr(N_GRID) = (radii(N_GRID) - radii(N_GRID - 1_ik)) / h
            do i = 2_ik, N_GRID - 1_ik
                dr(i) = (radii(i + 1_ik) - radii(i - 1_ik)) / (2.0_rk * h)
            end do
            do i = 1_ik, N_GRID
                theta = real(i - 1_ik, rk) * h
                integrand_v(i) = radii(i)**3 * sin(theta)
                integrand_s(i) = radii(i) * sin(theta) * sqrt(radii(i)**2 + dr(i)**2)
            end do
            vol  = (2.0_rk * PI_C / 3.0_rk) * trapz(integrand_v, h)
            surf = 2.0_rk * PI_C * trapz(integrand_s, h)

            call assert_true(vol > 0.0_rk .and. vol < 100.0_rk * sphere_vol, &
                    'breakdown scan: volume positive and finite')
            call assert_true(surf > 0.0_rk .and. surf < 1000.0_rk, &
                    'breakdown scan: surface positive and finite')
        end do

        call assert_true(n_valid >= 40_ik, 'breakdown scan: scan covered the valid range')
    end subroutine test_volume_surface_positive_inside

end program beta_param_breakdown_test
