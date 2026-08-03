!> Golden-independent physics properties: the COM root and volume conservation.
!!
!! Every reference number here is built test-side — own Gauss-Legendre rules,
!! own Legendre recurrence, own normalization constants. The suite touches no
!! library table and no stored expectation, so it fails only when the library's
!! two physics guarantees actually break.
program beta_param_property_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_gauss_legendre_quadrature_s, &
            compute_legendre_polynomials_s, &
            compute_spherical_harmonics_normalization_constants_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_resolve_shape_s, SHAPE_VALID

    implicit none

    integer(kind = ik), parameter :: N_PARAMS = 6_ik      !! 6-param cache, lambda = 1..6
    integer(kind = ik), parameter :: N_LEG    = 7_ik      !! P_0 .. P_6 in p(1..7)
    integer(kind = ik), parameter :: N_COM    = 512_ik    !! GL order for the COM reference
    integer(kind = ik), parameter :: N_VOL    = 128_ik    !! GL order for the volume reference
    integer(kind = ik), parameter :: N_THETA  = 32_ik     !! primary thetas for the COM cache

    !> CM_TOLERANCE, spec section 12: the library promises |z_cm| below this.
    real(kind = rk), parameter :: CM_TOLERANCE = 1.0e-5_rk

    real(kind = rk) :: com_x(N_COM), com_w(N_COM)
    real(kind = rk) :: vol_x(N_VOL), vol_w(N_VOL)
    real(kind = rk) :: norms(N_PARAMS)

    real(kind = rk), parameter :: SHAPE_A(N_PARAMS) = &
            [0.0_rk, 0.2_rk, 0.1_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: SHAPE_B(N_PARAMS) = &
            [0.1_rk, 0.3_rk, 0.15_rk, -0.05_rk, 0.05_rk, 0.02_rk]

    call compute_gauss_legendre_quadrature_s(N_COM, com_x, com_w)
    call compute_gauss_legendre_quadrature_s(N_VOL, vol_x, vol_w)
    call compute_spherical_harmonics_normalization_constants_s(norms, N_PARAMS)

    call check_com_property('shape A', SHAPE_A)
    call check_com_property('shape B', SHAPE_B)
    call check_volume_property('shape A', SHAPE_A)
    call check_volume_property('shape B', SHAPE_B)

    call test_summary()

contains

    !> R(x) = 1 + sum_k beta_con(k) P_k(x), evaluated with the test's own recurrence.
    function radius_at_f(beta_con, x) result(r)
        real(kind = rk), intent(in) :: beta_con(N_PARAMS)
        real(kind = rk), intent(in) :: x
        real(kind = rk) :: r
        real(kind = rk) :: p(N_LEG)
        integer(kind = ik) :: k

        call compute_legendre_polynomials_s(N_LEG, x, p)
        r = 1.0_rk
        do k = 1_ik, N_PARAMS
            r = r + beta_con(k) * p(k + 1_ik)
        end do
    end function radius_at_f

    !> z_cm of the shape defined by beta_con, from the test's own GL-512 rule:
    !! z_num = sum w x R^4, volume_integral = sum w R^3, vf = (2/volume_integral)^(1/3),
    !! z_cm = 3 z_num vf^3 / 8.
    function z_cm_of_f(beta_con) result(z_cm)
        real(kind = rk), intent(in) :: beta_con(N_PARAMS)
        real(kind = rk) :: z_cm
        real(kind = rk) :: r, volume_integral, z_num, vf
        integer(kind = ik) :: i

        volume_integral = 0.0_rk
        z_num           = 0.0_rk
        do i = 1_ik, N_COM
            r = radius_at_f(beta_con, com_x(i))
            volume_integral = volume_integral + com_w(i) * r**3
            z_num           = z_num + com_w(i) * com_x(i) * r**4
        end do
        vf   = (2.0_rk / volume_integral)**(1.0_rk / 3.0_rk)
        z_cm = 3.0_rk * z_num * vf**3 / 8.0_rk
    end function z_cm_of_f

    !> The corrected beta10 the library returns must put the centre of mass at
    !! the origin, measured with quadrature the library never sees.
    subroutine check_com_property(label, params)
        character(len = *), intent(in) :: label
        real(kind = rk),    intent(in) :: params(N_PARAMS)

        type(cache_t)      :: cache
        real(kind = rk)    :: thetas(N_THETA), beta_con(N_PARAMS)
        real(kind = rk)    :: corrected_beta10, r_north, r_south, volume_factor, z_cm
        integer(kind = ik) :: status, i

        ! Open uniform primary set: no poles, spacing irrelevant to the property.
        do i = 1_ik, N_THETA
            thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
        end do

        call cache_init_s(cache, N_PARAMS, thetas, .false., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'com ' // label // ': cache init')

        call cache_resolve_shape_s(cache, params, corrected_beta10, r_north, &
                r_south, volume_factor, status)
        call assert_int_eq(status, SHAPE_VALID, 'com ' // label // ': resolve')
        call assert_close(volume_factor, 1.0_rk, 0.0_rk, &
                'com ' // label // ': volume factor is 1 without conservation')

        ! Reference shape: library's corrected beta10, caller's beta2..beta6.
        do i = 1_ik, N_PARAMS
            beta_con(i) = params(i) * norms(i)
        end do
        beta_con(1) = corrected_beta10 * norms(1)

        z_cm = z_cm_of_f(beta_con)
        call assert_true(abs(z_cm) < CM_TOLERANCE, &
                'com ' // label // ': |z_cm(corrected)| < CM_TOLERANCE')

        call cache_free_s(cache)
    end subroutine check_com_property

    !> With volume conservation on, the scaled radii integrate to the unit-sphere
    !! volume: sum w R^3 = 2 over the same GL nodes the primary set was built from.
    subroutine check_volume_property(label, params)
        character(len = *), intent(in) :: label
        real(kind = rk),    intent(in) :: params(N_PARAMS)

        type(cache_t)      :: cache
        real(kind = rk)    :: thetas(N_VOL), radii(N_VOL), volume_integral
        integer(kind = ik) :: status, i

        ! Primary thetas = acos of the test's own GL-128 nodes, so the cache
        ! reports R exactly where the reference rule needs it. Nodes descend in
        ! x, so the thetas ascend; both are strictly interior, so none is polar.
        do i = 1_ik, N_VOL
            thetas(i) = acos(vol_x(i))
        end do

        call cache_init_s(cache, N_PARAMS, thetas, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'volume ' // label // ': cache init')

        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'volume ' // label // ': radius grid')

        volume_integral = 0.0_rk
        do i = 1_ik, N_VOL
            volume_integral = volume_integral + vol_w(i) * radii(i)**3
        end do
        call assert_close(volume_integral, 2.0_rk, 1.0e-12_rk, &
                'volume ' // label // ': scaled radii integrate to the sphere volume')

        call cache_free_s(cache)
    end subroutine check_volume_property

end program beta_param_property_test
