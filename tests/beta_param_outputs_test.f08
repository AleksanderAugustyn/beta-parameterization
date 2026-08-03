program beta_param_outputs_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_legendre_polynomials_s, &
            compute_spherical_harmonics_normalization_constants_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, cache_resolve_shape_s, &
            SHAPE_VALID, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE
    implicit none

    type(cache_t) :: cache
    integer(kind = ik) :: status, i
    real(kind = rk) :: thetas(16), params(4), radii(16), drs(16), bad(5)
    real(kind = rk) :: norms(4), p(5), expected, b10, rn, rs, vf

    do i = 1_ik, 16_ik
        thetas(i) = real(i, rk) * PI_C / 17.0_rk
    end do
    call compute_spherical_harmonics_normalization_constants_s(norms, 4_ik)
    params = [0.0_rk, 0.2_rk, 0.0_rk, 0.0_rk]

    ! no COM, no volume: pure table evaluation is exactly the Legendre sum
    call cache_init_s(cache, 4_ik, thetas, .false., .false., status)
    call cache_radius_grid_s(cache, params, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'radius grid ok')
    do i = 1_ik, 16_ik
        call compute_legendre_polynomials_s(5_ik, cos(thetas(i)), p)
        expected = 1.0_rk + params(2) * norms(2) * p(3)
        call assert_close(radii(i), expected, 1.0e-14_rk, 'R matches Legendre sum')
    end do

    call cache_radius_and_derivative_s(cache, params, radii, drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'radius+derivative ok')
    call assert_close(drs(8) + drs(9), 0.0_rk, 1.0e-12_rk, &
            'derivative antisymmetric for even-only shape')

    call cache_radius_grid_s(cache, params, bad(1:5), status)
    call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, 'bad buffer 104')
    call cache_free_s(cache)

    ! volume conservation scales outputs so that resolve reports the factor
    call cache_init_s(cache, 4_ik, thetas, .true., .false., status)
    call cache_resolve_shape_s(cache, params, b10, rn, rs, vf, status)
    call cache_radius_grid_s(cache, params, radii, status)
    call compute_legendre_polynomials_s(5_ik, cos(thetas(1)), p)
    expected = (1.0_rk + params(2) * norms(2) * p(3)) * vf
    call assert_close(radii(1), expected, 1.0e-13_rk, 'volume factor applied to radii')
    call cache_free_s(cache)
    call test_summary()
end program beta_param_outputs_test
