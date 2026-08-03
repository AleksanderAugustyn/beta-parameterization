program beta_param_resolve_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary
    use beta_parameterization_mod, only: &
            cache_t, cache_init_s, cache_free_s, cache_resolve_shape_s, &
            SHAPE_VALID, SHAPE_ERROR_CACHE_NOT_INITIALIZED, &
            SHAPE_ERROR_WRONG_PARAM_COUNT
    implicit none

    type(cache_t) :: cache
    integer(kind = ik) :: status, i
    real(kind = rk) :: thetas(64), params(4)
    real(kind = rk) :: b10, rn, rs, vf

    do i = 1_ik, 64_ik
        thetas(i) = real(i, rk) * PI_C / 65.0_rk
    end do

    call cache_resolve_shape_s(cache, [0.0_rk], b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_ERROR_CACHE_NOT_INITIALIZED, 'uninit cache rejected')

    call cache_init_s(cache, 4_ik, thetas, .true., .true., status)
    params = [0.05_rk, 0.20_rk, 0.10_rk, 0.02_rk]
    call cache_resolve_shape_s(cache, params, b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_VALID, 'asymmetric shape resolves')
    call assert_true(abs(b10 - params(1)) > 1.0e-6_rk, 'COM moved beta10')
    call assert_true(rn > 0.0_rk .and. rs > 0.0_rk, 'poles positive')
    call assert_true(vf > 0.0_rk .and. abs(vf - 1.0_rk) < 0.5_rk, 'volume factor sane')

    ! sphere: no COM motion, factor exactly 1 before scaling arithmetic
    call cache_resolve_shape_s(cache, [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], &
            b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_VALID, 'sphere resolves')
    call assert_close(b10, 0.0_rk, 1.0e-14_rk, 'sphere beta10 stays 0')
    call assert_close(vf, 1.0_rk, 1.0e-12_rk, 'sphere volume factor 1')
    call assert_close(rn, 1.0_rk, 1.0e-12_rk, 'sphere north pole 1')

    call cache_resolve_shape_s(cache, [0.0_rk, 0.2_rk], b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'exact n_params enforced')

    ! genuinely invalid shape: beta20 = -1.8 drives R(0) = 1 + C2*beta20 < 0
    ! (C2 ~ 0.631; -0.99 would NOT be invalid — R(0) ~ 0.38). Expected failure
    ! is 100/101 from validation, or 103 if the COM pass trips first — assert
    ! nonzero, not a specific code.
    call cache_resolve_shape_s(cache, [0.0_rk, -1.8_rk, 0.0_rk, 0.0_rk], &
            b10, rn, rs, vf, status)
    call assert_true(status /= SHAPE_VALID, 'invalid shape rejected')
    call assert_close(rn, 0.0_rk, 1.0e-15_rk, 'outputs zero-filled on failure')

    call cache_free_s(cache)
    call test_summary()
end program beta_param_resolve_test
