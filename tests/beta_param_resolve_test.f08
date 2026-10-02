program beta_param_resolve_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary
    use beta_parameterization_mod, only: &
            cache_t, cache_init_s, cache_free_s, cache_resolve_shape_s, &
            SHAPE_VALID, SHAPE_ERROR_CACHE_NOT_INITIALIZED, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_COM_NOT_CONVERGED
    implicit none

    type(cache_t) :: cache
    integer(kind = ik) :: status, i
    real(kind = rk) :: thetas(64), params(4), no_params(0)
    real(kind = rk) :: b10, rn, rs, vf

    do i = 1_ik, 64_ik
        thetas(i) = real(i, rk) * PI_C / 65.0_rk
    end do

    call cache_resolve_shape_s(cache, [0.0_rk], .true., .true., b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_ERROR_CACHE_NOT_INITIALIZED, 'uninit cache rejected')

    call cache_init_s(cache, 4_ik, thetas, status)
    params = [0.05_rk, 0.20_rk, 0.10_rk, 0.02_rk]
    call cache_resolve_shape_s(cache, params, .true., .true., b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_VALID, 'asymmetric shape resolves')
    call assert_true(abs(b10 - params(1)) > 1.0e-6_rk, 'COM moved beta10')
    call assert_true(rn > 0.0_rk .and. rs > 0.0_rk, 'poles positive')
    call assert_true(vf > 0.0_rk .and. abs(vf - 1.0_rk) < 0.5_rk, 'volume factor sane')

    ! options are per call: the same cache serves the raw resolve too
    call cache_resolve_shape_s(cache, params, .false., .false., b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_VALID, 'raw resolve on the same cache')
    call assert_close(b10, params(1), 0.0_rk, 'no COM: beta10 untouched')
    call assert_close(vf, 1.0_rk, 0.0_rk, 'no volume conservation: factor exactly 1')

    ! sphere: no COM motion, factor exactly 1 before scaling arithmetic
    call cache_resolve_shape_s(cache, [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], .true., .true., &
            b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_VALID, 'sphere resolves')
    call assert_close(b10, 0.0_rk, 1.0e-14_rk, 'sphere beta10 stays 0')
    call assert_close(vf, 1.0_rk, 1.0e-12_rk, 'sphere volume factor 1')
    call assert_close(rn, 1.0_rk, 1.0e-12_rk, 'sphere north pole 1')

    ! a short vector is accepted: missing trailing parameters are zero
    call cache_resolve_shape_s(cache, [0.0_rk, 0.2_rk], .true., .true., b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_VALID, 'short vector accepted')

    ! outside 1..max_params: wrong parameter count
    call cache_resolve_shape_s(cache, [0.0_rk, 0.2_rk, 0.0_rk, 0.0_rk, 0.0_rk], &
            .true., .true., b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'five parameters on a 4-cache')
    call cache_resolve_shape_s(cache, no_params, .true., .true., b10, rn, rs, vf, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'empty vector on a cache')

    ! genuinely invalid shape: beta20 = -1.8 drives R(0) = 1 + C2*beta20 < 0
    ! (C2 ~ 0.631; -0.99 would NOT be invalid — R(0) ~ 0.38).
    call cache_resolve_shape_s(cache, [0.0_rk, -1.8_rk, 0.0_rk, 0.0_rk], .true., .true., &
            b10, rn, rs, vf, status)
    call assert_int_eq(status, BETA_PARAM_ERROR_NORTH_POLE, 'invalid shape rejected with 100')
    call assert_close(b10, 0.0_rk, 0.0_rk, 'beta10 zero-filled on failure')
    call assert_close(rn, 0.0_rk, 0.0_rk, 'r_north zero-filled on failure')
    call assert_close(rs, 0.0_rk, 0.0_rk, 'r_south zero-filled on failure')
    call assert_close(vf, 0.0_rk, 0.0_rk, 'volume factor zero-filled on failure')

    ! COM guard: beta20 = -20 makes the volume integral negative. The Newton
    ! step must report non-convergence instead of taking the cube root of a
    ! negative number (a SIGFPE in Debug before the guard).
    call cache_resolve_shape_s(cache, [0.0_rk, -20.0_rk], .false., .true., &
            b10, rn, rs, vf, status)
    call assert_int_eq(status, BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            'negative volume integral with COM rejected with 103')
    call cache_resolve_shape_s(cache, [0.0_rk, -20.0_rk], .false., .false., &
            b10, rn, rs, vf, status)
    call assert_int_eq(status, BETA_PARAM_ERROR_NORTH_POLE, &
            'same vector without COM rejected with 100')

    call cache_free_s(cache)
    call test_summary()
end program beta_param_resolve_test
