!> Tier-1 standalone entries: no cache, no engine, same math as the cached path.
program beta_param_standalone_test
    use, intrinsic :: iso_fortran_env, only: int64
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_INVALID_INIT, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE
    implicit none

    type(cache_t) :: cache
    integer(kind = ik) :: status, i
    real(kind = rk) :: thetas(16), params(4)
    real(kind = rk) :: radii_sa(16), drs_sa(16), radii_c(16), drs_c(16)
    real(kind = rk) :: many_params(65), no_params(0)
    real(kind = rk) :: one_theta(1), one_radius(1)
    real(kind = rk) :: dr_bad(5)

    do i = 1_ik, 16_ik
        thetas(i) = real(i, rk) * PI_C / 17.0_rk
    end do
    params = [0.05_rk, 0.2_rk, 0.1_rk, 0.02_rk]

    !---------------------------------------------------------------------------
    ! 1. Standalone reproduces the cached pipeline bit for bit: same tables, same
    !    kernels, same order of operations, so the results must be identical.
    !---------------------------------------------------------------------------
    call cache_init_s(cache, 4_ik, thetas, .true., .true., status)
    call assert_int_eq(status, SHAPE_VALID, 'cached init')
    call cache_radius_and_derivative_s(cache, params, radii_c, drs_c, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached radius+derivative')

    call compute_radius_and_derivative_standalone_s(params, thetas, .true., .true., &
            radii_sa, drs_sa, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone radius+derivative')
    call assert_true(bits_equal_f(radii_sa, radii_c), &
            'standalone radii bit-identical to cached')
    call assert_true(bits_equal_f(drs_sa, drs_c), &
            'standalone derivatives bit-identical to cached')

    call cache_radius_grid_s(cache, params, radii_c, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached radius grid')
    call cache_free_s(cache)

    radii_sa(:) = 0.0_rk
    call compute_radius_grid_standalone_s(params, thetas, .true., .true., radii_sa, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone radius grid')
    call assert_true(bits_equal_f(radii_sa, radii_c), &
            'standalone grid bit-identical to cached')

    !---------------------------------------------------------------------------
    ! 2. Above the tier-1 parameter cap: status 1, outputs zero-filled.
    !---------------------------------------------------------------------------
    many_params(:) = 0.001_rk
    radii_sa(:)    = 1.0_rk
    call compute_radius_grid_standalone_s(many_params, thetas, .false., .false., &
            radii_sa, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, '65 params rejected with 1')
    call assert_close(maxval(abs(radii_sa)), 0.0_rk, 0.0_rk, 'too many params zero-fills radii')

    !---------------------------------------------------------------------------
    ! 3. No parameters at all: invalid init arguments.
    !---------------------------------------------------------------------------
    radii_sa(:) = 1.0_rk
    call compute_radius_grid_standalone_s(no_params, thetas, .false., .false., &
            radii_sa, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'empty params rejected with 5')
    call assert_close(maxval(abs(radii_sa)), 0.0_rk, 0.0_rk, 'empty params zero-fills radii')

    !---------------------------------------------------------------------------
    ! 4. Theta validation happens inside tables_init_s: one node is below the
    !    minimum grid size.
    !---------------------------------------------------------------------------
    one_theta(1)  = PI_C / 3.0_rk
    one_radius(1) = 1.0_rk
    call compute_radius_grid_standalone_s(params, one_theta, .false., .false., &
            one_radius, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'one-node theta set rejected with 3')
    call assert_close(maxval(abs(one_radius)), 0.0_rk, 0.0_rk, 'bad grid zero-fills radii')

    !---------------------------------------------------------------------------
    ! 5. The derivative variant checks both buffers; a failure zero-fills both.
    !---------------------------------------------------------------------------
    radii_sa(:) = 1.0_rk
    dr_bad(:)   = 1.0_rk
    call compute_radius_and_derivative_standalone_s(params, thetas, .false., .false., &
            radii_sa, dr_bad, status)
    call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            'wrong-length derivative buffer rejected with 104')
    call assert_close(maxval(abs(radii_sa)), 0.0_rk, 0.0_rk, 'bad buffer zero-fills radii')
    call assert_close(maxval(abs(dr_bad)), 0.0_rk, 0.0_rk, 'bad buffer zero-fills derivatives')

    call test_summary()

contains

    !> Exact equality on the bit patterns: -Wextra rejects `==` on reals, and
    !! a scalar transfer avoids the array temporary -Warray-temporaries flags.
    logical function bits_equal_f(a, b) result(ok)
        real(kind = rk), intent(in) :: a(:), b(:)
        integer(kind = ik) :: j
        ok = size(a, kind = ik) == size(b, kind = ik)
        if (.not. ok) return
        do j = 1_ik, size(a, kind = ik)
            if (transfer(a(j), 0_int64) /= transfer(b(j), 0_int64)) then
                ok = .false.
                return
            end if
        end do
    end function bits_equal_f
end program beta_param_standalone_test
