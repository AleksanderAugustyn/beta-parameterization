program beta_param_outputs_test
    use, intrinsic :: iso_fortran_env, only: int64
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_legendre_polynomials_s, &
            compute_spherical_harmonics_normalization_constants_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, cache_resolve_shape_s, &
            tables_t, tables_init_s, tables_free_s, &
            node_set_t, node_set_build_s, node_set_free_s, &
            cache_init_shared_s, cache_node_radius_and_derivative_s, &
            cache_radius_grid_unchecked_s, &
            SHAPE_VALID, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            BETA_PARAM_ERROR_NODE_SET_MISMATCH
    implicit none

    type(cache_t) :: cache
    type(tables_t), target :: tables
    type(tables_t) :: tables_small
    type(node_set_t) :: nodes, nodes_small
    integer(kind = ik) :: status, i
    ! bad/bad2 are two distinct undersized buffers: a call taking both radii and
    ! dr_dthetas needs separate actual arguments (two intent(out) dummies must
    ! not be aliased, F2018 15.5.2.13)
    real(kind = rk) :: thetas(16), params(4), radii(16), drs(16), bad(5), bad2(5)
    real(kind = rk) :: node_radii(16), node_drs(16)
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

    ! a node set over the primary thetas reproduces the primary output exactly:
    ! same params, same tables, same eval kernel, so the bits must agree
    call tables_init_s(tables, 4_ik, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'shared tables init')
    call node_set_build_s(nodes, tables, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'node set over the primary thetas')
    call cache_init_shared_s(cache, tables, 4_ik, .true., .false., status)
    call assert_int_eq(status, SHAPE_VALID, 'shared cache init')

    call cache_radius_and_derivative_s(cache, params, radii, drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'primary radius+derivative ok')
    call cache_node_radius_and_derivative_s(cache, nodes, params, node_radii, &
            node_drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'node radius+derivative ok')
    call assert_true(bits_equal_f(node_radii, radii), 'node radii bit-identical to primary')
    call assert_true(bits_equal_f(node_drs, drs), 'node derivatives bit-identical to primary')

    ! buffers are checked against the node count, before the node-set check
    call cache_node_radius_and_derivative_s(cache, nodes, params, bad(1:5), &
            bad2(1:5), status)
    call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            'node buffers sized to the node count')

    ! a node set whose max_l is below the cache parameter count is rejected
    call tables_init_s(tables_small, 2_ik, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'small tables init')
    call node_set_build_s(nodes_small, tables_small, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'small node set build')
    call cache_node_radius_and_derivative_s(cache, nodes_small, params, node_radii, &
            node_drs, status)
    call assert_int_eq(status, BETA_PARAM_ERROR_NODE_SET_MISMATCH, 'node set mismatch 106')
    call assert_close(maxval(abs(node_radii)), 0.0_rk, 0.0_rk, 'mismatch zero-fills radii')
    call assert_close(maxval(abs(node_drs)), 0.0_rk, 0.0_rk, 'mismatch zero-fills derivatives')

    call node_set_free_s(nodes_small)
    call tables_free_s(tables_small)
    call node_set_free_s(nodes)
    call cache_free_s(cache)
    call tables_free_s(tables)

    ! the unchecked path renders a shape validation would reject
    call cache_init_s(cache, 4_ik, thetas, .false., .false., status)
    call cache_radius_grid_unchecked_s(cache, [0.0_rk, -1.8_rk, 0.0_rk, 0.0_rk], &
            radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'unchecked path skips validation')
    call assert_true(minval(radii) < 0.0_rk, 'broken outline survives unchecked')

    call cache_radius_grid_unchecked_s(cache, params, bad(1:5), status)
    call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, 'unchecked bad buffer 104')
    call cache_free_s(cache)
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
end program beta_param_outputs_test
