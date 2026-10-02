program beta_param_outputs_test
    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_legendre_polynomials_s, &
            compute_spherical_harmonics_normalization_constants_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary, &
            bits_equal_f
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, cache_resolve_shape_s, &
            node_set_t, node_set_build_s, node_set_free_s, &
            cache_node_radius_and_derivative_s, cache_radius_grid_unchecked_s, &
            SHAPE_VALID, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            BETA_PARAM_ERROR_NODE_SET_MISMATCH, BETA_PARAM_ERROR_COM_NOT_CONVERGED
    implicit none

    !> IEEE754 -0.0; built from bits so that -fno-signed-zeros cannot fold it
    !! into +0.0 the way a `-0.0_rk` literal could.
    integer(kind = ikl), parameter :: NEG_ZERO_BITS = int(z'8000000000000000', kind = ikl)

    type(cache_t) :: cache, cache_small, cache_reversed
    !> volatile: forces the -0.0 bit patterns to memory under -fno-signed-zeros.
    real(kind = rk), volatile :: negative_zeros(4)
    real(kind = rk) :: thetas_reversed(16), radii_reversed(16), drs_reversed(16)
    type(node_set_t) :: nodes, nodes_small, nodes_unbuilt
    integer(kind = ik) :: status, i
    ! bad/bad2 are two distinct undersized buffers: a call taking both radii and
    ! dr_dthetas needs separate actual arguments (two intent(out) dummies must
    ! not be aliased, F2018 15.5.2.13)
    real(kind = rk) :: thetas(16), params(4), radii(16), drs(16), bad(5), bad2(5)
    real(kind = rk) :: node_radii(16), node_drs(16), none(0), none2(0)
    real(kind = rk) :: norms(4), p(5), expected, b10, rn, rs, vf

    do i = 1_ik, 16_ik
        thetas(i) = real(i, rk) * PI_C / 17.0_rk
    end do
    call compute_spherical_harmonics_normalization_constants_s(norms, 4_ik)
    params = [0.0_rk, 0.2_rk, 0.0_rk, 0.0_rk]

    call cache_init_s(cache, 4_ik, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'cache init')

    ! no COM, no volume: pure table evaluation is exactly the Legendre sum
    call cache_radius_grid_s(cache, params, .false., .false., radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'radius grid ok')
    do i = 1_ik, 16_ik
        call compute_legendre_polynomials_s(5_ik, cos(thetas(i)), p)
        expected = 1.0_rk + params(2) * norms(2) * p(3)
        call assert_close(radii(i), expected, 1.0e-14_rk, 'R matches Legendre sum')
    end do

    call cache_radius_and_derivative_s(cache, params, .false., .false., radii, drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'radius+derivative ok')
    call assert_close(drs(8) + drs(9), 0.0_rk, 1.0e-12_rk, &
            'derivative antisymmetric for even-only shape')

    call cache_radius_grid_s(cache, params, .false., .false., bad(1:5), status)
    call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, 'bad buffer 104')
    call assert_close(maxval(abs(bad)), 0.0_rk, 0.0_rk, 'bad buffer zero-filled')

    ! volume conservation is a per-call option on the same cache: it scales the
    ! outputs by the factor resolve reports
    call cache_resolve_shape_s(cache, params, .true., .false., b10, rn, rs, vf, status)
    call cache_radius_grid_s(cache, params, .true., .false., radii, status)
    call compute_legendre_polynomials_s(5_ik, cos(thetas(1)), p)
    expected = (1.0_rk + params(2) * norms(2) * p(3)) * vf
    call assert_close(radii(1), expected, 1.0e-13_rk, 'volume factor applied to radii')
    call assert_true(abs(vf - 1.0_rk) > 1.0e-6_rk, 'volume factor differs from 1')

    ! a node set over the primary thetas reproduces the primary output exactly:
    ! same params, same tables, same eval kernel, so the bits must agree
    call node_set_build_s(nodes, cache, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'node set over the primary thetas')
    call cache_radius_and_derivative_s(cache, params, .true., .false., radii, drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'primary radius+derivative ok')
    call cache_node_radius_and_derivative_s(cache, params, nodes, .true., .false., &
            node_radii, node_drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'node radius+derivative ok')
    call assert_true(bits_equal_f(node_radii, radii), 'node radii bit-identical to primary')
    call assert_true(bits_equal_f(node_drs, drs), 'node derivatives bit-identical to primary')

    ! buffers are checked against the node count, before the node-set check
    call cache_node_radius_and_derivative_s(cache, params, nodes, .true., .false., &
            bad(1:5), bad2(1:5), status)
    call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            'node buffers sized to the node count')

    ! an unbuilt node set has zero nodes: only zero-length buffers reach 106
    call cache_node_radius_and_derivative_s(cache, params, nodes_unbuilt, .true., .false., &
            none, none2, status)
    call assert_int_eq(status, BETA_PARAM_ERROR_NODE_SET_MISMATCH, 'unbuilt node set 106')

    ! a node set built for a smaller cache cannot serve this one
    call cache_init_s(cache_small, 2_ik, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'small cache init')
    call node_set_build_s(nodes_small, cache_small, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'small node set build')
    call cache_node_radius_and_derivative_s(cache, params, nodes_small, .true., .false., &
            node_radii, node_drs, status)
    call assert_int_eq(status, BETA_PARAM_ERROR_NODE_SET_MISMATCH, 'node set mismatch 106')
    call assert_close(maxval(abs(node_radii)), 0.0_rk, 0.0_rk, 'mismatch zero-fills radii')
    call assert_close(maxval(abs(node_drs)), 0.0_rk, 0.0_rk, 'mismatch zero-fills derivatives')

    ! a node set built for a larger cache is accepted by a smaller one
    call cache_node_radius_and_derivative_s(cache_small, params(1:2), nodes, .true., .false., &
            node_radii, node_drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'larger node set serves a smaller cache')
    call assert_true(bits_equal_f(node_radii, radii), &
            'larger node set on a smaller cache: same bits as the primary')

    call node_set_free_s(nodes_small)
    call cache_free_s(cache_small)
    call node_set_free_s(nodes)

    ! the unchecked path renders a shape validation would reject
    call cache_radius_grid_unchecked_s(cache, [0.0_rk, -1.8_rk, 0.0_rk, 0.0_rk], .false., &
            radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'unchecked path skips validation')
    call assert_true(minval(radii) < 0.0_rk, 'broken outline survives unchecked')

    call cache_radius_grid_unchecked_s(cache, params, .false., bad(1:5), status)
    call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, 'unchecked bad buffer 104')

    ! with COM, a non-positive volume integral is a COM failure even unchecked
    call cache_radius_grid_unchecked_s(cache, [0.0_rk, -20.0_rk], .true., radii, status)
    call assert_int_eq(status, BETA_PARAM_ERROR_COM_NOT_CONVERGED, 'unchecked COM failure 103')
    call assert_close(maxval(abs(radii)), 0.0_rk, 0.0_rk, 'unchecked COM failure zero-fills')

    ! trailing zeros are trimmed: short and zero-padded forms give the same bits
    call cache_radius_and_derivative_s(cache, [0.05_rk, 0.2_rk], .true., .true., &
            radii, drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'short vector ok')
    call cache_radius_and_derivative_s(cache, [0.05_rk, 0.2_rk, 0.0_rk, 0.0_rk], &
            .true., .true., node_radii, node_drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'zero-padded vector ok')
    call assert_true(bits_equal_f(node_radii, radii), 'short == zero-padded radii, bitwise')
    call assert_true(bits_equal_f(node_drs, drs), 'short == zero-padded derivatives, bitwise')

    ! only TRAILING zeros are trimmed: an interior zero keeps its place, so
    ! (0, 0, 0.1) is the pure l = 3 shape, not the l = 1 shape
    call cache_radius_grid_s(cache, [0.0_rk, 0.0_rk, 0.1_rk], .false., .false., radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'interior-zero vector ok')
    do i = 1_ik, 16_ik
        call compute_legendre_polynomials_s(5_ik, cos(thetas(i)), p)
        expected = 1.0_rk + 0.1_rk * norms(3) * p(4)
        call assert_close(radii(i), expected, 1.0e-14_rk, 'interior zero: R is the l = 3 sum')
    end do

    ! an all-zero vector of any length is the unit sphere, with the same bits
    ! as the one-element zero vector — +0.0 or -0.0
    call cache_radius_and_derivative_s(cache, [0.0_rk], .true., .true., radii, drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'one-element zero vector ok')
    call assert_close(radii(1), 1.0_rk, 1.0e-12_rk, 'zero vector is the unit sphere')
    call cache_radius_and_derivative_s(cache, [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], &
            .true., .true., node_radii, node_drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'four-element zero vector ok')
    call assert_true(bits_equal_f(node_radii, radii) .and. bits_equal_f(node_drs, drs), &
            'all-zero vector: same bits as the one-element zero vector')
    negative_zeros(:) = transfer(NEG_ZERO_BITS, 1.0_rk)
    call cache_radius_and_derivative_s(cache, negative_zeros, .true., .true., &
            node_radii, node_drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'all -0.0 vector ok')
    call assert_true(bits_equal_f(node_radii, radii) .and. bits_equal_f(node_drs, drs), &
            'all -0.0 vector: same bits as the one-element zero vector')

    ! the theta set need not be sorted: a reversed set gives the same shape at
    ! the same angles (R reversed, dR/dtheta reversed)
    call cache_radius_and_derivative_s(cache, [0.05_rk, 0.2_rk, 0.1_rk], .true., .true., &
            radii, drs, status)
    call assert_int_eq(status, SHAPE_VALID, 'forward theta set ok')
    thetas_reversed(:) = thetas(16:1:-1)
    call cache_init_s(cache_reversed, 4_ik, thetas_reversed, status)
    call assert_int_eq(status, SHAPE_VALID, 'cache over a descending theta set')
    call cache_radius_and_derivative_s(cache_reversed, [0.05_rk, 0.2_rk, 0.1_rk], &
            .true., .true., radii_reversed, drs_reversed, status)
    call assert_int_eq(status, SHAPE_VALID, 'reversed theta set ok')
    do i = 1_ik, 16_ik
        call assert_close(radii_reversed(i), radii(17_ik - i), 1.0e-14_rk, &
                'reversed thetas: same R at the same angle')
        call assert_close(drs_reversed(i), drs(17_ik - i), 1.0e-14_rk, &
                'reversed thetas: same dR/dtheta at the same angle')
    end do
    call cache_free_s(cache_reversed)

    call cache_free_s(cache)
    call test_summary()

end program beta_param_outputs_test
