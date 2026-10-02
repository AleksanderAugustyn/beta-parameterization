!> Contract family 4: boundaries. The parameter limit L = 64 at init, in a
!! cached call and in a one-shot call, and the theta-grid floor, must accept
!! the last legal value and reject the first illegal one with the documented
!! code — off-by-one on either side is a contract break, so both sides of every
!! edge are asserted.
program beta_param_boundary_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary, bits_equal_f
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_node_radius_and_derivative_s, &
            node_set_t, node_set_build_s, node_set_free_s, &
            compute_radius_grid_standalone_s, &
            SHAPE_MAX_PARAMS, MAX_BETA_PARAMS_LIMIT, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            BETA_PARAM_ERROR_POLE_NODE, BETA_PARAM_ERROR_NODE_SET_MISMATCH

    implicit none

    integer(kind = ik), parameter :: N_THETAS = 32_ik
    integer(kind = ik), parameter :: N_PARAMS_SMALL = 2_ik
    !> The contract's L = min(SHAPE_MAX_PARAMS, N_max).
    integer(kind = ik), parameter :: LIMIT = min(SHAPE_MAX_PARAMS, MAX_BETA_PARAMS_LIMIT)

    !> A valid eight-parameter shape.
    real(kind = rk), parameter :: PARAMS_8(8) = &
            [0.02_rk, 0.15_rk, 0.08_rk, -0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, -0.002_rk]

    !> The three ways a theta node can sit on a pole: exactly north, exactly
    !! south, and close enough to north that cos(theta) rounds to 1. The
    !! derivative recurrence divides by 1 - cos(theta)**2, so all three are
    !! rejected the same way.
    real(kind = rk), parameter :: POLE_THETAS(2, 3) = reshape( &
            [0.0_rk, 1.0_rk, PI_C, 1.0_rk, 1.0e-9_rk, 1.0_rk], [2, 3])
    character(len = 12), parameter :: POLE_NAMES(3) = &
            ['theta = 0   ', 'theta = pi  ', 'theta = 1e-9']

    real(kind = rk) :: thetas(N_THETAS)
    real(kind = rk) :: params_at(LIMIT), params_over(LIMIT + 1_ik)
    integer(kind = ik) :: i

    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do
    ! Alternating, decaying amplitudes: a valid shape at the full length.
    do i = 1_ik, LIMIT
        params_at(i) = 0.02_rk / real(i, rk)
        if (mod(i, 2_ik) == 0_ik) params_at(i) = -params_at(i)
    end do
    params_over(1:LIMIT)     = params_at
    params_over(LIMIT + 1_ik) = 1.0e-4_rk

    call run_limit_value_s()
    call run_init_limit_s()
    call run_cached_call_limit_s()
    call run_one_shot_limit_s()
    call run_grid_floor_s()
    call run_pole_thetas_s()
    call run_node_set_across_caches_s()

    call test_summary()

contains

    !> The limit is part of the published contract; everything below is written
    !! against this value.
    subroutine run_limit_value_s()
        call assert_int_eq(LIMIT, 64_ik, 'L == 64')
    end subroutine run_limit_value_s

    !> Init accepts max_params = L, rejects L + 1 with 1 and 0 with 5.
    subroutine run_init_limit_s()
        type(cache_t) :: cache
        integer(kind = ik) :: status

        call cache_init_s(cache, LIMIT, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'init at max_params = L')
        call cache_free_s(cache)

        call cache_init_s(cache, LIMIT + 1_ik, thetas, status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'init at L + 1 rejected with 1')

        call cache_init_s(cache, 0_ik, thetas, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'init at 0 rejected with 5')
    end subroutine run_init_limit_s

    !> A cached call accepts size(params) = max_params and rejects one more, or
    !! none, with 4. Checked on a small cache and on one at the limit, where
    !! perturbing the last parameter must change the output (no truncation).
    subroutine run_cached_call_limit_s()
        type(cache_t) :: cache
        real(kind = rk) :: radii(N_THETAS), radii_perturbed(N_THETAS)
        real(kind = rk) :: no_params(0), perturbed(LIMIT)
        integer(kind = ik) :: status

        call cache_init_s(cache, 8_ik, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'cache init at max_params = 8')
        call cache_radius_grid_s(cache, PARAMS_8, .true., .true., radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'call with size(params) = max_params')
        call cache_radius_grid_s(cache, params_at(1:9), .true., .true., radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, &
                'call with max_params + 1 rejected with 4')
        call cache_radius_grid_s(cache, no_params, .true., .true., radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, &
                'call with an empty vector rejected with 4')
        call cache_free_s(cache)

        call cache_init_s(cache, LIMIT, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'cache init at max_params = L')
        call cache_radius_grid_s(cache, params_at, .true., .true., radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'cached call with L parameters')
        perturbed = params_at
        perturbed(LIMIT) = perturbed(LIMIT) + 1.0e-3_rk
        call cache_radius_grid_s(cache, perturbed, .true., .true., radii_perturbed, status)
        call assert_int_eq(status, SHAPE_VALID, 'cached call with the last parameter perturbed')
        call assert_true(.not. bits_equal_f(radii, radii_perturbed), &
                'cached: parameter L changes the output')
        call cache_radius_grid_s(cache, params_over, .true., .true., radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, &
                'cached call with L + 1 parameters rejected with 4')
        call cache_free_s(cache)
    end subroutine run_cached_call_limit_s

    !> A one-shot call accepts L parameters, rejects L + 1 with 1 and an empty
    !! vector with 4; at L the last parameter changes the output.
    subroutine run_one_shot_limit_s()
        real(kind = rk) :: radii(N_THETAS), radii_perturbed(N_THETAS)
        real(kind = rk) :: no_params(0), perturbed(LIMIT)
        integer(kind = ik) :: status

        call compute_radius_grid_standalone_s(params_at, thetas, .true., .true., &
                radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'one-shot with L parameters')

        perturbed = params_at
        perturbed(LIMIT) = perturbed(LIMIT) + 1.0e-3_rk
        call compute_radius_grid_standalone_s(perturbed, thetas, .true., .true., &
                radii_perturbed, status)
        call assert_int_eq(status, SHAPE_VALID, 'one-shot with the last parameter perturbed')
        call assert_true(.not. bits_equal_f(radii, radii_perturbed), &
                'one-shot: parameter L changes the output')

        call compute_radius_grid_standalone_s(params_over, thetas, .true., .true., &
                radii, status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, &
                'one-shot with L + 1 parameters rejected with 1')

        call compute_radius_grid_standalone_s(no_params, thetas, .true., .true., &
                radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, &
                'one-shot with an empty vector rejected with 4')
    end subroutine run_one_shot_limit_s

    !> Two theta nodes are the fewest the Legendre tables can be built on; one is
    !! rejected with SHAPE_ERROR_INVALID_GRID at init, at node-set build and in
    !! the one-shot.
    subroutine run_grid_floor_s()
        type(cache_t) :: cache
        type(node_set_t) :: nodes
        real(kind = rk) :: thetas_2(2), thetas_1(1), radii_2(2), radii_1(1)
        integer(kind = ik) :: status

        thetas_2 = [0.7_rk, 2.0_rk]
        thetas_1 = [0.7_rk]

        call cache_init_s(cache, N_PARAMS_SMALL, thetas_1, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'cache init on one theta')

        call cache_init_s(cache, N_PARAMS_SMALL, thetas_2, status)
        call assert_int_eq(status, SHAPE_VALID, 'cache init on two thetas')

        call node_set_build_s(nodes, cache, thetas_2, status)
        call assert_int_eq(status, SHAPE_VALID, 'node set on two thetas')
        call node_set_build_s(nodes, cache, thetas_1, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'node set on one theta')
        call cache_free_s(cache)

        call compute_radius_grid_standalone_s(PARAMS_8(1:2), thetas_2, .true., .true., &
                radii_2, status)
        call assert_int_eq(status, SHAPE_VALID, 'one-shot on two thetas')
        call compute_radius_grid_standalone_s(PARAMS_8(1:2), thetas_1, .true., .true., &
                radii_1, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'one-shot on one theta')
    end subroutine run_grid_floor_s

    !> A polar node is rejected with BETA_PARAM_ERROR_POLE_NODE at init and at
    !! node-set build, including the node that is polar only after rounding.
    subroutine run_pole_thetas_s()
        type(cache_t) :: cache, good_cache
        type(node_set_t) :: nodes
        integer(kind = ik) :: k, status

        call cache_init_s(good_cache, N_PARAMS_SMALL, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'pole-test cache init')

        do k = 1_ik, 3_ik
            call cache_init_s(cache, N_PARAMS_SMALL, POLE_THETAS(:, k), status)
            call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, &
                    'cache init with ' // trim(POLE_NAMES(k)))

            call node_set_build_s(nodes, good_cache, POLE_THETAS(:, k), status)
            call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, &
                    'node set with ' // trim(POLE_NAMES(k)))
        end do
        call cache_free_s(good_cache)
    end subroutine run_pole_thetas_s

    !> A node set serves a cache only if it covers every vector that cache
    !! accepts. The verdict depends on the two objects, never on the vector:
    !! a short vector and its zero-padded forms get the same status.
    subroutine run_node_set_across_caches_s()
        type(cache_t) :: cache_3, cache_8
        type(node_set_t) :: nodes_3, nodes_8
        real(kind = rk) :: radii(N_THETAS), drs(N_THETAS)
        real(kind = rk) :: padded(8)
        integer(kind = ik) :: status

        padded(:)   = 0.0_rk
        padded(1:3) = PARAMS_8(1:3)

        call cache_init_s(cache_3, 3_ik, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'cache_3 init')
        call cache_init_s(cache_8, 8_ik, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'cache_8 init')
        call node_set_build_s(nodes_3, cache_3, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'nodes_3 build')
        call node_set_build_s(nodes_8, cache_8, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'nodes_8 build')

        ! small node set through the large cache: 106 whatever the vector length
        call cache_node_radius_and_derivative_s(cache_8, padded(1:3), nodes_3, &
                .true., .true., radii, drs, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_NODE_SET_MISMATCH, &
                'small node set, length-3 vector')
        call cache_node_radius_and_derivative_s(cache_8, padded(1:5), nodes_3, &
                .true., .true., radii, drs, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_NODE_SET_MISMATCH, &
                'small node set, vector padded to 5')
        call cache_node_radius_and_derivative_s(cache_8, padded, nodes_3, &
                .true., .true., radii, drs, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_NODE_SET_MISMATCH, &
                'small node set, vector padded to 8')

        ! large node set through the small cache: accepted
        call cache_node_radius_and_derivative_s(cache_3, padded(1:3), nodes_8, &
                .true., .true., radii, drs, status)
        call assert_int_eq(status, SHAPE_VALID, 'large node set through the small cache')

        call node_set_free_s(nodes_3)
        call node_set_free_s(nodes_8)
        call cache_free_s(cache_3)
        call cache_free_s(cache_8)
    end subroutine run_node_set_across_caches_s

end program beta_param_boundary_test
