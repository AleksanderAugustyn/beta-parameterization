!> Contract family 4: boundaries. The two parameter caps (cached tier 2 and
!! standalone tier 1) and the theta-grid floor must accept the last legal value
!! and reject the first illegal one with the documented code — off-by-one on
!! either side is a contract break, so both sides of every edge are asserted.
program beta_param_boundary_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_int_eq, test_summary
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, &
            tables_t, tables_init_s, tables_free_s, &
            compute_radius_grid_standalone_s, &
            SHAPE_CACHE_MAX_PARAMS, SHAPE_STANDALONE_MAX_PARAMS, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_INVALID_GRID, &
            BETA_PARAM_ERROR_POLE_NODE

    implicit none

    integer(kind = ik), parameter :: N_THETAS = 32_ik
    integer(kind = ik), parameter :: N_PARAMS_SMALL = 2_ik

    !> A valid eight-parameter shape, used to show the cap value is not merely
    !! accepted by the init but actually usable.
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
    integer(kind = ik) :: i

    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do

    call run_cap_mirrors_s()
    call run_cache_cap_s()
    call run_standalone_cap_s()
    call run_grid_floor_s()
    call run_pole_thetas_s()

    call test_summary()

contains

    !> The two caps are part of the published contract; everything below is
    !! written against these values.
    subroutine run_cap_mirrors_s()
        call assert_int_eq(SHAPE_CACHE_MAX_PARAMS, 8_ik, 'SHAPE_CACHE_MAX_PARAMS == 8')
        call assert_int_eq(SHAPE_STANDALONE_MAX_PARAMS, 64_ik, 'SHAPE_STANDALONE_MAX_PARAMS == 64')
    end subroutine run_cap_mirrors_s

    !> Tier 2 accepts exactly SHAPE_CACHE_MAX_PARAMS parameters and rejects one
    !! more with SHAPE_ERROR_TOO_MANY_PARAMS.
    subroutine run_cache_cap_s()
        type(cache_t) :: at_cap, over_cap
        real(kind = rk) :: radii(N_THETAS)
        integer(kind = ik) :: status

        call cache_init_s(at_cap, SHAPE_CACHE_MAX_PARAMS, thetas, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'cache init at the parameter cap')
        ! The cap is usable, not just accepted.
        call cache_radius_grid_s(at_cap, PARAMS_8, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'compute on a cache at the parameter cap')
        call cache_free_s(at_cap)

        call cache_init_s(over_cap, SHAPE_CACHE_MAX_PARAMS + 1_ik, thetas, &
                .true., .true., status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, &
                'cache init one parameter over the cap')
    end subroutine run_cache_cap_s

    !> Tier 1 accepts exactly SHAPE_STANDALONE_MAX_PARAMS parameters — eight
    !! times the cached cap — and rejects one more the same way.
    subroutine run_standalone_cap_s()
        real(kind = rk) :: params_at(SHAPE_STANDALONE_MAX_PARAMS)
        real(kind = rk) :: params_over(SHAPE_STANDALONE_MAX_PARAMS + 1_ik)
        real(kind = rk) :: radii(N_THETAS)
        integer(kind = ik) :: status

        params_at(:) = 0.0_rk
        params_at(2) = 0.1_rk
        call compute_radius_grid_standalone_s(params_at, thetas, .true., .true., &
                radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'standalone at the parameter cap')

        params_over(:) = 0.0_rk
        params_over(2) = 0.1_rk
        call compute_radius_grid_standalone_s(params_over, thetas, .true., .true., &
                radii, status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, &
                'standalone one parameter over the cap')
    end subroutine run_standalone_cap_s

    !> Two theta nodes are the fewest the Legendre tables can be built on; one is
    !! rejected with SHAPE_ERROR_INVALID_GRID. Both levels enforce the same
    !! floor, because a cache validates its grid through tables_init_s.
    subroutine run_grid_floor_s()
        type(tables_t) :: tables
        type(cache_t) :: cache
        real(kind = rk) :: thetas_2(2), thetas_1(1)
        integer(kind = ik) :: status

        thetas_2 = [0.7_rk, 2.0_rk]
        thetas_1 = [0.7_rk]

        call tables_init_s(tables, N_PARAMS_SMALL, thetas_2, status)
        call assert_int_eq(status, SHAPE_VALID, 'tables init on two thetas')
        call tables_free_s(tables)

        call cache_init_s(cache, N_PARAMS_SMALL, thetas_2, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'cache init on two thetas')
        call cache_free_s(cache)

        call tables_init_s(tables, N_PARAMS_SMALL, thetas_1, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'tables init on one theta')

        call cache_init_s(cache, N_PARAMS_SMALL, thetas_1, .true., .true., status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'cache init on one theta')
    end subroutine run_grid_floor_s

    !> A polar node is rejected with BETA_PARAM_ERROR_POLE_NODE at both levels,
    !! including the node that is polar only after rounding.
    subroutine run_pole_thetas_s()
        type(tables_t) :: tables
        type(cache_t) :: cache
        integer(kind = ik) :: k, status

        do k = 1_ik, 3_ik
            call tables_init_s(tables, N_PARAMS_SMALL, POLE_THETAS(:, k), status)
            call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, &
                    'tables init with ' // trim(POLE_NAMES(k)))

            call cache_init_s(cache, N_PARAMS_SMALL, POLE_THETAS(:, k), &
                    .true., .true., status)
            call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, &
                    'cache init with ' // trim(POLE_NAMES(k)))
        end do
    end subroutine run_pole_thetas_s

end program beta_param_boundary_test
