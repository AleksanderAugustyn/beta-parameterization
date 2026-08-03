!> Contract family 2: minimality. Every cached compute must recompute exactly
!! the intermediates its parameter diff invalidated and its `up_to` stage needs
!! — no more, no fewer. The engine's per-intermediate recompute counters are the
!! observable; this suite pins the whole counter table, checkpoint by checkpoint.
program beta_param_minimality_test
    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    use beta_parameterization_mod, only: cache_t, cache_free_s, &
            cache_init_shared_s, cache_resolve_shape_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_node_radius_and_derivative_s, cache_recompute_count_f, &
            tables_t, tables_init_s, tables_free_s, &
            node_set_t, node_set_build_s, node_set_free_s, &
            BETA_PARAM_I_RESOLVED, BETA_PARAM_I_MIN_RADIUS, BETA_PARAM_I_VOLUME, &
            BETA_PARAM_I_RADII, BETA_PARAM_I_DERIV, &
            SHAPE_VALID, SHAPE_ERROR_WRONG_PARAM_COUNT

    implicit none

    integer(kind = ik), parameter :: N_THETAS = 32_ik
    integer(kind = ik), parameter :: N_PARAMS = 6_ik
    integer(kind = ik), parameter :: N_NODES = 8_ik
    integer(kind = ik), parameter :: N_TRACKED = 5_ik
    real(kind = rk),    parameter :: DELTA = 0.013_rk

    !> IEEE754 -0.0; built from bits so that -fno-signed-zeros cannot fold it
    !! into +0.0 the way a `-0.0_rk` literal could.
    integer(kind = ikl), parameter :: NEG_ZERO_BITS = int(z'8000000000000000', kind = ikl)

    !> The five tracked intermediates, in dependency order.
    integer(kind = ik), parameter :: INTERMEDIATES(N_TRACKED) = &
            [BETA_PARAM_I_RESOLVED, BETA_PARAM_I_MIN_RADIUS, BETA_PARAM_I_VOLUME, &
             BETA_PARAM_I_RADII, BETA_PARAM_I_DERIV]
    character(len = 12), parameter :: INTERMEDIATE_NAMES(N_TRACKED) = &
            ['I_RESOLVED  ', 'I_MIN_RADIUS', 'I_VOLUME    ', 'I_RADII     ', 'I_DERIV     ']

    type(tables_t), target :: tables
    type(node_set_t) :: nodes
    real(kind = rk) :: thetas(N_THETAS), node_thetas(N_NODES)
    real(kind = rk) :: base_params(N_PARAMS)
    integer(kind = ik) :: i, status

    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do
    do i = 1_ik, N_NODES
        node_thetas(i) = real(2_ik * i - 1_ik, rk) * PI_C / real(2_ik * N_NODES, rk)
    end do
    base_params = [0.02_rk, 0.15_rk, 0.08_rk, -0.05_rk, 0.03_rk, 0.01_rk]

    ! The whole suite runs over one shared tables_t: the node-set path needs a
    ! tables_t the test can reach, and the counter semantics are identical to a
    ! private cache's.
    call tables_init_s(tables, N_PARAMS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'shared tables init')
    call node_set_build_s(nodes, tables, node_thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'node set build')

    call run_index_mirror_s()
    call run_counter_table_s()
    call run_cold_partial_stages_s()
    call run_negzero_counter_s()

    call node_set_free_s(nodes)
    call tables_free_s(tables)

    call test_summary()

contains

    !> The public mirrors must carry the values the library uses internally;
    !! everything below indexes the counters through them. A never-initialized
    !! cache answers 0 rather than reading an unconfigured engine.
    subroutine run_index_mirror_s()
        type(cache_t) :: fresh
        integer(kind = ik) :: want(N_TRACKED)
        want = 0_ik
        call check_counters_s(fresh, want, 'no init at all')
        call assert_int_eq(BETA_PARAM_I_RESOLVED, 1_ik, 'BETA_PARAM_I_RESOLVED == 1')
        call assert_int_eq(BETA_PARAM_I_MIN_RADIUS, 2_ik, 'BETA_PARAM_I_MIN_RADIUS == 2')
        call assert_int_eq(BETA_PARAM_I_VOLUME, 3_ik, 'BETA_PARAM_I_VOLUME == 3')
        call assert_int_eq(BETA_PARAM_I_RADII, 4_ik, 'BETA_PARAM_I_RADII == 4')
        call assert_int_eq(BETA_PARAM_I_DERIV, 5_ik, 'BETA_PARAM_I_DERIV == 5')
    end subroutine run_index_mirror_s

    !> The contract's counter table, in order, on one cache.
    subroutine run_counter_table_s()
        type(cache_t) :: cache
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: radii(N_THETAS), drs(N_THETAS)
        real(kind = rk) :: node_radii(N_NODES), node_drs(N_NODES)
        real(kind = rk) :: b10, rn, rs, vf
        integer(kind = ik) :: want(N_TRACKED)
        integer(kind = ik) :: j, status
        character(len = 32) :: label

        call cache_init_shared_s(cache, tables, N_PARAMS, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'shared cache init')

        ! 1. Cold radius grid: intermediates 1-4 only, the derivative stays cold.
        params = base_params
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'cold radius grid status')
        want = [1_ik, 1_ik, 1_ik, 1_ik, 0_ik]
        call check_counters_s(cache, want, 'cold radius grid')

        ! 2. Same vector again: nothing is recomputed.
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'repeat radius grid status')
        want = [1_ik, 1_ik, 1_ik, 1_ik, 0_ik]
        call check_counters_s(cache, want, 'repeat radius grid')

        ! 3. Same vector, derivative requested: only I_DERIV is computed.
        call cache_radius_and_derivative_s(cache, params, radii, drs, status)
        call assert_int_eq(status, SHAPE_VALID, 'warm derivative status')
        want = [1_ik, 1_ik, 1_ik, 1_ik, 1_ik]
        call check_counters_s(cache, want, 'warm derivative')
        call assert_true(maxval(abs(drs)) > 0.0_rk, &
                'derivatives are not identically zero, warm derivative')

        ! 4. Same vector, resolve: every stamp it needs is already valid.
        call cache_resolve_shape_s(cache, params, b10, rn, rs, vf, status)
        call assert_int_eq(status, SHAPE_VALID, 'warm resolve status')
        want = [1_ik, 1_ik, 1_ik, 1_ik, 1_ik]
        call check_counters_s(cache, want, 'warm resolve')

        ! 5. Every parameter in turn: all masks are full, so each perturbation
        !    invalidates all five intermediates and a derivative call recomputes
        !    all five. Checking one index only would let a wrong mask on any
        !    other parameter pass.
        do j = 1_ik, N_PARAMS
            write(label, '(A,I0)') 'perturbation of param ', j
            params(j) = params(j) + DELTA
            call cache_radius_and_derivative_s(cache, params, radii, drs, status)
            call assert_int_eq(status, SHAPE_VALID, 'status after ' // trim(label))
            want = 1_ik + j
            call check_counters_s(cache, want, trim(label))
            call assert_true(maxval(abs(drs)) > 0.0_rk, &
                    'derivatives are not identically zero after ' // trim(label))
        end do

        ! 6. The node path is uncached: intermediates 1-3 are already valid on
        !    this vector and the node evaluation itself is tracked by nothing.
        call cache_node_radius_and_derivative_s(cache, nodes, params, node_radii, &
                node_drs, status)
        call assert_int_eq(status, SHAPE_VALID, 'node call status')
        want = 1_ik + N_PARAMS
        call check_counters_s(cache, want, 'node call')
        call assert_true(maxval(abs(node_drs)) > 0.0_rk, &
                'node derivatives are not identically zero')

        ! 7. Precedence: a wrong parameter count outranks a wrong buffer size
        !    (status 4, not 104), and a rejected call recomputes nothing.
        call cache_radius_grid_s(cache, params(1:N_PARAMS - 1_ik), &
                radii(1:N_THETAS - 1_ik), status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, &
                'wrong param count outranks wrong buffer size')
        want = 1_ik + N_PARAMS
        call check_counters_s(cache, want, 'rejected wrong-count call')

        ! Out-of-range indices are answered, not trapped.
        call assert_int_eq(int(cache_recompute_count_f(cache, 0_ik), ik), 0_ik, &
                'counter of intermediate 0 is 0')
        call assert_int_eq(int(cache_recompute_count_f(cache, 99_ik), ik), 0_ik, &
                'counter of an out-of-range intermediate is 0')

        call cache_free_s(cache)

        ! A freed cache is uninitialized again: every counter reads 0.
        want = 0_ik
        call check_counters_s(cache, want, 'cache_free_s')
    end subroutine run_counter_table_s

    !> The two entry points that stop at I_VOLUME, each on a FRESH cache. Warm
    !! checkpoints cannot see over-computation here: on a cache where all five
    !! stamps are already valid, an ensure call raised to I_DERIV would recompute
    !! nothing and the counters would not move. Cold is the only state that
    !! shows where each entry point actually stops.
    subroutine run_cold_partial_stages_s()
        type(cache_t) :: resolve_cache, node_cache
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: node_radii(N_NODES), node_drs(N_NODES)
        real(kind = rk) :: b10, rn, rs, vf
        integer(kind = ik) :: want(N_TRACKED)
        integer(kind = ik) :: status

        params = base_params
        want = [1_ik, 1_ik, 1_ik, 0_ik, 0_ik]

        ! Cold resolve: intermediates 1-3 only; the radius and derivative tables
        ! stay cold because resolve never asks for them.
        call cache_init_shared_s(resolve_cache, tables, N_PARAMS, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'cold resolve cache init')
        call cache_resolve_shape_s(resolve_cache, params, b10, rn, rs, vf, status)
        call assert_int_eq(status, SHAPE_VALID, 'cold resolve status')
        call check_counters_s(resolve_cache, want, 'cold resolve')
        call cache_free_s(resolve_cache)

        ! Cold node call: the node evaluation is uncached, so it shares
        ! intermediates 1-3 with the grid path and computes nothing else.
        call cache_init_shared_s(node_cache, tables, N_PARAMS, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'cold node cache init')
        call cache_node_radius_and_derivative_s(node_cache, nodes, params, node_radii, &
                node_drs, status)
        call assert_int_eq(status, SHAPE_VALID, 'cold node status')
        call check_counters_s(node_cache, want, 'cold node call')
        call assert_true(maxval(abs(node_drs)) > 0.0_rk, &
                'cold node derivatives are not identically zero')
        call cache_free_s(node_cache)
    end subroutine run_cold_partial_stages_s

    !> Counter half of the -0.0 case: the parameter diff is on bit patterns, so
    !! replacing +0.0 by -0.0 is a change even though the two compare equal.
    subroutine run_negzero_counter_s()
        type(cache_t) :: cache
        !> volatile: under -fno-signed-zeros (Release -ffast-math) the compiler
        !! may elide a store of -0.0 over a location it knows holds +0.0, so the
        !! bit pattern would never reach the engine; volatile forces the store.
        real(kind = rk), volatile :: pv(N_PARAMS)
        real(kind = rk) :: radii(N_THETAS), drs(N_THETAS)
        integer(kind = ik) :: want(N_TRACKED)
        integer(kind = ik) :: k, status

        call cache_init_shared_s(cache, tables, N_PARAMS, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'negative-zero cache init')

        do k = 1_ik, N_PARAMS
            pv(k) = base_params(k)
        end do
        pv(1) = 0.0_rk
        call cache_radius_and_derivative_s(cache, pv, radii, drs, status)
        call assert_int_eq(status, SHAPE_VALID, 'negative-zero warm-up status')
        want = [1_ik, 1_ik, 1_ik, 1_ik, 1_ik]
        call check_counters_s(cache, want, 'negative-zero warm-up')

        pv(1) = transfer(NEG_ZERO_BITS, 1.0_rk)
        call assert_true(transfer(pv(1), 0_ikl) == NEG_ZERO_BITS, &
                'the -0.0 store reached memory')
        call cache_radius_and_derivative_s(cache, pv, radii, drs, status)
        call assert_int_eq(status, SHAPE_VALID, 'negative-zero incremental status')
        want = [2_ik, 2_ik, 2_ik, 2_ik, 2_ik]
        call check_counters_s(cache, want, 'negative-zero incremental call')

        call cache_free_s(cache)
    end subroutine run_negzero_counter_s

    !> Assert all five recompute counters at one checkpoint.
    subroutine check_counters_s(cache, expected, label)
        type(cache_t),      intent(in) :: cache
        integer(kind = ik), intent(in) :: expected(N_TRACKED)
        character(len = *), intent(in) :: label
        integer(kind = ik) :: k
        do k = 1_ik, N_TRACKED
            call assert_int_eq(int(cache_recompute_count_f(cache, INTERMEDIATES(k)), ik), &
                    expected(k), trim(INTERMEDIATE_NAMES(k)) // ' count after ' // label)
        end do
    end subroutine check_counters_s

end program beta_param_minimality_test
