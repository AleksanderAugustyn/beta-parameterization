!> Contract family 3: interleaving and recovery. A cache serves calls in any
!! order, on any sequence of parameter vectors, and a rejected call leaves no
!! trace — every output must match what a fresh cache would have produced from
!! that single call alone. All comparisons here are on IEEE754 bit patterns; no
!! tolerances appear anywhere in this suite.
program beta_param_interleave_test
    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    use beta_parameterization_mod, only: cache_t, cache_free_s, &
            cache_init_shared_s, cache_resolve_shape_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_node_radius_and_derivative_s, cache_radius_grid_unchecked_s, &
            cache_recompute_count_f, &
            tables_t, tables_init_s, tables_free_s, &
            node_set_t, node_set_build_s, node_set_free_s, &
            BETA_PARAM_I_RESOLVED, BETA_PARAM_I_MIN_RADIUS, BETA_PARAM_I_VOLUME, &
            BETA_PARAM_I_RADII, BETA_PARAM_I_DERIV, &
            BETA_PARAM_ERROR_NORTH_POLE, &
            SHAPE_VALID, SHAPE_ERROR_WRONG_PARAM_COUNT

    implicit none

    integer(kind = ik), parameter :: N_THETAS = 32_ik
    integer(kind = ik), parameter :: N_PARAMS = 4_ik
    integer(kind = ik), parameter :: N_NODES = 8_ik
    integer(kind = ik), parameter :: N_SCALARS = 4_ik
    integer(kind = ik), parameter :: N_STEPS = 8_ik
    integer(kind = ik), parameter :: N_TRACKED = 5_ik

    !> Two distinct valid vectors. A is the asymmetric pick the other contract
    !! suites use; B differs in every component so no diff mask is trivially empty.
    real(kind = rk), parameter :: PARAMS_A(N_PARAMS) = [0.05_rk, 0.20_rk, 0.10_rk, 0.02_rk]
    real(kind = rk), parameter :: PARAMS_B(N_PARAMS) = [-0.03_rk, 0.25_rk, 0.05_rk, -0.01_rk]

    !> beta2 = -1.8 drives R(theta = 0) negative, so the validity scan rejects it
    !! with BETA_PARAM_ERROR_NORTH_POLE. (A merely large NEGATIVE beta such as
    !! -0.99 is not invalid — the shape stays positive everywhere.)
    real(kind = rk), parameter :: PARAMS_INVALID(N_PARAMS) = &
            [0.0_rk, -1.8_rk, 0.0_rk, 0.0_rk]

    !> Which entry point each step calls, and which vector it calls it on. Every
    !! entry point is exercised on both vectors, and steps 4 and 5 repeat vector
    !! B across an entry-point change.
    integer(kind = ik), parameter :: STEP_ROUTINE(N_STEPS) = [1_ik, 2_ik, 3_ik, 4_ik, 1_ik, 2_ik, 3_ik, 4_ik]
    integer(kind = ik), parameter :: STEP_VECTOR(N_STEPS)  = [1_ik, 2_ik, 1_ik, 2_ik, 2_ik, 1_ik, 2_ik, 1_ik]

    character(len = 26), parameter :: ROUTINE_NAMES(4) = &
            ['radius_grid               ', 'radius_and_derivative     ', &
             'resolve_shape             ', 'node_radius_and_derivative']
    character(len = 1), parameter :: VECTOR_NAMES(2) = ['A', 'B']
    character(len = 4), parameter :: SCALAR_NAMES(N_SCALARS) = ['b10 ', 'rn  ', 'rs  ', 'vf  ']

    !> The five tracked intermediates, in dependency order.
    integer(kind = ik), parameter :: INTERMEDIATES(N_TRACKED) = &
            [BETA_PARAM_I_RESOLVED, BETA_PARAM_I_MIN_RADIUS, BETA_PARAM_I_VOLUME, &
             BETA_PARAM_I_RADII, BETA_PARAM_I_DERIV]
    character(len = 12), parameter :: INTERMEDIATE_NAMES(N_TRACKED) = &
            ['I_RESOLVED  ', 'I_MIN_RADIUS', 'I_VOLUME    ', 'I_RADII     ', 'I_DERIV     ']

    ! One shared tables_t for the whole suite: the node path needs a tables_t the
    ! test can reach, and every cache below is a fresh cache over it.
    type(tables_t), target :: tables
    type(node_set_t) :: nodes
    real(kind = rk) :: thetas(N_THETAS), node_thetas(N_NODES)
    real(kind = rk) :: params_ab(N_PARAMS, 2)

    ! Cold references, one column per vector. Each was produced by a fresh cache
    ! making exactly one call.
    real(kind = rk) :: cold_grid(N_THETAS, 2)
    real(kind = rk) :: cold_radii(N_THETAS, 2), cold_drs(N_THETAS, 2)
    real(kind = rk) :: cold_scalars(N_SCALARS, 2)
    real(kind = rk) :: cold_node_radii(N_NODES, 2), cold_node_drs(N_NODES, 2)

    integer(kind = ik) :: i, status

    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do
    do i = 1_ik, N_NODES
        node_thetas(i) = real(2_ik * i - 1_ik, rk) * PI_C / real(2_ik * N_NODES, rk)
    end do
    params_ab(:, 1) = PARAMS_A
    params_ab(:, 2) = PARAMS_B

    call tables_init_s(tables, N_PARAMS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'shared tables init')
    call node_set_build_s(nodes, tables, node_thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'node set build')

    call build_cold_references_s(1_ik)
    call build_cold_references_s(2_ik)

    call run_reference_distinctness_s()
    call run_alternation_s()
    call run_invalid_recovery_s()
    call run_wrong_length_recovery_s()
    call run_unchecked_then_checked_s()
    call run_node_before_primary_s()

    call node_set_free_s(nodes)
    call tables_free_s(tables)

    call test_summary()

contains

    !> The cold half of every comparison below: four fresh caches, one call each,
    !! so no reference is contaminated by another entry point's state.
    subroutine build_cold_references_s(iv)
        integer(kind = ik), intent(in) :: iv
        type(cache_t) :: grid_cache, deriv_cache, resolve_cache, node_cache
        integer(kind = ik) :: st
        character(len = 16) :: label

        write(label, '(A,A)') 'cold vector ', VECTOR_NAMES(iv)

        call cache_init_shared_s(grid_cache, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'grid cache init, ' // trim(label))
        call cache_radius_grid_s(grid_cache, params_ab(:, iv), cold_grid(:, iv), st)
        call assert_int_eq(st, SHAPE_VALID, 'radius_grid, ' // trim(label))
        call cache_free_s(grid_cache)

        call cache_init_shared_s(deriv_cache, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'derivative cache init, ' // trim(label))
        call cache_radius_and_derivative_s(deriv_cache, params_ab(:, iv), &
                cold_radii(:, iv), cold_drs(:, iv), st)
        call assert_int_eq(st, SHAPE_VALID, 'radius_and_derivative, ' // trim(label))
        call cache_free_s(deriv_cache)

        call cache_init_shared_s(resolve_cache, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'resolve cache init, ' // trim(label))
        call cache_resolve_shape_s(resolve_cache, params_ab(:, iv), cold_scalars(1, iv), &
                cold_scalars(2, iv), cold_scalars(3, iv), cold_scalars(4, iv), st)
        call assert_int_eq(st, SHAPE_VALID, 'resolve_shape, ' // trim(label))
        call cache_free_s(resolve_cache)

        call cache_init_shared_s(node_cache, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'node cache init, ' // trim(label))
        call cache_node_radius_and_derivative_s(node_cache, nodes, params_ab(:, iv), &
                cold_node_radii(:, iv), cold_node_drs(:, iv), st)
        call assert_int_eq(st, SHAPE_VALID, 'node_radius_and_derivative, ' // trim(label))
        call cache_free_s(node_cache)
    end subroutine build_cold_references_s

    !> Guard against a vacuous suite: every comparison below would also hold if
    !! A and B produced the same numbers, so the two vectors must be shown to
    !! separate every output first.
    subroutine run_reference_distinctness_s()
        call assert_true(.not. bits_equal_f(cold_grid(:, 1), cold_grid(:, 2)), &
                'cold radius grids of A and B differ')
        call assert_true(.not. bits_equal_f(cold_radii(:, 1), cold_radii(:, 2)), &
                'cold derivative-path radii of A and B differ')
        call assert_true(.not. bits_equal_f(cold_drs(:, 1), cold_drs(:, 2)), &
                'cold derivatives of A and B differ')
        call assert_true(.not. bits_equal_f(cold_scalars(1:1, 1), cold_scalars(1:1, 2)), &
                'cold b10 of A and B differ')
        call assert_true(.not. bits_equal_f(cold_scalars(2:2, 1), cold_scalars(2:2, 2)), &
                'cold r_north of A and B differ')
        call assert_true(.not. bits_equal_f(cold_node_radii(:, 1), cold_node_radii(:, 2)), &
                'cold node radii of A and B differ')
        call assert_true(.not. bits_equal_f(cold_node_drs(:, 1), cold_node_drs(:, 2)), &
                'cold node derivatives of A and B differ')
    end subroutine run_reference_distinctness_s

    !> One long interleaved sequence on a single cache: entry points cycle, the
    !! parameter vector alternates, and every output still equals its cold twin.
    subroutine run_alternation_s()
        type(cache_t) :: cache
        integer(kind = ik) :: s, st
        character(len = 64) :: label

        call cache_init_shared_s(cache, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'alternation cache init')

        do s = 1_ik, N_STEPS
            write(label, '(A,I0,A,A,A,A)') 'step ', s, ': ', &
                    trim(ROUTINE_NAMES(STEP_ROUTINE(s))), ' on ', VECTOR_NAMES(STEP_VECTOR(s))
            call call_and_compare_s(cache, STEP_ROUTINE(s), STEP_VECTOR(s), trim(label))
        end do

        call cache_free_s(cache)
    end subroutine run_alternation_s

    !> Call one entry point on one vector and assert every output against the
    !! cold reference for that (entry point, vector) pair.
    subroutine call_and_compare_s(cache, routine_id, iv, label)
        type(cache_t),      intent(inout) :: cache
        integer(kind = ik), intent(in)    :: routine_id, iv
        character(len = *), intent(in)    :: label
        real(kind = rk) :: radii(N_THETAS), drs(N_THETAS)
        real(kind = rk) :: node_radii(N_NODES), node_drs(N_NODES)
        real(kind = rk) :: scalars(N_SCALARS)
        integer(kind = ik) :: k, st

        select case (routine_id)
        case (1_ik)
            call cache_radius_grid_s(cache, params_ab(:, iv), radii, st)
            call assert_int_eq(st, SHAPE_VALID, 'status, ' // label)
            call assert_true(bits_equal_f(radii, cold_grid(:, iv)), &
                    'radii bitwise interleaved == cold, ' // label)
        case (2_ik)
            call cache_radius_and_derivative_s(cache, params_ab(:, iv), radii, drs, st)
            call assert_int_eq(st, SHAPE_VALID, 'status, ' // label)
            call assert_true(bits_equal_f(radii, cold_radii(:, iv)), &
                    'radii bitwise interleaved == cold, ' // label)
            call assert_true(bits_equal_f(drs, cold_drs(:, iv)), &
                    'derivatives bitwise interleaved == cold, ' // label)
        case (3_ik)
            call cache_resolve_shape_s(cache, params_ab(:, iv), scalars(1), scalars(2), &
                    scalars(3), scalars(4), st)
            call assert_int_eq(st, SHAPE_VALID, 'status, ' // label)
            do k = 1_ik, N_SCALARS
                call assert_true(bits_equal_f(scalars(k:k), cold_scalars(k:k, iv)), &
                        trim(SCALAR_NAMES(k)) // ' bitwise interleaved == cold, ' // label)
            end do
        case default
            call cache_node_radius_and_derivative_s(cache, nodes, params_ab(:, iv), &
                    node_radii, node_drs, st)
            call assert_int_eq(st, SHAPE_VALID, 'status, ' // label)
            call assert_true(bits_equal_f(node_radii, cold_node_radii(:, iv)), &
                    'node radii bitwise interleaved == cold, ' // label)
            call assert_true(bits_equal_f(node_drs, cold_node_drs(:, iv)), &
                    'node derivatives bitwise interleaved == cold, ' // label)
        end select
    end subroutine call_and_compare_s

    !> A shape the validity scan rejects must not poison the next call: the
    !! rejection returns the engine to cold, so the following valid call is a
    !! cold compute and must reproduce the cold reference bit for bit. Checked
    !! both from a fresh cache and from one already warm on the other vector.
    subroutine run_invalid_recovery_s()
        type(cache_t) :: fresh, warm
        real(kind = rk) :: radii(N_THETAS), drs(N_THETAS)
        integer(kind = ik) :: st

        call cache_init_shared_s(fresh, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'invalid-recovery fresh cache init')
        call cache_radius_grid_s(fresh, PARAMS_INVALID, radii, st)
        call assert_int_eq(st, BETA_PARAM_ERROR_NORTH_POLE, &
                'invalid vector rejected on a fresh cache')
        call cache_radius_grid_s(fresh, params_ab(:, 1), radii, st)
        call assert_int_eq(st, SHAPE_VALID, 'status after invalid, fresh cache')
        call assert_true(bits_equal_f(radii, cold_grid(:, 1)), &
                'radii bitwise after invalid == cold A, fresh cache')
        call cache_free_s(fresh)

        call cache_init_shared_s(warm, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'invalid-recovery warm cache init')
        call cache_radius_and_derivative_s(warm, params_ab(:, 2), radii, drs, st)
        call assert_int_eq(st, SHAPE_VALID, 'warm-up on B before the invalid call')
        call cache_radius_grid_s(warm, PARAMS_INVALID, radii, st)
        call assert_int_eq(st, BETA_PARAM_ERROR_NORTH_POLE, &
                'invalid vector rejected on a warm cache')
        call cache_radius_and_derivative_s(warm, params_ab(:, 1), radii, drs, st)
        call assert_int_eq(st, SHAPE_VALID, 'status after invalid, warm cache')
        call assert_true(bits_equal_f(radii, cold_radii(:, 1)), &
                'radii bitwise after invalid == cold A, warm cache')
        call assert_true(bits_equal_f(drs, cold_drs(:, 1)), &
                'derivatives bitwise after invalid == cold A, warm cache')
        call cache_free_s(warm)
    end subroutine run_invalid_recovery_s

    !> A wrong-length parameter vector is rejected with status 4 and, like every
    !! rejection, leaves the cache in a state the next valid call recovers from
    !! exactly.
    subroutine run_wrong_length_recovery_s()
        type(cache_t) :: fresh, warm
        real(kind = rk) :: radii(N_THETAS), drs(N_THETAS)
        integer(kind = ik) :: st

        call cache_init_shared_s(fresh, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'wrong-length fresh cache init')
        call cache_radius_grid_s(fresh, params_ab(1:N_PARAMS - 1_ik, 1), radii, st)
        call assert_int_eq(st, SHAPE_ERROR_WRONG_PARAM_COUNT, &
                'three-element vector rejected on a fresh cache')
        call cache_radius_grid_s(fresh, params_ab(:, 1), radii, st)
        call assert_int_eq(st, SHAPE_VALID, 'status after wrong length, fresh cache')
        call assert_true(bits_equal_f(radii, cold_grid(:, 1)), &
                'radii bitwise after wrong length == cold A, fresh cache')
        call cache_free_s(fresh)

        call cache_init_shared_s(warm, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'wrong-length warm cache init')
        call cache_radius_and_derivative_s(warm, params_ab(:, 2), radii, drs, st)
        call assert_int_eq(st, SHAPE_VALID, 'warm-up on B before the wrong-length call')
        call cache_radius_grid_s(warm, params_ab(1:N_PARAMS - 1_ik, 1), radii, st)
        call assert_int_eq(st, SHAPE_ERROR_WRONG_PARAM_COUNT, &
                'three-element vector rejected on a warm cache')
        call cache_radius_and_derivative_s(warm, params_ab(:, 1), radii, drs, st)
        call assert_int_eq(st, SHAPE_VALID, 'status after wrong length, warm cache')
        call assert_true(bits_equal_f(radii, cold_radii(:, 1)), &
                'radii bitwise after wrong length == cold A, warm cache')
        call assert_true(bits_equal_f(drs, cold_drs(:, 1)), &
                'derivatives bitwise after wrong length == cold A, warm cache')
        call cache_free_s(warm)
    end subroutine run_wrong_length_recovery_s

    !> Stamp-hole regression guard. The unchecked path stops at I_RESOLVED and
    !! leaves I_MIN_RADIUS, I_VOLUME and I_RADII cold; a checked call on the same
    !! vector must then fill exactly those holes and produce the cold answer. A
    !! stamp wrongly marked valid by the unchecked call would show up here as a
    !! wrong (unscaled) radius.
    !!
    !! The counter check on the same cache pins the unchecked path's own
    !! minimality: cold, it computes I_RESOLVED and nothing else.
    subroutine run_unchecked_then_checked_s()
        type(cache_t) :: cache
        real(kind = rk) :: radii(N_THETAS)
        integer(kind = ik) :: want(N_TRACKED)
        integer(kind = ik) :: k, st

        call cache_init_shared_s(cache, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'unchecked-then-checked cache init')

        call cache_radius_grid_unchecked_s(cache, params_ab(:, 1), radii, st)
        call assert_int_eq(st, SHAPE_VALID, 'cold unchecked call on A')
        want = [1_ik, 0_ik, 0_ik, 0_ik, 0_ik]
        do k = 1_ik, N_TRACKED
            call assert_int_eq(int(cache_recompute_count_f(cache, INTERMEDIATES(k)), ik), &
                    want(k), trim(INTERMEDIATE_NAMES(k)) // ' count after the cold unchecked call')
        end do

        call cache_radius_grid_s(cache, params_ab(:, 1), radii, st)
        call assert_int_eq(st, SHAPE_VALID, 'checked call after the unchecked call')
        call assert_true(bits_equal_f(radii, cold_grid(:, 1)), &
                'radii bitwise checked-after-unchecked == cold A')

        call cache_free_s(cache)
    end subroutine run_unchecked_then_checked_s

    !> Ordering guard: the node path shares intermediates 1-3 with the grid path
    !! but writes nothing into the cached grids. Evaluating nodes first must
    !! therefore leave a following primary-grid call exactly cold.
    subroutine run_node_before_primary_s()
        type(cache_t) :: cache
        real(kind = rk) :: node_radii(N_NODES), node_drs(N_NODES)
        real(kind = rk) :: radii(N_THETAS), drs(N_THETAS)
        integer(kind = ik) :: st

        call cache_init_shared_s(cache, tables, N_PARAMS, .true., .true., st)
        call assert_int_eq(st, SHAPE_VALID, 'node-before-primary cache init')

        call cache_node_radius_and_derivative_s(cache, nodes, params_ab(:, 1), &
                node_radii, node_drs, st)
        call assert_int_eq(st, SHAPE_VALID, 'node call before any primary call')
        call assert_true(bits_equal_f(node_radii, cold_node_radii(:, 1)), &
                'node radii bitwise node-first == cold A')
        call assert_true(bits_equal_f(node_drs, cold_node_drs(:, 1)), &
                'node derivatives bitwise node-first == cold A')

        call cache_radius_and_derivative_s(cache, params_ab(:, 1), radii, drs, st)
        call assert_int_eq(st, SHAPE_VALID, 'primary call after the node call')
        call assert_true(bits_equal_f(radii, cold_radii(:, 1)), &
                'radii bitwise primary-after-node == cold A')
        call assert_true(bits_equal_f(drs, cold_drs(:, 1)), &
                'derivatives bitwise primary-after-node == cold A')

        call cache_free_s(cache)
    end subroutine run_node_before_primary_s

    !> Exact equality on the bit patterns: -Wcompare-reals rejects `==` on
    !! reals, and a scalar transfer avoids the array temporary that
    !! -Warray-temporaries flags.
    logical function bits_equal_f(a, b) result(ok)
        real(kind = rk), intent(in) :: a(:), b(:)
        integer(kind = ik) :: j
        ok = size(a, kind = ik) == size(b, kind = ik)
        if (.not. ok) return
        do j = 1_ik, size(a, kind = ik)
            if (transfer(a(j), 0_ikl) /= transfer(b(j), 0_ikl)) then
                ok = .false.
                return
            end if
        end do
    end function bits_equal_f

end program beta_param_interleave_test
