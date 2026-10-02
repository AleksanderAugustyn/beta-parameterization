!> Contract family 2: statelessness. A cache carries nothing from one call to
!! the next: a repeated call, a call after a rejected call, and a call in the
!! middle of any mixed sequence all return what a fresh cache would have
!! returned from that single call alone. All comparisons here are on IEEE754
!! bit patterns; no tolerances appear anywhere in this suite.
program beta_param_statelessness_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary, bits_equal_f
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_resolve_shape_s, cache_radius_grid_s, &
            cache_radius_and_derivative_s, cache_radius_grid_unchecked_s, &
            cache_node_radius_and_derivative_s, &
            node_set_t, node_set_build_s, node_set_free_s, &
            SHAPE_VALID, SHAPE_ERROR_WRONG_PARAM_COUNT, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_SOUTH_POLE, &
            BETA_PARAM_ERROR_INTERIOR_NEGATIVE, BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, BETA_PARAM_ERROR_NODE_SET_MISMATCH

    implicit none

    integer(kind = ik), parameter :: N_THETAS   = 32_ik
    integer(kind = ik), parameter :: N_NODES    = 8_ik
    integer(kind = ik), parameter :: MAX_PARAMS = 8_ik
    integer(kind = ik), parameter :: N_ROUTINES = 5_ik
    integer(kind = ik), parameter :: N_VECTORS  = 3_ik
    integer(kind = ik), parameter :: N_REGIMES  = 4_ik
    !> Everything one call can return, flattened: radii, derivatives, node
    !! radii, node derivatives, four scalars. Unused slots stay zero.
    integer(kind = ik), parameter :: N_OUT = 2_ik * N_THETAS + 2_ik * N_NODES + 4_ik

    !> Three valid vectors of three different lengths.
    real(kind = rk), parameter :: VECTOR_A(4) = [0.05_rk, 0.20_rk, 0.10_rk, 0.02_rk]
    real(kind = rk), parameter :: VECTOR_B(8) = &
            [0.02_rk, 0.15_rk, 0.08_rk, -0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, -0.002_rk]
    real(kind = rk), parameter :: VECTOR_C(2) = [-0.03_rk, 0.25_rk]

    character(len = 26), parameter :: ROUTINE_NAMES(N_ROUTINES) = &
            ['radius_grid               ', 'radius_and_derivative     ', &
             'resolve_shape             ', 'node_radius_and_derivative', &
             'radius_grid_unchecked     ']

    type(cache_t) :: shared_cache
    !> Two valid node sets over different angles. Node evaluations alternate
    !! between them (see node_choice_f), so a table left over from the previous
    !! node call would show up as a bit mismatch.
    type(node_set_t) :: nodes_a, nodes_b
    real(kind = rk) :: thetas(N_THETAS), node_thetas_a(N_NODES), node_thetas_b(N_NODES)
    !> One fresh-cache, single-call reference per (routine, vector, regime).
    real(kind = rk) :: reference(N_OUT, N_ROUTINES, N_VECTORS, N_REGIMES)
    integer(kind = ik) :: i, status

    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do
    do i = 1_ik, N_NODES
        node_thetas_a(i) = real(2_ik * i - 1_ik, rk) * PI_C / real(2_ik * N_NODES, rk)
        node_thetas_b(i) = real(i, rk) * PI_C / real(N_NODES + 1_ik, rk)
    end do

    call cache_init_s(shared_cache, MAX_PARAMS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'shared cache init')
    call node_set_build_s(nodes_a, shared_cache, node_thetas_a, status)
    call assert_int_eq(status, SHAPE_VALID, 'node set A build')
    call node_set_build_s(nodes_b, shared_cache, node_thetas_b, status)
    call assert_int_eq(status, SHAPE_VALID, 'node set B build')

    call build_references_s()
    call run_references_distinct_s()
    call run_repeat_s()
    call run_after_rejection_s()
    call run_mixed_sequence_s()

    call node_set_free_s(nodes_a)
    call node_set_free_s(nodes_b)
    call cache_free_s(shared_cache)
    call test_summary()

contains

    !> Which node set a node evaluation of vector `v` in regime `g` uses:
    !! 1 = A, 2 = B. A fixed function of (v, g), so the reference and the
    !! compared call agree, while consecutive node calls in the mixed sequence
    !! switch sets.
    pure function node_choice_f(v, g) result(choice)
        integer(kind = ik), intent(in) :: v, g
        integer(kind = ik) :: choice
        choice = mod(v + g, 2_ik) + 1_ik
    end function node_choice_f

    !> One call of routine `r` on vector `v` in regime `g` against `cache`,
    !! with every output flattened into `out`.
    subroutine call_one_s(cache, r, v, g, out, status)
        type(cache_t),      intent(in)  :: cache
        integer(kind = ik), intent(in)  :: r, v, g
        real(kind = rk),    intent(out) :: out(N_OUT)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik), parameter :: I_RADII = 1_ik, I_DRS = I_RADII + N_THETAS, &
                I_NODE_RADII = I_DRS + N_THETAS, I_NODE_DRS = I_NODE_RADII + N_NODES, &
                I_SCALARS = I_NODE_DRS + N_NODES
        real(kind = rk) :: params(MAX_PARAMS)
        integer(kind = ik) :: n
        logical :: conserve_volume, apply_com

        conserve_volume = (g == 2_ik .or. g == 4_ik)
        apply_com       = (g == 3_ik .or. g == 4_ik)

        params(:) = 0.0_rk
        select case (v)
        case (1_ik)
            n = size(VECTOR_A, kind = ik)
            params(1:n) = VECTOR_A
        case (2_ik)
            n = size(VECTOR_B, kind = ik)
            params(1:n) = VECTOR_B
        case default
            n = size(VECTOR_C, kind = ik)
            params(1:n) = VECTOR_C
        end select

        out(:) = 0.0_rk
        select case (r)
        case (1_ik)
            call cache_radius_grid_s(cache, params(1:n), conserve_volume, apply_com, &
                    out(I_RADII:I_DRS - 1_ik), status)
        case (2_ik)
            call cache_radius_and_derivative_s(cache, params(1:n), conserve_volume, &
                    apply_com, out(I_RADII:I_DRS - 1_ik), &
                    out(I_DRS:I_NODE_RADII - 1_ik), status)
        case (3_ik)
            call cache_resolve_shape_s(cache, params(1:n), conserve_volume, apply_com, &
                    out(I_SCALARS), out(I_SCALARS + 1_ik), out(I_SCALARS + 2_ik), &
                    out(I_SCALARS + 3_ik), status)
        case (4_ik)
            if (node_choice_f(v, g) == 1_ik) then
                call cache_node_radius_and_derivative_s(cache, params(1:n), nodes_a, &
                        conserve_volume, apply_com, out(I_NODE_RADII:I_NODE_DRS - 1_ik), &
                        out(I_NODE_DRS:I_SCALARS - 1_ik), status)
            else
                call cache_node_radius_and_derivative_s(cache, params(1:n), nodes_b, &
                        conserve_volume, apply_com, out(I_NODE_RADII:I_NODE_DRS - 1_ik), &
                        out(I_NODE_DRS:I_SCALARS - 1_ik), status)
            end if
        case default
            call cache_radius_grid_unchecked_s(cache, params(1:n), apply_com, &
                    out(I_RADII:I_DRS - 1_ik), status)
        end select
    end subroutine call_one_s

    !> The reference half of every comparison: a fresh cache per entry, exactly
    !! one call on it.
    subroutine build_references_s()
        type(cache_t) :: fresh
        integer(kind = ik) :: r, v, g, st
        logical :: all_valid

        all_valid = .true.
        do g = 1_ik, N_REGIMES
            do v = 1_ik, N_VECTORS
                do r = 1_ik, N_ROUTINES
                    call cache_init_s(fresh, MAX_PARAMS, thetas, st)
                    all_valid = all_valid .and. (st == SHAPE_VALID)
                    call call_one_s(fresh, r, v, g, reference(:, r, v, g), st)
                    all_valid = all_valid .and. (st == SHAPE_VALID)
                    call cache_free_s(fresh)
                end do
            end do
        end do
        call assert_true(all_valid, 'every reference call returned SHAPE_VALID')
    end subroutine build_references_s

    !> Guard against a vacuous suite: the comparisons below would also hold if
    !! vectors or option combinations produced the same numbers.
    subroutine run_references_distinct_s()
        real(kind = rk) :: radii_a(N_NODES), drs_a(N_NODES), radii_b(N_NODES), drs_b(N_NODES)
        integer(kind = ik) :: st

        call assert_true(.not. bits_equal_f(reference(:, 1, 1, 4), reference(:, 1, 2, 4)), &
                'vectors A and B give different radius grids')
        call assert_true(.not. bits_equal_f(reference(:, 1, 1, 4), reference(:, 1, 3, 4)), &
                'vectors A and C give different radius grids')
        call assert_true(.not. bits_equal_f(reference(:, 1, 1, 1), reference(:, 1, 1, 2)), &
                'conserve_volume changes the radius grid')
        call assert_true(.not. bits_equal_f(reference(:, 1, 1, 1), reference(:, 1, 1, 3)), &
                'apply_com changes the radius grid')
        call assert_true(.not. bits_equal_f(reference(:, 5, 1, 1), reference(:, 5, 1, 3)), &
                'apply_com changes the unchecked radius grid')
        ! The two node sets must separate the node outputs, and the mixed
        ! sequence must reach both of them.
        call cache_node_radius_and_derivative_s(shared_cache, VECTOR_A, nodes_a, &
                .true., .true., radii_a, drs_a, st)
        call assert_int_eq(st, SHAPE_VALID, 'node set A evaluates')
        call cache_node_radius_and_derivative_s(shared_cache, VECTOR_A, nodes_b, &
                .true., .true., radii_b, drs_b, st)
        call assert_int_eq(st, SHAPE_VALID, 'node set B evaluates')
        call assert_true(.not. bits_equal_f(radii_a, radii_b), &
                'node sets A and B give different radii')
        call assert_true(.not. bits_equal_f(drs_a, drs_b), &
                'node sets A and B give different derivatives')
        call assert_int_eq(node_choice_f(1_ik, 1_ik), 1_ik, 'vector 1, regime 1 uses node set A')
        call assert_int_eq(node_choice_f(1_ik, 2_ik), 2_ik, 'vector 1, regime 2 uses node set B')
    end subroutine run_references_distinct_s

    !> The same call twice in a row on the shared cache.
    subroutine run_repeat_s()
        real(kind = rk) :: out(N_OUT)
        integer(kind = ik) :: r, pass, st

        do r = 1_ik, N_ROUTINES
            do pass = 1_ik, 2_ik
                call call_one_s(shared_cache, r, 1_ik, 4_ik, out, st)
                call assert_int_eq(st, SHAPE_VALID, &
                        'repeat status, ' // trim(ROUTINE_NAMES(r)))
                call assert_true(bits_equal_f(out, reference(:, r, 1_ik, 4_ik)), &
                        'repeated call == reference, ' // trim(ROUTINE_NAMES(r)))
            end do
        end do
    end subroutine run_repeat_s

    !> Every kind of rejected call on the shared cache, each followed by every
    !! routine on a valid vector: the rejection leaves no trace.
    subroutine run_after_rejection_s()
        type(cache_t) :: small_cache
        type(node_set_t) :: small_nodes
        real(kind = rk) :: radii(N_THETAS), bad(5), node_radii(N_NODES), node_drs(N_NODES)
        real(kind = rk) :: nine(9)
        integer(kind = ik) :: st

        nine(:) = 0.01_rk
        call cache_radius_grid_s(shared_cache, nine, .true., .true., radii, st)
        call assert_int_eq(st, SHAPE_ERROR_WRONG_PARAM_COUNT, 'nine parameters rejected')
        call check_all_routines_s('after wrong parameter count')

        call cache_radius_grid_s(shared_cache, VECTOR_A, .true., .true., bad, st)
        call assert_int_eq(st, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, 'bad buffer rejected')
        call check_all_routines_s('after bad buffer')

        call cache_init_s(small_cache, 2_ik, thetas, st)
        call assert_int_eq(st, SHAPE_VALID, 'small cache init')
        call node_set_build_s(small_nodes, small_cache, node_thetas_a, st)
        call assert_int_eq(st, SHAPE_VALID, 'small node set build')
        call cache_node_radius_and_derivative_s(shared_cache, VECTOR_A, small_nodes, &
                .true., .true., node_radii, node_drs, st)
        call assert_int_eq(st, BETA_PARAM_ERROR_NODE_SET_MISMATCH, 'small node set rejected')
        call check_all_routines_s('after node-set mismatch')
        call node_set_free_s(small_nodes)
        call cache_free_s(small_cache)

        call cache_radius_grid_s(shared_cache, [0.0_rk, -1.8_rk], .true., .true., radii, st)
        call assert_int_eq(st, BETA_PARAM_ERROR_NORTH_POLE, 'north-pole shape rejected')
        call check_all_routines_s('after north-pole rejection')

        call cache_radius_grid_s(shared_cache, [0.0_rk, 0.0_rk, 2.0_rk], .true., .true., &
                radii, st)
        call assert_int_eq(st, BETA_PARAM_ERROR_SOUTH_POLE, 'south-pole shape rejected')
        call check_all_routines_s('after south-pole rejection')

        call cache_radius_grid_s(shared_cache, [0.0_rk, 4.0_rk], .true., .true., radii, st)
        call assert_int_eq(st, BETA_PARAM_ERROR_INTERIOR_NEGATIVE, 'interior shape rejected')
        call check_all_routines_s('after interior rejection')

        call cache_radius_grid_s(shared_cache, [1.0e10_rk], .true., .true., radii, st)
        call assert_int_eq(st, BETA_PARAM_ERROR_COM_NOT_CONVERGED, 'COM failure rejected')
        call check_all_routines_s('after COM rejection')
    end subroutine run_after_rejection_s

    !> Every routine on vector A, regime TT, against its reference.
    subroutine check_all_routines_s(label)
        character(len = *), intent(in) :: label
        real(kind = rk) :: out(N_OUT)
        integer(kind = ik) :: r, st

        do r = 1_ik, N_ROUTINES
            call call_one_s(shared_cache, r, 1_ik, 4_ik, out, st)
            call assert_int_eq(st, SHAPE_VALID, &
                    'status ' // label // ', ' // trim(ROUTINE_NAMES(r)))
            call assert_true(bits_equal_f(out, reference(:, r, 1_ik, 4_ik)), &
                    'bits ' // label // ', ' // trim(ROUTINE_NAMES(r)))
        end do
    end subroutine check_all_routines_s

    !> One long sequence on the shared cache in which routine, vector (and so
    !! vector length) and option combination all change on every step. 5, 3 and
    !! 4 are pairwise coprime, so 60 steps visit each combination exactly once.
    subroutine run_mixed_sequence_s()
        real(kind = rk) :: out(N_OUT)
        integer(kind = ik) :: s, r, v, g, st
        character(len = 64) :: label

        do s = 0_ik, N_ROUTINES * N_VECTORS * N_REGIMES - 1_ik
            r = mod(s, N_ROUTINES) + 1_ik
            v = mod(s, N_VECTORS) + 1_ik
            g = mod(s, N_REGIMES) + 1_ik
            write(label, '(A,I0,A,A,A,I0,A,I0)') 'step ', s, ': ', &
                    trim(ROUTINE_NAMES(r)), ', vector ', v, ', regime ', g
            call call_one_s(shared_cache, r, v, g, out, st)
            call assert_int_eq(st, SHAPE_VALID, 'status, ' // trim(label))
            call assert_true(bits_equal_f(out, reference(:, r, v, g)), &
                    'mixed sequence == reference, ' // trim(label))
        end do
    end subroutine run_mixed_sequence_s

end program beta_param_statelessness_test
