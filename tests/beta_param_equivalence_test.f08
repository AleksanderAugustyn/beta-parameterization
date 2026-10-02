!> Contract family 1: equivalence. A one-shot call and a cached call agree —
!! same status, bitwise-identical outputs — for any cache whose max_params is
!! at least the vector length, and a short vector agrees with its zero-padded
!! form. Every comparison here is on the IEEE754 bit patterns; no tolerances
!! appear anywhere in this suite.
program beta_param_equivalence_test
    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary, &
            bits_equal_f, all_zero_f
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_resolve_shape_s, cache_radius_grid_s, &
            cache_radius_and_derivative_s, cache_radius_grid_unchecked_s, &
            cache_node_radius_and_derivative_s, &
            node_set_t, node_set_build_s, node_set_free_s, &
            compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, &
            SHAPE_MAX_PARAMS, SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_SOUTH_POLE, &
            BETA_PARAM_ERROR_INTERIOR_NEGATIVE, BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, BETA_PARAM_ERROR_POLE_NODE

    implicit none

    integer(kind = ik), parameter :: N_THETAS = 32_ik
    integer(kind = ik), parameter :: N_NODES  = 8_ik
    integer(kind = ik), parameter :: N_SHORT  = 3_ik
    integer(kind = ik), parameter :: N_PADDED = 8_ik

    !> IEEE754 -0.0; built from bits so that -fno-signed-zeros cannot fold it
    !! into +0.0 the way a `-0.0_rk` literal could.
    integer(kind = ikl), parameter :: NEG_ZERO_BITS = int(z'8000000000000000', kind = ikl)

    real(kind = rk), parameter :: PARAMS_1(1) = [0.15_rk]
    real(kind = rk), parameter :: PARAMS_4(4) = [0.05_rk, 0.20_rk, 0.10_rk, 0.02_rk]
    real(kind = rk), parameter :: PARAMS_8(8) = &
            [0.02_rk, 0.15_rk, 0.08_rk, -0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, -0.002_rk]
    real(kind = rk), parameter :: PARAMS_SHORT(N_SHORT) = [0.05_rk, 0.20_rk, 0.10_rk]

    !> One vector per rejection code. NORTH/SOUTH/INTERIOR fail the validity
    !! gate in every option combination; COM_103 and NEG_VOLUME fail the COM
    !! correction when it is on and a pole when it is off.
    real(kind = rk), parameter :: PARAMS_NORTH(4)      = [0.0_rk, -1.8_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: PARAMS_SOUTH(4)      = [0.0_rk, 0.0_rk, 2.0_rk, 0.0_rk]
    real(kind = rk), parameter :: PARAMS_INTERIOR(4)   = [0.0_rk, 4.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: PARAMS_COM_103(1)    = [1.0e10_rk]
    real(kind = rk), parameter :: PARAMS_NEG_VOLUME(2) = [0.0_rk, -20.0_rk]

    real(kind = rk) :: thetas(N_THETAS), node_thetas(N_NODES)
    real(kind = rk) :: params_64(SHAPE_MAX_PARAMS)
    integer(kind = ik) :: i

    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do
    do i = 1_ik, N_NODES
        node_thetas(i) = real(2_ik * i - 1_ik, rk) * PI_C / real(2_ik * N_NODES, rk)
    end do
    ! Alternating, decaying amplitudes: a valid 64-parameter shape.
    do i = 1_ik, SHAPE_MAX_PARAMS
        params_64(i) = 0.02_rk / real(i, rk)
        if (mod(i, 2_ik) == 0_ik) params_64(i) = -params_64(i)
    end do

    ! E1: one-shot == cached, every option combination, several cache sizes.
    call run_tiers_agree_s(PARAMS_1, 'n = 1')
    call run_tiers_agree_s(PARAMS_4, 'n = 4')
    call run_tiers_agree_s(PARAMS_8, 'n = 8')
    call run_tiers_agree_s(params_64, 'n = 64')

    ! E2 / E3: every rejection code, both tiers.
    call run_shape_rejections_s()
    call run_buffer_rejection_s()
    call run_theta_rejections_s()

    ! E4: short == zero-padded, +0.0 and -0.0 padding, every output.
    call run_padding_s(.false.)
    call run_padding_s(.true.)

    call test_summary()

contains

    !> (conserve_volume, apply_com) for regime 1..4: FF, TF, FT, TT.
    pure subroutine regime_flags_s(regime, conserve_volume, apply_com)
        integer(kind = ik), intent(in)  :: regime
        logical,            intent(out) :: conserve_volume, apply_com
        conserve_volume = (regime == 2_ik .or. regime == 4_ik)
        apply_com       = (regime == 3_ik .or. regime == 4_ik)
    end subroutine regime_flags_s

    !> E1. For one valid vector: both radius outputs agree between the one-shot
    !! tier and caches built with max_params = n, n + 3 and 64.
    subroutine run_tiers_agree_s(params, label)
        real(kind = rk),    intent(in) :: params(:)
        character(len = *), intent(in) :: label

        type(cache_t) :: cache
        real(kind = rk) :: grid_sa(N_THETAS), radii_sa(N_THETAS), drs_sa(N_THETAS)
        real(kind = rk) :: grid_c(N_THETAS), radii_c(N_THETAS), drs_c(N_THETAS)
        integer(kind = ik) :: sizes(3)
        integer(kind = ik) :: regime, k, n, status
        logical :: conserve_volume, apply_com
        character(len = 64) :: tag

        n = size(params, kind = ik)
        sizes(1) = n
        sizes(2) = min(n + 3_ik, SHAPE_MAX_PARAMS)
        sizes(3) = SHAPE_MAX_PARAMS

        do regime = 1_ik, 4_ik
            call regime_flags_s(regime, conserve_volume, apply_com)

            call compute_radius_grid_standalone_s(params, thetas, conserve_volume, &
                    apply_com, grid_sa, status)
            call assert_int_eq(status, SHAPE_VALID, 'one-shot grid valid, ' // label)
            call compute_radius_and_derivative_standalone_s(params, thetas, &
                    conserve_volume, apply_com, radii_sa, drs_sa, status)
            call assert_int_eq(status, SHAPE_VALID, 'one-shot derivative valid, ' // label)

            do k = 1_ik, 3_ik
                write(tag, '(A,A,I0,A,I0)') label, ', regime ', regime, &
                        ', max_params ', sizes(k)
                call cache_init_s(cache, sizes(k), thetas, status)
                call assert_int_eq(status, SHAPE_VALID, 'cache init, ' // trim(tag))

                call cache_radius_grid_s(cache, params, conserve_volume, apply_com, &
                        grid_c, status)
                call assert_int_eq(status, SHAPE_VALID, 'cached grid valid, ' // trim(tag))
                call assert_true(bits_equal_f(grid_c, grid_sa), &
                        'grid: one-shot == cached, ' // trim(tag))

                call cache_radius_and_derivative_s(cache, params, conserve_volume, &
                        apply_com, radii_c, drs_c, status)
                call assert_int_eq(status, SHAPE_VALID, 'cached derivative valid, ' // trim(tag))
                call assert_true(bits_equal_f(radii_c, radii_sa), &
                        'radii: one-shot == cached, ' // trim(tag))
                call assert_true(bits_equal_f(drs_c, drs_sa), &
                        'derivatives: one-shot == cached, ' // trim(tag))

                call cache_free_s(cache)
            end do
        end do
    end subroutine run_tiers_agree_s

    !> E2 / E3. Each shape rejection code comes back from both tiers with
    !! zero-filled outputs, in every option combination.
    subroutine run_shape_rejections_s()
        integer(kind = ik) :: regime
        logical :: conserve_volume, apply_com

        do regime = 1_ik, 4_ik
            call regime_flags_s(regime, conserve_volume, apply_com)
            call check_rejection_s(PARAMS_NORTH, conserve_volume, apply_com, &
                    BETA_PARAM_ERROR_NORTH_POLE, 'north pole')
            call check_rejection_s(PARAMS_SOUTH, conserve_volume, apply_com, &
                    BETA_PARAM_ERROR_SOUTH_POLE, 'south pole')
            call check_rejection_s(PARAMS_INTERIOR, conserve_volume, apply_com, &
                    BETA_PARAM_ERROR_INTERIOR_NEGATIVE, 'interior')
            if (apply_com) then
                call check_rejection_s(PARAMS_COM_103, conserve_volume, apply_com, &
                        BETA_PARAM_ERROR_COM_NOT_CONVERGED, 'COM not converged')
                call check_rejection_s(PARAMS_NEG_VOLUME, conserve_volume, apply_com, &
                        BETA_PARAM_ERROR_COM_NOT_CONVERGED, 'negative volume integral')
            else
                call check_rejection_s(PARAMS_COM_103, conserve_volume, apply_com, &
                        BETA_PARAM_ERROR_SOUTH_POLE, 'beta10 = 1e10 without COM')
                call check_rejection_s(PARAMS_NEG_VOLUME, conserve_volume, apply_com, &
                        BETA_PARAM_ERROR_NORTH_POLE, 'beta20 = -20 without COM')
            end if
        end do
    end subroutine run_shape_rejections_s

    !> One rejected vector, both tiers, both radius outputs: `expected` status
    !! and zero-filled buffers everywhere.
    subroutine check_rejection_s(params, conserve_volume, apply_com, expected, label)
        real(kind = rk),    intent(in) :: params(:)
        logical,            intent(in) :: conserve_volume, apply_com
        integer(kind = ik), intent(in) :: expected
        character(len = *), intent(in) :: label

        type(cache_t) :: cache
        real(kind = rk) :: radii(N_THETAS), drs(N_THETAS)
        integer(kind = ik) :: status

        call cache_init_s(cache, 8_ik, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'cache init, ' // label)

        radii(:) = 1.0_rk
        call cache_radius_grid_s(cache, params, conserve_volume, apply_com, radii, status)
        call assert_int_eq(status, expected, 'cached grid status, ' // label)
        call assert_true(all_zero_f(radii), 'cached grid zero-filled, ' // label)

        radii(:) = 1.0_rk
        drs(:)   = 1.0_rk
        call cache_radius_and_derivative_s(cache, params, conserve_volume, apply_com, &
                radii, drs, status)
        call assert_int_eq(status, expected, 'cached derivative status, ' // label)
        call assert_true(all_zero_f(radii) .and. all_zero_f(drs), &
                'cached derivative zero-filled, ' // label)
        call cache_free_s(cache)

        radii(:) = 1.0_rk
        call compute_radius_grid_standalone_s(params, thetas, conserve_volume, &
                apply_com, radii, status)
        call assert_int_eq(status, expected, 'one-shot grid status, ' // label)
        call assert_true(all_zero_f(radii), 'one-shot grid zero-filled, ' // label)

        radii(:) = 1.0_rk
        drs(:)   = 1.0_rk
        call compute_radius_and_derivative_standalone_s(params, thetas, &
                conserve_volume, apply_com, radii, drs, status)
        call assert_int_eq(status, expected, 'one-shot derivative status, ' // label)
        call assert_true(all_zero_f(radii) .and. all_zero_f(drs), &
                'one-shot derivative zero-filled, ' // label)
    end subroutine check_rejection_s

    !> E2. A wrong output buffer size is 104 in both tiers.
    subroutine run_buffer_rejection_s()
        type(cache_t) :: cache
        real(kind = rk) :: short_radii(N_THETAS - 1_ik)
        integer(kind = ik) :: status

        call cache_init_s(cache, 4_ik, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'buffer-rejection cache init')

        short_radii(:) = 1.0_rk
        call cache_radius_grid_s(cache, PARAMS_4, .true., .true., short_radii, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, 'cached bad buffer')
        call assert_true(all_zero_f(short_radii), 'cached bad buffer zero-filled')
        call cache_free_s(cache)

        short_radii(:) = 1.0_rk
        call compute_radius_grid_standalone_s(PARAMS_4, thetas, .true., .true., &
                short_radii, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, 'one-shot bad buffer')
        call assert_true(all_zero_f(short_radii), 'one-shot bad buffer zero-filled')
    end subroutine run_buffer_rejection_s

    !> E2. A bad theta set gets the same code from the cached tier (at init)
    !! and from the one-shot tier (at the call).
    subroutine run_theta_rejections_s()
        type(cache_t) :: cache
        real(kind = rk) :: one_theta(1), one_radius(1)
        real(kind = rk) :: pole_thetas(2), two_radii(2)
        integer(kind = ik) :: status

        one_theta = [0.7_rk]
        call cache_init_s(cache, 4_ik, one_theta, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'cached: one theta')
        call compute_radius_grid_standalone_s(PARAMS_4, one_theta, .false., .false., &
                one_radius, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'one-shot: one theta')

        pole_thetas = [0.0_rk, 1.0_rk]
        call cache_init_s(cache, 4_ik, pole_thetas, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, 'cached: pole theta')
        call compute_radius_grid_standalone_s(PARAMS_4, pole_thetas, .false., .false., &
                two_radii, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, 'one-shot: pole theta')
    end subroutine run_theta_rejections_s

    !> E4. PARAMS_SHORT against its forms padded to 5 and to 8 entries, on one
    !! max_params = 8 cache: all five cached outputs and both one-shot outputs,
    !! every option combination, every bit.
    !!
    !! @param[in] negative_zero  Pad with -0.0 instead of +0.0
    subroutine run_padding_s(negative_zero)
        logical, intent(in) :: negative_zero

        type(cache_t) :: cache
        type(node_set_t) :: nodes
        !> volatile: under -fno-signed-zeros (Release -ffast-math) the compiler
        !! may fold a -0.0 store into +0.0; volatile forces the bits to memory.
        real(kind = rk), volatile :: padded(N_PADDED)
        real(kind = rk) :: grid(N_THETAS, 2), radii(N_THETAS, 2), drs(N_THETAS, 2)
        real(kind = rk) :: unchecked(N_THETAS, 2), scalars(4, 2)
        real(kind = rk) :: node_radii(N_NODES, 2), node_drs(N_NODES, 2)
        real(kind = rk) :: grid_sa(N_THETAS, 2), radii_sa(N_THETAS, 2), drs_sa(N_THETAS, 2)
        integer(kind = ik) :: lengths(2)
        integer(kind = ik) :: regime, k, n, status
        logical :: conserve_volume, apply_com, all_valid
        character(len = 64) :: tag

        call cache_init_s(cache, N_PADDED, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'padding cache init')
        call node_set_build_s(nodes, cache, node_thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'padding node set build')

        padded(1:N_SHORT) = PARAMS_SHORT
        if (negative_zero) then
            padded(N_SHORT + 1_ik:) = transfer(NEG_ZERO_BITS, 1.0_rk)
            call assert_true(transfer(padded(N_PADDED), 0_ikl) == NEG_ZERO_BITS, &
                    'the -0.0 padding reached memory')
        else
            padded(N_SHORT + 1_ik:) = 0.0_rk
        end if
        lengths(1) = 5_ik
        lengths(2) = N_PADDED

        do regime = 1_ik, 4_ik
            call regime_flags_s(regime, conserve_volume, apply_com)
            all_valid = .true.

            ! Column 1: the short vector.
            call cache_radius_grid_s(cache, PARAMS_SHORT, conserve_volume, apply_com, &
                    grid(:, 1), status)
            all_valid = all_valid .and. (status == SHAPE_VALID)
            call cache_radius_and_derivative_s(cache, PARAMS_SHORT, conserve_volume, &
                    apply_com, radii(:, 1), drs(:, 1), status)
            all_valid = all_valid .and. (status == SHAPE_VALID)
            call cache_resolve_shape_s(cache, PARAMS_SHORT, conserve_volume, apply_com, &
                    scalars(1, 1), scalars(2, 1), scalars(3, 1), scalars(4, 1), status)
            all_valid = all_valid .and. (status == SHAPE_VALID)
            call cache_node_radius_and_derivative_s(cache, PARAMS_SHORT, nodes, &
                    conserve_volume, apply_com, node_radii(:, 1), node_drs(:, 1), status)
            all_valid = all_valid .and. (status == SHAPE_VALID)
            call cache_radius_grid_unchecked_s(cache, PARAMS_SHORT, apply_com, &
                    unchecked(:, 1), status)
            all_valid = all_valid .and. (status == SHAPE_VALID)
            call compute_radius_grid_standalone_s(PARAMS_SHORT, thetas, conserve_volume, &
                    apply_com, grid_sa(:, 1), status)
            all_valid = all_valid .and. (status == SHAPE_VALID)
            call compute_radius_and_derivative_standalone_s(PARAMS_SHORT, thetas, &
                    conserve_volume, apply_com, radii_sa(:, 1), drs_sa(:, 1), status)
            all_valid = all_valid .and. (status == SHAPE_VALID)

            ! Column 2: each padded form in turn.
            do k = 1_ik, 2_ik
                n = lengths(k)
                write(tag, '(A,I0,A,I0,A,L1)') 'regime ', regime, ', padded to ', n, &
                        ', -0.0: ', negative_zero

                call cache_radius_grid_s(cache, padded(1:n), conserve_volume, apply_com, &
                        grid(:, 2), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call cache_radius_and_derivative_s(cache, padded(1:n), conserve_volume, &
                        apply_com, radii(:, 2), drs(:, 2), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call cache_resolve_shape_s(cache, padded(1:n), conserve_volume, apply_com, &
                        scalars(1, 2), scalars(2, 2), scalars(3, 2), scalars(4, 2), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call cache_node_radius_and_derivative_s(cache, padded(1:n), nodes, &
                        conserve_volume, apply_com, node_radii(:, 2), node_drs(:, 2), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call cache_radius_grid_unchecked_s(cache, padded(1:n), apply_com, &
                        unchecked(:, 2), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call compute_radius_grid_standalone_s(padded(1:n), thetas, &
                        conserve_volume, apply_com, grid_sa(:, 2), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call compute_radius_and_derivative_standalone_s(padded(1:n), thetas, &
                        conserve_volume, apply_com, radii_sa(:, 2), drs_sa(:, 2), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)

                call assert_true(bits_equal_f(grid(:, 1), grid(:, 2)), &
                        'cached grid: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(radii(:, 1), radii(:, 2)), &
                        'cached radii: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(drs(:, 1), drs(:, 2)), &
                        'cached derivatives: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(scalars(:, 1), scalars(:, 2)), &
                        'resolve scalars: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(node_radii(:, 1), node_radii(:, 2)), &
                        'node radii: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(node_drs(:, 1), node_drs(:, 2)), &
                        'node derivatives: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(unchecked(:, 1), unchecked(:, 2)), &
                        'unchecked radii: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(grid_sa(:, 1), grid_sa(:, 2)), &
                        'one-shot grid: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(radii_sa(:, 1), radii_sa(:, 2)), &
                        'one-shot radii: short == padded, ' // trim(tag))
                call assert_true(bits_equal_f(drs_sa(:, 1), drs_sa(:, 2)), &
                        'one-shot derivatives: short == padded, ' // trim(tag))
            end do

            call assert_true(all_valid, 'every padding call returned SHAPE_VALID')
        end do

        call node_set_free_s(nodes)
        call cache_free_s(cache)
    end subroutine run_padding_s

end program beta_param_equivalence_test
