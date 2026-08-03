!> Contract family 1: an incremental compute (warm cache, one parameter
!! perturbed) must reproduce a cold compute (fresh cache, same init arguments,
!! same final parameter vector) bit for bit. Every comparison here is on the
!! IEEE754 bit patterns; no tolerances appear anywhere in this suite.
program beta_param_bitwise_test
    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_init_shared_s, cache_resolve_shape_s, &
            cache_radius_and_derivative_s, cache_radius_grid_unchecked_s, &
            cache_node_radius_and_derivative_s, &
            tables_t, tables_init_s, tables_free_s, &
            node_set_t, node_set_build_s, node_set_free_s, &
            SHAPE_VALID

    implicit none

    integer(kind = ik), parameter :: N_THETAS = 32_ik
    integer(kind = ik), parameter :: N_PARAMS = 6_ik
    integer(kind = ik), parameter :: N_NODES = 8_ik
    real(kind = rk),    parameter :: DELTA = 0.013_rk

    !> IEEE754 -0.0; built from bits so that -fno-signed-zeros cannot fold it
    !! into +0.0 the way a `-0.0_rk` literal could.
    integer(kind = ikl), parameter :: NEG_ZERO_BITS = int(z'8000000000000000', kind = ikl)

    real(kind = rk) :: thetas(N_THETAS), node_thetas(N_NODES)
    real(kind = rk) :: base_params(N_PARAMS)
    integer(kind = ik) :: i

    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do
    do i = 1_ik, N_NODES
        node_thetas(i) = real(2_ik * i - 1_ik, rk) * PI_C / real(2_ik * N_NODES, rk)
    end do
    base_params = [0.02_rk, 0.15_rk, 0.08_rk, -0.05_rk, 0.03_rk, 0.01_rk]

    call run_grid_sweep_s()
    call run_node_case_s()
    call run_unchecked_case_s()
    call run_negzero_case_s()

    call test_summary()

contains

    !> Every regime x every perturbed parameter: incremental == cold, bitwise,
    !! for the radii, the derivatives and all four resolve scalars.
    subroutine run_grid_sweep_s()
        type(cache_t) :: warm, cold
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: radii_inc(N_THETAS), drs_inc(N_THETAS)
        real(kind = rk) :: radii_cold(N_THETAS), drs_cold(N_THETAS)
        real(kind = rk) :: scalars_inc(4), scalars_cold(4)
        character(len = 4), parameter :: SCALAR_NAMES(4) = ['b10 ', 'rn  ', 'rs  ', 'vf  ']
        character(len = 48) :: label
        logical :: conserve_volume, apply_com, all_valid
        integer(kind = ik) :: regime, j, k, status

        all_valid = .true.

        do regime = 1_ik, 4_ik
            conserve_volume = (regime == 2_ik .or. regime == 4_ik)
            apply_com = (regime == 3_ik .or. regime == 4_ik)

            do j = 1_ik, N_PARAMS
                write(label, '(A,I0,A,I0)') 'regime ', regime, ', param ', j

                ! Warm the cache on the base vector, then perturb one parameter.
                call cache_init_s(warm, N_PARAMS, thetas, conserve_volume, apply_com, status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                params = base_params
                call cache_radius_and_derivative_s(warm, params, radii_inc, drs_inc, status)
                all_valid = all_valid .and. (status == SHAPE_VALID)

                params(j) = params(j) + DELTA
                call cache_radius_and_derivative_s(warm, params, radii_inc, drs_inc, status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call cache_resolve_shape_s(warm, params, scalars_inc(1), scalars_inc(2), &
                        scalars_inc(3), scalars_inc(4), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)

                ! Cold reference: identical init arguments, perturbed vector only.
                call cache_init_s(cold, N_PARAMS, thetas, conserve_volume, apply_com, status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call cache_radius_and_derivative_s(cold, params, radii_cold, drs_cold, status)
                all_valid = all_valid .and. (status == SHAPE_VALID)
                call cache_resolve_shape_s(cold, params, scalars_cold(1), scalars_cold(2), &
                        scalars_cold(3), scalars_cold(4), status)
                all_valid = all_valid .and. (status == SHAPE_VALID)

                call assert_true(bits_equal_f(radii_inc, radii_cold), &
                        'radii bitwise incremental == cold, ' // trim(label))
                call assert_true(bits_equal_f(drs_inc, drs_cold), &
                        'derivatives bitwise incremental == cold, ' // trim(label))
                do k = 1_ik, 4_ik
                    call assert_true(bits_equal_f(scalars_inc(k:k), scalars_cold(k:k)), &
                            trim(SCALAR_NAMES(k)) // ' bitwise incremental == cold, ' // trim(label))
                end do

                call cache_free_s(warm)
                call cache_free_s(cold)
            end do
        end do

        call assert_true(all_valid, 'every sweep call returned SHAPE_VALID')
    end subroutine run_grid_sweep_s

    !> The node-set path over shared tables obeys the same contract.
    subroutine run_node_case_s()
        type(tables_t), target :: tables
        type(node_set_t) :: nodes
        type(cache_t) :: warm, cold
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: radii_inc(N_NODES), drs_inc(N_NODES)
        real(kind = rk) :: radii_cold(N_NODES), drs_cold(N_NODES)
        integer(kind = ik) :: status

        call tables_init_s(tables, N_PARAMS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'shared tables init')
        call node_set_build_s(nodes, tables, node_thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'node set build')

        call cache_init_shared_s(warm, tables, N_PARAMS, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'warm shared cache init')
        params = base_params
        call cache_node_radius_and_derivative_s(warm, nodes, params, radii_inc, drs_inc, status)
        call assert_int_eq(status, SHAPE_VALID, 'node warm-up call')

        params(3) = params(3) + DELTA
        call cache_node_radius_and_derivative_s(warm, nodes, params, radii_inc, drs_inc, status)
        call assert_int_eq(status, SHAPE_VALID, 'node incremental call')

        call cache_init_shared_s(cold, tables, N_PARAMS, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'cold shared cache init')
        call cache_node_radius_and_derivative_s(cold, nodes, params, radii_cold, drs_cold, status)
        call assert_int_eq(status, SHAPE_VALID, 'node cold call')

        call assert_true(bits_equal_f(radii_inc, radii_cold), &
                'node radii bitwise incremental == cold')
        call assert_true(bits_equal_f(drs_inc, drs_cold), &
                'node derivatives bitwise incremental == cold')

        call cache_free_s(warm)
        call cache_free_s(cold)
        call node_set_free_s(nodes)
        call tables_free_s(tables)
    end subroutine run_node_case_s

    !> The unchecked path reuses the same cached state: warm checked call first,
    !! then unchecked must match a cold cache's unchecked call bit for bit,
    !! both on the warm vector and after a perturbation.
    subroutine run_unchecked_case_s()
        type(cache_t) :: warm, cold
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: radii_inc(N_THETAS), drs_inc(N_THETAS), radii_cold(N_THETAS)
        integer(kind = ik) :: status

        call cache_init_s(warm, N_PARAMS, thetas, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'unchecked warm cache init')
        call cache_init_s(cold, N_PARAMS, thetas, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'unchecked cold cache init')

        params = base_params
        call cache_radius_and_derivative_s(warm, params, radii_inc, drs_inc, status)
        call assert_int_eq(status, SHAPE_VALID, 'unchecked warm-up checked call')

        call cache_radius_grid_unchecked_s(warm, params, radii_inc, status)
        call assert_int_eq(status, SHAPE_VALID, 'unchecked call on the warm vector')
        call cache_radius_grid_unchecked_s(cold, params, radii_cold, status)
        call assert_int_eq(status, SHAPE_VALID, 'cold unchecked call')
        call assert_true(bits_equal_f(radii_inc, radii_cold), &
                'unchecked radii bitwise warm == cold, same vector')

        call cache_free_s(cold)

        params(5) = params(5) + DELTA
        call cache_radius_grid_unchecked_s(warm, params, radii_inc, status)
        call assert_int_eq(status, SHAPE_VALID, 'unchecked incremental call')
        call cache_init_s(cold, N_PARAMS, thetas, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'unchecked cold cache re-init')
        call cache_radius_grid_unchecked_s(cold, params, radii_cold, status)
        call assert_int_eq(status, SHAPE_VALID, 'cold unchecked call, perturbed vector')
        call assert_true(bits_equal_f(radii_inc, radii_cold), &
                'unchecked radii bitwise incremental == cold')

        call cache_free_s(warm)
        call cache_free_s(cold)
    end subroutine run_unchecked_case_s

    !> Bit semantics: replacing +0.0 by -0.0 is a parameter change even though
    !! the two compare equal as floating point. The incremental result must
    !! match a cold compute on (-0.0, 0.2).
    subroutine run_negzero_case_s()
        type(cache_t) :: warm, cold
        !> volatile: under -fno-signed-zeros (Release -ffast-math) the compiler
        !! may elide a store of -0.0 over a location it knows holds +0.0, so the
        !! bit pattern would never reach the engine; volatile forces the store.
        real(kind = rk), volatile :: pv(2), pv_cold(2)
        real(kind = rk) :: radii_inc(N_THETAS), drs_inc(N_THETAS)
        real(kind = rk) :: radii_cold(N_THETAS), drs_cold(N_THETAS)
        integer(kind = ik) :: status

        call cache_init_s(warm, 2_ik, thetas, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'negative-zero warm cache init')

        pv(1) = 0.0_rk
        pv(2) = 0.2_rk
        call cache_radius_and_derivative_s(warm, pv, radii_inc, drs_inc, status)
        call assert_int_eq(status, SHAPE_VALID, 'negative-zero warm-up call')

        pv(1) = transfer(NEG_ZERO_BITS, 1.0_rk)
        call assert_true(transfer(pv(1), 0_ikl) == NEG_ZERO_BITS, &
                'the -0.0 store reached memory')
        call cache_radius_and_derivative_s(warm, pv, radii_inc, drs_inc, status)
        call assert_int_eq(status, SHAPE_VALID, 'negative-zero incremental call')

        call cache_init_s(cold, 2_ik, thetas, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'negative-zero cold cache init')
        pv_cold(1) = transfer(NEG_ZERO_BITS, 1.0_rk)
        pv_cold(2) = 0.2_rk
        call cache_radius_and_derivative_s(cold, pv_cold, radii_cold, drs_cold, status)
        call assert_int_eq(status, SHAPE_VALID, 'negative-zero cold call')

        call assert_true(bits_equal_f(radii_inc, radii_cold), &
                'negative-zero radii bitwise incremental == cold')
        call assert_true(bits_equal_f(drs_inc, drs_cold), &
                'negative-zero derivatives bitwise incremental == cold')

        call cache_free_s(warm)
        call cache_free_s(cold)
    end subroutine run_negzero_case_s

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

end program beta_param_bitwise_test
