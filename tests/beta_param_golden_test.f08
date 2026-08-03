!> Golden baseline for the 3.0.0 pipeline (Newton COM, in-library volume).
!!
!! Six representative shapes x the two extreme regimes
!! (conserve_volume, apply_com) = (F,F) and (T,T). The literals were captured
!! from this library by `golden_capture` and pin cross-version drift: any
!! change to the resolve/volume/render path that moves a number by more than
!! 1e-15 relative shows up here.
!!
!! Bit-level identity between call paths is the bitwise suite's job; this suite
!! answers "does the library still produce the same physics as the baseline".
!! Regenerate with `./build/golden_capture` — and only with a reason.
program beta_param_golden_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_int_eq, assert_close, test_summary
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_resolve_shape_s, cache_radius_grid_s, SHAPE_VALID

    implicit none

    integer(kind = ik), parameter :: N_THETAS = 16_ik

    real(kind = rk)    :: thetas(N_THETAS)
    real(kind = rk)    :: radii(N_THETAS)
    real(kind = rk)    :: corrected_beta10, volume_factor
    integer(kind = ik) :: i

    ! Same open uniform grid the capture used: no node on a pole.
    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do

    ! S1: sphere
    call run_case_s([0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], .false., .false., 'S1 FF', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.00000000000000000E+00_rk, 1.0e-15_rk, 'S1 FF radii(1)')
    call assert_close(radii(16), 1.00000000000000000E+00_rk, 1.0e-15_rk, 'S1 FF radii(16)')
    call assert_close(corrected_beta10, 0.00000000000000000E+00_rk, 1.0e-15_rk, 'S1 FF corrected_beta10')
    call assert_close(volume_factor, 1.00000000000000000E+00_rk, 1.0e-15_rk, 'S1 FF volume_factor')
    call run_case_s([0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], .true., .true., 'S1 TT', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.00000000000002487E+00_rk, 1.0e-15_rk, 'S1 TT radii(1)')
    call assert_close(radii(16), 1.00000000000002487E+00_rk, 1.0e-15_rk, 'S1 TT radii(16)')
    call assert_close(corrected_beta10, 0.00000000000000000E+00_rk, 1.0e-15_rk, 'S1 TT corrected_beta10')
    call assert_close(volume_factor, 1.00000000000002487E+00_rk, 1.0e-15_rk, 'S1 TT volume_factor')

    ! S2: prolate
    call run_case_s([0.0_rk, 0.3_rk, 0.0_rk, 0.0_rk], .false., .false., 'S2 FF', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.17965097898173399E+00_rk, 1.0e-15_rk, 'S2 FF radii(1)')
    call assert_close(radii(16), 1.17965097898173399E+00_rk, 1.0e-15_rk, 'S2 FF radii(16)')
    call assert_close(corrected_beta10, 0.00000000000000000E+00_rk, 1.0e-15_rk, 'S2 FF corrected_beta10')
    call assert_close(volume_factor, 1.00000000000000000E+00_rk, 1.0e-15_rk, 'S2 FF volume_factor')
    call run_case_s([0.0_rk, 0.3_rk, 0.0_rk, 0.0_rk], .true., .true., 'S2 TT', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.17117341036512990E+00_rk, 1.0e-15_rk, 'S2 TT radii(1)')
    call assert_close(radii(16), 1.17117341036512990E+00_rk, 1.0e-15_rk, 'S2 TT radii(16)')
    call assert_close(corrected_beta10, 0.00000000000000000E+00_rk, 1.0e-15_rk, 'S2 TT corrected_beta10')
    call assert_close(volume_factor, 9.92813494187982815E-01_rk, 1.0e-15_rk, 'S2 TT volume_factor')

    ! S3: asymmetric
    call run_case_s([0.05_rk, 0.25_rk, 0.12_rk, 0.03_rk], .false., .false., 'S3 FF', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.27555852826677296E+00_rk, 1.0e-15_rk, 'S3 FF radii(1)')
    call assert_close(radii(16), 1.06631792857040297E+00_rk, 1.0e-15_rk, 'S3 FF radii(16)')
    call assert_close(corrected_beta10, 5.00000000000000028E-02_rk, 1.0e-15_rk, 'S3 FF corrected_beta10')
    call assert_close(volume_factor, 1.00000000000000000E+00_rk, 1.0e-15_rk, 'S3 FF volume_factor')
    call run_case_s([0.05_rk, 0.25_rk, 0.12_rk, 0.03_rk], .true., .true., 'S3 TT', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.23265776951246386E+00_rk, 1.0e-15_rk, 'S3 TT radii(1)')
    call assert_close(radii(16), 1.09448160290258545E+00_rk, 1.0e-15_rk, 'S3 TT radii(16)')
    call assert_close(corrected_beta10, -2.30708623724193915E-02_rk, 1.0e-15_rk, 'S3 TT corrected_beta10')
    call assert_close(volume_factor, 9.93707146942315656E-01_rk, 1.0e-15_rk, 'S3 TT volume_factor')

    ! S4: six-parameter mixed
    call run_case_s([0.02_rk, 0.15_rk, 0.08_rk, -0.05_rk, 0.03_rk, 0.01_rk], .false., .false., 'S4 FF', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.14593727573622317E+00_rk, 1.0e-15_rk, 'S4 FF radii(1)')
    call assert_close(radii(16), 9.76617097509943410E-01_rk, 1.0e-15_rk, 'S4 FF radii(16)')
    call assert_close(corrected_beta10, 2.00000000000000004E-02_rk, 1.0e-15_rk, 'S4 FF corrected_beta10')
    call assert_close(volume_factor, 1.00000000000000000E+00_rk, 1.0e-15_rk, 'S4 FF volume_factor')
    call run_case_s([0.02_rk, 0.15_rk, 0.08_rk, -0.05_rk, 0.03_rk, 0.01_rk], .true., .true., 'S4 TT', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.13105543631066063E+00_rk, 1.0e-15_rk, 'S4 TT radii(1)')
    call assert_close(radii(16), 9.86012148345387418E-01_rk, 1.0e-15_rk, 'S4 TT radii(16)')
    call assert_close(corrected_beta10, -4.88218027607693547E-03_rk, 1.0e-15_rk, 'S4 TT corrected_beta10')
    call assert_close(volume_factor, 9.97415006814771465E-01_rk, 1.0e-15_rk, 'S4 TT volume_factor')

    ! S5: eight-parameter mixed
    call run_case_s([0.02_rk, 0.2_rk, 0.1_rk, -0.04_rk, 0.03_rk, -0.02_rk, &
            0.01_rk, 0.005_rk], .false., .false., 'S5 FF', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.18492682251660542E+00_rk, 1.0e-15_rk, 'S5 FF radii(1)')
    call assert_close(radii(16), 9.76162497710678645E-01_rk, 1.0e-15_rk, 'S5 FF radii(16)')
    call assert_close(corrected_beta10, 2.00000000000000004E-02_rk, 1.0e-15_rk, 'S5 FF corrected_beta10')
    call assert_close(volume_factor, 1.00000000000000000E+00_rk, 1.0e-15_rk, 'S5 FF volume_factor')
    call run_case_s([0.02_rk, 0.2_rk, 0.1_rk, -0.04_rk, 0.03_rk, -0.02_rk, &
            0.01_rk, 0.005_rk], .true., .true., 'S5 TT', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.16565815802888451E+00_rk, 1.0e-15_rk, 'S5 TT radii(1)')
    call assert_close(radii(16), 9.86262088836088791E-01_rk, 1.0e-15_rk, 'S5 TT radii(16)')
    call assert_close(corrected_beta10, -9.77813399136218467E-03_rk, 1.0e-15_rk, 'S5 TT corrected_beta10')
    call assert_close(volume_factor, 9.95757198336741589E-01_rk, 1.0e-15_rk, 'S5 TT volume_factor')

    ! S6: strongly deformed, near the validity edge
    call run_case_s([0.0_rk, 0.35_rk, 0.25_rk, 0.1_rk], .false., .false., 'S6 FF', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.44828587213220206E+00_rk, 1.0e-15_rk, 'S6 FF radii(1)')
    call assert_close(radii(16), 1.11242694060612979E+00_rk, 1.0e-15_rk, 'S6 FF radii(16)')
    call assert_close(corrected_beta10, 0.00000000000000000E+00_rk, 1.0e-15_rk, 'S6 FF corrected_beta10')
    call assert_close(volume_factor, 1.00000000000000000E+00_rk, 1.0e-15_rk, 'S6 FF volume_factor')
    call run_case_s([0.0_rk, 0.35_rk, 0.25_rk, 0.1_rk], .true., .true., 'S6 TT', &
            radii, corrected_beta10, volume_factor)
    call assert_close(radii(1), 1.38953768578917458E+00_rk, 1.0e-15_rk, 'S6 TT radii(1)')
    call assert_close(radii(16), 1.13018457995286492E+00_rk, 1.0e-15_rk, 'S6 TT radii(16)')
    call assert_close(corrected_beta10, -7.52542571777818359E-02_rk, 1.0e-15_rk, 'S6 TT corrected_beta10')
    call assert_close(volume_factor, 9.83992524740617602E-01_rk, 1.0e-15_rk, 'S6 TT volume_factor')

    call test_summary()

contains

    !> One (shape, regime) case, through exactly the calls the capture used:
    !! init, resolve, radius grid — same order, same arguments.
    subroutine run_case_s(params, conserve_volume, apply_com, label, &
            out_radii, out_corrected_beta10, out_volume_factor)
        real(kind = rk),    intent(in)  :: params(:)
        logical,            intent(in)  :: conserve_volume, apply_com
        character(len = *), intent(in)  :: label
        real(kind = rk),    intent(out) :: out_radii(:)
        real(kind = rk),    intent(out) :: out_corrected_beta10, out_volume_factor

        type(cache_t)      :: cache
        real(kind = rk)    :: r_north, r_south
        integer(kind = ik) :: status

        call cache_init_s(cache, size(params, kind = ik), thetas, conserve_volume, &
                apply_com, status)
        call assert_int_eq(status, SHAPE_VALID, label // ' init')

        call cache_resolve_shape_s(cache, params, out_corrected_beta10, r_north, &
                r_south, out_volume_factor, status)
        call assert_int_eq(status, SHAPE_VALID, label // ' resolve')

        call cache_radius_grid_s(cache, params, out_radii, status)
        call assert_int_eq(status, SHAPE_VALID, label // ' radius grid')

        call cache_free_s(cache)
    end subroutine run_case_s

end program beta_param_golden_test
