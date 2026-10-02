!> Capture tool for the 3.0.0 golden baseline.
!!
!! Prints paste-ready `assert_close` lines for six representative shapes in the
!! two extreme regimes — (conserve_volume, apply_com) = (F,F) and (T,T). The
!! output goes verbatim into `beta_param_golden_test.f08`, which recomputes the
!! same quantities through the same public calls and compares.
!!
!! Not a test: no add_test, it only writes to stdout. It stops hard if any
!! captured run reports a status other than SHAPE_VALID — a golden may only be
!! taken from a shape the library accepts.
program golden_capture

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_resolve_shape_s, cache_radius_grid_s, SHAPE_VALID

    implicit none

    integer(kind = ik), parameter :: N_THETAS = 16_ik

    real(kind = rk) :: thetas(N_THETAS)
    integer(kind = ik) :: i

    ! Open uniform grid: no node sits on a pole, so the poles stay the job of
    ! resolve and never leak into the grid goldens.
    do i = 1_ik, N_THETAS
        thetas(i) = real(i, rk) * PI_C / real(N_THETAS + 1_ik, rk)
    end do

    call capture_shape_s(1_ik, [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk])
    call capture_shape_s(2_ik, [0.0_rk, 0.3_rk, 0.0_rk, 0.0_rk])
    call capture_shape_s(3_ik, [0.05_rk, 0.25_rk, 0.12_rk, 0.03_rk])
    call capture_shape_s(4_ik, [0.02_rk, 0.15_rk, 0.08_rk, -0.05_rk, 0.03_rk, 0.01_rk])
    call capture_shape_s(5_ik, [0.02_rk, 0.2_rk, 0.1_rk, -0.04_rk, 0.03_rk, -0.02_rk, &
            0.01_rk, 0.005_rk])
    call capture_shape_s(6_ik, [0.0_rk, 0.35_rk, 0.25_rk, 0.1_rk])

contains

    !> Both regimes for one shape.
    subroutine capture_shape_s(id, params)
        integer(kind = ik), intent(in) :: id
        real(kind = rk),    intent(in) :: params(:)

        call capture_regime_s(id, params, .false., .false., 'FF')
        call capture_regime_s(id, params, .true., .true., 'TT')
    end subroutine capture_shape_s

    !> One (shape, regime) block: init, resolve, radius grid, then four literals.
    subroutine capture_regime_s(id, params, conserve_volume, apply_com, tag)
        integer(kind = ik), intent(in) :: id
        real(kind = rk),    intent(in) :: params(:)
        logical,            intent(in) :: conserve_volume, apply_com
        character(len = *), intent(in) :: tag

        type(cache_t)      :: cache
        real(kind = rk)    :: radii(N_THETAS)
        real(kind = rk)    :: corrected_beta10, r_north, r_south, volume_factor
        integer(kind = ik) :: status
        character(len = 16) :: prefix

        write(prefix, '(A,I0,A,A)') 'S', id, ' ', tag

        call cache_init_s(cache, size(params, kind = ik), thetas, status)
        call require_valid_s(prefix, 'init', status)

        call cache_resolve_shape_s(cache, params, conserve_volume, apply_com, &
                corrected_beta10, r_north, r_south, volume_factor, status)
        call require_valid_s(prefix, 'resolve', status)

        call cache_radius_grid_s(cache, params, conserve_volume, apply_com, radii, status)
        call require_valid_s(prefix, 'radius grid', status)

        write(*, '(A,A,A,I0,A,L1,A,L1,A)') '    ! ', trim(prefix), ': n_params = ', &
                size(params, kind = ik), ', conserve_volume = ', conserve_volume, &
                ', apply_com = ', apply_com, ', status = 0'
        call emit_s('radii(1)',         radii(1),         trim(prefix))
        call emit_s('radii(16)',        radii(N_THETAS),  trim(prefix))
        call emit_s('corrected_beta10', corrected_beta10, trim(prefix))
        call emit_s('volume_factor',    volume_factor,    trim(prefix))

        call cache_free_s(cache)
    end subroutine capture_regime_s

    !> One paste-ready assertion, value at full double round-trip precision.
    subroutine emit_s(expr, value, prefix)
        character(len = *), intent(in) :: expr
        real(kind = rk),    intent(in) :: value
        character(len = *), intent(in) :: prefix

        character(len = 32) :: buffer

        write(buffer, '(ES25.17)') value
        write(*, '(A)') '    call assert_close(' // expr // ', ' // &
                trim(adjustl(buffer)) // '_rk, 1.0e-15_rk, ''' // &
                prefix // ' ' // expr // ''')'
    end subroutine emit_s

    !> A golden is only valid if the library accepted the shape.
    subroutine require_valid_s(prefix, stage, status)
        character(len = *), intent(in) :: prefix, stage
        integer(kind = ik), intent(in) :: status

        if (status /= SHAPE_VALID) then
            write(*, '(A,A,A,A,A,I0)') '    ! CAPTURE ABORTED: ', prefix, ' ', stage, &
                    ' returned status ', status
            error stop 1
        end if
    end subroutine require_valid_s

end program golden_capture
