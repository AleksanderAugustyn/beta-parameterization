!> Prints golden values for the four regression shapes as paste-ready Fortran.
!! Run AFTER all Task-3 fixes; output goes into beta_param_golden_test.f08.
program golden_capture

    use precision_utilities_mod, only: ik, rk
    use beta_parameterization_mod, only: cache_t, LEGENDRE_VALID
    use test_utils_mod, only: assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: IDX(7) = [1_ik, 31_ik, 61_ik, 91_ik, 121_ik, 151_ik, 181_ik]
    type(cache_t)        :: cache
    integer(kind = ik)   :: code
    character(len = 256) :: message

    call cache%init(8_ik, 181_ik, code, message)
    call assert_int_eq(code, LEGENDRE_VALID, 'capture: cache init')

    call capture('G1', [0.0_rk, 0.215_rk, 0.0_rk, 0.095_rk], .false.)
    call capture('G2', [0.0_rk, 0.85_rk, 0.35_rk, 0.18_rk, 0.05_rk, 0.02_rk], .true.)
    call capture('G3', [0.0_rk, -0.35_rk, 0.0_rk, 0.05_rk], .false.)
    call capture('G4', [0.0_rk, 0.40_rk, 0.20_rk, 0.10_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk], .true.)
    call test_summary()

contains

    subroutine capture(name, params, with_shift)
        character(len = *), intent(in) :: name
        real(kind = rk),    intent(in) :: params(:)
        logical,            intent(in) :: with_shift
        real(kind = rk)      :: radii(181), corrected
        integer(kind = ik)   :: c, i
        character(len = 256) :: msg

        if (with_shift) then
            call cache%compute_radius_grid_with_com_shift(params, radii, corrected, c, msg)
        else
            call cache%compute_radius_grid(params, radii, c, msg)
            corrected = 0.0_rk
        end if
        call assert_int_eq(c, LEGENDRE_VALID, name // ': valid')

        write(*, '(A,A,A)') 'real(kind = rk), parameter :: ', name, '_EXPECTED(7) = [ &'
        do i = 1_ik, 6_ik
            write(*, '(A,ES24.16E3,A)') '        ', radii(IDX(i)), '_rk, &'
        end do
        write(*, '(A,ES24.16E3,A)') '        ', radii(IDX(7)), '_rk]'
        if (with_shift) then
            write(*, '(A,A,A,ES24.16E3,A)') 'real(kind = rk), parameter :: ', name, &
                    '_CORRECTED_B10 = ', corrected, '_rk'
        end if
    end subroutine capture

end program golden_capture
