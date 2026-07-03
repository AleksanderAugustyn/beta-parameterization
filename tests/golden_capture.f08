!> Prints golden values for the four regression shapes as paste-ready Fortran.
!! Run AFTER all Task-3 fixes; output goes into beta_param_golden_test.f08.
program golden_capture

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use beta_parameterization_mod, only: cache_t, node_set_t, LEGENDRE_VALID
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

        call capture_node_set(name, params)
    end subroutine capture

    !> Node-set goldens: resolve + evaluate R and dR/dtheta at three fixed thetas.
    !! Layout: [corrected_beta10, r_north, r_south, R(1:3), dR(1:3)].
    subroutine capture_node_set(name, params)
        character(len = *), intent(in) :: name
        real(kind = rk),    intent(in) :: params(:)
        type(node_set_t)     :: node_set
        real(kind = rk)      :: thetas(3), beta_con(8), radii(3), dr(3)
        real(kind = rk)      :: corrected_beta10, r_north, r_south, values(9)
        integer(kind = ik)   :: c, i
        character(len = 256) :: msg

        thetas = [PI_C / 8.0_rk, PI_C / 2.0_rk, 7.0_rk * PI_C / 8.0_rk]
        call cache%build_node_set(thetas, node_set, c, msg)
        call assert_int_eq(c, LEGENDRE_VALID, name // ': node set built')
        call cache%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, c, msg)
        call assert_int_eq(c, LEGENDRE_VALID, name // ': resolved')
        call cache%compute_radius_and_derivative(beta_con, node_set, radii, dr, c, msg)
        call assert_int_eq(c, LEGENDRE_VALID, name // ': evaluated')

        values = [corrected_beta10, r_north, r_south, radii, dr]
        write(*, '(A,A,A)') 'real(kind = rk), parameter :: ', name, '_NODE_SET_EXPECTED(9) = [ &'
        do i = 1_ik, 8_ik
            write(*, '(A,ES24.16E3,A)') '        ', values(i), '_rk, &'
        end do
        write(*, '(A,ES24.16E3,A)') '        ', values(9), '_rk]'
    end subroutine capture_node_set

end program golden_capture
