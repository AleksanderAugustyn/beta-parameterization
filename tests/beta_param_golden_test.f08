!> Locked regression goldens for four representative shapes.
!! Values captured by tests/golden_capture.f08 at commit 8bfe973 (post Task-3 fixes).
program beta_param_golden_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use beta_parameterization_mod, only: cache_t, node_set_t, LEGENDRE_VALID
    use test_utils_mod, only: assert_int_eq, assert_close, test_summary

    implicit none

    integer(kind = ik), parameter :: IDX(7) = [1_ik, 31_ik, 61_ik, 91_ik, 121_ik, 151_ik, 181_ik]
    real(kind = rk),    parameter :: TOL = 1.0e-15_rk

    real(kind = rk), parameter :: G1_EXPECTED(7) = [ &
             1.2160153887141389E+000_rk, &
             1.0866457882160419E+000_rk, &
             9.5980794102974309E-001_rk, &
             9.6233969434154143E-001_rk, &
             9.5980794102974309E-001_rk, &
             1.0866457882160419E+000_rk, &
             1.2160153887141389E+000_rk]
    real(kind = rk), parameter :: G2_EXPECTED(7) = [ &
             1.9145931553115427E+000_rk, &
             1.3169048809475408E+000_rk, &
             7.3431444834952730E-001_rk, &
             7.8268444464280551E-001_rk, &
             1.0567285473303596E+000_rk, &
             1.3452258418380265E+000_rk, &
             1.5030848311140956E+000_rk]
    real(kind = rk), parameter :: G2_CORRECTED_B10 = -2.0926908316026530E-001_rk
    real(kind = rk), parameter :: G3_EXPECTED(7) = [ &
             8.2154012308931779E-001_rk, &
             8.6300792970435258E-001_rk, &
             1.0153653080975249E+000_rk, &
             1.1262548798756626E+000_rk, &
             1.0153653080975251E+000_rk, &
             8.6300792970435258E-001_rk, &
             8.2154012308931779E-001_rk]
    real(kind = rk), parameter :: G4_EXPECTED(7) = [ &
             1.5346514262194999E+000_rk, &
             1.1529734346542979E+000_rk, &
             8.7376808780044701E-001_rk, &
             9.0081230258281431E-001_rk, &
             1.0265221633106583E+000_rk, &
             1.1472278889699292E+000_rk, &
             1.1915473089293442E+000_rk]
    real(kind = rk), parameter :: G4_CORRECTED_B10 = -7.2500830051600282E-002_rk

    ! Node-set goldens (captured by golden_capture, v2.2.0 node-set API):
    ! layout [corrected_beta10, r_north, r_south, R(1:3), dR(1:3)] at
    ! thetas = [pi/8, pi/2, 7pi/8].
    real(kind = rk), parameter :: G1_NODE_SET_EXPECTED(9) = [ &
             0.0000000000000000E+000_rk, &
             1.2160153887141389E+000_rk, &
             1.2160153887141389E+000_rk, &
             1.1348983254477432E+000_rk, &
             9.6233969434154143E-001_rk, &
             1.1348983254477432E+000_rk, &
            -3.5524427545959703E-001_rk, &
             1.2009039485575869E-017_rk, &
             3.5524427545959714E-001_rk]
    real(kind = rk), parameter :: G2_NODE_SET_EXPECTED(9) = [ &
            -2.0926908316026530E-001_rk, &
             1.9145931553115429E+000_rk, &
             1.5030848311140959E+000_rk, &
             1.5366433412487048E+000_rk, &
             7.8268444464280551E-001_rk, &
             1.4071871239950646E+000_rk, &
            -1.6309751173882105E+000_rk, &
             4.0637180707528003E-001_rk, &
             4.3722024621518629E-001_rk]
    real(kind = rk), parameter :: G3_NODE_SET_EXPECTED(9) = [ &
             0.0000000000000000E+000_rk, &
             8.2154012308931779E-001_rk, &
             8.2154012308931779E-001_rk, &
             8.4302397766917225E-001_rk, &
             1.1262548798756626E+000_rk, &
             8.4302397766917225E-001_rk, &
             1.2290351730043024E-001_rk, &
             5.9988033154643641E-017_rk, &
            -1.2290351730043028E-001_rk]
    real(kind = rk), parameter :: G4_NODE_SET_EXPECTED(9) = [ &
            -7.2500830051600282E-002_rk, &
             1.5346514262195001E+000_rk, &
             1.1915473089293445E+000_rk, &
             1.2820856414035799E+000_rk, &
             9.0081230258281431E-001_rk, &
             1.1656292020106376E+000_rk, &
            -1.0077766770385752E+000_rk, &
             1.9551664231146185E-001_rk, &
             1.2297403623871962E-001_rk]

    type(cache_t)        :: cache
    integer(kind = ik)   :: code
    character(len = 256) :: message

    call cache%init(8_ik, 181_ik, code, message)
    call assert_int_eq(code, LEGENDRE_VALID, 'golden: cache init')

    call check('G1', [0.0_rk, 0.215_rk, 0.0_rk, 0.095_rk], .false., G1_EXPECTED, 0.0_rk)
    call check('G2', [0.0_rk, 0.85_rk, 0.35_rk, 0.18_rk, 0.05_rk, 0.02_rk], .true., &
            G2_EXPECTED, G2_CORRECTED_B10)
    call check('G3', [0.0_rk, -0.35_rk, 0.0_rk, 0.05_rk], .false., G3_EXPECTED, 0.0_rk)
    call check('G4', [0.0_rk, 0.40_rk, 0.20_rk, 0.10_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk], &
            .true., G4_EXPECTED, G4_CORRECTED_B10)

    call check_node_set('G1', [0.0_rk, 0.215_rk, 0.0_rk, 0.095_rk], G1_NODE_SET_EXPECTED)
    call check_node_set('G2', [0.0_rk, 0.85_rk, 0.35_rk, 0.18_rk, 0.05_rk, 0.02_rk], &
            G2_NODE_SET_EXPECTED)
    call check_node_set('G3', [0.0_rk, -0.35_rk, 0.0_rk, 0.05_rk], G3_NODE_SET_EXPECTED)
    call check_node_set('G4', [0.0_rk, 0.40_rk, 0.20_rk, 0.10_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk], &
            G4_NODE_SET_EXPECTED)
    call test_summary()

contains

    subroutine check(name, params, with_shift, expected, expected_b10)
        character(len = *), intent(in) :: name
        real(kind = rk),    intent(in) :: params(:), expected(7), expected_b10
        logical,            intent(in) :: with_shift
        real(kind = rk)      :: radii(181), corrected
        integer(kind = ik)   :: c, i
        character(len = 256) :: msg

        if (with_shift) then
            call cache%compute_radius_grid_with_com_shift(params, radii, corrected, c, msg)
        else
            call cache%compute_radius_grid(params, radii, c, msg)
        end if
        call assert_int_eq(c, LEGENDRE_VALID, name // ': valid')
        do i = 1_ik, 7_ik
            call assert_close(radii(IDX(i)), expected(i), TOL, name // ': golden radius')
        end do
        if (with_shift) call assert_close(corrected, expected_b10, TOL, name // ': golden corrected beta10')
    end subroutine check

    !> Locked node-set goldens: resolve + evaluate at thetas = [pi/8, pi/2, 7pi/8].
    subroutine check_node_set(name, params, expected)
        character(len = *), intent(in) :: name
        real(kind = rk),    intent(in) :: params(:), expected(9)
        type(node_set_t)     :: node_set
        real(kind = rk)      :: thetas(3), beta_con(8), radii(3), dr(3)
        real(kind = rk)      :: corrected_beta10, r_north, r_south
        integer(kind = ik)   :: c, i
        character(len = 256) :: msg

        thetas = [PI_C / 8.0_rk, PI_C / 2.0_rk, 7.0_rk * PI_C / 8.0_rk]
        call cache%build_node_set(thetas, node_set, c, msg)
        call assert_int_eq(c, LEGENDRE_VALID, name // ': node set built')
        call cache%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, c, msg)
        call assert_int_eq(c, LEGENDRE_VALID, name // ': resolved')
        call cache%compute_radius_and_derivative(beta_con, node_set, radii, dr, c, msg)
        call assert_int_eq(c, LEGENDRE_VALID, name // ': evaluated')

        call assert_close(corrected_beta10, expected(1), TOL, name // ': golden node-set corrected beta10')
        call assert_close(r_north, expected(2), TOL, name // ': golden r_north')
        call assert_close(r_south, expected(3), TOL, name // ': golden r_south')
        do i = 1_ik, 3_ik
            call assert_close(radii(i), expected(3 + i), TOL, name // ': golden node-set R')
            call assert_close(dr(i), expected(6 + i), TOL, name // ': golden node-set dR/dtheta')
        end do
    end subroutine check_node_set

end program beta_param_golden_test
