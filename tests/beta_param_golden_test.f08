!> Locked regression goldens for four representative shapes.
!! Values captured by tests/golden_capture.f08 at commit 8bfe973 (post Task-3 fixes).
program beta_param_golden_test

    use precision_utilities_mod, only: ik, rk
    use beta_parameterization_mod, only: cache_t, LEGENDRE_VALID
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

end program beta_param_golden_test
