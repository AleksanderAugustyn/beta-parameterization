!> Node-set API tests: derivative tables, node-set build, resolve/evaluate split.
program beta_param_node_set_test

    use precision_utilities_mod, only: ik, rk
    use beta_parameterization_workers_mod, only: precompute_legendre_table_s, &
            precompute_legendre_derivative_table_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary

    implicit none

    call test_derivative_table_analytic()
    call test_summary()

contains

    !> P_k prime against closed forms: P1=1, P2=3x, P3=(15x^2-3)/2, P4=(35x^3-15x)/2.
    subroutine test_derivative_table_analytic()
        real(kind = rk)    :: x(4), p_table(4, 5), dp_table(4, 5), expected
        integer(kind = ik) :: i

        x = [-0.9_rk, -0.3_rk, 0.2_rk, 0.7_rk]
        call precompute_legendre_table_s(x, 4_ik, p_table)
        call precompute_legendre_derivative_table_s(x, 4_ik, p_table, dp_table)

        do i = 1_ik, 4_ik
            call assert_close(dp_table(i, 1), 0.0_rk, 1.0e-15_rk, 'deriv table: dP0 == 0')
            call assert_close(dp_table(i, 2), 1.0_rk, 1.0e-15_rk, 'deriv table: dP1 == 1')
            expected = 3.0_rk * x(i)
            call assert_close(dp_table(i, 3), expected, 1.0e-14_rk, 'deriv table: dP2 == 3x')
            expected = (15.0_rk * x(i)**2 - 3.0_rk) / 2.0_rk
            call assert_close(dp_table(i, 4), expected, 1.0e-14_rk, 'deriv table: dP3')
            expected = (35.0_rk * x(i)**3 - 15.0_rk * x(i)) / 2.0_rk
            call assert_close(dp_table(i, 5), expected, 1.0e-14_rk, 'deriv table: dP4')
        end do
    end subroutine test_derivative_table_analytic

end program beta_param_node_set_test
