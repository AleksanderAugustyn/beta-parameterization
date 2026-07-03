!> Node-set API tests: derivative tables, node-set build, resolve/evaluate split.
program beta_param_node_set_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use beta_parameterization_workers_mod, only: precompute_legendre_table_s, &
            precompute_legendre_derivative_table_s, eval_radius_derivative_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary

    implicit none

    call test_derivative_table_analytic()
    call test_radius_derivative_eval()
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

    !> beta_con(2)-only shape: dR/dtheta = -sin(theta) * b * 3cos(theta) exactly.
    subroutine test_radius_derivative_eval()
        real(kind = rk)    :: thetas(3), x(3), sin_thetas(3)
        real(kind = rk)    :: p_table(3, 3), dp_table(3, 3), beta_con(2), dr(3), expected
        integer(kind = ik) :: i

        thetas = [PI_C / 6.0_rk, PI_C / 3.0_rk, 2.0_rk * PI_C / 3.0_rk]
        x = cos(thetas)
        sin_thetas = sin(thetas)
        call precompute_legendre_table_s(x, 2_ik, p_table)
        call precompute_legendre_derivative_table_s(x, 2_ik, p_table, dp_table)

        beta_con = [0.0_rk, 0.35_rk]
        call eval_radius_derivative_s(beta_con, dp_table, sin_thetas, dr)

        do i = 1_ik, 3_ik
            expected = -sin_thetas(i) * beta_con(2) * 3.0_rk * x(i)
            call assert_close(dr(i), expected, 1.0e-14_rk, 'radius derivative: b2 closed form')
        end do
    end subroutine test_radius_derivative_eval

end program beta_param_node_set_test
