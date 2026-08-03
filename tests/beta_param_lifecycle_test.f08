program beta_param_lifecycle_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    use beta_parameterization_mod, only: &
            tables_t, tables_init_s, tables_free_s, tables_max_l_f, tables_n_thetas_f, &
            SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_INVALID_INIT, &
            BETA_PARAM_ERROR_POLE_NODE
    implicit none

    call run_tables_tests_s()
    call test_summary()

contains

    subroutine run_tables_tests_s()
        type(tables_t) :: tables
        integer(kind = ik) :: status
        real(kind = rk) :: thetas4(4)
        integer(kind = ik) :: i

        do i = 1_ik, 4_ik
            thetas4(i) = real(i, rk) * PI_C / 5.0_rk
        end do

        call tables_init_s(tables, 8_ik, thetas4, status)
        call assert_int_eq(status, SHAPE_VALID, 'tables init ok')
        call assert_int_eq(tables_max_l_f(tables), 8_ik, 'max_l getter')
        call assert_int_eq(tables_n_thetas_f(tables), 4_ik, 'n_thetas getter')

        call tables_free_s(tables)
        call assert_int_eq(tables_max_l_f(tables), 0_ik, 'freed tables report 0')

        call tables_init_s(tables, 0_ik, thetas4, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'max_l 0 rejected')
        call tables_init_s(tables, 65_ik, thetas4, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'max_l 65 rejected')

        call tables_init_s(tables, 4_ik, thetas4(1:1), status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'single theta rejected')

        thetas4(2) = 0.0_rk
        call tables_init_s(tables, 4_ik, thetas4, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, 'exact pole rejected')
        thetas4(2) = 1.0e-9_rk   ! cos rounds to 1: rounded-cosine guard must fire
        call tables_init_s(tables, 4_ik, thetas4, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, 'rounded-cosine pole rejected')
    end subroutine run_tables_tests_s

end program beta_param_lifecycle_test
