program beta_param_lifecycle_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    use beta_parameterization_mod, only: &
            tables_t, tables_init_s, tables_free_s, tables_max_l_f, tables_n_thetas_f, &
            node_set_t, node_set_build_s, node_set_free_s, node_set_n_nodes_f, &
            SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
            BETA_PARAM_ERROR_POLE_NODE
    implicit none

    call run_tables_tests_s()
    call run_node_set_tests_s()
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

    subroutine run_node_set_tests_s()
        type(tables_t)   :: tables
        type(node_set_t) :: ns
        integer(kind = ik) :: status, i
        real(kind = rk) :: thetas6(6), nodes3(3)

        do i = 1_ik, 6_ik
            thetas6(i) = real(i, rk) * PI_C / 7.0_rk
        end do
        nodes3 = [0.3_rk, 1.1_rk, 2.6_rk]

        call node_set_build_s(ns, tables, nodes3, status)
        call assert_int_eq(status, SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
                'node set from uninitialized tables rejected')
        call assert_int_eq(node_set_n_nodes_f(ns), 0_ik, 'failed build leaves unbuilt')

        call tables_init_s(tables, 6_ik, thetas6, status)
        call node_set_build_s(ns, tables, nodes3, status)
        call assert_int_eq(status, SHAPE_VALID, 'node set builds')
        call assert_int_eq(node_set_n_nodes_f(ns), 3_ik, 'n_nodes getter')

        nodes3(2) = PI_C
        call node_set_build_s(ns, tables, nodes3, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, 'pole node rejected')

        call node_set_free_s(ns)
        call tables_free_s(tables)
        call assert_int_eq(node_set_n_nodes_f(ns), 0_ik, 'freed set reports 0')
    end subroutine run_node_set_tests_s

end program beta_param_lifecycle_test
