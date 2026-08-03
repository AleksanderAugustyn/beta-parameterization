program beta_param_lifecycle_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    use beta_parameterization_mod, only: &
            tables_t, tables_init_s, tables_free_s, tables_max_l_f, tables_n_thetas_f, &
            node_set_t, node_set_build_s, node_set_free_s, node_set_n_nodes_f, &
            cache_t, cache_init_s, cache_init_shared_s, cache_free_s, &
            cache_n_params_f, cache_n_thetas_f, cache_is_initialized_f, &
            SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, SHAPE_ERROR_TOO_MANY_PARAMS, &
            BETA_PARAM_ERROR_POLE_NODE
    implicit none

    call run_tables_tests_s()
    call run_node_set_tests_s()
    call run_cache_tests_s()
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

    subroutine run_cache_tests_s()
        type(tables_t), target :: tables   ! target REQUIRED for cache_init_shared_s
        type(cache_t)  :: cache
        integer(kind = ik) :: status, i
        real(kind = rk) :: thetas8(8)

        do i = 1_ik, 8_ik
            thetas8(i) = real(i, rk) * PI_C / 9.0_rk
        end do

        call cache_init_s(cache, 4_ik, thetas8, .true., .true., status)
        call assert_int_eq(status, SHAPE_VALID, 'private cache init ok')
        call assert_int_eq(cache_n_params_f(cache), 4_ik, 'n_params getter')
        call assert_int_eq(cache_n_thetas_f(cache), 8_ik, 'n_thetas getter')
        call assert_true(cache_is_initialized_f(cache), 'initialized getter')
        call cache_free_s(cache)
        call assert_true(.not. cache_is_initialized_f(cache), 'freed cache')

        call cache_init_s(cache, 9_ik, thetas8, .false., .false., status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'cap 8: 9 rejected')
        call cache_init_s(cache, 0_ik, thetas8, .false., .false., status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'n_params 0 rejected')
        call cache_init_s(cache, 4_ik, thetas8(1:1), .false., .false., status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'theta floor at cache init')

        call cache_init_shared_s(cache, tables, 4_ik, .false., .false., status)
        call assert_int_eq(status, SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
                'shared init needs initialized tables')
        call tables_init_s(tables, 3_ik, thetas8, status)
        call cache_init_shared_s(cache, tables, 4_ik, .false., .false., status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'n_params > tables max_l rejected')
        call cache_init_shared_s(cache, tables, 3_ik, .false., .false., status)
        call assert_int_eq(status, SHAPE_VALID, 'shared init ok')
        call assert_int_eq(cache_n_thetas_f(cache), 8_ik, 'shared cache sees primary set')
        call cache_free_s(cache)
        call tables_free_s(tables)
    end subroutine run_cache_tests_s

end program beta_param_lifecycle_test
