program beta_param_lifecycle_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    use beta_parameterization_mod, only: &
            cache_t, cache_init_s, cache_free_s, &
            cache_max_params_f, cache_n_thetas_f, cache_is_initialized_f, &
            node_set_t, node_set_build_s, node_set_free_s, node_set_n_nodes_f, &
            SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_TOO_MANY_PARAMS, &
            BETA_PARAM_ERROR_POLE_NODE
    implicit none

    call run_cache_tests_s()
    call run_node_set_tests_s()
    call test_summary()

contains

    subroutine run_cache_tests_s()
        type(cache_t) :: cache, copy
        integer(kind = ik) :: status
        real(kind = rk) :: thetas4(4)
        integer(kind = ik) :: i

        do i = 1_ik, 4_ik
            thetas4(i) = real(i, rk) * PI_C / 5.0_rk
        end do

        ! free of a never-initialized cache is a no-op
        call cache_free_s(cache)
        call assert_true(.not. cache_is_initialized_f(cache), 'free of an uninitialized cache')
        call assert_int_eq(cache_max_params_f(cache), 0_ik, 'uninitialized max_params is 0')
        call assert_int_eq(cache_n_thetas_f(cache), 0_ik, 'uninitialized n_thetas is 0')

        call cache_init_s(cache, 8_ik, thetas4, status)
        call assert_int_eq(status, SHAPE_VALID, 'cache init ok')
        call assert_int_eq(cache_max_params_f(cache), 8_ik, 'max_params getter')
        call assert_int_eq(cache_n_thetas_f(cache), 4_ik, 'n_thetas getter')
        call assert_true(cache_is_initialized_f(cache), 'initialized getter')

        ! re-init of a live cache without a free: intent(out) releases the old one
        call cache_init_s(cache, 3_ik, thetas4(1:2), status)
        call assert_int_eq(status, SHAPE_VALID, 're-init without free ok')
        call assert_int_eq(cache_max_params_f(cache), 3_ik, 're-init replaced max_params')
        call assert_int_eq(cache_n_thetas_f(cache), 2_ik, 're-init replaced n_thetas')

        ! intrinsic assignment is a deep copy: the copy survives the original
        copy = cache
        call cache_free_s(cache)
        call assert_true(.not. cache_is_initialized_f(cache), 'freed cache')
        call assert_int_eq(cache_max_params_f(cache), 0_ik, 'freed cache reports 0')
        call assert_true(cache_is_initialized_f(copy), 'copy outlives the original')
        call assert_int_eq(cache_max_params_f(copy), 3_ik, 'copy kept max_params')
        call cache_free_s(copy)

        call cache_init_s(cache, 0_ik, thetas4, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'max_params 0 rejected with 5')
        call assert_true(.not. cache_is_initialized_f(cache), 'failed init leaves uninitialized')
        call cache_init_s(cache, 65_ik, thetas4, status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'max_params 65 rejected with 1')
        call cache_init_s(cache, 64_ik, thetas4, status)
        call assert_int_eq(status, SHAPE_VALID, 'max_params 64 accepted')
        call cache_free_s(cache)

        call cache_init_s(cache, 4_ik, thetas4(1:1), status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'single theta rejected')

        thetas4(2) = 0.0_rk
        call cache_init_s(cache, 4_ik, thetas4, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, 'exact pole rejected')
        thetas4(2) = 1.0e-9_rk   ! cos rounds to 1: rounded-cosine guard must fire
        call cache_init_s(cache, 4_ik, thetas4, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, 'rounded-cosine pole rejected')
    end subroutine run_cache_tests_s

    subroutine run_node_set_tests_s()
        type(cache_t)    :: cache
        type(node_set_t) :: ns
        integer(kind = ik) :: status, i
        real(kind = rk) :: thetas6(6), nodes3(3)

        do i = 1_ik, 6_ik
            thetas6(i) = real(i, rk) * PI_C / 7.0_rk
        end do
        nodes3 = [0.3_rk, 1.1_rk, 2.6_rk]

        call node_set_build_s(ns, cache, nodes3, status)
        call assert_int_eq(status, SHAPE_ERROR_CACHE_NOT_INITIALIZED, &
                'node set from an uninitialized cache rejected with 2')
        call assert_int_eq(node_set_n_nodes_f(ns), 0_ik, 'failed build leaves unbuilt')

        call cache_init_s(cache, 6_ik, thetas6, status)
        call node_set_build_s(ns, cache, nodes3, status)
        call assert_int_eq(status, SHAPE_VALID, 'node set builds')
        call assert_int_eq(node_set_n_nodes_f(ns), 3_ik, 'n_nodes getter')

        call node_set_build_s(ns, cache, nodes3(1:1), status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'single node rejected with 3')
        call assert_int_eq(node_set_n_nodes_f(ns), 0_ik, 'rejected rebuild leaves unbuilt')

        nodes3(2) = PI_C
        call node_set_build_s(ns, cache, nodes3, status)
        call assert_int_eq(status, BETA_PARAM_ERROR_POLE_NODE, 'pole node rejected')

        ! contract lifetime rule: node sets are freed before their cache
        nodes3(2) = 1.1_rk
        call node_set_build_s(ns, cache, nodes3, status)
        call assert_int_eq(status, SHAPE_VALID, 'node set rebuilds after a rejection')
        call node_set_free_s(ns)
        call assert_int_eq(node_set_n_nodes_f(ns), 0_ik, 'freed set reports 0')
        call cache_free_s(cache)
    end subroutine run_node_set_tests_s

end program beta_param_lifecycle_test
