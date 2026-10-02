program beta_param_status_test
    use precision_utilities_mod, only: ik
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    use beta_parameterization_mod, only: &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_MAX_PARAMS, MAX_BETA_PARAMS_LIMIT, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_SOUTH_POLE, &
            BETA_PARAM_ERROR_INTERIOR_NEGATIVE, BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, BETA_PARAM_ERROR_POLE_NODE, &
            BETA_PARAM_ERROR_NODE_SET_MISMATCH, STATUS_MESSAGE_LEN, status_message_f
    implicit none

    character(len = STATUS_MESSAGE_LEN) :: msg

    call assert_int_eq(SHAPE_VALID, 0_ik, 'SHAPE_VALID is 0')
    call assert_int_eq(SHAPE_ERROR_TOO_MANY_PARAMS, 1_ik, 'shared 1 re-exported')
    call assert_int_eq(SHAPE_ERROR_CACHE_NOT_INITIALIZED, 2_ik, 'shared 2 re-exported')
    call assert_int_eq(SHAPE_ERROR_INVALID_GRID, 3_ik, 'shared 3 re-exported')
    call assert_int_eq(SHAPE_ERROR_WRONG_PARAM_COUNT, 4_ik, 'shared 4 re-exported')
    call assert_int_eq(SHAPE_ERROR_INVALID_INIT, 5_ik, 'shared 5 re-exported')
    call assert_int_eq(SHAPE_MAX_PARAMS, 64_ik, 'SHAPE_MAX_PARAMS re-exported')
    call assert_int_eq(MAX_BETA_PARAMS_LIMIT, 64_ik, 'library N_max is 64')
    call assert_int_eq(BETA_PARAM_ERROR_NORTH_POLE, 100_ik, 'north pole 100')
    call assert_int_eq(BETA_PARAM_ERROR_SOUTH_POLE, 101_ik, 'south pole 101')
    call assert_int_eq(BETA_PARAM_ERROR_INTERIOR_NEGATIVE, 102_ik, 'interior 102')
    call assert_int_eq(BETA_PARAM_ERROR_COM_NOT_CONVERGED, 103_ik, 'com 103')
    call assert_int_eq(BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, 104_ik, 'buffer 104')
    call assert_int_eq(BETA_PARAM_ERROR_POLE_NODE, 105_ik, 'pole node 105')
    call assert_int_eq(BETA_PARAM_ERROR_NODE_SET_MISMATCH, 106_ik, 'mismatch 106')

    msg = status_message_f(SHAPE_VALID)
    call assert_true(index(msg, 'valid') > 0, 'valid message mentions valid')
    msg = status_message_f(BETA_PARAM_ERROR_NORTH_POLE)
    call assert_true(index(msg, 'north') > 0, 'north message mentions north')
    msg = status_message_f(SHAPE_ERROR_WRONG_PARAM_COUNT)
    call assert_true(index(msg, 'max_params') > 0, 'code 4 message names max_params')
    msg = status_message_f(6_ik)
    call assert_true(index(msg, 'unknown') > 0, 'retired code 6 maps to unknown')
    msg = status_message_f(9999_ik)
    call assert_true(index(msg, 'unknown') > 0, 'unknown code maps to unknown')
    call test_summary()
end program beta_param_status_test
