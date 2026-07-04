!> Node-set API tests: derivative tables, node-set build, resolve/evaluate split.
program beta_param_node_set_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use beta_parameterization_workers_mod, only: precompute_legendre_table_s, &
            precompute_legendre_derivative_table_s, eval_radius_derivative_s
    use beta_parameterization_mod, only: cache_t, node_set_t, LEGENDRE_VALID, &
            LEGENDRE_ERROR_POLE_NODE, LEGENDRE_ERROR_INVALID_BUFFER_SIZE, &
            LEGENDRE_ERROR_INVALID_MAX_PARAMS, LEGENDRE_ERROR_EMPTY_PARAMS, &
            LEGENDRE_ERROR_INTERIOR_NEGATIVE, LEGENDRE_ERROR_NO_UNIFORM_GRID
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, test_summary

    implicit none

    call test_derivative_table_analytic()
    call test_radius_derivative_eval()
    call test_build_node_set()
    call test_resolve_shape()
    call test_resolve_shape_no_com()
    call test_evaluate_uniform_parity()
    call test_evaluate_fd_golden()
    call test_evaluate_errors()
    call test_node_set_only_cache()
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

    subroutine test_build_node_set()
        type(cache_t)        :: cache, cold_cache
        type(node_set_t)     :: node_set
        real(kind = rk)      :: thetas(8), pole_thetas(2), empty(0)
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call cache%init(8_ik, 181_ik, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'build: cache init')

        do i = 1_ik, 8_ik
            thetas(i) = real(i, rk) * PI_C / 9.0_rk
        end do
        call cache%build_node_set(thetas, node_set, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'build: valid mid thetas')
        call assert_true(node_set%is_built_get(), 'build: is_built set')
        call assert_int_eq(node_set%n_nodes_get(), 8_ik, 'build: n_nodes')

        pole_thetas = [0.0_rk, PI_C / 2.0_rk]
        call cache%build_node_set(pole_thetas, node_set, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_POLE_NODE, 'build: pole node rejected')

        call cache%build_node_set(empty, node_set, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_INVALID_BUFFER_SIZE, 'build: empty thetas rejected')

        call cold_cache%build_node_set(thetas, node_set, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_INVALID_MAX_PARAMS, 'build: uninitialized cache rejected')
    end subroutine test_build_node_set

    !> resolve_shape must reproduce the with_com_shift path: same corrected_beta10,
    !! and r_north/r_south equal to the uniform grid endpoints (theta = 0 and pi
    !! are on that grid; the library does not volume-scale).
    subroutine test_resolve_shape()
        type(cache_t)        :: cache
        real(kind = rk)      :: params(4), beta_con(8), radii(181)
        real(kind = rk)      :: corrected_beta10, corrected_ref, r_north, r_south
        real(kind = rk)      :: bad_beta_con(3), empty(0)
        integer(kind = ik)   :: code
        character(len = 256) :: message

        call cache%init(8_ik, 181_ik, code, message)
        params = [0.1_rk, 0.2_rk, 0.05_rk, 0.1_rk]  ! odd terms -> COM iteration active

        call cache%compute_radius_grid_with_com_shift(params, radii, corrected_ref, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'resolve: reference path valid')

        call cache%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'resolve: valid')
        call assert_close(corrected_beta10, corrected_ref, 1.0e-15_rk, 'resolve: corrected_beta10 parity')
        call assert_close(r_north, radii(1), 1.0e-15_rk, 'resolve: r_north == R(0)')
        call assert_close(r_south, radii(181), 1.0e-15_rk, 'resolve: r_south == R(pi)')

        call cache%resolve_shape(empty, beta_con, corrected_beta10, r_north, r_south, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_EMPTY_PARAMS, 'resolve: empty params')

        call cache%resolve_shape(params, bad_beta_con, corrected_beta10, r_north, r_south, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_INVALID_BUFFER_SIZE, 'resolve: beta_con size mismatch')
    end subroutine test_resolve_shape

    !> apply_com_correction = .false. must reproduce compute_radius_grid (no
    !! COM shift): same tables, same dot products, beta10 reported as the
    !! input value instead of the COM-corrected one.
    subroutine test_resolve_shape_no_com()
        type(cache_t)        :: cache
        type(node_set_t)     :: node_set
        real(kind = rk)      :: params(4), beta_con(8), radii_ref(181)
        real(kind = rk)      :: thetas(179), radii(179), dr(179)
        real(kind = rk)      :: corrected_beta10, r_north, r_south
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call cache%init(8_ik, 181_ik, code, message)
        params = [0.1_rk, 0.215_rk, 0.05_rk, 0.095_rk]  ! odd terms present but COM skipped

        call cache%compute_radius_grid(params, radii_ref, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'no-com: reference grid valid')

        do i = 1_ik, 179_ik
            thetas(i) = real(i, rk) * PI_C / 180.0_rk
        end do
        call cache%build_node_set(thetas, node_set, code, message)

        call cache%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, &
                code, message, apply_com_correction = .false.)
        call assert_int_eq(code, LEGENDRE_VALID, 'no-com: resolve valid')
        call assert_close(corrected_beta10, params(1), 0.0_rk, 'no-com: beta10 passthrough')

        call cache%compute_radius_and_derivative(beta_con, node_set, radii, dr, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'no-com: evaluate valid')

        do i = 1_ik, 179_ik
            call assert_close(radii(i), radii_ref(i + 1), 1.0e-15_rk, 'no-com: R matches compute_radius_grid')
        end do
        call assert_close(r_north, radii_ref(1), 1.0e-15_rk, 'no-com: r_north == R(0)')
        call assert_close(r_south, radii_ref(181), 1.0e-15_rk, 'no-com: r_south == R(pi)')
    end subroutine test_resolve_shape_no_com

    !> Node-set evaluation at the uniform grid's own thetas must reproduce
    !! compute_radius_grid_with_com_shift exactly (same beta_con, same tables,
    !! same summation worker).
    subroutine test_evaluate_uniform_parity()
        type(cache_t)        :: cache
        type(node_set_t)     :: node_set
        real(kind = rk)      :: params(4), beta_con(8), radii_ref(181)
        real(kind = rk)      :: thetas(179), radii(179), dr(179)
        real(kind = rk)      :: corrected_beta10, corrected_ref, r_north, r_south
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call cache%init(8_ik, 181_ik, code, message)
        params = [0.1_rk, 0.2_rk, 0.05_rk, 0.1_rk]

        call cache%compute_radius_grid_with_com_shift(params, radii_ref, corrected_ref, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'parity: reference valid')

        ! Interior uniform-grid thetas (skip the two poles)
        do i = 1_ik, 179_ik
            thetas(i) = real(i, rk) * PI_C / 180.0_rk
        end do
        call cache%build_node_set(thetas, node_set, code, message)
        call cache%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, code, message)
        call cache%compute_radius_and_derivative(beta_con, node_set, radii, dr, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'parity: evaluate valid')

        do i = 1_ik, 179_ik
            call assert_close(radii(i), radii_ref(i + 1), 1.0e-15_rk, 'parity: R matches uniform grid')
        end do
    end subroutine test_evaluate_uniform_parity

    !> dR/dtheta against 5-point central FD of node-set R values, tol 1e-9
    !! (spec gate), plus the beta2-only closed form at machine precision.
    subroutine test_evaluate_fd_golden()
        type(cache_t)        :: cache
        type(node_set_t)     :: node_set, stencil_set
        real(kind = rk), parameter :: h = 1.0e-3_rk
        real(kind = rk)      :: params(5), params_b2(2), beta_con(8)
        real(kind = rk)      :: thetas(5), stencil(20), r5(5), dr5(5), rs(20), drs(20)
        real(kind = rk)      :: corrected_beta10, r_north, r_south, fd, expected
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call cache%init(8_ik, 181_ik, code, message)

        thetas = [0.3_rk, 0.9_rk, PI_C / 2.0_rk, 2.2_rk, 2.9_rk]
        do i = 1_ik, 5_ik
            stencil(4 * i - 3) = thetas(i) - 2.0_rk * h
            stencil(4 * i - 2) = thetas(i) - h
            stencil(4 * i - 1) = thetas(i) + h
            stencil(4 * i)     = thetas(i) + 2.0_rk * h
        end do
        call cache%build_node_set(thetas, node_set, code, message)
        call cache%build_node_set(stencil, stencil_set, code, message)

        params = [0.15_rk, 0.25_rk, 0.1_rk, 0.05_rk, 0.02_rk]
        call cache%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, code, message)
        call cache%compute_radius_and_derivative(beta_con, node_set, r5, dr5, code, message)
        call cache%compute_radius_and_derivative(beta_con, stencil_set, rs, drs, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'fd: evaluate valid')

        do i = 1_ik, 5_ik
            fd = (rs(4 * i - 3) - 8.0_rk * rs(4 * i - 2) &
                    + 8.0_rk * rs(4 * i - 1) - rs(4 * i)) / (12.0_rk * h)
            call assert_close(dr5(i), fd, 1.0e-9_rk, 'fd: analytic dR/dtheta vs 5-point FD')
        end do

        ! beta2-only: dR/dtheta = -sin * beta_con(2) * 3cos, exactly
        params_b2 = [0.0_rk, 0.3_rk]
        call cache%resolve_shape(params_b2, beta_con, corrected_beta10, r_north, r_south, code, message)
        call cache%compute_radius_and_derivative(beta_con, node_set, r5, dr5, code, message)
        do i = 1_ik, 5_ik
            expected = -sin(thetas(i)) * beta_con(2) * 3.0_rk * cos(thetas(i))
            call assert_close(dr5(i), expected, 1.0e-14_rk, 'fd: b2 closed form end-to-end')
        end do
    end subroutine test_evaluate_fd_golden

    subroutine test_evaluate_errors()
        type(cache_t)        :: cache
        type(node_set_t)     :: node_set, unbuilt
        real(kind = rk)      :: params(2), beta_con(8), bad_beta_con(3)
        real(kind = rk)      :: thetas(3), radii(3), dr(3), bad_radii(5)
        real(kind = rk)      :: corrected_beta10, r_north, r_south
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call cache%init(8_ik, 181_ik, code, message)
        thetas = [0.5_rk, PI_C / 2.0_rk, 2.5_rk]
        call cache%build_node_set(thetas, node_set, code, message)

        ! Interior-negative: beta2 large enough that R(pi/2) = 1 - beta_con(2)/2 < 0,
        ! while both poles stay positive (even multipole, no COM shift).
        params = [0.0_rk, 3.5_rk]
        call cache%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'errors: extreme b2 resolves (poles positive)')
        call cache%compute_radius_and_derivative(beta_con, node_set, radii, dr, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_INTERIOR_NEGATIVE, 'errors: interior negative')
        do i = 1_ik, 3_ik
            call assert_close(radii(i), 0.0_rk, 1.0e-15_rk, 'errors: radii zero-filled')
            call assert_close(dr(i), 0.0_rk, 1.0e-15_rk, 'errors: dr zero-filled')
        end do

        params = [0.0_rk, 0.3_rk]
        call cache%resolve_shape(params, beta_con, corrected_beta10, r_north, r_south, code, message)
        call cache%compute_radius_and_derivative(bad_beta_con, node_set, radii, dr, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_INVALID_BUFFER_SIZE, 'errors: beta_con size mismatch')
        call cache%compute_radius_and_derivative(beta_con, node_set, bad_radii, dr, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_INVALID_BUFFER_SIZE, 'errors: radii size mismatch')
        call cache%compute_radius_and_derivative(beta_con, unbuilt, radii, dr, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_INVALID_BUFFER_SIZE, 'errors: unbuilt node set')
    end subroutine test_evaluate_errors

    !> A cache built without n_grid serves the node-set API identically to a
    !> uniform-grid cache; the legacy uniform entry points refuse cleanly.
    subroutine test_node_set_only_cache()
        type(cache_t)        :: full, lean
        type(node_set_t)     :: ns_full, ns_lean
        real(kind = rk)      :: thetas(5), params(4)
        real(kind = rk)      :: bc_full(8), bc_lean(8)
        real(kind = rk)      :: r_full(5), r_lean(5), dr_full(5), dr_lean(5)
        real(kind = rk)      :: radii_buf(181)
        real(kind = rk)      :: corr_full, corr_lean, rn_f, rs_f, rn_l, rs_l
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        call lean%init(8_ik, error_code = code, message = message)
        call assert_int_eq(code, LEGENDRE_VALID, 'lean init: VALID')
        call assert_int_eq(lean%n_grid_get(), 0_ik, 'lean init: n_grid == 0')

        call full%init(8_ik, 181_ik, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'full init: VALID')

        thetas = [0.3_rk, 0.9_rk, 1.5_rk, 2.1_rk, 2.7_rk]
        params = [0.0_rk, 0.85_rk, 0.35_rk, 0.18_rk]

        call lean%build_node_set(thetas, ns_lean, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'lean node set: built')
        call full%build_node_set(thetas, ns_full, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'full node set: built')

        call lean%resolve_shape(params, bc_lean, corr_lean, rn_l, rs_l, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'lean resolve: VALID')
        call full%resolve_shape(params, bc_full, corr_full, rn_f, rs_f, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'full resolve: VALID')
        call assert_close(corr_lean, corr_full, 0.0_rk, 'resolve: corrected_beta10 identical')
        call assert_close(rn_l, rn_f, 0.0_rk, 'resolve: r_north identical')
        call assert_close(rs_l, rs_f, 0.0_rk, 'resolve: r_south identical')

        call lean%compute_radius_and_derivative(bc_lean, ns_lean, r_lean, dr_lean, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'lean evaluate: VALID')
        call full%compute_radius_and_derivative(bc_full, ns_full, r_full, dr_full, code, message)
        call assert_int_eq(code, LEGENDRE_VALID, 'full evaluate: VALID')
        do i = 1_ik, 5_ik
            call assert_close(r_lean(i), r_full(i), 0.0_rk, 'evaluate: radii identical')
            call assert_close(dr_lean(i), dr_full(i), 0.0_rk, 'evaluate: dr identical')
        end do

        call lean%compute_radius_grid(params, radii_buf, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_NO_UNIFORM_GRID, &
                'lean radius_grid: NO_UNIFORM_GRID')
        call assert_true(len_trim(message) > 0, 'lean radius_grid: message non-empty')
        call lean%compute_radius_grid_with_com_shift(params, radii_buf, corr_lean, code, message)
        call assert_int_eq(code, LEGENDRE_ERROR_NO_UNIFORM_GRID, &
                'lean com radius_grid: NO_UNIFORM_GRID')
    end subroutine test_node_set_only_cache

end program beta_param_node_set_test
