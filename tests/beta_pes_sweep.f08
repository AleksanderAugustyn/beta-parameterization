!> PES-range sweep: resolve every shape on a regular 8D beta grid.
!!
!! Two runs over the same grid, both with volume conservation on:
!!   Run A - apply_com = .false., all 8 dimensions swept.
!!   Run B - apply_com = .true.,  dimension 1 (beta10) pinned to 0, since the
!!           COM correction determines beta10 itself.
!!
!! Every accepted shape is checked against an INDEPENDENT volume integral: the
!! sweep owns its Gauss-Legendre weights and evaluates SUM(w_i R(x_i)^3), which
!! must equal 2 (the unit-sphere value) to within 1e-12 relative. Nothing in the
!! check reads library internals.
!!
!! Not a ctest suite: the full grid is 1.55M shapes. Built as an executable and
!! driven by the `test_beta_pes_sweep` custom target. Invalid shapes are a
!! normal outcome and are only reported; the program fails (exit 1) only on a
!! volume mismatch or a status outside the permitted validation set.
!!
!! Usage: beta_pes_sweep [--step <value>]   (default 0.05)
program beta_pes_sweep

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_utilities_mod, only: compute_gauss_legendre_quadrature_s
    use beta_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_resolve_shape_s, cache_radius_grid_s, SHAPE_VALID, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_SOUTH_POLE, &
            BETA_PARAM_ERROR_INTERIOR_NEGATIVE, BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            BETA_PARAM_ERROR_POLE_NODE

    implicit none

    integer(kind = ik), parameter :: N_DIMS = 8_ik, N_GL = 128_ik
    real(kind = rk), parameter :: LO(N_DIMS) = &
            [-0.20_rk, 0.00_rk, 0.00_rk, -0.15_rk, 0.00_rk, -0.10_rk, 0.00_rk, -0.10_rk]
    real(kind = rk), parameter :: HI(N_DIMS) = &
            [ 0.20_rk, 0.40_rk, 0.25_rk,  0.20_rk, 0.15_rk,  0.10_rk, 0.15_rk,  0.10_rk]
    real(kind = rk) :: step
    ! per-dim counts: n(d) = nint((HI(d)-LO(d))/step) + 1  -> at step 0.05:
    ! 9 9 6 8 4 5 4 5 -> 1,555,200 points (8D) / 172,800 (7D, dim 1 pinned to 0)

    ! Volume of the unit sphere in the quadrature's own units: INT_{-1}^{1} dx.
    real(kind = rk), parameter :: TARGET_VOLUME = 2.0_rk
    real(kind = rk), parameter :: VOLUME_REL_TOL = 1.0e-12_rk

    step = 0.05_rk
    call parse_args_s(step)          ! accepts: --step <value>
    call run_sweep_s(apply_com = .false., step = step)   ! Run A: 8D
    call run_sweep_s(apply_com = .true.,  step = step)   ! Run B: 7D

contains

    !> Read the optional `--step <value>` argument; anything else is fatal.
    subroutine parse_args_s(step_out)
        real(kind = rk), intent(inout) :: step_out

        character(len = 64) :: arg, value_arg
        integer(kind = ik)  :: n_args, i, io_status

        n_args = int(command_argument_count(), ik)
        i = 1_ik
        do while (i <= n_args)
            call get_command_argument(int(i), arg)
            select case (trim(arg))
            case ('--step')
                if (i == n_args) call usage_stop_s('--step requires a value')
                call get_command_argument(int(i + 1_ik), value_arg)
                read(value_arg, *, iostat = io_status) step_out
                if (io_status /= 0) call usage_stop_s('--step value is not a number')
                if (step_out <= 0.0_rk) call usage_stop_s('--step value must be positive')
                i = i + 2_ik
            case default
                call usage_stop_s('unknown argument: ' // trim(arg))
            end select
        end do
    end subroutine parse_args_s

    !> Print the usage line with a reason, then stop hard.
    subroutine usage_stop_s(reason)
        character(len = *), intent(in) :: reason

        write(*, '(A,A)') 'error: ', reason
        write(*, '(A)')   'usage: beta_pes_sweep [--step <value>]   (default 0.05)'
        error stop 1
    end subroutine usage_stop_s

    !> One full grid sweep in one regime.
    subroutine run_sweep_s(apply_com, step)
        logical,         intent(in) :: apply_com
        real(kind = rk), intent(in) :: step

        type(cache_t)      :: cache
        real(kind = rk)    :: gl_x(N_GL), gl_w(N_GL), thetas(N_GL), radii(N_GL)
        real(kind = rk)    :: params(N_DIMS), first_bad_params(N_DIMS)
        real(kind = rk)    :: corrected_beta10, r_north, r_south, volume_factor
        real(kind = rk)    :: volume, first_bad_volume, wall_seconds, rate_hz
        integer(kind = ik) :: idx(N_DIMS), n(N_DIMS), d, status
        integer(kind = ikl) :: n_total, n_valid, n_volume_fail, n_usage_error
        integer(kind = ikl) :: n_north, n_south, n_interior, n_com, n_pole_node
        integer(kind = ikl) :: tick_start, tick_end, tick_rate
        logical             :: have_first_bad

        ! Quadrature: nodes in x = cos(theta), weights kept for the volume check.
        call compute_gauss_legendre_quadrature_s(N_GL, gl_x, gl_w)
        thetas(:) = acos(gl_x(:))

        call cache_init_s(cache, N_DIMS, thetas, .true., apply_com, status)
        if (status /= SHAPE_VALID) then
            write(*, '(A,I0)') 'error: cache_init_s failed with status ', status
            error stop 1
        end if

        n(:) = nint((HI(:) - LO(:)) / step, kind = ik) + 1_ik
        if (apply_com) n(1) = 1_ik      ! beta10 is set by the COM correction

        n_total = 1_ikl
        do d = 1_ik, N_DIMS
            n_total = n_total * int(n(d), ikl)
        end do

        n_valid = 0_ikl; n_volume_fail = 0_ikl; n_usage_error = 0_ikl
        n_north = 0_ikl; n_south = 0_ikl; n_interior = 0_ikl
        n_com = 0_ikl;   n_pole_node = 0_ikl
        have_first_bad = .false.
        first_bad_params(:) = 0.0_rk
        first_bad_volume = 0.0_rk

        idx(:) = 1_ik
        call system_clock(count = tick_start, count_rate = tick_rate)

        do
            do d = 1_ik, N_DIMS
                params(d) = LO(d) + real(idx(d) - 1_ik, rk) * step
            end do
            if (apply_com) params(1) = 0.0_rk

            call cache_resolve_shape_s(cache, params, corrected_beta10, r_north, &
                    r_south, volume_factor, status)

            select case (status)
            case (SHAPE_VALID)
                n_valid = n_valid + 1_ikl

                call cache_radius_grid_s(cache, params, radii, status)
                if (status /= SHAPE_VALID) then
                    ! A resolved shape has already passed every validation stage,
                    ! so only a misuse of the API can fail the grid call here.
                    n_usage_error = n_usage_error + 1_ikl
                else
                    volume = sum(gl_w(:) * radii(:)**3)
                    if (abs(volume - TARGET_VOLUME) > VOLUME_REL_TOL * TARGET_VOLUME) then
                        n_volume_fail = n_volume_fail + 1_ikl
                        if (.not. have_first_bad) then
                            have_first_bad = .true.
                            first_bad_params(:) = params(:)
                            first_bad_volume = volume
                        end if
                    end if
                end if
            case (BETA_PARAM_ERROR_NORTH_POLE)
                n_north = n_north + 1_ikl
            case (BETA_PARAM_ERROR_SOUTH_POLE)
                n_south = n_south + 1_ikl
            case (BETA_PARAM_ERROR_INTERIOR_NEGATIVE)
                n_interior = n_interior + 1_ikl
            case (BETA_PARAM_ERROR_COM_NOT_CONVERGED)
                n_com = n_com + 1_ikl
            case (BETA_PARAM_ERROR_POLE_NODE)
                ! Permitted shape-validation status; the sweep's node set cannot
                ! trigger it today, but it must not read as a sweep failure.
                n_pole_node = n_pole_node + 1_ikl
            case default
                n_usage_error = n_usage_error + 1_ikl
                if (n_usage_error == 1_ikl) then
                    write(*, '(A,I0,A)') 'error: unexpected status ', status, &
                            ' at params:'
                    write(*, '(A,8(1X,F8.4))') '   ', params(:)
                end if
            end select

            ! Odometer step: last dimension fastest.
            d = N_DIMS
            do
                idx(d) = idx(d) + 1_ik
                if (idx(d) <= n(d)) exit
                idx(d) = 1_ik
                d = d - 1_ik
                if (d < 1_ik) exit
            end do
            if (d < 1_ik) exit
        end do

        call system_clock(count = tick_end)
        call cache_free_s(cache)

        wall_seconds = real(tick_end - tick_start, rk) / real(tick_rate, rk)
        if (wall_seconds > 0.0_rk) then
            rate_hz = real(n_total, rk) / wall_seconds
        else
            rate_hz = 0.0_rk
        end if

        write(*, '(A)') repeat('-', 66)
        write(*, '(A,L1,A,I0,A,F8.5)') 'PES sweep  apply_com = ', apply_com, &
                '  dims = ', merge(N_DIMS - 1_ik, N_DIMS, apply_com), &
                '  step = ', step
        write(*, '(A)') repeat('-', 66)
        write(*, '(A,8(1X,I0))')  '  per-dim counts       :', n(:)
        write(*, '(A,I0)')        '  total points         : ', n_total
        write(*, '(A,I0,A,F7.3,A)') '  valid shapes         : ', n_valid, &
                '  (', 100.0_rk * real(n_valid, rk) / real(n_total, rk), ' %)'
        write(*, '(A,I0)')        '  status 100 north pole : ', n_north
        write(*, '(A,I0)')        '  status 101 south pole : ', n_south
        write(*, '(A,I0)')        '  status 102 interior   : ', n_interior
        write(*, '(A,I0)')        '  status 103 com        : ', n_com
        write(*, '(A,I0)')        '  status 105 pole node  : ', n_pole_node
        write(*, '(A,I0)')        '  unexpected statuses   : ', n_usage_error
        write(*, '(A,I0)')        '  volume-check failures : ', n_volume_fail
        write(*, '(A,F12.3)')     '  wall seconds          : ', wall_seconds
        write(*, '(A,F14.1)')     '  shapes / second       : ', rate_hz

        if (have_first_bad) then
            write(*, '(A,ES23.15)') '  first bad volume      : ', first_bad_volume
            write(*, '(A,8(1X,F8.4))') '  first bad params      :', first_bad_params(:)
        end if

        if (n_volume_fail > 0_ikl .or. n_usage_error > 0_ikl) then
            write(*, '(A)') 'SWEEP FAILED'
            error stop 1
        end if
    end subroutine run_sweep_s

end program beta_pes_sweep
