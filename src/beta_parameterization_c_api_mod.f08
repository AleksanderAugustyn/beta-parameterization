!> C-interop API for the beta parameterization library.
!!
!! Provides bind(c) procedures callable from C, C++, Python (ctypes/cffi),
!! and Rust. All arrays are caller-allocated with explicit sizes.
module beta_parameterization_c_api_mod

    use c_bindings_mod, only: ik_c, rk_c, c_char, c_null_char
    use precision_utilities_mod, only: ik, rk

    implicit none

    private

    public :: beta_param_compute_radius_grid

contains

    !> Compute the Legendre radius grid from deformation parameters.
    !!
    !! C signature:
    !! ```c
    !! void beta_param_compute_radius_grid(
    !!     const double *params, int n_params, int n_grid, int n_quad,
    !!     int apply_com_correction,
    !!     double *radii, double *corrected_beta10, int *is_valid,
    !!     int message_buf_len, char *message_buf);
    !! ```
    !!
    !! @param[in]  params                Deformation parameters β₁..βₙ
    !! @param[in]  n_params              Number of deformation parameters
    !! @param[in]  n_grid                Number of θ grid points
    !! @param[in]  n_quad                Number of GL quadrature points (0 = default)
    !! @param[in]  apply_com_correction  Apply COM correction (0=false, nonzero=true)
    !! @param[out] radii                 R(θ) values (size n_grid), caller-allocated
    !! @param[out] corrected_beta10      Corrected β₁₀ value
    !! @param[out] is_valid              1 if valid, 0 if invalid
    !! @param[in]  message_buf_len       Size of message_buf including null terminator
    !! @param[out] message_buf           Error message buffer, null-terminated
    subroutine beta_param_compute_radius_grid( &
            params, n_params, n_grid, n_quad, apply_com_correction, &
            radii, corrected_beta10, is_valid, &
            message_buf_len, message_buf) &
            bind(c, name='beta_param_compute_radius_grid')

        use beta_parameterization_mod, only: compute_legendre_radius_grid_standalone_s

        implicit none

        integer(kind = ik_c), intent(in), value :: n_params
        real(kind = rk_c), intent(in) :: params(n_params)
        integer(kind = ik_c), intent(in), value :: n_grid
        integer(kind = ik_c), intent(in), value :: n_quad
        integer(kind = ik_c), intent(in), value :: apply_com_correction
        real(kind = rk_c), intent(out) :: radii(n_grid)
        real(kind = rk_c), intent(out) :: corrected_beta10
        integer(kind = ik_c), intent(out) :: is_valid
        integer(kind = ik_c), intent(in), value :: message_buf_len
        character(kind = c_char), intent(out) :: message_buf(message_buf_len)

        ! Local variables
        real(kind = rk) :: f_params(n_params)
        real(kind = rk) :: f_radii(n_grid)
        real(kind = rk) :: f_corrected_beta10
        logical :: f_is_valid, f_apply_com
        character(len = 256) :: f_message
        integer(kind = ik) :: f_n_quad
        integer :: i, msg_len

        ! Convert C inputs to Fortran types
        f_params(:) = real(params(:), rk)
        f_apply_com = (apply_com_correction /= 0_ik_c)

        ! Call the standalone routine
        if (n_quad > 0_ik_c) then
            f_n_quad = int(n_quad, ik)
            call compute_legendre_radius_grid_standalone_s( &
                    f_params, int(n_grid, ik), &
                    f_radii, f_corrected_beta10, f_is_valid, f_message, &
                    f_apply_com, f_n_quad)
        else
            call compute_legendre_radius_grid_standalone_s( &
                    f_params, int(n_grid, ik), &
                    f_radii, f_corrected_beta10, f_is_valid, f_message, &
                    f_apply_com)
        end if

        ! Convert Fortran outputs to C types
        radii(:) = real(f_radii(:), rk_c)
        corrected_beta10 = real(f_corrected_beta10, rk_c)

        if (f_is_valid) then
            is_valid = 1_ik_c
        else
            is_valid = 0_ik_c
        end if

        ! Copy message string with null termination
        msg_len = len_trim(f_message)
        if (msg_len > message_buf_len - 1) msg_len = message_buf_len - 1

        do i = 1, msg_len
            message_buf(i) = f_message(i:i)
        end do
        message_buf(msg_len + 1) = c_null_char

    end subroutine beta_param_compute_radius_grid

end module beta_parameterization_c_api_mod
