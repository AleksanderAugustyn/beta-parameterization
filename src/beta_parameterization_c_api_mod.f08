!> C-interop API for the beta parameterization library.
!!
!! Three layers of procedures:
!!   - Cache lifecycle: create/destroy, opaque handle (`type(c_ptr)`).
!!   - Cache hot path:  compute radius grid (with/without COM shift).
!!   - Standalone:      single-shape one-off (builds + destroys a cache).
!!
!! All status codes match the Fortran `LEGENDRE_*` parameters and the C
!! `BETA_PARAM_*` macros — single source of truth.
module beta_parameterization_c_api_mod

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer, c_associated, c_null_ptr
    use c_bindings_mod, only: ik_c, rk_c, c_char, c_null_char
    use precision_utilities_mod, only: ik, rk
    use beta_parameterization_mod, only: &
            cache_t, &
            compute_radius_grid_standalone_s, &
            compute_radius_grid_standalone_with_com_shift_s, &
            LEGENDRE_VALID, LEGENDRE_ERROR_INVALID_MAX_PARAMS

    implicit none

    private

    public :: beta_param_cache_create
    public :: beta_param_cache_destroy
    public :: beta_param_cache_compute_radius_grid
    public :: beta_param_cache_compute_radius_grid_with_com_shift
    public :: beta_param_compute_radius_grid_standalone
    public :: beta_param_compute_radius_grid_standalone_with_com_shift

contains

    !===========================================================================
    ! INTERNAL HELPER — Fortran string → C buffer (null-terminated, truncated to fit)
    !===========================================================================

    pure subroutine marshal_message_to_c(f_message, c_buf, c_buf_len)
        character(len = *),         intent(in)  :: f_message
        integer(kind = ik_c),       intent(in)  :: c_buf_len
        character(kind = c_char),   intent(out) :: c_buf(c_buf_len)

        integer(kind = ik) :: i, msg_len, max_copy

        if (c_buf_len < 1_ik_c) return
        msg_len  = len_trim(f_message)
        max_copy = min(int(msg_len, ik), int(c_buf_len, ik) - 1_ik)
        do i = 1_ik, max_copy
            c_buf(i) = f_message(i:i)
        end do
        c_buf(max_copy + 1_ik) = c_null_char
    end subroutine marshal_message_to_c

    !===========================================================================
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Allocate a cache_t on the heap, initialize it, return an opaque handle.
    !! Returns c_null_ptr on init failure; reason written to message_buf.
    function beta_param_cache_create( &
            max_beta_params, n_grid, message_buf_len, message_buf) &
            result(handle) bind(c, name='beta_param_cache_create')

        integer(kind = ik_c),     intent(in), value :: max_beta_params
        integer(kind = ik_c),     intent(in), value :: n_grid
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        type(c_ptr) :: handle

        type(cache_t), pointer :: cache_ptr
        integer(kind = ik)     :: error_code
        character(len = 256)   :: f_message

        allocate(cache_ptr)
        call cache_ptr%init( &
                int(max_beta_params, ik), int(n_grid, ik), &
                error_code, f_message)

        if (error_code /= LEGENDRE_VALID) then
            deallocate(cache_ptr)
            handle = c_null_ptr
            call marshal_message_to_c(f_message, message_buf, message_buf_len)
            return
        end if

        handle = c_loc(cache_ptr)
        call marshal_message_to_c('', message_buf, message_buf_len)
    end function beta_param_cache_create

    !> Destroy a cache. Null-safe.
    subroutine beta_param_cache_destroy(handle) &
            bind(c, name='beta_param_cache_destroy')
        type(c_ptr), intent(in), value :: handle
        type(cache_t), pointer :: cache_ptr
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, cache_ptr)
        call cache_ptr%destroy()
        deallocate(cache_ptr)
    end subroutine beta_param_cache_destroy

    !===========================================================================
    ! CACHE HOT PATH
    !===========================================================================

    function beta_param_cache_compute_radius_grid( &
            handle, params, n_params, radii, message_buf_len, message_buf) &
            result(status) bind(c, name='beta_param_cache_compute_radius_grid')

        type(c_ptr),              intent(in), value :: handle
        integer(kind = ik_c),     intent(in), value :: n_params
        real(kind = rk_c),        intent(in)        :: params(n_params)
        real(kind = rk_c),        intent(out)       :: radii(*)
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        type(cache_t), pointer :: cache_ptr
        real(kind = rk), allocatable :: f_params(:), f_radii(:)
        integer(kind = ik)   :: error_code, n_grid_local
        character(len = 256) :: f_message

        if (.not. c_associated(handle)) then
            status = int(LEGENDRE_ERROR_INVALID_MAX_PARAMS, ik_c)
            call marshal_message_to_c('null cache handle', message_buf, message_buf_len)
            return
        end if

        call c_f_pointer(handle, cache_ptr)
        n_grid_local = cache_ptr%n_grid_get()

        allocate(f_params(int(n_params, ik)))
        allocate(f_radii(n_grid_local))
        f_params(:) = real(params(:), rk)

        call cache_ptr%compute_radius_grid(f_params, f_radii, error_code, f_message)

        radii(1:n_grid_local) = real(f_radii(:), rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
    end function beta_param_cache_compute_radius_grid

    function beta_param_cache_compute_radius_grid_with_com_shift( &
            handle, params, n_params, radii, corrected_beta10, &
            message_buf_len, message_buf) &
            result(status) &
            bind(c, name='beta_param_cache_compute_radius_grid_with_com_shift')

        type(c_ptr),              intent(in), value :: handle
        integer(kind = ik_c),     intent(in), value :: n_params
        real(kind = rk_c),        intent(in)        :: params(n_params)
        real(kind = rk_c),        intent(out)       :: radii(*)
        real(kind = rk_c),        intent(out)       :: corrected_beta10
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        type(cache_t), pointer :: cache_ptr
        real(kind = rk), allocatable :: f_params(:), f_radii(:)
        real(kind = rk)      :: f_corrected
        integer(kind = ik)   :: error_code, n_grid_local
        character(len = 256) :: f_message

        if (.not. c_associated(handle)) then
            status = int(LEGENDRE_ERROR_INVALID_MAX_PARAMS, ik_c)
            corrected_beta10 = 0.0_rk_c
            call marshal_message_to_c('null cache handle', message_buf, message_buf_len)
            return
        end if

        call c_f_pointer(handle, cache_ptr)
        n_grid_local = cache_ptr%n_grid_get()

        allocate(f_params(int(n_params, ik)))
        allocate(f_radii(n_grid_local))
        f_params(:) = real(params(:), rk)

        call cache_ptr%compute_radius_grid_with_com_shift( &
                f_params, f_radii, f_corrected, error_code, f_message)

        radii(1:n_grid_local) = real(f_radii(:), rk_c)
        corrected_beta10 = real(f_corrected, rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
    end function beta_param_cache_compute_radius_grid_with_com_shift

    !===========================================================================
    ! STANDALONE
    !===========================================================================

    function beta_param_compute_radius_grid_standalone( &
            params, n_params, n_grid, radii, message_buf_len, message_buf) &
            result(status) bind(c, name='beta_param_compute_radius_grid_standalone')

        integer(kind = ik_c),     intent(in), value :: n_params
        real(kind = rk_c),        intent(in)        :: params(n_params)
        integer(kind = ik_c),     intent(in), value :: n_grid
        real(kind = rk_c),        intent(out)       :: radii(n_grid)
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        real(kind = rk), allocatable :: f_params(:), f_radii(:)
        integer(kind = ik)   :: error_code
        character(len = 256) :: f_message

        allocate(f_params(int(n_params, ik)))
        allocate(f_radii(int(n_grid, ik)))
        f_params(:) = real(params(:), rk)

        call compute_radius_grid_standalone_s( &
                f_params, int(n_grid, ik), f_radii, error_code, f_message)

        radii(:) = real(f_radii(:), rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
    end function beta_param_compute_radius_grid_standalone

    function beta_param_compute_radius_grid_standalone_with_com_shift( &
            params, n_params, n_grid, radii, corrected_beta10, &
            message_buf_len, message_buf) &
            result(status) &
            bind(c, name='beta_param_compute_radius_grid_standalone_with_com_shift')

        integer(kind = ik_c),     intent(in), value :: n_params
        real(kind = rk_c),        intent(in)        :: params(n_params)
        integer(kind = ik_c),     intent(in), value :: n_grid
        real(kind = rk_c),        intent(out)       :: radii(n_grid)
        real(kind = rk_c),        intent(out)       :: corrected_beta10
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        real(kind = rk), allocatable :: f_params(:), f_radii(:)
        real(kind = rk)      :: f_corrected
        integer(kind = ik)   :: error_code
        character(len = 256) :: f_message

        allocate(f_params(int(n_params, ik)))
        allocate(f_radii(int(n_grid, ik)))
        f_params(:) = real(params(:), rk)

        call compute_radius_grid_standalone_with_com_shift_s( &
                f_params, int(n_grid, ik), f_radii, f_corrected, error_code, f_message)

        radii(:) = real(f_radii(:), rk_c)
        corrected_beta10 = real(f_corrected, rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
    end function beta_param_compute_radius_grid_standalone_with_com_shift

end module beta_parameterization_c_api_mod
