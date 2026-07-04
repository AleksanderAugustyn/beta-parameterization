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
            cache_t, node_set_t, &
            compute_radius_grid_standalone_s, &
            compute_radius_grid_standalone_with_com_shift_s, &
            LEGENDRE_VALID, LEGENDRE_ERROR_INVALID_MAX_PARAMS, &
            LEGENDRE_ERROR_INVALID_BUFFER_SIZE

    implicit none

    private

    public :: beta_param_cache_create
    public :: beta_param_cache_destroy
    public :: beta_param_cache_compute_radius_grid
    public :: beta_param_cache_compute_radius_grid_with_com_shift
    public :: beta_param_node_set_create
    public :: beta_param_node_set_destroy
    public :: beta_param_cache_resolve_shape
    public :: beta_param_cache_compute_radius_and_derivative
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
        allocate(f_radii(n_grid_local), source = 0.0_rk)
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
        allocate(f_radii(n_grid_local), source = 0.0_rk)
        f_params(:) = real(params(:), rk)
        f_corrected = 0.0_rk

        call cache_ptr%compute_radius_grid_with_com_shift( &
                f_params, f_radii, f_corrected, error_code, f_message)

        radii(1:n_grid_local) = real(f_radii(:), rk_c)
        corrected_beta10 = real(f_corrected, rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
    end function beta_param_cache_compute_radius_grid_with_com_shift

    !===========================================================================
    ! NODE SETS
    !===========================================================================

    !> Allocate a node_set_t on the heap, build it, return an opaque handle.
    !! Returns c_null_ptr on failure; reason written to message_buf.
    function beta_param_node_set_create( &
            cache_handle, thetas, n_thetas, message_buf_len, message_buf) &
            result(handle) bind(c, name='beta_param_node_set_create')

        type(c_ptr),              intent(in), value :: cache_handle
        integer(kind = ik_c),     intent(in), value :: n_thetas
        real(kind = rk_c),        intent(in)        :: thetas(n_thetas)
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        type(c_ptr) :: handle

        type(cache_t),    pointer :: cache_ptr
        type(node_set_t), pointer :: node_set_ptr
        real(kind = rk), allocatable :: f_thetas(:)
        integer(kind = ik)        :: error_code
        character(len = 256)      :: f_message

        handle = c_null_ptr
        if (.not. c_associated(cache_handle)) then
            call marshal_message_to_c('node_set_create: NULL cache handle', &
                    message_buf, message_buf_len)
            return
        end if
        call c_f_pointer(cache_handle, cache_ptr)

        allocate(f_thetas(int(n_thetas, ik)))
        f_thetas(:) = real(thetas(:), rk)

        allocate(node_set_ptr)
        call cache_ptr%build_node_set(f_thetas, node_set_ptr, error_code, f_message)
        if (error_code /= LEGENDRE_VALID) then
            deallocate(node_set_ptr)
            call marshal_message_to_c(f_message, message_buf, message_buf_len)
            return
        end if

        handle = c_loc(node_set_ptr)
        call marshal_message_to_c('', message_buf, message_buf_len)
    end function beta_param_node_set_create

    !> Destroy a node set. Null-safe.
    subroutine beta_param_node_set_destroy(handle) &
            bind(c, name='beta_param_node_set_destroy')
        type(c_ptr), intent(in), value :: handle
        type(node_set_t), pointer :: node_set_ptr
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, node_set_ptr)
        call node_set_ptr%destroy()
        deallocate(node_set_ptr)
    end subroutine beta_param_node_set_destroy

    function beta_param_cache_resolve_shape( &
            handle, params, n_params, beta_con, corrected_beta10, &
            r_north, r_south, apply_com_correction, message_buf_len, message_buf) &
            result(status) bind(c, name='beta_param_cache_resolve_shape')

        type(c_ptr),              intent(in), value :: handle
        integer(kind = ik_c),     intent(in), value :: n_params
        real(kind = rk_c),        intent(in)        :: params(n_params)
        real(kind = rk_c),        intent(out)       :: beta_con(*)
        real(kind = rk_c),        intent(out)       :: corrected_beta10
        real(kind = rk_c),        intent(out)       :: r_north
        real(kind = rk_c),        intent(out)       :: r_south
        integer(kind = ik_c),     intent(in), value :: apply_com_correction
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        type(cache_t), pointer :: cache_ptr
        real(kind = rk), allocatable :: f_params(:), f_beta_con(:)
        real(kind = rk)      :: f_corrected, f_r_north, f_r_south
        integer(kind = ik)   :: error_code, max_params_local
        character(len = 256) :: f_message

        corrected_beta10 = 0.0_rk_c
        r_north          = 0.0_rk_c
        r_south          = 0.0_rk_c

        if (.not. c_associated(handle)) then
            status = int(LEGENDRE_ERROR_INVALID_MAX_PARAMS, ik_c)
            call marshal_message_to_c('null cache handle', message_buf, message_buf_len)
            return
        end if

        call c_f_pointer(handle, cache_ptr)
        max_params_local = cache_ptr%max_beta_params_get()

        allocate(f_params(int(n_params, ik)))
        allocate(f_beta_con(max_params_local), source = 0.0_rk)
        f_params(:) = real(params(:), rk)

        call cache_ptr%resolve_shape( &
                f_params, f_beta_con, f_corrected, f_r_north, f_r_south, &
                error_code, f_message, &
                apply_com_correction = (apply_com_correction /= 0_ik_c))

        beta_con(1:max_params_local) = real(f_beta_con(:), rk_c)
        corrected_beta10 = real(f_corrected, rk_c)
        r_north          = real(f_r_north, rk_c)
        r_south          = real(f_r_south, rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
    end function beta_param_cache_resolve_shape

    function beta_param_cache_compute_radius_and_derivative( &
            handle, node_set_handle, beta_con, n_beta_con, &
            radii, dr_dtheta, n_nodes, message_buf_len, message_buf) &
            result(status) bind(c, name='beta_param_cache_compute_radius_and_derivative')

        type(c_ptr),              intent(in), value :: handle
        type(c_ptr),              intent(in), value :: node_set_handle
        integer(kind = ik_c),     intent(in), value :: n_beta_con
        real(kind = rk_c),        intent(in)        :: beta_con(n_beta_con)
        integer(kind = ik_c),     intent(in), value :: n_nodes
        real(kind = rk_c),        intent(out)       :: radii(n_nodes)
        real(kind = rk_c),        intent(out)       :: dr_dtheta(n_nodes)
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        type(cache_t),    pointer :: cache_ptr
        type(node_set_t), pointer :: node_set_ptr
        real(kind = rk), allocatable :: f_beta_con(:), f_radii(:), f_dr(:)
        integer(kind = ik)   :: error_code
        character(len = 256) :: f_message

        radii(:)     = 0.0_rk_c
        dr_dtheta(:) = 0.0_rk_c

        if (.not. c_associated(handle)) then
            status = int(LEGENDRE_ERROR_INVALID_MAX_PARAMS, ik_c)
            call marshal_message_to_c('null cache handle', message_buf, message_buf_len)
            return
        end if
        if (.not. c_associated(node_set_handle)) then
            status = int(LEGENDRE_ERROR_INVALID_BUFFER_SIZE, ik_c)
            call marshal_message_to_c('null node set handle', message_buf, message_buf_len)
            return
        end if

        call c_f_pointer(handle, cache_ptr)
        call c_f_pointer(node_set_handle, node_set_ptr)

        allocate(f_beta_con(int(n_beta_con, ik)))
        allocate(f_radii(int(n_nodes, ik)), source = 0.0_rk)
        allocate(f_dr(int(n_nodes, ik)), source = 0.0_rk)
        f_beta_con(:) = real(beta_con(:), rk)

        call cache_ptr%compute_radius_and_derivative( &
                f_beta_con, node_set_ptr, f_radii, f_dr, error_code, f_message)

        radii(:)     = real(f_radii(:), rk_c)
        dr_dtheta(:) = real(f_dr(:), rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
    end function beta_param_cache_compute_radius_and_derivative

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
        allocate(f_radii(int(n_grid, ik)), source = 0.0_rk)
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
        allocate(f_radii(int(n_grid, ik)), source = 0.0_rk)
        f_params(:) = real(params(:), rk)
        f_corrected = 0.0_rk

        call compute_radius_grid_standalone_with_com_shift_s( &
                f_params, int(n_grid, ik), f_radii, f_corrected, error_code, f_message)

        radii(:) = real(f_radii(:), rk_c)
        corrected_beta10 = real(f_corrected, rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
    end function beta_param_compute_radius_grid_standalone_with_com_shift

end module beta_parameterization_c_api_mod
