!> C-interop layer for the beta parameterization library (`beta_parameterization.h`).
!!
!! Three opaque handles, each a `c_loc` of a heap-allocated derived type:
!!   - `beta_param_tables_t`   -> `tables_t`   (shared, immutable after create)
!!   - `beta_param_node_set_t` -> `node_set_t` (shared, immutable after create)
!!   - `beta_param_cache_t`    -> `cache_t`    (THREAD-CONFINED, mutated by
!!                                              every compute call)
!!
!! Every `_create` returns a null handle on failure and reports the cause
!! through a trailing nullable `int* status` (an absent `optional` dummy when C
!! passes NULL). Every other entry returns the status code directly. No entry
!! point stops, allocates for the caller, or formats a message: diagnostics are
!! the static strings behind `beta_param_status_message`.
!!
!! This layer only marshals. Buffer-length and parameter-count contracts are
!! checked once, in `beta_parameterization_mod`, and their codes are passed
!! through untouched; the one check that lives here is the NULL-handle guard,
!! which the Fortran layer cannot see.
module beta_parameterization_c_api_mod

    use, intrinsic :: iso_c_binding, only: &
            c_int, c_double, c_char, c_ptr, c_loc, c_f_pointer, c_associated, &
            c_null_ptr, c_null_char
    use precision_utilities_mod, only: ik, rk
    use beta_parameterization_mod, only: &
            tables_t, cache_t, node_set_t, &
            tables_init_s, tables_free_s, &
            node_set_build_s, node_set_free_s, &
            cache_init_s, cache_init_shared_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_radius_grid_unchecked_s, cache_resolve_shape_s, &
            cache_node_radius_and_derivative_s, &
            compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, &
            STATUS_MESSAGE_LEN, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_SOUTH_POLE, &
            BETA_PARAM_ERROR_INTERIOR_NEGATIVE, &
            BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            BETA_PARAM_ERROR_POLE_NODE, BETA_PARAM_ERROR_NODE_SET_MISMATCH

    implicit none

    private

    public :: beta_param_status_message
    public :: beta_param_tables_create, beta_param_tables_destroy
    public :: beta_param_cache_create, beta_param_cache_create_shared
    public :: beta_param_cache_destroy
    public :: beta_param_node_set_create, beta_param_node_set_destroy
    public :: beta_param_cache_radius_grid
    public :: beta_param_cache_radius_and_derivative
    public :: beta_param_cache_radius_grid_unchecked
    public :: beta_param_cache_resolve_shape
    public :: beta_param_cache_node_radius_and_derivative
    public :: beta_param_radius_grid_standalone
    public :: beta_param_radius_and_derivative_standalone

    !---------------------------------------------------------------------------
    ! Static status strings
    !---------------------------------------------------------------------------
    !> One null-terminated string per status code, laid out as fixed columns of a
    !! module array that is initialized at compile time and never written again.
    !! `beta_param_status_message` hands out `c_loc` of a column, so the pointer
    !! is valid forever and every thread gets the same read-only bytes — a shared
    !! mutable buffer would make the message racy.
    !!
    !! The texts are the `status_message_f` texts: keep the two in sync (a C
    !! caller and a Fortran caller must not read different words for one code).
    integer(kind = ik), parameter :: MSG_LEN = STATUS_MESSAGE_LEN + 1_ik
    integer(kind = ik), parameter :: N_MSG   = 15_ik

    !> Codes in column order; column N_MSG is the unknown-code fallback and has
    !! no entry here.
    integer(kind = ik), parameter :: MSG_CODE(N_MSG - 1_ik) = [ &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_SOUTH_POLE, &
            BETA_PARAM_ERROR_INTERIOR_NEGATIVE, &
            BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            BETA_PARAM_ERROR_POLE_NODE, BETA_PARAM_ERROR_NODE_SET_MISMATCH]

    character(kind = c_char, len = MSG_LEN), parameter :: MSG_TEXT(N_MSG) = &
            [character(kind = c_char, len = MSG_LEN) :: &
                    'valid' // c_null_char, &
                    'too many parameters for this tier' // c_null_char, &
                    'cache not initialized' // c_null_char, &
                    'theta grid below minimum size (2)' // c_null_char, &
                    'params length differs from n_params' // c_null_char, &
                    'invalid init arguments' // c_null_char, &
                    'tables not initialized' // c_null_char, &
                    'north pole radius not positive' // c_null_char, &
                    'south pole radius not positive' // c_null_char, &
                    'interior radius not positive' // c_null_char, &
                    'COM correction did not converge' // c_null_char, &
                    'output buffer size mismatch' // c_null_char, &
                    'theta at or beyond a pole' // c_null_char, &
                    'node set unbuilt or max_l too small' // c_null_char, &
                    'unknown status code' // c_null_char]

    character(kind = c_char), save, target :: MSG_CHARS(MSG_LEN, N_MSG) = &
            reshape(transfer(MSG_TEXT, c_null_char, MSG_LEN * N_MSG), &
                    [MSG_LEN, N_MSG])

contains

    !===========================================================================
    ! DIAGNOSTICS
    !===========================================================================

    !> Pointer to the static, null-terminated description of `status`.
    !! Unknown codes get the fallback column. Never returns NULL.
    function beta_param_status_message(status) result(msg) &
            bind(c, name = 'beta_param_status_message')
        integer(c_int), value, intent(in) :: status
        type(c_ptr) :: msg
        integer(kind = ik) :: i, idx

        idx = N_MSG
        do i = 1_ik, N_MSG - 1_ik
            if (status == int(MSG_CODE(i), c_int)) then
                idx = i
                exit
            end if
        end do
        msg = c_loc(MSG_CHARS(1_ik, idx))
    end function beta_param_status_message

    !===========================================================================
    ! TABLES LIFECYCLE
    !===========================================================================

    !> Build the shared immutable level. Null handle on failure.
    function beta_param_tables_create(max_l, thetas, n_thetas, status) &
            result(handle) bind(c, name = 'beta_param_tables_create')
        integer(c_int), value, intent(in) :: max_l, n_thetas
        real(c_double), intent(in) :: thetas(n_thetas)
        integer(c_int), intent(out), optional :: status   ! NULL-able from C
        type(c_ptr) :: handle

        type(tables_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: thetas_f(n_thetas)

        handle = c_null_ptr
        st = SHAPE_ERROR_INVALID_INIT
        if (n_thetas >= 0_c_int) then
            thetas_f = real(thetas, rk)
            allocate(p)
            call tables_init_s(p, int(max_l, ik), thetas_f, st)
            if (st == SHAPE_VALID) then
                handle = c_loc(p)
            else
                deallocate(p)
            end if
        end if
        if (present(status)) status = int(st, c_int)
    end function beta_param_tables_create

    !> Release a tables handle. NULL-safe.
    subroutine beta_param_tables_destroy(handle) &
            bind(c, name = 'beta_param_tables_destroy')
        type(c_ptr), value, intent(in) :: handle
        type(tables_t), pointer :: p
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, p)
        call tables_free_s(p)
        deallocate(p)
    end subroutine beta_param_tables_destroy

    !===========================================================================
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Create a cache owning private tables (max_l = n_params). Null on failure.
    function beta_param_cache_create(n_params, thetas, n_thetas, &
            conserve_volume, apply_com, status) &
            result(handle) bind(c, name = 'beta_param_cache_create')
        integer(c_int), value, intent(in) :: n_params, n_thetas
        real(c_double), intent(in) :: thetas(n_thetas)
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        integer(c_int), intent(out), optional :: status
        type(c_ptr) :: handle

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: thetas_f(n_thetas)

        handle = c_null_ptr
        st = SHAPE_ERROR_INVALID_INIT
        if (n_thetas >= 0_c_int) then
            thetas_f = real(thetas, rk)
            allocate(p)
            call cache_init_s(p, int(n_params, ik), thetas_f, &
                    conserve_volume /= 0_c_int, apply_com /= 0_c_int, st)
            if (st == SHAPE_VALID) then
                handle = c_loc(p)
            else
                call cache_free_s(p)   ! a failed init may still own tables
                deallocate(p)
            end if
        end if
        if (present(status)) status = int(st, c_int)
    end function beta_param_cache_create

    !> Create a cache bound to caller-owned shared tables, which must outlive it.
    !! Null on failure.
    function beta_param_cache_create_shared(tables, n_params, &
            conserve_volume, apply_com, status) &
            result(handle) bind(c, name = 'beta_param_cache_create_shared')
        type(c_ptr), value, intent(in) :: tables
        integer(c_int), value, intent(in) :: n_params
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        integer(c_int), intent(out), optional :: status
        type(c_ptr) :: handle

        type(cache_t),  pointer :: p
        type(tables_t), pointer :: tp
        integer(kind = ik) :: st

        handle = c_null_ptr
        st = SHAPE_ERROR_TABLES_NOT_INITIALIZED
        if (c_associated(tables)) then
            call c_f_pointer(tables, tp)
            allocate(p)
            call cache_init_shared_s(p, tp, int(n_params, ik), &
                    conserve_volume /= 0_c_int, apply_com /= 0_c_int, st)
            if (st == SHAPE_VALID) then
                handle = c_loc(p)
            else
                call cache_free_s(p)
                deallocate(p)
            end if
        end if
        if (present(status)) status = int(st, c_int)
    end function beta_param_cache_create_shared

    !> Release a cache. NULL-safe. Shared tables are left to their owner.
    subroutine beta_param_cache_destroy(handle) &
            bind(c, name = 'beta_param_cache_destroy')
        type(c_ptr), value, intent(in) :: handle
        type(cache_t), pointer :: p
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, p)
        call cache_free_s(p)
        deallocate(p)
    end subroutine beta_param_cache_destroy

    !===========================================================================
    ! NODE-SET LIFECYCLE
    !===========================================================================

    !> Build an extra evaluation set against shared tables. Null on failure.
    function beta_param_node_set_create(tables, thetas, n_thetas, status) &
            result(handle) bind(c, name = 'beta_param_node_set_create')
        type(c_ptr), value, intent(in) :: tables
        integer(c_int), value, intent(in) :: n_thetas
        real(c_double), intent(in) :: thetas(n_thetas)
        integer(c_int), intent(out), optional :: status
        type(c_ptr) :: handle

        type(node_set_t), pointer :: p
        type(tables_t),   pointer :: tp
        integer(kind = ik) :: st
        real(kind = rk) :: thetas_f(n_thetas)

        handle = c_null_ptr
        st = SHAPE_ERROR_TABLES_NOT_INITIALIZED
        if (c_associated(tables)) then
            call c_f_pointer(tables, tp)
            st = SHAPE_ERROR_INVALID_GRID
            if (n_thetas >= 0_c_int) then
                thetas_f = real(thetas, rk)
                allocate(p)
                call node_set_build_s(p, tp, thetas_f, st)
                if (st == SHAPE_VALID) then
                    handle = c_loc(p)
                else
                    call node_set_free_s(p)
                    deallocate(p)
                end if
            end if
        end if
        if (present(status)) status = int(st, c_int)
    end function beta_param_node_set_create

    !> Release a node set. NULL-safe.
    subroutine beta_param_node_set_destroy(handle) &
            bind(c, name = 'beta_param_node_set_destroy')
        type(c_ptr), value, intent(in) :: handle
        type(node_set_t), pointer :: p
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, p)
        call node_set_free_s(p)
        deallocate(p)
    end subroutine beta_param_node_set_destroy

    !===========================================================================
    ! CACHED COMPUTES
    !===========================================================================

    !> R at the cache's primary thetas.
    function beta_param_cache_radius_grid(cache, params, n_params, radii, n_radii) &
            result(status) bind(c, name = 'beta_param_cache_radius_grid')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_radii
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_radii)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(n_params), radii_f(n_radii)

        radii = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_radius_grid_s(p, params_f, radii_f, st)
        radii = real(radii_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_radius_grid

    !> R and dR/dtheta at the cache's primary thetas.
    function beta_param_cache_radius_and_derivative(cache, params, n_params, &
            radii, dr_dthetas, n_radii) &
            result(status) bind(c, name = 'beta_param_cache_radius_and_derivative')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_radii
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_radii), dr_dthetas(n_radii)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(n_params), radii_f(n_radii), dr_f(n_radii)

        radii      = 0.0_c_double
        dr_dthetas = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_radius_and_derivative_s(p, params_f, radii_f, dr_f, st)
        radii      = real(radii_f, c_double)
        dr_dthetas = real(dr_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_radius_and_derivative

    !> R at the primary thetas with no validation gates and no volume scaling.
    function beta_param_cache_radius_grid_unchecked(cache, params, n_params, &
            radii, n_radii) &
            result(status) bind(c, name = 'beta_param_cache_radius_grid_unchecked')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_radii
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_radii)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(n_params), radii_f(n_radii)

        radii = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_radius_grid_unchecked_s(p, params_f, radii_f, st)
        radii = real(radii_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_radius_grid_unchecked

    !> COM-corrected beta10, analytic polar radii and the applied volume factor.
    function beta_param_cache_resolve_shape(cache, params, n_params, &
            corrected_beta10, r_north, r_south, volume_factor) &
            result(status) bind(c, name = 'beta_param_cache_resolve_shape')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: corrected_beta10, r_north, r_south, volume_factor
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(n_params)
        real(kind = rk) :: beta10_f, r_north_f, r_south_f, volume_f

        corrected_beta10 = 0.0_c_double
        r_north          = 0.0_c_double
        r_south          = 0.0_c_double
        volume_factor    = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_resolve_shape_s(p, params_f, beta10_f, r_north_f, r_south_f, &
                volume_f, st)
        corrected_beta10 = real(beta10_f, c_double)
        r_north          = real(r_north_f, c_double)
        r_south          = real(r_south_f, c_double)
        volume_factor    = real(volume_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_resolve_shape

    !> R and dR/dtheta at a node set's thetas (uncached evaluation).
    function beta_param_cache_node_radius_and_derivative(cache, node_set, &
            params, n_params, radii, dr_dthetas, n_nodes) &
            result(status) bind(c, name = 'beta_param_cache_node_radius_and_derivative')
        type(c_ptr), value, intent(in) :: cache, node_set
        integer(c_int), value, intent(in) :: n_params, n_nodes
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_nodes), dr_dthetas(n_nodes)
        integer(c_int) :: status

        type(cache_t),    pointer :: p
        type(node_set_t), pointer :: ns
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(n_params), radii_f(n_nodes), dr_f(n_nodes)

        radii      = 0.0_c_double
        dr_dthetas = 0.0_c_double
        if (.not. c_associated(cache) .or. .not. c_associated(node_set)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        call c_f_pointer(node_set, ns)
        params_f = real(params, rk)
        call cache_node_radius_and_derivative_s(p, ns, params_f, radii_f, dr_f, st)
        radii      = real(radii_f, c_double)
        dr_dthetas = real(dr_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_node_radius_and_derivative

    !===========================================================================
    ! STANDALONE COMPUTES (tier 1)
    !===========================================================================

    !> One-shot R at caller thetas; builds and discards its own tables.
    function beta_param_radius_grid_standalone(params, n_params, thetas, n_thetas, &
            conserve_volume, apply_com, radii) &
            result(status) bind(c, name = 'beta_param_radius_grid_standalone')
        integer(c_int), value, intent(in) :: n_params, n_thetas
        real(c_double), intent(in) :: params(n_params), thetas(n_thetas)
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        real(c_double), intent(out) :: radii(n_thetas)
        integer(c_int) :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(n_params), thetas_f(n_thetas), radii_f(n_thetas)

        radii = 0.0_c_double
        if (n_params < 0_c_int .or. n_thetas < 0_c_int) then
            status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        thetas_f = real(thetas, rk)
        call compute_radius_grid_standalone_s(params_f, thetas_f, &
                conserve_volume /= 0_c_int, apply_com /= 0_c_int, radii_f, st)
        radii = real(radii_f, c_double)
        status = int(st, c_int)
    end function beta_param_radius_grid_standalone

    !> One-shot R and dR/dtheta at caller thetas.
    function beta_param_radius_and_derivative_standalone(params, n_params, &
            thetas, n_thetas, conserve_volume, apply_com, radii, dr_dthetas) &
            result(status) bind(c, name = 'beta_param_radius_and_derivative_standalone')
        integer(c_int), value, intent(in) :: n_params, n_thetas
        real(c_double), intent(in) :: params(n_params), thetas(n_thetas)
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        real(c_double), intent(out) :: radii(n_thetas), dr_dthetas(n_thetas)
        integer(c_int) :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(n_params), thetas_f(n_thetas)
        real(kind = rk) :: radii_f(n_thetas), dr_f(n_thetas)

        radii      = 0.0_c_double
        dr_dthetas = 0.0_c_double
        if (n_params < 0_c_int .or. n_thetas < 0_c_int) then
            status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        thetas_f = real(thetas, rk)
        call compute_radius_and_derivative_standalone_s(params_f, thetas_f, &
                conserve_volume /= 0_c_int, apply_com /= 0_c_int, radii_f, dr_f, st)
        radii      = real(radii_f, c_double)
        dr_dthetas = real(dr_f, c_double)
        status = int(st, c_int)
    end function beta_param_radius_and_derivative_standalone

end module beta_parameterization_c_api_mod
