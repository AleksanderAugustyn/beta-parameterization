!> C-interop layer for the beta parameterization library (`beta_parameterization.h`).
!!
!! Two opaque handles, each a `c_loc` of a heap-allocated derived type:
!!   - `beta_param_cache_t`    -> `cache_t`    (immutable after create)
!!   - `beta_param_node_set_t` -> `node_set_t` (immutable after create)
!!
!! Both may be shared across threads: every compute takes its handles read-only.
!!
!! Every `_create` returns a null handle on failure and reports the cause
!! through a trailing nullable `int* status` (an absent `optional` dummy when C
!! passes NULL). Every other entry returns the status code directly. No entry
!! point stops, allocates for the caller, or formats a message: diagnostics are
!! the static strings behind `beta_param_status_message`.
!!
!! This layer only marshals. Buffer-length and parameter-count contracts are
!! checked once, in `beta_parameterization_mod`, and their codes are passed
!! through untouched. The checks that live here are the ones the Fortran layer
!! cannot see: the NULL-handle guard, the negative-count clamp, and the
!! allocation guard described next.
!!
!! ## Marshalling buffers are HEAP, not automatic
!!
!! An automatic array is allocated on procedure entry, BEFORE the stated size
!! reaches the Fortran tier that would reject it. Under Release
!! (`-fstack-arrays`) a large wrong size argument — the wrong variable passed
!! as `n_radii` — was therefore a stack overflow instead of a status code.
!! Every caller-sized marshalling buffer here is `allocatable` with an explicit
!! `allocate(..., stat = ...)`, so an outsized request either reaches the
!! Fortran tier (and is rejected there) or fails allocation recoverably.
!! Allocation failure maps to the code of the size argument implicated:
!! `n_radii` / `n_nodes` -> 104; `n_thetas` -> 3; `n_params` -> 4 in a cached
!! call and 1 in a one-shot call (a count that cannot be allocated exceeds 64).
!! A one-shot call has no output size argument — its buffers are sized by
!! `n_thetas` — so every buffer failure there is a 3. Negative counts are
!! clamped to zero and judged by the Fortran tier.
!!
!! ## A stated size must be the ACTUAL buffer extent
!!
!! The size arguments are a contract, not a bound this layer can verify: an
!! output dummy is declared with the caller's stated extent, so this layer
!! zero-fills exactly that many elements. A stated size LARGER than the
!! caller's real buffer is undefined behaviour that no check here can catch.
module beta_parameterization_c_api_mod

    use, intrinsic :: iso_c_binding, only: &
            c_int, c_double, c_char, c_ptr, c_loc, c_f_pointer, c_associated, &
            c_null_ptr, c_null_char
    use precision_utilities_mod, only: ik, rk
    use beta_parameterization_mod, only: &
            cache_t, node_set_t, &
            cache_init_s, cache_free_s, &
            node_set_build_s, node_set_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_radius_grid_unchecked_s, cache_resolve_shape_s, &
            cache_node_radius_and_derivative_s, &
            compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, &
            STATUS_MESSAGE_LEN, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_SOUTH_POLE, &
            BETA_PARAM_ERROR_INTERIOR_NEGATIVE, &
            BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            BETA_PARAM_ERROR_POLE_NODE, BETA_PARAM_ERROR_NODE_SET_MISMATCH

    implicit none

    private

    public :: beta_param_status_message
    public :: beta_param_cache_create, beta_param_cache_destroy
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
    integer(kind = ik), parameter :: N_MSG   = 14_ik

    !> Codes in column order; column N_MSG is the unknown-code fallback and has
    !! no entry here.
    integer(kind = ik), parameter :: MSG_CODE(N_MSG - 1_ik) = [ &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            BETA_PARAM_ERROR_NORTH_POLE, BETA_PARAM_ERROR_SOUTH_POLE, &
            BETA_PARAM_ERROR_INTERIOR_NEGATIVE, &
            BETA_PARAM_ERROR_COM_NOT_CONVERGED, &
            BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, &
            BETA_PARAM_ERROR_POLE_NODE, BETA_PARAM_ERROR_NODE_SET_MISMATCH]

    character(kind = c_char, len = MSG_LEN), parameter :: MSG_TEXT(N_MSG) = &
            [character(kind = c_char, len = MSG_LEN) :: &
                    'valid' // c_null_char, &
                    'too many parameters' // c_null_char, &
                    'cache not initialized' // c_null_char, &
                    'theta grid below minimum size (2)' // c_null_char, &
                    'params length outside 1..max_params' // c_null_char, &
                    'invalid init arguments' // c_null_char, &
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

    !> A C count as a Fortran extent: negative counts become zero, so the
    !! Fortran tier sees an empty array and rejects it with its own code.
    pure function extent_f(n) result(extent)
        integer(c_int), intent(in) :: n
        integer(kind = ik) :: extent
        extent = max(int(n, ik), 0_ik)
    end function extent_f

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
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Build the read-only cache. Null handle on failure.
    function beta_param_cache_create(max_params, thetas, n_thetas, status) &
            result(handle) bind(c, name = 'beta_param_cache_create')
        integer(c_int), value, intent(in) :: max_params, n_thetas
        real(c_double), intent(in) :: thetas(n_thetas)
        integer(c_int), intent(out), optional :: status   ! NULL-able from C
        type(c_ptr) :: handle

        type(cache_t), pointer :: p
        real(kind = rk), allocatable :: thetas_f(:)
        integer(kind = ik) :: st
        integer :: alloc_stat

        handle = c_null_ptr
        st = SHAPE_ERROR_INVALID_GRID          ! thetas cannot be marshalled
        allocate(thetas_f(extent_f(n_thetas)), stat = alloc_stat)
        if (alloc_stat == 0) then
            thetas_f(:) = real(thetas, rk)
            st = SHAPE_ERROR_INVALID_INIT      ! the handle object cannot be allocated
            allocate(p, stat = alloc_stat)
            if (alloc_stat == 0) then
                call cache_init_s(p, int(max_params, ik), thetas_f, st)
                if (st == SHAPE_VALID) then
                    handle = c_loc(p)
                else
                    deallocate(p)
                end if
            end if
        end if
        if (present(status)) status = int(st, c_int)
    end function beta_param_cache_create

    !> Release a cache. NULL-safe.
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

    !> Build an extra evaluation set from a cache. Null on failure.
    function beta_param_node_set_create(cache, thetas, n_thetas, status) &
            result(handle) bind(c, name = 'beta_param_node_set_create')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_thetas
        real(c_double), intent(in) :: thetas(n_thetas)
        integer(c_int), intent(out), optional :: status
        type(c_ptr) :: handle

        type(node_set_t), pointer :: p
        type(cache_t),    pointer :: cp
        real(kind = rk), allocatable :: thetas_f(:)
        integer(kind = ik) :: st
        integer :: alloc_stat

        handle = c_null_ptr
        st = SHAPE_ERROR_CACHE_NOT_INITIALIZED
        if (c_associated(cache)) then
            call c_f_pointer(cache, cp)
            st = SHAPE_ERROR_INVALID_GRID
            allocate(thetas_f(extent_f(n_thetas)), stat = alloc_stat)
            if (alloc_stat == 0) then
                thetas_f(:) = real(thetas, rk)
                st = SHAPE_ERROR_INVALID_INIT
                allocate(p, stat = alloc_stat)
                if (alloc_stat == 0) then
                    call node_set_build_s(p, cp, thetas_f, st)
                    if (st == SHAPE_VALID) then
                        handle = c_loc(p)
                    else
                        deallocate(p)
                    end if
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
    function beta_param_cache_radius_grid(cache, params, n_params, &
            conserve_volume, apply_com, radii, n_radii) &
            result(status) bind(c, name = 'beta_param_cache_radius_grid')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_radii
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_radii)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        real(kind = rk), allocatable :: params_f(:), radii_f(:)
        integer(kind = ik) :: st
        integer :: alloc_stat

        radii = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        allocate(params_f(extent_f(n_params)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_WRONG_PARAM_COUNT, c_int)
            return
        end if
        allocate(radii_f(extent_f(n_radii)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, c_int)
            return
        end if
        params_f(:) = real(params, rk)
        call cache_radius_grid_s(p, params_f, conserve_volume /= 0_c_int, &
                apply_com /= 0_c_int, radii_f, st)
        radii = real(radii_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_radius_grid

    !> R and dR/dtheta at the cache's primary thetas.
    function beta_param_cache_radius_and_derivative(cache, params, n_params, &
            conserve_volume, apply_com, radii, dr_dthetas, n_radii) &
            result(status) bind(c, name = 'beta_param_cache_radius_and_derivative')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_radii
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_radii), dr_dthetas(n_radii)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        real(kind = rk), allocatable :: params_f(:), radii_f(:), dr_f(:)
        integer(kind = ik) :: st
        integer :: alloc_stat

        radii      = 0.0_c_double
        dr_dthetas = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        allocate(params_f(extent_f(n_params)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_WRONG_PARAM_COUNT, c_int)
            return
        end if
        allocate(radii_f(extent_f(n_radii)), dr_f(extent_f(n_radii)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, c_int)
            return
        end if
        params_f(:) = real(params, rk)
        call cache_radius_and_derivative_s(p, params_f, conserve_volume /= 0_c_int, &
                apply_com /= 0_c_int, radii_f, dr_f, st)
        radii      = real(radii_f, c_double)
        dr_dthetas = real(dr_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_radius_and_derivative

    !> R at the primary thetas with no validation gates and no volume scaling.
    function beta_param_cache_radius_grid_unchecked(cache, params, n_params, &
            apply_com, radii, n_radii) &
            result(status) bind(c, name = 'beta_param_cache_radius_grid_unchecked')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_radii
        integer(c_int), value, intent(in) :: apply_com
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_radii)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        real(kind = rk), allocatable :: params_f(:), radii_f(:)
        integer(kind = ik) :: st
        integer :: alloc_stat

        radii = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        allocate(params_f(extent_f(n_params)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_WRONG_PARAM_COUNT, c_int)
            return
        end if
        allocate(radii_f(extent_f(n_radii)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, c_int)
            return
        end if
        params_f(:) = real(params, rk)
        call cache_radius_grid_unchecked_s(p, params_f, apply_com /= 0_c_int, radii_f, st)
        radii = real(radii_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_radius_grid_unchecked

    !> COM-corrected beta10, analytic polar radii and the applied volume factor.
    function beta_param_cache_resolve_shape(cache, params, n_params, &
            conserve_volume, apply_com, &
            corrected_beta10, r_north, r_south, volume_factor) &
            result(status) bind(c, name = 'beta_param_cache_resolve_shape')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: corrected_beta10, r_north, r_south, volume_factor
        integer(c_int) :: status

        type(cache_t), pointer :: p
        real(kind = rk), allocatable :: params_f(:)
        real(kind = rk) :: beta10_f, r_north_f, r_south_f, volume_f
        integer(kind = ik) :: st
        integer :: alloc_stat

        corrected_beta10 = 0.0_c_double
        r_north          = 0.0_c_double
        r_south          = 0.0_c_double
        volume_factor    = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        allocate(params_f(extent_f(n_params)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_WRONG_PARAM_COUNT, c_int)
            return
        end if
        params_f(:) = real(params, rk)
        call cache_resolve_shape_s(p, params_f, conserve_volume /= 0_c_int, &
                apply_com /= 0_c_int, beta10_f, r_north_f, r_south_f, volume_f, st)
        corrected_beta10 = real(beta10_f, c_double)
        r_north          = real(r_north_f, c_double)
        r_south          = real(r_south_f, c_double)
        volume_factor    = real(volume_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_resolve_shape

    !> R and dR/dtheta at a node set's thetas.
    function beta_param_cache_node_radius_and_derivative(cache, params, n_params, &
            node_set, conserve_volume, apply_com, radii, dr_dthetas, n_nodes) &
            result(status) bind(c, name = 'beta_param_cache_node_radius_and_derivative')
        type(c_ptr), value, intent(in) :: cache, node_set
        integer(c_int), value, intent(in) :: n_params, n_nodes
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_nodes), dr_dthetas(n_nodes)
        integer(c_int) :: status

        type(cache_t),    pointer :: p
        type(node_set_t), pointer :: ns
        real(kind = rk), allocatable :: params_f(:), radii_f(:), dr_f(:)
        integer(kind = ik) :: st
        integer :: alloc_stat

        radii      = 0.0_c_double
        dr_dthetas = 0.0_c_double
        if (.not. c_associated(cache) .or. .not. c_associated(node_set)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        call c_f_pointer(node_set, ns)
        allocate(params_f(extent_f(n_params)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_WRONG_PARAM_COUNT, c_int)
            return
        end if
        allocate(radii_f(extent_f(n_nodes)), dr_f(extent_f(n_nodes)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(BETA_PARAM_ERROR_INVALID_BUFFER_SIZE, c_int)
            return
        end if
        params_f(:) = real(params, rk)
        call cache_node_radius_and_derivative_s(p, params_f, ns, &
                conserve_volume /= 0_c_int, apply_com /= 0_c_int, radii_f, dr_f, st)
        radii      = real(radii_f, c_double)
        dr_dthetas = real(dr_f, c_double)
        status = int(st, c_int)
    end function beta_param_cache_node_radius_and_derivative

    !===========================================================================
    ! ONE-SHOT COMPUTES (tier 1)
    !===========================================================================

    !> One-shot R at caller thetas; builds and discards its own cache.
    function beta_param_radius_grid_standalone(params, n_params, thetas, n_thetas, &
            conserve_volume, apply_com, radii) &
            result(status) bind(c, name = 'beta_param_radius_grid_standalone')
        integer(c_int), value, intent(in) :: n_params, n_thetas
        real(c_double), intent(in) :: params(n_params), thetas(n_thetas)
        integer(c_int), value, intent(in) :: conserve_volume, apply_com
        real(c_double), intent(out) :: radii(n_thetas)
        integer(c_int) :: status

        real(kind = rk), allocatable :: params_f(:), thetas_f(:), radii_f(:)
        integer(kind = ik) :: st
        integer :: alloc_stat

        radii = 0.0_c_double
        allocate(params_f(extent_f(n_params)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_TOO_MANY_PARAMS, c_int)
            return
        end if
        allocate(thetas_f(extent_f(n_thetas)), radii_f(extent_f(n_thetas)), &
                stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_INVALID_GRID, c_int)
            return
        end if
        params_f(:) = real(params, rk)
        thetas_f(:) = real(thetas, rk)
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

        real(kind = rk), allocatable :: params_f(:), thetas_f(:), radii_f(:), dr_f(:)
        integer(kind = ik) :: st
        integer :: alloc_stat

        radii      = 0.0_c_double
        dr_dthetas = 0.0_c_double
        allocate(params_f(extent_f(n_params)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_TOO_MANY_PARAMS, c_int)
            return
        end if
        allocate(thetas_f(extent_f(n_thetas)), radii_f(extent_f(n_thetas)), &
                dr_f(extent_f(n_thetas)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = int(SHAPE_ERROR_INVALID_GRID, c_int)
            return
        end if
        params_f(:) = real(params, rk)
        thetas_f(:) = real(thetas, rk)
        call compute_radius_and_derivative_standalone_s(params_f, thetas_f, &
                conserve_volume /= 0_c_int, apply_com /= 0_c_int, radii_f, dr_f, st)
        radii      = real(radii_f, c_double)
        dr_dthetas = real(dr_f, c_double)
        status = int(st, c_int)
    end function beta_param_radius_and_derivative_standalone

end module beta_parameterization_c_api_mod
