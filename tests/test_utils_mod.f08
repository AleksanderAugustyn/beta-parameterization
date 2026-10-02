!> Minimal assertion helpers for the beta-parameterization test suites.
module test_utils_mod

    use precision_utilities_mod, only: ik, ikl, rk

    implicit none

    private

    public :: assert_true, assert_int_eq, assert_close, test_summary
    public :: bits_equal_f, all_zero_f

    integer(kind = ik) :: n_pass = 0_ik
    integer(kind = ik) :: n_fail = 0_ik

contains

    subroutine assert_true(cond, label)
        logical,            intent(in) :: cond
        character(len = *), intent(in) :: label
        if (cond) then
            n_pass = n_pass + 1_ik
        else
            n_fail = n_fail + 1_ik
            write(*, '(A,A)') 'FAIL: ', label
        end if
    end subroutine assert_true

    subroutine assert_int_eq(got, want, label)
        integer(kind = ik), intent(in) :: got, want
        character(len = *), intent(in) :: label
        if (got == want) then
            n_pass = n_pass + 1_ik
        else
            n_fail = n_fail + 1_ik
            write(*, '(A,A,A,I0,A,I0)') 'FAIL: ', label, ' — got ', got, ', want ', want
        end if
    end subroutine assert_int_eq

    !> Mixed absolute/relative closeness: |got - want| <= tol * max(1, |want|).
    subroutine assert_close(got, want, tol, label)
        real(kind = rk),    intent(in) :: got, want, tol
        character(len = *), intent(in) :: label
        if (abs(got - want) <= tol * max(1.0_rk, abs(want))) then
            n_pass = n_pass + 1_ik
        else
            n_fail = n_fail + 1_ik
            write(*, '(A,A,A,ES23.16,A,ES23.16)') 'FAIL: ', label, ' — got ', got, ', want ', want
        end if
    end subroutine assert_close

    !> Exact equality on the IEEE754 bit patterns, so +0.0 and -0.0 differ.
    !! -Wcompare-reals rejects `==` on reals, and a scalar transfer avoids the
    !! array temporary that -Warray-temporaries flags.
    pure function bits_equal_f(a, b) result(ok)
        real(kind = rk), intent(in) :: a(:), b(:)
        logical :: ok
        integer(kind = ik) :: j
        ok = size(a, kind = ik) == size(b, kind = ik)
        if (.not. ok) return
        do j = 1_ik, size(a, kind = ik)
            if (transfer(a(j), 0_ikl) /= transfer(b(j), 0_ikl)) then
                ok = .false.
                return
            end if
        end do
    end function bits_equal_f

    !> .true. when every element is +0.0 or -0.0: the zero-fill a rejected
    !! call must leave behind.
    pure function all_zero_f(a) result(ok)
        real(kind = rk), intent(in) :: a(:)
        logical :: ok
        integer(kind = ik) :: j
        ok = .true.
        do j = 1_ik, size(a, kind = ik)
            if (abs(a(j)) > 0.0_rk) then
                ok = .false.
                return
            end if
        end do
    end function all_zero_f

    subroutine test_summary()
        write(*, '(A,I0,A,I0,A)') 'Tests: ', n_pass, ' passed, ', n_fail, ' failed.'
        if (n_fail > 0_ik) error stop 1
    end subroutine test_summary

end module test_utils_mod
