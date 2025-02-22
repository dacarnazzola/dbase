module array_utils
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, c_bool
implicit none
private

    interface swap
        module procedure :: swap_i32
        module procedure :: swap_i64
        module procedure :: swap_sp
        module procedure :: swap_dp
        module procedure :: swap_c_bool
    end interface swap

    interface flip
        module procedure :: flip_i32
        module procedure :: flip_i64
        module procedure :: flip_sp
        module procedure :: flip_dp
        module procedure :: flip_c_bool
    end interface flip

    public :: swap, flip

contains

    pure elemental subroutine swap_i32(a, b)
        integer(kind=i32), intent(inout) :: a, b
        integer(kind=i32) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_i32

    pure elemental subroutine swap_i64(a, b)
        integer(kind=i64), intent(inout) :: a, b
        integer(kind=i64) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_i64

    pure elemental subroutine swap_sp(a, b)
        real(kind=sp), intent(inout) :: a, b
        real(kind=sp) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_sp

    pure elemental subroutine swap_dp(a, b)
        real(kind=dp), intent(inout) :: a, b
        real(kind=dp) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_dp

    pure elemental subroutine swap_c_bool(a, b)
        logical(kind=c_bool), intent(inout) :: a, b
        logical(kind=c_bool) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap_c_bool

    pure subroutine flip_i32(x)
        integer(kind=i32), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_i32

    pure subroutine flip_i64(x)
        integer(kind=i64), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_i64

    pure subroutine flip_sp(x)
        real(kind=sp), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_sp

    pure subroutine flip_dp(x)
        real(kind=dp), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_dp

    pure subroutine flip_c_bool(x)
        logical(kind=c_bool), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        n = size(x, kind=i64)
        do i=1_i64,n/2_i64
            call swap(x(i), x(n-i+1_i64))
        end do
    end subroutine flip_c_bool

end module array_utils
