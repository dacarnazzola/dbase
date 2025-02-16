module sorting
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, c_bool
use, non_intrinsic :: constants, only: i32_vec_len
implicit none
private
    
    interface is_sorted
        module procedure :: is_sorted_i32
        module procedure :: is_sorted_i64
        module procedure :: is_sorted_sp
        module procedure :: is_sorted_dp
    end interface is_sorted

    interface sort
        module procedure :: qsort_i32
        module procedure :: qsort2_i32
    end interface sort

    public :: is_sorted, sort

contains

    pure function is_sorted_i32(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        logical(kind=c_bool) :: val
        integer(kind=i64) :: n, i
        val = logical(.true., kind=c_bool)
        n = size(x, kind=i64)
        do i=2_i64,n
            if (x(i) < x(i-1_i64)) then
                val = logical(.false., kind=c_bool)
                return
            end if
        end do
    end function is_sorted_i32

    pure function is_sorted_i64(x) result(val)
        integer(kind=i64), intent(in) :: x(:)
        logical(kind=c_bool) :: val
        integer(kind=i64) :: n, i
        val = logical(.true., kind=c_bool)
        n = size(x, kind=i64)
        do i=2_i64,n
            if (x(i) < x(i-1_i64)) then
                val = logical(.false., kind=c_bool)
                return
            end if
        end do
    end function is_sorted_i64

    pure function is_sorted_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        logical(kind=c_bool) :: val
        integer(kind=i64) :: n, i
        val = logical(.true., kind=c_bool)
        n = size(x, kind=i64)
        do i=2_i64,n
            if (x(i) < x(i-1_i64)) then
                val = logical(.false., kind=c_bool)
                return
            end if
        end do
    end function is_sorted_sp

    pure function is_sorted_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        logical(kind=c_bool) :: val
        integer(kind=i64) :: n, i
        val = logical(.true., kind=c_bool)
        n = size(x, kind=i64)
        do i=2_i64,n
            if (x(i) < x(i-1_i64)) then
                val = logical(.false., kind=c_bool)
                return
            end if
        end do
    end function is_sorted_dp

    pure recursive subroutine qsort_i32(x)
        integer(kind=i32), intent(inout) :: x(:)
        integer(kind=i64) :: n, i, ii
        integer(kind=i32) :: u, v, w
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            call isort_i32(x)
        else
            u = x(1_i64)
            v = x(n/2_i64)
            w = x(n)
            if ((v > u) .neqv. (v > w)) then !! v is median
                x(1_i64) = v
                x(n/2_i64) = u
            else if ((w > u) .neqv. (w > v)) then !! w is median
                x(1_i64) = w
                x(n) = u
            end if !! u is median
            ii = 1_i64
            do i=2_i64,n
                if (x(i) < x(1_i64)) then
                    ii = ii + 1_i64
                    u = x(i)
                    x(i) = x(ii)
                    x(ii) = u
                end if
            end do
            u = x(ii)
            x(ii) = x(1_i64)
            x(1_i64) = u
            call qsort_i32(x(1_i64:ii-1_i64))
            call qsort_i32(x(ii+1_i64:n))
        end if
    end subroutine qsort_i32

    pure recursive subroutine qsort2_i32(x, xi)
        integer(kind=i32), intent(inout) :: x(:), xi(size(x, kind=i64))
        integer(kind=i64) :: n, i, ii
        integer(kind=i32) :: u, v, w, uxi, vxi, wxi
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            call isort2_i32(x, xi)
        else
            u = x(1_i64)
            v = x(n/2_i64)
            w = x(n)
            uxi = xi(1_i64)
            vxi = xi(n/2_i64)
            wxi = xi(n)
            if ((v > u) .neqv. (v > w)) then !! v is median
                x(1_i64) = v
                x(n/2_i64) = u
                xi(1_i64) = vxi
                xi(n/2_i64) = uxi
            else if ((w > u) .neqv. (w > v)) then !! w is median
                x(1_i64) = w
                x(n) = u
                xi(1_i64) = wxi
                xi(n) = uxi
            end if !! u is median
            ii = 1_i64
            do i=2_i64,n
                if (x(i) < x(1_i64)) then
                    ii = ii + 1_i64
                    u = x(i)
                    x(i) = x(ii)
                    x(ii) = u
                    u = xi(i)
                    xi(i) = xi(ii)
                    xi(ii) = u
                end if
            end do
            u = x(ii)
            x(ii) = x(1_i64)
            x(1_i64) = u
            u = xi(ii)
            xi(ii) = xi(1_i64)
            xi(1_i64) = u
            call qsort2_i32(x(1_i64:ii-1_i64), xi(1_i64:ii-1_i64))
            call qsort2_i32(x(ii+1_i64:n), xi(ii+1_i64:n))
        end if
    end subroutine qsort2_i32

    pure subroutine isort_i32(x)
        integer(kind=i32), intent(inout) :: x(:)
        integer(kind=i64) :: n, i, j
        integer(kind=i32) :: temp
        n = size(x, kind=i64)
        do i=2_i64,n
            if (x(i) < x(i-1_i64)) then
                temp = x(i)
                loop_j: do j=i,2_i64,-1_i64
                    if (temp < x(j-1_i64)) then
                        x(j) = x(j-1_i64)
                    else
                        exit loop_j
                    end if
                end do loop_j
                x(j) = temp
            end if
        end do
    end subroutine isort_i32

    pure subroutine isort2_i32(x, xi)
        integer(kind=i32), intent(inout) :: x(:), xi(size(x, kind=i64))
        integer(kind=i64) :: n, i, j
        integer(kind=i32) :: temp, tempi
        n = size(x, kind=i64)
        do i=2_i64,n
            if (x(i) < x(i-1_i64)) then
                temp = x(i)
                tempi = xi(i)
                loop_j: do j=i,2_i64,-1_i64
                    if (temp < x(j-1_i64)) then
                        x(j) = x(j-1_i64)
                        xi(j) = xi(j-1_i64)
                    else
                        exit loop_j
                    end if
                end do loop_j
                x(j) = temp
                xi(j) = tempi
            end if
        end do
    end subroutine isort2_i32

end module sorting
