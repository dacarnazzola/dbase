module sorting
use, non_intrinsic :: kinds, only: i32, i64, c_bool
implicit none
private
    
    interface is_sorted
        module procedure :: is_sorted_i32
    end interface is_sorted

    public :: is_sorted, qsort, isort

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

    pure recursive function qsort(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        if (n > 256_i64) then
            val = [qsort(pack(x, x < x(1_i64))), &
                   x(1_i64), &
                   qsort(pack(x(2_i64:n), x(2_i64:n) >= x(1_i64)))]
        else
            val = isort(x)
        end if
    end function qsort

    pure function isort(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, i, j
        val = x
        n = size(val, kind=i64)
        loop_i: do i=2_i64,n
            if (val(i) < val(i-1_i64)) then
                do j=1_i64,i-1_i64
                    if (val(i) < val(j)) then
                        val(1_i64:i) = [val(1_i64:j-1_i64), val(i), val(j:i-1_i64)]
                        cycle loop_i
                    end if
                end do
            end if 
        end do loop_i
    end function isort

end module sorting
