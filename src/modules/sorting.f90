module sorting
use, non_intrinsic :: kinds, only: i32, i64, c_bool
implicit none
private
    
    interface is_sorted
        module procedure :: is_sorted_i32
    end interface is_sorted

    public :: is_sorted, qsort, isort, isort2, isort3, qsort2, msort

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
        if (n > 2048_i64) then
            val = [qsort(pack(x, x < x(1_i64))), &
                   x(1_i64), &
                   qsort(pack(x(2_i64:n), x(2_i64:n) >= x(1_i64)))]
        else
            val = isort2(x)
        end if
    end function qsort

    pure recursive function qsort2(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, u, v, w
        n = size(x, kind=i64)
        if (n <= 2048_i64) then
            val = isort2(x)
        else
            val = x
            u = val(1)
            v = val(n/2)
            w = val(n)
            if ((v > u) .neqv. (v > w)) then !! v is median
                val(1) = v
                val(n/2) = u
            else if ((w > u) .neqv. (w > v)) then !! w is median
                val(1) = w
                val(n) = u
            end if !! u is median
            val = [qsort2(pack(val, val < val(1_i64))), &
                   val(1_i64), &
                   qsort2(pack(val(2_i64:n), val(2_i64:n) >= val(1_i64)))]
        end if
    end function qsort2

    pure recursive function msort(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, nl, nr, il, ir, ii, left(size(x, kind=i64)/2_i64), right(size(x, kind=i64)/2_i64+1_i64)
        n = size(x, kind=i64)
        if (n <= 512_i64) then
            val = isort3(x)
        else
            nl = n/2_i64
            nr = n - nl
            left = msort(x(1:nl))
            right = msort(x(nl+1:n))
            il = 1
            ir = 1
            ii = 1
            do while ((il <= nl) .and. (ir <= nr))
                if (left(il) <= right(ir)) then
                    val(ii) = left(il)
                    il = il + 1
                else
                    val(ii) = right(ir)
                    ir = ir + 1
                end if
                ii = ii + 1
            end do
            do while (il <= nl)
                val(ii) = left(il)
                il = il + 1
                ii = ii + 1
            end do
            do while (ir <= nr)
                val(ii) = right(ir)
                ir = ir + 1
                ii = ii + 1
            end do
        end if
    end function msort

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

    pure function isort2(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, i, j
        integer(kind=i32) :: temp
        val = x
        n = size(val, kind=i64)
        do i=2_i64,n
            j = i
            temp = val(j)
            do while ((j>1_i64).and.(temp<val(j-1_i64)))
                val(j) = val(j-1_i64)
                j = j - 1_i64
            end do
            val(j) = temp
        end do
    end function isort2

    pure function isort3(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, i, j
        integer(kind=i32) :: temp
        val = x
        n = size(val, kind=i64)
        do i=2_i64,n
            if (val(i) < val(i-1_i64)) then
                temp = val(i)
                do j=i,2_i64,-1_i64
                    if (temp < val(j-1_i64)) then
                        val(j) = val(j-1_i64)
                    else
                        exit
                    end if
                end do
                val(j) = temp
            end if
        end do
    end function isort3

end module sorting
