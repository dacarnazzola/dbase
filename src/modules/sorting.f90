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

    public :: is_sorted, qsort, isort, isort2, isort3, qsort2, msort, msort2, msort3, qsort3, qsorts, qsorts2

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

    pure recursive function qsort(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        if (n > i32_vec_len) then
            val = [qsort(pack(x, x < x(1_i64))), &
                   x(1_i64), &
                   qsort(pack(x(2_i64:n), x(2_i64:n) >= x(1_i64)))]
        else
            val = isort3(x)
        end if
    end function qsort

    pure recursive function qsort2(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n
        integer(kind=i32) :: u, v, w
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            val = isort3(x)
        else
            val = x
            u = val(1)
            v = val(n/2)
            w = val(n)
            if ((v > u) .neqv. (v > w)) then !! v is median
                val(1_i64) = v
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

    pure recursive function qsort3(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, i
        integer(kind=i32) :: u, v, w
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            val = isort3(x)
        else
            val = x
            u = val(1)
            v = val(n/2)
            w = val(n)
            if ((v > u) .neqv. (v > w)) then !! v is median
                val(1_i64) = v
                val(n/2_i64) = u
            else if ((w > u) .neqv. (w > v)) then !! w is median
                val(1) = w
                val(n) = u
            end if !! u is median
            v = 1_i64
            do i=2_i64,n
                if (val(i) < val(1_i64)) then
                    v = v + 1_i32
                    u = val(i)
                    val(i) = val(v)
                    val(v) = u
                end if
            end do
            u = val(v)
            val(v) = val(1_i64)
            val(1_i64) = u
            val(1_i64:v-1_i64) = qsort3(val(1_i64:v-1_i64))
            val(v+1_i64:n) = qsort3(val(v+1_i64:n))
        end if
    end function qsort3

    pure recursive subroutine qsorts2(x, xi)
        integer(kind=i32), intent(inout) :: x(:), xi(size(x, kind=i64))
        integer(kind=i64) :: n, i, ii
        integer(kind=i32) :: u, v, w, uxi, vxi, wxi
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            call isorts2(x, xi)
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
            call qsorts2(x(1_i64:ii-1_i64), xi(1_i64:ii-1_i64))
            call qsorts2(x(ii+1_i64:n), xi(ii+1_i64:n))
        end if
    end subroutine qsorts2

    pure recursive subroutine qsorts(x)
        integer(kind=i32), intent(inout) :: x(:)
        integer(kind=i64) :: n, i
        integer(kind=i32) :: u, v, w
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            call isorts(x)
        else
            u = x(1)
            v = x(n/2)
            w = x(n)
            if ((v > u) .neqv. (v > w)) then !! v is median
                x(1) = v
                x(n/2) = u
            else if ((w > u) .neqv. (w > v)) then !! w is median
                x(1) = w
                x(n) = u
            end if !! u is median
            v = 1_i32
            do i=2_i64,n
                if (x(i) < x(1_i64)) then
                    v = v + 1_i32
                    u = x(i)
                    x(i) = x(v)
                    x(v) = u
                end if
            end do
            u = x(v)
            x(v) = x(1_i64)
            x(1_i64) = u
            call qsorts(x(1_i64:v-1_i64))
            call qsorts(x(v+1_i64:n))
        end if
    end subroutine qsorts

    pure recursive function msort(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, nl, nr, il, ir, ii
        integer(kind=i32) :: left(size(x, kind=i64)/2_i64), right(size(x, kind=i64)/2_i64+1_i64)
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            val = isort3(x)
        else
            nl = n/2_i64
            nr = n - nl
            left = msort(x(1_i64:nl))
            right = msort(x(nl+1_i64:n))
            il = 1_i64
            ir = 1_i64
            ii = 1_i64
            do while ((il <= nl) .and. (ir <= nr))
                if (left(il) <= right(ir)) then
                    val(ii) = left(il)
                    il = il + 1_i64
                else
                    val(ii) = right(ir)
                    ir = ir + 1_i64
                end if
                ii = ii + 1_i64
            end do
            do while (il <= nl)
                val(ii) = left(il)
                il = il + 1_i64
                ii = ii + 1_i64
            end do
            do while (ir <= nr)
                val(ii) = right(ir)
                ir = ir + 1_i64
                ii = ii + 1_i64
            end do
        end if
    end function msort

    pure recursive function msort2(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, nl, nr, il, ir, ii, ii_presorted
        integer(kind=i32) :: left(size(x, kind=i64)/2_i64), right(size(x, kind=i64)/2_i64+1_i64)
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            val = isort3(x)
        else
            ii_presorted = 1_i64
            do il=2_i64,n
                if (x(il) < x(il-1_i64)) then !! out of place
                    ii_presorted = il - 1_i64
                    exit
                end if
            end do
            if (il >= (n - 1_i64)) then
                val = isort3(x)
                return
            end if
            nl = (n-ii_presorted+1_i64)/2_i64
            nr = (n-ii_presorted+1_i64) - nl
            left = msort2(x(ii_presorted:ii_presorted+nl-1_i64))
            right = msort2(x(ii_presorted+nl:n))
            il = 1_i64
            ir = 1_i64
            ii = ii_presorted
            do while ((il <= nl) .and. (ir <= nr))
                if (left(il) <= right(ir)) then
                    val(ii) = left(il)
                    il = il + 1_i64
                else
                    val(ii) = right(ir)
                    ir = ir + 1_i64
                end if
                ii = ii + 1_i64
            end do
            do while (il <= nl)
                val(ii) = left(il)
                il = il + 1_i64
                ii = ii + 1_i64
            end do
            do while (ir <= nr)
                val(ii) = right(ir)
                ir = ir + 1_i64
                ii = ii + 1_i64
            end do
            val = isort3(val)
        end if
    end function msort2

    pure recursive function msort3(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val(size(x, kind=i64))
        integer(kind=i64) :: n, nl, nr, il, ir, ii
        integer(kind=i32) :: left(size(x, kind=i64)/2_i64), right(size(x, kind=i64)/2_i64+1_i64)
        n = size(x, kind=i64)
        if (n <= i32_vec_len) then
            val = isort3(x)
        else
            nl = n/2_i64
            nr = n - nl
            left = msort3(x(1_i64:nl))
            right = msort3(x(nl+1_i64:n))
            il = 1_i64
            ir = 1_i64
            do ii=1_i64,n
                if ((il <= nl) .and. (ir <= nr)) then
                    if (left(il) <= right(ir)) then
                        val(ii) = left(il)
                        il = il + 1_i64
                    else
                        val(ii) = right(ir)
                        ir = ir + 1_i64
                    end if
                else if (ir <= nr) then
                    val(ii:n) = right(ir:nr)
                    return
                else
                    val(ii:n) = left(il:nl)
                    return
                end if
            end do
        end if
    end function msort3

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

    pure subroutine isorts2(x, xi)
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
    end subroutine isorts2

    pure subroutine isorts(x)
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
    end subroutine isorts

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
                loop_j: do j=i,2_i64,-1_i64
                    if (temp < val(j-1_i64)) then
                        val(j) = val(j-1_i64)
                    else
                        exit loop_j
                    end if
                end do loop_j
                val(j) = temp
            end if
        end do
    end function isort3

end module sorting
