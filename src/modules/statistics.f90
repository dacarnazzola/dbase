module statistics
use, non_intrinsic :: kinds, only: i32, i64, sp, dp
use, non_intrinsic :: constants, only: i32_vec_len, i64_vec_len, sp_vec_len, dp_vec_len
implicit none
private

    interface dsum
        module procedure :: dsum_i32
        module procedure :: dsum_i64
        module procedure :: dsum_sp
        module procedure :: dsum_dp
    end interface dsum


    public :: dsum

contains

    pure function dsum_i32(x) result(val)
        integer(kind=i32), intent(in) :: x(:)
        integer(kind=i32) :: val
        integer(kind=i32) :: accumulators(i32_vec_len)
        integer(kind=i64) :: n, i
        val = 0_i32
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 0_i32
            do i = 1_i64, n - sp_vec_len + 1_i64, sp_vec_len
                accumulators = accumulators + x(i:i+sp_vec_len-1_i64)
            end do
            i = n - (n/sp_vec_len)*sp_vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + x(n-i+1_i64:n)
            do i = 1_i64, sp_vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function dsum_i32

    pure function dsum_i64(x) result(val)
        integer(kind=i64), intent(in) :: x(:)
        integer(kind=i64) :: val
        integer(kind=i64) :: n, i, accumulators(i64_vec_len)
        val = 0_i64
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 0_i64
            do i = 1_i64, n - sp_vec_len + 1_i64, sp_vec_len
                accumulators = accumulators + x(i:i+sp_vec_len-1_i64)
            end do
            i = n - (n/sp_vec_len)*sp_vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + x(n-i+1_i64:n)
            do i = 1_i64, sp_vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function dsum_i64

    pure function dsum_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp) :: val
        real(kind=sp) :: accumulators(sp_vec_len)
        integer(kind=i64) :: n, i
        val = 0.0_sp
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 0.0_sp
            do i = 1_i64, n - sp_vec_len + 1_i64, sp_vec_len
                accumulators = accumulators + x(i:i+sp_vec_len-1_i64)
            end do
            i = n - (n/sp_vec_len)*sp_vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + x(n-i+1_i64:n)
            do i = 1_i64, sp_vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function dsum_sp

    pure function dsum_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: val
        real(kind=dp) :: accumulators(dp_vec_len)
        integer(kind=i64) :: n, i
        val = 0.0_dp
        n = size(x, kind=i64)
        if (n > 0_i64) then
            accumulators = 0.0_dp
            do i = 1_i64, n - sp_vec_len + 1_i64, sp_vec_len
                accumulators = accumulators + x(i:i+sp_vec_len-1_i64)
            end do
            i = n - (n/sp_vec_len)*sp_vec_len
            accumulators(1_i64:i) = accumulators(1_i64:i) + x(n-i+1_i64:n)
            do i = 1_i64, sp_vec_len, 2_i64
                val = val + accumulators(i) + accumulators(i + 1_i64)
            end do
        end if
    end function dsum_dp

end module statistics
