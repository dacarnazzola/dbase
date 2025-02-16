module statistics
use, non_intrinsic :: kinds, only: i32, i64, sp, dp, c_bool
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: constants, only: i32_vec_len, i64_vec_len, sp_vec_len, dp_vec_len, twopi_sp, twopi_dp
implicit none
private

    interface dsum
        module procedure :: dsum_i32
        module procedure :: dsum_i64
        module procedure :: dsum_sp
        module procedure :: dsum_dp
    end interface dsum

    interface avg
        module procedure :: avg_sp
        module procedure :: avg_dp
    end interface avg

    interface std
        module procedure :: std_sp
        module procedure :: std_dp
    end interface std

    interface normal_pdf
        module procedure :: normal_pdf_sp
        module procedure :: normal_pdf_dp
    end interface normal_pdf

    public :: dsum, avg, std, normal_pdf

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

    pure function avg_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp) :: val
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: avg function invalid for array with length < 1')
        val = dsum(x)/real(n, kind=sp)
    end function avg_sp

    pure function avg_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: val
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 1_i64, kind=c_bool), &
                                   'module STATISTICS :: avg function invalid for array with length < 1')
        val = dsum(x)/real(n, kind=dp)
    end function avg_dp

    pure function std_sp(x) result(val)
        real(kind=sp), intent(in) :: x(:)
        real(kind=sp) :: val
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: std function invalid for array with length < 2')
        val = sqrt(dsum((x - avg(x))**2_i32)/real(n - 1_i64, kind=sp))
    end function std_sp

    pure function std_dp(x) result(val)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: val
        integer(kind=i64) :: n
        n = size(x, kind=i64)
        call debug_error_condition(logical(n < 2_i64, kind=c_bool), &
                                   'module STATISTICS :: std function invalid for array with length < 2')
        val = sqrt(dsum((x - avg(x))**2_i32)/real(n - 1_i64, kind=dp))
    end function std_dp

    pure elemental function normal_pdf_sp(x, mu, sig) result(val)
        real(kind=sp), intent(in) :: x, mu, sig
        real(kind=sp) :: val
        real(kind=sp) :: sig2
        call debug_error_condition(nearly(sig, 0.0_sp), &
                                   'module STATISTICS :: normal_pdf function invalid for input with sig == 0.0')
        sig2 = sig*sig
        val = exp(-(x - mu)**2_i32/(2.0_sp*sig2))/sqrt(twopi_sp*sig2)
    end function normal_pdf_sp

    pure elemental function normal_pdf_dp(x, mu, sig) result(val)
        real(kind=dp), intent(in) :: x, mu, sig
        real(kind=dp) :: val
        real(kind=dp) :: sig2
        call debug_error_condition(nearly(sig, 0.0_dp), &
                                   'module STATISTICS :: normal_pdf function invalid for input with sig == 0.0')
        sig2 = sig*sig
        val = exp(-(x - mu)**2_i32/(2.0_dp*sig2))/sqrt(twopi_dp*sig2)
    end function normal_pdf_dp

end module statistics
