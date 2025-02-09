module random
use, non_intrinsic :: kinds, only: i32, i64, sp, dp
use, non_intrinsic :: constants, only: twopi_sp, twopi_dp
implicit none
private

    interface random_uniform
        module procedure :: random_uniform_i32
        module procedure :: random_uniform_i64
        module procedure :: random_uniform_sp
        module procedure :: random_uniform_dp
    end interface random_uniform

    interface random_normal
        module procedure :: random_normal_sp
        module procedure :: random_normal_dp
    end interface random_normal

    public :: random_uniform, random_normal

contains

    impure subroutine random_uniform_i32(vout, vmin, vmax)
        integer(kind=i32), intent(out) :: vout(:)
        integer(kind=i32), intent(in) :: vmin, vmax
        real(kind=sp) :: work(size(vout, kind=i64))
        call random_number(work)
        vout = int(work*real(vmax - vmin + 1_i64, kind=sp), kind=i32) + vmin
    end subroutine random_uniform_i32

    impure subroutine random_uniform_i64(vout, vmin, vmax)
        integer(kind=i64), intent(out) :: vout(:)
        integer(kind=i64), intent(in) :: vmin, vmax
        real(kind=sp) :: work(size(vout, kind=i64))
        call random_number(work)
        vout = int(work*real(vmax - vmin + 1_i64, kind=sp), kind=i64) + vmin
    end subroutine random_uniform_i64

    impure subroutine random_uniform_sp(vout, vmin, vmax)
        real(kind=sp), intent(out) :: vout(:)
        real(kind=sp), intent(in) :: vmin, vmax
        call random_number(vout)
        vout = vout*(vmax - vmin) + vmin
    end subroutine random_uniform_sp

    impure subroutine random_uniform_dp(vout, vmin, vmax)
        real(kind=dp), intent(out) :: vout(:)
        real(kind=dp), intent(in) :: vmin, vmax
        call random_number(vout)
        vout = vout*(vmax - vmin) + vmin
    end subroutine random_uniform_dp

    impure subroutine random_normal_sp(vout, mu, sig)
        real(kind=sp), intent(out) :: vout(:)
        real(kind=sp), intent(in) :: mu, sig
        real(kind=sp) :: u((size(vout,kind=i64)+1_i64)/2_i64), r((size(vout,kind=i64)+1_i64)/2_i64)
        integer(kind=i64) :: n
        call random_number(r)
        r = 1.0_sp - r
        r = sqrt(-2.0_sp*log(r))
        call random_number(u)
        n = size(vout, kind=i64)
        vout(1_i64:(n/2_i64)) = mu + sig*r(1_i64:(n/2_i64))*cos(twopi_sp*u(1_i64:(n/2_i64)))
        vout((n/2_i64+1):n) = mu + sig*r*sin(twopi_sp*u)
    end subroutine random_normal_sp

    impure subroutine random_normal_dp(vout, mu, sig)
        real(kind=dp), intent(out) :: vout(:)
        real(kind=dp), intent(in) :: mu, sig
        real(kind=dp) :: u((size(vout,kind=i64)+1_i64)/2_i64), r((size(vout,kind=i64)+1_i64)/2_i64)
        integer(kind=i64) :: n
        call random_number(r)
        r = 1.0_dp - r
        r = sqrt(-2.0_dp*log(r))
        call random_number(u)
        n = size(vout, kind=i64)
        vout(1_i64:(n/2_i64)) = mu + sig*r(1_i64:(n/2_i64))*cos(twopi_dp*u(1_i64:(n/2_i64)))
        vout((n/2_i64+1):n) = mu + sig*r*sin(twopi_dp*u)
    end subroutine random_normal_dp

end module random
