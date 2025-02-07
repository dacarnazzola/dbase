module random
use, non_intrinsic :: kinds, only: i32, i64, sp, dp
implicit none
private

    interface random_uniform
        module procedure :: random_uniform_i32
        module procedure :: random_uniform_i64
        module procedure :: random_uniform_sp
        module procedure :: random_uniform_dp
    end interface random_uniform

    public :: random_uniform

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

end module random
