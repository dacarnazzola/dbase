module math_utils
use, non_intrinsic :: kinds, only: i32, i64, sp, dp
implicit none
private

    interface clamp
        module procedure :: clamp_i32
        module procedure :: clamp_i64
        module procedure :: clamp_sp
        module procedure :: clamp_dp
    end interface clamp

    public :: clamp

contains

    pure elemental subroutine clamp_i32(x, xlo, xhi)
        integer(kind=i32), intent(inout) :: x
        integer(kind=i32), intent(in) :: xlo, xhi
        x = min(max(x, xlo), xhi)
    end subroutine clamp_i32

    pure elemental subroutine clamp_i64(x, xlo, xhi)
        integer(kind=i64), intent(inout) :: x
        integer(kind=i64), intent(in) :: xlo, xhi
        x = min(max(x, xlo), xhi)
    end subroutine clamp_i64

    pure elemental subroutine clamp_sp(x, xlo, xhi)
        real(kind=sp), intent(inout) :: x
        real(kind=sp), intent(in) :: xlo, xhi
        x = min(max(x, xlo), xhi)
    end subroutine clamp_sp

    pure elemental subroutine clamp_dp(x, xlo, xhi)
        real(kind=dp), intent(inout) :: x
        real(kind=dp), intent(in) :: xlo, xhi
        x = min(max(x, xlo), xhi)
    end subroutine clamp_dp

end module math_utils
