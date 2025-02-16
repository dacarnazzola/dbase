module vector_math
use, non_intrinsic :: kinds, only: sp, dp
use, non_intrinsic :: statistics, only: dsum
use, non_intrinsic :: system, only: debug_error_condition, nearly
implicit none
private

    interface vmag2
        module procedure :: vmag2_sp
        module procedure :: vmag2_dp
    end interface vmag2

    interface vmag3
        module procedure :: vmag3_sp
        module procedure :: vmag3_dp
    end interface vmag3

    interface vmag
        module procedure :: vmag_sp
        module procedure :: vmag_dp
    end interface vmag

    interface vunit
        module procedure :: vunit_sp
        module procedure :: vunit_dp
    end interface vunit

    public :: vmag2, vmag3, vmag, vunit

contains

    pure function vmag2_sp(v) result(val)
        real(kind=sp), intent(in) :: v(2)
        real(kind=sp) :: val
        val = sqrt(v(1)*v(1) + v(2)*v(2))
    end function vmag2_sp

    pure function vmag2_dp(v) result(val)
        real(kind=dp), intent(in) :: v(2)
        real(kind=dp) :: val
        val = sqrt(v(1)*v(1) + v(2)*v(2))
    end function vmag2_dp

    pure function vmag3_sp(v) result(val)
        real(kind=sp), intent(in) :: v(3)
        real(kind=sp) :: val
        val = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
    end function vmag3_sp

    pure function vmag3_dp(v) result(val)
        real(kind=dp), intent(in) :: v(3)
        real(kind=dp) :: val
        val = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
    end function vmag3_dp

    pure function vmag_sp(v) result(val)
        real(kind=sp), intent(in) :: v(:)
        real(kind=sp) :: val
        val = sqrt(dsum(v*v))
    end function vmag_sp

    pure function vmag_dp(v) result(val)
        real(kind=dp), intent(in) :: v(:)
        real(kind=dp) :: val
        val = sqrt(dsum(v*v))
    end function vmag_dp

    pure subroutine vunit_sp(v, vhat)
        real(kind=sp), intent(in) :: v(:)
        real(kind=sp), intent(out) :: vhat(:)
        real(kind=sp) :: mag
        mag = vmag(v)
        call debug_error_condition(nearly(mag, 0.0_sp), &
                                   'module VECTOR_MATH :: vunit subroutine invalid for vectors with magnitude near 0.0')
        vhat = v/mag
    end subroutine vunit_sp

    pure subroutine vunit_dp(v, vhat)
        real(kind=dp), intent(in) :: v(:)
        real(kind=dp), intent(out) :: vhat(:)
        real(kind=dp) :: mag
        mag = vmag(v)
        call debug_error_condition(nearly(mag, 0.0_dp), &
                                   'module VECTOR_MATH :: vunit subroutine invalid for vectors with magnitude near 0.0')
        vhat = v/vmag(v)
    end subroutine vunit_dp

end module vector_math
