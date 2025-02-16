module vector_math
use, non_intrinsic :: kinds, only: sp
use, non_intrinsic :: statistics, only: dsum
implicit none
private

    interface vmag2
        module procedure :: vmag2_sp
    end interface vmag2

    interface vmag3
        module procedure :: vmag3_sp
    end interface vmag3

    interface vmag
        module procedure :: vmag_sp
    end interface vmag

    public :: vmag2, vmag3, vmag

contains

    pure function vmag2_sp(v) result(val)
        real(kind=sp), intent(in) :: v(2)
        real(kind=sp) :: val
        val = sqrt(v(1)*v(1) + v(2)*v(2))
    end function vmag2_sp

    pure function vmag3_sp(v) result(val)
        real(kind=sp), intent(in) :: v(3)
        real(kind=sp) :: val
        val = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
    end function vmag3_sp

    pure function vmag_sp(v) result(val)
        real(kind=sp), intent(in) :: v(:)
        real(kind=sp) :: val
        val = sqrt(dsum(v*v))
    end function vmag_sp

end module vector_math
