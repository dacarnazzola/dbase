module constants
use, non_intrinsic :: kinds, only: sp, dp, c_bool
implicit none
private

#ifdef __COMPILE_FOR_DEBUG__
    logical(kind=c_bool), parameter :: debug = .true.
#else
    logical(kind=c_bool), parameter :: debug = .false.
#endif

    real(kind=sp), parameter :: pi_sp = acos(-1.0_sp)
    real(kind=dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: debug, pi_sp, pi_dp

end module constants
