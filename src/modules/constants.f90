module constants
use, non_intrinsic :: kinds, only: i32, sp, dp, c_bool
implicit none
private

#ifdef __COMPILE_FOR_DEBUG__
    logical(kind=c_bool), parameter :: debug = .true.
#else
    logical(kind=c_bool), parameter :: debug = .false.
#endif

    integer(kind=i32), parameter :: i32_vec_len = 64_i32
    integer(kind=i32), parameter :: i64_vec_len = 64_i32
    integer(kind=i32), parameter :: sp_vec_len  = 64_i32
    integer(kind=i32), parameter :: dp_vec_len  = 64_i32

    real(kind=sp), parameter :: pi_sp = acos(-1.0_sp)
    real(kind=dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: debug, pi_sp, pi_dp, i32_vec_len, i64_vec_len, sp_vec_len, dp_vec_len

end module constants
