program test_statistics
use, non_intrinsic :: kinds, only: stdout, i32, i64, sp, dp
use, non_intrinsic :: statistics, only: dsum
use, non_intrinsic :: timing, only: timer, tic, toc, get_elapsed
implicit none

    integer(kind=i32), parameter :: n = 2_i32**20_i32 - 1_i32

    integer(kind=i32) :: xi32(n), yi32
    integer(kind=i64) :: xi64(n), yi64
    real(kind=sp) :: xsp(n), ysp
    real(kind=dp) :: xdp(n), ydp
    type(timer) :: bench

    xi32 = 1_i32
    call tic(bench)
    yi32 = dsum(xi32)
    call toc(bench)
    write(unit=stdout, fmt='(a,i0,a,i0,a,f0.6,a)') 'TEST_STATISTICS :: dsum(xi32[',n,']): ',yi32,' (',get_elapsed(bench),' sec)'

    xi64 = 1_i64
    call tic(bench)
    yi64 = dsum(xi64)
    call toc(bench)
    write(unit=stdout, fmt='(a,i0,a,i0,a,f0.6,a)') 'TEST_STATISTICS :: dsum(xi64[',n,']): ',yi64,' (',get_elapsed(bench),' sec)'

    xsp = 1.0_sp
    call tic(bench)
    ysp = dsum(xsp)
    call toc(bench)
    write(unit=stdout, fmt='(a,i0,a,f0.1,a,f0.6,a)') 'TEST_STATISTICS :: dsum(xsp[',n,']): ',ysp,' (',get_elapsed(bench),' sec)'

    xdp = 1.0_dp
    call tic(bench)
    ydp = dsum(xdp)
    call toc(bench)
    write(unit=stdout, fmt='(a,i0,a,f0.1,a,f0.6,a)') 'TEST_STATISTICS :: dsum(xdp[',n,']): ',ydp,' (',get_elapsed(bench),' sec)'

end program test_statistics
