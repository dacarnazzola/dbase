program test_debug
use, non_intrinsic :: kinds, only: stdout
use, non_intrinsic :: constants, only: debug
use, non_intrinsic :: timing, only: timer, tic, toc, get_elapsed
implicit none

    type(timer) :: bench

    call tic(bench)
    write(unit=stdout, fmt='(a,l1)', advance='no') 'TEST_DEBUG :: status of logical parameter DEBUG = ',debug
    call toc(bench)
    write(unit=stdout, fmt='(a,e22.15,a)') ' (',get_elapsed(bench),' sec)'

    write(unit=stdout, fmt='(a)', advance='no') 'TEST_DEBUG :: timing raw tic/toc ... '
    call tic(bench)
    call toc(bench)
    write(unit=stdout, fmt='(a,e22.15,a)') '(',get_elapsed(bench),' sec)'

end program test_debug
