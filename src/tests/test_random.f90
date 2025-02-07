program test_random
use, non_intrinsic :: kinds, only: stdout, i32, i64, sp, dp
use, non_intrinsic :: timing, only: timer, tic, toc, get_elapsed
use, non_intrinsic :: random, only: random_uniform
implicit none
    
    integer(kind=i32), parameter :: trials = 1000_i32

    integer(kind=i32) :: ni, n, i
    integer(kind=i32), allocatable :: xi32(:)
    integer(kind=i64), allocatable :: xi64(:)
    real(kind=sp), allocatable :: xsp(:)
    real(kind=dp), allocatable :: xdp(:)
    real(kind=dp) :: times(4)
    type(timer) :: bench(4)

    do ni=8_i32,20_i32,4_i32
        n = 2_i32**ni
        if (allocated(xi32)) deallocate(xi32)
        if (allocated(xi64)) deallocate(xi64)
        if (allocated(xsp)) deallocate(xsp)
        if (allocated(xdp)) deallocate(xdp)
        allocate(xi32(n), xi64(n), xsp(n), xdp(n))
        times = 0.0_dp
        do i=1_i32,trials

            call tic(bench(1))
            call random_uniform(xi32, 1_i32, 6_i32)
            call toc(bench(1))

            call tic(bench(2))
            call random_uniform(xi64, 1_i64, 6_i64)
            call toc(bench(2))

            call tic(bench(3))
            call random_uniform(xsp, 1.0_sp, 6.0_sp)
            call toc(bench(3))

            call tic(bench(4))
            call random_uniform(xdp, 1.0_dp, 6.0_dp)
            call toc(bench(4))

            times = times + get_elapsed(bench)

            if (i == trials) then
                write(unit=stdout, fmt='(a,i0,a)') 'TEST_RANDOM :: generating 2**',ni,' elements'
                write(unit=stdout, fmt='(a16,i7,a5,i3,a,f0.6,a)') 'TEST_RANDOM :: xi32(',n,'): ',xi32(n),' ( ',times(1),' sec)'
                write(unit=stdout, fmt='(a16,i7,a5,i3,a,f0.6,a)') 'TEST_RANDOM :: xi64(',n,'): ',xi64(n),' ( ',times(2),' sec)'
                write(unit=stdout, fmt='(a16,i7,a5,f3.1,a,f0.6,a)') 'TEST_RANDOM :: xsp(',n,'): ',xsp(n),' ( ',times(3),' sec)'
                write(unit=stdout, fmt='(a16,i7,a5,f3.1,a,f0.6,a,/)') 'TEST_RANDOM :: xdp(',n,'): ',xdp(n),' ( ',times(4),' sec)'
            end if

        end do
    end do

end program test_random
