program test_sorting
use, non_intrinsic :: kinds, only: stdout, i32, i64, c_bool, dp
use, non_intrinsic :: random, only: random_uniform
use, non_intrinsic :: sorting, only: is_sorted, sort
use, non_intrinsic :: timing, only: timer, tic, toc, get_elapsed
implicit none

    integer(kind=i64), parameter :: n = 32_i64
    integer(kind=i32), parameter :: k_max = 1000_i32

    integer(kind=i32) :: arr(n), arr_dupe(n), i, k, arri(n)
    logical(kind=c_bool) :: sorted
    integer(kind=i32), allocatable :: x(:), ns(:), x2(:), xi(:), xi2(:)
    type(timer) :: bench(2)
    real(kind=dp) :: times(size(bench))

    arri = [(i, i=1,n)]
    call random_uniform(arr, -100_i32, 100_i32)
    arr_dupe = arr
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING ::      arr, sorted=',sorted,' :: ',arr

    call sort(arr)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: SORT arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    call sort(arr, arri)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: SORT arr, sorted=',sorted,' :: ',arr

    if (allocated(ns)) deallocate(ns); allocate(ns(0))
    do i=1_i32,6_i32
        ns = [ns, 10_i32**i]
    end do
    do i=0_i32,20_i32
        ns = [ns, 2_i32**i]
    end do
    call sort(ns)
    write(unit=stdout, fmt='(a,*(i0," "))') 'TEST_SORTING :: ns: ',ns

    write(unit=stdout, fmt='(a)') 'TEST_SORTING :: benchmark qsort vs isort'
    write(unit=stdout, fmt='(*(a22))') 'N Elements |','SORT Time |','SORT2 Time |'
    do i=1_i32,size(ns, kind=i32)
        if (allocated(x)) deallocate(x)
        if (allocated(x2)) deallocate(x2)
        if (allocated(xi)) deallocate(xi)
        if (allocated(xi2)) deallocate(xi2)
        allocate(x(ns(i)), x2(ns(i)), xi(ns(i)), xi2(ns(i)))
        xi = [(k, k=1,ns(i))]
        xi2 = xi
        times = 0.0_dp
        do k=1_i32,k_max
            call random_uniform(x, 0_i32, huge(1_i32))
            x2 = x
            call tic(bench(1))
            call sort(x)
            call toc(bench(1))
            if (.not.is_sorted(x)) error stop 'x not sorted after SORT (values)'
            x = x2
            xi = xi2
            call tic(bench(2))
            call sort(x, xi)
            call toc(bench(2))
            if (.not.is_sorted(x) .or. any(x /= x2(xi))) error stop 'x not sorted after SORT (values and indices)'
            times = times + get_elapsed(bench)
        end do
        times = times/real(k_max, kind=dp)
        write(unit=stdout, fmt='(i22,*(e22.15))') ns(i),times
    end do

    write(unit=stdout, fmt='(a,2(i0," "))') 'COMPLETE !! ',minval(x),maxval(x)

end program test_sorting
