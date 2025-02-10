program test_sorting
use, non_intrinsic :: kinds, only: stdout, i32, i64, c_bool
use, non_intrinsic :: random, only: random_uniform
use, non_intrinsic :: sorting, only: is_sorted, qsort, isort
use, non_intrinsic :: timing, only: timer, tic, toc, get_elapsed
implicit none

    integer(kind=i64), parameter :: n = 10_i64

    integer(kind=i32) :: arr(n), arr_dupe(n), i
    logical(kind=c_bool) :: sorted
    integer(kind=i32), allocatable :: x(:), ns(:), x2(:)
    type(timer) :: bench(2)

    call random_uniform(arr, -100_i32, 100_i32)
    arr_dupe = arr
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = isort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = qsort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    if (allocated(ns)) deallocate(ns); allocate(ns(0))
    do i=1_i32,6_i32
        ns = [ns, 10_i32**i]
    end do
    do i=0_i32,20_i32
        ns = [ns, 2_i32**i]
    end do
    ns = qsort(ns)
    write(unit=stdout, fmt='(a,*(i0," "))') 'TEST_SORTING :: ns: ',ns

    write(unit=stdout, fmt='(a)') 'TEST_SORTING :: benchmark qsort vs isort'
    write(unit=stdout, fmt='(3(a13))') 'N Elements','| QSORT Time','| ISORT Time'
    do i=1_i32,size(ns, kind=i32)
        if (allocated(x)) deallocate(x)
        if (allocated(x2)) deallocate(x2)
        allocate(x(ns(i)), x2(ns(i)))
        call random_uniform(x, 0_i32, huge(1_i32))
        x2 = x
        call tic(bench(1))
        x = qsort(x2)
        call toc(bench(1))
        if (.not.is_sorted(x)) error stop 'x not sorted after QSORT'
        call tic(bench(2))
        x = isort(x2)
        call toc(bench(2))
        if (.not.is_sorted(x)) error stop 'x not sorted after ISORT'
        write(unit=stdout, fmt='(i13,e13.6,e13.6)') ns(i),get_elapsed(bench)
    end do

    write(unit=stdout, fmt='(a,2(i0," "))') 'COMPLETE !! ',minval(x),maxval(x)

end program test_sorting
