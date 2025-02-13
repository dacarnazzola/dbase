program test_sorting
use, non_intrinsic :: kinds, only: stdout, i32, i64, c_bool
use, non_intrinsic :: random, only: random_uniform
use, non_intrinsic :: sorting, only: is_sorted, qsort, isort, isort2, isort3, qsort2, msort
use, non_intrinsic :: timing, only: timer, tic, toc, get_elapsed
implicit none

    integer(kind=i64), parameter :: n = 10_i64

    integer(kind=i32) :: arr(n), arr_dupe(n), i
    logical(kind=c_bool) :: sorted
    integer(kind=i32), allocatable :: x(:), ns(:), x2(:)
    type(timer) :: bench(6)

    call random_uniform(arr, -100_i32, 100_i32)
    arr_dupe = arr
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING ::        arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = isort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: ISORT  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = isort2(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: ISORT2 arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = isort3(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: ISORT3 arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = qsort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: QSORT  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = qsort2(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: QSORT2 arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = msort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: MSORT  arr, sorted=',sorted,' :: ',arr
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
    write(unit=stdout, fmt='(*(a16))') 'N Elements','| QSORT Time','| ISORT Time','| ISORT2 Time', &
                                                    '| ISORT3 Time','| QSORT2 Time','| MSORT Time'
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
        call tic(bench(3))
        x = isort2(x2)
        call toc(bench(3))
        if (.not.is_sorted(x)) error stop 'x not sorted after ISORT2'
        call tic(bench(4))
        x = isort3(x2)
        call toc(bench(4))
        if (.not.is_sorted(x)) error stop 'x not sorted after ISORT3'
        call tic(bench(5))
        x = qsort2(x2)
        call toc(bench(5))
        if (.not.is_sorted(x)) error stop 'x not sorted after QSORT2'
        call tic(bench(6))
        x = msort(x2)
        call toc(bench(6))
        if (.not.is_sorted(x)) error stop 'x not sorted after MSORT'
        write(unit=stdout, fmt='(i16,*(e16.6))') ns(i),get_elapsed(bench)
    end do

    write(unit=stdout, fmt='(a,2(i0," "))') 'COMPLETE !! ',minval(x),maxval(x)

end program test_sorting
