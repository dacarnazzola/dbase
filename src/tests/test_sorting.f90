program test_sorting
use, non_intrinsic :: kinds, only: stdout, i32, i64, c_bool, dp
use, non_intrinsic :: random, only: random_uniform
use, non_intrinsic :: sorting, only: is_sorted, qsort, isort, isort2, isort3, qsort2, msort, msort2, msort3, qsort3, qsorts, qsorts2
use, non_intrinsic :: timing, only: timer, tic, toc, get_elapsed
implicit none

    integer(kind=i64), parameter :: n = 32_i64
    integer(kind=i32), parameter :: k_max = 100_i32

    integer(kind=i32) :: arr(n), arr_dupe(n), i, k, arri(n)
    logical(kind=c_bool) :: sorted
    integer(kind=i32), allocatable :: x(:), ns(:), x2(:), xi(:), xi2(:)
    type(timer) :: bench(8)
    real(kind=dp) :: times(size(bench))

    arri = [(i, i=1,n)]
    call random_uniform(arr, -100_i32, 100_i32)
    arr_dupe = arr
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING ::         arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = isort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: ISORT   arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = isort2(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: ISORT2  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = isort3(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: ISORT3  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = qsort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: QSORT   arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = qsort2(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: QSORT2  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = qsort3(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: QSORT3  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = msort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: MSORT   arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = msort2(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: MSORT2  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = msort3(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: MSORT3  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    call qsorts(arr)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: QSORTS  arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    call qsorts2(arr, arri)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'TEST_SORTING :: QSORTS2 arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    if (allocated(ns)) deallocate(ns); allocate(ns(0))
    do i=1_i32,5_i32
        do k=1_i32,9_i32
            ns = [ns, k*10_i32**i]
        end do
    end do
    do i=0_i32,20_i32
        ns = [ns, 2_i32**i]
    end do
    ns = msort(ns)
    write(unit=stdout, fmt='(a,*(i0," "))') 'TEST_SORTING :: ns: ',ns

    write(unit=stdout, fmt='(a)') 'TEST_SORTING :: benchmark qsort vs isort'
    write(unit=stdout, fmt='(*(a16))') 'N Elements','| QSORT Time','| QSORT2 Time','| QSORT3 Time', &
                                                    '| MSORT Time','| MSORT2 Time','| MSORT3 Time', &
                                                    '| QSORTS Time','| QSORTS2 Time'
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
            x = qsort(x)
            call toc(bench(1))
            if (.not.is_sorted(x)) error stop 'x not sorted after QSORT'
            x = x2
            call tic(bench(2))
            x = qsort2(x)
            call toc(bench(2))
            if (.not.is_sorted(x)) error stop 'x not sorted after QSORT2'
            x = x2
            call tic(bench(3))
            x = qsort3(x)
            call toc(bench(3))
            if (.not.is_sorted(x)) error stop 'x not sorted after QSORT3'
            x = x2
            call tic(bench(4))
            x = msort(x)
            call toc(bench(4))
            if (.not.is_sorted(x)) error stop 'x not sorted after MSORT'
            x = x2
            call tic(bench(5))
            x = msort2(x)
            call toc(bench(5))
            if (.not.is_sorted(x)) error stop 'x not sorted after MSORT2'
            x = x2
            call tic(bench(6))
            x = msort3(x)
            call toc(bench(6))
            if (.not.is_sorted(x)) error stop 'x not sorted after MSORT3'
            x = x2
            call tic(bench(7))
            call qsorts(x)
            call toc(bench(7))
            if (.not.is_sorted(x)) error stop 'x not sorted after QSORTS'
            x = x2
            xi = xi2
            call tic(bench(8))
            call qsorts2(x, xi)
            call toc(bench(8))
            if (.not.is_sorted(x) .or. any(x /= x2(xi))) error stop 'x not sorted after QSORTS2'
            times = times + get_elapsed(bench)
        end do
        times = times/real(k_max, kind=dp)
        times = times/minval(times)
        write(unit=stdout, fmt='(i16,*(f16.2))') ns(i),times
    end do

    write(unit=stdout, fmt='(a,2(i0," "))') 'COMPLETE !! ',minval(x),maxval(x)

end program test_sorting
