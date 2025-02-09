program test_sorting
use, non_intrinsic :: kinds, only: stdout, i32, i64, c_bool
use, non_intrinsic :: random, only: random_uniform
use, non_intrinsic :: sorting, only: is_sorted, qsort, isort
implicit none

    integer(kind=i64), parameter :: n = 10_i64

    integer(kind=i32) :: arr(n), arr_dupe(n)
    logical(kind=c_bool) :: sorted

    call random_uniform(arr, -100_i32, 100_i32)
    arr_dupe = arr
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = isort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

    arr = qsort(arr_dupe)
    sorted = is_sorted(arr)
    write(unit=stdout, fmt='(a,l1,a,*(i4," "))') 'arr, sorted=',sorted,' :: ',arr
    arr = arr_dupe

end program test_sorting
