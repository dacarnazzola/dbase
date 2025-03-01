module matrix_math
use, non_intrinsic :: kinds, only: i64, sp, dp, c_bool
use, non_intrinsic :: system, only: debug_error_condition
implicit none
private

    interface chol
        module procedure :: chol_sp
        module procedure :: chol_dp
    end interface chol

    public :: chol

contains

    pure subroutine chol_sp(A, L)
        real(kind=sp), intent(in) :: A(:,:)
        real(kind=sp), intent(out) :: L(size(A, dim=1, kind=i64),size(A, dim=2, kind=i64))
        integer(kind=i64) :: n, i, j, k
        real(kind=sp) :: sum_of_squares
        call debug_error_condition(logical(size(A, dim=1, kind=i64) /= size(A, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine chol requires A is a square (n x n) matrix')
        L = 0.0_sp
        n = size(A, dim=1, kind=i64)
        do i=1_i64,n
            sum_of_squares = A(i,i)
            do j=1_i64,i-1_i64
                sum_of_squares = sum_of_squares - L(i,j)**2
            end do
            L(i,i) = sqrt(sum_of_squares)
            do j=i+1_i64,n
                sum_of_squares = A(j,i)
                do k=1_i64,i-1_i64
                    sum_of_squares = sum_of_squares - L(j,k)*L(i,k)
                end do
                L(j,i) = sum_of_squares/L(i,i)
            end do
        end do
    end subroutine chol_sp

    pure subroutine chol_dp(A, L)
        real(kind=dp), intent(in) :: A(:,:)
        real(kind=dp), intent(out) :: L(size(A, dim=1, kind=i64),size(A, dim=2, kind=i64))
        integer(kind=i64) :: n, i, j, k
        real(kind=dp) :: sum_of_squares
        call debug_error_condition(logical(size(A, dim=1, kind=i64) /= size(A, dim=1, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine chol requires A is a square (n x n) matrix')
        L = 0.0_dp
        n = size(A, dim=1, kind=i64)
        do i=1_i64,n
            sum_of_squares = A(i,i)
            do j=1_i64,i-1_i64
                sum_of_squares = sum_of_squares - L(i,j)**2
            end do
            L(i,i) = sqrt(sum_of_squares)
            do j=i+1_i64,n
                sum_of_squares = A(j,i)
                do k=1_i64,i-1_i64
                    sum_of_squares = sum_of_squares - L(j,k)*L(i,k)
                end do
                L(j,i) = sum_of_squares/L(i,i)
            end do
        end do
    end subroutine chol_dp

end module matrix_math
