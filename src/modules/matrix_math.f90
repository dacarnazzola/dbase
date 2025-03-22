module matrix_math
use, non_intrinsic :: kinds, only: i64, sp, dp, c_bool
use, non_intrinsic :: constants, only: eps_sp, eps_dp
use, non_intrinsic :: system, only: debug_error_condition
implicit none
private

    interface chol
        module procedure :: chol_sp
        module procedure :: chol_dp
    end interface chol

    interface qr
    end interface qr

    public :: chol, qr

contains

    pure subroutine chol_sp(A, L)
        real(kind=sp), intent(in) :: A(:,:)
        real(kind=sp), intent(out) :: L(size(A, dim=1, kind=i64),size(A, dim=2, kind=i64))
        integer(kind=i64) :: n, i, j, k
        real(kind=sp) :: sum_of_squares, tol
        call debug_error_condition(logical(size(A, dim=1, kind=i64) /= size(A, dim=2, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine chol requires A is a square (n x n) matrix')
        tol = max(eps_sp, eps_sp*maxval(abs(A)))
        L = 0.0_sp
        n = size(A, dim=1, kind=i64)
        do i=1_i64,n
            sum_of_squares = A(i,i)
            do j=1_i64,i-1_i64
                sum_of_squares = sum_of_squares - L(i,j)**2
            end do
            if (sum_of_squares > tol) then
                L(i,i) = sqrt(sum_of_squares)
                do j=i+1_i64,n
                    sum_of_squares = A(j,i)
                    do k=1_i64,i-1_i64
                        sum_of_squares = sum_of_squares - L(j,k)*L(i,k)
                    end do
                    L(j,i) = sum_of_squares/L(i,i)
                end do
            else if (sum_of_squares > -tol) then
                L(:,i) = 0.0_sp
            else
                call debug_error_condition(logical(.true., kind=c_bool), &
                                           'module MATRIX_MATH :: subroutine chol requires A be positive semi-definite')
            end if
        end do
    end subroutine chol_sp

    pure subroutine chol_dp(A, L)
        real(kind=dp), intent(in) :: A(:,:)
        real(kind=dp), intent(out) :: L(size(A, dim=1, kind=i64),size(A, dim=2, kind=i64))
        integer(kind=i64) :: n, i, j, k
        real(kind=dp) :: sum_of_squares, tol
        call debug_error_condition(logical(size(A, dim=1, kind=i64) /= size(A, dim=2, kind=i64), kind=c_bool), &
                                   'module MATRIX_MATH :: subroutine chol requires A is a square (n x n) matrix')
        tol = max(eps_dp, eps_dp*maxval(abs(A)))
        L = 0.0_dp
        n = size(A, dim=1, kind=i64)
        do i=1_i64,n
            sum_of_squares = A(i,i)
            do j=1_i64,i-1_i64
                sum_of_squares = sum_of_squares - L(i,j)**2
            end do
            if (sum_of_squares > tol) then
                L(i,i) = sqrt(sum_of_squares)
                do j=i+1_i64,n
                    sum_of_squares = A(j,i)
                    do k=1_i64,i-1_i64
                        sum_of_squares = sum_of_squares - L(j,k)*L(i,k)
                    end do
                    L(j,i) = sum_of_squares/L(i,i)
                end do
            else if (sum_of_squares > -tol) then
                L(:,i) = 0.0_dp
            else
                call debug_error_condition(logical(.true., kind=c_bool), &
                                           'module MATRIX_MATH :: subroutine chol requires A be positive semi-definite')
            end if
        end do
    end subroutine chol_dp

    pure subroutine qr(A, tau)
        real(kind=dp), intent(inout) :: A(:,:)
        real(kind=dp), intent(out) :: tau(min(size(A, dim=1, kind=i64), size(A, dim=2, kind=i64)))
        real(kind=dp) :: x(size(A, dim=1, kind=i64)), s, alpha, u(size(A, dim=1, kind=i64))
    end subroutine qr

end module matrix_math
