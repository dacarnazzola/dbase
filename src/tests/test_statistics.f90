program test_statistics
use, non_intrinsic :: kinds, only: stdout, i32, i64, sp, dp
use, non_intrinsic :: statistics, only: dsum, avg, std, normal_pdf
use, non_intrinsic :: random, only: random_normal
implicit none

    integer(kind=i32), parameter :: n = 2_i32**20_i32, val = 1_i32, n2 = 10_i32*n

    integer(kind=i32) :: xi32(n), i
    integer(kind=i64) :: xi64(n)
    real(kind=sp) :: xsp(n), xsp2(n2), isp, pdfi, pdf_sum
    real(kind=dp) :: xdp(n), xdp2(n2)

    write(unit=stdout, fmt='(a,i0,a,i0)') 'TEST_STATISTICS :: ',n,' element arrays initialized to ',val

    xi32 = int(val, kind=i32)
    write(unit=stdout, fmt='(a,i0)') 'TEST_STATISTICS :: xi32 sum ',dsum(xi32)

    xi64 = int(val, kind=i64)
    write(unit=stdout, fmt='(a,i0)') 'TEST_STATISTICS :: xi64 sum ',dsum(xi64)

    xsp = real(val, kind=sp)
    write(unit=stdout, fmt='(a,3(" ",e13.6))') 'TEST_STATISTICS :: xsp sum/avg/std',dsum(xsp),avg(xsp),std(xsp)

    xdp = real(val, kind=dp)
    write(unit=stdout, fmt='(a,3(" ",e22.15))') 'TEST_STATISTICS :: xdp sum/avg/std',dsum(xdp),avg(xdp),std(xdp)
    
    write(unit=stdout, fmt='(a,i0,a)') 'TEST_STATISTICS :: ',n2,' element arrays initialized with random_normal'

    call random_normal(xsp2, 19.93_sp, 8.31_sp)
    write(unit=stdout, fmt='(a," ",e13.6,2(" ",f0.2))') 'TEST_STATISTICS :: xsp2 sum/avg/std',dsum(xsp2),avg(xsp2),std(xsp2)

    call random_normal(xdp2, 19.94_dp, 5.30_dp)
    write(unit=stdout, fmt='(a," ",e22.15,2(" ",f0.2))') 'TEST_STATISTICS :: xdp2 sum/avg/std',dsum(xdp2),avg(xdp2),std(xdp2)

    pdf_sum = 0.0_sp
    do i=-6_i32,6_i32
        isp = real(i, kind=sp)
        pdfi = normal_pdf(isp, 0.0_sp, 1.0_sp)
        pdf_sum = pdf_sum + pdfi
        write(unit=stdout, fmt='(a,f0.1,a,f0.6)') 'TEST_STATISTICS :: normal_pdf(',isp,', 0.0_sp, 1.0_sp): ',pdfi
    end do
    write(unit=stdout, fmt='(a,f0.6)') 'TEST STATISTICS :: sum of normal_pdf (expect 1.0): ',pdf_sum

end program test_statistics
