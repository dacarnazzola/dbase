module pf
use, non_intrinsic :: kinds, only: stdout, dp
implicit none
private

    public :: stdout, dp

end module pf


program ex_particle_filter
use, non_intrinsic :: pf
implicit none

    real(dp) :: truth(6), update(6,6), t, dt

    truth = [607600.0_dp, 120000.0_dp, -4.0_dp*1125.0_dp, 0.0_dp, 0.0_dp, -32.2_dp]
    dt = 1.0_dp
    update = transpose(reshape([1.0_dp, 0.0_dp,     dt, 0.0_dp, 0.5_dp*dt**2,       0.0_dp, &
                                0.0_dp, 1.0_dp, 0.0_dp,     dt,       0.0_dp, 0.5_dp*dt**2, &
                                0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp,           dt,       0.0_dp, &
                                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp,       0.0_dp,           dt, &
                                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,       1.0_dp,       0.0_dp, &
                                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,       0.0_dp,       1.0_dp], [6, 6]))
    t = 0.0_dp
    do
        truth = matmul(update, truth)
        if (truth(2) > 0.0_dp) then
            t = t + dt
            write(stdout,'(7(a,f0.1))') 't: ',t, &
                                        'x: ',truth(1),', y: ',truth(2), &
                                        ', vx: ',truth(3),', vy: ',truth(4), &
                                        ', ax: ',truth(5),', ay: ',truth(6)
        else
            exit
        end if
    end do

end program ex_particle_filter
