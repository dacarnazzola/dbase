module pf
use, non_intrinsic :: kinds, only: stdout, i64, dp, c_bool
use, non_intrinsic :: constants, only: twopi_dp, deg2rad_dp
use, non_intrinsic :: random, only: random_uniform, random_normal
use, non_intrinsic :: vector_math, only: vmag2, vunit2, vdot
implicit none

    real(kind=dp), parameter :: min_range = 1.0_dp, max_range = 300.0_dp*6076.0_dp, &
                                default_range = max_range, default_range_sig = default_range/3.0_dp, &
                                default_range_rate = -1000.0_dp, default_range_rate_sig = 1000.0_dp, &
                                min_accel = -9.0_dp*32.2_dp, max_accel = 9.0_dp*32.2_dp

contains

    pure subroutine gen_meas(obs_state, tgt_state, obs_dim_mask, meas)
        real(kind=dp), intent(in) :: obs_state(6), tgt_state(6)
        logical(kind=c_bool), intent(in) :: obs_dim_mask(3)
        real(kind=dp), intent(out) :: meas(3)
        real(kind=dp) :: dx(2), dvx(2), dxhat(2)
        dx = tgt_state(1:2) - obs_state(1:2)
        if (obs_dim_mask(1)) then !! measurement includes range
            meas(1) = vmag2(dx)
        else
            meas(1) = -1.0_dp
        end if
        if (obs_dim_mask(2)) then !! measurement includes angle
            meas(2) = atan2(dx(2), dx(1))
        else
            meas(2) = -1.0_dp
        end if
        dvx = tgt_state(3:4) - obs_state(3:4)
        call vunit2(dx, dxhat)
        if (obs_dim_mask(3)) then !! measurement includes range-rate
            meas(3) = vdot(dvx, dxhat)
        else
            meas(3) = -1.0_dp
        end if
    end subroutine gen_meas

    impure subroutine init_particles(obs_state, obs_dim_mask, meas, meas_sig, num_particles, particles)
        real(kind=dp), intent(in) :: obs_state(6)
        logical(kind=c_bool), intent(in) :: obs_dim_mask(3)
        real(kind=dp), intent(in) :: meas(3), meas_sig(3)
        integer(kind=i64), intent(in) :: num_particles
        real(kind=dp), intent(out) :: particles(num_particles,6)
        real(kind=dp) :: rand_range(num_particles), rand_angle(num_particles), rand_range_rate(num_particles), &
                         los_obs2tgt(2)
        integer(kind=i64) :: i
        if (obs_dim_mask(1)) then !! measurement includes range
            call random_normal(rand_range, meas(1), meas_sig(1))
        else
            call random_normal(rand_range, default_range, default_range_sig)
        end if
        if (obs_dim_mask(2)) then !! measurement includes angle
            call random_normal(rand_angle, meas(2), meas_sig(2))
        else
            error stop 'lack of angle measurement not supported'
        end if
        rand_angle = mod(rand_angle, twopi_dp)
        if (obs_dim_mask(3)) then !! measurement includes range-rate
            call random_normal(rand_range_rate, meas(3), meas_sig(3))
        else
            call random_normal(rand_range_rate, default_range_rate, default_range_rate_sig)
        end if
        particles(:,1) = rand_range*cos(rand_angle)
        particles(:,2) = rand_range*sin(rand_angle)
        do i=1_i64,num_particles
            los_obs2tgt = particles(i,1:2) - obs_state(1:2)
            call vunit2(los_obs2tgt)
            particles(i,3:4) = rand_range_rate(i)*los_obs2tgt
        end do
        particles(:,5:6) = 0.0_dp !! initialize with 0.0 acceleration
    end subroutine init_particles

end module pf


program ex_particle_filter
use, non_intrinsic :: pf
implicit none

    integer(i64), parameter :: num_particles = 10000_i64
    real(dp), parameter :: obs_state(6) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
                           meas_sig(3) = [6076.0_dp, 1.0_dp*deg2rad_dp, 100.0_dp]
    logical(c_bool), parameter :: obs_dim_mask(3) = [logical(.true., kind=c_bool), &
                                                     logical(.true., kind=c_bool), &
                                                     logical(.true., kind=c_bool)]

    real(dp) :: truth(6), update(6,6), t, dt, particles(num_particles,6), meas(3)

    truth = [607600.0_dp, 120000.0_dp, -4.0_dp*1125.0_dp, 0.0_dp, 0.0_dp, -32.2_dp]
    dt = 1.0_dp
    update = transpose(reshape([1.0_dp, 0.0_dp,     dt, 0.0_dp, 0.5_dp*dt**2,       0.0_dp, &
                                0.0_dp, 1.0_dp, 0.0_dp,     dt,       0.0_dp, 0.5_dp*dt**2, &
                                0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp,           dt,       0.0_dp, &
                                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp,       0.0_dp,           dt, &
                                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,       1.0_dp,       0.0_dp, &
                                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,       0.0_dp,       1.0_dp], [6, 6]))
    t = 0.0_dp
    write(stdout,'(7(a,f0.1))') 't: ',t, &
                                ', x: ',truth(1),', y: ',truth(2), &
                                ', vx: ',truth(3),', vy: ',truth(4), &
                                ', ax: ',truth(5),', ay: ',truth(6)

    call gen_meas(obs_state=obs_state, tgt_state=truth, obs_dim_mask=obs_dim_mask, meas=meas)
    write(*,*) 'TRUTH meas R:',meas(1),', ang: ',meas(2),', R-rt: ',meas(3)
    call init_particles(obs_state=obs_state, &
                        obs_dim_mask=obs_dim_mask, &
                        meas=meas, &
                        meas_sig=meas_sig, &
                        num_particles=num_particles, &
                        particles=particles)
!    do
!        truth = matmul(update, truth)
!        if (truth(2) > 0.0_dp) then
!            t = t + dt
!            write(stdout,'(7(a,f0.1))') 't: ',t, &
!                                        ', x: ',truth(1),', y: ',truth(2), &
!                                        ', vx: ',truth(3),', vy: ',truth(4), &
!                                        ', ax: ',truth(5),', ay: ',truth(6)
!        else
!            exit
!        end if
!    end do


end program ex_particle_filter
