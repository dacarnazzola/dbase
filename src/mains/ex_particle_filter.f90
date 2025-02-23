module pf
use, non_intrinsic :: kinds, only: stdout, i64, dp, c_bool
use, non_intrinsic :: constants, only: twopi_dp, deg2rad_dp
use, non_intrinsic :: random, only: random_uniform, random_normal
use, non_intrinsic :: vector_math, only: vmag, vunit, vdot
use, non_intrinsic :: statistics, only: mv_normal_pdf, normalize
implicit none

    real(kind=dp), parameter :: min_range = 1.0_dp, max_range = 300.0_dp*6076.0_dp, &
                                default_range = max_range, default_range_sig = default_range/3.0_dp, &
                                default_range_rate = -1000.0_dp, default_range_rate_sig = 1000.0_dp, &
                                min_accel = -9.0_dp*32.2_dp, max_accel = 9.0_dp*32.2_dp

contains

    pure subroutine cart2pol(obs_state, tgt_state, obs_dim_mask, meas)
        real(kind=dp), intent(in) :: obs_state(6), tgt_state(6)
        logical(kind=c_bool), intent(in) :: obs_dim_mask(3)
        real(kind=dp), intent(out) :: meas(3)
        real(kind=dp) :: dx(2), dvx(2), dxhat(2)
        dx = tgt_state(1:2) - obs_state(1:2)
        if (obs_dim_mask(1)) then !! measurement includes range
            meas(1) = vmag(dx)
        else
            meas(1) = -1.0_dp
        end if
        if (obs_dim_mask(2)) then !! measurement includes angle
            meas(2) = atan2(dx(2), dx(1))
        else
            meas(2) = -1.0_dp
        end if
        dvx = tgt_state(3:4) - obs_state(3:4)
        call vunit(dx, dxhat)
        if (obs_dim_mask(3)) then !! measurement includes range-rate
            meas(3) = vdot(dvx, dxhat)
        else
            meas(3) = -1.0_dp
        end if
    end subroutine cart2pol

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
            call vunit(los_obs2tgt)
            particles(i,3:4) = rand_range_rate(i)*los_obs2tgt
        end do
        particles(:,5:6) = 0.0_dp !! initialize with 0.0 acceleration
    end subroutine init_particles

    pure subroutine particles_cart2pol(num_particles, obs_state, particles_cart, obs_dim_mask, particles_pol)
        integer(kind=i64), intent(in) :: num_particles
        real(kind=dp), intent(in) :: obs_state(6), particles_cart(num_particles,6)
        logical(kind=c_bool), intent(in) :: obs_dim_mask(3)
        real(kind=dp), intent(out) :: particles_pol(3,num_particles)
        integer(kind=i64) :: i
        do i=1_i64,num_particles
            call cart2pol(obs_state, particles_cart(i,:), obs_dim_mask, particles_pol(:,i))
        end do
    end subroutine particles_cart2pol

    pure subroutine calc_weights(meas, meas_sig, num_particles, particles, weights)
        real(kind=dp), intent(in) :: meas(:), meas_sig(size(meas,kind=i64))
        integer(kind=i64), intent(in) :: num_particles
        real(kind=dp), intent(in) :: particles(size(meas,kind=i64),num_particles)
        real(kind=dp), intent(inout) :: weights(num_particles)
        real(kind=dp) :: pdf_vals(num_particles)
        call mv_normal_pdf(pdf_vals, size(meas, kind=i64), particles, meas, meas_sig)
        weights = weights*pdf_vals + 1.0e-6_dp
        call normalize(weights)
    end subroutine calc_weights

    pure function neff(weights) result(val)
        real(kind=dp), intent(in) :: weights(:)
        real(kind=dp) :: val
        val = 1.0_dp/vdot(weights, weights)
    end function neff

end module pf


program ex_particle_filter
use, non_intrinsic :: pf
implicit none

    integer(i64), parameter :: num_particles = 10000_i64
    real(dp), parameter :: obs_state(6) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
!                           meas_sig(3) = [6076.0_dp, 1.0_dp*deg2rad_dp, 100.0_dp]
                           meas_sig(3) = [10.0_dp, 0.1_dp*deg2rad_dp, 1.0_dp]
    logical(c_bool), parameter :: obs_dim_mask(3) = [logical(.true., kind=c_bool), &
                                                     logical(.true., kind=c_bool), &
                                                     logical(.true., kind=c_bool)]

    real(dp) :: truth(6), update(6,6), t, dt, particles(num_particles,6), meas(3), particles_pol(3,num_particles), &
                weights(num_particles), best_weight, worst_weight
    integer(i64) :: i, best_weight_ii, worst_weight_ii

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

    call cart2pol(obs_state=obs_state, tgt_state=truth, obs_dim_mask=obs_dim_mask, meas=meas)
    write(*,'(a,f0.1,a,f0.6,a,f0.1)') 'TRUTH meas R:',meas(1),', ang: ',meas(2),', R-rt: ',meas(3)
    call init_particles(obs_state=obs_state, &
                        obs_dim_mask=obs_dim_mask, &
                        meas=meas, &
                        meas_sig=meas_sig, &
                        num_particles=num_particles, &
                        particles=particles)
    call particles_cart2pol(num_particles=num_particles, &
                            obs_state=obs_state, &
                            particles_cart=particles, &
                            obs_dim_mask=obs_dim_mask, &
                            particles_pol=particles_pol)
    weights = 1.0_dp/real(num_particles, kind=dp)
    call calc_weights(meas, meas_sig, num_particles, particles_pol, weights)
    if (num_particles < 10_i64) then
        do i=1_i64,num_particles
            write(*,'(a,i0,a,f0.1,a,f0.6,a,f0.1,a,e22.15)') 'particle ',i,' R: ',particles_pol(1,i), &
                                                                          ', ang: ',particles_pol(2,i), &
                                                                          ', R-rt: ',particles_pol(3,i), &
                                                                          ' -- ',weights(i)
        end do
    end if
    write(*,'(a,f0.1)') 'neff: ',neff(weights)
    best_weight = maxval(weights)
    best_weight_ii = findloc(weights, best_weight, dim=1)
    write(*,'(a,f0.6,a,i0)') 'best weight ',best_weight,' at ',best_weight_ii
    write(*,'(a,i0,a,f0.1,a,f0.6,a,f0.1,a,e22.15)') 'particle ',best_weight_ii,' R: ',particles_pol(1,best_weight_ii), &
                                                                  ', ang: ',particles_pol(2,best_weight_ii), &
                                                                  ', R-rt: ',particles_pol(3,best_weight_ii), &
                                                                  ' -- ',weights(best_weight_ii)
    worst_weight = minval(weights)
    worst_weight_ii = findloc(weights, worst_weight, dim=1)
    write(*,'(a,f0.6,a,i0)') 'worst weight ',worst_weight,' at ',worst_weight_ii
    write(*,'(a,i0,a,f0.1,a,f0.6,a,f0.1,a,e22.15)') 'particle ',worst_weight_ii,' R: ',particles_pol(1,worst_weight_ii), &
                                                                  ', ang: ',particles_pol(2,worst_weight_ii), &
                                                                  ', R-rt: ',particles_pol(3,worst_weight_ii), &
                                                                  ' -- ',weights(worst_weight_ii)
    write(*,'(a,e22.15)') 'difference worst - best: ',best_weight - worst_weight

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
