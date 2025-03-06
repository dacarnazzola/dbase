module pf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: pi_dp, twopi_dp, deg2rad_dp, rad2deg_dp, eps_dp
use, non_intrinsic :: random, only: random_normal, random_uniform
use, non_intrinsic :: statistics, only: dsum, avg, std, cumsum, cov, normalize
use, non_intrinsic :: vector_math, only: vdot
use, non_intrinsic :: system, only: debug_error_condition, nearly
implicit none
private

    logical, parameter :: debug = .false.
    
    real(dp), parameter :: global_minimum_range = 1.0_dp
    real(dp), parameter :: nmi2ft = 1852.0_dp*100.0_dp/2.54_dp/12.0_dp
    real(dp), parameter :: ft2nmi = 12.0_dp*2.54_dp/100.0_dp/1852.0_dp
    real(dp), parameter :: g = 32.2_dp

    interface
        impure subroutine resample_cartesian_particles(n, cartesian_particles, weights)
        import dp
        implicit none
            integer, intent(in) :: n
            real(dp), intent(inout) :: cartesian_particles(:,:), weights(:)
        end
    end interface

    public :: debug, dp, nmi2ft, ft2nmi, deg2rad_dp, rad2deg_dp, g, pi_dp, twopi_dp, &
              perfect_cart2pol, generate_measurements, initialize_particles, convert_particles_cart2pol, &
              generate_state_estimate, calculate_weights, apply_roughening, &
              resample_cartesian_particles, resample_systematic, resample_multinomial, resample_residual, resample_stratified, &
              fly_constant_acceleration, propagate_particles, apply_constraints, &
              rmse, neff, &
              avg, std, nearly

contains

    pure subroutine perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
        real(dp), intent(in) :: obs_cart(:), tgt_cart(size(obs_cart))
        real(dp), intent(out) :: tgt_pol(3)
        call cart2pol_inner(obs_cart, tgt_cart, 0.0_dp, 0.0_dp, 0.0_dp, tgt_pol)
    end

    pure subroutine cart2pol_inner(obs_cart, tgt_cart, err_rng, err_ang, err_rngrt, tgt_pol)
        real(dp), intent(in) :: obs_cart(:), tgt_cart(size(obs_cart)), err_rng, err_ang, err_rngrt
        real(dp), intent(out) :: tgt_pol(3)
        real(dp) :: rng, loshat(2), losv(2)
        rng = sqrt((tgt_cart(1) - obs_cart(1))**2 + (tgt_cart(2) - obs_cart(2))**2)
        tgt_pol(1) = max(global_minimum_range, rng + err_rng)
        loshat = (tgt_cart(1:2) - obs_cart(1:2))/rng
        tgt_pol(2) = mod(atan2(loshat(2), loshat(1)) + err_ang + pi_dp, twopi_dp) - pi_dp
        losv = tgt_cart(3:4) - obs_cart(3:4)
        tgt_pol(3) = losv(1)*loshat(1) + losv(2)*loshat(2) + err_rngrt
    end

    impure subroutine generate_measurements(obs_cart, tgt_cart, meas_sig, no_meas_max_rng, no_meas_max_spd, n, polar_measurements)
        real(dp), intent(in) :: obs_cart(:), tgt_cart(size(obs_cart)), meas_sig(3), no_meas_max_rng, no_meas_max_spd
        integer, intent(in) :: n
        real(dp), intent(out) :: polar_measurements(3,n)
        real(dp) :: err_rng(n), err_ang(n), err_rngrt(n), true_rng, true_ang, losv(2), loshat(2), true_rngrt
        integer :: i
        if (meas_sig(1) > 0) then
            call random_normal(err_rng, 0.0_dp, meas_sig(1))
        else
            true_rng = max(1.0e-6_dp, sqrt((tgt_cart(1) - obs_cart(1))**2 + (tgt_cart(2) - obs_cart(2))**2))
            call random_uniform(err_rng, global_minimum_range - true_rng, no_meas_max_rng - true_rng)
        end if
        if (meas_sig(2) > 0) then
            call random_normal(err_ang, 0.0_dp, meas_sig(2))
        else
            true_ang = atan2(tgt_cart(2) - obs_cart(2), tgt_cart(1) - obs_cart(1))
            call random_uniform(err_ang, -pi_dp - true_ang, pi_dp - true_ang)
        end if
        if (meas_sig(3) > 0) then
            call random_normal(err_rngrt, 0.0_dp, meas_sig(3))
        else
            losv = tgt_cart(3:4) - obs_cart(3:4)
            loshat = (tgt_cart(1:2) - obs_cart(1:2)) / &
                     max(1.0e-6_dp, sqrt((tgt_cart(1) - obs_cart(1))**2 + (tgt_cart(2) - obs_cart(2))**2))
            true_rngrt = losv(1)*loshat(1) + losv(2)*loshat(2)
            call random_uniform(err_rngrt, -no_meas_max_spd - true_rngrt, no_meas_max_spd - true_rngrt)
        end if
        do concurrent (i=1:n)
            call cart2pol_inner(obs_cart, tgt_cart, err_rng(i), err_ang(i), err_rngrt(i), polar_measurements(:,i))
        end do
    end

    impure subroutine resample_measurements(meas, meas_sig, no_meas_max_rng, no_meas_max_spd, n, polar_measurements)
        real(dp), intent(in) :: meas(3), meas_sig(3), no_meas_max_rng, no_meas_max_spd
        integer, intent(in) :: n
        real(dp), intent(out) :: polar_measurements(3,n)
        if (meas_sig(1) > 0) then
            call random_normal(polar_measurements(1,:), meas(1), meas_sig(1))
        else
            call random_uniform(polar_measurements(1,:), global_minimum_range, no_meas_max_rng)
        end if
        if (meas_sig(2) > 0) then
            call random_normal(polar_measurements(2,:), meas(2), meas_sig(2))
        else
            call random_uniform(polar_measurements(2,:), -pi_dp, pi_dp)
        end if
        if (meas_sig(3) > 0) then
            call random_normal(polar_measurements(3,:), meas(3), meas_sig(3))
        else
            call random_uniform(polar_measurements(3,:), -no_meas_max_spd, no_meas_max_spd)
        end if
    end

    impure subroutine initialize_particles(obs_cart, meas, meas_sig, tgt_max_rng, tgt_max_spd, n, cartesian_particles, weights)
        real(dp), intent(in) :: obs_cart(:), meas(3), meas_sig(3), tgt_max_rng, tgt_max_spd
        integer, intent(in) :: n
        real(dp), intent(out) :: cartesian_particles(size(obs_cart),n), weights(n)
        real(dp) :: polar_measurements(3,n), tgt_spd_scale(n), ax(n), ay(n)
        integer :: i
        call resample_measurements(meas, meas_sig, tgt_max_rng, tgt_max_spd, n, polar_measurements)
        call random_uniform(tgt_spd_scale, 0.0_dp, 1.0_dp)
        !! assume 6 cartesian states correspond to [x, y, vx, vy, ax, ay]
        if (size(obs_cart) == 6) then
            !! sample acceleration ax and ay from N(0,1)
            call random_normal(ax, 0.0_dp, 1.0_dp*g)
            call random_normal(ay, 0.0_dp, 1.0_dp*g)
        end if
        do concurrent (i=1:n)
            call pol2cart_inner(obs_cart, tgt_max_rng, tgt_max_spd, tgt_spd_scale(i), &
                                polar_measurements(:,i), cartesian_particles(:,i))
            if (size(obs_cart) == 6) then
                cartesian_particles(5,i) = ax(i)
                cartesian_particles(6,i) = ay(i)
            end if
        end do
        weights = 1.0_dp/real(n,dp)
    end

    pure subroutine pol2cart_inner(obs_cart, tgt_max_rng, tgt_max_spd, tgt_spd_scale, tgt_pol, tgt_cart)
        real(dp), intent(in) :: obs_cart(:), tgt_max_rng, tgt_max_spd, tgt_spd_scale, tgt_pol(3)
        real(dp), intent(out) :: tgt_cart(size(obs_cart))
        real(dp) :: cos_ang, sin_ang, tgt_rng, tgt_min_spd, tgt_spd, tgt_vt_mag, tgt_final_spd_scale
        integer :: tgt_vt_dir
        cos_ang = cos(tgt_pol(2))
        sin_ang = sin(tgt_pol(2))
        tgt_rng = min(tgt_max_rng, tgt_pol(1))
        tgt_cart(1) = obs_cart(1) + tgt_rng*cos_ang
        tgt_cart(2) = obs_cart(2) + tgt_rng*sin_ang
        tgt_cart(3) = obs_cart(3) + tgt_pol(3)*cos_ang
        tgt_cart(4) = obs_cart(4) + tgt_pol(3)*sin_ang
        tgt_min_spd = min(tgt_max_spd, sqrt(tgt_cart(3)**2 + tgt_cart(4)**2))
        tgt_spd = tgt_spd_scale*(tgt_max_spd - tgt_min_spd) + tgt_min_spd
        tgt_vt_mag = sqrt(max(0.0_dp, tgt_spd**2 - tgt_min_spd**2))
        tgt_vt_dir = (-1)**mod(floor(tgt_spd_scale*100000000.0_dp), 2)
        tgt_cart(3) = tgt_cart(3) - tgt_vt_dir*tgt_vt_mag*sin_ang
        tgt_cart(4) = tgt_cart(4) + tgt_vt_dir*tgt_vt_mag*cos_ang
        tgt_final_spd_scale = max(1.0_dp, sqrt(tgt_cart(3)**2 + tgt_cart(4)**2)/tgt_spd)
        tgt_cart(3:4) = tgt_cart(3:4)/tgt_final_spd_scale
    end

    pure subroutine convert_particles_cart2pol(obs_cart, n, cartesian_particles, polar_particles)
        real(dp), intent(in) :: obs_cart(:) 
        integer, intent(in) :: n
        real(dp), intent(in) :: cartesian_particles(size(obs_cart),n)
        real(dp), intent(out) :: polar_particles(3,n)
        integer :: i
        do concurrent (i=1:n)
            call perfect_cart2pol(obs_cart, cartesian_particles(:,i), polar_particles(:,i))
        end do
    end

    pure subroutine generate_state_estimate(state_estimate, n, particles, weights)
        real(dp), intent(out) :: state_estimate(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: particles(size(state_estimate),n), weights(n)
        integer :: i
        do concurrent (i=1:size(state_estimate))
            state_estimate(i) = vdot(particles(i,:), weights)
        end do
    end

    pure subroutine calculate_weights(tgt_pol, meas_sig, n, polar_particles, weights)
        real(dp), intent(in) :: tgt_pol(3), meas_sig(3)
        integer, intent(in) :: n
        real(dp), intent(in) :: polar_particles(3,n)
        real(dp), intent(inout) :: weights(n)
        real(dp) :: err_fac(3)
        integer :: i
        do concurrent (i=1:n)
            err_fac = 0.0_dp
            if (meas_sig(1) > 0) err_fac(1) = (polar_particles(1,i) - tgt_pol(1))**2/meas_sig(1)**2
            if (meas_sig(2) > 0) err_fac(2) = (mod(polar_particles(2,i) - tgt_pol(2) + pi_dp, twopi_dp) - pi_dp)**2/meas_sig(2)**2
            if (meas_sig(3) > 0) err_fac(3) = (polar_particles(3,i) - tgt_pol(3))**2/meas_sig(3)**2
            weights(i) = weights(i)*exp(-0.5_dp*(err_fac(1) + err_fac(2) + err_fac(3))) + eps_dp
        end do
        weights = weights/dsum(weights)
    end

    impure subroutine resample_systematic(n, cartesian_particles, weights)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:)
        real(dp) :: u(1), inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2))
        integer :: i, j
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        call cumsum(weights)
        weights(n) = 1.0_dp
        inv_n = 1.0_dp/real(n, kind=dp)
        call random_uniform(u, 0.0_dp, inv_n)
        j = 1
        do i=1,n
            do while (u(1) > weights(j))
                j = j + 1
            end do
            new_particles(:,i) = cartesian_particles(:,j)
            u = u + inv_n
        end do
        cartesian_particles = new_particles
        weights = inv_n
    end

    impure subroutine resample_multinomial(n, cartesian_particles, weights)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:)
        real(dp) :: u(n), inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2))
        integer :: i, j
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        call cumsum(weights)
        weights(n) = 1.0_dp
        call random_uniform(u, 0.0_dp, 1.0_dp)
        do i=1,n
            j = 1
            do while (u(i) > weights(j))
                j = j + 1
            end do
            new_particles(:,i) = cartesian_particles(:,j)
        end do
        cartesian_particles = new_particles
        inv_n = 1.0_dp/real(n, kind=dp)
        weights = inv_n
    end

    impure subroutine resample_residual(n, cartesian_particles, weights)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:)
        real(dp) :: inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2)), resamples(n), u(n)
        integer :: i, resamples_int, j, k
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        resamples = weights*n
        j = 1
        do i=1,n
            resamples_int = floor(resamples(i))
            do k=1,resamples_int
                new_particles(:,j) = cartesian_particles(:,i)
                j = j + 1
            end do
        end do
        resamples = resamples - real(floor(resamples), kind=dp)
        if (dsum(resamples) > 0.0_dp) then
            call normalize(resamples)
            call cumsum(resamples)
            resamples(n) = 1.0_dp
            call random_uniform(u(1:n-j+1), 0.0_dp, 1.0_dp)
            do i=1,n-j+1
                k = 1
                do while (u(i) > resamples(k))
                    k = k + 1
                end do
                new_particles(:,i+j-1) = cartesian_particles(:,k)
            end do
        end if
        cartesian_particles = new_particles
        inv_n = 1.0_dp/real(n, kind=dp)
        weights = inv_n
    end

    impure subroutine resample_stratified(n, cartesian_particles, weights)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:)
        real(dp) :: u(n), inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2))
        integer :: i, j
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        inv_n = 1.0_dp/real(n, kind=dp)
        call random_uniform(u, 0.0_dp, 1.0_dp)
        do i=1,n
            u(i) = (u(i) + real(i-1,dp))*inv_n
        end do
        call cumsum(weights)
        weights(n) = 1.0_dp
        j = 1
        do i=1,n
            do while (u(i) > weights(j))
                j = j + 1
            end do
            new_particles(:,i) = cartesian_particles(:,j)
        end do
        cartesian_particles = new_particles
        weights = inv_n
    end

    impure subroutine apply_roughening(z_scale, n, particles)
        real(dp), intent(in) :: z_scale
        integer, intent(in) :: n
        real(dp), intent(inout) :: particles(:,:)
        real(dp) :: particles_cov(size(particles,dim=1),size(particles,dim=1)), &
                    mu_zero(size(particles,dim=1)), z(size(particles,dim=1),size(particles,dim=2))
        integer :: i
        call debug_error_condition(size(particles, dim=2) /= n, 'mismatch in particles shape')
        call cov(particles, particles_cov)
        mu_zero = 0.0_dp
        call random_normal(z, mu_zero, particles_cov)
        do concurrent (i=1:n)
            particles(:,i) = particles(:,i) + z(:,i)*z_scale
        end do
    end

    pure subroutine fly_constant_acceleration(cart6_state, dt, max_spd)
        real(dp), intent(inout) :: cart6_state(:)
        real(dp), intent(in) :: dt, max_spd
        real(dp) :: spd_scale
        call debug_error_condition(size(cart6_state) /= 6, 'fly_constant_acceleration assumes state [x, y, vx, vy, ax, ay]')
        cart6_state(3:4) = cart6_state(3:4) + dt*cart6_state(5:6)
        spd_scale = max(1.0_dp, sqrt(cart6_state(3)**2 + cart6_state(4)**2)/max_spd)
        cart6_state(3:4) = cart6_state(3:4)/spd_scale
        cart6_state(1:2) = cart6_state(1:2) + dt*cart6_state(3:4)
    end

    impure subroutine propagate_particles(dt, jerk_sig, max_spd, max_acc, n, cartesian_particles)
        real(dp), intent(in) :: dt, jerk_sig, max_spd, max_acc
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:)
        real(dp) :: jx(n), jy(n), acc_scale
        integer :: i
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call random_normal(jx, 0.0_dp, jerk_sig)
        call random_normal(jy, 0.0_dp, jerk_sig)
        do concurrent (i=1:n)
            cartesian_particles(5,i) = cartesian_particles(5,i) + jx(i)
            cartesian_particles(6,i) = cartesian_particles(6,i) + jy(i)
            acc_scale = max(1.0_dp, sqrt(cartesian_particles(5,i)**2 + cartesian_particles(6,i)**2)/max_acc)
            cartesian_particles(5:6,i) = cartesian_particles(5:6,i)/acc_scale
            call fly_constant_acceleration(cartesian_particles(:,i), dt, max_spd)
        end do
    end

    pure subroutine apply_constraints(obs_cart, max_rng, max_spd, max_acc, n, cartesian_particles, polar_particles)
        real(dp), intent(in) :: obs_cart(:), max_rng, max_spd, max_acc
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:)
        real(dp), intent(out) :: polar_particles(:,:)
        real(dp) :: los(2), ang, spd_scale, acc_scale
        integer :: i
        call debug_error_condition(size(cartesian_particles, dim=1) /= size(obs_cart), 'mismatch in obs_cart vs particle dims')
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(polar_particles, dim=2) /= n, 'mismatch in polar_particles shape')
        do concurrent (i=1:n)
            los = cartesian_particles(1:2,i) - obs_cart(1:2)
            if ((los(1)**2 + los(2)**2) > max_rng**2) then
                ang = atan2(los(2), los(1))
                cartesian_particles(1,i) = obs_cart(1) + max_rng*cos(ang)
                cartesian_particles(2,i) = obs_cart(2) + max_rng*sin(ang)
            end if
            spd_scale = max(1.0_dp, sqrt(cartesian_particles(3,i)**2 + cartesian_particles(4,i)**2)/max_spd)
            cartesian_particles(3:4,i) = cartesian_particles(3:4,i)/spd_scale
            acc_scale = max(1.0_dp, sqrt(cartesian_particles(5,i)**2 + cartesian_particles(6,i)**2)/max_acc)
            cartesian_particles(5:6,i) = cartesian_particles(5:6,i)/acc_scale
        end do
        call convert_particles_cart2pol(obs_cart, n, cartesian_particles, polar_particles)
    end

    pure function rmse(predicted, observed) result(val)
        real(dp), intent(in) :: predicted(:), observed
        real(dp) :: val
        real(dp) :: diff2(size(predicted))
        diff2 = (predicted - observed)**2
        val = sqrt(avg(diff2))
    end

    pure function neff(weights) result(val)
        real(dp), intent(in) :: weights(:)
        real(dp) :: val
        real(dp) :: weights2(size(weights))
        weights2 = weights**2
        val = 1.0_dp/dsum(weights2)
    end

end module pf


program ex_particle_filter
use, non_intrinsic :: pf
implicit none

    integer, parameter :: max_trials = 128
    real(dp), parameter :: max_rng   = 500.0_dp*nmi2ft, &
                           max_spd   = 10000.0_dp, &
                           max_acc   = 9.0_dp*g, &
                           dt        = 1.0_dp, &
                           meas_sig(3) = [10.0_dp*nmi2ft, 5.0_dp*deg2rad_dp, 200.0_dp] ! poor measurements
!                           meas_sig(3) = [100.0_dp, 0.1_dp*deg2rad_dp, 10.0_dp] ! standard measurements
!                           meas_sig(3) = [1.0_dp, 0.001_dp, 1.0_dp] ! exquisite measurements
!                           meas_sig(3) = [-1.0_dp, 0.1_dp*deg2rad_dp, -1.0_dp] ! bearing only

    procedure(resample_cartesian_particles), pointer :: resample_subroutine
    integer :: dr, tmax_int, spd, rough_fac_int, neff_pass_int, num_particles, trial, num_particles_ii, init_meas_sig, &
               resampling_method, jerk_sig_pow
    real(dp) :: obs_cart(6), tgt_cart(6), tgt_pol(3), t, tmax, rough_fac, neff_thresh, neff_pass, neff0, &
                meas(3), est_cart(6), est_pol(3), meas_sig_fac, jerk_sig, &
                x_err(max_trials), y_err(max_trials), pos_err(max_trials), &
                vx_err(max_trials), vy_err(max_trials), vel_err(max_trials), &
                ax_err(max_trials), ay_err(max_trials), acc_err(max_trials), &
                rng_err(max_trials), ang_err(max_trials), rngrt_err(max_trials)
    real(dp), allocatable :: cartesian_particles(:,:), weights(:), polar_particles(:,:)
    character(len=1024) :: fmtstr, resample_strategy

    write(831,'(a)') 'num_particles,init_meas_sig,resampling_method,jerk_sig,rough_fac,neff_thresh,'// &
                     'tmax_sec,downrange_nmi,spd_mach,'// &
                     'rmse_rng,rmse_ang,rmse_rngrt,rmse_x,rmse_y,rmse_vx,rmse_vy,rmse_ax,rmse_ay,rmse_pos,rmse_vel,rmse_acc'
    do num_particles_ii=1,11,1
        num_particles = 10**(num_particles_ii/10)*1000*mod(num_particles_ii,10)
        if (num_particles == 0) cycle
        if (allocated(cartesian_particles)) deallocate(cartesian_particles)
        if (allocated(weights)) deallocate(weights)
        if (allocated(polar_particles)) deallocate(polar_particles)
        allocate(cartesian_particles(size(obs_cart),num_particles), weights(num_particles), polar_particles(3,num_particles))
    do init_meas_sig=1,6 !! multiply meas_sig by init_meas_sig for particle initialization
        meas_sig_fac = real(init_meas_sig, kind=dp)
    do jerk_sig_pow=-1,-4,-1 !! sigma used for random jerk in propagation step                                        
        jerk_sig = 10.0_dp**jerk_sig_pow*g                                                                            
    do neff_pass_int=5,100,5 !! respample when Neff is 5-100% of num_particles                                        
        neff_thresh = real(neff_pass_int, kind=dp)/100.0_dp                                                           
        neff_pass = neff_thresh*num_particles                                                                         
        neff0 = num_particles                                                                                         
    do rough_fac_int=5,100,5 !! scaling Neff**(-1/(d+4))
        rough_fac = real(rough_fac_int, kind=dp)/100.0_dp
    do resampling_method=1,4 !! simple select case on method used
        select case (resampling_method)
            case (1)
                resample_strategy = 'systematic'
                resample_subroutine => resample_systematic
            case (2)
                resample_strategy = 'multinomial'
                resample_subroutine => resample_multinomial
            case (3)
                resample_strategy = 'residual'
                resample_subroutine => resample_residual
            case (4)
                resample_strategy = 'stratified'
                resample_subroutine => resample_stratified
            case default
                error stop 'not implemented yet'
        end select
    do tmax_int=3,10 !! performance should increase with tmax until target passes max_spd
        tmax = real(tmax_int, kind=dp)
        t = 0.0_dp
    do dr=250,250 !! NMI
    do spd=3,3 !! ~Mach
    !$omp parallel do default(firstprivate) shared(x_err, y_err, vx_err, vy_err, ax_err, ay_err, pos_err, vel_err, acc_err, &
    !$omp&                                         rng_err, ang_err, rngrt_err)
    do trial=1,max_trials
        !! reset starting positions for each run
        obs_cart = 0.0_dp
        tgt_cart = 0.0_dp !! initialize all state components to zero (0.0)
        tgt_cart(1) = dr*nmi2ft !! x position [ft]
        tgt_cart(3) = -spd*cos(45*deg2rad_dp) !! vx velocity [ft/sec]
        tgt_cart(4) = spd*sin(45*deg2rad_dp) !! vy velocity [ft/sec]
        tgt_cart(6) = -32.2_dp !! ay acceleration [ft/sec**2]
        call perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
        call generate_measurements(obs_cart, tgt_cart, meas_sig, max_rng, max_spd, 1, meas)
        call initialize_particles(obs_cart, meas, meas_sig_fac*meas_sig, max_rng, max_spd, &
                                  num_particles, cartesian_particles, weights)
        t = 0.0_dp
        do while (t < tmax)
            !! advance time and target independent of anything else
            t = t + dt
            call fly_constant_acceleration(tgt_cart, dt, huge(max_spd))
            call fly_constant_acceleration(obs_cart, dt, huge(max_spd))
            !! update tgt_pol
            call perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
            !! generate new measurement from observer perspective
            call generate_measurements(obs_cart, tgt_cart, meas_sig, max_rng, max_spd, 1, meas)
            !! propagate current particles
            call propagate_particles(dt, jerk_sig, max_spd, max_acc, num_particles, cartesian_particles)
            !! apply constraints and convert particles from cartesian to polar representation
            call apply_constraints(obs_cart, max_rng, max_spd, max_acc, num_particles, cartesian_particles, polar_particles)
            !! calculate weights, resample if Neff is below neff_pass (acceptable percentage of original num_particles)
            call calculate_weights(meas, meas_sig, num_particles, polar_particles, weights)
            neff0 = neff(weights)
            if (neff0 < neff_pass) then
                call resample_subroutine(num_particles, cartesian_particles, weights)
                call apply_roughening(rough_fac*neff0**(-1.0_dp/(size(cartesian_particles,dim=1)+4)), &
                                      num_particles, cartesian_particles)
            end if
        end do
        call generate_state_estimate(est_pol, num_particles, polar_particles, weights)
        rng_err(trial) = tgt_pol(1) - est_pol(1)
        ang_err(trial) = mod(tgt_pol(2) - est_pol(2) + pi_dp, twopi_dp) - pi_dp
        rngrt_err(trial) = tgt_pol(3) - est_pol(3)
        call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
        x_err(trial) = tgt_cart(1) - est_cart(1)
        y_err(trial) = tgt_cart(2) - est_cart(2)
        vx_err(trial) = tgt_cart(3) - est_cart(3)
        vy_err(trial) = tgt_cart(4) - est_cart(4)
        ax_err(trial) = tgt_cart(5) - est_cart(5)
        ay_err(trial) = tgt_cart(6) - est_cart(6)
        pos_err(trial) = sqrt((tgt_cart(1) - est_cart(1))**2 + (tgt_cart(2) - est_cart(2))**2)
        vel_err(trial) = sqrt((tgt_cart(3) - est_cart(3))**2 + (tgt_cart(4) - est_cart(4))**2)
        acc_err(trial) = sqrt((tgt_cart(5) - est_cart(5))**2 + (tgt_cart(6) - est_cart(6))**2)
    end do
    !$omp end parallel do
        !! trials done
        fmtstr = '(2(i0,","),a,",",e13.6,",",2(f0.4,","),f0.1,",",2(i0,","),11(e13.6,","),e13.6)'
        if (all(nearly(pos_err, pos_err))) then
            write(831,fmt=trim(fmtstr)) num_particles,init_meas_sig,trim(resample_strategy),jerk_sig,rough_fac,neff_thresh, &
                                  tmax,dr,spd, &
                                  rmse(rng_err,0.0_dp),rmse(ang_err,0.0_dp),rmse(rngrt_err,0.0_dp), &
                                  rmse(x_err,0.0_dp),rmse(y_err,0.0_dp), &
                                  rmse(vx_err,0.0_dp),rmse(vy_err,0.0_dp), &
                                  rmse(ax_err,0.0_dp),rmse(ay_err,0.0_dp), &
                                  rmse(pos_err,0.0_dp),rmse(vel_err,0.0_dp),rmse(acc_err,0.0_dp)
        else
            error stop 'nan should not happen'
        end if
        fmtstr = '(2(a,i0),a,a11,a,e13.6,2(a,f0.2),a,f0.1,a,e13.6)'
        write(*,fmtstr) 'particles: ',num_particles,', init_meas_sig: ',init_meas_sig, &
                        ', ',trim(resample_strategy), &
                        ', jerk_sig: ',jerk_sig, &
                        ', rough_fac: ',rough_fac,', neff_thresh: ',neff_thresh, &
                        ', tmax: ',tmax, &
                        ', rmse pos_err: ',rmse(pos_err,0.0_dp)
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do

end program ex_particle_filter
