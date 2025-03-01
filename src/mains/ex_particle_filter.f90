module pf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: pi_dp, twopi_dp, deg2rad_dp, rad2deg_dp
use, non_intrinsic :: random, only: random_normal, random_uniform
use, non_intrinsic :: statistics, only: dsum, avg, std, cumsum
use, non_intrinsic :: vector_math, only: vdot
use, non_intrinsic :: system, only: debug_error_condition
implicit none
private

    logical, parameter :: debug = .false.
    
    real(dp), parameter :: global_minimum_range = 1.0_dp
    real(dp), parameter :: nmi2ft = 1852.0_dp*100.0_dp/2.54_dp/12.0_dp
    real(dp), parameter :: ft2nmi = 12.0_dp*2.54_dp/100.0_dp/1852.0_dp
    real(dp), parameter :: g = 32.2_dp

    public :: debug, dp, nmi2ft, ft2nmi, deg2rad_dp, rad2deg_dp, g, &
              perfect_cart2pol, generate_measurements, initialize_particles, convert_particles_cart2pol, &
              generate_state_estimate, calculate_weights, resample_systematic, &
              fly_constant_acceleration, propagate_particles, &
              rmse, neff, print_summary, &
              avg, std

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

    impure subroutine initialize_particles(obs_cart, tgt_cart, meas_sig, tgt_max_rng, tgt_max_spd, n, cartesian_particles, weights)
        real(dp), intent(in) :: obs_cart(:), tgt_cart(size(obs_cart)), meas_sig(3), tgt_max_rng, tgt_max_spd
        integer, intent(in) :: n
        real(dp), intent(out) :: cartesian_particles(size(obs_cart),n), weights(n)
        real(dp) :: polar_measurements(3,n), tgt_spd_scale(n), ax(n), ay(n)
        integer :: i
        call generate_measurements(obs_cart, tgt_cart, meas_sig, tgt_max_rng, tgt_max_spd, n, polar_measurements)
        call random_uniform(tgt_spd_scale, 0.0_dp, 1.0_dp)
        call random_normal(ax, 0.0_dp, 1.0_dp*g)
        call random_normal(ay, 0.0_dp, 1.0_dp*g)
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
            if (meas_sig(2) > 0) err_fac(2) = (polar_particles(2,i) - tgt_pol(2))**2/meas_sig(2)**2
            if (meas_sig(3) > 0) err_fac(3) = (polar_particles(3,i) - tgt_pol(3))**2/meas_sig(3)**2
            weights(i) = weights(i)*exp(-0.5*(err_fac(1) + err_fac(2) + err_fac(3)))
        end do
        weights = weights/dsum(weights)
    end

    impure subroutine resample_systematic(n, cartesian_particles, weights)
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:), weights(:)
        real(dp) :: weights_cdf(n), u(1), inv_n, new_particles(size(cartesian_particles,dim=1),size(cartesian_particles,dim=2))
        integer :: i, j
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call debug_error_condition(size(weights) /= n, 'size of weights array needs to match n')
        call cumsum(weights, weights_cdf)
        inv_n = 1.0_dp/real(n, kind=dp)
        call random_uniform(u, 0.0_dp, inv_n)
        j = 1
        do i=1,n
            do while (u(1) > weights_cdf(j))
                j = j + 1
                call debug_error_condition(j > n, 'uhh?')
            end do
            new_particles(:,i) = cartesian_particles(:,j)
            u = u + inv_n
        end do
        cartesian_particles = new_particles
        weights = inv_n
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

    impure subroutine propagate_particles(dt, jerk_sig, max_spd, n, cartesian_particles)
        real(dp), intent(in) :: dt, jerk_sig, max_spd
        integer, intent(in) :: n
        real(dp), intent(inout) :: cartesian_particles(:,:)
        real(dp) :: jx(n), jy(n)
        integer :: i
        call debug_error_condition(size(cartesian_particles, dim=2) /= n, 'mismatch in cartesian_particles shape')
        call random_normal(jx, 0.0_dp, jerk_sig)
        call random_normal(jy, 0.0_dp, jerk_sig)
        do concurrent (i=1:n)
            cartesian_particles(5,i) = cartesian_particles(5,i) + jx(i)
            cartesian_particles(6,i) = cartesian_particles(6,i) + jy(i)
            call fly_constant_acceleration(cartesian_particles(:,i), dt, max_spd)
        end do
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

    impure subroutine print_summary(msg, obs_cart, tgt_cart, n, cartesian_particles, weights)
        character(len=*), intent(in) :: msg
        real(dp), intent(in) :: obs_cart(:), tgt_cart(size(obs_cart))
        integer, intent(in) :: n
        real(dp), intent(in) :: cartesian_particles(size(obs_cart),n), weights(n)
        real(dp) :: tgt_pol(3), tgt_spd, tgt_acc, cart_spd(n), cart_acc(n), &
                    polar_particles(3,n), est_cart(size(obs_cart)), est_pol(3)
        call perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
        call convert_particles_cart2pol(obs_cart, n, cartesian_particles, polar_particles)
        call generate_state_estimate(est_cart, n, cartesian_particles, weights)
        call generate_state_estimate(est_pol, n, polar_particles, weights)
        tgt_spd = sqrt(tgt_cart(3)**2 + tgt_cart(4)**2)
        tgt_acc = sqrt(tgt_cart(5)**2 + tgt_cart(6)**2)
        cart_spd = sqrt(sum(cartesian_particles(3:4,:)**2, dim=1))
        cart_acc = sqrt(sum(cartesian_particles(5:6,:)**2, dim=1))
        write(*,'(a)') repeat('=', 32)
        write(*,'(a,f0.1)') msg//' :: Neff: ',neff(weights)
        write(*,'(5a26)') 'Quantity / ','Truth / ','Particle EST / ','Particle STD / ','Particle RMSE / '
        write(*,'(a26,3(f23.4,a3),f23.4)') 'range [NMI]: ',tgt_pol(1)*ft2nmi,' / ', &
                                                           est_pol(1)*ft2nmi,' / ', &
                                                           std(polar_particles(1,:))*ft2nmi,' / ', &
                                                           rmse(polar_particles(1,:), tgt_pol(1))*ft2nmi
        write(*,'(a26,3(f23.4,a3),f23.4)') 'bearing angle [DEG]: ',tgt_pol(2)*rad2deg_dp,' / ', &
                                                                   est_pol(2)*rad2deg_dp,' / ', &
                                                                   std(polar_particles(2,:))*rad2deg_dp,' / ', &
                                                                   rmse(polar_particles(2,:), tgt_pol(2))*rad2deg_dp
        write(*,'(a26,3(f23.4,a3),f23.4)') 'range-rate [FT/SEC]: ',tgt_pol(3),' / ', &
                                                                   est_pol(3),' / ', &
                                                                   std(polar_particles(3,:)),' / ', &
                                                                   rmse(polar_particles(3,:), tgt_pol(3))
        write(*,'(a26,3(f23.4,a3),f23.4)') 'x [NMI]: ',tgt_cart(1)*ft2nmi,' / ', &
                                                       est_cart(1)*ft2nmi,' / ', &
                                                       std(cartesian_particles(1,:))*ft2nmi,' / ', &
                                                       rmse(cartesian_particles(1,:), tgt_cart(1))*ft2nmi
        write(*,'(a26,3(f23.4,a3),f23.4)') 'y [NMI]: ',tgt_cart(2)*ft2nmi,' / ', &
                                                       est_cart(2)*ft2nmi,' / ', &
                                                       std(cartesian_particles(2,:))*ft2nmi,' / ', &
                                                       rmse(cartesian_particles(2,:), tgt_cart(2))*ft2nmi
        write(*,'(a26,3(f23.4,a3),f23.4)') 'vx [FT/SEC]: ',tgt_cart(3),' / ', &
                                                           est_cart(3),' / ', &
                                                           std(cartesian_particles(3,:)),' / ', &
                                                           rmse(cartesian_particles(3,:), tgt_cart(3))
        write(*,'(a26,3(f23.4,a3),f23.4)') 'vy [FT/SEC]: ',tgt_cart(4),' / ', &
                                                           est_cart(4),' / ', &
                                                           std(cartesian_particles(4,:)),' / ', &
                                                           rmse(cartesian_particles(4,:), tgt_cart(4))
        write(*,'(a26,3(f23.4,a3),f23.4)') 'speed [FT/SEC]: ',tgt_spd,' / ', &
                                                              vdot(cart_spd, weights),' / ', &
                                                              std(cart_spd),' / ', &
                                                              rmse(cart_spd, tgt_spd)
        write(*,'(a26,3(f23.4,a3),f23.4)') 'ax [FT/SEC**2]: ',tgt_cart(5),' / ', &
                                                              est_cart(5),' / ', &
                                                              std(cartesian_particles(5,:)),' / ', &
                                                              rmse(cartesian_particles(5,:), tgt_cart(5))
        write(*,'(a26,3(f23.4,a3),f23.4)') 'ay [FT/SEC**2]: ',tgt_cart(6),' / ', &
                                                              est_cart(6),' / ', &
                                                              std(cartesian_particles(6,:)),' / ', &
                                                              rmse(cartesian_particles(6,:), tgt_cart(6))
        write(*,'(a26,3(f23.4,a3),f23.4)') 'acceleration [FT/SEC**2]: ',tgt_acc,' / ', &
                                                                        vdot(cart_acc, weights),' / ', &
                                                                        std(cart_acc),' / ', &
                                                                        rmse(cart_acc, tgt_acc)
        write(*,'(a)') repeat('=', 32)
    end

end module pf


program ex_particle_filter
use, non_intrinsic :: pf
implicit none
    
    real(dp), parameter :: max_rng = 500.0_dp*nmi2ft, &
                           max_spd = 10000.0_dp, &
                           jerk_sig = 0.1_dp*g, &
                           obs_cart(6) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
                           dt = 1.0_dp, &
                           meas_sig(3) = [10.0_dp*nmi2ft, 1.0_dp*deg2rad_dp, 100.0_dp] ! poor measurements
!                           meas_sig(3) = [100.0_dp, 0.1_dp*deg2rad_dp, 10.0_dp] ! standard measurements
!                           meas_sig(3) = [1.0_dp, 0.001_dp, 1.0_dp] ! exquisite measurements

    integer :: ang, spd, num_particles        
    real(dp) :: tgt_cart(6), tgt_pol(3), t, meas(3), est_cart(6), pos_err, vel_err, acc_err
    real(dp), allocatable :: cartesian_particles(:,:), weights(:), polar_particles(:,:)

    do ang=15,15,5
        do spd=3000,3000
            do num_particles=10000000,10000000
                !! reset target for each run
                tgt_cart = 0.0_dp !! initialize all state components to zero (0.0)
                tgt_cart(1) = 100.0_dp*nmi2ft*cos(real(ang,dp)*deg2rad_dp) !! x position [ft]
                tgt_cart(2) = 100.0_dp*nmi2ft*sin(real(ang,dp)*deg2rad_dp) !! y position [ft]
                tgt_cart(3) = -real(spd,dp) !! vx velocity [ft/sec]
                tgt_cart(4) = 1000.0_dp !! vy velocity [ft/sec]
                tgt_cart(6) = -32.2_dp !! ay acceleration [ft/sec**2]
                call perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
                if (allocated(cartesian_particles)) deallocate(cartesian_particles)
                if (allocated(weights)) deallocate(weights)
                if (allocated(polar_particles)) deallocate(polar_particles)
                allocate(cartesian_particles(size(obs_cart),num_particles), &
                         weights(num_particles), &
                         polar_particles(3,num_particles))
                call initialize_particles(obs_cart, tgt_cart, meas_sig, max_rng, max_spd, &
                                          num_particles, cartesian_particles, weights)
                call convert_particles_cart2pol(obs_cart, num_particles, cartesian_particles, polar_particles)
                call calculate_weights(tgt_pol, meas_sig, num_particles, polar_particles, weights)
                call print_summary('Initialization Complete',obs_cart, tgt_cart, num_particles, cartesian_particles, weights)
                !! reset weights to 1/num_particles
                weights = 1.0_dp/real(num_particles,dp)
                t = 0.0_dp
                do while (tgt_cart(2) > 0.0_dp)
                    !! advance time and target independent of anything else
                    t = t + dt
                    call fly_constant_acceleration(tgt_cart, dt, max_spd)
                    !! generate new measurement from observer perspective
                    call generate_measurements(obs_cart, tgt_cart, meas_sig, max_rng, max_spd, 1, meas)
                    !! propagate current particles
                    call propagate_particles(dt, jerk_sig, max_spd, num_particles, cartesian_particles)
                    !! convert cartesian particles to polar representation and refine weights accordingly
                    call convert_particles_cart2pol(obs_cart, num_particles, cartesian_particles, polar_particles)
                    call calculate_weights(meas, meas_sig, num_particles, polar_particles, weights)
                    !! resample if Neff < 50% of num_particles
                    if (neff(weights) < 0.5_dp*real(num_particles,dp)) then
                        call resample_systematic(num_particles, cartesian_particles, weights)
                        if (debug) write(*,*) '--> WEIGHTS RESAMPLED'
                    end if
                    if (debug) then
                        call print_summary('tracking update',obs_cart, tgt_cart, num_particles, cartesian_particles, weights)
                    else
                        call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
                        pos_err = sqrt((tgt_cart(1) - est_cart(1))**2 + (tgt_cart(2) - est_cart(2))**2)
                        vel_err = sqrt((tgt_cart(3) - est_cart(3))**2 + (tgt_cart(4) - est_cart(4))**2)
                        acc_err = sqrt((tgt_cart(5) - est_cart(5))**2 + (tgt_cart(6) - est_cart(6))**2)
                        write(*,'(a,f0.1,a,f10.1,a,3(f0.1,a))') 't: ',t,', neff: ',neff(weights), &
                                                                ', pos_err: ',pos_err,' ft'// &
                                                                ', vel_err: ',vel_err,' ft/sec'// &
                                                                ', acc_err: ',acc_err,' ft/sec**2'
                    end if
                end do
                write(*,'(a,f0.1,a,f0.3,a,f0.3,a,f0.1,a,f0.1,a)') 'MEAS_SIG: ', &                      !! a
                                                                   meas_sig(1),' ft (', &              !! f0.1,a
                                                                   meas_sig(1)*ft2nmi,' NMI), ', &     !! f0.3,a
                                                                   meas_sig(2)*rad2deg_dp,' deg (', &  !! f0.3,a
                                                                   meas_sig(2)*1000.0_dp,' mrad), ', & !! f0.1,a
                                                                   meas_sig(3),' ft/sec'               !! f0.1,a
                write(*,'(7(a,f0.1),a)') '         TRUTH t: ',t, &
                                         ' :: x: ',tgt_cart(1)*ft2nmi,' NMI'// &
                                         ', y: ',tgt_cart(2)*ft2nmi,' NMI'//  &
                                         ', vx: ',tgt_cart(3),' ft/sec'// &
                                         ', vy: ',tgt_cart(4),' ft/sec'// &
                                         ', ax: ',tgt_cart(5),' ft/sec**2'// &
                                         ', ay: ',tgt_cart(6),' ft/sec**2'
                call generate_state_estimate(est_cart, num_particles, cartesian_particles, weights)
                write(*,'(7(a,f0.1),a)') 'STATE ESTIMATE t: ',t, &
                                         ' :: x: ',est_cart(1)*ft2nmi,' NMI'// &
                                         ', y: ',est_cart(2)*ft2nmi,' NMI'//  &
                                         ', vx: ',est_cart(3),' ft/sec'// &
                                         ', vy: ',est_cart(4),' ft/sec'// &
                                         ', ax: ',est_cart(5),' ft/sec**2'// &
                                         ', ay: ',est_cart(6),' ft/sec**2'
                pos_err = sqrt((tgt_cart(1) - est_cart(1))**2 + (tgt_cart(2) - est_cart(2))**2)
                vel_err = sqrt((tgt_cart(3) - est_cart(3))**2 + (tgt_cart(4) - est_cart(4))**2)
                acc_err = sqrt((tgt_cart(5) - est_cart(5))**2 + (tgt_cart(6) - est_cart(6))**2)
                write(*,'(a,f0.1,a,f10.1,a,3(f0.1,a))') '               t: ',t,', neff: ',neff(weights), &
                                                                     ', pos_err: ',pos_err,' ft'// &
                                                                     ', vel_err: ',vel_err,' ft/sec'// &
                                                                     ', acc_err: ',acc_err,' ft/sec**2'
            end do
        end do
    end do

end program ex_particle_filter
