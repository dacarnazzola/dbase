module pf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: pi_dp, twopi_dp, deg2rad_dp, rad2deg_dp
use, non_intrinsic :: random, only: random_normal, random_uniform
use, non_intrinsic :: statistics, only: dsum, avg, std
use, non_intrinsic :: vector_math, only: vdot
implicit none
private
    
    real(dp), parameter :: global_minimum_range = 1.0_dp
    real(dp), parameter :: nmi2ft = 1852.0_dp*100.0_dp/2.54_dp/12.0_dp
    real(dp), parameter :: ft2nmi = 12.0_dp*2.54_dp/100.0_dp/1852.0_dp

    public :: dp, nmi2ft, ft2nmi, deg2rad_dp, rad2deg_dp, &
              perfect_cart2pol, generate_measurements, initialize_particles, convert_particles_cart2pol, &
              generate_state_estimate, calculate_weights, &
              rmse, neff, &
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
        real(dp) :: polar_measurements(3,n), tgt_spd_scale(n)
        integer :: i
        call generate_measurements(obs_cart, tgt_cart, meas_sig, tgt_max_rng, tgt_max_spd, n, polar_measurements)
        call random_uniform(tgt_spd_scale, 0.0_dp, 1.0_dp)
        do concurrent (i=1:n)
            call pol2cart_inner(obs_cart, tgt_max_rng, tgt_max_spd, tgt_spd_scale(i), &
                                polar_measurements(:,i), cartesian_particles(:,i))
            if (size(obs_cart) == 6) cartesian_particles(5:6,i) = 0.0_dp
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

    pure subroutine generate_state_estimate(state_estimate, n, cartesian_particles, weights)
        real(dp), intent(out) :: state_estimate(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: cartesian_particles(size(state_estimate),n), weights(n)
        integer :: i
        do concurrent (i=1:size(state_estimate))
            state_estimate(i) = vdot(cartesian_particles(i,:), weights)
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
    
    real(dp), parameter :: max_rng = 500.0_dp*nmi2ft, &
                           max_spd = 10000.0_dp, &
                           obs_cart(6) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
                           meas_sig(3) = [10.0_dp*nmi2ft, 1.0_dp*deg2rad_dp, 100.0_dp] ! poor measurements
!                           meas_sig(3) = [100.0_dp, 0.1_dp*deg2rad_dp, 10.0_dp] ! standard measurements
!                           meas_sig(3) = [1.0_dp, 0.001_dp, 1.0_dp] ! exquisite measurements

    integer :: ang, spd, num_particles        
    real(dp) :: tgt_cart(6), tgt_pol(3)
    real(dp), allocatable :: cartesian_particles(:,:), weights(:), cart_spd(:), polar_particles(:,:)

    tgt_cart = 0.0_dp
    do ang=0,90,30
        tgt_cart(1) = 100.0_dp*nmi2ft*cos(real(ang,dp)*deg2rad_dp)
        tgt_cart(2) = 100.0_dp*nmi2ft*sin(real(ang,dp)*deg2rad_dp)
        do spd=900,900
            tgt_cart(3) = -real(spd,dp)
            call perfect_cart2pol(obs_cart, tgt_cart, tgt_pol)
            do num_particles=10000000,10000000
                if (allocated(cartesian_particles)) deallocate(cartesian_particles)
                if (allocated(weights)) deallocate(weights)
                if (allocated(cart_spd)) deallocate(cart_spd)
                if (allocated(polar_particles)) deallocate(polar_particles)
                allocate(cartesian_particles(size(obs_cart),num_particles), &
                         weights(num_particles), &
                         cart_spd(num_particles), &
                         polar_particles(3,num_particles))
                call initialize_particles(obs_cart, tgt_cart, meas_sig, max_rng, max_spd, &
                                          num_particles, cartesian_particles, weights)
                cart_spd = sqrt(sum(cartesian_particles(3:4,:)**2, dim=1))
                call convert_particles_cart2pol(obs_cart, num_particles, cartesian_particles, polar_particles)
                write(*,'(5a21)') 'Quantity / ','Truth / ','Particle AVG / ','Particle STD / ','Particle RMSE / '
                write(*,'(a21,3(f18.4,a3),f18.4)') 'range [NMI]: ',tgt_pol(1)*ft2nmi,' / ', &
                                                                   avg(polar_particles(1,:))*ft2nmi,' / ', &
                                                                   std(polar_particles(1,:))*ft2nmi,' / ', &
                                                                   rmse(polar_particles(1,:), tgt_pol(1))*ft2nmi
                write(*,'(a21,3(f18.4,a3),f18.4)') 'bearing angle [DEG]: ',tgt_pol(2)*rad2deg_dp,' / ', &
                                                                           avg(polar_particles(2,:))*rad2deg_dp,' / ', &
                                                                           std(polar_particles(2,:))*rad2deg_dp,' / ', &
                                                                           rmse(polar_particles(2,:), tgt_pol(2))*rad2deg_dp
                write(*,'(a21,3(f18.4,a3),f18.4)') 'range-rate [FT/SEC]: ',tgt_pol(3),' / ', &
                                                                           avg(polar_particles(3,:)),' / ', &
                                                                           std(polar_particles(3,:)),' / ', &
                                                                           rmse(polar_particles(3,:), tgt_pol(3))
                write(*,'(a21,3(f18.4,a3),f18.4)') 'x [NMI]: ',tgt_cart(1)*ft2nmi,' / ', &
                                                               avg(cartesian_particles(1,:))*ft2nmi,' / ', &
                                                               std(cartesian_particles(1,:))*ft2nmi,' / ', &
                                                               rmse(cartesian_particles(1,:), tgt_cart(1))*ft2nmi
                write(*,'(a21,3(f18.4,a3),f18.4)') 'y [NMI]: ',tgt_cart(2)*ft2nmi,' / ', &
                                                               avg(cartesian_particles(2,:))*ft2nmi,' / ', &
                                                               std(cartesian_particles(2,:))*ft2nmi,' / ', &
                                                               rmse(cartesian_particles(2,:), tgt_cart(2))*ft2nmi
                write(*,'(a21,3(f18.4,a3),f18.4)') 'vx [FT/SEC]: ',tgt_cart(3),' / ', &
                                                                   avg(cartesian_particles(3,:)),' / ', &
                                                                   std(cartesian_particles(3,:)),' / ', &
                                                                   rmse(cartesian_particles(3,:), tgt_cart(3))
                write(*,'(a21,3(f18.4,a3),f18.4)') 'vy [FT/SEC]: ',tgt_cart(4),' / ', &
                                                                   avg(cartesian_particles(4,:)),' / ', &
                                                                   std(cartesian_particles(4,:)),' / ', &
                                                                   rmse(cartesian_particles(4,:), tgt_cart(4))
                write(*,'(a21,3(f18.4,a3),f18.4)') 'speed [FT/SEC]: ',real(spd,dp),' / ', &
                                                                      avg(cart_spd),' / ', &
                                                                      std(cart_spd),' / ', &
                                                                      rmse(cart_spd, real(spd,dp))
                write(*,'(a)') ''
            end do
        end do
    end do

end program ex_particle_filter
