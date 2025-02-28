module pf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: pi_dp, twopi_dp
use, non_intrinsic :: random, only: random_normal, random_uniform
implicit none
    
    real(dp), parameter :: global_minimum_range = 1.0_dp

    public :: generate_measurements, initialize_particles

contains

    impure subroutine generate_measurements(obs_cart, tgt_cart, meas_sig, no_meas_max_rng, no_meas_max_spd, n, polar_measurements)
        real(dp), intent(in) :: obs_cart(4), tgt_cart(4), meas_sig(3), no_meas_max_rng, no_meas_max_spd
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

    pure subroutine cart2pol_inner(obs_cart, tgt_cart, err_rng, err_ang, err_rngrt, tgt_pol)
        real(dp), intent(in) :: obs_cart(4), tgt_cart(4), err_rng, err_ang, err_rngrt
        real(dp), intent(out) :: tgt_pol(3)
        real(dp) :: rng, loshat(2), losv(2)
        rng = sqrt((tgt_cart(1) - obs_cart(1))**2 + (tgt_cart(2) - obs_cart(2))**2)
        tgt_pol(1) = max(global_minimum_range, rng + err_rng)
        loshat = (tgt_cart(1:2) - obs_cart(1:2))/rng
        tgt_pol(2) = mod(atan2(loshat(2), loshat(1)) + err_ang + pi_dp, twopi_dp) - pi_dp
        losv = tgt_cart(3:4) - obs_cart(3:4)
        tgt_pol(3) = losv(1)*loshat(1) + losv(2)*loshat(2) + err_rngrt
    end

    impure subroutine initialize_particles(obs_cart, tgt_cart, meas_sig, tgt_max_rng, tgt_max_spd, n, cartesian_particles)
        real(dp), intent(in) :: obs_cart(4), tgt_cart(4), meas_sig(3), tgt_max_rng, tgt_max_spd
        integer, intent(in) :: n
        real(dp), intent(out) :: cartesian_particles(4,n)
        real(dp) :: polar_measurements(3,n), tgt_spd_scale(n)
        integer :: i
        call generate_measurements(obs_cart, tgt_cart, meas_sig, tgt_max_rng, tgt_max_spd, n, polar_measurements)
        call random_uniform(tgt_spd_scale, 0.0_dp, 1.0_dp)
        do concurrent (i=1:n)
            call pol2cart_inner(obs_cart, tgt_max_rng, tgt_max_spd, tgt_spd_scale(i), &
                                polar_measurements(:,i), cartesian_particles(:,i))
        end do
    end

    pure subroutine pol2cart_inner(obs_cart, tgt_max_rng, tgt_max_spd, tgt_spd_scale, tgt_pol, tgt_cart)
        real(dp), intent(in) :: obs_cart(4), tgt_max_rng, tgt_max_spd, tgt_spd_scale, tgt_pol(3)
        real(dp), intent(out) :: tgt_cart(4)
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

end module pf


program ex_particle_filter
use, non_intrinsic :: pf
implicit none
end program ex_particle_filter
