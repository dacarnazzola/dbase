module ukf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: nmi2ft_dp
use, non_intrinsic :: system, only: debug_error_condition
implicit none
private

    real(dp), parameter :: g = 32.2_dp

    type :: sr_ukf_type
        !! default values for unscented transform given below, alternate defaults commented out
        real(dp) :: ut_kappa = 0.0_dp ! 3.0_dp - real(size(state_estimate),dp)
        real(dp) :: ut_alpha = 1.0_dp ! 1.0_dp
        real(dp) :: ut_beta = 2.0_dp
        !! unscented transform lambda is calculated once during initialization based on kappa, alpha, and beta
        real(dp) :: ut_lambda
        !! state_estimate defaults if no measurement information is provided
        real(dp) :: default_range = 500.0_dp*nmi2ft_dp
        real(dp) :: default_velocity = 0.0_dp
        real(dp) :: default_acceleration = 0.0_dp
        !! state_estimate maximum values
        real(dp) :: maximum_velocity = 10000.0_dp
        real(dp) :: maximum_acceleration = 9.0_dp*g
        !! state_estimate and covariance_square_root have no sensible default values
        real(dp) :: state_estimate(6)
        real(dp) :: covariance_square_root(6,6)
    end type sr_ukf_type

    public :: sr_ukf_type, initialize_sr_ukf

contains

    pure subroutine initialize_sr_ukf(obs, meas, meas_sig, filter, k, a, b, def_rng, def_vel, def_acc, max_vel, max_acc)
        real(dp), intent(in) :: obs(6), meas(4), meas_sig(4)
        type(sr_ukf_type), intent(out) :: filter
        real(dp), intent(in), optional :: k, a, b, def_rng, def_vel, def_acc, max_vel, max_acc
        real(dp) :: trk_n, rng_use, cos_ang, sin_ang, vt, trk_spd
        call debug_error_condition(meas_sig(2) <= 0.0_dp, &
                                   'track initialization not implemented for measurements lacking angle information')
        ! override any filter defaults with provided optional values
        if (present(k)) filter%ut_kappa = k
        if (present(a)) filter%ut_alpha = a
        if (present(b)) filter%ut_beta = b
        trk_n = real(size(filter%state_estimate), kind=dp)
        filter%ut_lambda = filter%ut_alpha**2 * (trk_n + filter%ut_kappa) - trk_n
        if (present(def_rng)) filter%default_range = def_rng
        if (present(def_vel)) filter%default_velocity = def_vel
        if (present(def_acc)) filter%default_acceleration = def_acc
        if (present(max_vel)) filter%maximum_velocity = max_vel
        if (present(max_acc)) filter%maximum_acceleration = max_acc
        call debug_error_condition(filter%default_velocity > filter%maximum_velocity, &
                                   'invalid velocity initialization, default > maximum ILLEGAL')
        call debug_error_condition(filter%default_acceleration > filter%maximum_acceleration, &
                                   'invalid acceleration initialization, default > maximum ILLEGAL')
        ! initialize track estimate based on observer state and measurement
        cos_ang = cos(meas(2))
        sin_ang = sin(meas(2))
        if (meas_sig(1) > 0.0_dp) then !! measurement contains range information
            rng_use = meas(1)
        else
            rng_use = filter%default_range
        end if
        filter%state_estimate(1) = obs(1) + rng_use*cos_ang
        filter%state_estimate(2) = obs(2) + rng_use*sin_ang
        if (meas_sig(3) > 0.0_dp) then !! measurement contains range rate information
            filter%state_estimate(3) = obs(3) + meas(3)*cos_ang
            filter%state_estimate(4) = obs(4) + meas(3)*sin_ang 
        else
            filter%state_estimate(3) = filter%default_velocity
            filter%state_estimate(4) = filter%default_velocity
        end if
        if (meas_sig(4) > 0.0_dp) then !! measurement contains angle rate information
            vt = rng_use*meas(4)
            filter%state_estimate(3) = filter%state_estimate(3) - vt*sin_ang
            filter%state_estimate(4) = filter%state_estimate(4) + vt*cos_ang
        end if
        !! make sure state_estimate does not violate maximum_velocity constraint
        trk_spd = sqrt(filter%state_estimate(3)**2 + filter%state_estimate(4)**2)
        if (trk_spd > filter%maximum_velocity) then
            if ((meas_sig(1) > 0.0_dp) .or. (meas_sig(4) < 0.0_dp)) then
                !! simply scale estimated velocity if range is measured or there is no angle rate contribution
                filter%state_estimate(3:4) = filter%state_estimate(3:4)/(trk_spd/filter%maximum_velocity)
            else if ((meas_sig(1) < 0.0_dp) .and. (meas_sig(4) > 0.0_dp)) then!! measurement LACKS range information and CONTAINS angle rate information
                !! first, remove current angle rate contribution to velocity
                filter%state_estimate(3) = filter%state_estimate(3) + vt*sin_ang
                filter%state_estimate(4) = filter%state_estimate(4) - vt*cos_ang
                !! next, calculate the maximum tangential velocity and corrected state estimate velocity components
                trk_spd = sqrt(filter%state_estimate(3)**2 + filter%state_estimate(4)**2)
                vt = sqrt(max(0.0_dp, filter%maximum_velocity**2 - trk_spd**2))
                filter%state_estimate(3) = filter%state_estimate(3) - vt*sin_ang
                filter%state_estimate(4) = filter%state_estimate(4) + vt*cos_ang
                !! next, calculate rng_use based on the maximum tangential velocity and measured angle rate
                rng_use = vt/meas(4)
                !! next, reposition state estimate position components based on updated rng_use
                filter%state_estimate(1) = obs(1) + rng_use*cos_ang
                filter%state_estimate(2) = obs(2) + rng_use*sin_ang
            else
                error stop 'logic error, this block should be unreachable, check combinations of meas_sig'
            end if
        end if
        ! assume no acceleration
        filter%state_estimate(5) = filter%default_acceleration
        filter%state_estimate(6) = filter%default_acceleration
    end

end module ukf


program ex_ukf
use, non_intrinsic :: ukf
implicit none
end program ex_ukf
