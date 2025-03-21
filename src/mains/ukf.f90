module ukf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: nmi2ft_dp
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: vector_math, only: vmag
use, non_intrinsic :: matrix_math, only: chol
implicit none
private

    real(dp), parameter :: g = 32.2_dp
    real(dp), parameter :: global_minimum_range = 1.0e-6_dp

    type :: sr_ukf_type
        !! default values for unscented transform given below, alternate defaults commented out
        real(dp) :: ut_kappa = 0.0_dp ! 3.0_dp - real(size(state_estimate),dp)
        real(dp) :: ut_alpha = 1.0_dp ! 1.0_dp
        real(dp) :: ut_beta = 2.0_dp
        !! unscented transform lambda is calculated once during initialization based on kappa, alpha, and beta
        real(dp) :: ut_lambda
        !! state_estimate defaults if no measurement information is provided
        real(dp) :: default_range = 500.0_dp*nmi2ft_dp
        real(dp) :: default_velocity(2) = 0.0_dp
        real(dp) :: default_acceleration(2) = 0.0_dp
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
        real(dp), intent(in), optional :: k, a, b, def_rng, def_vel(2), def_acc(2), max_vel, max_acc
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
        call debug_error_condition(vmag(filter%default_velocity) > filter%maximum_velocity, &
                                   'invalid velocity initialization, magnitude of default > maximum ILLEGAL')
        call debug_error_condition(vmag(filter%default_acceleration) > filter%maximum_acceleration, &
                                   'invalid acceleration initialization, magnitude of default > maximum ILLEGAL')
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
        if ((meas_sig(3) > 0.0_dp) .and. (meas_sig(4) > 0.0_dp)) then !! measurement contains range rate AND angle rate information
            vt = rng_use*meas(4)
            filter%state_estimate(3) = obs(3) + meas(3)*cos_ang - vt*sin_ang
            filter%state_estimate(4) = obs(4) + meas(3)*sin_ang + vt*cos_ang 
        else if (meas_sig(3) > 0.0_dp) then !! measurement contains ONLY range rate information
            filter%state_estimate(3) = obs(3) + meas(3)*cos_ang
            filter%state_estimate(4) = obs(4) + meas(3)*sin_ang
        else if (meas_sig(4) > 0.0_dp) then !! measurement contains ONLY angle rate information
            vt = rng_use*meas(4)
            filter%state_estimate(3) = obs(3) - vt*sin_ang
            filter%state_estimate(4) = obs(4) + vt*cos_ang 
        else !! measurement contains NO velocity information
            filter%state_estimate(3:4) = filter%default_velocity
        end if
        !! make sure state_estimate does not violate maximum_velocity constraint
        trk_spd = sqrt(filter%state_estimate(3)**2 + filter%state_estimate(4)**2)
        if (trk_spd > filter%maximum_velocity) then
            if ((meas_sig(1) > 0.0_dp) .or. (meas_sig(4) < 0.0_dp)) then !! measurement contains range OR lacks angle rate information
                !! simply scale estimated velocity if range is measured or there is no angle rate contribution
                filter%state_estimate(3:4) = filter%state_estimate(3:4)/(trk_spd/filter%maximum_velocity)
            else !! measurement LACKS range information and CONTAINS angle rate information
                !! first, remove current angle rate contribution to velocity
                filter%state_estimate(3) = filter%state_estimate(3) + vt*sin_ang
                filter%state_estimate(4) = filter%state_estimate(4) - vt*cos_ang
                !! next, calculate the maximum tangential velocity and corrected state estimate velocity components
                trk_spd = sqrt(filter%state_estimate(3)**2 + filter%state_estimate(4)**2)
                vt = sqrt(max(0.0_dp, filter%maximum_velocity**2 - trk_spd**2))
                filter%state_estimate(3) = filter%state_estimate(3) - vt*sin_ang
                filter%state_estimate(4) = filter%state_estimate(4) + vt*cos_ang
                !! next, calculate rng_use based on the maximum tangential velocity and measured angle rate
                if (.not.nearly(meas(4), 0.0_dp)) then
                    rng_use = vt/meas(4)
                    !! next, reposition state estimate position components based on updated rng_use
                    filter%state_estimate(1) = obs(1) + rng_use*cos_ang
                    filter%state_estimate(2) = obs(2) + rng_use*sin_ang
                end if
            end if
        end if
        ! assume no acceleration
        filter%state_estimate(5:6) = filter%default_acceleration
        ! initialize covariance as J * R * transpose(J), where J is the Jacobian of the state space and R is the measurement covariance
        call calculate_sqrt_jrjt(obs, filter%state_estimate, meas_sig, filter%covariance_square_root)
    end

    pure subroutine calculate_sqrt_jrjt(obs, tgt, meas_sig, sqrt_cov)
        real(dp), intent(in) :: obs(6), tgt(6), meas_sig(4)
        real(dp), intent(out) :: sqrt_cov(6,6)
        integer :: i, meas_dim
        real(dp) :: local_obs2tgt_pol(4), cos_ang, sin_ang, j(6,4), r(4,4), p(6,6)
        call cart2pol(obs, tgt, local_obs2tgt_pol)
        cos_ang = cos(local_obs2tgt_pol(2))
        sin_ang = sin(local_obs2tgt_pol(2))
        r = 0.0_dp
        meas_dim = 0
        do i=1,size(meas_sig)
            if (meas_sig(i) > 0.0_dp) then
                meas_dim = meas_dim + 1
                r(meas_dim,meas_dim) = meas_sig(i)
                select case (i)
                    case (1) ! d/d_range
                        j(1,meas_dim) = cos_ang
                        j(2,meas_dim) = sin_ang
                        j(3,meas_dim) = -local_obs2tgt_pol(4)*sin_ang
                        j(4,meas_dim) = local_obs2tgt_pol(4)*cos_ang
                        j(5,meas_dim) = 0.0_dp
                        j(6,meas_dim) = 0.0_dp
                    case (2) ! d/d_angle
                        j(1,meas_dim) = -local_obs2tgt_pol(1)*sin_ang
                        j(2,meas_dim) = local_obs2tgt_pol(1)*cos_ang
                        j(3,meas_dim) = -local_obs2tgt_pol(3)*sin_ang - local_obs2tgt_pol(1)*local_obs2tgt_pol(4)*cos_ang
                        j(4,meas_dim) = local_obs2tgt_pol(3)*cos_ang - local_obs2tgt_pol(1)*local_obs2tgt_pol(4)*sin_ang
                        j(5,meas_dim) = 0.0_dp
                        j(6,meas_dim) = 0.0_dp
                    case (3) ! d/d_range_rate
                        j(1,meas_dim) = 0.0_dp
                        j(2,meas_dim) = 0.0_dp
                        j(3,meas_dim) = cos_ang
                        j(4,meas_dim) = sin_ang
                        j(5,meas_dim) = 0.0_dp
                        j(6,meas_dim) = 0.0_dp
                    case (4) ! d/d_angle_rate
                        j(1,meas_dim) = 0.0_dp
                        j(2,meas_dim) = 0.0_dp
                        j(3,meas_dim) = -local_obs2tgt_pol(1)*sin_ang
                        j(4,meas_dim) = local_obs2tgt_pol(1)*cos_ang
                        j(5,meas_dim) = 0.0_dp
                        j(6,meas_dim) = 0.0_dp
                    case default
                        error stop 'this block should be unreachable'
                end select
            end if
        end do
        p = matmul(j, matmul(r, transpose(j))) ! P = J * R * transpose(J)
        call chol(p, sqrt_cov) ! S = sqrt(P)
    end

    pure subroutine cart2pol(obs, tgt, pol)
        real(dp), intent(in) :: obs(6), tgt(6)
        real(dp), intent(out) :: pol(4)
        real(dp) :: los(2), losv(2)
        los = tgt(1:2) - obs(1:2)
        pol(1) = max(global_minimum_range, sqrt(los(1)**2 + los(2)**2))
        pol(2) = atan2(los(2), los(1))
        losv = tgt(3:4) - obs(3:4)
        if (pol(1) > global_minimum_range) then
            pol(3) = (losv(1)*los(1) + losv(2)*los(2))/pol(1)
            pol(4) = (losv(2)*los(1) - losv(1)*los(2))/(pol(1)**2)
        else
            pol(3:4) = 0.0_dp
        end if
    end

end module ukf


program ex_ukf
use, non_intrinsic :: ukf
implicit none
end program ex_ukf
