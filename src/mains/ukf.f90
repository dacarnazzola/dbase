module ukf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: nmi2ft_dp
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: vector_math, only: vmag
use, non_intrinsic :: matrix_math, only: chol, qr
use, non_intrinsic :: array_utils, only: broadcast_sub
implicit none
private

    real(dp), parameter :: g = 32.2_dp
    real(dp), parameter :: global_minimum_range = 1.0e-6_dp

    type :: sr_ukf_type
        !! an uninitialized sr_ukf_type will have negative negative time, default initialization will set this to 0.0
        real(dp) :: state_estimate_time = -1.0_dp
        !! model process noise should handle maneuvers up to process_noise ft/sec**3
        real(dp) :: process_noise = 3.0_dp*g
        !! state_estimate maximum values
        real(dp) :: maximum_velocity = 10000.0_dp
        real(dp) :: maximum_acceleration = 9.0_dp*g
        !! unscented transform weights are calculated once during initialization based on kappa, alpha, and beta
        real(dp) :: wx_1
        real(dp) :: wc_1
        real(dp) :: w_2_2n1
        real(dp) :: ut_gamma
        !! state_estimate and covariance_square_root have no sensible default values
        real(dp) :: state_estimate(6)
        real(dp) :: covariance_square_root(6,6)
    end type sr_ukf_type

    interface
        pure subroutine dynamics_model(state, dt, opt_data)
        import dp
        implicit none
            real(dp), intent(inout) :: state(:)
            real(dp), intent(in) :: dt
            real(dp), intent(in), optional :: opt_data(:)
        end
    end interface

    public :: sr_ukf_type, initialize_sr_ukf, filter_time_update, dynamics_ballistic

contains

    pure subroutine initialize_sr_ukf(obs, meas, meas_sig, filter, t, noise, max_vel, def_vel, max_acc, def_acc, k, a, b, def_rng)
        real(dp), intent(in) :: obs(6), meas(4), meas_sig(4)
        type(sr_ukf_type), intent(out) :: filter
        real(dp), intent(in), optional :: t, noise, max_vel, max_acc, k, a, b, def_rng, def_vel(2), def_acc(2)
        real(dp) :: trk_n, rng_use, cos_ang, sin_ang, v_t, trk_spd, jrjt(6,6), dim_var, ut_kappa, ut_alpha, ut_beta, ut_lambda, &
                    default_range, default_velocity(2), default_acceleration(2)
        call debug_error_condition(meas_sig(2) <= 0.0_dp, &
                                   'track initialization not implemented for measurements lacking angle information')
        !! override any filter defaults with provided optional values
        if (present(t)) then
            filter%state_estimate_time = t
            call debug_error_condition(filter%state_estimate_time < 0.0_dp, 'state_estimate_time must be >= 0.0')
        else
            filter%state_estimate_time = 0.0_dp
        end if
        if (present(noise)) then
            filter%process_noise = noise
            call debug_error_condition(filter%process_noise <= 0.0_dp, 'process_noise must be >= 0.0')
        end if
        if (present(max_vel)) filter%maximum_velocity = max_vel
        if (present(def_vel)) then
            default_velocity = def_vel
        else
            default_velocity = 0.0_dp
        end if
        call debug_error_condition(vmag(default_velocity) > filter%maximum_velocity, &
                                   'magnitude of default_velocity must be <= maximum_velocity')
        if (present(max_acc)) filter%maximum_acceleration = max_acc
        if (present(def_acc)) then
            default_acceleration = def_acc
        else
            default_acceleration = 0.0_dp
        end if
        call debug_error_condition(vmag(default_acceleration) > filter%maximum_acceleration, &
                                   'magnitude of default_acceleration must be <= maximum_acceleration')
        if (present(k)) then
            ut_kappa = k
        else
            ut_kappa = 0.0_dp !! alternate, 3.0_dp - real(size(filter%state_estimate), kind=dp)
        end if
        if (present(a)) then
            ut_alpha = a
        else
            ut_alpha = 1.0_dp !! alternate, 1.0_dp
        end if
        if (present(b)) then
            ut_beta = b
        else
            ut_beta = 2.0_dp !! 2.0_dp is ideal for Gaussian distribution
        end if
        if (present(def_rng)) then
            default_range = def_rng
        else
            default_range = 500.0_dp*nmi2ft_dp
        end if
        trk_n = size(filter%state_estimate, kind=dp)
        ut_lambda = ut_alpha**2 * (trk_n + ut_kappa) - trk_n
        filter%wx_1 = ut_lambda/(trk_n + ut_lambda)
        filter%wc_1 = filter%wx_1 + 1.0_dp - ut_alpha**2 + ut_beta
        call debug_error_condition(filter%wc_1 < 0.0_dp, 'negative weights for covariance matrix not implemented')
        filter%wc_1 = sqrt(filter%wc_1) !! wc_1 only ever used under square root, requiring the check above - compute once here
        filter%w_2_2n1 = 0.5_dp/(trk_n + ut_lambda)
        filter%ut_gamma = sqrt(trk_n + ut_lambda)
        !! initialize track estimate based on observer state and measurement
        cos_ang = cos(meas(2))
        sin_ang = sin(meas(2))
        if (meas_sig(1) > 0.0_dp) then !! measurement contains range information
            rng_use = meas(1)
        else
            rng_use = default_range
        end if
        filter%state_estimate(1) = obs(1) + rng_use*cos_ang
        filter%state_estimate(2) = obs(2) + rng_use*sin_ang
        if ((meas_sig(3) > 0.0_dp) .and. (meas_sig(4) > 0.0_dp)) then !! measurement contains range rate AND angle rate information
            v_t = rng_use*meas(4)
            filter%state_estimate(3) = obs(3) + meas(3)*cos_ang - v_t*sin_ang
            filter%state_estimate(4) = obs(4) + meas(3)*sin_ang + v_t*cos_ang 
        else if (meas_sig(3) > 0.0_dp) then !! measurement contains ONLY range rate information
            filter%state_estimate(3) = obs(3) + meas(3)*cos_ang
            filter%state_estimate(4) = obs(4) + meas(3)*sin_ang
        else if (meas_sig(4) > 0.0_dp) then !! measurement contains ONLY angle rate information
            v_t = rng_use*meas(4)
            filter%state_estimate(3) = obs(3) - v_t*sin_ang
            filter%state_estimate(4) = obs(4) + v_t*cos_ang 
        else !! measurement contains NO velocity information
            filter%state_estimate(3:4) = default_velocity
        end if
        !! make sure state_estimate does not violate maximum_velocity constraint
        trk_spd = sqrt(filter%state_estimate(3)**2 + filter%state_estimate(4)**2)
        if (trk_spd > filter%maximum_velocity) then
            if ((meas_sig(1) > 0.0_dp) .or. (meas_sig(4) < 0.0_dp)) then !! measurement contains range OR lacks angle rate information
                !! simply scale estimated velocity if range is measured or there is no angle rate contribution
                filter%state_estimate(3:4) = filter%state_estimate(3:4)/(trk_spd/filter%maximum_velocity)
            else !! measurement LACKS range information and CONTAINS angle rate information
                !! first, remove current angle rate contribution to velocity
                filter%state_estimate(3) = filter%state_estimate(3) + v_t*sin_ang
                filter%state_estimate(4) = filter%state_estimate(4) - v_t*cos_ang
                !! next, calculate the maximum tangential velocity and corrected state estimate velocity components
                trk_spd = sqrt(filter%state_estimate(3)**2 + filter%state_estimate(4)**2)
                v_t = sqrt(max(0.0_dp, filter%maximum_velocity**2 - trk_spd**2))
                filter%state_estimate(3) = filter%state_estimate(3) - v_t*sin_ang
                filter%state_estimate(4) = filter%state_estimate(4) + v_t*cos_ang
                !! next, calculate rng_use based on the maximum tangential velocity and measured angle rate
                if (.not.nearly(meas(4), 0.0_dp)) then
                    rng_use = v_t/meas(4)
                    !! next, reposition state estimate position components based on updated rng_use
                    filter%state_estimate(1) = obs(1) + rng_use*cos_ang
                    filter%state_estimate(2) = obs(2) + rng_use*sin_ang
                end if
            end if
        end if
        !! assume no acceleration
        filter%state_estimate(5:6) = default_acceleration
        !! initialize covariance as J * R * transpose(J), where J is the Jacobian of the state space and R is the measurement covariance
        call generate_jrjt(obs, filter%state_estimate, meas_sig, jrjt)
        if (meas_sig(1) < 0.0_dp) then
            dim_var = (0.5_dp*rng_use)**2/3.0_dp
            jrjt(1,1) = jrjt(1,1) + dim_var*cos_ang**2
            jrjt(2,1) = jrjt(2,1) + dim_var*cos_ang*sin_ang
            jrjt(1,2) = jrjt(1,2) + dim_var*cos_ang*sin_ang
            jrjt(2,2) = jrjt(2,2) + dim_var*sin_ang**2
        end if
        if ((meas_sig(3) < 0.0_dp) .and. (meas_sig(4) < 0.0_dp)) then
            dim_var = filter%maximum_velocity**2/3.0_dp
            jrjt(3,3) = jrjt(3,3) + dim_var
            jrjt(4,4) = jrjt(4,4) + dim_var
        else if (meas_sig(3) < 0.0_dp) then !! missing range rate information, add radial velocity uncertainty
            dim_var = filter%maximum_velocity**2/3.0_dp
            jrjt(3,3) = jrjt(3,3) + dim_var*cos_ang**2
            jrjt(4,3) = jrjt(4,3) + dim_var*cos_ang*sin_ang
            jrjt(3,4) = jrjt(3,4) + dim_var*cos_ang*sin_ang
            jrjt(4,4) = jrjt(4,4) + dim_var*sin_ang**2
        else if (meas_sig(4) < 0.0_dp) then !! missing angle rate information, add tangential velocity uncertainty
            dim_var = filter%maximum_velocity**2/3.0_dp
            jrjt(3,3) = jrjt(3,3) + dim_var*sin_ang**2
            jrjt(4,3) = jrjt(4,3) - dim_var*cos_ang*sin_ang
            jrjt(3,4) = jrjt(3,4) - dim_var*cos_ang*sin_ang
            jrjt(4,4) = jrjt(4,4) + dim_var*cos_ang**2
        end if
        dim_var = filter%maximum_acceleration**2/3.0_dp
        jrjt(5,5) = jrjt(5,5) + dim_var
        jrjt(6,6) = jrjt(6,6) + dim_var
        call chol(jrjt, filter%covariance_square_root)
    end

    pure subroutine generate_jrjt(obs, tgt, meas_sig, jrjt)
        real(dp), intent(in) :: obs(6), tgt(6), meas_sig(4)
        real(dp), intent(out) :: jrjt(6,6)
        integer :: i, meas_dim
        real(dp) :: local_obs2tgt_pol(4), cos_ang, sin_ang, j(6,4), r(4,4)
        call cart2pol(obs, tgt, local_obs2tgt_pol)
        cos_ang = cos(local_obs2tgt_pol(2))
        sin_ang = sin(local_obs2tgt_pol(2))
        j = 0.0_dp
        r = 0.0_dp
        meas_dim = 0
        do i=1,size(meas_sig)
            if (meas_sig(i) > 0.0_dp) then
                meas_dim = meas_dim + 1
                r(meas_dim,meas_dim) = meas_sig(i)
                select case (i)
                    case (1) !! d/d_range
                        j(1,meas_dim) = cos_ang
                        j(2,meas_dim) = sin_ang
                        j(3,meas_dim) = -local_obs2tgt_pol(4)*sin_ang
                        j(4,meas_dim) = local_obs2tgt_pol(4)*cos_ang
                    case (2) !! d/d_angle
                        j(1,meas_dim) = -local_obs2tgt_pol(1)*sin_ang
                        j(2,meas_dim) = local_obs2tgt_pol(1)*cos_ang
                        j(3,meas_dim) = -local_obs2tgt_pol(3)*sin_ang - local_obs2tgt_pol(1)*local_obs2tgt_pol(4)*cos_ang
                        j(4,meas_dim) = local_obs2tgt_pol(3)*cos_ang - local_obs2tgt_pol(1)*local_obs2tgt_pol(4)*sin_ang
                    case (3) !! d/d_range_rate
                        j(3,meas_dim) = cos_ang
                        j(4,meas_dim) = sin_ang
                    case (4) !! d/d_angle_rate
                        j(3,meas_dim) = -local_obs2tgt_pol(1)*sin_ang
                        j(4,meas_dim) = local_obs2tgt_pol(1)*cos_ang
                end select
            end if
        end do
        jrjt = matmul(j(:,1:meas_dim), matmul(r(1:meas_dim,1:meas_dim), transpose(j(:,1:meas_dim)))) !! P = J * R * transpose(J)
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

    pure subroutine filter_time_update(filter, t, state_dynamics_model)
        type(sr_ukf_type), intent(inout) :: filter
        real(dp), intent(in) :: t
        procedure(dynamics_model) :: state_dynamics_model
        real(dp) :: dt, sigma_points(size(filter%state_estimate),2*size(filter%state_estimate)+1), &
                    amat(size(filter%state_estimate),3*size(filter%state_estimate)+1), sqrt_wc, &
                    unused_tau(min(size(sigma_points,dim=1),size(sigma_points,dim=2))), &
                    amat_t(size(amat,dim=2),size(amat,dim=1))
        integer :: i
        call debug_error_condition(filter%state_estimate_time < 0.0_dp, 'filter is not initialized')
        dt = t - filter%state_estimate_time
        call debug_error_condition(dt < 0.0_dp, 'negative dt time update not implemented')
        filter%state_estimate_time = t
        if (nearly(dt, 0.0_dp)) return !! return early if dt is approximately 0.0
        !! generate sigma points
        call generate_sigma_points(filter%state_estimate, filter%covariance_square_root, filter%ut_gamma, sigma_points)
        !! propagate sigma points through system dynamics
        do concurrent (i=1:size(sigma_points,dim=2))
            call state_dynamics_model(sigma_points(:,i), dt, [filter%maximum_velocity])
        end do
        !! calculate sigma point weighted average state estimate
        filter%state_estimate = filter%wx_1*sigma_points(:,1)
        do i=1,size(sigma_points,dim=1)
            filter%state_estimate = filter%state_estimate + filter%w_2_2n1*(sigma_points(:,i+1) + &
                                                                            sigma_points(:,i+1+size(sigma_points,dim=1)))
        end do
        !! center sigma points (sigma_points - state_estimate)
        call broadcast_sub(sigma_points, filter%state_estimate)
        !! form A matrix (amat) as [sqrt(w)*centered_sigma_points, sqrt(Q)]
        amat(:,1) = filter%wc_1*sigma_points(:,1)
        sqrt_wc = sqrt(filter%w_2_2n1)
        do concurrent (i=2:size(sigma_points,dim=2))
            amat(:,i) = sqrt_wc*sigma_points(:,i)
        end do
        call generate_sqrt_process_noise(filter%process_noise, dt, amat(:,size(sigma_points,dim=2)+1:size(amat,dim=2)))
        !! calculate new covariance_square_root
        amat_t = transpose(amat)
        call qr(amat_t, unused_tau)
        call extract_rt(amat_t(1:6,1:6), filter%covariance_square_root)
    end

    pure subroutine extract_rt(r_upper_triangle, rt)
        real(dp), intent(in) :: r_upper_triangle(6,6)
        real(dp), intent(out) :: rt(6,6)
        integer :: i
        rt = 0.0_dp
        do concurrent (i=1:6)
            rt(i:6,i) = r_upper_triangle(i,i:6)
        end do
    end

    pure subroutine generate_sigma_points(state, sqrt_cov, ut_gamma, sigma_points)
        real(dp), intent(in) :: state(:), sqrt_cov(size(state),size(state)), ut_gamma
        real(dp), intent(out) :: sigma_points(size(state),2*size(state)+1)
        integer :: trk_n, i
        trk_n = size(state)
        sigma_points(:,1) = state
        do i=1,trk_n
            sigma_points(:,i+1) = state + ut_gamma*sqrt_cov(:,i)
            sigma_points(:,i+1+trk_n) = state - ut_gamma*sqrt_cov(:,i)
        end do
    end

    pure subroutine dynamics_ballistic(state, dt, opt_data)
        real(dp), intent(inout) :: state(6)
        real(dp), intent(in) :: dt
        real(dp), intent(in), optional :: opt_data(:)
        real(dp) :: vel_max, spd0
        if (present(opt_data)) then
            call debug_error_condition(size(opt_data) /= 1, 'dynamics_ballistic can only accept size 1 optional data')
            vel_max = opt_data(1)
        else
            vel_max = huge(1.0_dp)
        end if
        state(5) = 0.0_dp
        state(6) = -1.0_dp*g
        state(3:4) = state(3:4) + dt*state(5:6)
        spd0 = sqrt(state(3)**2 + state(4)**2)
        if (spd0 > vel_max) state(3:4) = state(3:4)/(spd0/vel_max)
        state(1:2) = state(1:2) + dt*state(3:4)
    end

    pure subroutine generate_sqrt_process_noise(noise, dt, sqrt_q)
        real(dp), intent(in) :: noise, dt
        real(dp), intent(out) :: sqrt_q(:,:)
        real(dp) :: qmat(6,6), dt2, dt3, dt4, dt5, g2, g3, g6, g8, g20
        call debug_error_condition((size(sqrt_q, dim=1) /= 6) .or. (size(sqrt_q, dim=2) /= 6), &
                                   'process noise matrix (sqrt_q) must be 6x6')
        !! compute constants
        dt2 = dt**2
        dt3 = dt**3
        dt4 = dt**4
        dt5 = dt**5
        g2 = noise/2.0_dp
        g3 = noise/3.0_dp
        g6 = noise/6.0_dp
        g8 = noise/8.0_dp
        g20 = noise/20.0_dp
        !! initialize all of Q to 0.0
        qmat = 0.0_dp
        !! fill in Q matrix column by column
        qmat(1,1) = g20*dt5
        qmat(3,1) = g8*dt4
        qmat(5,1) = g6*dt3
        qmat(2,2) = qmat(1,1)
        qmat(4,2) = qmat(3,1)
        qmat(6,2) = qmat(5,1)
        qmat(1,3) = g8*dt4
        qmat(3,3) = g3*dt3
        qmat(5,3) = g2*dt2
        qmat(2,4) = qmat(1,3)
        qmat(4,4) = qmat(3,3)
        qmat(6,4) = qmat(5,3)
        qmat(1,5) = g6*dt3
        qmat(3,5) = g2*dt2
        qmat(5,5) = noise*dt
        qmat(2,6) = qmat(1,5)
        qmat(4,6) = qmat(3,5)
        qmat(6,6) = qmat(5,5)
        !! convert to square root form via Cholesky decomposition
        call chol(qmat, sqrt_q)
    end

end module ukf


program ex_ukf
use, non_intrinsic :: ukf
implicit none
end program ex_ukf
