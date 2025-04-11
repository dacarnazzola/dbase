module ukf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: nmi2ft_dp, deg2rad_dp, eps_dp, rad2deg_dp
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: vector_math, only: vmag, vdot
use, non_intrinsic :: matrix_math, only: chol, qr, forward_substitution, backward_substitution
use, non_intrinsic :: array_utils, only: broadcast_sub
use, non_intrinsic :: random, only: random_normal
use, non_intrinsic :: statistics, only: avg
use, non_intrinsic :: sorting, only: sort
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
        !! measurement update will iterate up to max_iterations to improve state_estimate
        integer :: max_iterations = 1
        !! state_estimate and covariance_square_root have no sensible default values
        real(dp) :: state_estimate(6)
        real(dp) :: covariance_square_root(6,6)
    end type sr_ukf_type

    interface
        pure subroutine dynamics_model(state, dt, opt_data)
        import dp
        implicit none
            real(dp), intent(inout) :: state(6)
            real(dp), intent(in) :: dt
            real(dp), intent(in) :: opt_data(:)
        end
    end interface

    public :: dp, nmi2ft_dp, deg2rad_dp, g, sr_ukf_type, initialize_sr_ukf, dynamics_ballistic, filter_measurement_update, &
              generate_observation, print_status, print_matrix, dynamics_constant_velocity, dynamics_constant_acceleration, &
              dynamics_model, rmse, sort, dump_states, dump_summary

contains

    pure subroutine initialize_sr_ukf(obs, meas, meas_sig, filter, t, noise, max_vel, def_vel, max_acc, def_acc, k, a, b, def_rng, &
                                      max_iterations)
        real(dp), intent(in) :: obs(6), meas(4), meas_sig(4)
        type(sr_ukf_type), intent(out) :: filter
        real(dp), intent(in), optional :: t, noise, max_vel, max_acc, k, a, b, def_rng, def_vel(2), def_acc(2)
        integer, intent(in), optional :: max_iterations
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
        if (present(max_iterations)) filter%max_iterations = max_iterations
        trk_n = 6.0_dp
        ut_lambda = ut_alpha**2 * (trk_n + ut_kappa) - trk_n
        filter%wx_1 = ut_lambda/(trk_n + ut_lambda)
        filter%wc_1 = filter%wx_1 + 1.0_dp - ut_alpha**2 + ut_beta
        call debug_error_condition(filter%wc_1 < 0.0_dp, 'negative weights for covariance matrix not implemented')
        filter%wc_1 = sqrt(filter%wc_1) !! wc_1 primarily used under square root, so compute that here (cross correlation will square this term)
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
        if (.false.) then
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
        else
            jrjt = 0.0_dp
            dim_var = (0.5_dp*default_range)**2/3.0_dp
            jrjt(1,1) = dim_var
            jrjt(2,2) = dim_var
            dim_var = filter%maximum_velocity**2/3.0_dp
            jrjt(3,3) = dim_var
            jrjt(4,4) = dim_var
        end if
        dim_var = filter%maximum_acceleration**2/3.0_dp
        jrjt(5,5) = jrjt(5,5) + dim_var
        jrjt(6,6) = jrjt(6,6) + dim_var
        call chol(jrjt, filter%covariance_square_root)
        call iterate_p0(obs, meas, meas_sig, filter)
    end

    pure subroutine iterate_p0(obs, meas, meas_sig, filter)
        real(dp), intent(in) :: obs(6), meas(4), meas_sig(4)
        type(sr_ukf_type), intent(inout) :: filter
        real(dp) :: sigma_points(6,13), sigma_points_meas(4,13), pred_meas(4), innovation(4), amat(4,19), amat_t(19,4), sqrt_wc, &
                    square_root_measurement_covariance(4,4), cross_correlation(6,4), p_xz_transpose(4,6), temp(4,6), &
                    s_z_transpose(4,4), k_t(4,6), kalman_gain(6,4), p(6,6), pzz_plus_r(4,4), y(4), chi2_0, chi2_now, chi2_diff
        integer :: i, meas_dim, meas_ii(4), meas_index, iteration
        !! return early if measurement provides no information
        if (count(meas_sig > 0.0_dp) < 1) then
            return
        else
        !! otherwise, collect valid measurement dimensions
            meas_dim = 0
            meas_ii = -1
            do i=1,4
                if (meas_sig(i) > 0.0_dp) then
                    meas_dim = meas_dim + 1
                    meas_ii(meas_dim) = i
                end if
            end do
        end if
        chi2_0 = huge(1.0_dp)
        do iteration=1,filter%max_iterations/2
            !! generate sigma points
            call generate_sigma_points(filter%state_estimate, filter%covariance_square_root, filter%ut_gamma, sigma_points)
            !! convert sigma_points to measurement space
            do concurrent (i=1:13)
                call cart2pol(obs, sigma_points(:,i), sigma_points_meas(:,i))
            end do
            !! center original state estimate dimension sigma_points for use later
            call broadcast_sub(sigma_points, filter%state_estimate)
            !! calculate sigma point weighted average predicted measurement
            pred_meas = filter%wx_1*sigma_points_meas(:,1)
            do i=1,6
                pred_meas = pred_meas + filter%w_2_2n1*(sigma_points_meas(:,i+1) + sigma_points_meas(:,i+7))
            end do
            !! calculate innovation
            innovation = 0.0_dp
            do concurrent (i=1:meas_dim)
                meas_index = meas_ii(i)
                if (meas_index == 2) then !! meas_dim(meas_ii(i)) is ANGLE, need to properly wrap -pi/+pi
                    innovation(meas_index) = atan2(sin(meas(meas_index) - pred_meas(meas_index)), &
                                                   cos(meas(meas_index) - pred_meas(meas_index)))
                else
                    innovation(meas_index) = meas(meas_index) - pred_meas(meas_index)
                end if
            end do
            !! center measurement space sigma points
            do concurrent (i=1:13)
                sigma_points_meas(:,i) = sigma_points_meas(:,i) - pred_meas
                if (meas_sig(2) > 0.0_dp) sigma_points_meas(2,i) = atan2(sin(sigma_points_meas(2,i)), cos(sigma_points_meas(2,i)))
            end do
            !! form amat as [sqrt(w)*centered_measurement_sigma_points, sqrt(R)] where R is measurement covariance (independent)
            amat = 0.0_dp
            amat(1:meas_dim,1) = filter%wc_1*sigma_points_meas(meas_ii(1:meas_dim),1)
            sqrt_wc = sqrt(filter%w_2_2n1)
            do concurrent (i=2:13)
                amat(1:meas_dim,i) = sqrt_wc*sigma_points_meas(meas_ii(1:meas_dim),i)
            end do
            do concurrent (i=1:meas_dim)
                amat(i,13+i) = meas_sig(meas_ii(i))
            end do
            !! calculate square root of measurement covariance
            amat_t = transpose(amat)
            call qr(amat_t(1:13+meas_dim,1:meas_dim))
            square_root_measurement_covariance = 0.0_dp
            call extract_rt(amat_t(1:meas_dim,1:meas_dim), square_root_measurement_covariance(1:meas_dim,1:meas_dim))
            !! calculate cross correlation
            cross_correlation = 0.0_dp
            cross_correlation(:,1:meas_dim) = filter%wc_1**2*matmul(sigma_points(:,1:1), &
                                                                    transpose(sigma_points_meas(meas_ii(1:meas_dim),1:1)))
            do i=2,13
                cross_correlation(:,1:meas_dim) = cross_correlation(:,1:meas_dim) + &
                                                  filter%w_2_2n1*matmul(sigma_points(:,i:i), &
                                                                        transpose(sigma_points_meas(meas_ii(1:meas_dim),i:i)))
            end do
            !! calculate Kalman gain
            p_xz_transpose(1:meas_dim,:) = transpose(cross_correlation(:,1:meas_dim))
            temp = 0.0_dp !! (meas_dim x 6) to be filled in
            do concurrent (i=1:6)
                call forward_substitution(L=square_root_measurement_covariance(1:meas_dim,1:meas_dim), & !! sqrt measurement cov
                                          b=p_xz_transpose(1:meas_dim,i), &                              !! cross covariance for each state dimension
                                          x=temp(1:meas_dim,i))                                          !! temporary solution vector
            end do
            s_z_transpose(1:meas_dim,1:meas_dim) = transpose(square_root_measurement_covariance(1:meas_dim,1:meas_dim))
            k_t = 0.0_dp !! transpose of Kalman gain, solved for using temp from above
            do concurrent (i=1:6)
                call backward_substitution(U=s_z_transpose(1:meas_dim,1:meas_dim), & !! transpose(Sz)
                                           b=temp(1:meas_dim,i), &                   !! columns of temp, solved for above
                                           x=k_t(1:meas_dim,i))                      !! columns of transpose(Kalman gain) --> rows of Kalman gain
            end do
            kalman_gain(:,1:meas_dim) = transpose(k_t(1:meas_dim,:))
            !! update filter%covariance_square_root
            call reform_cov(filter%covariance_square_root, p)
            call reform_cov(square_root_measurement_covariance(1:meas_dim,1:meas_dim), pzz_plus_r(1:meas_dim,1:meas_dim))
            p = p - matmul(kalman_gain(:,1:meas_dim), matmul(pzz_plus_r(1:meas_dim,1:meas_dim), k_t(1:meas_dim,:)))
            call chol(p, filter%covariance_square_root)
            !! solve sqrt(Pzz)*y = v for y
            call forward_substitution(L=square_root_measurement_covariance(1:meas_dim,1:meas_dim), &
                                      b=innovation(meas_ii(1:meas_dim)), &
                                      x=y(1:meas_dim))
            !! compute ||y||**2 = Chi**2 - monitor chi2 until it stops changing and is less than threshold
            chi2_now = vdot(y(1:meas_dim), y(1:meas_dim))
            chi2_diff = abs(chi2_0 - chi2_now)
            if (chi2_diff < 0.01_dp) then
                select case (meas_dim)
                    case (1)
                        if (chi2_now < 3.841_dp) exit
                    case (2)
                        if (chi2_now < 5.991_dp) exit
                    case (3)
                        if (chi2_now < 7.815_dp) exit
                    case (4)
                        if (chi2_now < 9.488_dp) exit
                    case default
                        error stop 'only implemented Chi**2 critical values for up to 4 degrees of freedom'
                end select
            else
                chi2_0 = chi2_now
            end if
        end do
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
        do i=1,4
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

    impure subroutine filter_measurement_update(obs, meas, meas_sig, filter, t, state_dynamics_model, opt_data)
        real(dp), intent(in) :: obs(6), meas(4), meas_sig(4)
        type(sr_ukf_type), intent(inout) :: filter
        real(dp), intent(in) :: t
        procedure(dynamics_model) :: state_dynamics_model
        real(dp), intent(in) :: opt_data(:)
        real(dp) :: sigma_points(6,13), sigma_points_meas(4,13), pred_meas(4), innovation(4), amat(4,19), amat_t(19,4), sqrt_wc, &
                    square_root_measurement_covariance(4,4), cross_correlation(6,4), p_xz_transpose(4,6), temp(4,6), &
                    s_z_transpose(4,4), k_t(4,6), kalman_gain(6,4), p(6,6), pzz_plus_r(4,4), y(4), chi2_0, chi2_now, chi2_diff
        integer :: i, meas_dim, meas_ii(4), meas_index, iteration, r_alpha
        !! update filter to current time
        call filter_time_update(filter, t, state_dynamics_model, opt_data)
        !! return early if measurement provides no information
        if (count(meas_sig > 0.0_dp) < 1) then
            return
        else
        !! otherwise, collect valid measurement dimensions
            meas_dim = 0
            meas_ii = -1
            do i=1,4
                if (meas_sig(i) > 0.0_dp) then
                    meas_dim = meas_dim + 1
                    meas_ii(meas_dim) = i
                end if
            end do
        end if
        chi2_0 = huge(1.0_dp)
        do iteration=1,filter%max_iterations
            !! generate sigma points
            call generate_sigma_points(filter%state_estimate, filter%covariance_square_root, filter%ut_gamma, sigma_points)
            !! convert sigma_points to measurement space
            do concurrent (i=1:13)
                call cart2pol(obs, sigma_points(:,i), sigma_points_meas(:,i))
            end do
            !! center original state estimate dimension sigma_points for use later
            call broadcast_sub(sigma_points, filter%state_estimate)
            !! calculate sigma point weighted average predicted measurement
            pred_meas = filter%wx_1*sigma_points_meas(:,1)
            do i=1,6
                pred_meas = pred_meas + filter%w_2_2n1*(sigma_points_meas(:,i+1) + sigma_points_meas(:,i+7))
            end do
            !! calculate innovation
            innovation = 0.0_dp
            do concurrent (i=1:meas_dim)
                meas_index = meas_ii(i)
                if (meas_index == 2) then !! meas_dim(meas_ii(i)) is ANGLE, need to properly wrap -pi/+pi
                    innovation(meas_index) = atan2(sin(meas(meas_index) - pred_meas(meas_index)), &
                                                   cos(meas(meas_index) - pred_meas(meas_index)))
                else
                    innovation(meas_index) = meas(meas_index) - pred_meas(meas_index)
                end if
            end do
            !! center measurement space sigma points
            do concurrent (i=1:13)
                sigma_points_meas(:,i) = sigma_points_meas(:,i) - pred_meas
                if (meas_sig(2) > 0.0_dp) sigma_points_meas(2,i) = atan2(sin(sigma_points_meas(2,i)), cos(sigma_points_meas(2,i)))
            end do
            !! form amat as [sqrt(w)*centered_measurement_sigma_points, sqrt(R)] where R is measurement covariance (independent)
            amat = 0.0_dp
            amat(1:meas_dim,1) = filter%wc_1*sigma_points_meas(meas_ii(1:meas_dim),1)
            sqrt_wc = sqrt(filter%w_2_2n1)
            do concurrent (i=2:13)
                amat(1:meas_dim,i) = sqrt_wc*sigma_points_meas(meas_ii(1:meas_dim),i)
            end do
            r_alpha = 1.0_dp + (iteration - 1)*0.1_dp
            do concurrent (i=1:meas_dim)
                amat(i,13+i) = meas_sig(meas_ii(i))*r_alpha
            end do
            !! calculate square root of measurement covariance
            amat_t = transpose(amat)
            call qr(amat_t(1:13+meas_dim,1:meas_dim))
            square_root_measurement_covariance = 0.0_dp
            call extract_rt(amat_t(1:meas_dim,1:meas_dim), square_root_measurement_covariance(1:meas_dim,1:meas_dim))
            !! calculate cross correlation
            cross_correlation = 0.0_dp
            cross_correlation(:,1:meas_dim) = filter%wc_1**2*matmul(sigma_points(:,1:1), &
                                                                    transpose(sigma_points_meas(meas_ii(1:meas_dim),1:1)))
            do i=2,13
                cross_correlation(:,1:meas_dim) = cross_correlation(:,1:meas_dim) + &
                                                  filter%w_2_2n1*matmul(sigma_points(:,i:i), &
                                                                        transpose(sigma_points_meas(meas_ii(1:meas_dim),i:i)))
            end do
            !! calculate Kalman gain
            p_xz_transpose(1:meas_dim,:) = transpose(cross_correlation(:,1:meas_dim))
            temp = 0.0_dp !! (meas_dim x 6) to be filled in
            do concurrent (i=1:6)
                call forward_substitution(L=square_root_measurement_covariance(1:meas_dim,1:meas_dim), & !! sqrt measurement cov
                                          b=p_xz_transpose(1:meas_dim,i), &                              !! cross covariance for each state dimension
                                          x=temp(1:meas_dim,i))                                          !! temporary solution vector
            end do
            s_z_transpose(1:meas_dim,1:meas_dim) = transpose(square_root_measurement_covariance(1:meas_dim,1:meas_dim))
            k_t = 0.0_dp !! transpose of Kalman gain, solved for using temp from above
            do concurrent (i=1:6)
                call backward_substitution(U=s_z_transpose(1:meas_dim,1:meas_dim), & !! transpose(Sz)
                                           b=temp(1:meas_dim,i), &                   !! columns of temp, solved for above
                                           x=k_t(1:meas_dim,i))                      !! columns of transpose(Kalman gain) --> rows of Kalman gain
            end do
            kalman_gain(:,1:meas_dim) = transpose(k_t(1:meas_dim,:))
            !! update filter%state_estimate
            filter%state_estimate = filter%state_estimate + matmul(kalman_gain(:,1:meas_dim), innovation(meas_ii(1:meas_dim)))
            !! solve sqrt(Pzz)*y = v for y
            call forward_substitution(L=square_root_measurement_covariance(1:meas_dim,1:meas_dim), &
                                      b=innovation(meas_ii(1:meas_dim)), &
                                      x=y(1:meas_dim))
            !! compute ||y||**2 = Chi**2 - monitor chi2 until it stops changing and is less than threshold
            chi2_now = vdot(y(1:meas_dim), y(1:meas_dim))
            chi2_diff = abs(chi2_0 - chi2_now)
            if (chi2_diff < 0.01_dp) then
                select case (meas_dim)
                    case (1)
                        if (chi2_now < 3.841_dp) exit
                    case (2)
                        if (chi2_now < 5.991_dp) exit
                    case (3)
                        if (chi2_now < 7.815_dp) exit
                    case (4)
                        if (chi2_now < 9.488_dp) exit
                    case default
                        error stop 'only implemented Chi**2 critical values for up to 4 degrees of freedom'
                end select
            else
                chi2_0 = chi2_now
            end if
        end do
        !! update filter%covariance_square_root
        call reform_cov(filter%covariance_square_root, p)
        call reform_cov(square_root_measurement_covariance(1:meas_dim,1:meas_dim), pzz_plus_r(1:meas_dim,1:meas_dim))
        p = p - matmul(kalman_gain(:,1:meas_dim), matmul(pzz_plus_r(1:meas_dim,1:meas_dim), k_t(1:meas_dim,:)))
        call chol(p, filter%covariance_square_root)
    end

    impure subroutine print_matrix(msg, mat)
        character(len=*), intent(in) :: msg
        real(dp), intent(in) :: mat(:,:)
        integer :: i
        write(*,'(a)') msg
        do i=1,size(mat,dim=1)
           write(*,'(*(e13.6," "))') mat(i,:) 
        end do
    end

    pure subroutine reform_cov(sqrt_cov, cov)
        real(dp), intent(in) :: sqrt_cov(:,:)
        real(dp), intent(out) :: cov(:,:)
        integer :: n, i, j
        call debug_error_condition(size(sqrt_cov,dim=1) /= size(sqrt_cov,dim=2), 'sqrt_cov must be square')
        call debug_error_condition(size(cov,dim=1) /= size(sqrt_cov,dim=1), 'cov dim 1 must match sqrt_cov dim 1')
        call debug_error_condition(size(cov,dim=2) /= size(sqrt_cov,dim=2), 'cov dim 1 must match sqrt_cov dim 1')
        cov = matmul(sqrt_cov, transpose(sqrt_cov))
        n = size(cov, dim=1)
        do i=1,n
            do j=i+1,n
                cov(j,i) = cov(i,j)
            end do
            cov(i,i) = cov(i,i) + eps_dp
        end do
    end

    pure subroutine filter_time_update(filter, t, state_dynamics_model, opt_data)
        type(sr_ukf_type), intent(inout) :: filter
        real(dp), intent(in) :: t
        procedure(dynamics_model) :: state_dynamics_model
        real(dp), intent(in) :: opt_data(:)
        real(dp) :: dt, sigma_points(6,13), amat(6,19), sqrt_wc, amat_t(19,6), r(6,6)
        integer :: i
        call debug_error_condition(filter%state_estimate_time < 0.0_dp, 'filter is not initialized')
        dt = t - filter%state_estimate_time
        call debug_error_condition(dt < 0.0_dp, 'negative dt time update not implemented')
        filter%state_estimate_time = t
        if (nearly(dt, 0.0_dp)) return !! return early if dt is approximately 0.0
        !! generate sigma points
        call generate_sigma_points(filter%state_estimate, filter%covariance_square_root, filter%ut_gamma, sigma_points)
        !! propagate sigma points through system dynamics
        do concurrent (i=1:13)
            call state_dynamics_model(sigma_points(:,i), dt, opt_data)
        end do
        !! calculate sigma point weighted average state estimate
        filter%state_estimate = filter%wx_1*sigma_points(:,1)
        do i=1,6
            filter%state_estimate = filter%state_estimate + filter%w_2_2n1*(sigma_points(:,i+1) + sigma_points(:,i+7))
        end do
        !! center sigma points (sigma_points - state_estimate)
        call broadcast_sub(sigma_points, filter%state_estimate)
        !! form A matrix (amat) as [sqrt(w)*centered_sigma_points, sqrt(Q)]
        amat(:,1) = filter%wc_1*sigma_points(:,1)
        sqrt_wc = sqrt(filter%w_2_2n1)
        do concurrent (i=2:13)
            amat(:,i) = sqrt_wc*sigma_points(:,i)
        end do
        call generate_sqrt_process_noise(filter%process_noise, dt, amat(:,14:19))
        !! calculate new covariance_square_root
        amat_t = transpose(amat)
        call qr(amat_t)
        r = amat_t(1:6,1:6)
        call extract_rt(r, filter%covariance_square_root)
    end

    pure subroutine extract_rt(r, r_upper_triangle_transpose)
        real(dp), intent(in) :: r(:,:)
        real(dp), intent(out) :: r_upper_triangle_transpose(size(r,dim=1),size(r,dim=2))
        integer :: n, i
        call debug_error_condition(size(r,dim=1) /= size(r,dim=2), 'R matrix must be square')
        n = size(r, dim=1)
        r_upper_triangle_transpose = 0.0_dp
        do concurrent (i=1:n)
            r_upper_triangle_transpose(i:n,i) = r(i,i:n)
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

    pure subroutine dynamics_constant_velocity(state, dt, opt_data)
        real(dp), intent(inout) :: state(6)
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: opt_data(:)
        call debug_error_condition(size(opt_data) /= 0, 'dynamics_constant_velocity does not accept optional data')
        state(1:2) = state(1:2) + dt*state(3:4)
    end

    pure subroutine dynamics_constant_acceleration(state, dt, opt_data)
        real(dp), intent(inout) :: state(6)
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: opt_data(:)
        real(dp) :: vel_max, spd0, acc_max, acc0
        call debug_error_condition((size(opt_data) /= 0) .and. (size(opt_data) /= 2), &
                                   'dynamics_constant_acceleration can only accept size 2 optional data')
        if (size(opt_data) == 2) then
            vel_max = opt_data(1)
            acc_max = opt_data(2)
        else
            vel_max = huge(1.0_dp)
            acc_max = huge(1.0_dp)
        end if
        acc0 = sqrt(state(5)**2 + state(6)**2)
        if (acc0 > acc_max) state(5:6) = state(5:6)/(acc0/acc_max)
        state(3:4) = state(3:4) + dt*state(5:6)
        spd0 = sqrt(state(3)**2 + state(4)**2)
        if (spd0 > vel_max) state(3:4) = state(3:4)/(spd0/vel_max)
        state(1:2) = state(1:2) + dt*state(3:4)
    end

    pure subroutine dynamics_ballistic(state, dt, opt_data)
        real(dp), intent(inout) :: state(6)
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: opt_data(:)
        real(dp) :: vel_max, spd0
        call debug_error_condition((size(opt_data) /= 0) .and. (size(opt_data) /= 1), &
                                   'dynamics_ballistic can only accept size 1 optional data')
        if (size(opt_data) == 1) then
            vel_max = opt_data(1)
        else
            vel_max = huge(1.0_dp)
        end if
        state(4) = state(4) - dt*g
        spd0 = sqrt(state(3)**2 + state(4)**2)
        if (spd0 > vel_max) state(3:4) = state(3:4)/(spd0/vel_max)
        state(1:2) = state(1:2) + dt*state(3:4)
    end

    pure subroutine generate_sqrt_process_noise(noise, dt, sqrt_q)
        real(dp), intent(in) :: noise, dt
        real(dp), intent(out) :: sqrt_q(6,6)
        real(dp) :: qmat(6,6), dt2, dt3, dt4, dt5, g2, g3, g6, g8, g20
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

    impure subroutine generate_observation(obs, tgt, meas, meas_sig)
        real(dp), intent(in) :: obs(6), tgt(6)
        real(dp), intent(out) :: meas(4)
        real(dp), intent(in), optional :: meas_sig(4)
        real(dp) :: z(1)
        integer :: i
        call cart2pol(obs, tgt, meas)
        if (present(meas_sig)) then
            do i=1,4
                if (meas_sig(i) > 0.0_dp) then
                    call random_normal(z, 0.0_dp, meas_sig(i))
                    meas(i) = meas(i) + z(1)
                else
                    meas(i) = huge(1.0_dp)
                end if
            end do
        end if
    end

    impure subroutine print_status(sim_t, obs, tgt, filter, extra)
        real(dp), intent(in) :: sim_t, obs(6), tgt(6)
        type(sr_ukf_type), intent(in) :: filter
        logical, intent(in), optional :: extra
        real(dp) :: p(6,6), sigma_points(6,13), sigma_points_meas(4,13), pred_meas(4), amat(4,13), amat_t(13,4), &
                    square_root_measurement_covariance(4,4), sqrt_wc, pz(4,4)
        integer :: i
        write(*,'(a)') repeat('=', 32)
        write(*,'(a,2(f0.4,a))') 'Simulation Time: ',sim_t,' seconds, Filter Time: ',filter%state_estimate_time,' seconds'
        write(*,'(6(a,e13.6))') '       Observer State :: x: ',obs(1),', y: ',obs(2), &
                               ', vx: ',obs(3),', vy: ',obs(4),', ax: ',obs(5),', ay: ',obs(6)
        write(*,'(6(a,e13.6))') '         Target State :: x: ',tgt(1),', y: ',tgt(2), &
                               ', vx: ',tgt(3),', vy: ',tgt(4),', ax: ',tgt(5),', ay: ',tgt(6)
        write(*,'(6(a,e13.6))') 'Filter State Estimate :: x: ',filter%state_estimate(1),', y: ',filter%state_estimate(2), &
                               ', vx: ',filter%state_estimate(3),', vy: ',filter%state_estimate(4), &
                               ', ax: ',filter%state_estimate(5),', ay: ',filter%state_estimate(6)
        call reform_cov(filter%covariance_square_root, p)
        write(*,'(6(a,e13.6))') '         Track Sigmas :: x: ',sqrt(p(1,1)),', y: ',sqrt(p(2,2)), &
                                ', vx: ',sqrt(p(3,3)),', vy: ',sqrt(p(4,4)),', ax: ',sqrt(p(5,5)),', ay: ',sqrt(p(6,6))
        call generate_sigma_points(filter%state_estimate, filter%covariance_square_root, filter%ut_gamma, sigma_points)
        do concurrent (i=1:13)
            call cart2pol(obs, sigma_points(:,i), sigma_points_meas(:,i))
        end do
        !! calculate sigma point weighted average predicted measurement
        pred_meas = filter%wx_1*sigma_points_meas(:,1)
        do i=1,6
            pred_meas = pred_meas + filter%w_2_2n1*(sigma_points_meas(:,i+1) + sigma_points_meas(:,i+7))
        end do
        !! center measurement space sigma points
        do concurrent (i=1:13)
            sigma_points_meas(:,i) = sigma_points_meas(:,i) - pred_meas
            sigma_points_meas(2,i) = atan2(sin(sigma_points_meas(2,i)), cos(sigma_points_meas(2,i)))
        end do
        !! form amat as [sqrt(w)*centered_measurement_sigma_points]
        amat = 0.0_dp
        amat(:,1) = filter%wc_1*sigma_points_meas(:,1)
        sqrt_wc = sqrt(filter%w_2_2n1)
        do concurrent (i=2:13)
            amat(:,i) = sqrt_wc*sigma_points_meas(:,i)
        end do
        !! calculate square root of measurement covariance
        amat_t = transpose(amat)
        call qr(amat_t)
        square_root_measurement_covariance = 0.0_dp
        call extract_rt(amat_t(1:4,:), square_root_measurement_covariance)
        call reform_cov(square_root_measurement_covariance, pz)
        write(*,'(4(a,f0.2))') '   Polar Track Sigmas :: range: ',sqrt(pz(1,1)),', angle: ',sqrt(pz(2,2))*rad2deg_dp, &
                               ', range-rate: ',sqrt(pz(3,3)),', angle-rate: ',sqrt(pz(4,4))*rad2deg_dp
        if (present(extra)) then
            if (extra) then
                do i=1,4
                    write(*,'(a,i0,a,4(e13.6," "))') 'Measurement-Space Covariance (row ',i,'): ',pz(i,:)
                end do
                do i=1,6
                    write(*,'(a,i0,a,6(e13.6," "))') 'Full Track Covariance (row ',i,'): ',p(i,:)
                end do
            end if
        end if
        write(*,'(a,/)') repeat('=', 32)
    end

    pure function rmse(predicted, observed) result(val)
        real(dp), intent(in) :: predicted(:), observed
        real(dp) :: val
        real(dp) :: diff2(size(predicted))
        diff2 = (predicted - observed)**2
        val = sqrt(avg(diff2))
    end

    impure subroutine dump_states(fid, model_ii, t, obs, tgt, trk)
        integer, intent(in) :: fid, model_ii
        real(dp), intent(in) :: t, obs(6), tgt(6), trk(6)
        real(dp) :: trk_err(6), tru_pol(4), trk_pol(4), trk_err_pol(4), pos_vel_acc_err(3)
        trk_err = abs(trk - tgt)
        pos_vel_acc_err(1) = sqrt(trk_err(1)**2 + trk_err(2))
        pos_vel_acc_err(2) = sqrt(trk_err(3)**2 + trk_err(4))
        pos_vel_acc_err(3) = sqrt(trk_err(5)**2 + trk_err(6))
        call cart2pol(obs, tgt, tru_pol)
        call cart2pol(obs, trk, trk_pol)
        trk_err_pol = abs(trk_pol - tru_pol)
        write(unit=fid, fmt='(i0,40(",",e13.6))') model_ii, t, obs, tgt, trk, trk_err, pos_vel_acc_err, &
                                                  tru_pol, trk_pol, trk_err_pol
    end

    impure subroutine dump_summary(fid, reps, model_ii, filter_noise, init_sig_scale, init_sig, meas_sig, se_err, se_err_pol, &
                                   ut_alpha, ut_lambda, ut_kappa, max_iterations)
        integer, intent(in) :: fid, reps, model_ii, max_iterations
        real(dp), intent(in) :: filter_noise, init_sig_scale, init_sig(4), meas_sig(4), se_err(:,:), se_err_pol(:,:), ut_alpha, &
                                ut_lambda, ut_kappa
        real(dp) :: rmse_state(size(se_err,dim=1)), rmse_pol(size(se_err_pol,dim=1))
        integer :: i
        call debug_error_condition(size(se_err,dim=1) /= 9, 'se_err should have 9 rows')
        call debug_error_condition(size(se_err_pol,dim=1) /= 4, 'se_err_pol should have 4 rows')
        call debug_error_condition(size(se_err,dim=2) /= size(se_err_pol,dim=2), 'se_err and se_err_pol must have same columns')
        do concurrent (i=1:size(rmse_state))
            rmse_state(i) = rmse(se_err(i,:), 0.0_dp)
        end do
        do concurrent (i=1:size(rmse_pol))
            rmse_pol(i) = rmse(se_err_pol(i,:), 0.0_dp)
        end do
        write(unit=fid, fmt='(i0,",",i0,26(",",e13.6),",",i0)') reps, model_ii, filter_noise, init_sig_scale, init_sig, &
                                                                meas_sig, rmse_state, rmse_pol, ut_alpha, ut_lambda, ut_kappa, &
                                                                max_iterations
    end

end module ukf


program ex_sr_ukf
use, non_intrinsic :: ukf
implicit none

    logical, parameter :: debug(*) = [.false., & !! 1, per-observation
                                      .false., & !! 2, per-trial
                                      .false., & !! 3, all-trial summary
                                      .false., & !! 4, per-observation output to file
                                      .false., & !! 5, per-trial output to file
                                      .true.] !! 6, all-trial summary output to file
    real(dp), parameter :: dt = 1.0_dp
    integer, parameter :: max_trials = 1024
    real(dp), parameter :: meas_sig1_list(*) = [-1.0_dp, 1.0_dp, 10.0_dp, 100.0_dp, nmi2ft_dp, 10.0_dp*nmi2ft_dp]
    real(dp), parameter :: meas_sig2_list(*) = [-1.0_dp, 0.01_dp*deg2rad_dp, 0.1_dp*deg2rad_dp, deg2rad_dp, 5.0_dp*deg2rad_dp]
    real(dp), parameter :: meas_sig3_list(*) = [-1.0_dp, 0.1_dp, 10.0_dp, 100.0_dp, 200.0_dp, 1000.0_dp]
    real(dp), parameter :: meas_sig4_list(*) = [-1.0_dp, 0.001_dp*deg2rad_dp, 0.01_dp*deg2rad_dp, 0.1_dp*deg2rad_dp]
    real(dp), parameter :: init_scale_list(*) = [1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, 3.5_dp, 4.0_dp, 4.5_dp, 5.0_dp, 5.5_dp, &
                                                 6.0_dp, 6.5_dp, 7.0_dp, 7.5_dp, 8.0_dp, 8.5_dp, 9.0_dp]
    real(dp), parameter :: noise_list(*) = [0.01_dp/g, 0.01_dp, 0.1_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 9.0_dp]*g
    real(dp), parameter :: ut_alpha_list(*) = [0.001_dp, 0.002_dp, 0.003_dp, 0.004_dp, 0.005_dp, 0.006_dp, 0.007_dp, 0.008_dp, &
                                               0.009_dp, 0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp, 0.06_dp, 0.07_dp, 0.08_dp, &
                                               0.09_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp, 0.7_dp, 0.8_dp, 0.9_dp, &
                                               1.0_dp]
    real(dp), parameter :: ut_lambda_list(*) = [0.0_dp, 0.001_dp, 0.01_dp, 0.1_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
                                                7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp]
    integer, parameter :: max_iterations_list(*) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    type(sr_ukf_type) :: filter
    real(dp) :: obs(6), tgt(6), meas(4), meas_sig(4), init_sig(4), t, no_opt_data(0), se_err(9,max_trials), tgt0(6), tgt_pol(4), &
                se_err_pol(4,max_trials), trk_pol(4), noise, init_sig_scale, ut_alpha, ut_lambda, ut_kappa
    procedure(dynamics_model), pointer :: state_dynamics_model
    character(len=:), allocatable :: state_dynamics_model_name
    real(dp), allocatable :: opt_data(:)
    integer :: state_model_ii, trial_ii, i, in_run_stats_fid, end_run_stats_fid, meas_sig1_ii, meas_sig2_ii, meas_sig3_ii, &
               meas_sig4_ii, sr_ukf_summary_fid, noise_ii, init_sig1_ii, init_sig2_ii, init_sig3_ii, init_sig4_ii, init_scale_ii, &
               ut_alpha_ii, ut_lambda_ii, max_iterations_ii, max_iterations

    !! observer initial conditions
    obs = 0.0_dp

    if (debug(4)) then
        open(newunit=in_run_stats_fid, file='/valinor/in-run-stats.csv', action='write')
        write(unit=in_run_stats_fid, fmt='(a)') 'state_model_ii,t,obs_x,obs_y,obs_vx,obs_vy,obs_ax,obs_ay,tgt_x,tgt_y,tgt_vx,'// &
                                                'tgt_vy,tgt_ax,tgt_ay,trk_x,trk_y,trk_vx,trk_vy,trk_ax,trk_ay,trk_err_x,'// &
                                                'trk_err_y,trk_err_vx,trk_err_vy,trk_err_ax,trk_err_ay,trk_pos_err,trk_vel_err,'// &
                                                'trk_acc_err,rng,ang,rng_rt,ang_rt,trk_rng,trk_ang,trk_rng_rt,trk_ang_rt,'// &
                                                'trk_rng_err,trk_ang_err,trk_rng_rt_err,trk_ang_rt_err'
    end if
    if (debug(5)) then
        open(newunit=end_run_stats_fid, file='/valinor/end-run-stats.csv', action='write')
        write(unit=end_run_stats_fid, fmt='(a)') 'state_model_ii,t,obs_x,obs_y,obs_vx,obs_vy,obs_ax,obs_ay,tgt_x,tgt_y,tgt_vx,'// &
                                                 'tgt_vy,tgt_ax,tgt_ay,trk_x,trk_y,trk_vx,trk_vy,trk_ax,trk_ay,trk_err_x,'// &
                                                 'trk_err_y,trk_err_vx,trk_err_vy,trk_err_ax,trk_err_ay,trk_pos_err,'// &
                                                 'trk_vel_err,trk_acc_err,rng,ang,rng_rt,ang_rt,trk_rng,trk_ang,trk_rng_rt,'// &
                                                 'trk_ang_rt,trk_rng_err,trk_ang_err,trk_rng_rt_err,trk_ang_rt_err'
    end if
    if (debug(6)) then
        open(newunit=sr_ukf_summary_fid, file='/valinor/sr-ukf-summary.csv', action='write')
        write(unit=sr_ukf_summary_fid, fmt='(a)') 'reps,state_model_ii,process_noise,init_sig_scale,init_sig_rng,init_sig_ang,'// &
                                                  'init_sig_rng_rt,init_sig_ang_rt,sig_rng,sig_ang,sig_rng_rt,sig_ang_rt,'// &
                                                  'rmse_x,rmse_y,rmse_vx,rmse_vy,rmse_ax,rmse_ay,rmse_pos,rmse_vel,rmse_acc,'// &
                                                  'rmse_rng,rmse_ang,rmse_rng_rt,rmse_ang_rt,ut_alpha,ut_lambda,ut_kappa,'// &
                                                  'max_iterations'
    end if

    do init_sig1_ii=6,6 ! 1,size(meas_sig1_list)
        init_sig(1) = meas_sig1_list(init_sig1_ii)
    do init_sig2_ii=5,5 ! 1,size(meas_sig2_list)
        init_sig(2) = meas_sig2_list(init_sig2_ii)
        if (init_sig(2) <= 0.0_dp) cycle
    do init_sig3_ii=6,6 ! 1,size(meas_sig3_list)
        init_sig(3) = meas_sig3_list(init_sig3_ii)
    do init_sig4_ii=1,1 ! 1,size(meas_sig4_list)
        init_sig(4) = meas_sig4_list(init_sig4_ii)
    do init_scale_ii=1,1 ! 1,size(init_scale_list)
        init_sig_scale = init_scale_list(init_scale_ii)
    do meas_sig1_ii=6,6 ! 2,size(meas_sig1_list)
        meas_sig(1) = meas_sig1_list(meas_sig1_ii)
    do meas_sig2_ii=5,5 ! 2,size(meas_sig2_list)
        meas_sig(2) = meas_sig2_list(meas_sig2_ii)
    do meas_sig3_ii=6,6 ! 2,size(meas_sig3_list)
        meas_sig(3) = meas_sig3_list(meas_sig3_ii)
    do meas_sig4_ii=1,1 ! 1,size(meas_sig4_list)
        meas_sig(4) = meas_sig4_list(meas_sig4_ii)
    do noise_ii=4,4 ! 1,size(noise_list)
        noise = noise_list(noise_ii)
    do ut_alpha_ii=28,28 ! 1,size(ut_alpha_list)
        ut_alpha = ut_alpha_list(ut_alpha_ii)
    do ut_lambda_ii=1,1 ! 1,size(ut_lambda_list)
        ut_lambda = ut_lambda_list(ut_lambda_ii)
        ut_kappa = (6.0_dp + ut_lambda)/ut_alpha**2 - 6.0_dp
    do max_iterations_ii=1,1 ! 1,size(max_iterations_list)
        max_iterations = max_iterations_list(max_iterations_ii)

    do state_model_ii=3,3
        if (allocated(opt_data)) deallocate(opt_data)
        select case (state_model_ii)
            case (1)
                state_dynamics_model => dynamics_constant_velocity
                state_dynamics_model_name = 'Constant-Velocity'
                allocate(opt_data(0))
            case (2)
                state_dynamics_model => dynamics_constant_acceleration
                state_dynamics_model_name = 'Constant-Acceleration'
                opt_data = [filter%maximum_velocity, filter%maximum_acceleration]
            case (3)
                state_dynamics_model => dynamics_ballistic
                state_dynamics_model_name = 'Ballistic'
                opt_data = [filter%maximum_velocity]
            case default
                error stop 'no state_dynamics_model defined'
        end select

    !$omp parallel do default(firstprivate) shared(se_err, se_err_pol)
    do trial_ii=1,max_trials

    !! target initial conditions
    tgt = [200.0_dp*nmi2ft_dp, 1.0_dp, &
           -3000.0_dp*cos(45.0_dp*deg2rad_dp), 3000.0_dp*sin(45.0_dp*deg2rad_dp), &
           0.0_dp, -1.0_dp*g]
    tgt0 = tgt
    !! generate initial measurement and initialize Square Root Unscented Kalman Filter
    t = 0.0_dp
!    call generate_observation(obs, tgt, meas, init_sig)
!    call initialize_sr_ukf(obs, meas, init_sig_scale*init_sig, filter, noise=noise, k=ut_kappa, a=ut_alpha)
    call generate_observation(obs, tgt, meas, meas_sig)
    call initialize_sr_ukf(obs, meas, init_sig_scale*meas_sig, filter, noise=noise, k=ut_kappa, a=ut_alpha, max_iterations=max_iterations)
    if (debug(1)) call print_status(t, obs, tgt, filter)

    do while (tgt(2) > 0.0_dp)
        !! evolve simulation time
        t = t + dt
        call dynamics_ballistic(tgt, dt, no_opt_data)
        !! incorporate new measurement
        call generate_observation(obs, tgt, meas, meas_sig)
        call filter_measurement_update(obs, meas, meas_sig, filter, t, state_dynamics_model, opt_data)
        if (debug(1)) call print_status(t, obs, tgt, filter)
        if (debug(4)) call dump_states(in_run_stats_fid, state_model_ii, t, obs, tgt, filter%state_estimate)
    end do
    if (debug(5)) call dump_states(end_run_stats_fid, state_model_ii, t, obs, tgt, filter%state_estimate)

    se_err(1:6,trial_ii) = filter%state_estimate - tgt
    se_err(7,trial_ii) = sqrt(se_err(1,trial_ii)**2 + se_err(2,trial_ii)**2) !! position error
    se_err(8,trial_ii) = sqrt(se_err(3,trial_ii)**2 + se_err(4,trial_ii)**2) !! velocity error
    se_err(9,trial_ii) = sqrt(se_err(5,trial_ii)**2 + se_err(6,trial_ii)**2) !! acceleration error
    call generate_observation(obs, tgt, tgt_pol)
    call generate_observation(obs, filter%state_estimate, trk_pol)
    se_err_pol(:,trial_ii) = trk_pol - tgt_pol
    if (debug(2)) then
        call print_status(t, obs, tgt, filter, .true.)
        write(*,'(a)') repeat('=', 32)
        write(*,'(a)') 'State Dynamics Model '//state_dynamics_model_name
        write(*,'(a,e13.6,a,f0.1,a)') 'Filter Noise: ',filter%process_noise,' ft/sec**2 (',filter%process_noise/g,' g)'
        write(*,'(6(a,e13.6))') 'Target Initial Conditions :: x: ',tgt0(1),', y: ',tgt0(2), &
                                ', vx: ',tgt0(3),', vy: ',tgt0(4),', ax: ',tgt0(5),', ay: ',tgt0(6)
        write(*,'(6(a,e13.6))') '              Track Error :: x: ',se_err(1,trial_ii), ', y: ',se_err(2,trial_ii), &
                                                          ', vx: ',se_err(3,trial_ii),', vy: ',se_err(4,trial_ii), &
                                                          ', ax: ',se_err(5,trial_ii),', ay: ',se_err(6,trial_ii)
        write(*,'(a,/)') repeat('=', 32)
    end if

    end do !! trial_ii
    !$omp end parallel do
    flush(in_run_stats_fid)
    flush(end_run_stats_fid)

    if (debug(6)) then
        call dump_summary(sr_ukf_summary_fid, max_trials, state_model_ii, noise, init_sig_scale, init_sig, &
                          meas_sig, se_err, se_err_pol, ut_alpha, ut_lambda, ut_kappa, max_iterations)
        flush(sr_ukf_summary_fid)
    end if
    if (debug(3)) then
        write(*,'(a)') repeat('=', 32)
        write(*,'(6(a,e13.6))') '               Track RMSE :: x: ',rmse(se_err(1,:),0.0_dp), ', y: ',rmse(se_err(2,:),0.0_dp), &
                                                          ', vx: ',rmse(se_err(3,:),0.0_dp),', vy: ',rmse(se_err(4,:),0.0_dp), &
                                                          ', ax: ',rmse(se_err(5,:),0.0_dp),', ay: ',rmse(se_err(6,:),0.0_dp)
        write(*,'(3(a,e13.6))') '       Track Summary RMSE :: pos: ',rmse(se_err(7,:), 0.0_dp), &
                                                           ', vel: ',rmse(se_err(8,:), 0.0_dp), &
                                                           ', acc: ',rmse(se_err(9,:), 0.0_dp)
        se_err = abs(se_err)
        do concurrent (i=7:9)
            call sort(se_err(i,:))
        end do
        do i=10,90,10
            trial_ii = min(max_trials, max(1, floor(i*max_trials/100.0_dp)))
            write(*,'(a,i0,3(a,e13.6))') 'percentile(',i,') RMSE :: pos: ',se_err(7,trial_ii), &
                                                                 ', vel: ',se_err(8,trial_ii), &
                                                                 ', acc: ',se_err(9,trial_ii)
        end do
        write(*,'(a,/)') repeat('=', 32)
    end if

    end do !! max_iterations_ii
    end do !! ut_lambda_ii
    end do !! ut_alpha_ii
    end do !! state_model_ii
    end do !! noise_ii
    end do !! meas_sig4_ii
    end do !! meas_sig3_ii
    end do !! meas_sig2_ii
    end do !! meas_sig1_ii
    end do !! init_scale_ii
    end do !! init_sig4_ii
    end do !! init_sig3_ii
    end do !! init_sig2_ii
    end do !! init_sig1_ii

    if (debug(4)) close(in_run_stats_fid)
    if (debug(5)) close(end_run_stats_fid)
    if (debug(6)) close(sr_ukf_summary_fid)

end program ex_sr_ukf
