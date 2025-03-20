module ukf
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: nmi2ft_dp
implicit none
private

    real(dp), parameter :: default_range = 500.0_dp*nmi2ft_dp

    type :: sr_ukf_type
        real(dp) :: k, a, y, b
        real(dp), allocatable :: trk_est(:), sqrt_cov(:,:)
    end type sr_ukf_type

    type :: cart2d_6
        real(dp) :: x, y, vx, vy, ax, ay
    end type

    type :: pol2d_4
        real(dp) :: range, bearing_angle, range_rate, bearing_angle_rate
        logical :: have_range, have_bearing_angle, have_range_rate, have_bearing_angle_rate
    end type

    public :: sr_ukf_type, cart2d_6, pol2d_4, initialize_sr_ukf

contains

    pure elemental subroutine initialize_sr_ukf(obs, meas, trk_est_dim, filter, k, a, b)
        type(cart2d_6), intent(in) :: obs
        type(pol2d_4), intent(in) :: meas
        integer, intent(in) :: trk_est_dim
        type(sr_ukf_type), intent(out) :: filter
        real(dp), intent(in), optional :: k, a, b
        real(dp) :: rng_use, cos_ang, sin_ang, vt
        ! set unscented transform parameters kappa (k), alpha (a), and beta (b)
        if (present(k) .and. present(a) .and. present(b)) then
            filter%k = k
            filter%a = a
            filter%b = b
        else
            filter%k = 0.0_dp
            filter%a = 1.0_dp
            filter%b = 2.0_dp
            ! alternatively consider k = 3.0_dp - trk_est_dim, a = 1.0_dp, b = 2.0_dp
        end if
        ! parameters determine lambda (y)
        filter%y = filter%a**2 * (trk_est_dim + filter%k) - trk_est_dim
        ! allocate track estimate trk_est and covariance matrix square root sqrt_cov
        allocate(filter%trk_est(trk_est_dim), filter%sqrt_cov(trk_est_dim,trk_est_dim))
        ! initialize track estimate based on observer state and measurement
        if (meas%have_range) then
            rng_use = meas%range
        else
            rng_use = default_range
        end if
        if (meas%have_bearing_angle) then
            cos_ang = cos(meas%bearing_angle)
            sin_ang = sin(meas%bearing_angle)
            filter%trk_est(1) = obs%x + rng_use*cos_ang
            filter%trk_est(2) = obs%y + rng_use*sin_ang
        else
            error stop 'track initialization not possible without angle information'
        end if
        if (meas%have_range_rate) then
            filter%trk_est(3) = obs%vx + meas%range_rate*cos_ang
            filter%trk_est(4) = obs%vy + meas%range_rate*sin_ang
        else
            filter%trk_est(3) = 0.0_dp
            filter%trk_est(4) = 0.0_dp
        end if
        if (meas%have_bearing_angle_rate) then
            vt = rng_use*meas%bearing_angle_rate
            filter%trk_est(3) = filter%trk_est(3) - vt*sin_ang
            filter%trk_est(4) = filter%trk_est(4) + vt*cos_ang
        end if
        ! assume no acceleration
        filter%trk_est(5) = 0.0_dp
        filter%trk_est(6) = 0.0_dp
    end

end module ukf


program ex_ukf
use, non_intrinsic :: ukf
implicit none
end program ex_ukf
