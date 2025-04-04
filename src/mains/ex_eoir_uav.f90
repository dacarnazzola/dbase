module eoir_uav_mod
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: nmi2ft_dp, ft2nmi_dp, deg2rad_dp, km2nmi_dp
use, non_intrinsic :: system, only: debug_error_condition, nearly
use, non_intrinsic :: math_utils, only: linspace
use, non_intrinsic :: statistics, only: normalize_min_max
implicit none
private

    real(dp), parameter :: R_air = 1716.46_dp !! https://en.wikipedia.org/wiki/Gas_constant#Specific_gas_constant

    public :: dp, km2nmi_dp, endurance_hr, cost_m, leg_hr, calc_patches, linspace, normalize_min_max

contains

    pure elemental function endurance_hr(mach, alt_kft) result(val)
        real(dp), intent(in) :: mach, alt_kft
        real(dp) :: val
        val = -18.75_dp*mach**2 + 8.0893_dp*mach + 0.01_dp*alt_kft**2 + 0.05_dp*alt_kft + 9.2105_dp
    end function endurance_hr

    pure elemental function cost_m(mach, alt_kft) result(val)
        real(dp), intent(in) :: mach, alt_kft
        real(dp) :: val
        val = 50.0_dp*mach**2 - 35.0_dp*mach + 0.03_dp*alt_kft**2 - 0.2_dp*alt_kft + 11.0_dp
    end function cost_m

    pure elemental function imperfect_gamma(T_R) result(val) !! https://www.grc.nasa.gov/www/BGH/realspec.html
        real(dp), intent(in) :: T_R
        real(dp) :: val
        real(dp) :: Theta_T_R
        call debug_error_condition(T_R <= 0.0_dp, 'module EOIR_UAV_MOD :: funciton imperfect_gamma invalid for T_R <= 0.0')
        Theta_T_R = 5500.0_dp/T_R
        val = 1.0_dp + 0.4_dp/(1.0_dp + 0.4_dp*(Theta_T_R**2*exp(Theta_T_R)/(exp(Theta_T_R)-1.0_dp)**2))
    end function imperfect_gamma

    pure elemental function t_alt(alt_ft) result(val) !! https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/atmos.html
        real(dp), intent(in) :: alt_ft
        real(dp) :: val
        if (alt_ft < 36152.0_dp) then
            val = 59.0_dp - 0.00356_dp*alt_ft
        else if(alt_ft > 82345.0_dp) then
            val = -205.05_dp + 0.00164_dp*alt_ft
        else !! 36152.0 <= alt_ft <= 82345
            val = -70.0_dp
        end if
        val = val + 459.67_dp
    end function t_alt

    pure elemental function mach1(alt_ft) result(val)
        real(dp), intent(in) :: alt_ft
        real(dp) :: val
        real(dp) :: T_R_at_alt, y_at_alt
        T_R_at_alt = t_alt(alt_ft)
        y_at_alt = imperfect_gamma(T_R_at_alt)
        val = sqrt(y_at_alt*R_air*T_R_at_alt)
    end function mach1

    pure elemental function leg_hr(leg_nmi, mach, alt_kft) result(val)
        real(dp), intent(in) :: leg_nmi, mach, alt_kft
        real(dp) :: val
        val = leg_nmi*nmi2ft_dp/(mach1(alt_kft*1000.0_dp)*mach)/3600.0_dp
    end function leg_hr

    pure elemental subroutine calc_patches(mach, alt_kft, fov_az_hw_deg, fov_el_hw_deg, az_center_deg, el_center_deg, &
                                           mission_width_nmi, mission_length_nmi, &
                                           patch_area_nmi2, mission_patches, total_route_nmi, total_route_hr)
        real(dp), intent(in) :: mach, alt_kft, fov_az_hw_deg, fov_el_hw_deg, az_center_deg, el_center_deg, &
                                mission_width_nmi, mission_length_nmi
        real(dp), intent(out) :: patch_area_nmi2
        integer, intent(out) :: mission_patches
        real(dp), intent(out) :: total_route_nmi, total_route_hr
        real(dp) :: patch_width_nmi, patch_length_nmi
        integer :: nrow, ncol
        call debug_error_condition(.not.(nearly(az_center_deg, 0.0_dp) .and. nearly(el_center_deg, 0.0_dp)), &
                                   'module EOIR_UAV_MOD :: function patch_area_nmi2 only implemented for pure lookdown case')
        patch_width_nmi = alt_kft*1000.0_dp*ft2nmi_dp*tan(fov_az_hw_deg*deg2rad_dp)*2.0_dp
        patch_length_nmi = alt_kft*1000.0_dp*ft2nmi_dp*tan(fov_el_hw_deg*deg2rad_dp)*2.0_dp
        patch_area_nmi2 = patch_length_nmi*patch_width_nmi ! area = length * width
        nrow = ceiling(mission_length_nmi/patch_length_nmi) ! round up for 100% coverage in length
        ncol = ceiling(mission_width_nmi/patch_width_nmi) ! round up for 100% coverage in width
        mission_patches = nrow*ncol ! area = length * width
        total_route_nmi = min(nrow*(ncol - 1)*patch_width_nmi + (nrow - 1)*patch_length_nmi, &
                              ncol*(nrow - 1)*patch_length_nmi + (ncol - 1)*patch_width_nmi) + &
                              0.5_dp*(patch_width_nmi + patch_length_nmi) !! add half of the patch length and width for entry/exit to pattern
        total_route_hr = leg_hr(total_route_nmi, mach, alt_kft)
    end subroutine calc_patches

end module eoir_uav_mod


program main
use, non_intrinsic :: eoir_uav_mod
implicit none

    integer, parameter :: mach_list_size = 11, alt_list_size = 21
    real(dp), parameter :: ingress_nmi = 100.0_dp*km2nmi_dp, egress_nmi = 100.0_dp*km2nmi_dp, &
                           mission_width_nmi = 100.0_dp*km2nmi_dp, mission_length_nmi = 100.0_dp*km2nmi_dp, &
                           fov_deg_list(*) = [15.0_dp, 30.0_dp, 60.0_dp], sensor_cost_m_list(*) = [0.05_dp, 1.0_dp, 10.0_dp], &
                           min_mach = 0.40_dp, max_mach = 0.90_dp, min_alt_kft = 5.0_dp, max_alt_kft = 25.0_dp

    integer :: mach_ii, alt_ii, fov_ii, mission_patches, fid, c_ii
    real(dp) :: mach, alt_kft, fov_deg, total_endurance_hr, ingress_time_hr, egress_time_hr, mission_time_hr, &
                mission_patch_nmi2, airframe_cost_m, sensor_cost_m, total_cost_m, mission_route_nmi, mission_route_hr, &
                mission_reps, mission_reps_per_cost_m, mission_reps_per_hour, mission_reps_per_hour_per_cost_m, &
                mission_one_rep_time_hr, mission_one_rep_time_hr_per_cost_m, &
                mach_list(mach_list_size), alt_list(alt_list_size), results(29,mach_list_size*alt_list_size*size(fov_deg_list)), &
                absolute_performance_weight, efficiency_weight, total_weight
    character(len=128) :: fmt_str

    if (command_argument_count() == 2) then
        call get_command_argument(1, fmt_str)
        read(fmt_str,*) absolute_performance_weight
        call get_command_argument(2, fmt_str)
        read(fmt_str,*) efficiency_weight
    else
        absolute_performance_weight = 0.50_dp
        efficiency_weight = 0.50_dp
    end if
    total_weight = absolute_performance_weight + efficiency_weight

    call linspace(mach_list, min_mach, max_mach)
    call linspace(alt_list, min_alt_kft, max_alt_kft)

    !$omp parallel do default(firstprivate) shared(results)
    do mach_ii=1,mach_list_size
        mach = mach_list(mach_ii)
        do alt_ii=1,alt_list_size
            alt_kft = alt_list(alt_ii)
            total_endurance_hr = endurance_hr(mach, alt_kft)
            ingress_time_hr = leg_hr(ingress_nmi, mach, alt_kft)
            egress_time_hr = leg_hr(egress_nmi, mach, alt_kft)
            mission_time_hr = total_endurance_hr - ingress_time_hr - egress_time_hr
            airframe_cost_m = cost_m(mach, alt_kft)
            do fov_ii=1,size(fov_deg_list)
                c_ii = fov_ii + (alt_ii - 1)*size(fov_deg_list) + (mach_ii - 1)*alt_list_size*size(fov_deg_list)
                fov_deg = fov_deg_list(fov_ii)
                sensor_cost_m = sensor_cost_m_list(fov_ii)
                total_cost_m = airframe_cost_m + sensor_cost_m
                call calc_patches(mach, alt_kft, fov_deg*0.5_dp, fov_deg*0.5_dp, 0.0_dp, 0.0_dp, &
                                  mission_width_nmi, mission_length_nmi, &
                                  mission_patch_nmi2, mission_patches, mission_route_nmi, mission_route_hr)
                mission_reps = mission_time_hr/mission_route_hr
                mission_reps_per_cost_m = mission_reps/total_cost_m
                mission_reps_per_hour = 1.0_dp/mission_route_hr
                mission_reps_per_hour_per_cost_m = mission_reps_per_hour/total_cost_m
                mission_one_rep_time_hr = ingress_time_hr + mission_route_hr + egress_time_hr
                mission_one_rep_time_hr_per_cost_m = mission_one_rep_time_hr/total_cost_m
                !! fill results
                results( 1,c_ii) = mach
                results( 2,c_ii) = alt_kft
                results( 3,c_ii) = fov_deg
                results( 4,c_ii) = ingress_time_hr
                results( 5,c_ii) = egress_time_hr
                results( 6,c_ii) = mission_time_hr
                results( 7,c_ii) = total_endurance_hr
                results( 8,c_ii) = mission_patch_nmi2
                results( 9,c_ii) = real(mission_patches, kind=dp)
                results(10,c_ii) = mission_route_nmi
                results(11,c_ii) = mission_route_hr
                results(12,c_ii) = airframe_cost_m
                results(13,c_ii) = sensor_cost_m
                results(14,c_ii) = total_cost_m
                results(15,c_ii) = mission_reps
                results(16,c_ii) = mission_reps_per_cost_m
                results(17,c_ii) = mission_reps_per_hour
                results(18,c_ii) = mission_reps_per_hour_per_cost_m
                results(19,c_ii) = mission_one_rep_time_hr
                results(20,c_ii) = mission_one_rep_time_hr_per_cost_m
            end do
        end do
    end do
    !$omp end parallel do

    call normalize_min_max(results(15,:), results(21,:)) !! higher = better
    call normalize_min_max(results(16,:), results(22,:)) !! higher = better
    call normalize_min_max(results(17,:), results(23,:)) !! higher = better
    call normalize_min_max(results(18,:), results(24,:)) !! higher = better
    call normalize_min_max(results(19,:), results(25,:)); results(25,:) = 1.0_dp - results(25,:) !! lower = better
    call normalize_min_max(results(20,:), results(26,:)); results(26,:) = 1.0_dp - results(26,:) !! lower = better

    results(27,:) = (results(21,:) + results(23,:) + results(25,:))/3.0_dp !! average absolute scores
    results(28,:) = (results(22,:) + results(24,:) + results(26,:))/3.0_dp !! average efficiency scores
    results(29,:) = (absolute_performance_weight*results(27,:) + efficiency_weight*results(28,:))/total_weight !! weighted average total scores

    open(newunit=fid, file='/valinor/eoir-uav-trade-study.csv', action='write')
    write(fid,'(a,f0.4,a,f0.4,a)') 'mach,alt_kft,sensor_fov_deg,'// &
                     'ingress_time_hr,egress_time_hr,mission_time_hr,total_endurance_hr,'// &
                     'survey_patch_nmi2,survey_patches_per_mission,mission_route_nmi,mission_route_hr,'// &
                     'airframe_cost_m,sensor_cost_m,total_cost_m,'// &
                     'mission_max_reps,mission_max_reps_per_cost_m,'// &
                     'mission_reps_per_hour,mission_reps_per_hour_per_cost_m,'// &
                     'mission_one_rep_time_hr,mission_one_rep_time_hr_per_cost_m,'// &
                     'norm_max_reps,norm_max_reps_per_cost_m,norm_reps_per_hour,norm_reps_per_hour_per_cost_m,'// &
                     'norm_one_rep_time_hr,norm_one_rep_time_hr_per_cost_m,'// &
                     'average_absolute_score',absolute_performance_weight,',average_efficiency_score',efficiency_weight, &
                     ',average_total_score'
    write(fmt_str,'(a,i0,a)') '(e22.15,',size(results,dim=2)-1,'(",",e22.15))'
    do c_ii=1,size(results,dim=2)
        write(unit=fid, fmt=fmt_str) results(:,c_ii)
    end do
    close(fid)

end program main
