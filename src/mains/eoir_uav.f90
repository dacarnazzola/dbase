module eoir_uav_mod
use, non_intrinsic :: kinds, only: dp
use, non_intrinsic :: constants, only: nmi2ft_dp, ft2nmi_dp, deg2rad_dp, km2nmi_dp
use, non_intrinsic :: system, only: debug_error_condition, nearly
implicit none
private

    real(dp), parameter :: R_air = 1716.46_dp !! https://en.wikipedia.org/wiki/Gas_constant#Specific_gas_constant

    public :: dp, km2nmi_dp, endurance_hr, cost_m, leg_hr, calc_patches

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

    real(dp), parameter :: ingress_nmi = 100.0_dp*km2nmi_dp, egress_nmi = 100.0_dp*km2nmi_dp, &
                           mission_width_nmi = 100.0_dp*km2nmi_dp, mission_length_nmi = 100.0_dp*km2nmi_dp, &
                           fov_deg_list(*) = [15.0_dp, 30.0_dp, 60.0_dp], sensor_cost_m_list(*) = [0.05_dp, 1.0_dp, 10.0_dp]

    integer :: mach_ii, alt_ii, fov_ii, mission_patches
    real(dp) :: mach, alt_kft, fov_deg, total_endurance_hr, ingress_time_hr, egress_time_hr, mission_time_hr, &
                mission_patch_nmi2, airframe_cost_m, sensor_cost_m, total_cost_m, mission_route_nmi, mission_route_hr, &
                mission_reps, mission_reps_per_cost_m, missions_per_hour, missions_per_hour_per_cost_m

    write(*,'(a)') 'mach,alt_kft,sensor_fov_deg,'// &
                   'ingress_hr,egress_hr,mission_hr,total_endurance_hr,'// &
                   'survey_patch_nmi2,survey_patches_per_mission,mission_route_nmi,mission_route_hr,'// &
                   'airframe_cost_m,sensor_cost_m,total_cost_m,'// &
                   'mission_max_reps,mission_max_reps_per_million,missions_per_hour,missions_per_hour_per_million'
    do mach_ii=40,90
        mach = real(mach_ii,dp)/100.0_dp
        do alt_ii=5,25
            alt_kft = real(alt_ii,dp)
            total_endurance_hr = endurance_hr(mach, alt_kft)
            ingress_time_hr = leg_hr(ingress_nmi, mach, alt_kft)
            egress_time_hr = leg_hr(egress_nmi, mach, alt_kft)
            mission_time_hr = total_endurance_hr - ingress_time_hr - egress_time_hr
            airframe_cost_m = cost_m(mach, alt_kft)
            do fov_ii=1,size(fov_deg_list)
                fov_deg = fov_deg_list(fov_ii)
                sensor_cost_m = sensor_cost_m_list(fov_ii)
                total_cost_m = airframe_cost_m + sensor_cost_m
                call calc_patches(mach, alt_kft, fov_deg*0.5_dp, fov_deg*0.5_dp, 0.0_dp, 0.0_dp, &
                                  mission_width_nmi, mission_length_nmi, &
                                  mission_patch_nmi2, mission_patches, mission_route_nmi, mission_route_hr)
                mission_reps = mission_time_hr/mission_route_hr
                mission_reps_per_cost_m = mission_reps/total_cost_m
                missions_per_hour = 1.0_dp/mission_route_hr
                missions_per_hour_per_cost_m = missions_per_hour/total_cost_m
                write(*,'(f0.2,2(",",f0.0),5(",",e13.6),",",i0,9(",",e13.6))') &
                                                        mach, alt_kft, fov_deg, &
                                                        ingress_time_hr, egress_time_hr, mission_time_hr, total_endurance_hr, &
                                                        mission_patch_nmi2, mission_patches, mission_route_nmi, mission_route_hr, &
                                                        airframe_cost_m, sensor_cost_m, total_cost_m, &
                                                        mission_reps, mission_reps_per_cost_m, &
                                                        missions_per_hour, missions_per_hour_per_cost_m
            end do
        end do
    end do

end program main
