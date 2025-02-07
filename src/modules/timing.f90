module timing
use, non_intrinsic :: kinds, only: i64, dp
use, non_intrinsic :: constants, only: debug
implicit none
private

    type :: timer
        private
        integer(kind=i64) :: count1
        integer(kind=i64) :: count_rate = 0.0_i64
        integer(kind=i64) :: count2
        real(kind=dp) :: elapsed = 0.0_dp
    end type timer

    public :: timer, tic, toc, get_elapsed

contains

    impure elemental subroutine tic(t)
        type(timer), intent(out) :: t
        call system_clock(count=t%count1, count_rate=t%count_rate)
    end subroutine tic

    impure elemental subroutine toc(t)
        type(timer), intent(inout) :: t
        call system_clock(count=t%count2)
        if (debug) then
            if (t%count_rate <= 0_i64) then
                error stop 'module TIMING :: must call tic(timer) before toc(timer)'
            end if
        end if
        t%elapsed = real(max(t%count2 - t%count1, 1_i64), kind=dp)/real(t%count_rate, kind=dp)
    end subroutine toc

    pure elemental function get_elapsed(t) result(elapsed_sec)
        type(timer), intent(in) :: t
        real(kind=dp) :: elapsed_sec
        if (debug) then
            if (t%elapsed <= 0.0_dp) then
                error stop 'module TIMING :: must call tic(timer) and toc(timer) before get_elapsed(timer)'
            end if
        end if
        elapsed_sec = t%elapsed
    end function get_elapsed

end module timing
