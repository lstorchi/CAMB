subroutine transferdata (statein, bessein, &
    s_tau0, s_chi0, s_curvature_radius, s_tau_start_redshiftwindows, &
    s_flat, s_closed, s_num_redshiftwindows, s_num_extra_redshiftwindows, &
    s_npoints, s_lowest, s_highest, b_npoints, b_lowest, b_highest)

    use results

    implicit none

    class(CAMBdata), intent(in) :: statein
    class(TRanges), intent(in) :: bessein

    double precision, intent(inout) :: s_tau0, s_chi0, s_curvature_radius, &
        s_tau_start_redshiftwindows
    logical, intent(inout) :: s_flat, s_closed
    integer, intent(inout) :: s_num_redshiftwindows, s_num_extra_redshiftwindows, &
        s_npoints, b_npoints
    double precision, intent(inout) :: s_lowest, s_highest, b_lowest, b_highest

    print *, "transferdata"

    s_tau0 = statein%tau0
    s_chi0 = statein%chi0
    s_flat = statein%flat
    s_closed = statein%closed
    s_curvature_radius = statein%curvature_radius
    s_num_redshiftwindows = statein%num_redshiftwindows
    s_num_extra_redshiftwindows = statein%num_extra_redshiftwindows
    s_npoints = statein%TimeSteps%npoints
    s_lowest = statein%TimeSteps%Lowest
    s_highest = statein%TimeSteps%Highest
    s_tau_start_redshiftwindows = statein%ThermoData%tau_start_redshiftwindows

    b_npoints = bessein%npoints
    b_highest = bessein%Highest
    b_lowest = bessein%Lowest

    print *, "s_tau0 = ", s_tau0
    print *, "s_chi0 = ", s_chi0
    print *, "s_flat = ", s_flat
    print *, "s_closed = ", s_closed
    print *, "s_curvature_radius = ", s_curvature_radius
    print *, "s_num_redshiftwindows = ", s_num_redshiftwindows
    print *, "s_num_extra_redshiftwindows = ", s_num_extra_redshiftwindows
    print *, "s_npoints = ", s_npoints
    print *, "s_lowest = ", s_lowest
    print *, "s_highest = ", s_highest
    print *, "s_tau_start_redshiftwindows = ", s_tau_start_redshiftwindows
    
    print *, "b_npoints = ", b_npoints
    print *, "b_lowest = ", b_lowest
    print *, "b_highest = ", b_highest

end subroutine transferdata

function statbesseindexof (count, R, npoints, Highest, tau) 
#ifdef USEACC
!$ACC ROUTINE SEQ
#endif
    use RangeUtils
    !statein%TimeSteps%IndexOf  RangeUtils.f90 procedure :: IndexOf => TRanges_IndexOf
    ! to test it compare respect to State.IndexOf 
    integer :: statbesseindexof
    double precision, intent(in) :: tau
    integer , intent(in) :: count
    type(TRange), intent(in) :: R(count)
    integer, intent(in) :: npoints
    double precision, intent(in) :: Highest
    type(TRange), pointer :: AReg
    integer :: i

    !print *, "I am in simple IndexOf "
    !print *, "  npoints: ", npoints
    !print *, "  count  : ", count
    !print *, "  tau    : ", tau 
    !print *, "  Hihest : ", Highest
    !print *, "   1st Low        :", R(1)%Low
    !print *, "   1st High       :", R(1)%High 
    !print *, "   1st delta      :", R(1)%delta
    !print *, "   1st start_index:", R(1)%start_index

    statbesseindexof=0
    do i=1, count
        associate(AReg => R(i))
            if (tau < AReg%High .and. tau >= AReg%Low) then
                if (AReg%IsLog) then
                    statbesseindexof = AReg%start_index + int(log(tau / AReg%Low) / AReg%delta)
                else
                    statbesseindexof = AReg%start_index + int((tau - AReg%Low) / AReg%delta)
                end if
                return
            end if
        end associate
    end do
    
    if (tau >= Highest) then
        statbesseindexof = npoints
    else
        stop
    end if
 
    return

end function statbesseindexof

function staterofchi (flat, closed, chi) 
#ifdef USEACC
!$ACC ROUTINE SEQ
#endif
    !statein%rofChi
    logical, intent(in) :: flat, closed
    double precision , intent(in) :: chi
    double precision :: staterofchi

    if (flat) then
        staterofchi=chi
    else if (closed) then
        staterofchi=sin(chi)
    else
        staterofchi=sinh(chi)
    endif
    
    return

end function staterofchi
