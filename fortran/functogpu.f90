subroutine transferdata (statein, &
    s_tau0, s_chi0, s_curvature_radius, s_tau_start_redshiftwindows, &
    s_flat, s_closed, s_num_redshiftwindows, s_num_extra_redshiftwindows, &
    s_npoints, s_lowest, s_highest)

    use results

    implicit none

    class(CAMBdata), intent(in) :: statein
    !class(TRanges), intent(in) :: bessein

    double precision, intent(inout) :: s_tau0, s_chi0, s_curvature_radius, &
        s_tau_start_redshiftwindows
    logical, intent(inout) :: s_flat, s_closed
    integer, intent(inout) :: s_num_redshiftwindows, s_num_extra_redshiftwindows, &
        s_npoints
    double precision, intent(inout) :: s_lowest, s_highest

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

    print *, "tau0 = ", s_tau0
    print *, "chi0 = ", s_chi0
    print *, "flat = ", s_flat
    print *, "closed = ", s_closed
    print *, "curvature_radius = ", s_curvature_radius
    print *, "num_redshiftwindows = ", s_num_redshiftwindows
    print *, "num_extra_redshiftwindows = ", s_num_extra_redshiftwindows
    print *, "npoints = ", s_npoints
    print *, "lowest = ", s_lowest
    print *, "highest = ", s_highest
    print *, "tau_start_redshiftwindows = ", s_tau_start_redshiftwindows

end subroutine transferdata

function statindexof (tau)

    !statein%TimeSteps%IndexOf  RangeUtils.f90 procedure :: IndexOf => TRanges_IndexOf
    ! to test it compare respect to State.IndexOf 

    integer :: statindexof
    double precision, intent(in) :: tau
    integer :: pointstep, i
    !  importante type(TRange), allocatable :: R(:) 
    !   that is  type(TRange), allocatable :: R(:)

    pointstep=0
    !do i=1, this%count
    !    associate(AReg => this%R(i)) here I need to check  
    !        if (tau < AReg%High .and. tau >= AReg%Low) then
    !            if (AReg%IsLog) then
    !                pointstep = AReg%start_index + int(log(tau / AReg%Low) / AReg%delta)
    !            else
    !                pointstep = AReg%start_index + int((tau - AReg%Low) / AReg%delta)
    !            end if
    !            return
    !        end if
    !    end associate
    !end do
    
    !if (tau >= this%Highest) then
    !    pointstep = this%npoints
    !else
    !    print *, "tau=", tau, ",this%Highest=", this%Highest
    !    call MpiStop('TRanges_IndexOf: value out of range')
    !end if
 
    return

end function statindexof

function besseindexof (tau)

    !bessein%IndexOf  RangeUtils.f90 

    integer :: besseindexof
    double precision, intent(in) :: tau
    
    besseindexof = 0
    
    return

end function besseindexof

function staterofchi (flat, closed, chi) 
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