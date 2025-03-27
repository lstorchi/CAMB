subroutine transferdata (statein, bessein, &
    s_tau0, s_chi0, s_curvature_radius, s_tau_start_redshiftwindows, &
    s_flat, s_closed, s_num_redshiftwindows, s_num_extra_redshiftwindows, &
    s_npoints, s_lowest, s_highest)

    use results

    class(CAMBdata), intent(in) :: statein
    class(TRanges), intent(in) :: bessein

    double precision, intent(inout) :: s_tau0, s_chi0, s_curvature_radius, &
        s_tau_start_redshiftwindows
    logical, intent(inout) :: s_flat, s_closed
    integer, intent(inout) :: s_num_redshiftwindows, s_num_extra_redshiftwindows, &
        s_npoints
    double precision, intent(inout) :: s_lowest, s_highest

    print *, "transferdata"

    s_tau0 = statein%tau0
    s_chi0 = statein%chi0
    s_falt = statein%flat
    s_closed = statein%closed
    s_curvature_radius = statein%curvature_radius
    s_num_redshiftwindows = statein%num_redshiftwindows
    s_num_extra_redshiftwindows = statein%num_extra_redshiftwindows
    s_npoints = statein%TimeSteps%npoints
    s_lowest = statein%TimeSteps%Lowest
    s_highest = statein%TimeSteps%Highest
    s_tau_start_redshiftwindows = statein%ThermoData%tau_start_redshiftwindows

    print *, "tau0 = ", tau0
    print *, "chi0 = ", chi0
    print *, "flat = ", flat
    print *, "closed = ", closed
    print *, "curvature_radius = ", curvature_radius
    print *, "num_redshiftwindows = ", num_redshiftwindows
    print *, "num_extra_redshiftwindows = ", num_extra_redshiftwindows
    print *, "npoints = ", npoints
    print *, "lowest = ", lowest
    print *, "highest = ", highest
    print *, "tau_start_redshiftwindows = ", tau_start_redshiftwindows

end subroutine transferdata

function statindexof ()

    !statein%TimeSteps%IndexOf  RangeUtils.f90 procedure :: IndexOf => TRanges_IndexOf

    integer :: statindexof
    
    statindexof = 0
    
    return

end function statindexof

function besseindexof ()

    !bessein%IndexOf  RangeUtils.f90 

    integer :: besseindexof
    
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