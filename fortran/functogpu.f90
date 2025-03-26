subroutine transferdata (statein, bessein)

    use results

    class(CAMBdata) :: statein
    class(TRanges) :: bessein

    real(dl) :: tau0, chi0, curvature_radius, tau_start_redshiftwindows
    logical :: flat, closed
    integer :: num_redshiftwindows, num_extra_redshiftwindows, &
        npoints
    double precision, dimension(:), allocatable :: points, dpoints
    double precision, dimension(:), allocatable :: b_points
    double precision :: lowest, highest

    tau0 = statein%tau0
    chi0 = statein%chi0
    falt = statein%flat
    closed = statein%closed
    curvature_radius = statein%curvature_radius
    num_redshiftwindows = statein%num_redshiftwindows
    num_extra_redshiftwindows = statein%num_extra_redshiftwindows
    npoints = statein%TimeSteps%npoints
    lowest = statein%TimeSteps%Lowest
    highest = statein%TimeSteps%Highest
    tau_start_redshiftwindows = statein%ThermoData%tau_start_redshiftwindows

    ! to copy the data 
    if (allocated(statein%TimeSteps%points)) then
        allocate(points( size(statein%TimeSteps%points) ))
        points = statein%TimeSteps%points
    end if
    if (allocated(statein%TimeSteps%dpoints)) then
        allocate(dpoints( size(statein%TimeSteps%dpoints) ))
        dpoints = statein%TimeSteps%dpoints
    end if
    if (allocated(statein%TimeSteps%points)) then
        allocate(b_points( size(statein%TimeSteps%points) ))
        b_points = statein%TimeSteps%points
    end if

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