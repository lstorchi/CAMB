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
    !statein%TimeSteps%points
    !statein%TimeSteps%dpoints 
    !bessein%points

    !functions or procedures 
    !bessein%IndexOf  RangeUtils.f90 
    !statein%rofChi
    !statein%TimeSteps%IndexOf  RangeUtils.f90 procedure :: IndexOf => TRanges_IndexOf

end subroutine transferdata

function statindexof ()

    integer :: statindexof
    
    statindexof = 0
    
    return

end function statindexof
