subroutine transferdata (statein, bessein)

    use results

    class(CAMBdata) :: statein
    class(TRange) :: bessein

    real(dl) :: tau0, chi0, curvature_radius
    logical :: flat, closed
    integer :: num_redshiftwindows, num_extra_redshiftwindows

    tau0 = statein%tau0
    chi0 = statein%chi0
    falt = statein%flat
    closed = statein%closed
    curvature_radius = statein%curvature_radius
    num_redshiftwindows = statein%num_redshiftwindows
    num_extra_redshiftwindows = statein%num_extra_redshiftwindows

    !statein%TimeSteps%npoints
    !statein%TimeSteps%points
    !statein%TimeSteps%IndexOf 
    !statein%TimeSteps%dpoints 
    !statein%TimeSteps%Highest 
    !statein%TimeSteps%Lowest 
    
    !statein%ThermoData%tau_start_redshiftwindows

    !bessein%IndexOf 
    !bessein%points

    !functions or procedures 
    !statein%rofChi

end subroutine transferdata

function statindexof ()

    integer :: statindexof
    
    statindexof = 0
    
    return

end function statindexof
