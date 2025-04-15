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
!$ACC ROUTINE
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

function UseLimberGPU(l, CPin)
#ifdef USEACC
        !$ACC ROUTINE
#endif
        !Calculate lensing potential power using Limber rather than j_l integration
        !even when sources calculated as part of temperature calculation
        !(Limber better on small scales unless step sizes made much smaller)
        !This affects speed, esp. of non-flat case
        use model

        logical :: UseLimberGPU
        integer l
        Type(CAMBParams) :: CPin
    
        !note increasing non-limber is not neccessarily more accurate unless AccuracyBoost much higher
        !use **0.5 to at least give some sensitivity to Limber effects
        !Could be lower but care with phi-T correlation at lower L
        if (CPin%SourceTerms%limber_windows) then
            UseLimberGPU = l >= CPin%SourceTerms%limber_phi_lmin
        else
            UseLimberGPU = l > 400 * (CPin%Accuracy%AccuracyBoost * CPin%Accuracy%LimberBoost)** 0.5
        end if
    
end function UseLimberGPU

subroutine SourceToTransfers(datasbin, &
    ThisCT, q_ix,  ThisSourcesin, CPin, ScaledSrcin, &
    ddScaledSrcin, max_etak_tensorin, max_etak_vectorin, &
    WantLateTimein, max_etak_scalarin, full_bessel_integrationin, &
    do_bispectrumin, max_bessels_l_indexin, IV, &
    xlimfracin, xlimminin, ajlin, ajlprin)
#ifdef USEACC
!$acc routine seq
#endif
    use CAMBmain
    use results
    use RangeUtils

    implicit none

    real(dl) :: xlimfracin, xlimminin
    real(dl), dimension(:,:), allocatable :: ajlin, ajlprin
    type(ClTransferData), target :: ThisCT 
    Type(TTimeSources) :: ThisSourcesin
    integer :: q_ix, max_bessels_l_indexin
    Type(CAMBParams) :: CPin
    type(IntegrationVars) :: IV
    real(dl), dimension(:,:,:) :: ScaledSrcin
    real(dl), dimension(:,:,:) :: ddScaledSrcin
    real(dl) :: max_etak_tensorin, max_etak_vectorin, max_etak_scalarin
    logical :: WantLateTimein
    logical :: full_bessel_integrationin, do_bispectrumin
    type(datastatebessel) :: datasbin   

    !call IntegrationVars_Init(IV, datasbin)
    ! to avoid a call 

    !print *, "allocated ajlin: ", allocated(ajlin)
    !print *, "allocated ajlprin: ", allocated(ajlprin)

    IV%Source_q(1,:)=0
    IV%Source_q(datasbin%s_npoints,:) = 0
    IV%Source_q(datasbin%s_npoints-1,:) = 0

    IV%q_ix = q_ix
    IV%q = ThisCT%q%points(q_ix)
    IV%dq= ThisCT%q%dpoints(q_ix)

    call InterpolateSources(IV, ThisSourcesin, CPin, ScaledSrcin, ddScaledSrcin, &
      max_etak_tensorin, max_etak_vectorin, WantLateTimein, max_etak_scalarin, &
      datasbin)

    !call DoSourceIntegration(IV, ThisCT, CPin, ThisSourcesin, &
    !        full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
    !        datasbin,xlimfracin,xlimminin,ajlin,ajlprin)

end subroutine SourceToTransfers

subroutine InterpolateSources(IV, ThisSourcesin, CPin, ScaledSrcin, &
    ddScaledSrcin, max_etak_tensorin, max_etak_vectorin, WantLateTimein, &
    max_etak_scalarin, datasbin)

    use CAMBmain
    use results

    implicit none
    integer i,khi,klo, step
    real(dl) xf,b0,ho,a0,ho2o6,a03,b03
    type(IntegrationVars) IV
    Type(CAMBParams) :: CPin
    Type(TTimeSources) :: ThisSourcesin
    real(dl), dimension(:,:,:) :: ScaledSrcin
    real(dl), dimension(:,:,:) :: ddScaledSrcin
    real(dl) :: max_etak_tensorin, max_etak_vectorin, max_etak_scalarin
    logical :: WantLateTimein
    type(datastatebessel) :: datasbin

    !     finding position of k in table Evolve_q to do the interpolation.

    !Can't use the following in closed case because regions are not set up (only points)
    !           klo = min(ThisSourcesin%Evolve_q%npoints-1,ThisSourcesin%Evolve_q%IndexOf(IV%q))
    !This is a bit inefficient, but thread safe
    klo=1
    do while ((IV%q > ThisSourcesin%Evolve_q%points(klo+1)).and.(klo < (ThisSourcesin%Evolve_q%npoints-1)))
        klo=klo+1
    end do

    khi=klo+1

    ho=ThisSourcesin%Evolve_q%points(khi)-ThisSourcesin%Evolve_q%points(klo)
    a0=(ThisSourcesin%Evolve_q%points(khi)-IV%q)/ho
    b0=(IV%q-ThisSourcesin%Evolve_q%points(klo))/ho
    ho2o6 = ho**2/6
    a03=(a0**3-a0)
    b03=(b0**3-b0)
    IV%SourceSteps = 0

    !Interpolating the source as a function of time for the present
    !wavelength.
    step=2
    do i=2, datasbin%s_npoints
        xf=IV%q*(datasbin%s_tau0-datasbin%s_points(i))
        if (CPin%WantTensors) then
            if (IV%q*datasbin%s_points(i) < max_etak_tensorin.and. xf > 1.e-8_dl) then
                step=i
                IV%Source_q(i,:) =a0*ScaledSrcin(klo,:,i)+&
                    b0*ScaledSrcin(khi,:,i)+(a03 *ddScaledSrcin(klo,:,i)+ &
                    b03*ddScaledSrcin(khi,:,i)) *ho2o6
            else
                IV%Source_q(i,:) = 0
            end if
        end if
        if (CPin%WantVectors) then
            if (IV%q*datasbin%s_points(i) < max_etak_vectorin.and. xf > 1.e-8_dl) then
                step=i
                IV%Source_q(i,:) =a0*ScaledSrcin(klo,:,i) + b0*ScaledSrcin(khi,:,i)+(a03 *ddScaledSrcin(klo,:,i)+ &
                    b03*ddScaledSrcin(khi,:,i)) *ho2o6
            else
                IV%Source_q(i,:) = 0
            end if
        end if

        if (CPin%WantScalars) then
            if ((DebugEvolution .or. WantLateTimein .or. IV%q*datasbin%s_points(i) < max_etak_scalarin) &
                .and. xf > 1.e-8_dl) then
                step=i
                IV%Source_q(i,:) = a0 * ScaledSrcin(klo,:,i) +  b0 * ScaledSrcin(khi,:,i) + (a03*ddScaledSrcin(klo,:,i) + &
                    b03 * ddScaledSrcin(khi,:,i)) * ho2o6
            else
                IV%Source_q(i,:) = 0
            end if
        end if
    end do
    IV%SourceSteps = step

    if (.not.datasbin%s_flat) then
        do i=1, ThisSourcesin%SourceNum
            call spline_def(datasbin%s_points,IV%Source_q(:,i),datasbin%s_npoints,&
                IV%ddSource_q(:,i))
        end do
    end if
 
end subroutine InterpolateSources

subroutine DoSourceIntegration(IV, ThisCT, CPin, ThisSourcesin, &
    full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
    datasbin,xlimfracin,xlimminin,ajlin,ajlprin) !for particular wave number q
    use CAMBmain
    use precision
    use model
    use results

    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT    
    real(dl), dimension(:,:), allocatable, intent(inout) :: ajlin, ajlprin
    real(dl) xlimfracin, xlimminin
    integer j,ll,llmax, max_bessels_l_indexin 
    real(dl) nu
    real(dl) :: sixpibynu
    Type(CAMBParams) :: CPin
    Type(TTimeSources) :: ThisSourcesin
    logical :: full_bessel_integrationin, do_bispectrumin
    type(datastatebessel) :: datasbin
    double precision, external :: staterofchi

    nu=IV%q*datasbin%s_curvature_radius
    sixpibynu  = 6._dl*3.1415926535897932384626433832795_dl/nu

    if (datasbin%s_closed) then
        if (nu<20 .or. datasbin%s_tau0/datasbin%s_curvature_radius+sixpibynu > const_pi/2) then
            llmax=nint(nu)-1
        else
            ! beeing if flat shoud be chi itslef
            !print *, "no Chi :", Statein%tau0/Statein%curvature_radius + sixpibynu
            !print *, "   Chi :", Statein%rofChi(Statein%tau0/Statein%curvature_radius + sixpibynu)
            !llmax=nint(nu*Statein%rofChi(Statein%tau0/Statein%curvature_radius + sixpibynu))
            ! check if it is correct 
            !print *, "orig llmax: ", nint(nu*State%rofChi(State%tau0/State%curvature_radius + sixpibynu))
            llmax=nint(nu*staterofchi(datasbin%s_flat, datasbin%s_closed, &
                datasbin%s_tau0/datasbin%s_curvature_radius + sixpibynu))
            !print *, "new llmax: ", llmax
            !llmax=nint(nu*(datasbin%s_tau0/datasbin%s_curvature_radius + sixpibynu))
            llmax=min(llmax,nint(nu)-1)  !nu >= l+1
        end if
    else
        llmax=nint(nu*datasbin%s_chi0)
        if (llmax<15) then
            llmax=17 !AL Sept2010 changed from 15 to get l=16 smooth
        else
            !print *, "no Chi: ", Statein%tau0/Statein%curvature_radius + sixpibynu
            !print *, "   Chi: ", Statein%rofChi(Statein%tau0/Statein%curvature_radius + sixpibynu)
            !llmax=nint(nu*Statein%rofChi(Statein%tau0/Statein%curvature_radius + sixpibynu))
            !print *, "orig llmax: ", nint(nu*State%rofChi(State%tau0/State%curvature_radius + sixpibynu))
            llmax = nint(nu*staterofchi (datasbin%s_flat, datasbin%s_closed, & 
                     datasbin%s_tau0/datasbin%s_curvature_radius + sixpibynu))
            !print *, "new llmax: ", llmax
            !llmax=nint(nu*(datasbin%s_tau0/datasbin%s_curvature_radius + sixpibynu))
        end if
    end if

    if (datasbin%s_flat) then
        call DoFlatIntegration(IV,ThisCT, llmax, CPin, ThisSourcesin, &
          full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
          datasbin,xlimfracin,xlimminin,ajlin,ajlprin)
    else
        print * , "not yet fully ported"
        stop
        !do j=1,ThisCT%ls%nl
        !    ll=ThisCT%ls%l(j)
        !    if (ll>llmax) exit
        !    call IntegrateSourcesBessels(IV,ThisCT,j,ll,nu,Statein,CPin,ThisSourcesin)
        !end do !j loop
    end if

end subroutine DoSourceIntegration

    !flat source integration
subroutine DoFlatIntegration(IV, ThisCT, llmax, CPin, ThisSourcesin, &
    full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
    datasbin,xlimfracin,xlimminin,ajlin,ajlprin)
#ifdef USEACC
    !$ACC ROUTINE
#endif

    use CAMBmain
    use precision
    use model
    use results

    implicit none
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT 
    Type(TTimeSources) :: ThisSourcesin
    integer llmax
    integer j
    logical DoInt
    real(dl) xlimfracin, xlimminin
    real(dl), dimension(:,:), allocatable, intent(inout) :: ajlin, ajlprin
    real(dl) xlim,xlmax1
    real(dl) tmin, tmax
    real(dl) a2, J_l, aa(IV%SourceSteps), fac(IV%SourceSteps)
    real(dl) xf, sums(ThisSourcesin%SourceNum)
    real(dl) qmax_int
    integer bes_ix,n, bes_index(IV%SourceSteps)
    integer custom_source_off, s_ix
    integer nwin, max_bessels_l_indexin
    real(dl) :: BessIntBoost
    Type(CAMBParams) :: CPin
    logical :: full_bessel_integrationin, do_bispectrumin
    type(datastatebessel) :: datasbin

    INTERFACE
        FUNCTION statbesseindexof (count, R, npoints, Highest, tau)
#ifdef USEACC
            !$ACC ROUTINE 
#endif
            USE RangeUtils
            INTEGER :: statbesseindexof

            INTEGER, INTENT(IN) :: count            
            TYPE(TRange), INTENT(IN) :: R(count)    
            INTEGER, INTENT(IN) :: npoints         
            DOUBLE PRECISION, INTENT(IN) :: Highest 
            DOUBLE PRECISION, INTENT(IN) :: tau     
        END FUNCTION statbesseindexof
    END INTERFACE

    INTERFACE
        FUNCTION UseLimberGPU(l, CPin)
#ifdef USEACC
            !$ACC ROUTINE
#endif
            use model
            LOGICAL :: UseLimberGPU
            INTEGER :: l
            TYPE(CAMBParams) :: CPin
        END FUNCTION UseLimberGPU
    END INTERFACE

#ifdef COMPARISON
    integer :: tocompare
#endif
    integer :: startloopidx, endloopidx

    BessIntBoost = CPin%Accuracy%AccuracyBoost*CPin%Accuracy%BessIntBoost
    custom_source_off = datasbin%s_num_redshiftwindows + datasbin%s_num_extra_redshiftwindows + 4

    !     Find the position in the xx table for the x correponding to each
    !     timestep

    do j=1,IV%SourceSteps !Precompute arrays for this k
        xf=abs(IV%q*(datasbin%s_tau0-datasbin%s_points(j)))
#ifdef COMPARISON
        ! in case need to use a statein as input
        !bes_index(j)=BessRanges%IndexOf(xf)
        !tocompare = bes_index(j)
#endif
        bes_index(j)=statbesseindexof (datasbin%b_count, &
            datasbin%b_R, datasbin%b_npoints, datasbin%b_Highest, xf)
#ifdef COMPARISON        
        !if (tocompare /= bes_index(j)) then
        !    print *, "Error in Bessel index: ", tocompare, bes_index(j)
        !    stop
        !end if
#endif
        bes_ix= bes_index(j)
        fac(j)=datasbin%b_points(bes_ix+1)-datasbin%b_points(bes_ix)
        aa(j)=(datasbin%b_points(bes_ix+1)-xf)/fac(j)
        fac(j)=fac(j)**2*aa(j)/6
    end do
    !print *, "Done first indexof"

    do j=1,max_bessels_l_indexin
        if (ThisCT%ls%l(j) > llmax) return
        xlim=xlimfracin*ThisCT%ls%l(j)
        xlim=max(xlim,xlimminin)
        xlim=ThisCT%ls%l(j)-xlim
        if (full_bessel_integrationin .or. do_bispectrumin) then
            tmin = datasbin%s_points(2)
        else
            xlmax1=80*ThisCT%ls%l(j)*BessIntBoost
            if (datasbin%s_num_redshiftwindows>0 .and. CPin%WantScalars) then
                xlmax1=80*ThisCT%ls%l(j)*8*BessIntBoost !Have to be careful if sharp spikes due to late time sources
            end if
            tmin=datasbin%s_tau0-xlmax1/IV%q
            tmin=max(datasbin%s_points(2),tmin)
        end if
        tmax=datasbin%s_tau0-xlim/IV%q
        tmax=min(datasbin%s_tau0,tmax)
        tmin=max(datasbin%s_points(2),tmin)
        if (.not. CPin%Want_CMB .and. .not. CPin%Want_CMB_lensing) &
            tmin = max(tmin, datasbin%s_tau_start_redshiftwindows)


        if (tmax < datasbin%s_points(2)) exit
        sums = 0

        !As long as we sample the source well enough, it is sufficient to
        !interpolate the Bessel functions only

        if (ThisSourcesin%SourceNum==2) then
            !This is the innermost loop, so we separate the no lensing scalar case to optimize it
            startloopidx = statbesseindexof (datasbin%s_count, datasbin%s_R, &
                datasbin%s_npoints, datasbin%s_Highest, tmin)
            endloopidx = min(IV%SourceSteps,statbesseindexof (datasbin%s_count, &
                datasbin%s_R, datasbin%s_npoints, datasbin%s_Highest, tmax)) 
#ifdef COMPARISON
            ! in case need to use a statein as input
            !tocompare = State%TimeSteps%IndexOf(tmin)
            !if (tocompare /= startloopidx) then
            !    print *, "Error in State index: ", tocompare, startloopidx
            !    stop
            !end if
            !tocompare = min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
            !if (tocompare /= endloopidx) then
            !    print *, "Error in State index: ", tocompare, endloopidx
            !    stop
            !end if
#endif
            do n= startloopidx,endloopidx
                a2=aa(n)
                bes_ix=bes_index(n)

                J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
                    *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline

                J_l = J_l*datasbin%s_dpoints(n)
                sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                sums(2) = sums(2) + IV%Source_q(n,2)*J_l
            end do
        else
            qmax_int= max(850,ThisCT%ls%l(j))*3*BessIntBoost/datasbin%s_tau0*1.2
            DoInt = .not. CPin%WantScalars .or. IV%q < qmax_int
            !Do integral if any useful contribution to the CMB, or large scale effects

            if (DoInt) then
                if (CPin%CustomSources%num_custom_sources==0 .and. datasbin%s_num_redshiftwindows==0) then
                    startloopidx = statbesseindexof (datasbin%s_count, datasbin%s_R, &
                       datasbin%s_npoints, datasbin%s_Highest, tmin)
                    endloopidx = min(IV%SourceSteps,statbesseindexof (datasbin%s_count, &
                       datasbin%s_R, datasbin%s_npoints, datasbin%s_Highest, tmax))
#ifdef COMPARISON
                    ! in case need to use a statein as input
                    !tocompare = State%TimeSteps%IndexOf(tmin)
                    !if (tocompare /= startloopidx) then
                    !    print *, "Error in State index: ", tocompare, startloopidx
                    !    stop
                    !end if
                    !tocompare = min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                    !if (tocompare /= endloopidx) then
                    !    print *, "Error in State index: ", tocompare, endloopidx
                    !    stop
                    !end if
#endif
                    do n=startloopidx,endloopidx
                       !Full Bessel integration
                       a2=aa(n)
                       bes_ix=bes_index(n)

                       J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
                           *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline
                       J_l = J_l*datasbin%s_dpoints(n)

                       !The unwrapped form is faster
                       sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                       sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                       sums(3) = sums(3) + IV%Source_q(n,3)*J_l
                   end do
                else
                    if (datasbin%s_num_redshiftwindows>0) then
                        nwin = statbesseindexof(datasbin%s_count, datasbin%s_R, &
                            datasbin%s_npoints, datasbin%s_Highest, &
                            datasbin%s_tau_start_redshiftwindows)
#ifdef COMPARISON  
                        ! in case need to use a statein as input
                        !tocompare = State%TimeSteps%IndexOf(State%ThermoData%tau_start_redshiftwindows)
                        !if (tocompare /= nwin) then
                        !    print *, "Error in State index: ", tocompare, nwin
                        !    stop
                        !end if     
#endif
                    else
                        !nwin = State%TimeSteps%npoints+1
                        nwin = datasbin%s_npoints+1
                    end if
                    if (CPin%CustomSources%num_custom_sources==0) then
                        startloopidx = statbesseindexof (datasbin%s_count, datasbin%s_R, &
                            datasbin%s_npoints, datasbin%s_Highest, tmin)   
                        endloopidx = min(IV%SourceSteps,statbesseindexof (datasbin%s_count, &
                            datasbin%s_R, datasbin%s_npoints, datasbin%s_Highest, tmax))
#ifdef COMPARISON
                        ! in case need to use a statein as input
                        !tocompare = State%TimeSteps%IndexOf(tmin)
                        !if (tocompare /= startloopidx) then
                        !    print *, "Error in State index: ", tocompare, startloopidx
                        !    stop
                        !end if
                        !tocompare = min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                        !if (tocompare /= endloopidx) then
                        !    print *, "Error in State index: ", tocompare, endloopidx
                        !    stop
                        !end if
#endif                        
                        do n= startloopidx,endloopidx
                           a2=aa(n)
                           bes_ix=bes_index(n)

                           J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
                               *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline
                           J_l = J_l*datasbin%s_dpoints(n)

                           !The unwrapped form is faster
                           sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                           sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                           sums(3) = sums(3) + IV%Source_q(n,3)*J_l
                           if (n >= nwin) then
                               do s_ix = 4, ThisSourcesin%SourceNum
                                   sums(s_ix) = sums(s_ix) + IV%Source_q(n,s_ix)*J_l
                               end do
                           end if
                       end do
                    else
                        startloopidx = statbesseindexof (datasbin%s_count, datasbin%s_R, &
                            datasbin%s_npoints, datasbin%s_Highest, tmin)
                        endloopidx = min(IV%SourceSteps,statbesseindexof (datasbin%s_count, &
                            datasbin%s_R, datasbin%s_npoints, datasbin%s_Highest, tmax))
#ifdef COMPARISON
                        ! in case need to use a statein as input
                        !tocompare = State%TimeSteps%IndexOf(tmin)
                        !if (tocompare /= startloopidx) then
                        !    print *, "Error in State index: ", tocompare, startloopidx
                        !    stop
                        !end if
                        !tocompare = min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                        !if (tocompare /= endloopidx) then
                        !    print *, "Error in State index: ", tocompare, endloopidx
                        !    stop
                        !end if
#endif
                        do n=startloopidx,endloopidx
                           a2=aa(n)
                           bes_ix=bes_index(n)

                           J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
                               *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline
                           J_l = J_l*datasbin%s_dpoints(n)

                           !The unwrapped form is faster
                           sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                           sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                           sums(3) = sums(3) + IV%Source_q(n,3)*J_l
                           sums(custom_source_off) = sums(custom_source_off) +  IV%Source_q(n,custom_source_off)*J_l
                           if (n >= nwin) then
                               do s_ix = 4, ThisSourcesin%NonCustomSourceNum
                                   sums(s_ix) = sums(s_ix) + IV%Source_q(n,s_ix)*J_l
                               end do
                           end if
                           do s_ix = custom_source_off+1, custom_source_off+CPin%CustomSources%num_custom_sources -1
                               sums(s_ix) = sums(s_ix)  + IV%Source_q(n,s_ix)*J_l
                           end do
                       end do
                    end if
                end if
            end if
            if (.not. DoInt .or. UseLimberGPU(ThisCT%ls%l(j), CPin) .and. CPin%WantScalars) then
                !Limber approximation for small scale lensing (better than poor version of above integral)
                xf = datasbin%s_tau0-(ThisCT%ls%l(j)+0.5_dl)/IV%q
                if (xf < datasbin%s_highest .and. xf > datasbin%s_lowest) then
                    n=statbesseindexof (datasbin%s_count, datasbin%s_R, &
                        datasbin%s_npoints, datasbin%s_Highest, xf)
#ifdef COMPARISON
                    ! in case need to use a statein as input
                    !tocompare=State%TimeSteps%IndexOf(xf)
                    !if (tocompare /= n) then
                    !    print *, "Error in State index: ", tocompare, n
                    !    stop
                    !end if
#endif
                    !n=statindexof(xf)
                    xf= (xf-datasbin%s_points(n))/(datasbin%s_points(n+1)-datasbin%s_points(n))
                    sums(3) = (IV%Source_q(n,3)*(1-xf) + xf*IV%Source_q(n+1,3))*&
                        sqrt(const_pi/2/(ThisCT%ls%l(j)+0.5_dl))/IV%q
                else
                    sums(3)=0
                end if
            end if
            if (.not. DoInt .and. ThisSourcesin%NonCustomSourceNum>3) then
                if (any(ThisCT%limber_l_min(4:ThisSourcesin%NonCustomSourceNum)==0 .or. &
                    ThisCT%limber_l_min(4:ThisSourcesin%NonCustomSourceNum) > j)) then
                    !When CMB does not need integral but other sources do
                    startloopidx = statbesseindexof (datasbin%s_count, datasbin%s_R, &
                        datasbin%s_npoints, datasbin%s_Highest, &
                        datasbin%s_tau_start_redshiftwindows)
                    endloopidx = min(IV%SourceSteps,statbesseindexof (datasbin%s_count, &
                        datasbin%s_R, datasbin%s_npoints, datasbin%s_Highest, tmax))
#ifdef COMPARISON
                    ! in case need to use a statein as input
                    !tocompare = State%TimeSteps%IndexOf(State%ThermoData%tau_start_redshiftwindows)
                    !if (tocompare /= startloopidx) then
                    !    print *, "Error in State index: ", tocompare, startloopidx
                    !    stop
                    !end if
                    !tocompare = min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                    !if (tocompare /= endloopidx) then
                    !    print *, "Error in State index: ", tocompare, endloopidx
                    !    stop
                    !end if
#endif
                    do n=startloopidx,endloopidx
                        !Full Bessel integration
                        a2 = aa(n)
                        bes_ix = bes_index(n)

                        J_l = a2 * ajlin(bes_ix, j) + (1 - a2) * (ajlin(bes_ix + 1, j) -&
                            ((a2 + 1) * ajlprin(bes_ix, j) + (2 - a2) * &
                            ajlprin(bes_ix + 1, j)) * fac(n)) !cubic spline
                        J_l = J_l * datasbin%s_dpoints(n)

                        sums(4) = sums(4) + IV%Source_q(n, 4) * J_l
                        do s_ix = 5, ThisSourcesin%NonCustomSourceNum
                            sums(s_ix) = sums(s_ix) + IV%Source_q(n, s_ix) * J_l
                        end do
                    end do
                end if
            end if
        end if

        ThisCT%Delta_p_l_k(:,j,IV%q_ix) = ThisCT%Delta_p_l_k(:,j,IV%q_ix) + sums
    end do

end subroutine DoFlatIntegration
