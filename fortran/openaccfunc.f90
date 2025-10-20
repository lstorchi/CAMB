! START OPENACC 

#define  IVSQROWS 3600
#define  IVSQCOLS 4

#ifndef ONLYFLAT
subroutine spline_def_local (x,y,n,d2)
    !Low-level initialize spline arrays with default boundary conditions 
    integer, intent(in) :: n
    real(dl), intent(in) :: x(n), y(n)
    real(dl), intent(out) :: d2(n)
    real(dl) ::  d11, d1n
    real(dl) xp,qn,sig,un,xxdiv,u(n-1),d1l,d1r
    real(dl), parameter :: LOCALSPLINE_DANGLE=1.d30
    integer i

    d11 = LOCALSPLINE_DANGLE 
    d1n = LOCALSPLINE_DANGLE

    d1r= (y(2)-y(1))/(x(2)-x(1))
    if (d11==SPLINE_DANGLE) then
        d2(1)=0.d0
        u(1)=0.d0
    else
        d2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*(d1r-d11)
    endif

    do i=2,n-1
        d1l=d1r
        d1r=(y(i+1)-y(i))/(x(i+1)-x(i))
        xxdiv=1.d0/(x(i+1)-x(i-1))
        sig=(x(i)-x(i-1))*xxdiv
        xp=1.d0/(sig*d2(i-1)+2.d0)

        d2(i)=(sig-1.d0)*xp

        u(i)=(6.d0*(d1r-d1l)*xxdiv-sig*u(i-1))*xp
    end do
    d1l=d1r

    if (d1n==LOCALSPLINE_DANGLE) then
        qn=0.d0
        un=0.d0
    else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(d1n-d1l)
    endif

    d2(n)=(un-qn*u(n-1))/(qn*d2(n-1)+1.d0)
    do i=n-1,1,-1
        d2(i)=d2(i)*d2(i+1)+u(i)
    end do
end subroutine spline_def_local
#endif

function statbesseindexof (count, R, npoints, Highest, tau) 
#ifdef USEACC
!$acc routine seq
#endif
    use RangeUtils
    !statein%TimeSteps%IndexOf  RangeUtils.f90 procedure :: IndexOf => TRanges_IndexOf
    ! to test it compare respect to State.IndexOf 
    integer :: statbesseindexof
    double precision, value, intent(in) :: tau
    integer , value, intent(in) :: count
    type(TRange), intent(in) :: R(count)
    integer, intent(in) :: npoints
    double precision, intent(in) :: Highest
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

    statbesseindexof=1
    do i=1, count
        if (tau < R(i)%High .and. tau >= R(i)%Low) then
            if (R(i)%IsLog) then
                statbesseindexof = R(i)%start_index + int(log(tau / R(i)%Low) / R(i)%delta)
            else
                statbesseindexof = R(i)%start_index + int((tau - R(i)%Low) / R(i)%delta)
            end if
            return
        end if
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

function UseLimberGPU(l,  datasb)
#ifdef USEACC
!$ACC ROUTINE SEQ
#endif
        !Calculate lensing potential power using Limber rather than j_l integration
        !even when sources calculated as part of temperature calculation
        !(Limber better on small scales unless step sizes made much smaller)
        !This affects speed, esp. of non-flat case
        !use model

        logical :: UseLimberGPU
        integer l
        Type(datastatebessel) :: datasb
    
        !note increasing non-limber is not neccessarily more accurate unless AccuracyBoost much higher
        !use **0.5 to at least give some sensitivity to Limber effects
        !Could be lower but care with phi-T correlation at lower L
        if (datasb%cp_st_limber_windows) then
            UseLimberGPU = l >= datasb%cp_st_limber_phi_lmin
        else
            UseLimberGPU = l > 400 * (datasb%cp_accuracy_boost * &
                datasb%cp_accuracy_liber_boost)** 0.5
        end if
    
end function UseLimberGPU

subroutine SourceToTransfers(datasb, &
    ThisCT, q_ix,  ThisSourcesin, ScaledSrcin, &
    ddScaledSrcin, max_etak_tensorin, max_etak_vectorin, &
    WantLateTimein, max_etak_scalarin, full_bessel_integrationin, &
    do_bispectrumin, max_bessels_l_indexin, &
    xlimfracin, xlimminin, ajlin, ajlprin, DebugEvolutionin) 
    !IVSource_q)
 
#ifdef USEACC
!$acc routine vector 
!acc routine 
#endif
!    use CAMBmain
!    use results
!    use RangeUtils

    implicit none

    real(dl) :: xlimfracin, xlimminin
    real(dl), dimension(:,:), allocatable :: ajlin, ajlprin
    type(ClTransferData), target :: ThisCT 
    Type(TTimeSources) :: ThisSourcesin
    integer :: q_ix, max_bessels_l_indexin
    real(dl), dimension(:,:,:) :: ScaledSrcin
    real(dl), dimension(:,:,:) :: ddScaledSrcin
    real(dl) :: max_etak_tensorin, max_etak_vectorin, max_etak_scalarin
    logical :: WantLateTimein
    logical :: full_bessel_integrationin, do_bispectrumin, DebugEvolutionin
    type(datastatebessel) :: datasb   
    type(PrivateIdxs) :: privateindexes
    !double precision , allocatable, dimension(:,:) :: IVSource_q
    double precision , dimension(IVSQROWS,IVSQCOLS) :: IVSource_q
    !call IntegrationVars_Init(IV, datasb)
    ! to avoid a call 

    IVSource_q(1,:)=0
    IVSource_q(datasb%s_npoints,:) = 0
    IVSource_q(datasb%s_npoints-1,:) = 0

    privateindexes%iv_q_ix = q_ix
    privateindexes%iv_q = ThisCT%q%points(q_ix)
    privateindexes%iv_dq = ThisCT%q%dpoints(q_ix)

    call InterpolateSources(ThisSourcesin, ScaledSrcin, ddScaledSrcin, &
      max_etak_tensorin, max_etak_vectorin, WantLateTimein, max_etak_scalarin, &
      datasb, DebugEvolutionin, privateindexes, IVSource_q)

    call DoSourceIntegration(ThisCT, ThisSourcesin, &
            full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
            datasb,xlimfracin,xlimminin,ajlin,ajlprin, &
            privateindexes, IVSource_q)

end subroutine SourceToTransfers


subroutine InterpolateSources(ThisSourcesin, ScaledSrcin, &
    ddScaledSrcin, max_etak_tensorin, max_etak_vectorin, &
    WantLateTimein, max_etak_scalarin, datasb, DebugEvolutionin, &
    privateindexes, IVSource_q)
#ifdef USEACC
!$acc routine vector
#endif

!    use CAMBmain
!    use results

    implicit none
    integer i,khi,klo, step
    real(dl) xf,b0,ho,a0,ho2o6,a03,b03
    Type(TTimeSources) :: ThisSourcesin
    real(dl), dimension(:,:,:) :: ScaledSrcin
    real(dl), dimension(:,:,:) :: ddScaledSrcin
    real(dl) :: max_etak_tensorin, max_etak_vectorin, max_etak_scalarin
    logical :: WantLateTimein
    type(datastatebessel) :: datasb
    logical :: DebugEvolutionin
    type(PrivateIdxs) :: privateindexes
    integer :: ixunit, local_step
    !double precision , allocatable, dimension(:,:) :: IVSource_q
    double precision , dimension(IVSQROWS,IVSQCOLS) :: IVSource_q

#ifndef USEACC
    print *, "datasb%cp_want_tensors: ", datasb%cp_want_tensors
    print *, "datasb%cp_want_vectors: ", datasb%cp_want_vectors
    print *, "datasb%cp_want_scalars: ", datasb%cp_want_scalars
#endif
    
    !character(len=30) :: filename
    !     finding position of k in table Evolve_q to do the interpolation.

    !Can't use the following in closed case because regions are not set up (only points)
    !           klo = min(ThisSourcesin%Evolve_q%npoints-1,ThisSourcesin%Evolve_q%IndexOf(privateindexes%iv_q))
    !This is a bit inefficient, but thread safe
    klo=1
    do while ((privateindexes%iv_q > ThisSourcesin%Evolve_q%points(klo+1)).and.&
        (klo < (ThisSourcesin%Evolve_q%npoints-1)))
        klo=klo+1
    end do

    khi=klo+1

    ho=ThisSourcesin%Evolve_q%points(khi)-ThisSourcesin%Evolve_q%points(klo)
    a0=(ThisSourcesin%Evolve_q%points(khi)-privateindexes%iv_q)/ho
    b0=(privateindexes%iv_q-ThisSourcesin%Evolve_q%points(klo))/ho
    ho2o6 = ho**2/6
    a03=(a0**3-a0)
    b03=(b0**3-b0)
    ixunit = privateindexes%iv_q_ix
    privateindexes%iv_sourcessteps = 0

    !Interpolating the source as a function of time for the present
    !wavelength.
    local_step = 2
    step=2
    !$acc loop vector reduction(max:local_step)
    do i=2, datasb%s_npoints
        xf=privateindexes%iv_q*(datasb%s_tau0-datasb%s_points(i))
        
        if (datasb%cp_want_tensors) then
            if (privateindexes%iv_q*datasb%s_points(i) < max_etak_tensorin.and. xf > 1.e-8_dl) then
                IVSource_q(i,:) =a0*ScaledSrcin(klo,:,i)+&
                    b0*ScaledSrcin(khi,:,i)+(a03 *ddScaledSrcin(klo,:,i)+ &
                    b03*ddScaledSrcin(khi,:,i)) *ho2o6
                local_step = i
            else
                IVSource_q(i,:) = 0.0_dl
            end if
        end if

        if (datasb%cp_want_vectors) then
            if (privateindexes%iv_q*datasb%s_points(i) < max_etak_vectorin.and. &
              xf > 1.e-8_dl) then
                IVSource_q(i,:) = a0*ScaledSrcin(klo,:,i) + & 
                  b0*ScaledSrcin(khi,:,i)+(a03 * ddScaledSrcin(klo,:,i)+ &
                  b03*ddScaledSrcin(khi,:,i)) *ho2o6
                local_step = i
            else
                IVSource_q(i,:) = 0.0_dl
            end if
        end if

        if (datasb%cp_want_scalars) then
            if ((DebugEvolutionin .or. WantLateTimein .or. &
                privateindexes%iv_q*datasb%s_points(i) < max_etak_scalarin) &
                .and. xf > 1.e-8_dl) then
                IVSource_q(i,:) = a0 * ScaledSrcin(klo,:,i) +  & 
                  b0 * ScaledSrcin(khi,:,i) + (a03*ddScaledSrcin(klo,:,i) + &
                  b03 * ddScaledSrcin(khi,:,i)) * ho2o6
                local_step = i
            else
                IVSource_q(i,:) = 0.0_dl
            end if
        end if
    end do
    !$acc end loop
    
    step = local_step 
    privateindexes%iv_sourcessteps = step

#ifndef ONLYFLAT  
    ! only non flat 
    if (.not.datasb%s_flat) then
        do i=1, datasb%ttsources_sourcenum
            call spline_def_local(datasb%s_points,IVSource_q(:,i),datasb%s_npoints,&
                IVddSource_q(:,i))
        end do
    end if
#endif

end subroutine InterpolateSources

subroutine DoSourceIntegration(ThisCT, ThisSourcesin, &
    full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
    datasb, xlimfracin, xlimminin, ajlin, ajlprin, privateindexes, &
    IVSource_q) !for particular wave number q
#ifdef USEACC
!$acc routine vector
#endif

!    use CAMBmain
!    use precision
!    use model
!    use results

    Type(ClTransferData) :: ThisCT    
    real(dl), dimension(:,:), allocatable, intent(inout) :: ajlin, ajlprin
    real(dl) xlimfracin, xlimminin
    integer j,ll,llmax, max_bessels_l_indexin 
    real(dl) nu
    real(dl) :: sixpibynu
    Type(TTimeSources) :: ThisSourcesin
    logical :: full_bessel_integrationin, do_bispectrumin
    type(datastatebessel) :: datasb
    type(PrivateIdxs) :: privateindexes
    !double precision , allocatable, dimension(:,:) :: IVSource_q
    double precision , dimension(IVSQROWS,IVSQCOLS) :: IVSource_q

    nu=privateindexes%iv_q*datasb%s_curvature_radius
    sixpibynu  = 6._dl*3.1415926535897932384626433832795_dl/nu

#ifndef USEACC
   print *, " datasb%s_closed: ", datasb%s_closed
   print *, "   datasb%s_flat: ", datasb%s_flat
#endif

    if (datasb%s_closed) then
        if (nu<20 .or. datasb%s_tau0/datasb%s_curvature_radius+sixpibynu > const_pi/2) then
            llmax=nint(nu)-1
        else
            ! beeing if flat shoud be chi itslef
            !print *, "no Chi :", Statein%tau0/Statein%curvature_radius + sixpibynu
            !print *, "   Chi :", Statein%rofChi(Statein%tau0/Statein%curvature_radius + sixpibynu)
            !llmax=nint(nu*Statein%rofChi(Statein%tau0/Statein%curvature_radius + sixpibynu))
            ! check if it is correct 
            !print *, "orig llmax: ",iv_sourcessteps nint(nu*State%rofChi(State%tau0/State%curvature_radius + sixpibynu))
            llmax=nint(nu*staterofchi(datasb%s_flat, datasb%s_closed, &
                datasb%s_tau0/datasb%s_curvature_radius + sixpibynu))
            !print *, "new llmax: ", llmax
            !llmax=nint(nu*(datasb%s_tau0/datasb%s_curvature_radius + sixpibynu))
            llmax=min(llmax,nint(nu)-1)  !nu >= l+1
        end if
    else
        llmax=nint(nu*datasb%s_chi0)
        if (llmax<15) then
            llmax=17 !AL Sept2010 changed from 15 to get l=16 smooth
        else
            !print *, "no Chi: ", Statein%tau0/Statein%curvature_radius + sixpibynu
            !print *, "   Chi: ", Statein%rofChi(Statein%tau0/Statein%curvature_radius + sixpibynu)
            !llmax=nint(nu*Statein%rofChi(Statein%tau0/Statein%curvature_radius + sixpibynu))
            !print *, "orig llmax: ", nint(nu*State%rofChi(State%tau0/State%curvature_radius + sixpibynu))
            llmax = nint(nu*staterofchi (datasb%s_flat, datasb%s_closed, & 
                     datasb%s_tau0/datasb%s_curvature_radius + sixpibynu))
            !print *, "new llmax: ", llmax
            !llmax=nint(nu*(datasb%s_tau0/datasb%s_curvature_radius + sixpibynu))
        end if
    end if

    if (datasb%s_flat) then
        call DoFlatIntegration(ThisCT, llmax, ThisSourcesin, &
          full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
          datasb,xlimfracin,xlimminin,ajlin,ajlprin,privateindexes, IVSource_q)
    else
#ifdef ONLYFLAT
        print * , "not yet fully ported"
        stop
#else
        do j=1,ThisCT%ls%nl
            ll=ThisCT%ls%l(j)
            if (ll>llmax) exit
            call IntegrateSourcesBessels(IV,ThisCT,j,ll,nu,Statein,ThisSourcesin)
        end do !j loop
#endif
    end if

end subroutine DoSourceIntegration


subroutine DoFlatIntegration(ThisCT, llmax, ThisSourcesin, &
    full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
    datasb, xlimfracin, xlimminin, ajlin, ajlprin, privateindexes, &
    IVSource_q)
#ifdef USEACC
!$acc routine vector
!acc routine
#endif

!    use CAMBmain
!    use precision
!    use model
!    use results
#ifndef USEACC
    use omp_lib
#endif
    implicit none

    ! input 
    Type(ClTransferData) :: ThisCT 
    integer llmax
    Type(TTimeSources) :: ThisSourcesin
    logical :: full_bessel_integrationin, do_bispectrumin
    integer :: max_bessels_l_indexin
    type(datastatebessel) :: datasb
    real(dl) xlimfracin, xlimminin
    real(dl), dimension(:,:), allocatable, intent(inout) :: ajlin, ajlprin
    type(PrivateIdxs) :: privateindexes
    !double precision , allocatable, dimension(:,:) :: IVSource_q
    double precision , dimension(IVSQROWS,IVSQCOLS) :: IVSource_q

    ! local vars    
    integer j
    logical DoInt
    real(dl) xlim,xlmax1
    real(dl) tmin, tmax
    real(dl) a2, J_l, aa(privateindexes%iv_sourcessteps), fac(privateindexes%iv_sourcessteps)
    real(dl) xf, sums(datasb%ttsources_sourcenum)
    real(dl) qmax_int
    integer bes_ix,n, bes_index(privateindexes%iv_sourcessteps)
    integer custom_source_off, s_ix
    integer nwin
    real(dl) :: BessIntBoost
    real(dl) :: temp_sum1, temp_sum2, temp_sum3, temp_sum4, temp_sum5, &
        temp_sum6, temp_sum7, temp_sum8, temp_sum9, temp_sum10
#ifndef USEACC
    integer :: omp_thread_num, thread_id
#endif

    !integer :: tocompare
    integer :: startloopidx, endloopidx

    if (datasb%ttsources_sourcenum.gt.10) then
       ! nee to increase temp sum or find a different approach
       return
    endif

#ifndef USEACC
    omp_thread_num = omp_get_max_threads()
    if (omp_thread_num > 1) then
        thread_id = omp_get_thread_num()
    end if
#endif

#ifndef USEACC
    print *, "full_bessel_integrationin: ", full_bessel_integrationin
    print *, "   do_bispectrumin: ", do_bispectrumin
    print *, "   datasb%ttsources_sourcenum: ", datasb%ttsources_sourcenum
    print *, "  datasb%s_num_redshiftwindows: ", datasb%s_num_redshiftwindows
    print *, "  datasb%s_num_extra_redshiftwindows: ", datasb%s_num_extra_redshiftwindows
    print *, "  datasb%cp_custom_sources_nam_custom: ", datasb%cp_custom_sources_nam_custom
#endif

    BessIntBoost = datasb%cp_accuracy_boost*datasb%cp_accuracy_bessintboost
    custom_source_off = datasb%s_num_redshiftwindows + datasb%s_num_extra_redshiftwindows + 4

    !     Find the position in the xx table for the x correponding to each
    !     timestep
    !$acc loop vector
    do j=1,privateindexes%iv_sourcessteps !Precompute arrays for this k
        xf=abs(privateindexes%iv_q*(datasb%s_tau0-datasb%s_points(j)))
        ! in case need to use a statein as input
        !tocompare=BessRanges%IndexOf(xf)
        
        bes_index(j)=statbesseindexof (datasb%b_count, &
            datasb%b_R, datasb%b_npoints, datasb%b_Highest, xf)
        bes_ix= bes_index(j)
        !bes_ix=1

        fac(j)=datasb%b_points(bes_ix+1)-datasb%b_points(bes_ix)
        aa(j)=(datasb%b_points(bes_ix+1)-xf)/fac(j)
        fac(j)=fac(j)**2*aa(j)/6
    end do
    !print *, "Done first indexof"

    !$acc loop vector 
    do j=1,max_bessels_l_indexin
        if (ThisCT%ls%l(j) > llmax) return
        xlim=xlimfracin*ThisCT%ls%l(j)
        xlim=max(xlim,xlimminin)
        xlim=ThisCT%ls%l(j)-xlim
        if (full_bessel_integrationin .or. do_bispectrumin) then
            tmin = datasb%s_points(2)
        else
            xlmax1=80*ThisCT%ls%l(j)*BessIntBoost
            if (datasb%s_num_redshiftwindows>0 .and. datasb%cp_want_scalars) then
                xlmax1=80*ThisCT%ls%l(j)*8*BessIntBoost !Have to be careful if sharp spikes due to late time sources
            end if
            tmin=datasb%s_tau0-xlmax1/privateindexes%iv_q
            tmin=max(datasb%s_points(2),tmin)
        end if
        tmax=datasb%s_tau0-xlim/privateindexes%iv_q
        tmax=min(datasb%s_tau0,tmax)
        tmin=max(datasb%s_points(2),tmin)
        if (.not. datasb%cp_want_cmb .and. .not. datasb%cp_want_cmp_lensing) &
            tmin = max(tmin, datasb%s_tau_start_redshiftwindows)

        if (tmax < datasb%s_points(2)) exit
        !sums = 0

        temp_sum1 = 0.0_dl
        temp_sum2 = 0.0_dl
        temp_sum3 = 0.0_dl
        temp_sum4 = 0.0_dl 
        temp_sum5 = 0.0_dl
        temp_sum6 = 0.0_dl
        temp_sum7 = 0.0_dl
        temp_sum8 = 0.0_dl
        temp_sum9 = 0.0_dl
        temp_sum10 = 0.0_dl

        !As long as we sample the source well enough, it is sufficient to
        !interpolate the Bessel functions only

        if (datasb%ttsources_sourcenum==2) then
            !This is the innermost loop, so we separate the no lensing scalar case to optimize it
            startloopidx = statbesseindexof (datasb%s_count, datasb%s_R, &
                datasb%s_npoints, datasb%s_Highest, tmin)
            endloopidx = min(privateindexes%iv_sourcessteps,statbesseindexof (datasb%s_count, &
                datasb%s_R, datasb%s_npoints, datasb%s_Highest, tmax)) 
            ! in case need to use a statein as input
            !tocompare = State%TimeSteps%IndexOf(tmin)
            !tocompare = min(privateindexes%iv_sourcessteps,State%TimeSteps%IndexOf(tmax))
            do n= startloopidx,endloopidx
                a2=aa(n)
                bes_ix=bes_index(n)

                J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
                    *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline

                J_l = J_l*datasb%s_dpoints(n)
                !sums(1) = sums(1) + IVSource_q(n,1)*J_l
                temp_sum1 = temp_sum1 + IVSource_q(n,1)*J_l
                !sums(2) = sums(2) + IVSource_q(n,2)*J_l
                temp_sum2 = temp_sum2 + IVSource_q(n,2)*J_l
            end do
        else
            qmax_int= max(850,ThisCT%ls%l(j))*3*BessIntBoost/datasb%s_tau0*1.2
            DoInt = .not. datasb%cp_want_scalars .or. privateindexes%iv_q < qmax_int
            !Do integral if any useful contribution to the CMB, or large scale effects

            if (DoInt) then
                 if (datasb%cp_custom_sources_nam_custom==0 .and. & 
                    datasb%s_num_redshiftwindows==0) then
                    startloopidx = statbesseindexof (datasb%s_count, datasb%s_R, &
                       datasb%s_npoints, datasb%s_Highest, tmin)
                    endloopidx = min(privateindexes%iv_sourcessteps,statbesseindexof (datasb%s_count, &
                       datasb%s_R, datasb%s_npoints, datasb%s_Highest, tmax))
                    ! in case need to use a statein as input
                    !tocompare = State%TimeSteps%IndexOf(tmin)
                    !tocompare = min(privateindexes%iv_sourcessteps,State%TimeSteps%IndexOf(tmax))
                    do n=startloopidx,endloopidx
                       !Full Bessel integration
                       a2=aa(n)
                       bes_ix=bes_index(n)
                      
                       J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
                           *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline
                       J_l = J_l*datasb%s_dpoints(n)

                       !The unwrapped form is faster
                       !sums(1) = sums(1) + IVSource_q(n,1)*J_l
                       temp_sum1 = temp_sum1 + IVSource_q(n,1)*J_l
                       !sums(2) = sums(2) + IVSource_q(n,2)*J_l
                       temp_sum2 = temp_sum2 + IVSource_q(n,2)*J_l
                       !sums(3) = sums(3) + IVSource_q(n,3)*J_l
                       temp_sum3 = temp_sum3 + IVSource_q(n,3)*J_l
                   end do
                else
                    if (datasb%s_num_redshiftwindows>0) then
                        nwin = statbesseindexof(datasb%s_count, datasb%s_R, &
                            datasb%s_npoints, datasb%s_Highest, &
                            datasb%s_tau_start_redshiftwindows)
                        ! in case need to use a statein as input
                        !tocompare = State%TimeSteps%IndexOf(State%ThermoData%tau_start_redshiftwindows)
                    else
                        !nwin = State%TimeSteps%npoints+1
                        nwin = datasb%s_npoints+1
                    end if
                    if (datasb%cp_custom_sources_nam_custom==0) then
                        startloopidx = statbesseindexof (datasb%s_count, datasb%s_R, &
                            datasb%s_npoints, datasb%s_Highest, tmin)   
                        endloopidx = min(privateindexes%iv_sourcessteps,statbesseindexof (datasb%s_count, &
                            datasb%s_R, datasb%s_npoints, datasb%s_Highest, tmax))
                        ! in case need to use a statein as input
                        !tocompare = State%TimeSteps%IndexOf(tmin)
                        !tocompare = min(privateindexes%iv_sourcessteps,State%TimeSteps%IndexOf(tmax))
                        do n= startloopidx,endloopidx
                           a2=aa(n)
                           bes_ix=bes_index(n)

                           J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
                               *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline
                           J_l = J_l*datasb%s_dpoints(n)

                           !The unwrapped form is faster
                           !sums(1) = sums(1) + IVSource_q(n,1)*J_l
                           temp_sum1 = temp_sum1 + IVSource_q(n,1)*J_l
                           !sums(2) = sums(2) + IVSource_q(n,2)*J_l
                           temp_sum2 = temp_sum2 + IVSource_q(n,2)*J_l
                           !sums(3) = sums(3) + IVSource_q(n,3)*J_l
                           temp_sum3 = temp_sum3 + IVSource_q(n,3)*J_l
                           if (n >= nwin) then
                               !do s_ix = 4, datasb%ttsources_sourcenum
                               !    sums(s_ix) = sums(s_ix) + IVSource_q(n,s_ix)*J_l
                               !end do
                               if (datasb%ttsources_sourcenum.ge.4) &
                                 temp_sum4 = temp_sum4 + IVSource_q(n,4)*J_l
                               if (datasb%ttsources_sourcenum.ge.5) &
                                 temp_sum5 = temp_sum5 + IVSource_q(n,5)*J_l
                               if (datasb%ttsources_sourcenum.ge.6) &
                                 temp_sum6 = temp_sum6 + IVSource_q(n,6)*J_l
                               if (datasb%ttsources_sourcenum.ge.7) &
                                 temp_sum7 = temp_sum7 + IVSource_q(n,7)*J_l
                               if (datasb%ttsources_sourcenum.ge.8) &
                                 temp_sum8 = temp_sum8 + IVSource_q(n,8)*J_l
                               if (datasb%ttsources_sourcenum.ge.9) &
                                 temp_sum9 = temp_sum9 + IVSource_q(n,9)*J_l
                               if (datasb%ttsources_sourcenum.ge.10) &
                                 temp_sum10 = temp_sum10 + IVSource_q(n,10)*J_l
                           end if
                       end do
                    else
                        stop ! see TODO
                        
                        startloopidx = statbesseindexof (datasb%s_count, datasb%s_R, &
                            datasb%s_npoints, datasb%s_Highest, tmin)
                        endloopidx = min(privateindexes%iv_sourcessteps,statbesseindexof (datasb%s_count, &
                            datasb%s_R, datasb%s_npoints, datasb%s_Highest, tmax))
                        ! in case need to use a statein as input
                        !tocompare = State%TimeSteps%IndexOf(tmin)
                        !tocompare = min(privateindexes%iv_sourcessteps,State%TimeSteps%IndexOf(tmax))
                        do n=startloopidx,endloopidx
                           a2=aa(n)
                           bes_ix=bes_index(n)

                           J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
                              *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline
                           J_l = J_l*datasb%s_dpoints(n)

                           !The unwrapped form is faster
                           !sums(1) = sums(1) + IVSource_q(n,1)*J_l
                           temp_sum1 = temp_sum1 + IVSource_q(n,1)*J_l
                           !sums(2) = sums(2) + IVSource_q(n,2)*J_l
                           temp_sum2 = temp_sum2 + IVSource_q(n,1)*J_l
                           !sums(3) = sums(3) + IVSource_q(n,3)*J_l
                           temp_sum3 = temp_sum3 + IVSource_q(n,1)*J_l
                           
                           ! TODO
                           !sums(custom_source_off) = sums(custom_source_off) +  IVSource_q(n,custom_source_off)*J_l
                           !temp_sum15 = temp_sum15 +  IVSource_q(n,custom_source_off)*J_l

                           if (n >= nwin) then
                               !do s_ix = 4, datasb%ttsources_non_custom_sources_num
                               !    sums(s_ix) = sums(s_ix) + IVSource_q(n,s_ix)*J_l
                               !end do
                           end if
                           !do s_ix = custom_source_off+1, custom_source_off+datasb%cp_custom_sources_nam_custom -1
                           !    sums(s_ix) = sums(s_ix)  + IVSource_q(n,s_ix)*J_l
                           !end do


                       end do
                    end if
                end if
            end if
            if (.not. DoInt .or. UseLimberGPU(ThisCT%ls%l(j), datasb) .and. datasb%cp_want_scalars) then
                !Limber approximation for small scale lensing (better than poor version of above integral)
                xf = datasb%s_tau0-(ThisCT%ls%l(j)+0.5_dl)/privateindexes%iv_q
                if (xf < datasb%s_highest .and. xf > datasb%s_lowest) then
                    n=statbesseindexof (datasb%s_count, datasb%s_R, &
                        datasb%s_npoints, datasb%s_Highest, xf)
                    ! in case need to use a statein as input
                    !tocompare=State%TimeSteps%IndexOf(xf)
                    !n=statindexof(xf)
                    xf= (xf-datasb%s_points(n))/(datasb%s_points(n+1)-datasb%s_points(n))
                    !sums(3) = (IVSource_q(n,3)*(1-xf) + xf*IVSource_q(n+1,3))*&
                    !    sqrt(const_pi/2/(ThisCT%ls%l(j)+0.5_dl))/privateindexes%iv_q
                    temp_sum3 = (IVSource_q(n,3)*(1-xf) + xf*IVSource_q(n+1,3))*&
                        sqrt(const_pi/2/(ThisCT%ls%l(j)+0.5_dl))/privateindexes%iv_q
                else
                    !sums(3)=0
                    temp_sum3 = 0
                end if
            end if
            if (.not. DoInt .and. datasb%ttsources_non_custom_sources_num>3) then
                if (any(ThisCT%limber_l_min(4:datasb%ttsources_non_custom_sources_num)==0 .or. &
                    ThisCT%limber_l_min(4:datasb%ttsources_non_custom_sources_num) > j)) then
                    !When CMB does not need integral but other sources do
                    startloopidx = statbesseindexof (datasb%s_count, datasb%s_R, &
                        datasb%s_npoints, datasb%s_Highest, &
                        datasb%s_tau_start_redshiftwindows)
                    endloopidx = min(privateindexes%iv_sourcessteps,statbesseindexof (datasb%s_count, &
                        datasb%s_R, datasb%s_npoints, datasb%s_Highest, tmax))
                    ! in case need to use a statein as input
                    !tocompare = State%TimeSteps%IndexOf(State%ThermoData%tau_start_redshiftwindows)
                    !tocompare = min(privateindexes%iv_sourcessteps,State%TimeSteps%IndexOf(tmax))
                    do n=startloopidx,endloopidx
                        !Full Bessel integration
                        a2 = aa(n)
                        bes_ix = bes_index(n)

                        J_l = a2 * ajlin(bes_ix, j) + (1 - a2) * (ajlin(bes_ix + 1, j) -&
                            ((a2 + 1) * ajlprin(bes_ix, j) + (2 - a2) * &
                            ajlprin(bes_ix + 1, j)) * fac(n)) !cubic spline
                        J_l = J_l * datasb%s_dpoints(n)

                        !sums(4) = sums(4) + IVSource_q(n, 4) * J_l
                        temp_sum4 = temp_sum4 + IVSource_q(n, 4) * J_l 
                        !do s_ix = 5, datasb%ttsources_non_custom_sources_num
                        !    sums(s_ix) = sums(s_ix) + IVSource_q(n, s_ix) * J_l
                        !end do
                        if (datasb%ttsources_sourcenum.ge.4) &
                          temp_sum4 = temp_sum4 + IVSource_q(n,4)*J_l
                        if (datasb%ttsources_sourcenum.ge.5) &
                          temp_sum5 = temp_sum5 + IVSource_q(n,5)*J_l
                        if (datasb%ttsources_sourcenum.ge.6) &
                          temp_sum6 = temp_sum6 + IVSource_q(n,6)*J_l
                        if (datasb%ttsources_sourcenum.ge.7) &
                          temp_sum7 = temp_sum7 + IVSource_q(n,7)*J_l
                        if (datasb%ttsources_sourcenum.ge.8) &
                          temp_sum8 = temp_sum8 + IVSource_q(n,8)*J_l
                        if (datasb%ttsources_sourcenum.ge.9) &
                          temp_sum9 = temp_sum9 + IVSource_q(n,9)*J_l
                        if (datasb%ttsources_sourcenum.ge.10) &
                          temp_sum10 = temp_sum10 + IVSource_q(n,10)*J_l
                    end do
                end if
            end if
        end if

        !ThisCT%Delta_p_l_k(:,j,privateindexes%iv_q_ix) = ThisCT%Delta_p_l_k(:,j,privateindexes%iv_q_ix) + sums
        ThisCT%Delta_p_l_k(1, j, privateindexes%iv_q_ix) = &
          ThisCT%Delta_p_l_k(1, j, privateindexes%iv_q_ix) + temp_sum1
        ThisCT%Delta_p_l_k(2, j, privateindexes%iv_q_ix) = &
          ThisCT%Delta_p_l_k(2, j, privateindexes%iv_q_ix) + temp_sum2
        ThisCT%Delta_p_l_k(3, j, privateindexes%iv_q_ix) = &
          ThisCT%Delta_p_l_k(3, j, privateindexes%iv_q_ix) + temp_sum3
        if (datasb%ttsources_sourcenum.ge.4) &
           ThisCT%Delta_p_l_k(4,j,privateindexes%iv_q_ix) = &
                ThisCT%Delta_p_l_k(4,j,privateindexes%iv_q_ix) + temp_sum4
        if (datasb%ttsources_sourcenum.ge.5) &
           ThisCT%Delta_p_l_k(5,j,privateindexes%iv_q_ix) = &
                ThisCT%Delta_p_l_k(5,j,privateindexes%iv_q_ix) + temp_sum5
        if (datasb%ttsources_sourcenum.ge.6) &
           ThisCT%Delta_p_l_k(6,j,privateindexes%iv_q_ix) = &
                ThisCT%Delta_p_l_k(6,j,privateindexes%iv_q_ix) + temp_sum5
        if (datasb%ttsources_sourcenum.ge.7) &
           ThisCT%Delta_p_l_k(7,j,privateindexes%iv_q_ix) = &
                ThisCT%Delta_p_l_k(7,j,privateindexes%iv_q_ix) + temp_sum5
        if (datasb%ttsources_sourcenum.ge.8) &
           ThisCT%Delta_p_l_k(8,j,privateindexes%iv_q_ix) = &
                ThisCT%Delta_p_l_k(8,j,privateindexes%iv_q_ix) + temp_sum5
        if (datasb%ttsources_sourcenum.ge.9) &
           ThisCT%Delta_p_l_k(9,j,privateindexes%iv_q_ix) = &
                ThisCT%Delta_p_l_k(9,j,privateindexes%iv_q_ix) + temp_sum5
        if (datasb%ttsources_sourcenum.ge.10) &
           ThisCT%Delta_p_l_k(10,j,privateindexes%iv_q_ix) = &
                ThisCT%Delta_p_l_k(10,j,privateindexes%iv_q_ix) + temp_sum5
     end do

end subroutine DoFlatIntegration

! END OPENACC 
