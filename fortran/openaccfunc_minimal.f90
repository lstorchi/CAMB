! START OPENACC

#define  IVSQROWS 3600
#define  IVSQCOLS 3

#define EXTRAVECTOR 1

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
!$acc routine vector
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
!$ACC ROUTINE
#endif
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
!$ACC ROUTINE
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
   integer :: ixunit
   !double precision , allocatable, dimension(:,:) :: IVSource_q
   double precision , dimension(IVSQROWS,IVSQCOLS) :: IVSource_q
#ifdef EXTRAVECTOR
   integer :: local_step
#endif

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

!#ifdef EXTRAVECTOR
!   local_step = 2 
!   step = 2
!   !$acc loop vector reduction(max:local_step)
!   do i=2, datasb%s_npoints
!     xf=privateindexes%iv_q*(datasb%s_tau0-datasb%s_points(i))
!     IVSource_q(i,:) = 0.0_dl ! Initialize to 0 first
!     if (datasb%cp_want_tensors) then
!       if (privateindexes%iv_q*datasb%s_points(i) < max_etak_tensorin .and. &
!         xf > 1.e-8_dl) then
!          IVSource_q(i,:) = a0*ScaledSrcin(klo,:,i)+&
!                           b0*ScaledSrcin(khi,:,i)+(a03 *ddScaledSrcin(klo,:,i)+ &
!                           b03*ddScaledSrcin(khi,:,i)) * ho2o6
!          local_step = i ! Update local_step if condition met
!        end if
!      end if
!      if (datasb%cp_want_scalars) then
!        if ((DebugEvolutionin .or. WantLateTimein .or. &
!            privateindexes%iv_q*datasb%s_points(i) < max_etak_scalarin) &
!            .and. xf > 1.e-8_dl) then
!          IVSource_q(i,:) = a0 * ScaledSrcin(klo,:,i) +  &
!                           b0 * ScaledSrcin(khi,:,i) + (a03*ddScaledSrcin(klo,:,i) + &
!                           b03 * ddScaledSrcin(khi,:,i)) * ho2o6
!          local_step = i ! Update local_step if condition met
!        end if
!      end if
!   end do
!   !$acc end loop 

!   step = local_step ! Assign the final max value to step
!   privateindexes%iv_sourcessteps = step
!#else
   step = 2
   do i=2, datasb%s_npoints
      xf=privateindexes%iv_q*(datasb%s_tau0-datasb%s_points(i))

      if (datasb%cp_want_tensors) then
         if (privateindexes%iv_q*datasb%s_points(i) < max_etak_tensorin.and. xf > 1.e-8_dl) then
            step=i
            IVSource_q(i,:) =a0*ScaledSrcin(klo,:,i)+&
               b0*ScaledSrcin(khi,:,i)+(a03 *ddScaledSrcin(klo,:,i)+ &
               b03*ddScaledSrcin(khi,:,i)) *ho2o6
         else
            IVSource_q(i,:) = 0.0_dl
         end if
      end if

      if (datasb%cp_want_scalars) then
         if ((DebugEvolutionin .or. WantLateTimein .or. &
            privateindexes%iv_q*datasb%s_points(i) < max_etak_scalarin) &
            .and. xf > 1.e-8_dl) then
            step=i
            IVSource_q(i,:) = a0 * ScaledSrcin(klo,:,i) +  &
               b0 * ScaledSrcin(khi,:,i) + (a03*ddScaledSrcin(klo,:,i) + &
               b03 * ddScaledSrcin(khi,:,i)) * ho2o6
         else
            IVSource_q(i,:) = 0.0_dl
         end if
      end if
   end do
   privateindexes%iv_sourcessteps = step
!#endif 

end subroutine InterpolateSources

subroutine DoSourceIntegration(ThisCT, ThisSourcesin, &
   full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
   datasb, xlimfracin, xlimminin, ajlin, ajlprin, privateindexes, &
   IVSource_q) !for particular wave number q
#ifdef USEACC
!$acc routine vector
!acc routine
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

   llmax=nint(nu*datasb%s_chi0)
   if (llmax<15) then
      llmax=17 !AL Sept2010 changed from 15 to get l=16 smooth
   else
      llmax = nint(nu*staterofchi (datasb%s_flat, datasb%s_closed, &
         datasb%s_tau0/datasb%s_curvature_radius + sixpibynu))
   end if

   call DoFlatIntegration(ThisCT, llmax, ThisSourcesin, &
      full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
      datasb,xlimfracin,xlimminin,ajlin,ajlprin,privateindexes, IVSource_q)

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
   real(dl) :: temp_sum1, temp_sum2, temp_sum3

   !integer :: tocompare
   integer :: startloopidx, endloopidx

   BessIntBoost = datasb%cp_accuracy_boost*datasb%cp_accuracy_bessintboost
   custom_source_off = datasb%s_num_redshiftwindows + datasb%s_num_extra_redshiftwindows + 4

   do j=1,privateindexes%iv_sourcessteps !Precompute arrays for this k
      xf=abs(privateindexes%iv_q*(datasb%s_tau0-datasb%s_points(j)))
      bes_index(j)=statbesseindexof (datasb%b_count, &
         datasb%b_R, datasb%b_npoints, datasb%b_Highest, xf)
      bes_ix= bes_index(j)

      fac(j)=datasb%b_points(bes_ix+1)-datasb%b_points(bes_ix)
      aa(j)=(datasb%b_points(bes_ix+1)-xf)/fac(j)
      fac(j)=fac(j)**2*aa(j)/6
   end do

#ifdef EXTRAVECTOR
   do j=1,max_bessels_l_indexin
     if (ThisCT%ls%l(j) > llmax) return
     xlim=xlimfracin*ThisCT%ls%l(j)
     xlim=max(xlim,xlimminin)
     xlim=ThisCT%ls%l(j)-xlim
     xlmax1=80*ThisCT%ls%l(j)*BessIntBoost
     tmin=datasb%s_tau0-xlmax1/privateindexes%iv_q
     tmin=max(datasb%s_points(2),tmin)
     tmax=datasb%s_tau0-xlim/privateindexes%iv_q
     tmax=min(datasb%s_tau0,tmax)
     tmin=max(datasb%s_points(2),tmin)
     if (tmax < datasb%s_points(2)) exit
     
     ! Initialize sums array and temporary scalar sums for each j iteration
     !sums = 0.0_dl
     temp_sum1 = 0.0_dl
     temp_sum2 = 0.0_dl
     temp_sum3 = 0.0_dl
   
     qmax_int= max(850,ThisCT%ls%l(j))*3*BessIntBoost/datasb%s_tau0*1.2
     DoInt = .not. datasb%cp_want_scalars .or. privateindexes%iv_q < qmax_int
   
     if (DoInt) then
       startloopidx = statbesseindexof (datasb%s_count, datasb%s_R, &
           datasb%s_npoints, datasb%s_Highest, tmin)
       endloopidx = min(privateindexes%iv_sourcessteps,statbesseindexof (datasb%s_count, &
           datasb%s_R, datasb%s_npoints, datasb%s_Highest, tmax))
   
       ! Apply reduction to scalar temporaries in the n-loop
       !$acc loop reduction(+:temp_sum1, temp_sum2, temp_sum3) 
       do n=startloopidx,endloopidx
         a2=aa(n)
         bes_ix=bes_index(n)
   
         J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
             *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline
         J_l = J_l*datasb%s_dpoints(n)
   
         temp_sum1 = temp_sum1 + IVSource_q(n,1)*J_l
         temp_sum2 = temp_sum2 + IVSource_q(n,2)*J_l
         temp_sum3 = temp_sum3 + IVSource_q(n,3)*J_l
       end do
       
       ! After the n-loop, update the sums array with the reduced scalar values
       sums(1) = temp_sum1
       sums(2) = temp_sum2
       sums(3) = temp_sum3
     end if
   
     ! This section updates sums(3) based on different logic.
     ! It's a direct assignment, not an accumulation within a parallel loop,
     ! so it should be fine with respect to the *reported* error.
     if (.not. DoInt .or. UseLimberGPU(ThisCT%ls%l(j), datasb) &
         .and. datasb%cp_want_scalars) then
       xf = datasb%s_tau0-(ThisCT%ls%l(j)+0.5_dl)/privateindexes%iv_q
       if (xf < datasb%s_highest .and. xf > datasb%s_lowest) then
         n=statbesseindexof (datasb%s_count, datasb%s_R, &
             datasb%s_npoints, datasb%s_Highest, xf)
         xf= (xf-datasb%s_points(n))/(datasb%s_points(n+1)-datasb%s_points(n))
         sums(3) = (IVSource_q(n,3)*(1-xf) + xf*IVSource_q(n+1,3))*&
             sqrt(const_pi/2/(ThisCT%ls%l(j)+0.5_dl))/privateindexes%iv_q
       else
         sums(3)=0.0_dl
       end if
     end if
   
     ! This final update might need !$acc atomic update if multiple 'j' iterations
     ! (gangs) could write to the same elements of ThisCT%Delta_p_l_k simultaneously.
     ! However, this is a separate issue from the error reported for 'sums(:)'.
     !$acc atomic update
     ThisCT%Delta_p_l_k(1, j, privateindexes%iv_q_ix) = ThisCT%Delta_p_l_k(1, j, privateindexes%iv_q_ix) + sums(1)
     !$acc atomic update
     ThisCT%Delta_p_l_k(2, j, privateindexes%iv_q_ix) = ThisCT%Delta_p_l_k(2, j, privateindexes%iv_q_ix) + sums(2)
     !$acc atomic update
     ThisCT%Delta_p_l_k(3, j, privateindexes%iv_q_ix) = ThisCT%Delta_p_l_k(3, j, privateindexes%iv_q_ix) + sums(3)
     !ThisCT%Delta_p_l_k(:,j,privateindexes%iv_q_ix) = ThisCT%Delta_p_l_k(:,j,privateindexes%iv_q_ix) + sums
   end do
#else
   do j=1,max_bessels_l_indexin
      if (ThisCT%ls%l(j) > llmax) return
      xlim=xlimfracin*ThisCT%ls%l(j)
      xlim=max(xlim,xlimminin)
      xlim=ThisCT%ls%l(j)-xlim
      xlmax1=80*ThisCT%ls%l(j)*BessIntBoost
      tmin=datasb%s_tau0-xlmax1/privateindexes%iv_q
      tmin=max(datasb%s_points(2),tmin)
      tmax=datasb%s_tau0-xlim/privateindexes%iv_q
      tmax=min(datasb%s_tau0,tmax)
      tmin=max(datasb%s_points(2),tmin)
      if (tmax < datasb%s_points(2)) exit
      sums = 0

      qmax_int= max(850,ThisCT%ls%l(j))*3*BessIntBoost/datasb%s_tau0*1.2
      DoInt = .not. datasb%cp_want_scalars .or. privateindexes%iv_q < qmax_int

      if (DoInt) then
         startloopidx = statbesseindexof (datasb%s_count, datasb%s_R, &
            datasb%s_npoints, datasb%s_Highest, tmin)
         endloopidx = min(privateindexes%iv_sourcessteps,statbesseindexof (datasb%s_count, &
            datasb%s_R, datasb%s_npoints, datasb%s_Highest, tmax))
         do n=startloopidx,endloopidx
            a2=aa(n)
            bes_ix=bes_index(n)

            J_l=a2*ajlin(bes_ix,j)+(1-a2)*(ajlin(bes_ix+1,j) - ((a2+1) &
               *ajlprin(bes_ix,j)+(2-a2)*ajlprin(bes_ix+1,j))* fac(n)) !cubic spline
            J_l = J_l*datasb%s_dpoints(n)

            sums(1) = sums(1) + IVSource_q(n,1)*J_l
            sums(2) = sums(2) + IVSource_q(n,2)*J_l
            sums(3) = sums(3) + IVSource_q(n,3)*J_l
         end do
      end if
      if (.not. DoInt .or. UseLimberGPU(ThisCT%ls%l(j), datasb) &
         .and. datasb%cp_want_scalars) then
         xf = datasb%s_tau0-(ThisCT%ls%l(j)+0.5_dl)/privateindexes%iv_q
         if (xf < datasb%s_highest .and. xf > datasb%s_lowest) then
            n=statbesseindexof (datasb%s_count, datasb%s_R, &
               datasb%s_npoints, datasb%s_Highest, xf)
            xf= (xf-datasb%s_points(n))/(datasb%s_points(n+1)-datasb%s_points(n))
            sums(3) = (IVSource_q(n,3)*(1-xf) + xf*IVSource_q(n+1,3))*&
               sqrt(const_pi/2/(ThisCT%ls%l(j)+0.5_dl))/privateindexes%iv_q
         else
            sums(3)=0
         end if
      end if
      
      ThisCT%Delta_p_l_k(:,j,privateindexes%iv_q_ix) = ThisCT%Delta_p_l_k(:,j,privateindexes%iv_q_ix) + sums
   end do
#endif

end subroutine DoFlatIntegration

! END OPENACC
