! START OPENACC

#define  IVSQROWS 3600
#define  IVSQCOLS 3

function statbesseindexof (count, R, npoints, Highest, tau)
#ifdef USEACC
!$acc routine 
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

subroutine SourceToTransfers(datasb, &
   ThisCT, q_ix,  ThisSourcesin, ScaledSrcin, &
   ddScaledSrcin, max_etak_tensorin, max_etak_vectorin, &
   WantLateTimein, max_etak_scalarin, full_bessel_integrationin, &
   do_bispectrumin, max_bessels_l_indexin, &
   xlimfracin, xlimminin, ajlin, ajlprin, DebugEvolutionin)
   !IVSource_q)
#ifdef USEACC
!acc routine vector 
!$acc routine 
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
   ! local data to avoid a call to IntegrationVars_Init
   integer i,khi,klo, step
   real(dl) xf,b0,ho,a0,ho2o6,a03,b03
   integer :: ixunit
   integer :: local_step
   !local data to avoid a call to DoSourceIntegration
   integer :: j,ll,llmax
   real(dl) nu
   real(dl) :: sixpibynu

   IVSource_q(1,:)=0
   IVSource_q(datasb%s_npoints,:) = 0
   IVSource_q(datasb%s_npoints-1,:) = 0

   privateindexes%iv_q_ix = q_ix
   privateindexes%iv_q = ThisCT%q%points(q_ix)
   privateindexes%iv_dq = ThisCT%q%dpoints(q_ix)

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

   local_step = 2 
   step = 2
   ! if GNU
   !acc parallel loop vector reduction(max:local_step)
   !$acc loop vector reduction(max:local_step)
   do i=2, datasb%s_npoints
     xf=privateindexes%iv_q*(datasb%s_tau0-datasb%s_points(i))
     IVSource_q(i,:) = 0.0_dl ! Initialize to 0 first
     if (datasb%cp_want_tensors) then
       if (privateindexes%iv_q*datasb%s_points(i) < max_etak_tensorin .and. &
         xf > 1.e-8_dl) then
          IVSource_q(i,:) = a0*ScaledSrcin(klo,:,i)+&
                           b0*ScaledSrcin(khi,:,i)+(a03 *ddScaledSrcin(klo,:,i)+ &
                           b03*ddScaledSrcin(khi,:,i)) * ho2o6
          local_step = i ! Update local_step if condition met
        end if
      end if
      if (datasb%cp_want_scalars) then
        if ((DebugEvolutionin .or. WantLateTimein .or. &
            privateindexes%iv_q*datasb%s_points(i) < max_etak_scalarin) &
            .and. xf > 1.e-8_dl) then
          IVSource_q(i,:) = a0 * ScaledSrcin(klo,:,i) +  &
                           b0 * ScaledSrcin(khi,:,i) + (a03*ddScaledSrcin(klo,:,i) + &
                           b03 * ddScaledSrcin(khi,:,i)) * ho2o6
          local_step = i ! Update local_step if condition met
        end if
      end if
   end do
   !$acc end loop
   ! if GNU
   !acc end parallel loop 

   step = local_step ! Assign the final max value to step
   privateindexes%iv_sourcessteps = step

   nu=privateindexes%iv_q*datasb%s_curvature_radius
   sixpibynu  = 6._dl*3.1415926535897932384626433832795_dl/nu

   llmax=nint(nu*datasb%s_chi0)
   if (llmax<15) then
      llmax=17 !AL Sept2010 changed from 15 to get l=16 smooth
   else
      llmax = nint(nu*(datasb%s_tau0/datasb%s_curvature_radius + sixpibynu))
   end if

   call DoFlatIntegration(ThisCT, llmax, ThisSourcesin, &
      full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
      datasb,xlimfracin,xlimminin,ajlin,ajlprin,privateindexes, IVSource_q)

end subroutine SourceToTransfers

subroutine DoFlatIntegration(ThisCT, llmax, ThisSourcesin, &
   full_bessel_integrationin, do_bispectrumin, max_bessels_l_indexin, &
   datasb, xlimfracin, xlimminin, ajlin, ajlprin, privateindexes, &
   IVSource_q)
#ifdef USEACC
!$acc routine 
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
   real(dl) xf
   real(dl) qmax_int
   integer bes_ix,n, bes_index(privateindexes%iv_sourcessteps)
   integer custom_source_off, s_ix
   integer nwin
   real(dl) :: BessIntBoost
   real(dl) :: temp_sum1, temp_sum2, temp_sum3
   integer :: startloopidx, endloopidx

   BessIntBoost = datasb%cp_accuracy_boost*datasb%cp_accuracy_bessintboost
   custom_source_off = datasb%s_num_redshiftwindows + datasb%s_num_extra_redshiftwindows + 4

   !$acc loop vector
   do j=1,privateindexes%iv_sourcessteps !Precompute arrays for this k
      xf=abs(privateindexes%iv_q*(datasb%s_tau0-datasb%s_points(j)))
      bes_index(j)=statbesseindexof (datasb%b_count, &
         datasb%b_R, datasb%b_npoints, datasb%b_Highest, xf)
      bes_ix= bes_index(j)

      fac(j)=datasb%b_points(bes_ix+1)-datasb%b_points(bes_ix)
      aa(j)=(datasb%b_points(bes_ix+1)-xf)/fac(j)
      fac(j)=fac(j)**2*aa(j)/6
   end do

   !$acc loop vector
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
       ! if uysiing GNU
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
       
     end if
   
     ! This section updates sums(3) based on different logic.
     ! It's a direct assignment, not an accumulation within a parallel loop,
     ! so it should be fine with respect to the *reported* error.
     if (.not. DoInt .or.&
          (ThisCT%ls%l(j) > 400 * (datasb%cp_accuracy_boost * &
            datasb%cp_accuracy_liber_boost)** 0.5) & ! was UseLimberGPU assuming datasb%cp_st_limber_windows F
         .and. datasb%cp_want_scalars) then
       xf = datasb%s_tau0-(ThisCT%ls%l(j)+0.5_dl)/privateindexes%iv_q
       if (xf < datasb%s_highest .and. xf > datasb%s_lowest) then
         n=statbesseindexof (datasb%s_count, datasb%s_R, &
             datasb%s_npoints, datasb%s_Highest, xf)
         xf= (xf-datasb%s_points(n))/(datasb%s_points(n+1)-datasb%s_points(n))
         temp_sum3 = (IVSource_q(n,3)*(1-xf) + xf*IVSource_q(n+1,3))*&
             sqrt(const_pi/2/(ThisCT%ls%l(j)+0.5_dl))/privateindexes%iv_q
       else
         temp_sum3 = 0.0_dl
       end if
     end if
   
     ThisCT%Delta_p_l_k(1, j, privateindexes%iv_q_ix) = &
       ThisCT%Delta_p_l_k(1, j, privateindexes%iv_q_ix) + temp_sum1
     ThisCT%Delta_p_l_k(2, j, privateindexes%iv_q_ix) = &
       ThisCT%Delta_p_l_k(2, j, privateindexes%iv_q_ix) + temp_sum2
     ThisCT%Delta_p_l_k(3, j, privateindexes%iv_q_ix) = & 
       ThisCT%Delta_p_l_k(3, j, privateindexes%iv_q_ix) + temp_sum3
     !ThisCT%Delta_p_l_k(:,j,privateindexes%iv_q_ix) = ThisCT%Delta_p_l_k(:,j,privateindexes%iv_q_ix) + sums
   end do

end subroutine DoFlatIntegration

! END OPENACC
