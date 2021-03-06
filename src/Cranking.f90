!-------------------------------------------------------------------------------
! Module that treats all the angular momentum cranking in the code.
!
!-------------------------------------------------------------------------------
! Note for anyone that wants to compare MOCCa to CR8
!   -> homegz in CR8 (and homegzb) are the hbar * omega from the 24Mg paper
!      WITH A MINUS SIGN!
!
!-------------------------------------------------------------------------------
module Cranking

  use CompilationInfo
  use Geninfo
  use SpwfStorage

  implicit none

  !-----------------------------------------------------------------------------
  ! Omega:
  !   Value called omega in the constraints. It is the lagrange multiplier.
  !
  !   (1) => Omega_x ; (2) => Omega_y ;  (3) => Omega_z
  !
  ! CrankValues:
  !   Values of J to use in the constraint.
  !
  ! CrankIntensityScale
  !   Scaling factor for automatically determined CrankIntensity (CrankType = 1)
  !
  ! CrankEnergy
  !   Energy associated with every cranking constraint.
  ! 
  ! CrankType
  !   Determines the type of constraint (for every direction separately)
  !   (0) => Linear constraint for fixed omega, read from data. (No cranking 
  !          when omega is 0).
  !   (1) => Cranking with a corrective step and updating of omega. 
  !
  !-----------------------------------------------------------------------------
  real(KIND=dp), public :: Omega(3)      = 0.0_dp, CrankValues(3)= 0.0_dp
  real(KIND=dp), public :: CrankEnergy(3)= 0.0_dp, OmegaSize     = 0.0_dp
  real(KIND=dp), public :: CrankIntensity(3)      = 0.0_dp
  real(KIND=dp), public :: CrankIntensityScale(3) = 1.0_dp
  real(KIND=dp), public :: Jtotal                 = 0.0_dp
  integer      , public :: CrankType(3)  = 0

  !-----------------------------------------------------------------------------
  ! Whether the cranking update is done with reference to 
  ! .false. => TotalAngMom     , calculated from the spwfs directly
  ! .true.  => TotalAngMom_dens, calculated by integrating the densities 
  logical :: crank_smooth = .false.
  !-----------------------------------------------------------------------------
  ! Whether or not to use the cranking info from file
  logical               :: ContinueCrank= .false.
  !-----------------------------------------------------------------------------
  ! Logical signalling other modules if there is a direction where we do other 
  ! steps.
  logical               :: AlternateCrank, RealignOmega = .false.
  !-----------------------------------------------------------------------------

contains

  subroutine ReadCrankingInfo()
    !---------------------------------------------------------------------------
    ! Subroutine for reading user input on the cranking variables.
    !
    !---------------------------------------------------------------------------
    real(KIND=dp) :: OmegaX=0.0_dp, OmegaY=0.0_dp, OmegaZ=0.0_dp
    real(KIND=dp) :: CrankX=0.0_dp, CrankY=0.0_dp, CrankZ=0.0_dp
    real(KIND=dp) :: IntensityX=0.0_dp,IntensityY=0.0_dp,IntensityZ=0.0_dp
    real(KIND=dp) :: ScaleIntensityX=1.0_dp,ScaleIntensityY=1.0_dp
    real(KIND=dp) :: ScaleIntensityZ=1.0_dp
    integer       :: CrankTypeX=0, CrankTypeY=0, CrankTypeZ=0

    NameList /Cranking/ OmegaX,OmegaY,OmegaZ,CrankX,CrankY,CrankZ,             &
    &                   ContinueCrank, CrankTypeX,CrankTypeY, CrankTypeZ,      &
    &                   OmegaSize, RealignOmega, Jtotal, IntensityX,           &
    &                   IntensityY, IntensityZ,                                &
    &                   ScaleIntensityX,ScaleIntensityY,ScaleIntensityZ,       &
    &                   crank_smooth

    read(unit=*, NML=Cranking)

    !Checking cranking constraints with symmetries.
    if(TRC.and. (CrankX.ne.0.0_dp .or. CrankY.ne.0.0_dp .or. CrankZ.ne.0.0_dp))&
      & call stp('MOCCa cannot crank when Time Reversal is conserved.')
    if(SC   .and. (CrankX.ne.0.0_dp)) then
      call stp('MOCCa cannot crank Jx when signature is conserved.')
    endif
    if(TSC .and. (CrankY.ne.0.0_dp)) then
      call stp('MOCCa cannot crank Jy when time simplex is conserved.')
    endif
    if(CrankTypeX.lt.0 .or. CrankTypeX .gt. 1) then
        call stp('Unrecognised cranktype in the X direction.')
    endif
    if(CrankTypeY.lt.0 .or. CrankTypeY .gt. 1) then
        call stp('Unrecognised cranktype in the Y direction.')
    endif
    if(CrankTypeZ.lt.0 .or. CrankTypeZ .gt. 1) then
        call stp('Unrecognised cranktype in the Z direction.')
    endif
    if(Jtotal .lt. 0.0_dp) then
        call stp('Jtotal constraint should be positive.')
    endif
    if(Jtotal .ne. 0.0_dp .and. SC) then
        call stp('Jtotal only usable when signature is broken at the moment.')
    endif

    !----------------- Assigning Constants based on Input ----------------------
    CrankValues    = (/ CrankX,CrankY,CrankZ/)
    Omega          = (/ OmegaX,OmegaY,OmegaZ/)
    CrankType      = (/ CrankTypeX, CrankTypeY, CrankTypeZ/)
    CrankIntensity = (/ IntensityX, IntensityY, IntensityZ/)
    CrankIntensityScale = (/ ScaleIntensityX, ScaleIntensityY, ScaleIntensityZ/)

    if(.not. ContinueCrank) then
        OmegaSize = sqrt((OmegaX**2 + OmegaY**2 + OmegaZ**2))
    endif
    
    if(any(CrankType.eq.1)) then
      AlternateCrank = .true.
    else
      AlternateCrank = .false.
    endif
    
    if(Jtotal.ne.0.0) then
        CrankType(1) = 1
        CrankType(2) = 0
        CrankType(3) = 1
    endif
    
    return
  end subroutine ReadCrankingInfo

  subroutine UpdateAM()
    !---------------------------------------------------------------------------
    ! This subroutine calculates the angular momentum in each direction for
    ! every Spwf, and sums it. It is also calculated by integration of the 
    ! current and spin density. 
    !---------------------------------------------------------------------------

    use Densities
    integer :: wave, i,P,S,it

    !Save the values of the previous iteration
    AngMomOld   = TotalAngMom
    OldJ2Total  = J2Total

    TotalAngMom = 0.0_dp ; J2Total = 0.0_dp ; AMIsoblock = 0.0_dp
    JTR = 0.0_dp         ; JTI = 0.0_dp    
    TotalAngMom_dens = 0.0_dp    

    do wave=1,nwt
      call HFBasis(wave)%DiagSpin()
      call HFBasis(wave)%DiagAng()
      ! When canonical basis is allocated we need to recalculate the angular
      ! momentum for that basis too.
      if(allocated(CanBasis)) then
        call  CanBasis(wave)%DiagAng()
      endif
    enddo

    !---------------------------------------------------------------------------
    !Only compute angular momentum when TimeReversal is not conserved
    if(.not. TRC) then
      do wave=1,nwt
        do i=1,3
          TotalAngMom(i) = TotalAngMom(i) + &
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%J(i)
          JTR(i)          = JTR(i)    +     & 
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%JTR(i)
          JTI(i)          = JTI(i)    +     & 
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%JTI(i)
          J2Total(i)     = J2Total(i) +     &
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%J2(i)

          S = (DensityBasis(wave)%GetSignature() + 3)/2
          P = (DensityBasis(wave)%GetParity()    + 3)/2
          it= (DensityBasis(wave)%GetIsospin()   + 3)/2

          AMIsoblock(S,P,it,i) = AMIsoblock(S,P,it,i) +                  &
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%J(i)
        enddo
      enddo
      !-------------------------------------------------------------------------
      ! And now we integrate the current density and spin density.
      do it=1,2
        ! Spin part
        TotalAngMom_dens(1) = TotalAngMom_dens(1) + &
        &                                    0.5 * sum(Density%vecs(:,:,:,1,it)) 
        TotalAngMom_dens(2) = TotalAngMom_dens(2) + & 
        &                                    0.5 * sum(Density%vecs(:,:,:,2,it)) 
        TotalAngMom_dens(3) = TotalAngMom_dens(3) + &
        &                                    0.5 * sum(Density%vecs(:,:,:,3,it)) 

        do i=1, nx*ny*nz
          TotalAngMom_dens(1)=TotalAngMom_dens(1) &
          & - Mesh3D(3,i,1,1) * Density%vecj(i,1,1,2,it) &
          & + Mesh3D(2,i,1,1) * Density%vecj(i,1,1,3,it)
     !    & + Mesh3D(3,i,1,1) * Density%vecj(i,1,1,2,it) &  ! bugfix MB 19/12/15
     !    & - Mesh3D(2,i,1,1) * Density%vecj(i,1,1,3,it)

          TotalAngMom_dens(2)=TotalAngMom_dens(2) &
          & - Mesh3D(1,i,1,1) * Density%vecj(i,1,1,3,it) &
          & + Mesh3D(3,i,1,1) * Density%vecj(i,1,1,1,it)
     !    & + Mesh3D(1,i,1,1) * Density%vecj(i,1,1,3,it) &  ! bugfix MB 19/12/15
     !    & - Mesh3D(3,i,1,1) * Density%vecj(i,1,1,1,it)

          TotalAngMom_dens(3)=TotalAngMom_dens(3) &
          & - Mesh3D(2,i,1,1) * Density%vecj(i,1,1,1,it) &
          & + Mesh3D(1,i,1,1) * Density%vecj(i,1,1,2,it)
        enddo
      enddo
      TotalAngMom_dens = TotalAngMom_dens * dv

      !-------------------------------------------------------------------------
      ! Artificially set total J_x and J_y to zero when dictated by symmetries.
      ! In that case the above summation is not \sum v^2_i <Psi_i|J_x|\Psi_i>
      ! but rather \sum v^2_i <Psi_i|J_x T |\Psi_i>
      if(SC)          TotalAngMom(1) = 0.0
      if(SC .or. TSC) TotalAngMom(2) = 0.0
      ! And the integrations of the spin and current densities are wrong because
      ! of the symmetries
      if(SC)          TotalAngMom_dens(1) = 0.0
      if(SC .or. TSC) TotalAngMom_dens(2) = 0.0
    endif

    do i=1,3
      CrankEnergy(i) =  - Omega(i) * TotalAngMom(i)
    enddo

  end subroutine UpdateAM

  function CrankAPot() result(CrankContribution)
    !---------------------------------------------------------------------------
    ! Function returning the contribution of a cranking constraint to the
    ! mean field potential A.
    !       APot => APot - hbar * f_cut(r) \vec{\omega} x \vec{r}
    !---------------------------------------------------------------------------
    use Moments, only: Cutoff
    use Mesh, only   : Mesh3D

    real(KIND=dp) :: CrankContribution(nx,ny,nz,3,2)
    integer       :: i,j,k,it
    
    CrankContribution = 0.0_dp

    ! MB: I added here the trivial if condition for non-zero LeviCivita, which
    ! reduces the cost of executing this to 6/27 of the original cost, modulo
    ! gains compiler optimisation might find in this. But as the function 
    ! LeviCivita() has to be called each time, I doubt this can be efficiently
    ! done anyway.
    do it=1,2
      do k=1,3
      do j=1,3
      do i=1,3
  !     if((i.ne.j).and.(j.ne.k).and.(k.ne.i)) &
            CrankContribution(:,:,:,i,it) = CrankContribution(:,:,:,i,it) -     &
          &     LeviCivita(i,j,k) * Omega(j) *  Mesh3D(k,:,:,:)*Cutoff(:,:,:,it)
     enddo
     enddo
     enddo
   enddo
   return
  end function CrankAPot

  function CrankDivAPot() result(CrankContribution)
    !---------------------------------------------------------------------------
    ! Function returning the contribution of a cranking constraint to the
    ! mean field potential A. Using that div.(fF) = grad(f).F + f div.F, one has
    !       DivAPot => DivAPot - hbar * Grad(f_cut(r)).\vec{\omega} x \vec{r}
    !---------------------------------------------------------------------------
    use Derivatives
    use Moments, only: Cutoff
    use Mesh, only   : Mesh3D
    use Densities

    real(KIND=dp) :: CrankContribution(nx,ny,nz,2)
    real(KIND=dp) :: GradCutoff(nx,ny,nz,3)
    integer       :: i,j,k,it
    
    CrankContribution = 0.0_dp

    do it=1,2
      GradCutoff(:,:,:,1) = &
        & DeriveX(Cutoff(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
      GradCutoff(:,:,:,2) = &
        & DeriveY(Cutoff(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
      GradCutoff(:,:,:,3) = &
        & DeriveZ(Cutoff(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
   
      CrankContribution(:,:,:,it) = CrankContribution(:,:,:,it)     &
          & - GradCutoff(:,:,:,1) * Omega(2) * Mesh3D(3,:,:,:)      &
          & + GradCutoff(:,:,:,1) * Omega(3) * Mesh3D(2,:,:,:)      &
          & - GradCutoff(:,:,:,2) * Omega(3) * Mesh3D(1,:,:,:)      &
          & + GradCutoff(:,:,:,2) * Omega(1) * Mesh3D(3,:,:,:)      &
          & - GradCutoff(:,:,:,3) * Omega(1) * Mesh3D(2,:,:,:)      &
          & + GradCutoff(:,:,:,3) * Omega(2) * Mesh3D(1,:,:,:)
    enddo

 !  print '(  "   k   z    ")'
 !  i = 1 ; j = 1
 !  do k=1,nz
 !    print '(i4,f7.3,8es16.8)',                                             &
 !     & k,Mesh3D(3,i,j,k),                                                  &
 !     & Density%Rho(i,j,k,1),Density%Rho(i,j,k,2),                          &
 !     & Cutoff(i,j,k,1),                                                    &
 !     & GradCutoff(i,j,k,1),GradCutoff(i,j,k,2),GradCutoff(i,j,k,3),        &
 !     & CrankContribution(i,j,k,1),CrankContribution(i,j,k,2) 
 !  enddo
 !  print '(" ")'

  end function CrankDivAPot

  function CrankSPot() result(CrankContribution)
   !----------------------------------------------------------------------------
   ! Function that returns the contribution of a cranking constraint to the S
   ! mean field potential.
   !       SPot => SPot - 1/2 * hbar * \vec{omega}
   !----------------------------------------------------------------------------
   use Moments, only: Cutoff

   real(KIND=dp) :: CrankContribution(nx,ny,nz,3,2)
   integer       :: i, it

   CrankContribution = 0.0_dp
   do it=1,2
     do i=1,3
        CrankContribution(:,:,:,i,it) =                                         &
        & CrankContribution(:,:,:,i,it) -1/2.0_dp*Omega(i)*Cutoff(:,:,:,it)
     enddo
   enddo
   return
  end function CrankSPot

  subroutine ReadjustCranking()
    !---------------------------------------------------------------------------
    ! Readjust the cranking constraints
    !---------------------------------------------------------------------------
    integer                 :: i,j
    real(KIND=dp),parameter :: d0 = 1.0_dp
    real(KIND=dp)           :: SizeJ, value

    OmegaSize = sqrt(sum(Omega**2))

    do i=1,3
        select case(CrankType(i))
        case(0)
          !---------------------------------------------------------------------
          ! Even for constant omega we readjust its direction if realignomega
          ! is active.
          if(realignomega) then
            Omega(i) = TotalAngMom(i) * OmegaSize  &
            & /sqrt(TotalAngMom(1)**2 + TotalAngMom(2)**2 + TotalAngMom(3)**2)
          endif
        case(1)
          !---------------------------------------------------------------------
          ! Judge the intensity of the cranking constraints if not specified by
          ! the user. CrankIntensityScale (default value = 1) rescales this 
          ! estimate if specified by the user.
          if(CrankIntensity(i) .eq. 0.0_dp) then
            CrankIntensity(i) =  CrankIntensityScale(i) / J2Total(i)
          endif
          ! Actual readjustment of the constraint           
          if(crank_smooth) then
            value = TotalAngMom_dens(i)
          else
            value = TotalAngMom(i)
          endif
          Omega(i) = Omega(i)- CrankIntensity(i)*(value-CrankValues(i))

          ! Debugging printout
!          print '(" ReadjustCranking ",i2,7f12.5,l3)',            &
!           & i, Omega(i), Jtotal, J2Total(i), CrankIntensity(i),  &
!           & TotalAngMom(i),CrankValues(i),                       &
!           & CrankIntensity(i)*(TotalAngMom(i) - CrankValues(i))
        end select
    enddo
  end subroutine ReadjustCranking

  subroutine PrintCranking
    !---------------------------------------------------------------------------
    ! Prints all kinds of information about the expectation value of the 
    ! angular momentum operator and all kinds of angles.
    !
    !---------------------------------------------------------------------------
    use Densities

 !  ! format restoration MB 19/12/15
 !  1 format (18('-'), ' Angular Momentum (hbar) ',17('-') )
 !  2 format (15x, 'Spwfs(*)', 4x, 'Densit.   ', 2x, 'Desired', 5x, 'Omega')
 ! 21 format (15x, 'Spwfs   ', 4x, 'Densit.(*)', 2x, 'Desired', 5x, 'Omega')
 !  6 format (2x,' _______________________________________________________' )
 !  3 format (3x,'J_',a1,'   ','|', 4f12.5 )
 ! 31 format (3x,'Size  |', 4f12.5)
 ! 32 format (1x,'ReJT',a1,'   ','|', 4f12.5 )
 ! 33 format (1x,'ImJT',a1,'   ','|', 4f12.5 )
 ! 34 format (2x,'|J|',a1,'   ','|', 4f12.5 )

    1 format (22('-'), ' Angular Momentum (hbar) ',23('-') )
    2 format (15x, 'Spwfs(*)  ',2x, 'Desired', 5x, 'Omega', 7x, 'Energy' 6x,'Densit. ')
   21 format (15x, 'Densit.(*)',2x, 'Desired', 5x, 'Omega', 7x, 'Energy' 6x,'Spwfs   ')
  211 format (15x, '(*) = used for the cranking constraints')
    3 format (3x,'J_',a1,'   ','|', 5f12.5 )
   31 format (3x,'Size  |', 3f12.5,12x,1f12.5)
   32 format (1x,'ReJT',a1,'   ','|', 5f12.5 )
   33 format (1x,'ImJT',a1,'   ','|', 5f12.5 )
   34 format (2x,'|J|' ,a1,'   ','|', 5f12.5 )
   
    4 format (3x,'Theta |', 3f12.5)
   41 format (3x,'Phi   |', 3f12.5)
    
    5 format (3x,'P R_z ','|',6x,'++',10x,'-+',10x,'+-',10x,'--')
    6 format (2x,' _______________________________________________________' )
   66 format (2x,' ___________________________________________________________________' )
    7 format (3x,a1,'_',a1,3x,'|',4f12.5)
   71 format (3x,'The_', a1' |'      , 4f12.5)
   72 format (3x,'Phi_', a1' |'      , 4f12.5)
   
    8 format (3x,'Open spin')
   10 format (3x,a1,1x,'|',3x,'|',4f12.5)
    
   11 format (3x, 'Attention: Omega gets realigned to J.') 
    
    integer :: i
    real(KIND=dp) :: theta(4), phi(4), J

    !---------------------------------------------------------------------------
    ! Print all info on the angular momentum
    print 1
    print *
    if(RealignOmega) then
        print 11 
    endif
    if(crank_smooth) then
      print 21
    else
      print 2
    endif
    print 66
    if(.not.SC) then
      ! format restoration MB 19/12/15
 !    print 3, 'x',TotalAngMom(1), TotalAngMom_dens(1), CrankValues(1), Omega(1) 
      if(crank_smooth) then
        print 3, 'x',TotalAngMom_dens(1), CrankValues(1), Omega(1), CrankEnergy(1), TotalAngMom     (1)
      else
        print 3, 'x',TotalAngMom     (1), CrankValues(1), Omega(1), CrankEnergy(1), TotalAngMom_dens(1)
      endif
    endif
    if(.not.TSC) then
      ! format restoration MB 19/12/15
 !    print 3, 'y',TotalAngMom(2), TotalAngMom_dens(2), CrankValues(2), Omega(2)
      if(crank_smooth) then
        print 3, 'y',TotalAngMom_dens(2), CrankValues(2), Omega(2), CrankEnergy(2), TotalAngMom     (2)
      else
        print 3, 'y',TotalAngMom     (2), CrankValues(2), Omega(2), CrankEnergy(2), TotalAngMom_dens(2)
      endif
    endif
    ! format restoration MB 19/12/15
 !  print 3, 'z',TotalAngMom(3), TotalAngMom_dens(3), CrankValues(3), Omega(3)
    if(crank_smooth) then
      print 3, 'z',TotalAngMom_dens(3), CrankValues(3), Omega(3), CrankEnergy(3), TotalAngMom     (3)
    else
      print 3, 'z',TotalAngMom     (3), CrankValues(3), Omega(3), CrankEnergy(3), TotalAngMom_dens(3)
    endif
    print 66
    ! format restoration MB 19/12/15
 !  print 31, sqrt(sum(totalangmom(1:3)**2)),sqrt(sum(totalangmom_dens(1:3)**2))&
 !  &       , sqrt(sum(crankvalues(1:3)**2)), sqrt(sum(omega(1:3)**2))
    if(crank_smooth) then
      print 31, sqrt(sum(totalangmom_dens(1:3)**2)) , sqrt(sum(crankvalues(1:3)**2)), &
    &           sqrt(sum(omega(1:3)**2))            , sqrt(sum(totalangmom(1:3)**2))
    else
      print 31, sqrt(sum(totalangmom(1:3)**2))      , sqrt(sum(crankvalues(1:3)**2)), &
    &           sqrt(sum(omega(1:3)**2))            , sqrt(sum(totalangmom_dens(1:3)**2))
    endif
    print *
    print 211
    print 6
    if(.not. SC) then
        if(TSC) then
            print 4, atan2(TotalAngMom(1), TotalAngMom(3)) * 180.0/pi, &
            &        atan2(crankvalues(1), crankvalues(3)) * 180.0/pi, &
            &        atan2(Omega(1), Omega(3)) * 180.0/pi
        else
            J = sqrt(sum(TotalAngMom**2))  
            
            print 4, acos(TotalAngMom(3)/J) * 180.0/pi, &
            &        acos(crankvalues(3)/J) * 180.0/pi, &
            &        acos(Omega(3)/sqrt(sum(Omega**2))) * 180.0/pi
        endif
    endif
    if(.not. SC .and. .not. TSC) then
        print 41, atan2(TotalAngMom(2), TotalAngMom(1)) * 180.0/pi, &
        &         atan2(crankvalues(2), crankvalues(1)) * 180.0/pi, &
        &         atan2(Omega(2), Omega(1)) * 180.0/pi
    endif

    print 5
    print 6
    print 32, 'x',JTR(1), 0.0, 0.0, 0.0
    print 32, 'y',JTR(2), 0.0, 0.0, 0.0
    print 32, 'z',JTR(3), 0.0, 0.0, 0.0
    print 6
    print 33, 'x',JTI(1), 0.0, 0.0, 0.0
    print 33, 'y',JTI(2), 0.0, 0.0, 0.0
    print 33, 'z',JTI(3), 0.0, 0.0, 0.0
    print 6
    print 31, sqrt(sum(JTR(:)**2 + JTI(:)**2)) ,0.0,0.0
    print 6
!    print 33, 'x',TotalAngMom(1) + JT(1), 0.0, 0.0, 0.0
!    print 33, 'y',TotalAngMom(2) + JT(2), 0.0, 0.0, 0.0
!    print 33, 'z',TotalAngMom(3) + JT(3), 0.0, 0.0, 0.0
!    print 6
!    print 31, sqrt(sum((TotalAngMom + JT)**2)), 0.0,0.0
    print 6
    print 34, 'x',sqrt(TotalAngMom(1)**2 + JTR(1)**2 + JTI(1)**2), 0.0, 0.0, 0.0
    print 34, 'y',sqrt(TotalAngMom(2)**2 + JTR(2)**2 + JTI(2)**2), 0.0, 0.0, 0.0
    print 34, 'z',sqrt(TotalAngMom(3)**2 + JTR(3)**2 + JTI(3)**2), 0.0, 0.0, 0.0
    print 6
    print 31, sqrt(sum((TotalAngMom**2 + JTR**2 + JTI**2))), 0.0,0.0
    print 6
    print *
    print 5
    print 6
    !---------------------------------------------------------------------------
    ! Print the individual contributions to the angular momentum vector of all
    ! isospin-parity-signature blocks
    if(.not. SC) then
      print 7, 'N','x', AMIsoblock(:,:,1,1)
      print 7, 'P','x', AMIsoblock(:,:,2,1)
      print *
    endif
    if((.not. TSC).and.(.not.SC)) then
      print 7, 'N','y', AMIsoblock(:,:,1,2)
      print 7, 'P','y', AMIsoblock(:,:,2,2)
      print *
    endif
    print 7, 'N','z', AMIsoblock(:,:,1,3)
    print 7, 'P','z', AMIsoblock(:,:,2,3)
    
    if(.not. SC) then
      if(TSC) then
          print *
          theta(1) = atan2(AMIsoblock(1,1,1,1), AMIsoblock(1,1,1,3))
          theta(2) = atan2(AMIsoblock(2,1,1,1), AMIsoblock(2,1,1,3))
          theta(3) = atan2(AMIsoblock(1,2,1,1), AMIsoblock(1,2,1,3))
          theta(4) = atan2(AMIsoblock(2,2,1,1), AMIsoblock(2,2,1,3))
          
          print 71, 'N', theta * 180.0/pi
          
          theta(1) = atan2(AMIsoblock(1,1,2,1), AMIsoblock(1,1,2,3))
          theta(2) = atan2(AMIsoblock(2,1,2,1), AMIsoblock(2,1,2,3))
          theta(3) = atan2(AMIsoblock(1,2,2,1), AMIsoblock(1,2,2,3))
          theta(4) = atan2(AMIsoblock(2,2,2,1), AMIsoblock(2,2,2,3))
          print 71, 'P', theta * 180.0/pi
      else
          print *
          theta(1) = acos(AMIsoblock(1,1,1,3)/sqrt(sum(AMIsoblock(1,1,1,:)**2)))
          theta(2) = acos(AMIsoblock(2,1,1,3)/sqrt(sum(AMIsoblock(2,1,1,:)**2)))
          theta(3) = acos(AMIsoblock(1,2,1,3)/sqrt(sum(AMIsoblock(1,2,1,:)**2)))
          theta(4) = acos(AMIsoblock(2,2,1,3)/sqrt(sum(AMIsoblock(2,2,1,:)**2)))
          print 71, 'N', theta * 180.0/pi
          theta(1) = acos(AMIsoblock(1,1,2,3)/sqrt(sum(AMIsoblock(1,1,2,:)**2)))
          theta(2) = acos(AMIsoblock(2,1,2,3)/sqrt(sum(AMIsoblock(2,1,2,:)**2)))
          theta(3) = acos(AMIsoblock(1,2,2,3)/sqrt(sum(AMIsoblock(1,2,2,:)**2)))
          theta(4) = acos(AMIsoblock(2,2,2,3)/sqrt(sum(AMIsoblock(2,2,2,:)**2)))
          print 71, 'P', theta * 180.0/pi
      endif
    endif
    if( .not. TSC .and. .not. SC) then
        print *
        phi(1) = atan2(AMIsoblock(1,1,1,2), AMIsoblock(1,1,1,1))
        phi(2) = atan2(AMIsoblock(2,1,1,2), AMIsoblock(2,1,1,1))
        phi(3) = atan2(AMIsoblock(1,2,1,2), AMIsoblock(1,2,1,1))
        phi(4) = atan2(AMIsoblock(2,2,1,2), AMIsoblock(2,2,1,1))
        print 72, 'N', phi * 180.0/pi
        
        phi(1) = atan2(AMIsoblock(1,1,2,2), AMIsoblock(1,1,2,1))
        phi(2) = atan2(AMIsoblock(2,1,2,2), AMIsoblock(2,1,2,1))
        phi(3) = atan2(AMIsoblock(1,2,2,2), AMIsoblock(1,2,2,1))
        phi(4) = atan2(AMIsoblock(2,2,2,2), AMIsoblock(2,2,2,1))
        print 72, 'P', phi * 180.0/pi
    endif
    print *
    print 8
    print 6
    if(.not. SC) then
        print 10, 'X', 0.5*sum(Density%vecs(:,:,:,1,1))*dv, &
        &              0.5*sum(Density%vecs(:,:,:,1,2))*dv
    endif
    if(.not. SC .and. .not. TSC) then
        print 10, 'Y', 0.5*sum(Density%vecs(:,:,:,2,1))*dv, &
        &              0.5*sum(Density%vecs(:,:,:,2,1))*dv
    endif
    print 10    , 'Z', 0.5*sum(Density%vecs(:,:,:,3,1))*dv, &
        &              0.5*sum(Density%vecs(:,:,:,3,2))*dv
    print 6
    !---------------------------------------------------------------------------
  end subroutine PrintCranking

  pure logical function ConverCranking(Prec) result(Converged)
  !-----------------------------------------------------------------------------
  ! Subroutine that checks the Cranking constraints for convergence.
  !-----------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Prec
    integer                   :: i

    Converged=.true.
    do i=1,3
      if(CrankType(i).gt.0) then
        if(abs(TotalAngMom(i) - CrankValues(i)).gt.Prec) then
          Converged=.false.
          return
        endif
      endif
    enddo

  end function ConverCranking

end module Cranking
