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
  !   Value called omega in the constraints. It is the lagrange multiplicator.
  !
  !   (1) => Omega_x ; (2) => Omega_y ;  (3) => Omega_z
  ! CrankValues:
  !   Values of J to use in the constraint.
  ! Jtotal
  !   Size of J if the constraint is only on the total size
  ! Crankdamp:
  !   Damping of the cranking potential.
  ! CrankIntensityScale
  !   Scaling factor for automatically determined CrankIntensity (CrankType = 3)
  ! CrankEnergy
  !   Energy associated with every cranking constraint.
  ! CrankType
  !   Determines the type of constraint
  !-----------------------------------------------------------------------------
  real(KIND=dp), public :: Omega(3)      = 0.0_dp, CrankValues(3)= 0.0_dp
  real(KIND=dp), public :: CrankReadj    = 1.0_dp, CrankDamp     = 0.95_dp
  real(KIND=dp), public :: CrankEnergy(3)= 0.0_dp, OmegaSize     = 0.0_dp
  real(KIND=dp), public :: CrankIntensity(3)      = 0.0_dp
  real(KIND=dp), public :: CrankIntensityScale(3) = 1.0_dp
  real(KIND=dp), public :: Jtotal                 = 0.0_dp
  integer      , public :: CrankType(3)  = 0
  !-----------------------------------------------------------------------------
  ! Whether or not to use the cranking info from file
  logical               :: ContinueCrank= .false.
  !-----------------------------------------------------------------------------
  real(KIND=dp)         :: CrankC0=0.8_dp
  !-----------------------------------------------------------------------------
  ! Logical signalling other modules whether or not to do Rutz correction steps
  logical               :: RutzCrank              , AlternateCrank
  logical               :: RealignOmega = .false.
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
    &                   CrankDamp,CrankReadj, ContinueCrank,                   &
    &                   CrankTypeX,CrankTypeY, CrankTypeZ, CrankC0,            &
    &                   OmegaSize, RealignOmega, Jtotal, IntensityX,           &
    &                   IntensityY, IntensityZ,                                &
    &                   ScaleIntensityX,ScaleIntensityY,ScaleIntensityZ

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
    if(CrankDamp.gt.1.0_dp .or. CrankDamp .lt.0.0_dp) then
      call stp('Crankdamp should be between 0 and 1.')
    endif
    if(CrankReadj.lt.0.0_dp) then
      call stp('CrankReadj should not be negative!')
    endif
    if(CrankTypeX.lt.0 .or. CrankTypeX .gt. 3) then
        call stp('Unrecognised cranktype in the X direction.')
    endif
    if(CrankTypeY.lt.0 .or. CrankTypeY .gt. 3) then
        call stp('Unrecognised cranktype in the Y direction.')
    endif
    if(CrankTypeZ.lt.0 .or. CrankTypeZ .gt. 3) then
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
    if(any(CrankType.eq.1) .and. CrankC0.ne.0.0_dp) then
      RutzCrank = .true.
    else
      RutzCrank = .false.
    endif
    
    if(any(CrankType.eq.3) .or. (Jtotal .ne. 0.0)) then
      AlternateCrank = .true.
    else
      AlternateCrank = .false.
    endif
    
    if(Jtotal.ne.0.0) then
        CrankType(1) = 3
        CrankType(2) = 0
        CrankType(3) = 3
    endif
    
    return
  end subroutine ReadCrankingInfo

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

  subroutine ReadjustCranking(Rutz)
    !---------------------------------------------------------------------------
    ! Readjust the Lagrange multipliers based on  Rutz' prescription.
    ! Rutz determines whether the Rutz or other constraints should be readjusted
    !---------------------------------------------------------------------------
    integer       :: i,j
    real(KIND=dp),parameter :: d0 = 1.0_dp
    logical, intent(in) :: Rutz
    real(KIND=dp) :: SizeJ

    OmegaSize = sqrt(sum(Omega**2))

    do i=1,3
        select case(CrankType(i))
        case(0)
          if(realignomega) then
            Omega(i) = TotalAngMom(i) * OmegaSize  &
            & /sqrt(TotalAngMom(1)**2 + TotalAngMom(2)**2 + TotalAngMom(3)**2)
          endif
        case(1)
          if(.not.Rutz) cycle
          Omega(i) = Omega(i) -                                                &
          & CrankReadj*(TotalAngMom(i) - AngMomOld(i))/(J2Total(i) + d0)
        case(2)
          if(Rutz) cycle
          Omega(i) = Omega(i) -                                                &                                     
          &     2*CrankIntensity(i)*Crankreadj*(TotalAngMom(i) - CrankValues(i))
        case(3)
     !    print '(" ReadjustCranking ",i2,7f12.5,l3)',   &
     !      & i,Omega(i),Jtotal,J2Total(i),CrankIntensity(i),TotalAngMom(i),CrankValues(i), &
     !      & CrankIntensity(i)*(TotalAngMom(i) - CrankValues(i)),Rutz

          if(.not. Rutz) cycle
          !---------------------------------------------------------------------
          ! Judge the intensity of the cranking constraints.
          ! CrankIntensityScale (default value = 1) rescales the estimate
          if(CrankIntensity(i) .eq. 0.0_dp) then
            CrankIntensity(i) = 2 * CrankIntensityScale(i) / J2Total(i)
          endif
          if(Jtotal.eq.0.0_dp) then
            Omega(i) = Omega(i) - CrankIntensity(i)*(TotalAngMom(i) - CrankValues(i))
          endif
     !    print '(" ReadjustCranking ",i2,1f12.5)',i,Omega(i)
        end select
    enddo
    
    if(Jtotal .ne. 0.0_dp .and. (.not. SC)) then
        SizeJ = sqrt(sum(TotalAngMom(1:3)**2))
        Omega(1) = Omega(1) - 2/(J2Total(1))*(1 - Jtotal/SizeJ)*(TotalAngMom(1))
        Omega(3) = Omega(3) - 2/(J2Total(3))*(1 - Jtotal/SizeJ)*(TotalAngMom(3))
    endif

  end subroutine ReadjustCranking

  subroutine PrintCranking
    !---------------------------------------------------------------------------
    ! Prints all kinds of information about the expectation value of the 
    ! angular momentum operator and all kinds of angles.
    !
    !---------------------------------------------------------------------------
    use Densities

    1 format (18('-'), ' Angular Momentum (hbar) ',17('-') )
    2 format (15x, 'Total', 7x, 'Desired', 5x, 'Omega', 7x, 'Energy')
    3 format (3x,'J_',a1,'   ','|', 4f12.5 )
   31 format (3x,'Size  |', 3f12.5)
   32 format (1x,'ReJT',a1,'   ','|', 4f12.5 )
   33 format (1x,'ImJT',a1,'   ','|', 4f12.5 )
   34 format (2x,'|J|',a1,'   ','|', 4f12.5 )
   
    4 format (3x,'Theta |', 3f12.5)
   41 format (3x,'Phi   |', 3f12.5)
    
    5 format (3x,'P R_z ','|',6x,'++',10x,'-+',10x,'+-',10x,'--')
    6 format (2x,' _______________________________________________________' )
    7 format (3x,a1,'_',a1,3x,'|',4f12.5)
   71 format (3x,'The_', a1' |'      , 4f12.5)
   72 format (3x,'Phi_', a1' |'      , 4f12.5)
   
    8 format (3x,'Open spin')
   10 format (3x,a1,1x,'|',3x,'|',4f12.5)
    
   11 format (3x, 'Attention: Omega gets realigned to J.') 
    
    integer :: i
    real(KIND=dp) :: theta(4), phi(4), J

    do i=1,3
      CrankEnergy(i) =  - Omega(i) * TotalAngMom(i)
    enddo
    !---------------------------------------------------------------------------
    ! Print all info on the angular momentum
    print 1
    print *
    if(RealignOmega) then
        print 11 
    endif
    print 2
    print 6
    if(.not.SC) then
      print 3, 'x',TotalAngMom(1), CrankValues(1), Omega(1), CrankEnergy(1)
    endif
    if(.not.TSC) then
      print 3, 'y',TotalAngMom(2), CrankValues(2), Omega(2), CrankEnergy(2)
    endif
    print 3, 'z',TotalAngMom(3), CrankValues(3), Omega(3), CrankEnergy(3)
    print 6
    print 31, sqrt(sum(totalangmom(1:3)**2)) , sqrt(sum(crankvalues(1:3)**2)), &
    &         sqrt(sum(omega(1:3)**2))
    print *
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
        print 10, 'X', 0.5*sum(Density%vecs(:,:,:,1,1))*dv, 0.5*sum(Density%vecs(:,:,:,1,2))*dv
    endif
    if(.not. SC .and. .not. TSC) then
        print 10, 'Y', 0.5*sum(Density%vecs(:,:,:,2,1))*dv, 0.5*sum(Density%vecs(:,:,:,2,1))*dv
    endif
    print 10    , 'Z', 0.5*sum(Density%vecs(:,:,:,3,1))*dv, 0.5*sum(Density%vecs(:,:,:,3,2))*dv
    print 6
    !---------------------------------------------------------------------------
  end subroutine PrintCranking

!   subroutine PrintCranking
!     !---------------------------------------------------------------------------
!     ! A subroutine that prints info on the cranking constraints.
!     ! Included in the output:
!     !   - Values for angular momentum
!     !   - Values for constraints (both readjusted & true)
!     !   - Energy of the constraints
!     !
!     !Also includes some more info on the constraints if this routine is called
!     !for the first time.
!     !---------------------------------------------------------------------------

!     1 format (18('-'), ' Angular Momentum (hbar) ',17('-') )
!     2 format (11x, 3x,'Jx',6x, 'Jy', 6x,'Jz')
!     3 format (3x,A8,3f8.3)
!     4 format ('Constraints :')
!     5 format ('Energy (MeV):')
!     6 format ('Cons. Type  :')
!     7 format ('  L. or Q.? :', 3(3x, A2, 3x))

!     9 format ('  Readjus.  :', 3x, f8.3)
!    10 format ('KermanOnishi:')
!    11 format ('  Omega x J :', 3x, f8.3)
!    12 format ('  (Di-Dj)*Lk:', 3x, f8.3)
!    13 format ('Angle mu/J  :', 3x, f8.3)
!    14 format (' Rutz C0    :', 3x, f8.3)
!   100 format (60("-"))

!     integer          :: i
!     logical          :: FirstTime=.true.
!     real(KIND=dp)    :: Angle

!     !No Cranking with time-reversal Symmetry
!     if(TRC) return

!     print 1
!     print 2
!     print 3,'        ', TotalAngMom

!     !First check for the presence of constraints. If none are present,
!     ! don't print them.
!     if(all(CrankType.eq.0)) return

!     print 4
!     if(CrankReadj.ne.0.0_dp) then
!       print 3,' Des.:', CrankValues
!     endif

!     print 3,' Omega :', Omega
!     print 5

!     do i=1,3
!       CrankEnergy(i) =  - Omega(i) * TotalAngMom(i)
!     enddo

!     print 3,'        ', CrankEnergy

!     !Printing some extra parameters if this is the first time.
!     if(FirstTime) then
!         print 9, CrankReadj
!         print 14, CrankC0
!         FirstTime=.false.
!     endif
!     !Checking the symmetries, if it is useful to calculate and print
!     ! the Kerman-Onishi conditions
!     !if(.not.SC .or. .not. TSC) then
!     !  call CalcKermanOnishi
!     !
!     !  print 10
!     !  print 11, KermanOnishiRHS
!     !  print 12, KermanOnishilHS
!     !endif

!     !Calculating the angle between desired and actual moment
!     !Angle = 0.0_dp
!     !do i=1,3
!     !  Angle = Angle + TrueCrank(i) * TotalAngMom(i)
!     !enddo
!     !Angle = Angle/(sqrt(sum(TrueCrank(1:3)**2) * sum(TotalAngMom(1:3)**2)))
!     !Angle = acos(Angle)
!     !
!     !print 13 , Angle/(2*pi) * 360.0_dp
!     !
!     !print 100

!     return
!   end subroutine PrintCranking

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
