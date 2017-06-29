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
  ! Crankdamp:
  !   Damping of the cranking potential.
  ! Crankdamp:
  !   Damping of the cranking potential.
  ! CrankEnergy
  !   Energy associated with every cranking constraint.
  ! CrankType
  !   Determines the type of constraint
  !-----------------------------------------------------------------------------
  real(KIND=dp), public :: Omega(3)      = 0.0_dp, CrankValues(3)= 0.0_dp
  real(KIND=dp), public :: CrankReadj    = 1.0_dp, CrankDamp     = 0.95_dp
  real(KIND=dp), public :: CrankEnergy(3)= 0.0_dp, OmegaSize     = 0.0_dp
  integer      , public :: CrankType(3)  = 0
  !-----------------------------------------------------------------------------
  ! Whether or not to use the cranking info from file
  logical               :: ContinueCrank= .false.
  !-----------------------------------------------------------------------------
  real(KIND=dp)         :: CrankC0=0.8_dp
  !-----------------------------------------------------------------------------
  ! Logical signalling other modules whether or not to do Rutz correction steps
  logical               :: RutzCrank
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
    integer       :: CrankTypeX=0, CrankTypeY=0, CrankTypeZ=0

    NameList /Cranking/ OmegaX,OmegaY,OmegaZ,CrankX,CrankY,CrankZ,             &
    &                   CrankDamp,CrankReadj, ContinueCrank,                   &
    &                   CrankTypeX,CrankTypeY, CrankTypeZ, CrankC0,            &
    &                   OmegaSize, RealignOmega

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
    if(CrankTypeX.lt.0 .or. CrankTypeX .gt. 2) then
        call stp('Unrecognised cranktype in the X direction.')
    endif
    if(CrankTypeY.lt.0 .or. CrankTypeY .gt. 2) then
        call stp('Unrecognised cranktype in the Y direction.')
    endif
    if(CrankTypeZ.lt.0 .or. CrankTypeZ .gt. 2) then
        call stp('Unrecognised cranktype in the Z direction.')
    endif

    !----------------- Assigning Constants based on Input ----------------------
    CrankValues = (/ CrankX,CrankY,CrankZ/)
    Omega       = (/ OmegaX,OmegaY,OmegaZ/)
    CrankType   = (/ CrankTypeX, CrankTypeY, CrankTypeZ/)

    OmegaSize = sqrt((OmegaX**2 + OmegaY**2 + OmegaZ**2))
    
    if(any(CrankType.eq.1) .and. CrankC0.ne.0.0_dp) then
      RutzCrank = .true.
    else
      RutzCrank = .false.
    endif

    return
  end subroutine ReadCrankingInfo

  function CrankAPot() result(CrankContribution)
    !---------------------------------------------------------------------------
    ! Function returning the contribution of a cranking constraint to the
    ! mean field potential A.
    !       APot => APot - hbar* \vec{\omega} x \vec{r}
    !---------------------------------------------------------------------------
    use Moments, only: Cutoff
    use Mesh, only   : Mesh3D

    real(KIND=dp) :: CrankContribution(nx,ny,nz,3,2)
    integer       :: i,j,k,it

    CrankContribution=0.0_dp

    do it=1,2
     do k=1,3
      do j=1,3
        do i=1,3
           CrankContribution(:,:,:,i,it) = CrankContribution(:,:,:,i,it) -     &
           & LeviCivita(i,j,k) * Omega(j) *  Mesh3D(k,:,:,:)*Cutoff(:,:,:,it)
        enddo
      enddo
     enddo
   enddo
   return
  end function CrankAPot

  function CrankSPot() result(CrankContribution)
   !----------------------------------------------------------------------------
   ! Function that returns the contribution ReadjustCrankingof a cranking constraint to the S
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
          Omega(i) = Omega(i) - 2*CrankC0*(TotalAngMom(i) - CrankValues(i))
        end select
    enddo

  end subroutine ReadjustCranking

  subroutine PrintCranking
    
    use Densities

    1 format (18('-'), ' Angular Momentum (hbar) ',17('-') )
    2 format (15x, 'Total', 4x, 'Desired', 4x, 'Omega', 4x, 'Energy')
    3 format (3x,'J_',a1,'   ','|', 4f10.3 )
   31 format (3x,'Size  |', 3f10.3)
   
    4 format (3x,'Theta |', 3f10.3)
   
    5 format (3x,'P R_z ','|',6x,'++',8x,'-+',8x,'+-',8x,'--')
    6 format (2x,' ___________________________________________________' )
    7 format (3x,a1,'_',a1,3x,'|',4f10.3)

    8 format (3x,'Open spin')
    9 format (2x,' ___________________________________________________' )
   10 format (3x,a1,1x,'|',3x,'|',4f10.3)
    
   11 format (3x, 'Attention: Omega gets realigned to J.') 
    
    integer :: i

    do i=1,3
      CrankEnergy(i) =  - Omega(i) * TotalAngMom(i)
    enddo

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
    print 31, sum(totalangmom(1:3)**2) , sum(crankvalues(1:3)**2),             &
    &         sum(omega(1:3)**2)
    if(.not. SC) then
        print 4, atan2(TotalAngMom(1), TotalAngMom(3)), atan2(crankvalues(1), crankvalues(3)) &
        &        atan2(Omega(1), Omega(3))
    endif

    print *
    print *
    print 5
    if(.not. SC) then
      print 6
      print 7, 'N','x', AMIsoblock(:,:,1,1)
      print 7, 'P','x', AMIsoblock(:,:,2,1)
    endif
    if(.not. TSC) then
      print 6
      print 7, 'N','y', AMIsoblock(:,:,1,2)
      print 7, 'P','y', AMIsoblock(:,:,2,2)
    endif

    print 6
    print 7, 'N','z', AMIsoblock(:,:,1,3)
    print 7, 'P','z', AMIsoblock(:,:,2,3)

    print *
    print 8
    print 9
    if(.not. SC) then
        print 10, 'X', 0.5*sum(Density%vecs(:,:,:,1,1))*dv, 0.5*sum(Density%vecs(:,:,:,1,2))*dv
    endif
    if(.not. SC .and. .not. TSC) then
        print 10, 'Y', 0.5*sum(Density%vecs(:,:,:,2,1))*dv, 0.5*sum(Density%vecs(:,:,:,2,1))*dv
    endif
    print 10    , 'Z', 0.5*sum(Density%vecs(:,:,:,3,1))*dv, 0.5*sum(Density%vecs(:,:,:,3,2))*dv
    print *

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
