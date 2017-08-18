module ImaginaryTime

  use CompilationInfo
  use Geninfo
  use MeanFields
  use Spwfstorage

  implicit none

  !-----------------------------------------------------------------------------
  !Procedure that determines the evolution of a Spwf under imaginary time.
  abstract interface
      subroutine Evolve_interface(Iteration)
        integer, intent(in)       :: iteration
      end subroutine
  end interface
  procedure(Evolve_Interface),pointer :: EvolveSpwf

contains

  function hPsi(Psi)
    !---------------------------------------------------------------------------
    ! Calculate the action of \hat{h} on a wavefunction.
    !---------------------------------------------------------------------------
    ! Depending on whether time-reversal is conserved or not, and depending on 
    ! the type of functional, different fields are used to construct the 
    ! single-particle hamiltonian.
    !
    ! Note that all of the signs have been absorbed into the definitions of the 
    ! actions, in order that all signs in this routine can be '+'.
    !---------------------------------------------------------------------------
    use WaveFunctions
    use Force 
    use Cranking

    type(Spwf)  , intent(in) :: Psi
    type(Spinor)             :: hPsi, U, B, S, A, W, C, D, Crank, DN2LO, X, T, P
    integer                  :: i

    hPsi = NewSpinor()
    
    !---------------------------------------------------------------------------
    ! Time-even part 
    B = ActionOfB(Psi) ! Field of tau (or tau_mn)
    U = ActionOfU(Psi) ! Field of rho
    W = ActionOfW(Psi) ! Spin-orbit density
   
    hPsi = B + U + W

    if(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp) then
        !--------------------------------------
        ! N2LO time-even fields 
        DN2LO = ActionOfDN2LO(Psi)   ! Field of Q density
        X     = ActionOfX(Psi)       ! Field of Vmn density
        T     = ActionOfImTField(Psi)! Field of Im Tmnk density
        hPsi = hpsi + DN2LO + X + T
    endif

    !---------------------------------------------------------------------------
    ! Time-odd part 
    if(.not.TRC) then
        S = ActionOfS(Psi)  ! Field of (small) s density
        A = ActionOfA(Psi)  ! Field of j density
        hPsi = hPsi + S + A
        
        if(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp) then
            !--------------------------------------
            ! Time-odd N2LO fields
            T  = ActionOfReTField(Psi) ! Field of Re Tmnk density
            P  = ActionOfPi(Psi)       ! Field of Pi density
            S  = ActionOfSN2LO(Psi)    ! Field of (big) S density
            hPsi = hPsi + T + P + S         
        endif

        if(B14.ne.0.0_dp .or. B15 .ne. 0.0_dp) then
            !--------------------------------------
            ! time-odd tensor fields
            C = ActionOfC(Psi)
            hPsi = hPsi + C
        endif

        if(B16.ne.0.0_dp .or. B17 .ne. 0.0_dp) then
            !--------------------------------------
            ! Other time-odd tensor field
            D = ActionOfD(Psi)
            hPsi = hPsi + D
        endif
        
    endif
  end function hPsi

  subroutine RutzCorrectionStep
    !---------------------------------------------------------------------------
    ! Subroutine that takes an extra corrective step to accelerate the
    ! convergence of constrained calculations.
    !---------------------------------------------------------------------------
    use Wavefunctions
    use Moments
    use Spwfstorage
    use Cranking

    real(KIND=dp) :: Correction(nx,ny,nz,2), O2(2), Update(nx,ny,nz,2)
    real(KIND=dp) :: CrankFactor(3), Value(2), Desired(2)
    integer       :: it, i,j, power
    type(Moment),pointer  :: Current, Extra
    type(Spinor)  :: QPsi, TempSpinor

    Current    => Root
    Correction = 0.0_dp
    !---------------------------------------------------------------------------
    ! We first construct the correction from all the moments that are Rutz-type
    ! constrained.
    do while(associated(Current%Next))
      Current => Current%Next
      if(Current%ConstraintType.ne.2) cycle

      !Ordinary constraints
      select case(Current%Isoswitch)
      case(1)
        !---------------------------------------------------------------------
        ! Constrain the sum of the proton & neutron parts
        Desired = Current%Constraint(1)

        power = 1
        O2    = sum(Current%Squared)
        Value = sum(Current%Value)

        if(Current%Total) then
          Value   = sum(CalculateTotalQl(Current%l))
          O2      = sum(Current%Squared)
        endif
        !endif
      case(2)
        !---------------------------------------------------------------------
        !Constrain proton & neutron contributions to different values
        Desired = Current%Constraint
        power = 1
        O2    = Current%Squared
        Value = Current%Value
      case(3)
        !---------------------------------------------------------------------
        ! Constrain the difference between proton and neutron
        ! Note the extra minus signs!
        Desired(1) =  Current%Constraint(1)
        Desired(2) = -Current%Constraint(2)

        power = 1
        O2    = sum(Current%Squared)
        Value(1) =   Current%Value(1) - Current%Value(2)
        Value(2) = - Value(1)
      end select

      !-------------------------------------------------------------------------
      !Calculate the update
      do it=1,2
        Update(:,:,:,it) = c0/2*(Value(it)  - Desired(it))/(O2(it) + d0)*      &
        &                  Current%SpherHarm
      enddo
      Correction = Correction + Update
    enddo
    call CompCutoff
    !---------------------------------------------------------------------------
    ! Then construct the correction for the cranking constraints.
    CrankFactor=0.0_dp
    if(RutzCrank) then
      do i=1,3
          if(CrankType(i).ne.1) cycle
          CrankFactor(i)=CrankC0*                                                &
          &              (TotalAngMom(i)-CrankValues(i))/(J2Total(i)+d0)
      enddo
    endif
    !---------------------------------------------------------------------------
    ! Then we update the spwf-functions
    do i=1,nwt
      QPsi = HFBasis(i)%GetValue()
      it   = (HFBasis(i)%GetIsospin() + 3)/2

      TempSpinor = NewSpinor()
      !Calculating the actions of the cranking correction
      do j=1,3
        if(CrankFactor(j).eq.0.0_dp) cycle
        TempSpinor = TempSpinor + CrankFactor(j)*AngMomOperator(HFBasis(i), j)
      enddo

      !Substituting the correction
      QPsi = QPsi - Correction(:,:,:,it)*Cutoff(:,:,:,it)*QPsi - TempSpinor
      call HFBasis(i)%SetGrid(QPsi)
    enddo
    ! Finally, orthonormalisation
    call Gramschmidt
  end subroutine RutzCorrectionStep
  
  subroutine AlternateStep
  !-----------------------------------------------------------------------------
  ! Subroutine performing one (or more) alternate step for the alternating
  ! constraints. The idea is a a simple gradient step in the direction of a
  ! satisfied constraint, meaning that the objective function being minimized is
  !
  !  c0 ( < O > - O_target )^2
  ! 
  ! and we update the single-particle wavefunctions according to 
  !
  ! psi = ( 1 - epsilon \hat{O} ) psi
  !  
  ! with
  !
  ! epsilon = 1/2 * ( <C> - C )/( < C^2 >)
  !
  ! where < C >^2 is the one-body part of the two-body operator C.
  !-----------------------------------------------------------------------------

   use Wavefunctions
   use Moments
   use Spwfstorage
   use Cranking

   real(KIND=dp) :: O2(2), Targ(2), Des(2), Value(2), CrankFactor(3)
   real(KIND=dp) :: multipole(nx,ny,nz,2), update(nx,ny,nz,2)
   integer       :: it, i,j, power
   type(Moment),pointer  :: Current
   type(Spinor)  :: QPsi, tempspinor

   Current    => Root
   multipole = 0.0_dp
   
   do while(associated(Current%Next))
    Current => Current%next
    
    if(Current%ConstraintType.ne.3) cycle
    select case(Current%Isoswitch)
    case(1)
        !-----------------------------------------------------------------------
        ! Constrain the sum of the proton & neutron parts       
        power = 1
        O2    = sum(Current%Squared)                    ! < C^2 >
        Value = sum(Current%Value)                      ! Current value of <C>
        Targ  = Current%Constraint(1)                   ! Current targeted value
        Des   = Current%TrueConstraint(1)               ! Desired final value
        
        !-----------------------------------------------------------------------
        !Calculate the update
        do it=1,2
            Update(:,:,:,it) = 0.5*(Value(it)  - Des(it))/(O2(it) + d0)*      &
            &                                                  Current%SpherHarm
        enddo
        multipole = multipole + Update
    case DEFAULT
        call stp('Alternating constraints do not recognize this type of constraint.')
    end select
   enddo
   !---------------------------------------------------------------------------
   ! Then construct the correction for the cranking constraints.
   CrankFactor=0.0_dp
   if(AlternateCrank) then
     do i=1,3
         if(CrankType(i).ne.3) cycle
         CrankFactor(i) = 0.5*(TotalAngMom(i)-CrankValues(i))/(J2Total(i)+d0)
     enddo
   endif
   
   !---------------------------------------------------------------------------
   ! With the update in hand, we update the spwfs
   call compcutoff()
    do i=1,nwt
      QPsi = HFBasis(i)%GetValue()
      it   = (HFBasis(i)%GetIsospin() + 3)/2
      
      TempSpinor = NewSpinor()
      !Calculating the actions of the cranking correction
      do j=1,3
        if(CrankFactor(j).eq.0.0_dp) cycle
        TempSpinor = TempSpinor + CrankFactor(j)*AngMomOperator(HFBasis(i), j)
      enddo
      
      !Substituting the correction
      QPsi = QPsi - multipole(:,:,:,it)*Cutoff(:,:,:,it)*QPsi - TempSpinor
      call HFBasis(i)%SetGrid(QPsi)
    enddo
   !---------------------------------------------------------------------------
   ! Finally, orthonormalisation
   call Gramschmidt
  end subroutine AlternateStep

  subroutine GradDesc(Iteration)
    !----------------------------------------------------------------
    ! Subroutine that performs gradient descent on the wavefunctions,
    ! also known as imaginary timestep.
    !----------------------------------------------------------------

    use GenInfo
    use Damping

    integer :: i
    real(KIND=dp)             :: SpEnergy, SpDispersion, Propfactor
    type(Spinor)              :: Current, ActionOfH, ActionOfH2
    integer, intent(in)       :: iteration

    type(Spwf)                :: TempWf
    Propfactor = dt/hbar

    if(TaylorOrder.eq.2) print *, 'Taylor'

    do i=1,nwt
      Current      = HFBasis(i)%GetValue()
      ActionofH    = hPsi(HFbasis(i))
      SpEnergy     = InproductSpinorReal(Current, ActionOfH)
      SpDispersion = InproductSpinorReal(ActionOfH,ActionOfH) - SpEnergy**2

      !Precondition the descent direction if InverseKineticDamping is active.
      if(InverseKineticDamping) then
        ActionOfH = ActionOfH - SpEnergy*Current
        ActionOfH = InverseKinetic(ActionOfH,E0,                        &
        &           HFBasis(i)%GetParity(), HFBasis(i)%GetSignature(),  &
        &           HFBasis(i)%GetTimeSimplex(), HFBasis(i)%GetIsospin())
      endif

      if(TaylorOrder.eq.2) then
          TempWF = NewWaveFunction(ActionOfH,HFBasis(i)%GetIsospin(),          &
          &      HFBasis(i)%GetTimeSimplex(), HFBasis(i)%GetParity(),          &
          &      HFBasis(i)%GetSignature(), HFBasis(i)%GetTimeReversal())
          call TempWF%CompDer()

          ActionOfH2 = hPsi(TempWF)
          if(InverseKineticDamping) then
            ActionOfH2 = ActionOfH2 - SpDispersion * Current
            ActionOfH2 = InverseKinetic(ActionOfH2, E0 ,HFBasis(i)%GetParity(),&
            & HFBasis(i)%GetSignature(), HFBasis(i)%GetTimeSimplex(),          &
            & HFBasis(i)%GetIsospin())
          endif

          ActionOfH = ActionofH - (propfactor * 0.5_dp) * ActionofH2
      endif

      ActionOFH    = - propfactor*ActionOfH
      ActionOfH    =   Current + ActionOfH
      call HFBasis(i)%SetGrid(ActionOfH)
      call HFBasis(i)%SetEnergy(SpEnergy)
      call HFBasis(i)%SetDispersion(SpDispersion)
    enddo

    !Orthonormalisation
    call Gramschmidt
  end subroutine GradDesc

  subroutine Nesterov(Iteration)
    !----------------------------------
    ! Nesterov optimal gradient method.
    !
    !----------------------------------
    use geninfo
    use Damping
    use Energy, only : totalenergy, oldenergy, nesterovsignal
    
    real(KIND=dp), save             :: Alpha, OldAlpha, f
    type(Spwf), save,allocatable    :: NesterovVectors(:)
    real(KIND=dp)                   :: propfactor,SpDispersion, q
    integer, intent(in)             :: Iteration
    integer                         :: i
    integer,allocatable,save        :: NIter(:)
    type(Spinor)                    :: x(nwt), hx(nwt), y(nwt), hy(nwt), temp(nwt)

    propfactor = dt/hbar

    if(Iteration.eq.1) then
      Alpha = 0 ; Oldalpha=0
      allocate(NesterovVectors(nwt)) ; NesterovVectors = HFBasis
    endif


    do i=1,nwt
      x(i)= HFBasis(i)%GetValue()
      !if(NesterovSignal.eq.1) then 
      !  y(i)= x(i)
      !else
        y(i)= NesterovVectors(i)%GetValue()
      !endif
      hy(i)  = hPsi(NesterovVectors(i))
      ! Do the Nesterov update
      call HFBasis(i)%SetGrid( y(i)  - propfactor * hy(i) )
    enddo

    !Orthonormalise the wavefunctions
    call Gramschmidt

    !Calculate the new nesterovalpha
    !if(NesterovSignal.eq.1) then
      !if (mod(Iteration,Restart) .eq. 0)  then
!          Alpha = 1 ; OldAlpha=1
!          print *, 'restarted', TotalEnergy, OldEnergy(1)
      !endif
!      do i=1,nwt
!          ! Calculate new Nesterov vectors
!          y(i) =HFBasis(i)%GetValue()
!          call NesterovVectors(i)%SetGrid(y(i))
!          call NesterovVectors(i)%CompDer
!      enddo      
    !endif
    Alpha = (1 + sqrt( 4 *OldAlpha**2 + 1))/2
    do i=1,nwt
        ! Calculate new Nesterov vectors
        y(i) = ((alpha + oldalpha - 1)/alpha)*HFBasis(i)%GetValue() +  ((1-Oldalpha))/alpha *x(i)
        call NesterovVectors(i)%SetGrid(y(i))
        call NesterovVectors(i)%CompDer
    enddo
    !endif

    !Update Ak
    OldAlpha = Alpha
  end subroutine Nesterov

end module Imaginarytime
