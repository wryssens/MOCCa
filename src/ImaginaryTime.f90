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
    use WaveFunctions
    use Force, only : B14, B15, B16, B17
    use Cranking

    type(Spwf)  , intent(in) :: Psi
    type(Spinor)             :: hPsi, U, B, S, A, W, C, D, Crank
    integer                  :: i

    B = ActionOfB(Psi)
    U = ActionOfU(Psi)
    W = ActionOfW(Psi)

    hPsi = NewSpinor()
    hPsi = B + U + W

    if(.not.TRC) then
        S = ActionOfS(Psi)
        A = ActionOfA(Psi)

        hPsi = hPsi + S
        hPsi = hPsi + A

        if(B14.ne.0.0_dp .or. B15 .ne. 0.0_dp) then
            C = ActionOfC(Psi)
            hPsi = hPsi + C
        endif

        if(B16.ne.0.0_dp .or. B17 .ne. 0.0_dp) then
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

      !Calculate the update
      do it=1,2
        Update(:,:,:,it) = c0*(Value(it)  - Desired(it))/(O2(it) + d0)*       &
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

    real(KIND=dp), save             :: Alpha, OldAlpha, f
    type(Spwf), save,allocatable    :: NesterovVectors(:)
    real(KIND=dp)                   :: propfactor,SpEnergy, SpDispersion, q
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
      y(i)= NesterovVectors(i)%GetValue()
      hy(i)  = hPsi(NesterovVectors(i))
      ! Do the Nesterov update
      call HFBasis(i)%SetGrid( y(i)  - propfactor * hy(i) )
    enddo

    !Orthonormalise the wavefunctions
    call Gramschmidt

    !Calculate the new nesterovalpha
    if(Restart.ne.0) then
      if (mod(Iteration,Restart) .eq. 0)  then
          Alpha = 0 ; OldAlpha=0
      endif
    endif

    Alpha = (1 + sqrt( 4 *OldAlpha**2 + 1))/2
    print *, alpha
    print *, (alpha + oldalpha - 1)/alpha,  ( 1 - Oldalpha)/alpha

    do i=1,nwt
        ! Calculate new Nesterov vectors
        y(i) = (alpha + oldalpha - 1)/alpha*HFBasis(i)%GetValue() +  ( 1 - Oldalpha)/alpha *x(i)
        call NesterovVectors(i)%SetGrid(y(i))
        call NesterovVectors(i)%CompDer
    enddo

    !Update Ak
    OldAlpha = Alpha
  end subroutine Nesterov

end module Imaginarytime
