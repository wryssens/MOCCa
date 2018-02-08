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

  !-----------------------------------------------------------------------------
  ! Spwf to estimate the maximum eigenvalue of h.
  type(Spwf), allocatable :: maxspwf
  !-----------------------------------------------------------------------------
  ! Array that stores the previous update.
  type(Spinor), allocatable :: updates(:)
  !-----------------------------------------------------------------------------
  ! Arrays containing the parameters governing the iterative process. 
  real(KIND=dp), allocatable ::  mom_estimate(:), dt_estimate(:)
  
contains

  subroutine IterativeEstimation(Iteration)
    !---------------------------------------------------------------------------
    ! Estimate optimum parameters 
    !   a) dt_estimate
    !   b) mom_estimate
    ! at every iteration to accelerate convergence. 
    !---------------------------------------------------------------------------
    
    use Pairing
    use HFB
    use force
    
    1 format (a20, 99f10.3)
   
    integer, intent(in) :: iteration
    
    integer             :: iter, estiter, i,j,it,P,S, ind(2,2), N,ii,jj
    real(KIND=dp)       :: timestep, mu, maxE, E, con
    real(KIND=dp)       :: kappa, minE, relE(2,2), compare, ground(2,2)
    type(Spinor)        :: actionofh, update
    
    !---------------------------------------------------------------------------
    ! Initialization
    if(.not.allocated(mom_estimate)) then
        allocate(mom_estimate(nwt)) ; mom_estimate = 0 
        allocate(dt_estimate(nwt)) ;  dt_estimate  = 0
    endif
    
    if(ParameterEstimation .eq. 0) then
        !-----------------------------------------------------------------------
        ! Take user-supplied values when not asked to do estimation.
        dt_estimate  = dt
        mom_estimate = momentum
        return
    endif

    if(Iteration.eq.1) then 
        !-----------------------------------------------------------------------
        ! Use starting parameters if this is the first iteration.
        dt_estimate  = dt
        mom_estimate = momentum
    endif
    
    !---------------------------------------------------------------------------
    ! Step 1: evolve the maxspwf in order to estimate the largest eigenvalue 
    !         on the mesh.
    estiter = 2
    if(Iteration .eq.1) then
        !-----------------------------------------------------------------------
        ! Initialize randomly at the start
        maxspwf = copywavefunction(HFBasis(1))
        call random_number(maxspwf%value%grid)
        call maxspwf%compnorm()
        maxspwf%value = 1.0/sqrt(maxspwf%norm) * maxspwf%value
        
        call maxspwf%compder()  
        
        if(t1n2.ne.0.0_dp .or. t2n2 .ne.0.0_dp) then
            call maxspwf%compsecondder()
        endif
    endif
    estiter = 100
    update=NewSpinor()
    update%grid=0.0

    do iter=1,estiter
        !-----------------------------------------------------------------------
        actionofh = hpsi(maxspwf)
        con       = maxE
        maxE      = InproductSpinorReal(maxspwf%value, actionofh)
        con       = con - maxE
        !-----------------------------------------------------------------------
        ! notice the sign, we are maximising instead of minimising.
        update        = dt_estimate(1)/hbar*(actionofh-maxE*maxspwf%value + mom_estimate(1)*update) 
        maxspwf%value = maxspwf%value + update
        !-----------------------------------------------------------------------
        ! Normalize
        call maxspwf%compnorm()
        maxspwf%value = 1.0/sqrt(maxspwf%norm) * maxspwf%value
        call maxspwf%compder()
        if(t1n2.ne.0.0_dp .or. t2n2 .ne.0.0_dp) then
            call maxspwf%compsecondder()
        endif

        !-----------------------------------------------------------------------
        ! Don't be to picky about convergence, within the order of an MeV is
        ! good enough.
        if(abs(con).lt. 1d-2) exit
    enddo
    !if(abs(con).gt. 1d-2) then
    !    call stp('unreliable upper bound on MaxE')
    !endif
    !---------------------------------------------------------------------------
    ! Step two: find the next excitation value of the next many-body state.
    select case(PairingType)
    
    case(0)
        !-----------------------------------------------------------------------
        ! The excitation energy is here simply the energy difference of the first 
        ! unoccupied state with the last occupied state.
        !-----------------------------------------------------------------------
        relE = 10000
        do it=1,2
            do P=1,2
                 do i=1,nwt
                    if(((HFBasis(i)%parity  + 3)/2) .ne. P ) cycle
                    if(((HFBasis(i)%isospin + 3)/2) .ne. it) cycle
                    if(abs(HFBasis(i)%occupation)   .lt.0.5) cycle
                    
                    do j=1,nwt
                        if(((HFBasis(j)%parity  + 3)/2) .ne. P ) cycle
                        if(((HFBasis(j)%isospin + 3)/2) .ne. it) cycle
                        if(abs(HFBasis(j)%occupation)   .gt.0.5) cycle
 
                        compare = - HFBasis(i)%energy  &
                        &         + HFBasis(j)%energy
        
                        if(compare .gt. 0.1) then
                          relE(P,it) = min(relE(P,it), compare)
                        endif
                     enddo
                 enddo
            enddo
        enddo
        relE = 2*relE 
!        relE = -10000 
!        do it=1,2
!            do P=1,2
!                do i=1,nwt
!                    if((HFBasis(i)%isospin +3)/2 .ne. it) cycle
!                    if((HFBasis(i)%parity  +3)/2 .ne. P ) cycle
!                    if( HFBasis(i)%occupation    .eq. 0)  cycle
!                    
!                    if(HFBasis(i)%energy .gt. relE(P,it)) then
!                        relE(P,it) = HFBasis(i)%energy
!                        ind(P,it)  = i
!                    endif
!                enddo
!            enddo
!        enddo
!        !-----------------------------------------------------------------------
!        ! Step three, estimate the eigenvalue right above
!        relE = 10000
!        do it=1,2
!          do P=1,2
!            do i=1,nwt
!                if((HFBasis(i)%isospin +3)/2 .ne. it) cycle
!                if((HFBasis(i)%parity  +3)/2 .ne. P ) cycle
!                if( HFBasis(i)%occupation    .eq. 1)  cycle
!                if(HFBasis(i)%energy - HFBasis(ind(P,it))%energy .le. 0) cycle                 

!                if(HFBasis(i)%energy-HFBasis(ind(P,it))%energy.lt.relE(P,it))then            
!                    relE(P,it) = HFBasis(i)%energy - HFBasis(ind(P,it))%energy
!                endif
!            enddo
!          enddo
!        enddo
!        !-----------------------------------------------------------------------
!        relE = relE*2
        !-----------------------------------------------------------------------
    case(1)
        !-----------------------------------------------------------------------
        ! When doing BCS, the excitation is the lowest available qp energy.
        relE = 10000
        do it=1,2
            do P=1,2
                 do i=1,nwt
                    if(((HFBasis(i)%parity  + 3)/2) .ne. P ) cycle
                    if(((HFBasis(i)%isospin + 3)/2) .ne. it) cycle
                    
                    do j=1,nwt
                        if(((HFBasis(j)%parity  + 3)/2) .ne. P ) cycle
                        if(((HFBasis(j)%isospin + 3)/2) .ne. it) cycle
                        
                        compare = HFBasis(i)%eqp  &
                        &       + HFBasis(j)%eqp

                        relE(P,it) = min(relE(P,it), compare)
                    enddo
                 enddo
            enddo
        enddo
        !-----------------------------------------------------------------------
    case(2)
        !-----------------------------------------------------------------------
        ! When doing HFB, there are some nasty pitfalls that can happen.
        ! 
        ! The configuration at this particular iteration need not be the lowest
        ! in this particular symmetry block.
        !   a) we are doing cr8-like blocking for the non-lowest qp
        !   b) the HFB solver has not yet found the lowest qp
        !   c) probably other things can happen on the way to convergence
        !
        ! For this reason, the minimum two-qp excitation within the symmetry
        ! block need not be positive! In that case we are not in the ground 
        ! state and need to figure out which one is 
        ! a) the ground state
        ! b) the first excited state 
        !    (since our configuration can be even higher)
        !-----------------------------------------------------------------------
        
        relE = 100000
        do it=1,Iindex
            do P=1,Pindex
                 N = blocksizes(P,it)
                 !--------------------------------------------------------------
                 ! Step one: find the minimum sum of two qps
                 do i=1,N
                    ii = HFBcolumns(i,P,it)
                    do j=1,N
                        jj = HFBcolumns(j,P,it)
                        
                        compare = QuasiEnergies(ii,P,it)  &
                        &       + QuasiEnergies(jj,P,it)
                        relE(P,it) = min(relE(P,it), compare)
                    enddo
                 enddo
                 !--------------------------------------------------------------
                 ! Step two: if the minimum energy is negative, we are not in 
                 ! the ground state currently.
                 !--------------------------------------------------------------
                 if(relE(P,it).lt.0) then
                    ! This is the new reference
                    ground(P,it) = relE(P,it)
                    
                    relE(P,it) = 10000                    
                    ! Find the next state
                    do i=1,N
                        ii = HFBcolumns(i,P,it)
                        do j=1,N
                            jj = HFBcolumns(j,P,it)

                            compare = QuasiEnergies(ii,P,it)  &
                            &       + QuasiEnergies(jj,P,it)
                            if(compare.gt. ground(P,it)) then
                                relE(P,it) = min(relE(P,it), compare)
                            endif
                        enddo
                    enddo
                    relE(P,it) = relE(P,it) - ground(P,it)
                 endif
            enddo
        enddo
    end select

    !---------------------------------------------------------------------------
    ! Have a failsafe for accidental degeneracies
    !where(relE .lt. 0.1) relE = 1
    maxE = maxE - HFBasis(1)%energy
    relE = relE/2
    !---------------------------------------------------------------------------
    ! Step four, estimate dt and the momentum. 
    kappa        = minval(relE)/maxE
    mom_estimate = ((sqrt(kappa) - 1)/(sqrt(kappa) + 1))**2

    dt_estimate = 4.0/(maxE+minval(relE)+2*sqrt(maxE*minval(relE)))*hbar*0.80

    !---------------------------------------------------------------------------
    ! Temporary printing.
    !print *, '----------------'
    !print *, ' MAXE:' , maxE
    !print *, ' dt   ' , dt_estimate(1)
    !print *, ' mom  ' , mom_estimate(1) 
    !print *, ' ind  ' , ind
    !print *, ' relE ' , relE 

  !-----------------------------------------------------------------------------
  end subroutine IterativeEstimation

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
  
  subroutine AlternateStep(AllConstraints)
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

   real(KIND=dp) :: O2(2), Targ(2), Des(2), Value(2), CrankFactor(3), sizeJ
   real(KIND=dp) :: multipole(nx,ny,nz,2), update(nx,ny,nz,2), total, fac
   integer       :: it, i,j, power
   type(Moment),pointer  :: Current
   type(Spinor)  :: QPsi, tempspinor
   logical, intent(in):: AllConstraints

   Current    => Root
   multipole = 0.0_dp
   call compcutoff()
   
   do while(associated(Current%Next))
    Current => Current%next
   
    if(.not. AllConstraints) then
        if(Current%ConstraintType.lt.3) cycle
    else
        if(Current%ConstraintType.eq.0) cycle
    endif
    
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
        if(.not. Current%total) then
            do it=1,2
                Update(:,:,:,it) = 0.5*(Value(it)  - Des(it))/(O2(it) + d0)*   &
                &                      Cutoff(:,:,:,it)*Current%SpherHarm
            enddo
        else
            total = sum(CalculateTotalQl(current%l))
            O2    = sum(CalculateTotalsquared(current%l))
            fac   = sqrt(current%l/(16*pi/5))
            if(Current%m .ne. 0) fac = 2 * fac      
            do it=1,2
                Update(:,:,:,it)=fac*Value(it)*(1 - Des(it)/total)/(O2(it))* &
                &                      Cutoff(:,:,:,it)*Current%SpherHarm
            enddo
        endif
        
        multipole = multipole + Update
    case DEFAULT
        call stp('Alternating constraints do not recognize this type of constraint.')
    end select
   enddo
   !---------------------------------------------------------------------------
   ! Then construct the correction for the cranking constraints.
   CrankFactor=0.0_dp
   if(AlternateCrank) then
     if(Jtotal.eq.0.0_dp) then        
         do i=1,3
             if(.not. AllConstraints) then
                if(CrankType(i).ne.3) cycle
             else
                if(CrankType(i).eq.0) cycle
             endif
             CrankFactor(i) = 0.5*(TotalAngMom(i)-CrankValues(i))/(J2Total(i)+d0)
         enddo
     else
        ! Constraints on total J
        ! Currently only working for J_x + J_z
        sizeJ = sqrt(sum(TotalAngMom(1:3)**2))
        CrankFactor(1) = 4*(1 - JTotal/sizeJ)/(J2Total(1)+d0) * TotalAngMom(1)
        CrankFactor(2) = 0
        CrankFactor(3) = 4*(1 - JTotal/sizeJ)/(J2Total(3)+d0) * TotalAngMom(3)
     endif
   endif
   
   !---------------------------------------------------------------------------
   ! With the update in hand, we update the spwfs
   do i=1,nwt
      QPsi = HFBasis(i)%GetValue()
      it   = (HFBasis(i)%GetIsospin() + 3)/2
      
      TempSpinor = NewSpinor()
      !-------------------------------------------------------------------------
      !Calculating the actions of the cranking correction
      do j=1,3
        if(CrankFactor(j).eq.0.0_dp) cycle
        TempSpinor = TempSpinor + CrankFactor(j)*AngMomOperator(HFBasis(i), j)
      enddo
      
      !Substituting the correction
      QPsi = QPsi - multipole(:,:,:,it)*QPsi - TempSpinor
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
    real(KIND=dp)             :: SpEnergy, SpDispersion, Propfactor, m
    type(Spinor)              :: Current, ActionOfH, ActionOfH2
    integer, intent(in)       :: iteration
    type(Spwf)                :: TempWf
    
    !---------------------------------------------------------------------------
    ! Initialization.
    if(.not. allocated(updates)) then
        allocate(updates(nwt))
        do i=1,nwt
            updates(i) = NewSpinor()
        enddo
    endif
    
    !---------------------------------------------------------------------------
    ! Try to find an appropriate set of parameters (dt, mu), possibly for every
    ! wavefunction. 
    call IterativeEstimation(iteration)
    
    do i=1,nwt
      Current      = HFBasis(i)%GetValue()
      ActionofH    = hPsi(HFbasis(i))
      SpEnergy     = InproductSpinorReal(Current, ActionOfH)
      SpDispersion = InproductSpinorReal(ActionOfH,ActionOfH) - SpEnergy**2
      !-------------------------------------------------------------------------
      !Precondition the descent direction if InverseKineticDamping is active.
      if(InverseKineticDamping) then
        ActionOfH = InverseKinetic(ActionOfH,E0,                        &
        &           HFBasis(i)%GetParity(), HFBasis(i)%GetSignature(),  &
        &           HFBasis(i)%GetTimeSimplex(), HFBasis(i)%GetIsospin())
      endif
      !-------------------------------------------------------------------------
      ! Add a heavy-ball term and save the update for next time
      ActionOfH = -(dt_estimate(i)/hbar)*ActionofH + mom_estimate(i)*updates(i)
      ! Save current value to calculate difference
      updates(i)= current
      !-------------------------------------------------------------------------
      ! Update and save      
      ActionOfH = Current + ActionOfH
    
      call HFBasis(i)%SetGrid(ActionOfH)
      call HFBasis(i)%SetEnergy(SpEnergy)
      call HFBasis(i)%SetDispersion(SpDispersion)
    enddo
    !---------------------------------------------------------------------------
    !Orthonormalisation
    call Gramschmidt
    !---------------------------------------------------------------------------
    ! Calculate the difference for next time (after orthonormalisation!)
    do i=1,nwt
        updates(i) = HFBasis(i)%Value - updates(i)
    enddo
  end subroutine GradDesc
    
  subroutine Nesterov(Iteration)
    !---------------------------------------------------------------------------
    ! Nesterov optimal gradient method.
    !
    !---------------------------------------------------------------------------
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
