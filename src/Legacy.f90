module legacy

  implicit none 

contains

subroutine Evolve_OLD(MaxIterations, iprint)
  !-----------------------------------------------------------------------------
  ! This subroutine commands the time evolution of the system.
  !
  ! The integer iprint is a flag: iprint=1 gives extra diagnostic printouts.
  !
  ! I (W.R.) supply this old routine WITH ABSOLUTELY NO GUARANTEE!
  !
  !-----------------------------------------------------------------------------
  use CompilationInfo
  use GenInfo
  use WaveFunctions
  use Derivatives
  use PotentialMixing
  use ImaginaryTime
  use SpwfStorage
  use Densities
  use Moments    
  use Pairing
  use Cranking
  use Cranking
  use DensityMixing
  use Energy 
  use HFB
  use Testing

  implicit none

  integer :: dummyHF(2,2,2,2) = 0

  !---------------------------------------------------------------------------
  ! Interface for the Summaryprinting routines, since FORTRAN is rather
  ! peculiar about routines in the same program...
  interface
    subroutine PrintSummary_v1(Iteration)
      integer, intent(in) :: iteration
    end subroutine PrintSummary_v1
    subroutine PrintSummary_v2(Iteration)
      integer, intent(in) :: iteration
    end subroutine PrintSummary_v2
  end interface

  procedure(PrintSummary_v1), pointer :: PrintSummary


  1 format ( 20x, "Iteration ", i4, 20x,/)
  2 format ( "Total Time Elapsed:", f8.3 )
  3 format (A11 , f8.3 , " or ", f8.3 , " %")
  4 format (60("-"))
  5 format (22('-'),"Congratulations!",22('-'),/,                            &
  &                " Your calculation converged, to the following accuracy: ")
  6 format ("Multipole moments changed for less than : ", es9.2)
  7 format ("Energy            changed for less than : ", es9.2)
  8 format ("Total dispersion of the occupied Spwfs  : ", es9.2)
  9 format ("Average dispersion of the occupied Spwfs: ", es9.2)
100 format (/,60('='),/, 18x,' START OF THE ITERATIVE PROCESS ' ,/, 60('='))
101 format (/,60('='),/, 18x,' **FINAL** Iteration ')
102 format (/,60('='),/, 18x,' PROJECTION ON FEASIBLE SUBSPACE ' ,/, 55('='))

  logical, intent(in) :: iprint
  integer, intent(in) :: MaxIterations
  integer             :: Iteration
  logical             :: Convergence
  logical, external   :: ConvergenceCheck
  integer             :: i,wave
  logical             :: AlternateCheck
  real(KIND=dp)       :: Canenergy

  !--------------------------------------------------------------------------
  ! Decide which summary to print
  PrintSummary => PrintSummary_v2
  
  ! Assign correct evolution operator for the iterations
  select case(IterType)
  case('IMTS')
    EvolveSpwf => GradDesc
  case('NEST')
    EvolveSpwf => Nesterov
  case DEFAULT
    call stp('Itertype is not recognised.')
  end select

  !Updating certain parameters of the wavefunctions that were not read from file
  do i=1,nwt
    call HFBasis(i)%SymmetryOperators()
  enddo
  
  if(SolvePairingStart) then
    if(FreezeOccupation .and. PairingType.eq.0) then
        call HFFill(DummyHF) ! Fill in the HF way anyway for the first iteration
    else
        call SolvePairing
    endif
  else
  ! Temporary, until MOCCa saves canonical basis to file.
    if(PairingType.eq.2) then
      call stp('Need to calculate Canonical basis at the start')
    endif
  endif
  if(t1n2.ne. 0.0_dp .or. t2n2.ne.0.0_dp)   then
    if(MAXFDORDER .ne. -1 .or. MAXFDLAPORDER .NE. -1) then
        call stp('It is irresponsible to calculate N2LO without lag derivatives.')
    endif
    call N2LODerive()
  else
    call DeriveAll()  
  endif
  !Update the Angular momentum variables
  call UpdateAm()
  !Make sure that there is no improper readjusting of cranking constraints
  AngMomOld = TotalAngMom

  ! Recalculate densities if needed
  ! This includes several scenarios:
  ! 1) If the user asked for it
  ! 2) if the input file was a CR8 file, the time-odd densities could then
  !    not be read from file.
  ! 3) if the input file was an EV8 file, some derivatives of densities could not
  !    be read.
  if(Recalc) then
    call UpdateDensities(0)
  endif
   
  ! Printing the parameters of this run.
  call PrintStartInfo
  !Loop over the iterations
  Iteration   = 0
  Convergence = .false.
  ! We calculate all of the multipole moments three times. This is strictly 
  ! wasteful, but only done at the zeroth iteration. This is to make sure
  ! the correct values are used in all readjustment formulas.
  call CalculateAllMoments(1)
  call CalculateAllMoments(1)
  call CalculateAllMoments(1)

  !-----------------------------------------------------------------------------
  !Construct the mean-field potentials
  call ConstructPotentials

  !Calculating  The Energy
  call CompEnergy()

  if(Pairingtype.eq.2) then
    do i=1,nwt
      CanEnergy = InproductSpinorReal(CanBasis(i)%GetValue(),hPsi(Canbasis(i)))
      call Canbasis(i)%SetEnergy(CanEnergy)
    enddo
  endif

  !Printing observables
  if(maxiterations .eq.0) print 101
  call PrintIterationInfo(0, .true.)

  !Checking whether we should do alternate constraint steps.
  AlternateCheck = CheckForAlternateConstraints() 
 
  if(MaxIterations.ne.0) then
        print 100
  endif
  !-----------------------------------------------------------------------------
  !Start of the Mean-field iterations
  iteration = 0
  do while((Iteration.lt.MaxIterations))
    !Incrementing the iteration number
    Iteration = Iteration + 1
    !---------------------------------------------------------------------------
    ! Calculate all the different potentials to construct the single particle
    ! hamiltonian h and then substitute every spwf |\Psi> by 1 - dt/hbar*h|\Psi>
    ! or more complicated algorithms in general.
    ! Note that the orthogonalisation is managed by this routine too.
    call EvolveSpwf(Iteration)
    
    !---------------------------------------------------------------------------
    !Pairing.
    call SolvePairing

    !---------------------------------------------------------------------------
    AlternateCheck = CheckForAlternateConstraints()
    
    !---------------------------------------------------------------------------
    ! When moments or angular moments are constrained according to K.Rutz'
    ! prescription, there needs to be a correction step, dependent on the
    ! density or angular moment obtained in  the meantime.
    if(AlternateCheck) then
      ! Save the value used to construct the mean-fields to DensityHistory
      call UpdateDensities(0,.true.) 
      call CalculateAllMoments(1) ! Save old values to history
      call ReadjustAllMoments()  ! Only readjust Rutz-style moments
      call ReadjustAllMoments()  ! Only readjust alternating
    endif
    if( AlternateCrank) then
        !Deriving all Spwf
        if(t1n2.ne. 0.0_dp .or. t2n2.ne.0.0_dp) then
            call N2LODerive()
        else
            call DeriveAll()    
        endif
        call updateAm()
        call ReadjustCranking()
    endif
    
    !---------------------------------------------------------------------------
!    Checking for the presence of alternating constraints
    if(AlternateCheck .or. AlternateCrank) then
        call AlternateStep()
        if(t1n2.ne. 0.0_dp .or. t2n2.ne.0.0_dp) then
          call N2LODerive()
        else
          call DeriveAll()
        endif
        call SolvePairing
    endif

    !---------------------------------------------------------------------------
!    Deriving all Spwf
    if(t1n2.ne. 0.0_dp .or. t2n2.ne.0.0_dp) then
      call N2LODerive()
    else
      call DeriveAll()
    endif
    !---------------------------------------------------------------------------
    !Calculate the expectation values of the symmetry operators and the spwf
    !energies of the canonical basis
    do i=1,nwt
      call HFBasis(i)%SymmetryOperators()
      if(Pairingtype.eq.2) then
        CanEnergy = InproductSpinorReal(CanBasis(i)%GetValue(), hPsi(Canbasis(i)))
        call Canbasis(i)%SetEnergy(CanEnergy)
        call CanBasis(i)%SymmetryOperators()
      endif
    enddo

    !Update all the angular momentum information. This can be skipped on
    !iterations without printout if there is no cranking.
    !if( any(CrankType.ne.0)  .or.  mod(Iteration, PrintIter).eq.0 ) then
    call UpdateAM()
    !endif

    !Print a summary
    if(mod(Iteration, PrintIter).ne.0) call PrintSummary(Iteration)

    !Computing the densities.
    if(AlternateCheck ) then
        ! If we did an intermediate step, don't use the current values in 
        ! Density as history.
        call UpdateDensities(1)
    else
        ! If we did not do an intermediate step, these are still fine.
        call UpdateDensities(0)
    endif
  
    !Smooth the densities using the difference between this and the previous
    ! iteration
    call MixDensities(Iteration)

    !See if some moments were temporary
    call TurnOffConstraints(iteration)

    !Calculating the expected values of all relevant Moments
    call CalculateAllMoments(1)

    !Readjust the constraints
    call ReadjustAllMoments()
    
    call ReadjustCranking()

    !Construct the mean-field potentials
    call ConstructPotentials

    !Calculating the Energy
    call CompEnergy
    
    !Checking for convergence
    Convergence = ConvergenceCheck()
  	if(Convergence) exit
    
    if(mod(Iteration, PrintIter).eq.0 .or. Iteration.eq.MaxIterations) then
      if((Iteration.eq.MaxIterations) .and. (.not. Maxiterations.eq.0)) print 101
      !Print info after selected amount of iterations.
      call PrintIterationInfo(Iteration, .true.)
    endif
  enddo
  
  !End of the mean-field iterations
  !-----------------------------------------------------------------------------
  if (Convergence .and. mod(Iteration, PrintIter).ne.0) then
    ! If convergence happens during iterations, it is possible
    ! that a printout has been skipped
    call printIterationInfo(Iteration, .true.)
  endif

  !Reanalysis of the result with Lagrange derivatives.
  call FinalIteration()

  if(Convergence) then
     print 5
     print 6, MomentPrec
     print 7, EnergyPrec
     print 8, TotalDispersion
     print 9, TotalDispersion/(protons+neutrons)
     print 4
  endif
  
!  if(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp) call testN2LOsph
  
  !Write densities to files.
  if(Pictures) then
    call PlotDensity()
!    call PlotCurrents(15,1)
  endif
end subroutine Evolve_OLD

end module legacy
