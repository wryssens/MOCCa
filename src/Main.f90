program MOCCa
    !---------------------------------------------------------------------------
    !            ; ,                                                           |
    !           ) ;( (                                                         |
    !          ( (  ) ;                                                        |
    !           ,-"""-.                                                        |
    !        ,-|`-...-'|                                                       |
    !       ((_|       |                                                       |
    !        `-\       /                                                       |
    !           `.___.'                                                        |
    ! Creation: 8th of April   , 2013.                                         |
    ! Wouter Ryssens                                                           |
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Some general comments and things to keep in mind
    !
    ! *) Output is currently roughly organized as envisioned by Wouter, and as 
    !    such might seem suboptimal to many people. It is however slated for 
    !    redesign as soon as Michael, Wouter and whoever has a strong opinion on
    !    the subject agree on a course to follow. 
    ! *) For now (November 2016) a small redesign is planned to make the output
    !    more grep/awk/python-analysable. 
    !---------------------------------------------------------------------------

    use GenInfo
    use InOutput, only    : Input, Output, PlotDensity
    use SpwfStorage, only : PrintSpwf, DensityBasis, HFBasis
    use Pairing, only     : PairingType
    use Testing
    use Transform
    implicit none

    real*8 :: time(2)
    integer :: i
    
 200 format (/,' ___________________________________________________________', &
     &       /,'|                                                          |', &
     &       /,'|  MOCCa                                                   |', &
     &       /,'|                                          Unstable        |', &
     &       /,'|  Copyright  P.-H. Heenen, M.Bender & W. Ryssens          |', &
     &       /,'|                 ; ,                                      |', &
     &       /,'|                ) ;( (                                    |', &
     &       /,'|               ( (  ) ;                                   |', &
     &       /,'|                ,-"""-.                                   |', &
     &       /,"|             ,-|`-...-'|                                  |", &
     &       /,'|            ((_|       |                                  |', &
     &       /,'|             `-\       /                                  |', &
     &       /,"|                `.___.'                                   |", &
     &       /,'|                                                          |', &
     &       /,'|__________________________________________________________|')

     print 200

     !--------------------------------------------------------------------------
     ! Reading the user input and the wavefunction input.
     call Input()    !Override is disabled
     !--------------------------------------------------------------------------
     !Canbasis needs to be properly initialised in order to start if doing HFB
     DensityBasis => HFBasis

     !--------------------------------------------------------------------------
     ! In testing mode.
     if(TestRun.eq.1) then
        print*, "MOCCa is entering test mode!"
        call stp('End of TestRun')
     endif
     !--------------------------------------------------------------------------
     ! Solve the mean-field equations.
     call Evolve(MaxIter, .true.)
     !--------------------------------------------------------------------------
     ! Write wavefunctions to file!
     call Output
     call stp('MOCCa exits correctly.')

end program MOCCa

subroutine Evolve(MaxIterations, iprint)
  !-----------------------------------------------------------------------------
  ! This subroutine commands the time evolution of the system.
  !
  ! The integer iprint is a flag: iprint=1 gives extra diagnostic printouts.
  !-----------------------------------------------------------------------------
  use CompilationInfo
  use GenInfo
  use WaveFunctions
  use Derivatives,only  : MaxFDOrder
  use MeanFields, only  : ConstructPotentials
  use ImaginaryTime
  use SpwfStorage
  use Densities, only   : UpdateDensities, DampingParam, Density, Recalc
  use Moments, only     : CalculateAllMoments, ReadjustAllMoments,             &
  &                       CheckForRutzMoments, TurnOffConstraints,             &
  &                       CheckForAlternateConstraints, totaldeviation,        &
  &                       ConstraintsProject, projectioniter, projectionprecision
  use Pairing, only     : SolvePairing, PairingType, SolvePairingStart
  use Cranking, only    : ReadjustCranking, CrankC0
  use DensityMixing,only: MixDensities
  use Energy 
  use Cranking, only    : RutzCrank, CrankType, AlternateCrank
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
  6 format ("Multipole moments changed for less than : ", e9.2)
  7 format ("Energy            changed for less than : ", e9.2)
  8 format ("Total dispersion of the occupied Spwfs  : ", e9.2)
  9 format ("Average dispersion of the occupied Spwfs: ", e9.2)
100 format (/,60('='),/, 18x,' START OF THE ITERATIVE PROCESS ' ,/, 60('='))
101 format (/,60('='),/, 18x,' **FINAL** Iteration ')
102 format (/,60('='),/, 18x,' PROJECTION ON FEASILBLE CONSTRAINTS ' ,/, 55('='))

  logical, intent(in) :: iprint
  integer, intent(in) :: MaxIterations
  integer             :: Iteration
  logical             :: Convergence
  logical, external   :: ConvergenceCheck
  integer             :: i,wave
  logical             :: RutzCheck, AlternateCheck
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
  call UpdateAm(.true.)
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
  
  if(Pairingtype.eq.2) then
    do i=1,nwt
      CanEnergy = InproductSpinorReal(CanBasis(i)%GetValue(),hPsi(Canbasis(i)))
      call Canbasis(i)%SetEnergy(CanEnergy)
    enddo
  endif
  
  !Calculating  The Energy
  call CompEnergy()
  
  !Printing observables
  call PrintIterationInfo(0, .true.)

  !-----------------------------------------------------------------------------
  ! If asked for, we project first on the subspace of valid constraints.
  if(ConstraintsProject) then
      do iteration = 1, ProjectionIter 
            !Apply Corrections and resolve pairing
            call AlternateStep(.true.)
            call SolvePairing
            call UpdateDensities(0)
            call updateAM(.false.)
            call CalculateAllMoments(1)
            print *, iteration, totaldeviation
            if(TotalDeviation .lt. ProjectionPrecision) exit
      enddo
  
      print 102
      !-------------------------------------------------------------------------
      ! Derive everything again
      if(t1n2.ne. 0.0_dp .or. t2n2.ne.0.0_dp)   then
        if(MAXFDORDER .ne. -1 .or. MAXFDLAPORDER .NE. -1) then
            call stp('It is irresponsible to calculate N2LO without lag derivatives.')
        endif
        call N2LODerive()
      else
        call DeriveAll()  
      endif
      call compEnergy()
      !Printing observables
      call PrintIterationInfo(0, .true.)
  endif
  !-----------------------------------------------------------------------------
  !Construct the mean-field potentials
  call ConstructPotentials
  
  if(Maxiterations.eq.0) then
        print 101
  endif

  !Checking for the presence of Rutz-Type constraints
  RutzCheck      = CheckForRutzMoments()
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

    ! Calculate all the different potentials to construct the single particle
    ! hamiltonian h and then substitute every spwf |\Psi> by 1 - dt/hbar*h|\Psi>
    ! or more complicated algorithms in general.
    ! Note that the orthogonalisation is managed by this routine too.
    call EvolveSpwf(Iteration)

    !Pairing.
    call SolvePairing
    !---------------------------------------------------------------------------
    !Checking for the presence of Rutz-Type constraints
    RutzCheck      = CheckForRutzMoments()
    AlternateCheck = CheckForAlternateConstraints()
    
    !---------------------------------------------------------------------------
    ! When moments or angular moments are constrained according to K.Rutz'
    ! prescription, there needs to be a correction step, dependent on the
    ! density or angular moment obtained in  the meantime.
    if(RutzCheck .or. AlternateCheck) then
      call UpdateDensities(0,.true.)
      call CalculateAllMoments(1) ! Save old values to history
      call ReadjustAllMoments(2)  ! Only readjust Rutz-style moments
    endif
    if(RutzCrank .or. AlternateCrank) then
        !Deriving all Spwf
        if(t1n2.ne. 0.0_dp .or. t2n2.ne.0.0_dp) then
            call N2LODerive()
        else
            call DeriveAll()
        endif
        call updateAm(.true.)
        call ReadjustCranking(.true.)
    endif
    
    if(RutzCheck .or. RutzCrank) then
        !Apply Corrections and resolve pairing
        call RutzCorrectionStep
        call SolvePairing
    endif
    !---------------------------------------------------------------------------
    !Checking for the presence of alternating constraints
    if(AlternateCheck .or. AlternateCrank) then
        call AlternateStep(.false.)
        call SolvePairing
    endif

    !---------------------------------------------------------------------------
    !Deriving all Spwf
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
    call UpdateAM(.false.)
    !endif

    !Print a summary
    if(mod(Iteration, PrintIter).ne.0) call PrintSummary(Iteration)

    !Computing the densities, without smoothing.
    call UpdateDensities(0)

    !Smooth the densities!
    call MixDensities(Iteration)

    !See if some moments were temporary
    call TurnOffConstraints(iteration)

    !Calculating the expected values of all relevant Moments
    call CalculateAllMoments(1)

    !Readjust the constraints
    call ReadjustAllMoments(1)
    call ReadjustAllMoments(3)
    
    call ReadjustCranking(.false.)

    !Construct the mean-field potentials
    call ConstructPotentials

    !Calculating the Energy
    call CompEnergy

    !Checking for convergence
    Convergence = ConvergenceCheck()
  	if(Convergence) exit

    if(mod(Iteration, PrintIter).eq.0 .or. Iteration.eq.MaxIterations) then
      if(Iteration.eq.MaxIterations) print 101
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

  !if(recalcN2LO) call N2LOAnalysis()

  !call writeN3LODensities(Density)

  if(Convergence) then
     print 5
     print 6, MomentPrec
     print 7, EnergyPrec
     print 8, TotalDispersion
     print 9, TotalDispersion/nwt
     print 4
  endif
  
  if(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp) call testN2LOsph
end subroutine Evolve

logical function ConvergenceCheck() result(Converged)
  !-----------------------------------------------------------------------------
  ! Aggregate function that combines the result from several modules
  ! regarding convergence.
  ! What should be converged:
  !       1) Multipole moments
  !       2) Angular Momentum(in the case of Time Reversal Breaking)
  !          (Not implemented yet!)
  !       3) Energies
  !       4) D2H of the spwf (not implemented yet)
  !       5) Fermi Energy( not implemented yet)
  !-----------------------------------------------------------------------------
  use GenInfo
  use Moments, only : ConverMultipoleAll
  use Energy, only  : ConverEnergy, TotalEnergy, PrintENergy
  use Pairing, only : ConverFermi
  use Cranking, only: ConverCranking

  implicit none

  Converged = .true.
  !Checking convergence of the multipole Degrees of freedom
  Converged = Converged .and. ConverMultipoleAll(MomentPrec)
  ! Checking convergence of the Energy
  Converged = Converged .and. ConverEnergy(EnergyPrec)
  !Checking convergence of the Fermi energy
  Converged = Converged .and. ConverFermi(PairingPrec)
  !Checking convergence of the angular momentum for cranked calculations
  Converged = Converged .and. ConverCranking(CrankPrec)

  if(TotalEnergy.ge.0.0_dp) then
    call PrintEnergy
    call stp('Positive total energy!', 'Total Energy', TotalEnergy)
  endif
  return
end function ConvergenceCheck

subroutine PrintStartInfo()
    !---------------------------------------------------------------------------
    ! A subroutine that prints all the relevant info that is available at
    ! the start of MOCCa.
    !---------------------------------------------------------------------------
    use InOutput, only : PrintInput,PrintFileInfo
    use Coulomb, only  : PrintCoulombInfo
    use Force,   only  : PrintForce
    use Moments, only  : PrintConstraints, PrintMomentInfo
    use Cranking,only  : PrintCranking
    use Pairing, only  : PrintPairingInfo
    use Spwfstorage, only: UpdateAM, DeriveAll
    use geninfo

    implicit none

    1 format(60('-'))

    !Printing compilation and runtime input.
    call PrintInput
    ! Print info on the force that is used.
    call PrintForce
    !Print info on the pairing interaction
    call PrintPairingInfo
    !Print info on the Coulomb Parameters
    call PrintCoulombInfo
    !Print info on the Moment parameters
    call printMomentInfo
    !Print info on the cranking
    call DeriveAll()
    ! Print info from the file
    call PrintFileInfo

    call UpdateAm(.true.)
    if(.not.TRC) then
      call PrintCranking
    endif

end subroutine PrintStartInfo

subroutine PrintIterationInfo(Iteration, PrintAll)
  !-----------------------------------------------------------------------------
  ! Wrapper subroutine for all the printing routines.
  !
  !-----------------------------------------------------------------------------
  use Spwfstorage, only : PrintSpwf
  use Moments,     only : PrintAllMoments
  use Cranking,    only : PrintCranking
  use Pairing,     only : PrintPairing, PairingType, Fermi
  use Energy,      only : PrintEnergy, PairingEnergy
  use HFB, only         : PrintQP
  use geninfo

  implicit none

  integer, intent(in) :: Iteration
  logical, intent(in) :: PrintAll

  1 format (/,60('='))
  2 format (' Iteration ', i5)

  print 1
  print 2, Iteration

  call printSpwf(PairingType, Fermi, PrintAll)
  if(PairingType.eq.2) call PrintQp()
  call PrintAllMoments
  if(.not. TRC) then
        call PrintCranking
  endif
  call PrintPairing
  call PrintEnergy

end subroutine PrintIterationInfo

subroutine FinalIteration()
  !-----------------------------------------------------------------------------
  ! Subroutine that does some final manipulations, though not a real iteration!
  ! 1) Recalculate densities with Lagrange derivatives.
  ! 2) Recalculate Energies with Lagrange derivatives.
  !-----------------------------------------------------------------------------
  use Derivatives
  use Spwfstorage, only : DeriveAll, DensityBasis, nwt
  use Densities,   only : UpdateDensities
  use Pairing,     only : SolvePairing
  use InOutput,    only : Pictures, PlotDensity
  use Energy
    
  implicit none

  integer               :: i

  ! Telling the derivatives module to use Lagrangian derivatives (but not
  ! changing the Coulomb derivations...)
  call AssignFDCoefs(-1, -1, CoulombLapOrder)

  !Deriving wavefunctions, either canonical or HF basis
  do i=1,nwt
    call DensityBasis(i)%CompDer()
  enddo

  !Recalculating densities without damping and/or Mixing
  call UpdateDensities(0)
   
  !Computing Energies
  call CompEnergy
  
  !Printing Energy (with Lagrange mention...)
  call PrintEnergy(.true.)

  ! Recalculating the finite difference derivatives, in order
  ! to write the densities of the final iteration to file
  call AssignFDCoefs(MaxFDOrder, MaxFDLapOrder, CoulombLapOrder)
  do i=1,nwt
    call DensityBasis(i)%CompDer()
  enddo

  !Write densities to files.
  if(Pictures) call PlotDensity()

  call UpdateDensities(0)

end subroutine FinalIteration

subroutine PrintSummary_v2(Iteration)
  !-----------------------------------------------------------------------------
  ! This routine prints a summary of the iteration.
  ! Multilined, as I like to have plenty of info.
  !
  !-----------------------------------------------------------------------------

  use CompilationInfo
  use Energy, only: TotalEnergy, OldEnergy, Routhian, OldRouthian
  use Spwfstorage, only : TotalAngMom
  use Densities,   only : DensityChange
  use Pairing, only : Fermi
  use Cranking
  use Moments
  use GenInfo

  implicit none

  integer, intent(in)  :: Iteration
  type(Moment),pointer :: Current
  real(KIND=dp):: Q20, Q22, dQ20, dQ22, Dispersion, Q30, Q32, dQ30, dQ32, Q10, dQ10
  integer      :: i

  1 format('----------- Iteration ', i6,'--------------')
  2 format("Summ.     E=" f10.3, "   dE=",e10.3)
 21 format("Summ.     R=" f10.3, "   dR=",e10.3)
 22 format("Summ. Fermi=" f10.3, "      ",f10.3) 
  
  3 format("Summ.   Q20=" f10.3, "  Q22=",f10.3)
  4 format("Summ.  dQ20=" e10.3, " dQ22=",e10.3)
 31 format("Summ.   Q10=" f10.3)
 41 format("Summ.  dQ10=" e10.3)
 32 format("Summ.   Q30=" f10.3, "  Q32=",f10.3)
 42 format("Summ.  dQ30=" e10.3, " dQ32=",e10.3)
  
  5 format('Summ.    Jz=',f10.3, '  OmZ=',f10.3)
 51 format('Summ.    Jx=',f10.3, '  OmX=',f10.3)
  6 format('Summ.   dH2=',e10.3)
  ! Printing energy
  print 1, Iteration
  print 2, TotalEnergy,abs((TotalEnergy - OldEnergy(1))/TotalEnergy)
  print 21, Routhian, abs((Routhian - OldRouthian(1))/Routhian)
  print 22, Fermi
  !------------------------------------------------------------------
  ! Note that we can best look back two values due when there is a
  ! Rutz constraint active.
  Current => FindMoment(2,0,.false.) ; Q20 = sum(Current%Value) ; dQ20 = Q20 - sum(Current%OldValue(:,2))
  Current => FindMoment(2,2,.false.) ; Q22 = sum(Current%Value) ; dQ22 = Q22 - sum(Current%OldValue(:,2))

  if(.not. PC) then
  	Current => FindMoment(1,0,.false.) ; Q10 = sum(Current%Value) ; dQ10 = Q10 - sum(Current%OldValue(:,2))
 	Current => FindMoment(3,0,.false.) ; Q30 = sum(Current%Value) ; dQ30 = Q30 - sum(Current%OldValue(:,2))
  	Current => FindMoment(3,2,.false.) ; Q32 = sum(Current%Value) ; dQ32 = Q32 - sum(Current%OldValue(:,2))
  endif  

  if(.not.PC) then
	print 31, Q10
        print 41, dQ10
  endif

  print 3, Q20,Q22
  print 4, dQ20, dQ22

  if(.not.PC) then
	print 32,  Q30, dQ30
        print 42, dQ30, dQ32
  endif


  if(.not.TRC) then
    if((.not. TRC) .and. (.not.SC)) then
        print 51, TotalAngMom(1), Omega(1)
    endif
    print 5, TotalAngMom(3), Omega(3)        
  endif

  !------------------------------------------------------------------------------
  ! Calculate sum of dispersions
  dispersion = 0.0_dp
  do i=1,nwt
    dispersion = dispersion + HFBasis(i)%Occupation*HFBasis(i)%Dispersion
  enddo

  print 6, Dispersion

end subroutine PrintSummary_v2

subroutine PrintSummary_v1(Iteration)
	!-----------------------------------------------------------------------------
	! This subroutine prints a summary of the iteration.
	!
	! Included:
	!	- Relative change in energy
	!	- Relative (or absolute) change in constrained moments
	! - Angular momentum for cranked directions
	!-----------------------------------------------------------------------------
	use CompilationInfo
	use Energy, only: TotalEnergy, OldEnergy
	use Spwfstorage, only : TotalAngMom
	use Densities,   only : DensityChange
	use Cranking
	use Moments
	use GenInfo

	implicit none

	1 format ("Summ: dE=",  e10.3, ' ','dRho=',e10.3, ' ', a80 )
	2 format (' ',a1, 'Q', 2i1,'=',e10.3, ' ')
	3 format (' <',a2,'>=',e9.2, ' ')
    4 format (' dH2=',f10.7)

	integer, intent(in) :: Iteration
	type(Moment), pointer :: Current
	character(40)       :: MomString, CrankString
	character(len=18)   :: Temp=''
	character(len=1)    :: R
	integer             :: i
	real(KIND=dp)       :: MomentDiff,Dispersion
	character(len=2)    :: DirectionName(4) = (/ "Jx", "Jy", "Jz", "J2" /)

	Current => Root
	MomString=''
	!-----------------------------------------------------------------------------
	!Find Constrained Moments
	do while(associated(Current%Next))
		Current => Current%Next
		if(Current%ConstraintType .ne. 0) then
			MomentDiff=abs(sum(Current%Value) - sum(Current%OldValue(:,1)))
			if(sum(Current%Value) .gt. 1d-3) then
			  MomentDiff = MomentDiff/abs(sum(Current%Value))
		  endif
			R='R'
			if(Current%Impart) R='I'

			write(Temp, 2),R, Current%l, Current%m, MomentDiff
			MomString = adjustl(trim(MomString)//Temp)
		endif
	enddo

	!-----------------------------------------------------------------------------
	!Find cranked directions
	Temp=''
	CrankString=''
	do i=1,3
	  if(Omega(i).ne.0.0_dp) then
	    write(Temp,fmt=3), DirectionName(i), TotalAngMom(i)
	    CrankString = adjustl(trim(CrankString)//Temp)
	  endif
	enddo
  !------------------------------------------------------------------------------
  ! Calculate sum of dispersions
  dispersion = 0.0_dp
  do i=1,nwt
    dispersion = dispersion + HFBasis(i)%Occupation*HFBasis(i)%Dispersion
  enddo

  print 1, abs((TotalEnergy - OldEnergy(1))/TotalEnergy), DensityChange,       &
  &        adjustl(adjustr(MomString)//adjustl(CrankString))
	nullify(Current)
end subroutine PrintSummary_v1

