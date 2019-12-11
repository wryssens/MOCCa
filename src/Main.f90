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
    use Legacy
    use Transform
    
    implicit none

    real*8 :: time(2)
    integer :: i
    
 200 format (/,' ___________________________________________________________', &
     &       /,'|                                                          |', &
     &       /,'|  MOCCa                                                   |', &
     &       /,'|                                                          |', &
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

 299 format (  ' ----------------------------------------------------------')
 300 format (  ' Version information from GIT' )
 301 format (  '  VERSION1') ! Git commit
 302 format (  '  VERSION2') ! Author of commit
 303 format (  '  VERSION3') ! Date

     print 200
     
     print 299
     print 300
     print 301
     print 302
     print 303
     print 299

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
  ! This subroutine commands the evolution of the system.
  !
  ! The integer iprint is a flag: iprint=1 gives extra diagnostic printouts.
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
  use DensityMixing
  use Energy 
  use HFB
  use Testing

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

  integer             :: dummyHF(2,2,2,2) = 0
  logical, intent(in) :: iprint
  integer, intent(in) :: MaxIterations
  integer             :: Iteration
  logical             :: Convergence, AlternateCheck
  logical, external   :: ConvergenceCheck
  integer             :: i,wave
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

  !--------------------------------------------------------------------------
  !Updating certain parameters of the wavefunctions that were not read from file
  do i=1,nwt
    call HFBasis(i)%SymmetryOperators()
  enddo

  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------
  ! Derive all the spwfs to start the calculation
  if(t1n2.ne. 0.0_dp .or. t2n2.ne.0.0_dp)   then
    if(MAXFDORDER .ne. -1 .or. MAXFDLAPORDER .NE. -1) then
      call stp('It is irresponsible to calculate N2LO without lag derivatives.')
    endif
    call N2LODerive()
  else
    call DeriveAll()  
  endif

  !-----------------------------------------------------------------------------
  !Update the Angular momentum variables
  call UpdateAm()
  !Make sure that there is no improper readjusting of cranking constraints
  AngMomOld = TotalAngMom
  !-----------------------------------------------------------------------------
  ! Recalculate densities if needed
  ! This includes several scenarios:
  ! 1) If the user asked for it
  ! 2) if the input file was a CR8 file, the time-odd densities could then
  !    not be read from file.
  ! 3) if the input file was an EV8 file, some derivatives of densities could 
  !    not be read.
  if(Recalc) then
    call UpdateDensities(0)
  endif
  ! We calculate all of the multipole moments three times. This is strictly 
  ! wasteful, but only done at the zeroth iteration. This is to make sure
  ! the correct values are used in all readjustment formulas.
  call CalculateAllMoments(1)
  call CalculateAllMoments(1)
  call CalculateAllMoments(1)
  !Construct the mean-field potentials
  call ConstructPotentials
  !Calculating  The Energy
  call CompEnergy()
  ! Explicitly calculate <h> for the canonical basis
  if(Pairingtype.eq.2) then
    do i=1,nwt
      CanEnergy = InproductSpinorReal(CanBasis(i)%GetValue(),hPsi(Canbasis(i)))
      call Canbasis(i)%SetEnergy(CanEnergy)
    enddo
  endif

  !-----------------------------------------------------------------------------
  ! Printing the parameters of this run.
  call PrintStartInfo

  !Loop over the iterations
  Iteration   = 0
  Convergence = .false.

  !Printing observables
  if(maxiterations .eq.0) print 101
  call PrintIterationInfo(0, .true.)

  !Checking for the presence of constraints with extra steps
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
    ! Take a preliminary step to improve the constraints
    if(AlternateCheck .or. AlternateCrank) then
      call AlternateStep()
    endif
    !---------------------------------------------------------------------------
    ! Evolve the spwfs with the chosen strategy.
    ! Note that this includes orthonormalization.
    call EvolveSpwf(Iteration)    
    !---------------------------------------------------------------------------
    ! Solve the pairing subproblem in the new HFbasis.
    call SolvePairing
    !---------------------------------------------------------------------------
    ! Checking for the presence of alternate constraints, they may have been 
    ! switched off!
    AlternateCheck = CheckForAlternateConstraints()
    !----------------------------------------------------------------------------
    !Deriving all Spwfs
    if(t1n2.ne. 0.0_dp .or. t2n2.ne.0.0_dp) then
        call N2LODerive()
    else
        call DeriveAll()    
    endif
    do i=1,nwt
      call HFBasis(i)%SymmetryOperators()
      if(Pairingtype.eq.2) then
        CanEnergy = InproductSpinorReal(CanBasis(i)%GetValue(), hPsi(Canbasis(i)))
        call Canbasis(i)%SetEnergy(CanEnergy)
        call CanBasis(i)%SymmetryOperators()
      endif
    enddo

    !----------------------------------------------------------------------------
    ! Update all calculated quantities
    call UpdateDensities(0)
    ! Multipole moments 
    call CalculateAllMoments(1) 
    call ReadjustAllMoments()
    ! Angular momentum
    call updateAm()
    call ReadjustCranking()
    ! Smooth the densities using the difference between this and the previous
    ! iteration, if the user asked for it
    call MixDensities(Iteration)
    !See if some constraints were temporary
    call TurnOffConstraints(iteration)
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
    !Print a summary
    if(mod(Iteration, PrintIter).ne.0) call PrintSummary(Iteration)
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
  endif

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
  use InOutput, only: Pictures, PlotDensity
  use Force, only   : IsTrapped    
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
  if( TotalEnergy.ge.0.0_dp .and. ( .not. IsTrapped) ) then
    call PrintIterationInfo(-1, .true.)
    !Write densities to files.
    if(Pictures) then
       call PlotDensity()
    endif
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
    use Cranking,only  : PrintCranking, UpdateAM
    use Pairing, only  : PrintPairingInfo
    use Spwfstorage, only:  DeriveAll
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

    call UpdateAm()
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
  use Coulomb
  use ImaginaryTime 
  use geninfo

  implicit none

  integer, intent(in) :: Iteration
  logical, intent(in) :: PrintAll

  1 format (/,60('='))
  2 format (' Iteration ', i5, /, ' (Average) dt = ', f10.5, ' mu =', f10.5)

  print 1
  if(Iteration.eq.0) then 
    print 2, Iteration
  else
    print 2, Iteration, dt, momentum
  endif

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
  use Spwfstorage, only : DeriveAll, DensityBasis, nwt, NumberParityCounting
  use Densities,   only : UpdateDensities
  use Pairing,     only : SolvePairing
  use InOutput,    only : Pictures, PlotDensity, PlotCurrents
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

  ! Count the quantum numbers of the single-particle states
  if((.not. PC) .and. (.not. SC) .and. (.not. TRC))call NumberParityCounting

  ! Recalculating the finite difference derivatives, in order
  ! to write the densities of the final iteration to file
  call AssignFDCoefs(MaxFDOrder, MaxFDLapOrder, CoulombLapOrder)
  do i=1,nwt
    call DensityBasis(i)%CompDer()
  enddo

  call UpdateDensities(0)

end subroutine FinalIteration

subroutine PrintSummary_v2(Iteration)
  !-----------------------------------------------------------------------------
  ! This routine prints a summary of the iteration.
  ! Multilined, as I like to have plenty of info.
  !
  !-----------------------------------------------------------------------------

  use CompilationInfo
  use Energy, only: TotalEnergy, OldEnergy, Routhian, OldRouthian, LNEnergy
  use Spwfstorage, only : TotalAngMom
  use Densities,   only : DensityChange, Density
  use Pairing, only : Fermi, LNLambda, Lipkin
  use Cranking
  use Moments
  use GenInfo

  implicit none

  integer, intent(in)  :: Iteration
  type(Moment),pointer :: Current
  real(KIND=dp) :: Q20, Q22, dQ20, dQ22, Dispersion, Q30, Q32, dQ30, dQ32, Q10, dQ10
  real(KIND=dp), save :: oldln2(2), diff(2), oldlne(2), oldlambda(2), oldce
  integer       :: i

  1 format('----------- Iteration ', i6,'--------------')
  2 format("Summ.     E=" f10.3, "   dE=",es10.3)
 21 format("Summ.     R=" f10.3, "   dR=",es10.3)
 22 format("Summ. Fermi=" f10.5, "      ",f10.5) 
222 format("Summ.    dF="es10.3, "      ",es10.3)
 23 format("Summ.    L2=" f10.5, "      ",f10.5) 
231 format("Summ.   dL2=" es10.3, "      ",es10.3) 
232 format("Summ.  dEL2=" es10.3, "      ", es10.3)  

  3 format("Summ.   Q20=" f10.3, "  Q22=",f10.3)
  4 format("Summ.  dQ20=" es10.3, " dQ22=",es10.3)
 31 format("Summ.   Q10=" f10.3)
 41 format("Summ.  dQ10=" es10.3)
 32 format("Summ.   Q30=" f10.3, "  Q32=",f10.3)
 42 format("Summ.  dQ30=" es10.3, " dQ32=",es10.3)
  
 43 format("Summ.   dCE=" es10.3 )

  5 format('Summ.    Jz=',f10.6, '  OmZ=',f10.6)
 51 format('Summ.    Jx=',f10.6, '  OmX=',f10.6)
 52 format('Summ. ReJxT=',f10.6)
  6 format('Summ.   dH2=',es10.3)
  ! Printing energy
  print 1, Iteration
  print 2, TotalEnergy,abs((TotalEnergy - OldEnergy(1))/TotalEnergy)
  print 21, Routhian, abs((Routhian - OldRouthian(1))/Routhian)
  print 22, Fermi
  print 222, Fermi - oldlambda
  if ( Lipkin ) then
     print 23, LNLambda
     diff = LNLambda - oldln2
     print 231,diff
     print 232, abs(LNEnergy-oldlne)/abs(TotalEnergy)
     oldln2    = lnlambda
     oldlne    = lnenergy
  endif
  oldlambda = Fermi

  !-----------------------------------------------------------------------------
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
	print 32,  Q30, Q32
        print 42, dQ30, dQ32
  endif

  print 43, abs(oldce - sum(ConstraintEnergy*Density%Rho)*dv )/abs(TotalEnergy)
  oldce = sum(ConstraintEnergy*Density%Rho)*dv 

  if(.not.TRC) then
    if((.not. TRC) .and. (.not.SC)) then
        print 51, TotalAngMom(1), Omega(1)
    else
        print 52,JTR(1)
    endif
    print 5, TotalAngMom(3), Omega(3)
  endif

  !------------------------------------------------------------------------------
  ! Calculate sum of dispersions
  dispersion = 0.0_dp
  do i=1,nwt
    dispersion = dispersion + HFBasis(i)%Occupation*HFBasis(i)%Dispersion
  enddo

  print 6, abs(Dispersion)
  ! Saving
  TotalDispersion = Dispersion

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
  dispersion = dispersion/(neutrons + protons)

  print 1, abs((TotalEnergy - OldEnergy(1))/TotalEnergy), DensityChange,       &
  &        adjustl(adjustr(MomString)//adjustl(CrankString))
	nullify(Current)
end subroutine PrintSummary_v1

