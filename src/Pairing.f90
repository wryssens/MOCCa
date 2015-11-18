module Pairing
!-------------------------------------------------------------------------------
! Module that is able to solve the pairing problem.
!
!-------------------------------------------------------------------------------
! NOTES
! 1) This module works 'alone' in the sense that other modules reference it only
!    very sparingly. This module only supplies
!    a) SolvePairing pointer routine that actually solves the pairing problem.
!    b) Occupation probabilities and quasi-particle energies in either the
!       HF-basis or canonical basis. They are stored in the respective modules,
!       but need to be set in this module.
! 2) This module uses separate modules to actually solve the pairing problem.
!    In this sense, this module is only used to drive the other modules.
! 3) This contains the basic HF solver too, for simplicity.
! 4) Same remarks regarding precision of conversion with complex numbers as in 
!    the HFB module.s
! 5) The BCS and HFB module both have several routines that take arguments
!    that are not strictly needed. They are only needed in order to be able
!    to abstract the pairing solver and be able to point some procedure 
!    pointers.
!-------------------------------------------------------------------------------

  use CompilationInfo
  use Geninfo
  use SpwfStorage
  use BCS
  use HFB
  use PairingInteraction
  use LipkinNogami

  implicit none

  !-----------------------------------------------------------------------------
  ! Integer determines the type of pairing that is active
  !  0  =>  HF 
  !  1  =>  BCS
  !  2  =>  HFB
  integer :: PairingType=0
  !-----------------------------------------------------------------------------
  ! String for inputting pairingtype, for ease of use.
  ! 'HF ' => Pairingtype 0
  ! 'BCS' => PairingType 1
  ! 'HFB' => PairingType 2
  character(len=3) :: Type='   '
  !-----------------------------------------------------------------------------
  ! Different pointers to important parts of solving the pairing equations.
  procedure(BCSGaps), pointer                     :: GetGaps          !Pairing Gaps
  procedure(BCSEnergy), pointer                   :: CompPairingEnergy!Calculate energy
  procedure(BCSPairingField), pointer             :: GetPairingFields !Pairing fields
  procedure(BCSFindFermiEnergy), pointer          :: FindFermiEnergy  !Fermi Energy
  procedure(BCSOccupations), pointer              :: GetOccupations   !Occupations
  procedure(HFBComputePairingDensity), pointer    :: CompPDensity     
  !-----------------------------------------------------------------------------
  ! Universal properties
  !-----------------------------------------------------------------------------
  ! Dispersion of the particle numbers
  ! < N^2 > - < N >^2
  real(KIND=dp) :: PairingDisp(2) = 0.0_dp
  ! Energy associated with pairing
  real(KIND=dp) :: PairingEnergy(2)=0.0_dp
  ! Fermi levels and their history
  real(KIND=dp) :: Fermi(2)=-10.0_dp, FermiHistory(2,7)=0.0_dp
  ! Freeze the occupation or not.
  logical       :: FreezeOccupation=.false.  
  !-----------------------------------------------------------------------------
  ! Pairing gaps for the different pairing models. 
  !-----------------------------------------------------------------------------
  ! Indices for this quantity are:
  ! 1) wavefunction index
  ! 2) wavefunction index in the case of HFB, otherwise dimension 1
  ! 3) Parity index in case of HFB, otherwise dimension 1
  ! 4) Isospin index in case of HFB, otherwise dimension 1
  complex(KIND=dp), allocatable :: Delta(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! Pairing fields, used to calculate the gaps.
  ! Defined as a position-dependent Delta.
  ! Delta(x,y,z) = 
  !    Sum_{i,j}     2*fi*fj*kappa_{ij} * s * <r,s;r,-s|V|ij>
  ! which is general for HFB and reduces to the following for BCS:
  ! Delta(x,y,z) = 
  !    Sum_{i,ibar}  2*fi*fibar*u_i*v_i * s * <r,s;r,-s|V|ij>
  !-----------------------------------------------------------------------------
  complex(KIND=dp), allocatable  :: PairingField(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! Modified Lipkin-Nogami pairing field and pairing gaps
  complex(KIND=dp), allocatable  :: PairingFieldLN(:,:,:,:), DeltaLN(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! Integer for the maximum number of iterations for the readjustment of lambda
  integer                     :: MaxLambdaIter=1
  !-----------------------------------------------------------------------------
  ! Precision on the Fermi energy
  real(KIND=dp), parameter    :: FermiPrec=1d-6
  !-----------------------------------------------------------------------------
  ! Maximum iterations for the pairing solver inside one mean-field iteration.
  integer                     :: PairingIter=1
  !-----------------------------------------------------------------------------
  ! Value of the pairing gaps for constant pairing gap calculations and logical
  ! whether to use them or not.
  real(KIND=dp)               :: Gaps(2) = 0.0_dp
  logical                     :: ConstantGap=.false.
  !-----------------------------------------------------------------------------
  ! Lipkin-Nogami switch and Lipkin-nogami multiplier Lambda_2.
  logical                     :: Lipkin =.false.
  real(KIND=dp)               :: LNLambda(2)=0.0_dp
  !-----------------------------------------------------------------------------
  ! Logical that determines whether a reinitialisation of KappaHFB is done 
  ! for HFB calculations.
  logical                     :: GuessKappa=.false.
  !-----------------------------------------------------------------------------
  ! Whether or not the pairing model is solved before the start of iterations.
  ! Thus, if .false., the occupations from file are used.
  logical                     :: SolvePairingStart=.true.
  !-----------------------------------------------------------------------------
  ! Integer that tells MOCCa what number of blocked qps to look for.
  integer                     :: Block=0
  !------------------------------------------------------
  ! Temporary array for the pairing density rho~, until a 
  ! more general implementation is available.
  complex(KIND=dp),allocatable :: PairDensity(:,:,:,:)
  ! Response pairing density
  real(KIND=dp),allocatable    :: chi(:,:,:,:)
  ! Adjoint Response pairing density
  real(KIND=dp),allocatable    :: chistar(:,:,:,:)

contains
 
  subroutine ReadPairingInfo
    !---------------------------------------------------------------------------
    ! Subroutine that reads info on the pairing interaction from the namelist
    ! Pairing. It then checks as far as possible the input on consistency errors
    ! and allocates the necessary variables.
    !---------------------------------------------------------------------------
     real(KIND=dp) :: PairingNeutron=-100.0_dp, PairingProton=-100.0_dp
     real(KIND=dp) :: CutProton=-100.0_dp, CutNeutron=-100.0_dp
     real(KIND=dp) :: AlphaProton=-100.0_dp, AlphaNeutron=-100.0_dp
     real(KIND=dp) :: Protongap=0.0_dp, NeutronGap=0.0_dp
     integer       :: i
     logical       :: SemiBCS=.false., SemiBCSNeutron=.false.
     logical       :: SemiBCSProton=.false.
     character(len=3) :: UpperType

  
     NameList /Pairing/ PairingNeutron, PairingProton, CutProton , RhoSat,     &
     &                  alphaProton, AlphaNeutron,                             &
     &                  Type, CutNeutron, FreezeOccupation, PairingIter,       &
     &                  HFBMix, NeutronGap, ProtonGap, ConstantGap,            &
     &                  Lipkin, LNFraction, GuessKappa, PairingMu, CutType,    &
     &                  MaxLambdaIter, SemiBCS, SemiBCSNeutron, SemiBCSProton, &
     &                  QPinHFBasis, SolvePairingStart,QPPrintWindow, Block,   &
     &                  FermiSolver

     read(unit=*, NML=Pairing)
  
     PairingStrength(1) = PairingNeutron;  PairingStrength(2) = PairingProton   
     PairingCut(1)      = CutNeutron    ;  PairingCut(2)      = CutProton
     Alpha(1)           = AlphaNeutron  ;  Alpha(2)           = AlphaProton
     Gaps(1)            = NeutronGap    ;  Gaps(2)            = ProtonGap
    
     call to_upper(Type, Type)
     select case (Type)
     case ('HF ')     
      PairingType = 0
     case ('BCS')
      PairingType = 1
     case ('HFB')
      PairingType = 2
     case DEFAULT
      PairingType = 0 !Default to HF
     end select
     call to_upper(FermiSolver,FermiSolver)
     
     !-----------------------------------------------------------------------------
     !Checking input
     if(PairingStrength(2).eq.-100.0_dp) then 
      if(PairingStrength(1).ne.-100.0_dp) then
        PairingStrength(2) = PairingStrength(1)
      else
        PairingStrength=0.0_dp
      endif
     endif
     if(PairingStrength(1).eq.-100.0_dp) then 
      if(PairingStrength(2).ne.-100.0_dp) then
        PairingStrength(1) = PairingStrength(2)
      endif
     endif
     if(PairingCut(2).eq.-100.0_dp) then 
      if(PairingCut(1).ne.-100.0_dp) then
        PairingCut(2) = PairingCut(1)
      else
        PairingCut = 5.0_dp
      endif
     endif
     if(PairingCut(1).eq.-100.0_dp) then 
      if(PairingCut(2).ne.-100.0_dp) then
        PairingCut(1) = PairingCut(2)
      endif
     endif
      if(Alpha(2).eq.-100.0_dp) then 
      if(Alpha(1).ne.-100.0_dp) then
        Alpha(2) = Alpha(1)
      else
        Alpha = 0.0_dp
      endif
     endif
     if(Alpha(1).eq.-100.0_dp) then 
      if(Alpha(2).ne.-100.0_dp) then
        Alpha(1) = Alpha(2)
      endif
     endif

     if(PairingType.ne.2 .and. Block .ne. 0) then
      call stp('Cannot block quasiparticles when not doing HFB calculations.')
     endif

     if(Pairingtype.lt.0 .or. PairingType.gt.2) then
      call stp('Non-defined pairing treatment chosen. Choose either 0,1 or 2.',&
      &        'PairingType', PairingType)
     endif
     if(RhoSat .lt. 0.0_dp ) then
      call stp('Invalid value for RhoSat. ', 'RhoSat', Rhosat)
     endif
     if(PairingCut(1).le. 0.0_dp)  then
      call stp('Choose a positive value for neutron cutoff.',                  &
      &        'PairingCut',PairingCut(1))
     endif
     if(PairingCut(2).le. 0.0_dp)  then
      call stp('Choose a positive value for proton cutoff.',                   &
      &        'PairingCut',PairingCut(2))
     endif
     if(alpha(1) .lt. 0.0_dp .or. alpha(1) .gt. 1.0_dp) then
      call stp('You should probably choose alpha for neutrons between 0 and 1.'&
      &       ,'alpha', alpha(1))
     endif
     if(alpha(2) .lt. 0.0_dp .or. alpha(2) .gt. 1.0_dp) then
      call stp('You should probably choose alpha for protons between 0 and 1.' &
      &       ,'alpha', alpha(2))
     endif
     if( PairingStrength(1) .lt. 0.0_dp ) then
      call stp('Neutron pairing strength is not positive.',                    &
      &        ' Neutron pairing strength' , PairingStrength(1))
     endif
     if( PairingStrength(2) .lt. 0.0_dp ) then
      call stp('Proton pairing strength is not positive.',                     &
      &        ' Proton pairing strength' , PairingStrength(1))
     endif
     if(.not. TRC .and. PairingType .eq. 1) then
      !Don't allow BCS when time-reversal is not conserved!
      call stp("MOCCa can't do BCS when Time Reversal is not conserved!")
     endif
     if(ConstantGap .and. Lipkin) then
        call stp('Impossible to enforce constant pairing gaps when'            &
        &         // ' Lipkin-Nogami is enabled.')
     endif
     if(Block.lt.0) then
        call stp('Number of blocked quasiparticles should not be smaller than' &
          &    //' zero.')
     endif
     !-----------------------------------------------------------------------------
     !--------------------------------------------------------------------------
     ! Determine the cutoff type.
     select case(Cuttype)    
     case(1) 
       PairingCutoff => SymmetricFermi
     case(2)
       PairingCutoff => CosineCut
     case DEFAULT
      call stp('Non-recognised pairing cutoff type.')
     end select
     !-----------------------------------------------------------------------
     ! Determine the sizes of the matrices.
     ! Note that this parameter is necessary, independent of whether HFB
     ! is activated or not, since there is always a guess for Kappa that is
     ! written to file.
     if(TRC) then
        !Double the size if there is Time-reversal present
        HFBSize = 2*nwt
     else
        HFBSize = nwt
     endif
     !--------------------------------------------------------------------------
     !Assigning the correct subroutine to SolvePairing
     nullify(CompPairingEnergy); nullify(GetOccupations) ; nullify(GetGaps) 
     nullify(GetPairingFields) ; nullify(FindFermiEnergy)
     select case (PairingType) 
      case(0)
        !-----------------------------------------------------------------------
        ! Hartree-Fock
        CompPairingEnergy => HFEnergy
        !-----------------------------------------------------------------------
        ! Double-check neutron & proton number.
        if((int(Neutrons) - Neutrons).ne.0.0_dp) then
          call stp('Number of neutrons should be integer for HF calculations.')
        endif
        if((int(Protons) - Protons).ne.0.0_dp) then
          call stp('Number of protons should be integer for HF calculations.')
        endif
        
      case(1)
        if(Lipkin) call stp('Lipkin-Nogami is not implemented for BCS yet.')
        !-----------------------------------------------------------------------
        ! BCS Pairing!
        CompPairingEnergy => BCSEnergy
        GetPairingFields  => BCSPairingField
        GetOccupations    => BCSOccupations
        GetGaps           => BCSGaps
        FindFermiEnergy   => BCSFindFermiEnergy
        !-----------------------------------------------------------------------
        allocate(Delta(nwt,1,1,1)) ; allocate(PairingField(nx,ny,nz,2))
        Delta =0.0_dp ; PairingField = 0.0_dp        
        !Put the default value to 50 for BCS calculations
        if(PairingIter.eq.1) then
            PairingIter=50
        endif
      case(2)
        !-----------------------------------------------------------------------
        ! HFB pairing!
        !-----------------------------------------------------------------------
        ! Parsing SemiBCS input
        if(SemiBCS) then
          BCSinHFB = .true.
        elseif(SemiBCSNeutron) then
          BCSinHFB(1) = .true.
        elseif(SemiBCSProton) then
          BCSinHFB(2) = .false.
        else
          BCSinHFB = .false.
        endif      
        !Make sure the density basis points to the correct spot.
        DensityBasis      => CanBasis
        GetPairingFields  => HFBPairingField
        GetOccupations    => HFBOccupations
        GetGaps           => HFBGaps!_TIMEREV!HFBGaps
        CompPDensity      => HFBComputePairingDensity

        if(trim(FermiSolver).eq.'BROYDEN') then
          FindFermiEnergy   => HFBFindFermiEnergyBroyden
        elseif(trim(FermiSolver).eq.'BISECTION') then
          FindFermiEnergy   => HFBFindFermiEnergyBisection
        else
          call stp('Unknown Fermi solver requested.')
        endif

        CompPairingEnergy => HFBEnergy
        !-----------------------------------------------------------------------
        !Prepare the indices in the module
        allocate(Delta(HFBSize,HFBSize,2,2)); allocate(PairingField(nx,ny,nz,2))
        Delta =0.0_dp                               ;  PairingField = 0.0_dp
        allocate(RhoHFB(HFBSize,HFBSize,2,2))       ;  RhoHFB       = 0.0_dp
        allocate(KappaHFB(HFBSize,HFBSize,2,2))     ;  KappaHFB     = 0.0_dp
        allocate(OldRhoHFB(HFBSize,HFBSize,2,2))    ;  OldRhoHFB       = 0.0_dp
        allocate(OldKappaHFB(HFBSize,HFBSize,2,2))  ;  OldKappaHFB     = 0.0_dp     
        !-----------------------------------------------------------------------
        ! When Lipkin-Nogami is present, make space for the modified pairing
        ! gaps.
        if(Lipkin) then
          allocate(DeltaLN(HFBSize,HFBSize,2,2))   ; DeltaLN = 0.0_dp
          allocate(PairingFieldLN(nx,ny,nz,2))     ; PairingFieldLN = 0.0_dp
          allocate(PairDensity(nx,ny,nz,2), Chi(nx,ny,nz,2), ChiStar(nx,ny,nz,2))
        endif
        !Note that we will store the complete U & V matrices     
        allocate(U(HFBSize,2*HFBSize,2,2)) ; allocate(V(HFBSize,2*HFBSize,2,2))
        U = 0.0_dp ; V = 0.0_dp
        allocate(HFBHamil(2*HFBSize,2*HFBSize,2,2)) ; HFBHamil=0.0_dp
        !-----------------------------------------------------------------------
        ! By default, consider even-even nuclei. 
        ! The routine readblocking will change the numberparity if necessary.
        HFBNumberParity=1        
        !-----------------------------------------------------------------------
        !Reset the canonical wavefunctions
        allocate(CanBasis(nwt))
        do i=1,nwt
            call CanBasis(i)%Resetwf()
        enddo     
        !-----------------------------------------------------------------------       
        !Setting numerical default parameters
        if(PairingIter.eq.1) then
          PairingIter=1
        endif
        if(MaxLambdaIter.eq.1) then
          MaxLambdaIter=50
        endif
        !-----------------------------------------------------------------------
        ! Calling the readBlocking routine when the user wants to do blocking
        if(Block.ne.0) call ReadBlockingInfo(Block)
     end select
  
  end subroutine ReadPairingInfo

  subroutine PrintPairingInfo
    !---------------------------------------------------------------------------
    ! Subroutine prints pairing info to STDOUT at the start of the program.
    !---------------------------------------------------------------------------
    
    1  format (22('-'), 'Pairing Parameters', 21('-'))
    2  format (A26) 
    3  format (33x, ' N ',7x, ' P ') 
    4  format ('  Pairing Strength (MeV fm^3)', 2x,f8.3,2x,f8.3)
    5  format ('  Pairing Cutoff   (Mev)     ', 2x,f8.3,2x,f8.3)
    6  format ('  Pairing Cut width(MeV)     ', 2x,f8.3,2x,f8.3)
    7  format ('  Alpha                      ', 2x,f8.3,2x,f8.3)
    8  format ('  Attention: configuration is frozen.')
    9  format ('  HFBMixing                  ', 2x, f8.3)
   10  format ('  BCS-like: only Delta_{i,ibar} not zero.' )
   100 format ('  BCS-like Protons : only Delta_{i,ibar} not zero.' )
   101 format ('  BCS-like Neutrons: only Delta_{i,ibar} not zero.' ) 
   11  format ('  Constant Gap activated. ', /,                                 &
   &          '   Neutron Gap: ' , f8.3   , /,                                  &
   &          '   Proton  Gap: ' , f8.3   )
   12  format ('  Lipkin-Nogami prescription active.' )
   13  format ('  Cosine cutoff used. ')
   14  format ('  Symmetric Fermi cutoff used.')
   15  format ('  Maximum number of times the FermiSolver is called: ', i3)
   16  format ('  LNFraction on the hamiltonian:' , f8.3)
   17  format ('  Fermisolver used: ', a9)

    print 1
    select case(PairingType)
      case(0)
        print 2, 'Hartree-Fock; No pairing. '
        if(FreezeOccupation) print 8 
        return !Done printing when dealing with hartree-fock.
      case(1)
        print 2, 'BCS with zero-range force.'
        if(ConstantGap) print 11, Gaps
      case(2)
        print 2, 'HFB with zero-range force.'
        if(all(BCSinHFB)) then
          print 10
        elseif(BCSinHFB(1)) then
          print 101
        elseif(BCSinHFB(2)) then
          print 100
        endif
    end select
    
    select case(CutType)
    case(1) 
      print 14
    case(2) 
      print 13
    end select
    
    print 3
    print 4, PairingStrength
    print 5, PairingCut
    print 6, PairingMu
    print 7, Alpha
    
    if(PairingType.eq.2) print 9, HFBMix
    if(Lipkin)           print 12
    
    print *
    print 15 , PairingIter
    print 16 , LNFraction
    if(PairingType.eq.2) print 17, FermiSolver
    
    print * 

    if(allocated(QPExcitations)) call PrintBlocking

  end subroutine PrintPairingInfo

  subroutine PrintPairing
    !---------------------------------------------------------------------------
    ! Subroutine that prints intermediate info on the pairing
    !
    !---------------------------------------------------------------------------
  
    use Densities, only : Density
  
    1 format (27('-'), 'Pairing', 26('-'))
    2 format (25x, ' N ',7x, ' P ') 
    3 format (' Fermi Level (MeV) ',2x,f10.5,2x,f10.5)
    4 format (' Particles         ',2x,f10.5,2x,f10.5)
    5 format (' Dispersion        ',2x,f10.5,2x,f10.5)
    6 format (' Lipkin-N. Lambda_2',2x,f10.5,2x,f10.5)
    7 format (60('-'))
    
    select case(PairingType)   
      case (0) 
        !Don't print anything when doing HF
        return 
      case (1)
        !Print all that is below for BCS pairing
    end select
    
    print 1
    print 2
    print 3, Fermi
    print 4, sum(Density%Rho(:,:,:,1))*dv, sum(Density%Rho(:,:,:,2))*dv
    print 5, PairingDisp
    if(Lipkin) print 6, LNLambda
    print 7
  end subroutine PrintPairing

  function HFEnergy( Delta) result(E)
    !---------------------------------------------------------------------------
    ! Returns the pairing energy of a HF calculation, namely 0.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: E(2) 
    complex(KIND=dp), allocatable,intent(in) :: Delta(:,:,:,:)
  
    E = 0.0_dp
  
    return  
  end function HFEnergy
  
  subroutine HFFill()
    !---------------------------------------------------------------------------
    ! This subroutine finds the orbitals with the lowest single particle 
    ! energy and fills them, after sorting all the levels.
    ! Of course, nothing happens if the user has asked to freeze the occupation.
    !---------------------------------------------------------------------------
    integer :: i, n,p, ProtonUpperBound, NeutronUpperBound, order(nwt),j
    
    if(FreezeOccupation) return
    
    n=0; p=0
    !---------------------------------------------------------------------------
    !Setting all occupation numbers to Zero
    do i=1,nwt
            call HFBasis(i)%SetOcc(0.0_dp)                     
    enddo
    !---------------------------------------------------------------------------
    ! We can distribute the neutrons and protons in pairs in the case of
    ! Time Reversal Invariance
    if(TRC) then
      ProtonUpperBound  = floor(Protons/2)
      NeutronUpperBound = floor(Neutrons/2)
    else
      ProtonUpperBound  = floor(Protons)
      NeutronUpperBound = floor(Neutrons)
    endif
    !---------------------------------------------------------------------------
    !Finding the order of the spwfs, in terms of energy
    order = OrderSpwfs(.false.)
    !---------------------------------------------------------------------------
    !Filling in the lowest Proton orbitals. This is easy, since we know 
    !the order of Spwfs.
    i=1
    do while(p.lt.ProtonUpperBound .and. i.le.nwt)
       j = order(i)

      if(HFBasis(j)%GetOcc().eq.0.0_dp .and. HFBasis(j)%GetIsospin().eq.1) then
      !We've found an unoccupied Proton orbital
          call HFBasis(j)%SetOcc(1.0_dp)
          p = p + 1
      endif       
      i = i + 1     
    enddo
    !---------------------------------------------------------------------------
    !Filling in the lowest Neutron orbitals. This is easy,since we know the
    ! order of Spwfs.
    i=1
    do while(n.lt.NeutronUpperBound .and. i.le.nwt)
      j = order(i)
      if(HFBasis(j)%GetOcc().eq.0.0_dp .and. HFBasis(j)%GetIsospin().eq.-1) then
        !We've found an unoccupied Neutron orbital
        call HFBasis(j)%SetOcc(1.0_dp)
        n = n + 1
      endif       
      i = i + 1     
    enddo
    !---------------------------------------------------------------------------
    !Double all occupation Numbers in the case of Time Reversal Symmetry
    if(TRC) then
       do i=1,nwt
         call SetOcc(HFBasis(i),2.0_dp * HFBasis(i)%GetOcc())
       enddo
    endif
    !---------------------------------------------------------------------------
    !Putting the FermiEnergy equal to the energy of the highest occupied state
    do i=nwt,1,-1      
      j = order(i)               
      FermiEnergy = HFBasis(j)%GetEnergy()
      
      if(HfBasis(j)%GetOcc().ne.0.0_dp) exit
    enddo
    return
  end subroutine HFFill

  subroutine HFFermi
  !-----------------------------------------------------------------------------
  ! Find the fermi energy of a HF calculation. Trivial, since it is just the 
  ! energy of the highest occupied level.
  !-----------------------------------------------------------------------------
    integer       :: i, it
    real(KIND=dp) :: E
  
    Fermi = -1000000.0_dp
    do i=1,nwt
        it = (HFBasis(i)%GetIsospin() + 3)/2
        E  = HFBasis(i)%GetEnergy()
        if(HFBasis(i)%GetOcc().gt.0.0_dp.and.Fermi(it) .lt. E)then
            Fermi(it) = E
        endif
    enddo
  
  end subroutine HFFermi

  subroutine SolvePairing
  !-----------------------------------------------------------------------------
  ! Driver routine for the solving of the pairing equations, whatever the model
  ! chosen is exactly.
  !
  ! |--------- Iterate until the Fermi energy is stationary
  ! |          (Note that this iteration is not strictly necessary)
  ! |   1) Determine the pairing fields Delta(x,y,z)
  ! |   2) Determine the pairing gaps Delta_{ij} using the pairing fields.
  ! | |-------  Iterate until the particle number is correct
  ! | |         and the second Kamlah moment is stationary (when LN is active)
  ! | | 3) Guess Fermi energy (& second Kamlah moment when doing LN)
  ! | | 4) Find number of particles for current guess
  ! | | 5) Readjust Lambda ( & second Kamlah moment )
  ! | |-------
  ! |   6) Transform to canonical basis in case of HFB
  ! |   7) Find occupation probabilities for the Density basis
  ! |     (HFBasis for BCS, Canonical for HFB)
  ! |---------
  !-----------------------------------------------------------------------------
    1 format(' Warning! ',/,                                       &
    &        ' Pairing calculations did not converge.',/,          &
    &        ' Pairing iterations:  ', i4,/,                       &
    &        ' Old Fermi Energy:    ', 2f10.6,/,                   &
    &        ' New Fermi Energy:    ', 2f10.6)
    
    real(KIND=dp)    :: OldFermi(2), Particles(2)
    integer          :: outeriter
    !logical          :: ForceBisection=.false.
  
    if(PairingType.eq.0) then
     call HFFill()
     call HFFermi()
     return
    endif

    Particles(1) = Neutrons ; Particles(2) = Protons
    !Initialise the HFB indices properly
    if(.not.allocated(blockindices) .and. PairingType.eq.2) then
      call PrepareHFBModule
    endif
    !---------------------------------------------------------------------------
    ! Reinitialise Kappa when asked for
    if(GuessKappa .and. PairingType.eq.2) then 
      call WriteOutKappa(PairingType)
      GuessKappa = .false.
    endif
    !---------------------------------------------------------------------------
    ! Make sure the HFBHamiltonian at iteration 0 can be constructed
    if(PairingType.eq.2) then
      if(all(OldRhoHFB.eq.0.0_dp)) then
        OldRhoHFB   = RhoHFB
        OldKappaHFB = KappaHFB
      endif
    endif
    if(PairingType.eq.2) DensityBasis => CanBasis
    
    if(PairingType.eq.1 .and. maxval(abs(Delta)).lt.1d-5) Delta = 1.0_dp
    if(Lipkin)  where(LNLambda  .eq.0.0_dp) LNLambda= 0.1_dp
    if(.not.Lipkin) LNLambda= 0.0_dp
    !---------------------------------------------------------------------------
    ! Outer iterations: iterating gaps and pairingfield.
    do outerIter=1,PairingIter
            
        call ComputePairingCutoffs(Fermi)
        if(.not.ConstantGap) then
          !Don't get new pairingfields if the gaps are constant
          call GetPairingFields(PairingField,PairingFieldLN,Delta)
        endif
        call GetGaps(Delta,DeltaLN,PairingField,PairingFieldLN,Gaps,ConstantGap)
           
        OldFermi = Fermi           
        !-----------------------------------------------------------------------
        ! Get the correct Fermi energy and LN parameter          
        call FindFermiEnergy(Fermi,LNLambda,Delta,DeltaLN,Lipkin,FermiPrec)
        call ComputePairingCutoffs(Fermi)
        !-----------------------------------------------------------------------
        ! Get the occupations in the appropriate basis and as a side effect
        ! calculate the dispersion in the particle number.
        call GetOccupations(Fermi,Delta,LNLambda,PairingDisp)  
                
        if(all(abs(Fermi - OldFermi).lt.FermiPrec)) then
            exit
        elseif(outerIter.eq.PairingIter .and. .not. PairingIter.eq.1) then
            print 1, PairingIter, OldFermi, Fermi
        endif        
    enddo 
    !---------------------------------------------------------------------------
    call SetDeltas
    ! Compute pairingdensity for the calculation of the next Lipkin parameter
    if(Lipkin) then
      !call CompPDensity(PairDensity)
      !LNLambda = CalcLambda2()
    endif
  end subroutine SolvePairing
 
  subroutine SetDeltas
  !-----------------------------------------------------------------------------
  ! Routine that assigns all the deltas to the wavefunctions for printing
  ! purposes.
  !-----------------------------------------------------------------------------
  use HFB, only : blockindices, blocksizes
  
  integer :: i, partner(1), P, it,jj,jjj, check(1), ii, iii
  
  select case (PairingType)
    case(0)
        !HF calculation, no Deltas
    case(1) 
        !BCS calculation
        do i=1,nwt
            call HFBasis(i)%SetDelta(Delta(i,1,1,1))
        enddo        
    case(2) 
        !HFB calculation
        do it=1,Iindex
          do P=1,Pindex
            do i=1,blocksizes(P,it)
              ii   = blockindices(i,P,it)
              iii  = mod(ii-1,nwt)+1

              partner = maxloc(abs(DBLE(Delta(:,i,P,it))))            
              jj  = blockindices(Partner(1),P,it)
              jjj = mod(jj-1,nwt)+1
              call HFBasis(iii)%SetDelta(Delta(partner(1),i,P,it))
              call HFBasis(iii)%SetPairPartner(jjj)
            enddo
          enddo
        enddo
  end select
  
  end subroutine SetDeltas

  pure logical function ConverFermi(Prec) result(Converged)
  !-----------------------------------------------------------------------------
  ! Subroutine that checks the Fermi level for convergence across mean-field 
  ! iterations. 
  ! It checks whether
  !  abs(Fermi(:,Iteration) - Fermi(:,Iteration-n))< FermiPrec
  ! for n =1,7
  !-----------------------------------------------------------------------------
  integer                   :: iter, it
  real(KIND=dp), intent(in) :: Prec
  Converged = .true.
  do iter=1,7
    do it=1,2
      if(abs(Fermi(it) - FermiHistory(it,iter)).gt.Prec) then
        Converged=.false.
        return
      endif
    enddo
  enddo
  
  end function ConverFermi

end module Pairing
