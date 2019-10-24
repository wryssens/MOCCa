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
  use HartreeFock
  use GradientHFB

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
  !-----------------------------------------------------------------------------
  ! Universal properties
  !-----------------------------------------------------------------------------
  ! Dispersion of the particle numbers
  ! < N^2 > - < N >^2
  real(KIND=dp) :: PairingDisp(2) = 0.0_dp
  ! Energy associated with pairing obtained summing up weithed matrix elements
  ! of Delta_kl
  real(KIND=dp) :: PairingEnergy(2) = 0.0_dp
  ! Fermi levels and their history
  real(KIND=dp) :: Fermi(2) = -10.0_dp, FermiHistory(2,7) = 0.0_dp
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
  !matrix of overlaps < psi | T phi > in the HF basis used to identify
  ! "partner" states
  !-----------------------------------------------------------------------------
  real(KIND=dp), allocatable :: overlapT2(:,:,:,:)
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
  ! LNFix determines whether or not LNLambda is held fixed.
  logical                     :: Lipkin =.false.
  real(KIND=dp)               :: LNLambda(2)=0.0_dp, LNFix(2)=0.0_dp
  real(KIND=dp)               :: DN2(2) = 0.0_dp
  logical                     :: ConstrainDispersion=.false.
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
  !-----------------------------------------------------------------------------
  ! If MOCCa should look for a particular HF configuration
  logical       :: HFConfig=.false.

  !-----------------------------------------------------------------------------
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
     real(KIND=dp) :: CutProton     =-100.0_dp, CutNeutron   =-100.0_dp
     real(KIND=dp) :: AlphaProton   =-100.0_dp, AlphaNeutron =-100.0_dp
     real(KIND=dp) :: GammaProton   =-100.0_dp, GammaNeutron =-100.0_dp
     real(KIND=dp) :: Protongap     =   0.0_dp, NeutronGap   =   0.0_dp
     real(KIND=dp) :: LNFixN        =   0.0_dp, LNFixP       =   0.0_dp
     real(KIND=dp) :: DN2P          =   0.0_dp, DN2N         =   0.0_dp
     real(KIND=dp) :: PairingGradientNeutron       = -100.0_dp
     real(KIND=dp) :: PairingGradientProton        = -100.0_dp
     real(KIND=dp) :: PairingStabCutNeutron        = -100.0_dp
     real(KIND=dp) :: PairingStabCutProton         = -100.0_dp
     real(KIND=dp) :: PairingStabFacNeutron        = -100.0_dp
     real(KIND=dp) :: PairingStabFacProton         = -100.0_dp
     real(KIND=dp) :: HFBGaugeNeutron = 0.0_dp, HFBGaugeProton = 0.0_dp
     integer       :: i
     logical       :: SemiBCS=.false., SemiBCSNeutron=.false.
     logical       :: SemiBCSProton=.false.
     character(len=3) :: UpperType

     NameList /Pairing/ PairingNeutron, PairingProton, CutProton , RhoSat,     &
     &                  AlphaProton, AlphaNeutron,                             &
     &                  GammaProton, GammaNeutron, PairingContributionUpot,    &
     &                  PairingGradientNeutron, PairingGradientProton,         &
     &                  PairingStabCutNeutron, PairingStabCutProton ,          &
     &                  PairingStabFacNeutron, PairingStabFacProton ,          &
     &                  Type, CutNeutron, FreezeOccupation, PairingIter,       &
     &                  HFBMix, NeutronGap, ProtonGap, ConstantGap,            &
     &                  Lipkin, LNFraction, GuessKappa, PairingMu, CutType,    &
     &                  SemiBCS, SemiBCSNeutron, SemiBCSProton,                &
     &                  QPinHFBasis, SolvePairingStart,QPPrintWindow, Block,   &
     &                  FermiSolver, HFBIter,                                  &
     &                  HFBGaugeNeutron, HFBGaugeProton,                       &
     &                  HFConfig, LNFixN, LNFixP,                              &
     &                  DN2P, DN2N, ConstrainDispersion, HFBlock, HFBreduce,   &
     &                  Blockconsistent, aliyangle, fermimomentum, PfSolver,   &
     &                  AllowFuzzyNumber

     read(unit=*, NML=Pairing)

     PairingStrength(1) = PairingNeutron;  PairingStrength(2) = PairingProton
     PairingGradientStrength(1) = PairingGradientNeutron
     PairingGradientStrength(2) = PairingGradientProton
     PairingStabCut(1)  = PairingStabCutNeutron
     PairingStabCut(2)  = PairingStabCutProton
     StabilisingGapFactor(1) = PairingStabFacNeutron
     StabilisingGapFactor(2) = PairingStabFacProton
     PairingCut(1)      = CutNeutron    ;  PairingCut(2)      = CutProton
     Alpha(1)           = AlphaNeutron  ;  Alpha(2)           = AlphaProton
     DDExp(1)           = GammaNeutron  ;  DDExp(2)           = GammaProton
     Gaps(1)            = NeutronGap    ;  Gaps(2)            = ProtonGap
     LNFIx(1)           = LNFIXN        ;  LNFIX(2)           = LNFIXP
     DN2(1)             = DN2N          ;  DN2(2)             = DN2P
     HFBGauge(1)      = HFBGaugeNeutron ;  HFBGauge(2)      = HFBGaugeProton
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
     if(DDExp(2).eq.-100.0_dp) then
      if(DDExp(1).ne.-100.0_dp) then
        DDexp(2) = DDExp(1)
      else
        DDexp = 1.0_dp
      endif
     endif
     if(DDExp(1).eq.-100.0_dp) then
      if(DDExp(2).ne.-100.0_dp) then
        DDExp(1) = DDExp(2)
      endif
     endif
     !-----------------------------------------------------------------------
     ! use a pairing gradient term a la Fayans?
     !-----------------------------------------------------------------------
     if(PairingGradientStrength(2).eq.-100.0_dp) then
      if(PairingGradientStrength(1).ne.-100.0_dp) then
        PairingGradientStrength(2) = PairingGradientStrength(1)
      else
        PairingGradientStrength = 0.0_dp
      endif
     endif
     if(PairingGradientStrength(1).eq.-100.0_dp) then
      if(PairingGradientStrength(2).ne.-100.0_dp) then
        PairingGradientStrength(1) = PairingGradientStrength(2)
      endif
     endif
     !-----------------------------------------------------------------------
     ! Use Erler's stabilised functional?
     !-----------------------------------------------------------------------
     if(PairingStabCut(2).eq.-100.0_dp) then
      if(PairingStabCut(1).ne.-100.0_dp) then
        PairingStabCut(2) = PairingStabCut(1)
      else
        PairingStabCut = 0.0_dp
      endif
     endif
     if(PairingStabCut(1).eq.-100.0_dp) then
      if(PairingStabCut(2).ne.-100.0_dp) then
        PairingStabCut(1) = PairingStabCut(2)
      endif
     endif
     !-----------------------------------------------------------------------
     ! Use initial factor for Erler's stabilised functional?
     ! Setting this to specific values sometimes helps to stabilise the 
     ! initial phase of a code run.
     !-----------------------------------------------------------------------
     if(StabilisingGapFactor(2).eq.-100.0_dp) then
      if(StabilisingGapFactor(1).ne.-100.0_dp) then
        StabilisingGapFactor(2) = StabilisingGapFactor(1)
      else
        StabilisingGapFactor = 0.0_dp
      endif
     endif
     if(StabilisingGapFactor(1).eq.-100.0_dp) then
      if(StabilisingGapFactor(2).ne.-100.0_dp) then
        StabilisingGapFactor(1) = StabilisingGapFactor(2)
      endif
     endif

     !-----------------------------------------------------------------------
     ! Sanity checks
     !-----------------------------------------------------------------------
     if ( PairingType .ne. 2 .and. any( abs(PairingStabCut) .ne. 0.0 ) ) then
       print '(/," Stabilised EDF only implemented for HFB. Stopping.",/)'
       stop
     endif
     if ( Lipkin .and. any( abs(PairingStabCut) .ne. 0.0 ) ) then
       print '(/," LN and stabilised EDF are incompatible. Stopping.",/)'
       stop
     endif

     !-----------------------------------------------------------------------
     ! set flag when pair EDF is of old-school ULB form. The code then 
     ! bypasses some parts designed for the more complete form.
     !-----------------------------------------------------------------------
     if (  PairingGradientStrength(1) .eq. 0.0_dp .and. &
         & PairingGradientStrength(2) .eq. 0.0_dp .and. &
         & PairingStabCut(1)          .eq. 0.0_dp .and. &
         & PairingStabCut(2)          .eq. 0.0_dp .and. &
         & DDExp(1)                   .eq. 1.0_dp .and. &
         & DDExp(2)                   .eq. 1.0_dp ) then
       PairingULB = .true.
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
     if(any(DN2.lt.0.0_dp)) then
       call stp('Dispersion needs to be larger than 0.')
     endif
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
        !See if the user asked for a specific configuration
        if(HFConfig) then
          call ReadHFConfig()
          HFFill => PickHFConfig
        else
          HFFill => NaiveFill
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
        GetGaps           => HFBGaps

        !Decide on the diagonalization routine
        if(SC) then
          DiagonaliseHFBHamiltonian => DiagonaliseHFBHamiltonian_Signature
          LNCr8                     => LNCr8_sig
        else
          DiagonaliseHFBHamiltonian => DiagonaliseHFBHamiltonian_NoSignature
          lncr8                     => LNCr8_nosig
        endif

        ! Decide on the Fermi solver. In addition, if a default choice for
        ! HFBIter was taken, adapt it to the solver.
        if(trim(FermiSolver).eq.'BROYDEN') then
          FindFermiEnergy   => HFBFindFermiEnergyBroyden
          if(HFBIter.eq.-1) HFBIter=50
        elseif(trim(FermiSolver).eq.'GRADIENT') then
          FindFermiEnergy   => HFBFermiGradient
          if(HFBIter.eq.-1) HFBIter=500
        elseif(trim(FermiSolver).eq.'BISECTION') then
          if (Lipkin) then 
            FindFermiEnergy   => HFBLNFindFermiEnergyBisection
          else
            FindFermiEnergy   => HFBFindFermiEnergyBisection
          endif
          if(HFBIter.eq.-1) HFBIter=50
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
        allocate(overlapT2(HFBSize,HFBSize,2,2));

        !-----------------------------------------------------------------------
        ! When Lipkin-Nogami is present, make space for the modified pairing ga
        ! MB: as far as I can see, PairDensity, Chi, ChiStar are never
        ! calculated or used ...
        if(Lipkin) then
          allocate(DeltaLN(HFBSize,HFBSize,2,2))   ; DeltaLN = 0.0_dp
          allocate(PairingFieldLN(nx,ny,nz,2))     ; PairingFieldLN = 0.0_dp
          if (.not.allocated(PairDensity)) allocate(PairDensity(nx,ny,nz,2))
          allocate(Chi(nx,ny,nz,2), ChiStar(nx,ny,nz,2))
        endif
        !Note that we will store the complete U & V matrices
        allocate(U(HFBSize,2*HFBSize,2,2)) ; allocate(V(HFBSize,2*HFBSize,2,2))
        U = 0.0_dp ; V = 0.0_dp
        allocate(HFBHamil(2*HFBSize,2*HFBSize,2,2)) ; HFBHamil=0.0_dp
        !-----------------------------------------------------------------------
        ! the contribution of the pairing EDF to the single-particle-potential
        ! is asked for, in which case we need to sum up the pairi densities
        if ( PairingType .eq. 2) then
          if (.not.allocated(PairDensity)) allocate(PairDensity(nx,ny,nz,2))
        endif
        !-----------------------------------------------------------------------
        ! By default, consider even-even nuclei.
        ! The routine readblocking will change the numberparity if necessary.
        HFBNumberParity=1
        !-----------------------------------------------------------------------
        !Reset the canonical wavefunctions
        allocate(CanBasis(nwt)) ; allocate(CanDeriv(nwt))
        CanDeriv = 0
        do i=1,nwt
            call CanBasis(i)%Resetwf()
        enddo
        !-----------------------------------------------------------------------
        !Setting numerical default parameters
        if(PairingIter.eq.1) then
          PairingIter=1
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

    1  format (22('-'), 'Pairing Parameters', 22('-'))
    2  format (A26)
    3  format (45x, ' N ',8x, ' P ')
    4  format ('  Pairing Strength (MeV fm^3)          ', 1x,f10.3,1x,f10.3)
    5  format ('  Pairing Cutoff   (Mev)               ', 1x,f10.3,1x,f10.3)
    6  format ('  Pairing Cut width(MeV)               ', 1x,f10.3,1x,f10.3)
    7  format ('  Alpha                                ', 1x,f10.3,1x,f10.3)
    77 format ('  Gamma                                ', 1x,f10.3,1x,f10.3)
    78 format ('  Pairing Gradient Strength (MeV fm^5) ', 1x,f10.3,1x,f10.3)
    79 format ('  Pairing Stabilising Cut (MeV)        ', 1x,f10.3,1x,f10.3)
    80 format ('  Initial Stabilising Factor set to    ', 1x,f10.3,1x,f10.3)
    81 format ('  add contribution from pairing EDF to Upot')
    8  format ('  Attention: configuration is frozen.')
    9  format ('  HFBMixing                            ', 1x, f10.3)
   10  format ('  BCS-like: only Delta_{i,ibar} not zero.' )
   100 format ('  BCS-like Protons : only Delta_{i,ibar} not zero.' )
   101 format ('  BCS-like Neutrons: only Delta_{i,ibar} not zero.' )
   11  format ('  Constant Gap activated. ', /,                                 &
   &          '   Neutron Gap: ' , f8.3   , /,                                  &
   &          '   Proton  Gap: ' , f8.3   )
   12  format ('  Lipkin-Nogami prescription active.' )
   121 format ('  Constant \lambda_2                   ', 2x,f9.3,2x,f9.3)
   122 format ('  Dispersion constrained               ', 2x,f9.3,2x,f9.3)
   13  format ('  Cosine cutoff used. ')
   14  format ('  Symmetric Fermi cutoff used.')
   15  format ('  Fermisolver iterations (external):   ', i5)
   16  format ('  Fermisolver iterations (internal):   ', i5)
   17  format ('  LNFraction on the HFB hamiltonian:   ', f9.3)
   18  format ('  Gauge of the HFB hamiltonian:        ', 2x,f9.3,2x,f9.3)
   19  format ('  Fermisolver used:                    ', a9)
   20  format ('  HFBreduce is ON.')

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

    print  3
    print  4, PairingStrength
    print  5, PairingCut
    print  6, PairingMu
    print  7, Alpha
    print 77, DDExp
    print 78, PairingGradientStrength
    if ( any(PairingStabCut .ne. 0.0 ) ) print 79, PairingStabCut
    if (PairingContributionUpot) print 81

    if(PairingType.eq.2) print 9, HFBMix
    if(Lipkin)           print 12
    if(any(LNFIX.ne.0.0)) print 121, LNFIX
    if(ConstrainDispersion) print 122, DN2
    print *
    if(PairingType .eq. 2) print 19, FermiSolver
    if(HFBReduce) print 20
    print 15 , PairingIter
    if(PairingType .eq. 2) print 16, HFBIter
    if(Lipkin)             print 17, LNFraction
    if(PairingType .eq. 2) print 18, HFBGauge


    print *

    if(allocated(QPExcitations)) call PrintBlocking

  end subroutine PrintPairingInfo

  subroutine PrintPairing
    !---------------------------------------------------------------------------
    ! Subroutine that prints intermediate info on the pairing.
    !
    !---------------------------------------------------------------------------

    use Densities, only : Density

    1 format (/, 1x,35('-'), ' Pairing ', 36('-'))
    2 format (25x, ' N ',7x, ' P ')
    3 format (' Fermi Level (MeV) ',2x,f10.5,2x,f10.5)
    4 format (' Particles         ',2x,f10.5,2x,f10.5)
   51 format (' Constrained ')
    5 format (' Dispersion        ',2x,f10.5,2x,f10.5)
    6 format (' Lambda_2          ',2x,f10.5,2x,f10.5)
   61 format (' Lambda_2 (Constr.)',2x,f10.5,2x,f10.5 )
    7 format (1x,75('-'),/)
    print 1

    select case(PairingType)
      case (0)
        !-------------------------------------------------------------------
        !Print Configuration counting
        call PrintHFConfiguration
        return
      case (1,2)
        !-------------------------------------------------------------------
        print 2
        print 3, Fermi
        print 4, sum(Density%Rho(:,:,:,1))*dv, sum(Density%Rho(:,:,:,2))*dv

        if(ConstrainDispersion) then
          print 51
        endif
        print 5, PairingDisp
        if(Lipkin) then
          print 6, LNLambda
        else if(any(LNFIX.ne.0.0_dp)) then
          print 61, LNLambda
        endif
        !--------------------------------------------------------------------
    end select
    !------------------------------------------------------------------------
    ! Print the number parity of the HFB vacuum
    if(PairingType.eq.2) call PrintNumberParities

    if(PairingType.eq.2) call PrintHFBconvergence
    !------------------------------------------------------------------------
    print 7
    !------------------------------------------------------------------------
    select case(PairingType)
      case (2)
        call CompCutPairDensity
        call CompPairingEDF      (PairDensity)
        call CompAverageGaps     (PairDensity)
        call PrintPairingProfiles(PairDensity)
        call PrintPairingEnergies
    end select
    print 7

  end subroutine PrintPairing

  subroutine PrintHFConfiguration
    !---------------------------------------------------------------------------
    ! Simple subroutine to count in what kind of Hartree-Fock configuration
    ! we are currently sitting.
    !---------------------------------------------------------------------------

    1 format (' Hartree-Fock Configuration          ')
    2 format (' ------------------------------------')
    3 format ('  P, Rz      + +   + -   - +   - - ')
    4 format ('Neutron | ', 4i6)
    5 format ('Proton  | ', 4i6)

    integer :: N(2,2,2), it, P, S,i

    print 1
    print 2
    print 3
    print 2

    ! Count
    N = 0
    do i=1,nwt
        if(HFBasis(i)%GetOcc() .lt. 0.5_dp) cycle
        it = (HFBasis(i)%GetIsospin()  + 3)/2
        P  = (HFBasis(i)%GetParity()   + 3)/2
        S  = (HFBasis(i)%GetSignature()+ 3)/2

        N(P,S,it) = N(P,S,it) + 1
        if(TRC) N(P,mod(S,2)+1,it) = N(P,mod(S,2)+1,it) + 1
    enddo

    print 4, N(2,2,1), N(2,1,1),N(1,2,1), N(1,1,1)
    print 5, N(2,2,2), N(2,1,2),N(1,2,2), N(1,1,2)
    print 2

  end subroutine PrintHFConfiguration

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
    
        
    use HFB, only: blocksizes
    
    1 format('---------------------------------------',/,   &
    &        ' Warning! ',/,                                &
    &        ' Pairing calculations did not converge.',/,   &
    &        ' Pairing iterations:  ', i4,/,                &
    &        ' Old Fermi:    ', 2f12.7,/,                   &
    &        ' New Fermi:    ', 2f12.7,/,                   &
    &        '---------------------------------------')

    real(KIND=dp)    :: OldFermi(2), Particles(2)
    integer          :: outeriter,j


    if(PairingType.eq.0) then
     if(FreezeOccupation) return
     call HFFill(HFConfiguration)
     !if(HFConfig) FreezeOccupation=.true.
     return
    endif

    Particles(1) = Neutrons ; Particles(2) = Protons
    !Initialise the HFB indices properly
    if(.not.allocated(blockindices) .and. PairingType.eq.2) then
      call PrepareHFBModule
      !--------------------------------------------------------------------------
      ! Initialize Stabilising factor for Erler's stabilised pairing EDF when
      ! needed. During the first call in a program run, the pairing energy is not
      ! yet known, such that the StabilisingFactor will be set to a small value
      ! in the subroutine.
      if ( any(PairingStabCut .ne. 0.0 ) ) then
        call CompStabilisingFactor
      endif
    endif
    !---------------------------------------------------------------------------
    ! Reinitialise Kappa when asked for
    if(GuessKappa .and. PairingType.eq.2) then
      call GuessHFBMatrices(PairingType,Fermi)
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
    if(.not.Lipkin) then
      if(any(LNFix.ne.0.0_dp)) then
        LNLambda = LNFix
      else if(.not. ConstrainDispersion) then
        LNLambda= 0.0_dp
      endif
    endif

    call CompDensityFactor

    !---------------------------------------------------------------------------
    ! calculate matrix of overlaps between psi and T psi' in HF the basis
    if(PairingType.eq.2) then
      call HFBoverlapT2(overlapT2)
      call SetDeltas
    endif
    
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
        call FindFermiEnergy(Fermi,LNLambda,Delta,DeltaLN,Lipkin,DN2,          &
        &                                         ConstrainDispersion,FermiPrec)
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

    !---------------------------------------------------------------------------
    ! the pairing interaction asks for the contribution to the single-particle
    ! Hamiltonian. Requires to calculate pair densities (which is avoided for
    ! the rest of the treatment of pairing)
    if ( PairingType .eq. 2) then
      call CompCutPairDensity
      if ( PairingContributionUpot ) then
        call CompUpotPairingContribution(PairDensity)
      endif
    endif

    !---------------------------------------------------------------------------
    ! calculate Stabilising factor for Erler's stabilised pairing EDF to be
    ! used in the next iteration.
    ! NOTE: if there has not been pairing already, the stabilisation will not
    ! restart it on its own.
    if ( PairingType.eq.2 .and. any(PairingStabCut .ne. 0.0 ) ) then
       call CompPairingEDF(PairDensity)
       call CompStabilisingFactor
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
              ! scan overlap matrix |< psi | T phi >|^2 (if available), 
              ! which also works when pairing breaks down.
              ! Note that maxloc returns an array (here of rank 1), hence 
              ! partner(1) has to be an array too.
              ! Note that overlapT2 is the absolute square of the (complex) 
              ! overlap; hence it's positive
              if (allocated(overlapT2)) then
                partner = maxloc(overlapT2(:,i,P,it))
              else
                partner = maxloc(abs(DBLE(Delta(:,i,P,it))))
              endif
              jj  = blockindices(Partner(1),P,it)
              jjj = mod(jj-1,nwt)+1
              call HFBasis(iii)%SetDelta(Delta(partner(1),i,P,it))
              call HFBasis(iii)%SetPairPartner(jjj)
              ! print'(" SetDeltas ",10i4)',it,P,i,ii,iii,partner(1),jj,jjj
            enddo
          enddo
        enddo
  end select

  end subroutine SetDeltas

  subroutine CompCutPairDensity
    !---------------------------------------------------------------------------
    ! this subroutine is somehow doubling part of what is done in HFB through
    ! summing matrix elements. At time being, it is only called when the 
    ! contribution of the pairing EDF to the single-particle potential is needed
    ! Not entirely sure if the overall sign is correct. As all true observables
    ! are bilinear in the pair density, this is not a real problem, but the
    ! average gaps are dependent on its sign. Compared to the HFB module, there
    ! is an additional factor -1 below such that the trace of rho~ is positive.
    !---------------------------------------------------------------------------
    complex(KIND=dp)  :: ActionOfPairing(nx,ny,nz)
    real(KIND=dp)     :: factor
    integer           :: i, it, j, ii, sig1, sig2, jj, k, P, iii, jjj, l
    real(KIND=dp)     :: Cutoff(2)
    type(Spinor)      :: Temp(2)
    logical           :: TR

    if ( .not.allocated(PairDensity)) call stp ('CompCutPairDensity !')

    PairDensity = (0.0_dp,0.0_dp)

    ! using the skew-symmetry of the contributions, the sum is just over
    ! half the combinations
    do it=1,Iindex
    do P=1,Pindex
      do j=1,blocksizes(P,it)
        jj        = Blockindices(j,P,it)
        jjj       = mod(jj-1,nwt)+1
        sig1      = HFBasis(jjj)%GetSignature()
        Cutoff(1) = PCutoffs(jjj)
        if(Cutoff(1).lt.HFBNumCut) cycle
        if(TRC .and. jj .ne. jjj) then
          sig1 = - sig1
        endif
        do i=j+1,blocksizes(P,it)
          ii   = Blockindices(i,P,it)
          iii  = mod(ii-1,nwt)+1
          sig2 = HFBasis(iii)%GetSignature()
          if( TRC .and. ii .ne. iii) then
            sig2 = - sig2
          endif
          !Selection on signature quantum number; the loop structure already
          !selects on both parity and isospin.
          if(sig1 .ne. -sig2 ) cycle
          Cutoff(2) = PCutoffs(iii)

          ! Save some CPU cycles
          if(Cutoff(1)*Cutoff(2)*abs(KappaHFB(i,j,P,it)) .lt. HFBNumCut) cycle

          if( TRC .and. ( ii.ne. iii  .or. jj.ne.jjj ) .and. & 
              &  .not. (ii .ne.  iii .and. jj .ne. jjj )) then
            TR = .true.
          else
            ! Should only happen when Time-reversal is conserved, but signature
            ! isn't.
            TR = .false.
          endif

          ! Note that this does automatically include a Time-reversal operator
          ! when appropriate
          ! the factor 2 compensates for summing just over half of the 
          ! (symmetric) combinations of i and j
          ! the "-" sign assures that tr rho~ is usually positive.
          ActionOfPairing = GetPairDensity(HFBasis(jjj)%Value,HFBasis(iii)%Value,TR)
          factor = -2.0_dp * Cutoff(1)*Cutoff(2)*DBLE(KappaHFB(i,j,P,it))
          PairDensity(:,:,:,it) = PairDensity(:,:,:,it) + factor*ActionOfPairing(:,:,:)
        enddo
      enddo
    enddo
    enddo
  ! print '(" norm pair density ",4f12.6)',dv*sum(PairDensity(:,:,:,1)), &
  !                                 &      dv*sum(PairDensity(:,:,:,2))
  end subroutine CompCutPairDensity

  function GetUpotPairingContribution() result(UpotPC)
  !-----------------------------------------------------------------------------
  ! provide contribution from pairing EDF to single-particle potential Upot
  !-----------------------------------------------------------------------------
    real(KIND=dp) :: UpotPC(nx,ny,nz,2)
 
    UpotPC = 0.0_dp

    ! nothing to do - what am I doing here ...?
    if ( .not.PairingContributionUpot ) return

    ! there is no PairingContribution in a HF calculation ...
    if ( PairingType .ne. 2 ) return

    if ( .not.allocated(UpotPairingContribution) ) call stp('GetUpotPairingContribution!')

    UpotPC = UpotPairingContribution

  end function GetUpotPairingContribution

  function GetUpotPairingRearrangement() result(UpotPR)
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
     real(KIND=dp) :: UpotPR

    ! nothing to do - what am I doing here ...?
    if ( .not.PairingContributionUpot ) then
      UpotPR = 0.0_dp
    else
      UpotPR = UpotPairingRearrangementEnergy
    endif

  end function GetUpotPairingRearrangement

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

  if (pairingtype.eq.0) return !Don't do this check for HF calculations
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
