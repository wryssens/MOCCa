module HFB
!-------------------------------------------------------------------------------
! Module that contains everything needed to solve the HFB equations.
!-------------------------------------------------------------------------------
! Some general notes on technicalities:
! *) When conserving time-reversal not all spwfs are represented in memory.
!    the wavefunctions are actually stored in memory. However, for easier coding
!    we still construct explicitly the entire HFB matrices, even if strictly not
!    needed. This will save us severe trouble in 'index juggling'.
!    (Mostly because this symmetry is very easy to understand and enforce when
!    signature conservation is present, but very hard when it is not)
!    An important thing to understand is the labelling in this case:
!        ii = mod(i-1,nwt) + 1
!   which means that ii = i if i <= nwt and ii = i - nwt if i > nwt. This is
!   (I hope) more efficient than if-clauses in deeply nested loops.
!
! *) This module uses a lot of complex numbers. It is important to remember the
!    proper expressions to convert numbers from complex form to real and
!    and imaginary parts and vice versa while maintaining double precision.
!    (Naive use of cmplx, real and aimag will convert to single precision!)
!
!    R = DBLE ( C )
!    I = DIMAG( C )
!    C = DCMPLX( R, I)
!
!    Note that different compilers sometimes offer other intrinsics,
!    like gfortran offers `RealPart' and 'ImagPart', but this is not portable
!    at all.
!------------------------------------------------------------------------------

  use CompilationInfo
  use Geninfo
  use PairingInteraction
  use SpwfStorage
  use Spinors

  implicit none

  !-----------------------------------------------------------------------------
  ! Real that saves the number of particles. Useful because it can be iterated.
  real(KIND=dp), save :: Particles(2)
  !-----------------------------------------------------------------------------
  ! HFBSize controls the dimensions of all matrices related to HFB. It is
  ! just nwt when time-reversal is broken, but 2*nwt when it is conserved.
  integer :: HFBSize=0
  !-----------------------------------------------------------------------------
  ! Numerical cutoff to determine what components in the canonical
  ! transformation to ignore.
  real(KIND=dp),parameter :: HFBNumCut=1d-12
  !-----------------------------------------------------------------------------
  ! HFB density and anomalous density matrix note that these are complex
  ! in general. The old matrices are the matrices at previous iterations
  ! for smoothing purposes.
  ! Third and fourth index are for parity & isospin blocks.
  complex(KIND=dp), allocatable  ::    RhoHFB(:,:,:,:),    KappaHFB(:,:,:,:)
  complex(KIND=dp), allocatable  :: OldRhoHFB(:,:,:,:), OldKappaHFB(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! U & V matrices that diagonalise the HFB Hamiltonian
  ! CRITICALLY IMPORT NOTE:
  ! These store ALL of the eigenvectors of the HFBhamiltonian, so INCLUDING
  ! occupied and unoccupied levels. In this sense the U and V in this code
  ! DO NOT immediately correspond to U and V in most texts. In fact the U and V
  ! in texts are one half of these matrices.
  ! Third and fourth index are for parity & isospin blocks.
  complex(KIND=dp), allocatable :: U(:,:,:,:), V(:,:,:,:)
  complex(KIND=dp), allocatable :: OldU(:,:,:,:), OldV(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! HFB Hamiltonian
  ! Matrix of the form
  !      (   h      Delta  )
  !  H = (                 )
  !      (-Delta^*    -h^* )
  ! Third and fourth index are for parity & isospin blocks.
  complex(KIND=dp), allocatable :: HFBHamil(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! This integer stores all of the columns we want to take from the U and V
  ! to construct Density and anomalous density matrix.
   integer, allocatable :: HFBColumns(:,:,:), oldColumns(:,:,:)
   integer, allocatable :: NonBlockedHFBColumns(:,:,:)
  !-----------------------------------------------------------------------------
  ! Logical. Check all kinds of relations that should hold for correct HFB
  ! calculations. Only to be used for debugging purposes.
  logical,parameter :: HFBCheck=.true.
  !-----------------------------------------------------------------------------
  ! Real that determines how much to damp the RhoHFB and KappaHFB matrices.
  real(KIND=dp) :: HFBMix=0.1_dp
  !-----------------------------------------------------------------------------
  ! Switch to simulate BCS with HFB. Mostly useful for debugging purposes.
  ! It forces the code to have only diagonal elements in Delta. (Diagonal,
  ! in the sense between functions and their time-inverse.) Should not be
  ! activated when Time-Reversal is not conserved!
  ! ( I don't understand how this is correctly done in CR8. Sure, one can invent
  !   some sense of `diagonal', but how this is physically justified is beyond
  !   me.)
  logical :: BCSinHFB(2)=.false.
  !-----------------------------------------------------------------------------
  ! Fraction of the Lipkin-Nogami term to put into the pairing part of the
  ! HFB Hamiltonian, versus the part to put in the mean-field part.
  ! See the discussion in the PhD Thesis of Gall.
  real(KIND=dp) :: LNFraction=1.0_dp
  real(KIND=dp) :: HFBGauge(2) = 0.0_dp
  !----------------------------------------------------------------------------
  !Quasiparticle Energies: note that all of these get stored, even the occupied
  ! ones.
  ! Also the signatures of the quasiparticles are stored, since they are quite
  ! a bore to get everytime they are needed.
  real(KIND=dp), allocatable :: QuasiEnergies(:,:,:), QuasiSignatures(:,:,:)
  real(KIND=dp), allocatable :: Quasisimplex(:,:,:)
  !-----------------------------------------------------------------------------
  ! Maximum number of iterations that the Fermi solver can take
  ! Default is -1 so that MOCCa can detect if it has been read.
  integer :: HFBIter = -1
  !-----------------------------------------------------------------------------
  ! Logical that traces if quasiparticles excitations are defined in the HF or
  ! in the canonical basis.
  logical :: QPinHFBasis=.false.
  !-----------------------------------------------------------------------------
  ! Indices and other properties of the quasiparticle states that we want to 
  ! block:
  ! Positive QPExcitations(i) is the index of the single-particle statexs in 
  !   the HF basis that is supposed to be the dominant components of the blocked
  !   qp state i. For QPFilled(i) = 1 the U matrix is scanned for the largest
  !   overlap with this HF state, for PFilled(i) = 0 the V matrix is scanned for
  !   the largest overlap with this HF state instead. If the norm of the 
  ! Negative QPExcitations(i) indicates that the lowest qp state
  !   in some parity/isospin block is to be blocked.
  ! From these indices, the code determines some derived quantities that are
  ! useful at various places:
  ! QPParities(i) indicate the parity block in which the HF state with index
  !   QPExcitations(i) can be found
  ! QPIsospins(i) indicates the isospin block in which the HF state with index
  !   QPExcitations(i) can be found
  ! QPSignatures (seems to be obsolete?)
  ! QPBlockind(i) is the index of the HF state QPExcitations(i) in the 
  !   respective parity/isospin block it can be found in.
  ! QPFilled(i) is the information if the HF state that dominates this QP is to 
  !   be (almost) filled (QPFilled(i) = 1) or to be (almost) empty 
  !   (QPFilled(i) = 0). The default value is 1, which is what was implied in
  !   early code versions.
  !   Note: for decorrelated-pair 2qp states, fill = 0 and fill = 1 yield the
  !   same state.
  ! QPBlockPartnerInd(i) contains (if defined on input) the HF index of the 
  !   partnerstate in the HF basis that is to be blocked when the V matrix is 
  !   scanned. When not defined on input, the code scans for the HF state whose
  !   time-reversed has the largest overlap with the HF state with index
  !   QPExcitations(i) instead.
  ! QPisinDCP(i) is a flag that indicates that this QP is in a decorrelated 
  !   pair for which the selection of blocked qps is done with a separate 
  !   dedicated piece of code (which at time being seems to deliver exactly
  !   the same results as the standard procedure).
  ! QPBlockedBefore(i) QP indeces of previously blocked QPs. Used to avoid 
  ! blocking again an already blocked QP (might happen for decorrelated pair
  ! states)
  ! QPBlockLog(i)
  !-----------------------------------------------------------------------------
  integer,allocatable :: QPExcitations(:), QPblockind(:), QPFilled(:)
  integer,allocatable :: QPPartners(:), QPBlockPartnerInd(:)
  integer,allocatable :: QPParities(:),QPIsospins(:), QPSignatures(:)
  integer,allocatable :: QPBlockedBefore(:),QPisinDCP(:)
  character(len=1),allocatable :: QPBlockLog(:)
  !-----------------------------------------------------------------------------
  ! Logical that switches on an experimental handling of unusual situations
  ! that may be encountered by the bisection solver for the fermi energy in
  ! the HFB (w/o LN) case.
  logical :: HandleBisectionEmergencies = .false.
  !-----------------------------------------------------------------------------
  ! Logical that allows the code to continueHFB calculations with the bisection
  ! solver slightly wrong particle number instead of forcing the Fermi energy 
  ! to change a lot.
  logical :: AllowFuzzyNumber=.false.
  !-----------------------------------------------------------------------------
  ! Angle for the alirotation in the y-direction in degrees. Use with care.
  real(KIND=dp), allocatable :: aliyangle(:)
  !-----------------------------------------------------------------------------
  ! When using the gradient solver, block before doing the gradient procedure
  ! or after.
  logical :: BlockConsistent=.true.
  !-----------------------------------------------------------------------------
  ! Integer that controls whether the program looks for even or odd HFB states.
  ! This needs to be stored for every parity-isospin block.
  integer :: HFBNumberParity(2,2) = 0
  !-----------------------------------------------------------------------------
  ! Pointers for the parity-isospin blocks, to facilitate index juggling.
  ! Indices are parity & isospin respectively. In addition, these pointers can
  ! nicely be used to simplify the storage case for time-reversal conserving
  ! calculations.
  ! The array blocksize carries the total number of wavefunctions in each
  ! block, while blockindices carries the index in HFBasis.
  ! Pindex and Iindex count the number of parity and isospin blocks.
  integer,allocatable :: blocksizes(:,:), blockindices(:,:,:)
  integer             :: Pindex, Iindex
  !----------------------------------------------------------------------------
  ! Canonical transformation and occupation numbers, both obtained by
  ! diagonalising the density matrix.
  complex(KIND=dp), allocatable   :: CanTransfo(:,:,:,:)
  real(KIND=dp), allocatable      :: Occupations(:,:,:)
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! 'Distance in energy' for the printing of quasiparticle excitations
  real(Kind=dp) :: QPPrintWindow=100.0_dp
  !-----------------------------------------------------------------------------
  ! Whether we want to use Broyden (fast and not 100% guaranteed) or Bisection
  ! method to solve for the Fermi energy. Note that in HFB without LN, 
  ! Bisection is handled with the fast, reliable and stable Brent's method,
  ! while HFB+LN uses a more basic and less efficient, but also stable
  ! bisection.
  character(len=12) :: FermiSolver='Broyden'
  !-----------------------------------------------------------------------------
  logical :: FermiNonConvergenceHistory(7) = .false.
  !-----------------------------------------------------------------------------
  ! Whether or not to enable momentum in the gradient solver.
  ! HIGHLY EXPERIMENTAL
  logical :: Fermimomentum = .false.
  !-----------------------------------------------------------------------------
  ! Size of the change in rho and kappa in the last iteration.
  real(KIND=dp) :: drho(2,2), dkappa(2,2)
  !-----------------------------------------------------------------------------
  ! Storing the pfaffian of the HFBHamiltonian in every subspace
  ! QPswapped keeps track of swapped qps in the blocks for diagnostics
  logical       :: PfSolver = .false.
  real(KIND=dp) :: pf(2,2) = 0, oldpf(2,2) = 0 , refpf(2,2) = 0
  integer       :: QPswapped(2,2,2) = 0
  integer       :: NP(2,2)    = 0
  integer       :: NPVac(2,2) = 0
  !-----------------------------------------------------------------------------
  ! Procedure pointer for the diagonalisation of the HFBhamiltonian.
  ! Either with or without signature conservation.
  abstract interface
    function LNCr8_interface(Delta, DeltaLN, flag) result(LNLambda)
      import :: dp
      integer, intent(inout) :: flag(2)
      complex(KIND=dp), intent(in), allocatable :: Delta(:,:,:,:)
      complex(KIND=dp), intent(in), allocatable :: DeltaLN(:,:,:,:)
      real(KIND=dp)                             :: LNLambda(2)
    end function
  end interface
  abstract interface
        subroutine diag_interface()
        end Subroutine
  end interface

  procedure(LNCR8_interface), pointer  :: LNCr8
  procedure(Diag_interface), pointer :: DiagonaliseHFBHamiltonian
  procedure(HFBNumberofParticles_ordinary), pointer :: HFBNumberofParticles
contains
  
  subroutine PrepareHFBModule
    !---------------------------------------------------------------------------
    ! This subroutine nicely prepares the module at the start of the program.
    !---------------------------------------------------------------------------
    integer :: i, P, it, ii

    Pindex = 1 ; if(PC) Pindex = 2
    Iindex = 1 ; if(IC) Iindex = 2

    if(.not.allocated(Blocksizes))   allocate( blocksizes(Pindex,Iindex))
    if(.not.allocated(Blockindices)) allocate( blockindices(HFBsize,Pindex,Iindex))

    ! Counts the number of wavefunctions in every parity-isospin block
    blocksizes=1 ; blockindices = 0
    do i=1,HFBSize
      ii=mod(i-1,nwt)+1
      P = (HFBasis(ii)%GetParity() +3)/2
      it= (HFBasis(ii)%GetIsospin()+3)/2
      blockindices(blocksizes(P,it),P,it) = i
      blocksizes(P,it) = blocksizes(P,it) + 1
    enddo
    Blocksizes = Blocksizes - 1
    !Quick sanity check
    if(sum(Blocksizes).ne.HFBSize) then
      print *, Blocksizes, nwt, sum(blocksizes)
      call stp('Counting error in PrepareHFBModule.')
    endif

    if(.not.allocated(CanTransfo)) then
      allocate(CanTransfo(HFBSize,HFBSize,Pindex,Iindex)) ; CanTransfo     = 0.0_dp
    endif
    if(.not.allocated(Occupations)) then
      allocate(Occupations(HFbsize,Pindex,Iindex))        ; Occupations    = 0.0_dp
    endif
    if(.not.allocated(QuasiEnergies)) then
      allocate(QuasiEnergies(2*HFBSize,Pindex,Iindex))    ; QuasiEnergies  = 0.0_dp
    endif
    if(.not.allocated(QuasiSignatures)) then
      allocate(QuasiSignatures(2*HFBSize,Pindex,Iindex))  ; QuasiSignatures= 0.0_dp
    endif
    if(.not.allocated(HFBColumns)) then
      allocate(HFBColumns(HFBSize,Pindex,Iindex))     ; HFBColumns     = 0
    endif
    
    HFBNumberofParticles => HFBNumberofParticles_ordinary 
    
  end subroutine PrepareHFBModule

  function HFBEnergy(Delta) result(Energy)
    !---------------------------------------------------------------------------
    ! Placeholder function for calculating the HFB energy.
    ! Note that this formula introduces an extra cutoff factor compared to R&S.
    !---------------------------------------------------------------------------
    ! Note that this only makes sense if the pairing EDF is strictly bilinear
    ! in kappa. The stabilized pairing EDF a la Erler et al is not, therefore
    ! one has to correct for double counting at the end.
    !---------------------------------------------------------------------------
    real(KIND=dp)             :: Energy(2)
    integer                   :: it,P,i,j
    complex(KIND=dp), allocatable,intent(in)  :: Delta(:,:,:,:)

    Energy = 0.0_dp
    do it=1,Iindex
      do P=1,Pindex
        do i=1,blocksizes(P,it)
          do j=1,blocksizes(P,it)
            Energy(it) = Energy(it) - DBLE(Delta(i,j,P,it)*KappaHFB(i,j,P,it))
          enddo
        enddo
      enddo
    enddo
    Energy = - 0.5_dp * Energy

    !---------------------------------------------------------------------------
    ! correct for different factor in gaps and energy when using the stabilised
    ! pairing EDF of Erler et al.
    !---------------------------------------------------------------------------
    ! Assuming a pairing functional that is linear in (kappa kappa*),
    ! the non-stabilised gaps and pair energy are related by
    !   Delta  = d E_pair / d kappa*
    !   E_pair = sum kappa* Delta 
    ! The stabilised quantities are
    !   Delta^s  = [ 1 + StabilisingGapFactor ] Delta 
    !   E_pair^s = [ 1 - StabilisingGapFactor ] E_pair
    !            = [ 1 - StabilisingGapFactor ] sum kappa* Delta 
    ! where StabilisingGapFactor = PairingStabCut^2 / E_pair^2 is a global
    ! state-independent factor. Therefore
    !              [ 1 - StabilisingGapFactor ]
    !   E_pair^s = ---------------------------- sum kappa* Delta^s
    !              [ 1 + StabilisingGapFactor ]
    !---------------------------------------------------------------------------
    do it=1,Iindex
      if ( abs(StabilisingGapFactor(it)) .ne. 0.0 ) then
         ! print '(" HFBEnergy: for ",i2," correct energy with ",3f16.8)',it, &
         !    & Energy(it), ( 1.0_dp - StabilisingGapFactor(it) ) / ( 1.0_dp + StabilisingGapFactor(it) ), &
         !    & Energy(it)*( 1.0_dp - StabilisingGapFactor(it) ) / ( 1.0_dp + StabilisingGapFactor(it) )
         Energy(it) =  Energy(it) * ( 1.0_dp - StabilisingGapFactor(it) ) &
               &                  / ( 1.0_dp + StabilisingGapFactor(it) )
      endif
    enddo
  end function HFBEnergy

  subroutine HFBPairingField(Field,FieldLN, Delta)
    !---------------------------------------------------------------------------
    ! Computes the pairing field due to a HFB model pairing.
    !  Delta(x,y,z) =
    !  Sum_{i,j}  0.5*fi*fj * kappa_{ij} * s * <r,s;r,-s|V|ij>
    ! When Lipkin-Nogami is activated, this routine also sums the modified
    ! pairing fields needed for it.
    ! They correspond to the following modified formula:
    !  Delta_LN(x,y,z) =
    !  Sum_{i,j} 0.5*fi*fj *Gamma_{ij}*s*<r,s;r,-s|V|ij>
    !  where Gamma = (1 - 2 Rho) * Kappa
    !---------------------------------------------------------------------------
    use Spinors
    use PairingInteraction

    complex(KIND=dp), allocatable, intent(inout) :: Field(:,:,:,:)
    complex(KIND=dp), allocatable, intent(inout) :: FieldLN(:,:,:,:)
    complex(KIND=dp), allocatable, intent(inout) :: Delta(:,:,:,:)
    complex(KIND=dp)                             :: ActionOfPairing(nx,ny,nz)
    complex(KIND=dp),allocatable                 :: Gamma(:,:,:,:)
    real(KIND=dp):: factor, factorLN, StabFac
    integer      :: i, it, j, ii, sig1, sig2, jj,k,P, iii, jjj,l
    real(KIND=dp):: Cutoff(2)
    type(Spinor) :: Temp(2)
    logical      :: TR

    !Make sure the indices re correctly initiated
    if(.not.allocated(blocksizes)) then
      call PrepareHFBModule
    endif

    Field  = 0.0_dp
    if(.not.allocated(Gamma)) allocate(Gamma(HFBSize, HFBSize,Pindex,Iindex))
    if(allocated(FieldLN)) then
      FieldLN = 0.0_dp ; Gamma = 0.0_dp
      !Compute gamma if needed
      ! Gamma =(1 - 2 * Rho)Kappa
      do it=1,Iindex
        do P=1,Pindex
          do j=1,blocksizes(P,it)
            do i=1,blocksizes(P,it)
              Gamma(i,j,P,it) = KappaHFB(i,j,P,it)
              do k=1,blocksizes(P,it)
                Gamma(i,j,P,it) = Gamma(i,j,P,it)-2*KappaHFB(k,j,P,it)*RhoHFB(i,k,P,it)
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
    !--------------------------------------------------------------------------
    !Note that the formula for the pairing field is symmetric in the change
    ! (i,j) => (j,i) since kappa and ActionOfPairing are both antisymmetric.
    ! Thus we only calculate half of the contributions and multiply by a hidden
    ! factor 2.
    ! (Note that the i=j contributions are not calculated since they are zero
    ! anyway.)
    if(allocated(FieldLN)) then
      do it=1,Iindex
        do P=1,Pindex
          do j=1,blocksizes(P,it)
            jj        = Blockindices(j,P,it)
            jjj       = mod(jj-1,nwt)+1
            sig1      = HFBasis(jjj)%GetSignature()
            !Temp(1)   = HFBasis(jjj)%GetValue()
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

              if(TRC .and. ( ii.ne. iii  .or. jj.ne.jjj ) .and. .not. (ii .ne. iii .and. jj .ne. jjj )) then
                TR = .true.
              else
                ! Should only happen when Time-reversal is conserved, but signature
                ! isn't.
                TR = .false.
              endif

              ActionOfPairing = GetPairDensity(HFBasis(jjj)%Value,HFBasis(iii)%Value,TR)
              factor    = Cutoff(1)*Cutoff(2)*DBLE(KappaHFB(i,j,P,it))
              factorLN  = Cutoff(1)*Cutoff(2)*DBLE(Gamma(i,j,P,it))
              do l=1,nx*ny*nz
                Field(l,1,1,it)   = Field(l,1,1,it)   - factor  *ActionOfPairing(l,1,1)
                FieldLN(l,1,1,it) = FieldLN(l,1,1,it) - factorLN*ActionOfPairing(l,1,1)
              enddo
            enddo
          enddo
        enddo
      enddo
      do it=1,Iindex
        do i=1,nx*ny*nz
          FieldLN(i,1,1,it) = FieldLN(i,1,1,it) * DensityFactor(i,1,1,it)
          Field(i,1,1,it)   = Field(i,1,1,it)   * DensityFactor(i,1,1,it)
        enddo
      enddo
    else
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

              if(TRC .and. ( ii.ne. iii  .or. jj.ne.jjj ) .and. .not. (ii .ne. iii .and. jj .ne. jjj )) then
                TR = .true.
              else
                ! Should only happen when Time-reversal is conserved, but signature
                ! isn't.
                TR = .false.
              endif

              ! Note that this does automatically include a Time-reversal operator
              ! when appropriate
              ActionOfPairing = GetPairDensity(HFBasis(jjj)%Value,HFBasis(iii)%Value,TR)
              factor          = Cutoff(1)*Cutoff(2)*DBLE(KappaHFB(i,j,P,it))
              do l=1,nx*ny*nz
                Field(l,1,1,it) = Field(l,1,1,it) - factor*ActionOfPairing(l,1,1)
              enddo
            enddo
          enddo
        enddo
      enddo
      
      !-------------------------------------------------------------------------
      ! Multiply with stabilising factor of Erler's  stabilised EDF "StabFac"
      ! (which is 1 for non-stabilized EDF's)
      !-------------------------------------------------------------------------
      do it=1,Iindex
        StabFac = 1.0_dp + StabilisingGapFactor(it)
        do i=1,nx*ny*nz
          Field(i,1,1,it) = Field(i,1,1,it) * DensityFactor(i,1,1,it) * StabFac
        enddo
      enddo
    endif

  end subroutine HFBPairingField

  subroutine HFBGaps(Delta,DeltaLN,PairingField,PairingFieldLN,Gaps,ConstantGap)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the Delta_{i j} for the HFB model.
    ! The formula (see Ring & Schuck, page 254 eq. (7.41))
    ! is the following:
    !
    ! Delta_{i, j} = \int dr <i,j | Delta(r) | r >
    !
    ! When Lipkin-Nogami is activated, this routine also sums the modified gaps,
    ! which are just analogous to the normal gaps:
    ! Delta_{LN,i,j} = \int dr <i,j | Delta(r)_{LN} | r >
    !---------------------------------------------------------------------------
    ! Gaps & Constantgap are not used in this subroutine, but are present to
    ! be compatible with the BCS GetGaps subroutine header.
    !---------------------------------------------------------------------------
    !
    ! NOTE THAT THIS ROUTINE IS NOT READY FOR TIME-SIMPLEX BREAKING AS DELTA
    ! CAN BE COMPLEX IN THAT CASE!
    !---------------------------------------------------------------------------

    use Spinors

    complex(KIND=dp), allocatable,intent(inout) :: Delta(:,:,:,:)
    complex(KIND=dp), allocatable,intent(inout) :: DeltaLN(:,:,:,:)
    complex(KIND=dp), allocatable,intent(in)    :: Pairingfield(:,:,:,:)
    complex(KIND=dp), allocatable,intent(in)    :: PairingfieldLN(:,:,:,:)
    real(KIND=dp), intent(in)                   :: Gaps(2)
    logical,intent(in)                          :: ConstantGap
    integer                                     :: i, it,j, ii, jj, sig1,sig2, P
    integer                                     :: k
    integer                                     :: minj, maxj, iii, jjj
    type(Spinor)                                :: Psi1, Psi2
    real(KIND=dp)                               :: TempReal(nx,ny,nz)
    real(KIND=dp)                               :: TempIm(nx,ny,nz)
    complex(KIND=dp)                            :: Temp(nx,ny,nz)
    real(KIND=dp)                               :: Cutoff(2)
    logical                                     :: TR

    if(ConstantGap) call stp('Trying to do constant gap pairing in HFB!')

    Delta     = 0.0_dp  
    if(allocated(DeltaLN)) then
      DeltaLN    = 0.0_dp
    endif

    do it=1,Iindex
      do P=1,Pindex
        do i=1,Blocksizes(P,it)
          ii   = Blockindices(i,P,it)
          iii  = mod(ii-1,nwt)+1
          sig1 = HFBasis(iii)%GetSignature()
          if(TRC .and. ii.ne.iii) then
            sig1 = - sig1
          endif

          Cutoff(1) = PCutoffs(iii)

          if(.not. BCSinHFB(it)) then
            !Delta is antisymmetric, only compute the upper-diagonal part.
            !(And Thus Delta_{ii} = 0)
            minj=i+1; maxj=Blocksizes(P,it)
          else
            ! Only compute the 'diagonal' part of Delta_{i, ibar}
            if(Blockindices(i,P,it).le.Blocksizes(P,it)/2) then
                minj = i + Blocksizes(P,it)/2 ; maxj = i + Blocksizes(P,it)/2
            else
                minj = i - Blocksizes(P,it)/2 ; maxj = i - Blocksizes(P,it)/2
            endif
          endif
          do j=minj,maxj
            jj   = Blockindices(j,P,it)
            jjj  = mod(jj-1,nwt)+1
            sig2 = HFBasis(jjj)%GetSignature()
            if(TRC .and. jj.ne.jjj) sig2 = -sig2
            !-------------------------------------------------------------------
            !Selection on quantum numbers:
            ! only explicitly on signature, parity and isospin are taken care of
            ! by the loop structure
            if(sig1.ne.-sig2)  cycle
            Cutoff(2) = PCutoffs(jjj)

            if(Cutoff(1)*Cutoff(2) .lt. HFBNumCut) cycle

            if(TRC .and. ( ii.ne. iii  .or. jj.ne.jjj )  &
            &  .and. &
            &   .not. (ii .ne. iii .and. jj .ne. jjj )) then
              TR = .true.
            else
              ! Should only happen when Time-reversal is conserved, but signature
              ! isn't.
              TR = .false.
            endif
            Temp = GetPairDensity(HFBasis(iii)%Value,HFBasis(jjj)%Value,TR)
            !-------------------------------------------------------------------
            ! Don't forget the complex conjugate here!
            Delta(i,j,P,it) =   dv*Cutoff(1)*Cutoff(2)*                        &
            &            (sum(   DBLE(Temp)  * DBLE( PairingField(:,:,:,it)))  &
            &            +sum(   DIMAG(Temp) * DIMAG(PairingField(:,:,:,it))))
            !Delta is antisymmetric
            Delta(j,i,P,it) = - Delta(i,j,P,it)
            if(allocated(DeltaLN)) then
                DeltaLN(i,j,P,it) =   dv*Cutoff(1)*Cutoff(2)*                  &
                &   (sum(   DBLE(TEMP) * DBLE (PairingFieldLN(:,:,:,it)))      &
                &   +sum(   DIMAG(TEMP)* DIMAG(PairingFieldLN(:,:,:,it))))
                DeltaLN(j,i,P,it) =  - DeltaLN(i,j,P,it)
            endif
          enddo
        enddo
      enddo
    enddo

    !---------------------------------------------------------------------------
  end subroutine HFBGaps

  subroutine HFBGaps_TIMEREV(Delta,DeltaLN,PairingField,PairingFieldLN,Gaps,ConstantGap)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the Delta_{i j} for the HFB model.
    ! The formula (see Ring & Schuck, page 254 eq. (7.41))
    ! is the following:
    !
    ! Delta_{i, j} = \int dr <i,j | Delta(r) | r >
    !
    ! When Lipkin-Nogami is activated, this routine also sums the modified gaps,
    ! which are just analogous to the normal gaps:
    ! Delta_{LN,i,j} = \int dr <i,j | Delta(r)_{LN} | r >
    !---------------------------------------------------------------------------
    ! Gaps & Constantgap are not used in this subroutine, but are present to
    ! be compatible with the BCS GetGaps subroutine header.
    !---------------------------------------------------------------------------

    use Spinors

    complex(KIND=dp), intent(inout) :: Delta(HFBSize,HFBSize,2,2)
    complex(KIND=dp), intent(inout) :: DeltaLN(HFBSize,HFBSize,2,2)
    complex(KIND=dp), intent(in)    :: Pairingfield(nx,ny,nz,2)
    complex(KIND=dp), intent(in)    :: PairingfieldLN(nx,ny,nz,2)
    real(KIND=dp), intent(in)                   :: Gaps(2)
    logical,intent(in)                          :: ConstantGap
    integer                                     :: i, it,j, ii, jj, sig1,sig2, P
    integer                                     :: minj, maxj, iii, jjj
    type(Spinor)                                :: Psi1, Psi2
    real(KIND=dp)                               :: TempReal(nx,ny,nz)
    real(KIND=dp)                               :: TempIm(nx,ny,nz)
    complex(KIND=dp)                            :: Temp(nx,ny,nz)
    real(KIND=dp)                               :: Cutoff(2)
    integer, save                               :: C = 0

    Delta = 0.0_dp ; DeltaLN = 0.0_dp

    do it=1,Iindex
      do P=1,Pindex
        ! On the left, only take wavefunctions of signature +1
        do i=1,Blocksizes(P,it)/2
          ii   = Blockindices(i,P,it)
          iii  = mod(ii-1,nwt)+1

          sig1 = HFBasis(iii)%GetSignature()
          Psi1 = HFBasis(iii)%GetValue()

          Cutoff(1) = PCutoffs(iii)

          do j=Blocksizes(P,it)/2+1, Blocksizes(P,it)
            !On the irght, only take wavefunctions of signature -1

            jj   =   Blockindices(j,P,it)
            jjj  =   mod(jj-1,nwt)+1
            sig2 = - HFBasis(jjj)%GetSignature()

            Cutoff(2) = PCutoffs(jjj)

            if(Cutoff(1)*Cutoff(2) .lt. HFBNumCut) cycle

            Psi2 = HFBasis(jjj)%GetValue()

            !TempReal = - sum(Psi1%Grid * Psi2%Grid)
            TempReal =                                                         &
            &             - Psi1%Grid(:,:,:,3,1) * Psi2%Grid(:,:,:,3,1)        &
            &             - Psi1%Grid(:,:,:,4,1) * Psi2%Grid(:,:,:,4,1)        &
            &             - Psi1%Grid(:,:,:,1,1) * Psi2%Grid(:,:,:,1,1)        &
            &             - Psi1%Grid(:,:,:,2,1) * Psi2%Grid(:,:,:,2,1)

            ! Attention for the extra minus sign: the bra < i,j | indicates an
            ! extra complex conjugation.
            TempIm   =                                                         &
            &             - Psi1%Grid(:,:,:,4,1) * Psi2%Grid(:,:,:,3,1)        &
            &             - Psi1%Grid(:,:,:,3,1) * Psi2%Grid(:,:,:,4,1)        &
            &             - Psi1%Grid(:,:,:,2,1) * Psi2%Grid(:,:,:,1,1)        &
            &             - Psi1%Grid(:,:,:,1,1) * Psi2%Grid(:,:,:,2,1)

            ! Only valid when TimeSimplex is conserved!
            Delta(i,j,P,it) =   dv*Cutoff(1)*Cutoff(2)*                    &
            &            (sum(   TempReal * DBLE(PairingField(:,:,:,it)))  &
            &            -sum(   TempIm   * DIMAG(PairingField(:,:,:,it))))
            Delta(j,i,P,it) = - Delta(i,j,P,it)

            !Only valid when TimeSimplex is conserved
            DeltaLN(i,j,P,it) =   dv*Cutoff(1)*Cutoff(2)*              &
            &   (sum(   TempReal * DBLE (PairingFieldLN(:,:,:,it)))    &
            &   -sum(   TempIm   * DIMAG(PairingFieldLN(:,:,:,it))))
            DeltaLN(j,i,P,it) =  - DeltaLN(i,j,P,it)
          enddo
        enddo
      enddo
    enddo
  end subroutine HFBGaps_TIMEREV

  subroutine HFBoverlapT2(overlapT2)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the matrix of absolute squares of the overlap 
    ! between a state and the time-reversed of another state in order to 
    ! identify partner states in the HF basis. 
    !---------------------------------------------------------------------------
    ! Uses the same storage scheme as Delta(:,:,:,:) in order to be usable
    ! in the same way.
    !---------------------------------------------------------------------------
    ! Note that < psi | T phi > =   < T^+ T psi | T phi >
    !                           =   < T^+ T phi | T^+ psi >
    !                           =   < phi | T^+ psi >
    !                           = - < phi | T psi >
    ! is a complex skew-symmetric matrix (use relations from appendix A of 
    ! my thesis). For the purpose of this routine, however, it is sufficient 
    ! to calculate the real matrix | < psi | T phi > |^2.
    !---------------------------------------------------------------------------
    ! NOTE THAT THIS ROUTINE IS NOT READY FOR TIME-SIMPLEX BREAKING AS DELTA
    ! CAN BE COMPLEX IN THAT CASE!
    !---------------------------------------------------------------------------
    ! MB 19/08/30: for time-reversal-invariant states, this fails to correctly 
    ! compute the overlaps. I don't understand what the code is actually doing,
    ! as only half of the states are in storage, and the 
    ! itself
    !---------------------------------------------------------------------------
    integer                                     :: i, ii , iii , j , jj , jjj
    integer                                     :: it , P
    integer                                     :: sig1, sig2
    integer                                     :: minj, maxj
    real(KIND=dp)                               :: reo , imo
    real(KIND=dp), allocatable,intent(inout)    :: overlapT2(:,:,:,:)
    logical                                     :: TR

    if(.not.allocated(overlapT2)) allocate ( overlapT2(HFBSize,HFBSize,2,2) )

    overlapT2 = 0.0_dp

    ! explicit distinction between cases with and without time-reversal. 
    ! It seems that I (MB) cannot properly handle the clever index juggling
    ! done elsewhere.
    if (TRC) then
      do it=1,Iindex
      do P=1,Pindex
        do i=1,Blocksizes(P,it)/2
          j = Blocksizes(P,it)/2 + i 
          overlapT2(i,j,P,it) = 1.0_dp
          ! use symmetry of the square of overlaps to determine the other half
          overlapT2(j,i,P,it) = overlapT2(i,j,P,it)
        enddo
      enddo
      enddo
    else
      do it=1,Iindex
      do P=1,Pindex
        do i=1,Blocksizes(P,it)
          ii   = Blockindices(i,P,it)
          iii  = mod(ii-1,nwt)+1

          minj=i+1 
          maxj=Blocksizes(P,it)
          do j=minj,maxj
            jj   = Blockindices(j,P,it)
            jjj  = mod(jj-1,nwt)+1

            !-------------------------------------------------------------------
            reo = InproductSpinorReal(HFBasis(iii)%GetValue(),  &
              &                       TimeReverse(HFBasis(jjj)%GetValue()))
            imo = 0.0_dp
            if(.not.TSC) then
            imo = InproductSpinorImaginary(HFBasis(iii)%GetValue(), &
              &                            TimeReverse(HFBasis(jjj)%GetValue()))
            endif
            overlapT2(i,j,P,it) = reo*reo + imo*imo
            ! use symmetry of the square of overlaps to determine the other half
            overlapT2(j,i,P,it) = overlapT2(i,j,P,it)
          enddo
        enddo
      enddo
      enddo
    endif

    !---------------------------------------------------------------------------
  end subroutine HFBoverlapT2

  function HFBNumberofParticles_ordinary(Lambda, Delta, LNLambda) result(N)
    !---------------------------------------------------------------------
    ! Constructs the HFB state ( consisting of RhoHFB and KappaHFB) and
    ! returns the number of particles present.
    !---------------------------------------------------------------------
    real(Kind=dp), intent(in)                :: Lambda(2),LNLambda(2)
    complex(KIND=dp), allocatable,intent(in) :: Delta(:,:,:,:)
    real(KIND=dp), allocatable               :: GaugeEnergies(:,:,:)

    real(KIND=dp)    :: N(2), E, MinE
    complex(KIND=dp) :: Overlaps(nwt,nwt)
    integer          :: it, P, i,j,k,S, ind(2,2), minind(2,2), s1, s2

    call ConstructHFBHamiltonian(Lambda, Delta, LNLambda, HFBGauge)
    !---------------------------------------------------------------------------
    ! Calculate the pfaffian of the hamiltonian in the Majorana representation
    pf        = Pfaffian_HFBHamil()

    !---------------------------------------------------------------------------
    ! reinitialize tracking of swapped QPs
    QPswapped(:,:,:) = 0

    !---------------------------------------------------------------------------
    call DiagonaliseHFBHamiltonian
    !---------------------------------------------------------------------------
    ! Determining which eigenvectors of the HFB Hamiltonian to use
    HFBColumns  = 0
    
    if(PfSolver) then
      !-------------------------------------------------------------------------
      ! Use the number parity of the present HFB state to decide.
      ! The Pfaffian of the HFB Hamiltonian seems to be unreliable for this
      ! when pairing breaks down. It is nevertheless calculated. 
      !-------------------------------------------------------------------------
      ! This requires further future investigation of what actually goes on.
      !-------------------------------------------------------------------------
      ! MB: the use of the Pfaffian to follwo number parity is presently 
      ! (October 2019) a mess, which needs future cleanup after some testing 
      ! why the actual calculation of the pfaffian sometimes fails. 
      ! At present, PfSolver = .true. means that the code follows the number 
      ! parity by actually calculating it in the canonical basis from counting
      ! eigenvalues of 1 and then overriding signs of the pfaffians pf, refpf,
      ! and oldpf in a way that works but became untractable. This needs
      ! some cleaning as what is presently done has nothing to do with Pfaffians
      ! anymore.
      ! MB's working hypothesis is that the pfaffian of the HFB Hamiltonian 
      ! looses its meaning when pairing breaks down.
      !-------------------------------------------------------------------------
      ! From a practical point of view it is worth noting that after thousands
      ! of signature-breaking blocked & cranked calculations I could not
      ! identify a single instance where what is presently done failed. It
      ! works, but should be coded in a more transparent manner.
      !-------------------------------------------------------------------------
      do it=1,Iindex
        do P=1,Pindex
          ! First find all of the positive energy ones
          minE = 2000000
          ind(P,it) = 1
          S = blocksizes(P,it)

          !-------------------------------------------------------------------
          ! MB 19/07/03: it seems that the lowest QP with positive E_qp is 
          ! always the first one at S+1 anyway, at least when HFBGauge = 0, 
          ! which seems natural when diagonalizing a Hamiltonian with pairs
          ! of eigenvalues at +/- E with a routine that returns eigenvectors 
          ! sorted by eigenvalue
          !-------------------------------------------------------------------
   !      do i=1,2*S
   !         if(QuasiEnergies(i,P,it) .gt. 0.0_dp) then
   !            HFBColumns(ind(P,it),P,it) = i
   !            if (QuasiEnergies(i,P,it) .lt. minE) then
   !              minE = QuasiEnergies(i,P,it)
   !              minind(P,it) = ind(P,it)
   !            endif
   !            print '(" scan ",2i2,2i5,1f9.4)',it,P,i,ind(P,it),QuasiEnergies(i,P,it)
   !            ind(P,it) = ind(P,it) + 1
   !         else
   !            print '(" scan ",2i2,1i5,"    *",1f9.4)',it,P,i,QuasiEnergies(i,P,it)
   !         endif
   !      enddo

          !-------------------------------------------------------------------
          ! MB 19/08/22: When HFBGauge is larger than the lowest quasiparticle
          ! energy, looking for the QP with positive E_qp becomes a bug. For
          ! HFBGauge << 0, less than S states will be kept, with the lowest 
          ! QP states missing, and for HFBGauge >> 0, more than S states will 
          ! be kept, with states that should be discarded kept as lowest QPs.
          ! For HFBGauge << 0, we can trust in energy ordering and know that
          ! the empty (in the Dirac sense) QP states are those with indices 
          ! between S+1 and 2S. For HFBGauge >> 0, I do in fact not see how 
          ! to identify the relevant states other than by calculating the
          ! expectation value of the generalized density matrix (which is 
          ! not yet done).
          !-------------------------------------------------------------------
   !      if ( abs(HFBGauge(it)) .gt. 0.001 ) then
   !        do i=1,S
   !          print '(" scan ",2i2,1i5,"    *",1f9.4)',it,P,i,QuasiEnergies(i,P,it)
   !        enddo
   !      endif
          do i=S+1,2*S
            HFBColumns(ind(P,it),P,it) = i
            if (QuasiEnergies(i,P,it) .lt. minE) then
              minE = QuasiEnergies(i,P,it)
              minind(P,it) = ind(P,it)
            endif
   !        if ( abs(HFBGauge(it)) .gt. 0.001 ) then
   !          print '(" scan ",2i2,2i5,1f9.4)',it,P,i,ind(P,it),QuasiEnergies(i,P,it)
   !        endif
            ind(P,it) = ind(P,it) + 1
          enddo
          
          !-------------------------------------------------------------------
          ! estimate number parities of the non-manipulated vacuum 
          !-------------------------------------------------------------------
          call EstimateNumberParitiesVac(P,it)
          ! print '(/," it = ",i2," P = ",i2," BM ",3i5,5x,5i5)',it,P, &
          !      &  HFBNumberParity(P,it),(-1)**NP(P,it),NP(P,it),     &
          !      &  (-1)**NPVac(P,it),NPVac(P,it)

          !-------------------------------------------------------------------
          ! oldpf not yet set (meaning this is the first call of this routine) 
          ! so set it to pf. Its correct sign will be determined next.
          !-------------------------------------------------------------------
          if( oldpf(P,it) .eq. 0 ) then
            oldpf(P,it) = pf(P,it)
          endif

          !-------------------------------------------------------------------
          ! refpf will be a reference pfaffian of a state with correct
          ! number parity.
          !-------------------------------------------------------------------
          refpf(P,it) = pf(P,it)
         
          !-------------------------------------------------------------------
          ! Before an update made in June 2019 by MB, 
          ! the number parity of the non-blocked qp vacuum is followed with 
          ! the help of the sign of the HFB Hamiltonian in Majorana 
          ! representation in each of the isospin-parity blocks. It is 
          ! calculated after initialisation, but there is nothing that 
          ! guarantees a priori that at this stage the reference value 
          ! corresponds to an even particle number in a given block. 
          ! Therefore, when all necessary information is available
          ! the code checks if the estimate of number parity from counting 
          ! eigenvalues of rho in a given block is even or not. When it is 
          ! not, the sign of refpf(P,it) for the block concerned is changed.
          !-------------------------------------------------------------------
          ! Note: during SCF iteration 0, neither refpf(:,:) nor NP(:,:) are
          ! known. refpf(:,:) is therefore set when the HFB solver is called
          ! for the first time in SCF iteration 1. At this moment,  NP(:,:) 
          ! is non-zero for the first time. 
          !-------------------------------------------------------------------
          ! As it turned out, the sign of the pfaffian sometimes also changes
          ! without changing number parity (mainly, possibly exclusively 
          ! when pairing breaks down). The sign of the reference pfaffian is 
          ! changed when the code detects that this has happened. In the end, 
          ! the correct number parity is now detected by actually calculating
          ! it in the canonical basis (which is not 100% reliable), in spite
          ! of the variable "refpf" being used to trace the evolution of
          ! number parity. This is admittedly confusing and requires further
          ! updating whenever it is understood when and why the pfaffian 
          ! method fails.
          !-------------------------------------------------------------------
          if ( (-1)**NPVac(P,it) .ne. 1 ) then
 !          print '(/," Incorrect number parity of nonblocked HFB reference vacuum detected",  & 
 !                 &/," for it = ",i2," P = ",i2,                                              &
 !                 &  " number parity : ",i2," obtained from NPVac = ",1i5)',                  &
 !                 &   it,P,(-1)**NPVac(P,it),NPVac(P,it)
 !          print '(" change sign of refpf= ",1es16.6,/)',refpf(P,it)
            refpf(P,it) = -refpf(P,it)
          endif

          s1 = int(pf(p,it)/abs(pf(p,it)))
          s2 = int(refpf(p,it)/abs(refpf(p,it)))
          !-------------------------------------------------------------------
          ! Something changed, so swap the lowest qp in this block.
          ! QPswapped(P,it,2) is the index of the one swapped away,
          ! QPswapped(P,it,1) is the index of the one swapped into the vacuum.
          !-------------------------------------------------------------------
          if( s1 .ne. s2 ) then
            QPswapped(P,it,2) = HFBColumns(minind(P,it),P,it)
            HFBColumns(minind(P,it),P,it) = 2*S-HFBColumns(minind(P,it),P,it)+1
            QPswapped(P,it,1) = HFBColumns(minind(P,it),P,it)
          ! print '(" swap lowest qp for it = ",i1," P = ",i1," old = ",i5," new ",i5)',it,P,QPswapped(P,it,2),QPswapped(P,it,1)
          endif

          !-------------------------------------------------------------------
          ! estimate number parities of the manipulated vacuum 
          !-------------------------------------------------------------------
          call EstimateNumberParitiesVac(P,it)

 !          print '(" it = ",i2," P = ",i2," AM ",3i5,5x,2i5,/)',it,P, &
 !               &  HFBNumberParity(P,it),(-1)**NP(P,it),NP(P,it),     &
 !               &  (-1)**NPVac(P,it),NPVac(P,it)
        enddo
      enddo
    else
      !-------------------------------------------------------------------------
      ! Do as in CR8: take only the first half of the matrix, and conjugate
      ! the first quarter.
      if(SC) then
          do it=1,Iindex
              do P=1,Pindex
                  S = blocksizes(P,it)
                  do i=1,S/2
                      HFBColumns(i,P,it) = 2*S - i + 1
                  enddo
                  do i=S/2+1,S
                      HFBColumns(i,P,it) = i
                  enddo
              enddo
          enddo
      else
          !---------------------------------------------------------------------
          ! If signature is not conserved, just take all the positive
          ! qp energies and pray it will work.
          ind = 1
          do it=1,Iindex
          do P=1,Pindex
            S = blocksizes(P,it)
            do i=1,2*S
              if (QuasiEnergies(i,P,it) .gt. 0.0_dp) then
                HFBColumns(ind(P,it),P,it) = i
                ind(P,it) = ind(P,it) + 1
               endif
            enddo
          enddo
          enddo
      endif
    endif

    !---------------------------------------------------------------------------
    ! Block some quasiparticles
    if(allocated(qpexcitations)) then
      call BlockQuasiParticles(Lambda)
    endif
    call constructRhoHFB(HFBColumns)

    ! MB 19/07/10 Occupations(:,:,:) have not been recalculated at this point ...
    do it=1,Iindex
      do P=1,Pindex
        ind(P,it) = 0
        do i=1,Blocksizes(P,it)
          if(abs(Occupations(i,P,it) - 1) .lt. 1d-10) then
            ind(P,it) = ind(P,it) + 1
          endif
        enddo
      enddo
    enddo

    call constructKappaHFB(HFBColumns)
    !---------------------------------------------------------------------------
    ! Calculate the number of particles in this configuration and return.
    ! The Fermi energy can then use this to get adjusted.
    N = 0.0_dp
    do it=1,Iindex
      do P=1,Pindex
        do i=1,blocksizes(P,it)
          N(it) = N(it) + real(RhoHFB(i,i,P,it))
        enddo
      enddo
    enddo

  end function HFBNumberofParticles_ordinary

  subroutine HFBFindFermiEnergyBisection(Fermi,L2,Delta,DeltaLN,Lipkin,DN2,    &
    &                                                  ConstrainDispersion,Prec)
  !-----------------------------------------------------------------------------
  ! New bisection routine targeting pure HFB case. The previous solver that can 
  ! also handle LN is now HFBLNFindFermiEnergyBisection) and is still called
  ! for HFB+LN. 
  !-----------------------------------------------------------------------------
  ! The search for Fermi energies is now done in two steps.
  ! 
  ! The first step consists in bracketing the Fermi energy, i.e. finding an
  ! interval [A,B] with A < eps_F < B with N(A) - N_0 < 0 and N(B) - N_0 > 0
  !
  ! The algorithm is as follows:
  !
  ! 1) calculate <N> for Fermi energy eps_o found in storage
  !
  ! 2) if N(eps_o) = N_0 : don't look further, set A = B with N(A) = N(B) = N_0
  !    if N(eps_o) < N_0 : A = eps_o, look for B = eps_o + teststep
  !    if N(eps_o) > N_0 : B = eps_o, look for A = eps_o - teststep
  !
  ! with teststep = 0.01 MeV (see below)
  !
  ! 3) if N(A)*N(B) < 0 : Fermi energy is bracketed by [A,B], don't look further
  !    if N(A)*N(B) > 0 : Fermi energy is not bracketed yet by [A,B]
  ! 
  ! The stepsize of further updates grows with the number of failed updates. 
  ! A reasonable choice seems to be 
  !
  !            step(iter) = teststep * iter = 0.01 MeV * iter , 
  !
  ! leading to step(1) = 0.01, step(2) = 0.03, step(3) = 0.06, step(4) = 0.1, 
  ! step(5) = 0.21, ...  step(76) = 30.03 MeV (then the routine stops the code).
  ! Except during the initial phase, the Fermi energy typically changes on the
  ! few keV scale, such that the first step usually brackets the zero.
  !      if B is moving boundary : A = B , B = eps_o + teststep * iter
  !      if A is moving boundary : B = A , A = eps_o - teststep * iter
  ! repeat step 3 until convergence. 
  !-----------------------------------------------------------------------------
  ! Note: Reducing the first step to 0.005 very often requires two bracketing 
  ! steps, which is why I have put 0.01 as starting point for time being.
  ! Note: starting with eps_o +/- 0.005, one reaches only eps_o +/- 15.02 
  ! in 76 iterations
  !-----------------------------------------------------------------------------
  ! The second step uses Brent's method that switches between the secant method
  ! and inverse quadratic interpolation. This is done by call of the subroutine 
  ! FermiBrent()
  !-----------------------------------------------------------------------------
  ! The routine can also be called for active LN. In that case, it calculates
  ! the lambda2 once at the call of the routine from the present pairing fields
  ! without further updates during the iterations of the Fermi energy. With 
  ! this, the lambda2 will slowly converge with the HF iterations. As the 
  ! lambda2 usually will change only by small amounts anyway, this even might 
  ! stabilise the iterations.
  !-----------------------------------------------------------------------------
  ! In the present context of MOCCa this routine is unnecessarily complicated by
  ! having to calculate proton and neutron number simultaneously with 
  ! HFBNumberofParticles(), as it anticipated proton-neutron mixing.
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
  complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
  complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
  logical, intent(in)                       :: Lipkin,COnstrainDispersion
  real(KIND=dp), intent(in)                 :: Prec, DN2(2)

  real(KIND=dp) :: N(2), Num(2)
  real(KIND=dp) :: InitialBracket(2,2), FA(2), FB(2)
  real(KIND=dp) :: A(2), B(2)
  logical       :: NotFound(2)
  logical       :: Success
  integer       :: iter, FailCount, flag(2), it, i
  integer       :: idir(2) = 0 , idirsig(2) = 1

  !-----------------------------------------------------------------------------
  ! calculate the lambda2(it) from present occupations and pairing fields in 
  ! case of Lipkin-Nogami.
  !-----------------------------------------------------------------------------
  if (Lipkin) then
    flag = 0
    L2 = LNCR8(Delta, DeltaLN, flag)
  endif

  !-----------------------------------------------------------------------------
  ! The targeted particle numbers
  !-----------------------------------------------------------------------------
  Particles(1) = Neutrons
  Particles(2) = Protons

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  do it=1,2
    if ( Particles(it) .lt. 0.1_dp ) Fermi(it) = -100.0_dp
  enddo

  !-----------------------------------------------------------------------------
  ! Check if initial guess is already good enough. If so, this will be the 
  ! only call of the HFB solver.
  !-----------------------------------------------------------------------------
  N = HFBNumberofParticles(Fermi, Delta, L2 ) - Particles

  ! Initialise oldpf with present pf (first call only)
  ! NOTE: the number parity might be incorrect for the initialisation, such that
  ! one has to double-check.
  if( all(oldpf .eq. 0)) then
    oldpf = Pfaffian_HFBHamil()   
  endif

  if (all(abs(N) .lt. Prec)) return

  ! here would start the loop over iterations of lambda_2 in the LN case if
  ! an actual iteration is needed to improve convergence.

    Success  = .false.
    !---------------------------------------------------------------------------
    ! Use present Fermi energy as starting point and check the direction 
    ! where the zero of <N>-N0 can be expected.
    ! If <N>-N0 <  0, search at higher values.
    ! If <N>-N0 >= 0, search at lower  values.
    ! Initialize InitialBracket(it,1) = A, InitialBracket(it,2) = B with present
    ! Fermi energy.
    ! "dir" is the label of the InitialBracket(it,idir) that has to be moved,
    ! "idirsig" is the sign of steps needed to go into that direction.
    !---------------------------------------------------------------------------
    do it=1,2
      InitialBracket(it,:) = Fermi(it) 
      if ( N(it) .lt. 0.0_dp ) then 
        idir   (it) =  2
        idirsig(it) =  1
      else 
        idir   (it) =  1
        idirsig(it) = -1
      endif
    enddo

    !---------------------------------------------------------------------------
    ! Try to find a boundary that brackets the Fermi energy in the direction 
    ! into which the Fermi energy has to be changed.
    !---------------------------------------------------------------------------
    FailCount = -1
    NotFound(:) = .true. 

    do while(.not. Success)
      Success  = .true.
      !-------------------------------------------------------------------------
      ! check if particle number for one species takes already the targeted 
      ! value. If so, do not look further. If the particle numbers were correct
      ! for both isospins, the routine would already have been left above.
      ! Note: this test is only relevant for the first iteration.
      !-------------------------------------------------------------------------
      do it=1,2
        if ( abs(N(it)) .lt. Prec ) NotFound(it) = .false. 
      enddo
      FailCount = FailCount + 1
      !-------------------------------------------------------------------------
      ! update moving boundary and recalculate particle numbers at both.
      !-------------------------------------------------------------------------
      ! If the particle number is already correct for the starting point (as 
      ! might happen in the HF case or for near-converged states), NotFound(it) 
      ! is already .false. such that the moving boundary will not be updated, 
      ! leading to A = B with N(A) = N(B) approx N_0 for this isospin.
      !-------------------------------------------------------------------------
      do it=1,2
        if ( NotFound(it) ) InitialBracket(it,idir(it)) &
          & = InitialBracket(it,idir(it)) + idirsig(it) * 0.01_dp*(FailCount+1)
      enddo
      FA(:) = HFBNumberofParticles(InitialBracket(:,1),Delta,L2) - Particles(:)
      FB(:) = HFBNumberofParticles(InitialBracket(:,2),Delta,L2) - Particles(:)
      Num(:) = FB(:) + Particles(:)

      !========================================================================
      ! check if N(epsilon_F) is a monotonically growing function.
      ! It should be, but the code faces cases where it is not.
      !------------------------------------------------------------------------
      ! NOTE: That there is a region where N(epsilon_F) is monotonically 
      ! decreasing does not automatically mean that there are three (or even 
      ! more) zeros of N(epsilon_F) - N_0. Many of the examples I used to
      ! test workarounds still had one zero only, with a wide plateau either
      ! slightly above or below N(epsilon_F) - N_0 = 0.  
      !------------------------------------------------------------------------
      do it=1,2
        if ( FB(it) .lt. FA(it) ) then 
          !===================================================================
          ! we are in a state of emergency for this isospin!
          !===================================================================
          print '(" HFBFindFermiEnergyBisection: Warning N(eps_F) decreases for it = ",i2)',it
          print '(" A = ",f13.8," FA = ",f14.8," B = ",f13.8," FB = ",f14.8)', &
               & InitialBracket(it,1),FA(it)+Particles(it), &
               & InitialBracket(it,2),FB(it)+Particles(it)
          if ( HandleBisectionEmergencies ) then
            !------------------------------------------------------------------
            ! Make a guess about suitable action
            !------------------------------------------------------------------
            A(:) = InitialBracket(:,1)
            B(:) = InitialBracket(:,2)
            call HFBFindFermiEnergyBisectionEmergency(Fermi,L2,Delta,DeltaLN, &
                                         &  Lipkin,DN2,          &
                                         &  ConstrainDispersion,Prec, &
                                         &  A,B,FA,FB,it,NotFound)
            InitialBracket(:,1) = A(:)
            InitialBracket(:,2) = B(:)
            print '(" New interval set to ",4f13.8)', &
                    &    InitialBracket(it,1),FA(it), &
                    &    InitialBracket(it,2),FB(it)
          endif
        endif
        !-----------------------------------------------------------------------
        ! check if root is bracketed for isospin it after the update
        ! As the function might be constant, this has to be .le., not .lt.
        !-----------------------------------------------------------------------
        if( NotFound(it) ) then 
          if( FA(it)*FB(it) .le. 0.0_dp ) then 
            NotFound(it) = .false.
          endif
        endif
      enddo
      !-------------------------------------------------------------------------
      ! Check if root is bracketed for both isospins. If so, loop will be left.
      !-------------------------------------------------------------------------
      if ( any(NotFound) )  Success = .false.

      !-------------------------------------------------------------------------
      ! diagnostic printing for convergence analysis (usually commented out)
      !-------------------------------------------------------------------------
      ! For a convergence analysis, uncomment this and the similar statement in
      ! the subroutine FermiBrent().
      !-------------------------------------------------------------------------
!       print '(" Bracketing ",i4,1l2,2(1l2,2(f13.8,es16.7),f14.8))', &
!            & FailCount,Success,                                     &
!            & NotFound(1),InitialBracket(1,1),FA(1),                 &
!            &             InitialBracket(1,2),FB(1),Num(1),          &
!            & NotFound(2),InitialBracket(2,1),FA(2),                 &
!            &             InitialBracket(2,2),FB(2),Num(2)
      !-------------------------------------------------------------------------
      ! code failure (Fermi energy has changed by 30 MeV)
      !-------------------------------------------------------------------------
    ! if (Failcount .gt. 76) then
      if (Failcount .gt. 130) then
        print '(/," A(n) = ", f13.8, " FA(n) = ",1es12.4,             &
             &    " B(n) = ", f13.8, " FB(n) = ",1es12.4,             &
             &    " A(p) = ", f13.8, " FA(p) = ",1es12.4,             &
             &    " B(p) = ", f13.8, " FB(p) = ",1es12.4)',           &
             &  InitialBracket(1,1),FA(1),InitialBracket(1,2),FB(1),  &
             &  InitialBracket(2,1),FA(2),InitialBracket(2,2),FB(2) 
        call stp('HFBFindFermiEnergyBisection: Search for InitialBracket failed.')
      endif
      !-------------------------------------------------------------------------
      ! when zero is not bracketed, let the kept boundary follow the moving one 
      ! before doing the next step
      !-------------------------------------------------------------------------
      do it=1,2
        if( NotFound(it) ) then 
          InitialBracket(it,3-idir(it)) = InitialBracket(it,idir(it))
        endif
      enddo
    enddo

    !---------------------------------------------------------------------------
    ! call solver for Brent's method fo find the Fermi energy.
    !---------------------------------------------------------------------------
    ! that there is nothing to be iterated further for one isospin is 
    ! communicated to FermiBrent() as A(it) = B(it).
    ! This is particularly relevant for the breakdown of HFB to HF, where 
    ! N(Lambda) - N_0 is zero over a large interval, a case that cannot be 
    ! handled by FermiBrent().
    !---------------------------------------------------------------------------
    call FermiBrent(InitialBracket(:,1),InitialBracket(:,2),FA,FB,HFBIter, &
                 & Delta,L2,Prec,Fermi,N)

  end subroutine HFBFindFermiEnergyBisection

  subroutine HFBFindFermiEnergyBisectionEmergency(Fermi,L2,Delta,DeltaLN, &
                                         &  Lipkin,DN2,          &
                                         &  ConstrainDispersion,Prec, & 
                                         &  A,B,FA,FB,it,NotFound)
  !----------------------------------------------------------------------------
  ! Various optional measures taken when there is a problem with finding roots
  ! of N(e_fermi) - N0 in the subroutine HFBFindFermiEnergyBisection because 
  ! the function is decreasing over some short interval with the possible
  ! consequence of the function having several zeros. The call of this routine 
  ! is activated by the flag HandleBisectionEmergencies. It makes a guess
  ! about the most plausible root that leads to the targeted configuration.
  !----------------------------------------------------------------------------
  ! At time being, this subroutine always prints some diagnostics of what is
  ! happening and what it is doing about it.
  !----------------------------------------------------------------------------
  real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
  complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
  complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
  logical, intent(in)                       :: Lipkin,COnstrainDispersion
  real(KIND=dp), intent(in)                 :: Prec, DN2(2)
  integer, intent(in)                       :: it
  real(KIND=dp) :: N(2), Num(2)
  real(KIND=dp) :: C(2), D(2), E(2), F(2), FC(2), FD(2), FE(2), FF(2)
  real(KIND=dp) :: XMINA, MINA, XZA(2), ZA(2), XMAXB, MAXB, XZB(2), ZB(2)
  real(KIND=dp) :: Estep = 0.05_dp
  real(KIND=dp) :: A(2), B(2), FA(2), FB(2)
  real(KIND=dp) :: eFermiFuzzy 
  logical       :: NotFound(2)
  logical       :: GuessFuzzyFermi
  integer       :: i, icase

          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! old Fermi energy as a reminder
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          print '(" Old  : ",i3,3x,3f13.8)',it,Fermi(it),N(it)

          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! check what happens below in steps of Estep = 0.05 MeV until the
          ! slope is correct and N(epsilon_F) - N_0 is negative.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          XMAXB  = A (it)
          MAXB   = 0.0_dp
          XZB(:) = 0.0_dp
          ZB (:) = 0.0_dp
          XMINA  = B (it)
          MINA   = 0.0_dp
          XZA(:) = 0.0_dp
          ZA (:) = 0.0_dp

          D (:)  = A (:)
          FD(:)  = FA(:)
          print '(" Below: ",2i3,12f9.5)',it,0,A(it),FA(it), &
                                           & XMAXB,MAXB,XZB(:),XMINA,MINA,XZA(:)
          do i=1,21
            C   (it) = D  (it) - i*Estep
            C (3-it) = D(3-it)
            FC(:)    = HFBNumberofParticles(C, Delta, L2) - Particles(:)
            if ( (FC(it) - FA(it)) .gt. MAXB ) then
              XMAXB = C (it) 
              MAXB  = FC(it) - FA(it)
            endif
            if ( FA(it) .lt. 0.0_dp )  then 
              if ( FA(it)*FC(it) .lt. 0.0_dp ) then ! we crossed zero 
                XZB(1) = C (it)
                ZB (1) = FC(it)
              endif
              if ( FA(it)*FC(it) .gt. 0.0_dp .and. &
                  ZB(1) .ne. 0.0_dp .and. ZB(2) .eq. 0.0_dp ) then
                XZB(2) = C (it)
                ZB (2) = FC(it)
              endif
            endif
            print '(" Below: ",2i3,12f9.5)',it,i,C(it),FC(it), &
                                           & XMAXB,MAXB,XZB(:),XMINA,MINA,XZA(:)
            if ( FC(it) .lt. 0 .and. FD(it) .gt. FC(it) ) exit
            D (:) = C (:)
            FD(:) = FC(:)
          enddo
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! check what happens above in steps of Estep = 0.05 MeV until the
          ! slope is correct and N(epsilon_F) - N_0 is positive.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          E (:)  = B (:)
          FE(:)  = FB(:)
          print '(" ABove: ",2i3,12f9.5)',it,0,B(it),FB(it), &
                                           & XMAXB,MAXB,XZB(:),XMINA,MINA,XZA(:)
          do i=1,21
            F (  it) = E  (it) + i*Estep
            F (3-it) = E(3-it)
            FF(:)    = HFBNumberofParticles(F, Delta, L2) - Particles(:)
            if ( (FF(it) - FB(it)) .lt. MINA ) then
              XMINA = F (it)
              MINA  = FF(it) - FB(it)
            endif
            if ( FB(it) .gt. 0.0_dp ) then
              if ( FB(it)*FF(it) .lt. 0.0_dp ) then
                XZA(1) = F (it)
                ZA (1) = FF(it)
              endif
              if ( FB(it)*FF(it) .gt. 0.0_dp .and. & 
                  & ZA(1) .ne. 0.0_dp .and. ZA(2) .eq. 0.0_dp ) then
                XZA(2) = F (it)
                ZA (2) = FF(it)
              endif
            endif
            print '(" ABove: ",2i3,12f9.5)',it,i,F(it),FF(it), &
                                           & XMAXB,MAXB,XZB(:),XMINA,MINA,XZA(:)
            if ( FF(it) .gt. 0 .and. FF(it) .gt. FE(it) ) exit
            E (:) = F (:)
            FE(:) = FF(:)
          enddo
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! at this point:
          !   Fermi : old Fermi energy
          !   N     : N(Fermi)-N0
          !   MAXB  : largest  value of N(eps)-N0 found below A
          !   XMAXB : energy giving MAXB
          !   XZB(1): energy where N(eps)-N0 crossed 0 in wrong direction below 
          !   MINA  : smallest value of N(eps)-N0 found above B
          !   XMINA : energy giving MINA
          !   XZA(1): energy where N(eps)-N0 crossed 0 in wrong direction above
          !   the energies are C < A < B < F
          !   D = C + Estep, E = F - Estep
          !   FA < FB (which indicates a problem)
          !   FC < 0 with correct slope
          !   FF > 0 with correct slope
          ! possible cases:
          !   1.) FC < 0, FA < 0, FB < 0, FF > 0, FC*FB > 0, FA*FF < 0
          !   2.) FC < 0, FA < 0, FB < 0, FF > 0, FC*FB > 0, FA*FF < 0, XZB
          !   3.) FC < 0, FA > 0, FB < 0, FF > 0, FC*FB > 0, FA*FF > 0
          !   4.) FC < 0, FA > 0, FB > 0, FF > 0, FC*FB < 0, FA*FF > 0 
          !   5.) FC < 0, FA > 0, FB > 0, FF > 0, FC*FB < 0, FA*FF > 0, XZA
          ! In cases 1 and 4 there is one zero, in cases 2, 3, 5 there are three
          ! zeros.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          GuessFuzzyFermi = .false.
          eFermiFuzzy     = 10000.0_dp
          icase           = 0
          if ( FC(it)*FB(it) .gt. 0 .and. FA(it)*FF(it) .lt. 0 ) then
            if ( abs(XZB(1)) .gt. 0 ) then 
              icase = 2
            else
              icase = 1
            endif
          endif
          if ( FC(it)*FB(it) .gt. 0 .and. FA(it)*FF(it) .gt. 0 ) icase = 3
          if ( FC(it)*FB(it) .lt. 0 .and. FA(it)*FF(it) .gt. 0 ) then
            if ( abs(XZA(1)) .gt. 0 ) then 
              icase = 5
            else
              icase = 4
            endif
          endif
          !--------------------------------------------------------------------
          ! How to get out of this?
          ! case 1) only one zero above B
          !     1a) it is nearby old eps_F, continue bracketing in [E,F]
          !     1b) If it is not, take emergency measure (see below)
          ! case 2) two zeros below A and one zero above B. 
          !     2a) the one above is near old eps_F, continue in [E,F]
          !     2b) the lowest one is near old eps_F, continue in [C,D]
          ! case 3) 3 zeros: one below A, one in [A,B] one above B
          !     3a) if (FF-FB) << (FA-FC), continue bracketing in [E,F]
          !     3b) if (FF-FB) >> (FA-FC), continue bracketing in [C,D]
          !     3c) if it is neither, take emergency measure (see below)
          ! case 4) only one zero below A
          !     4a) the one below is near old eps_F, continue in [C,D]
          !     4b) If it is not, take emergency measure (see below)
          ! case 5) one zero below A and two zeros above B
          !     5a) the one below is near old eps_F, continue in [C,D]
          !     5b) the highest one is near old eps_F, continue in [E,F]
          ! For time being I don't consider taking the zero on the negative 
          ! slope found in cases 2, 3, and 5
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! As an emergency measure there are two possibilities:
          !     a) determine eps_F at all cost by continuing bracketing on the
          !        correct side. 
          !     b) make an educated guess for a "fuzzy" Fermi energy that
          !        gives a not too much incorrect particle number
          ! This case seems to emerge mainly when some qp is blocked, which
          ! might seriously impact the eps_F-dependence of particle number.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! Note: for time being, the information about the height of the
          ! "barriers" or "troughs" separating the zeros is not used for the 
          ! decision, as for all test cases I looked at distance of zeros 
          ! were clearly correlated with the height of the barriers. Future
          ! use cases might need to consider them as well.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! Note: the present treatment of the cases could be simplified, as
          ! the zero on the downsloping part of the curve is always excluded.
          ! Not clear if this is the final word, so I kept the complicated
          ! detailed distinction of cases.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          select case(icase)
          case (1)
            ! one zero above: C < D < A < B < E < eps < F 
            print '(" case 1 ",8f13.8)',abs(F(it)-Fermi(it)),MINA,MAXB
            if ( F(it) .lt. Fermi(it) ) then
              ! the descending N(eps) poses no problem as the old Fermi energy 
              ! was even further above than the region where everything is 
              ! normal, so continue there (this would also be found without 
              ! calling this subroutine)
              A(it)=E(it) ; B(it)=F(it) ; FA(it)=FE(it) ; FB(it)=FF(it)
            else
              ! old fermi energy was below the upper limit
              if ( abs(F(it)-Fermi(it)) .lt. 5.0_dp * Estep ) then
                ! no large jump, so continue bracketing 
                A(it)=E(it) ; B(it)=F(it) ; FA(it)=FE(it) ; FB(it)=FF(it)
              else
                ! large jump
                if ( AllowFuzzyNumber ) then
                  ! use Fuzzy particle number
                  if (  abs(XMAXB-Fermi(it))  .le.  5.0_dp * Estep ) then
                    ! use estimate for energy giving smallest N above [A,B]
                    eFermiFuzzy = XMAXB + 0.5_dp*(XMAXB-A(it))
                  else
                    GuessFuzzyFermi = .true. ! (which is the case already)
                  endif
                else
                  ! continue bracketing at all cost
                  A(it)=E(it) ; B(it)=F(it) ; FA(it)=FE(it) ; FB(it)=FF(it)
                endif
              endif
            endif

          case (2)
            ! two zeros below, one zero above: 
            ! FA and FB are negative.
            ! C < D < eps1 < eps2 < A < B < E < eps3 < F
            ! (eps1 might also be in [C,D])
            ! At time being, the intermediate one on the wrong slope
            ! is avoided when AllowFuzzyNumber is active.
            print '(" case 2 ",8f13.8)', &
              &  abs(F(it)-Fermi(it)),abs(XZB(1)-Fermi(it)), &
              &  abs(XZB(2)-Fermi(it)),MINA
            print '(" zeros near :",8f13.8)', & 
               &  C(it),D(it),XZB(2),XZB(1),A(it),B(it),E(it),F(it)

            if ( abs(F(it)-Fermi(it)) .lt.  abs(C(it)-Fermi(it)) ) then
              A(it) = E(it) ; B(it) = F(it) ; FA(it) = FE(it) ; FB(it) = FF(it)
            else
              A(it) = C(it) ; B(it) = D(it) ; FA(it) = FC(it) ; FB(it) = FD(it)
            endif

          case (3)
            ! three zeros
            ! C < eps1 < D < A < eps2 < B < E < eps3 < F
            print '(" case 3 ",8f13.8)', &
            &       abs(F(it)-Fermi(it)),MAXB,abs(F(it)-Fermi(it)),MINA
            print '(" zeros near :",8f13.8)', & 
               &  C(it),D(it),XZB(1),A(it),B(it),XZA(1),E(it),F(it)

            if ( abs(F(it)-Fermi(it)) .lt.  abs(C(it)-Fermi(it)) ) then
              A(it) = E(it) ; B(it) = F(it) ; FA(it) = FE(it) ; FB(it) = FF(it)
            else
              A(it) = C(it) ; B(it) = D(it) ; FA(it) = FC(it) ; FB(it) = FD(it)
            endif

          case (4)
            ! one zero below
            print '(" case 4 ",8f13.8)',abs(C(it)-Fermi(it)),MAXB,MINA
            if ( C(it) .gt. Fermi(it) ) then
              ! the descending N(eps) poses no problem as the old Fermi energy 
              ! was even further below than the region where everything is 
              ! normal, so continue there (this would also be found without 
              ! calling this subroutine)
              A(it)=C(it) ; B(it)=D(it) ; FA(it)=FC(it) ; FB(it)=FD(it)
            else
              ! old fermi energy was above the lower limit
              if ( abs(C(it)-Fermi(it)) .lt. 4.0_dp * Estep ) then
                ! no large jump, so continue bracketing 
                A(it)=C(it) ; B(it)=D(it) ; FA(it)=FC(it) ; FB(it)=FD(it)
              else
                ! large jump
                if ( AllowFuzzyNumber ) then
                  ! use Fuzzy particle number
                  if (  abs(XMINA-Fermi(it))  .le.  5.0_dp * Estep ) then
                    ! use estimate for energy giving smallest N above [A,B]
                    eFermiFuzzy = B(it) + 0.5_dp*(XMINA-B(it))
                  else
                    ! we are lost, so make a guess
                    GuessFuzzyFermi = .true. ! (which is the case already)
                  endif
                else
                  ! continue bracketing at all cost
                  A(it)=C(it) ; B(it)=D(it) ; FA(it)=FC(it) ; FB(it)=FD(it)
                endif
              endif
            endif

          case (5)
            ! one zero below, two zeros above: 
            ! C < eps1 < D < A < B < eps2 < eps3 < E < F
            ! (eps 2  and eps3 might also be in [E,F])
            ! At time being, the intermediate one on the wrong slope
            ! is avoided when AllowFuzzyNumber is active.
            print '(" case 5 ",8f13.8)', &
              &  abs(C(it)-Fermi(it)),abs(XZA(1)-Fermi(it)), &
              &  abs(XZA(2)-Fermi(it)),MINA
            print '(" zeros near :",8f13.8)', &
               &  C(it),D(it),A(it),B(it),XZA(1),XZA(2),E(it),F(it)

            if ( abs(F(it)-Fermi(it)) .lt.  abs(C(it)-Fermi(it)) ) then
              A(it) = E(it) ; B(it) = F(it) ; FA(it) = FE(it) ; FB(it) = FF(it)
            else
              A(it) = C(it) ; B(it) = D(it) ; FA(it) = FC(it) ; FB(it) = FD(it)
            endif

          case default
            if ( AllowFuzzyNumber ) GuessFuzzyFermi = .true.

          end select

          !--------------------------------------------------------------------
          ! Use a point half way between the boundaries as emergency Fermi 
          ! energy. Taking the average stirs some changes in occupations such 
          ! that the code has a chance to find a better solution next time. 
          ! Using the Fermi energy that gives the smallest deviation might 
          ! lock the code at some value of the Fermi energy that does not
          ! change anymore.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The manipulations set the boundaries such that the subsequent
          ! Brent search will not iterate for this isospin.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! Note: particle number will be incorrect in the densities used in 
          ! the following HF iteration.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          if ( GuessFuzzyFermi .and. AllowFuzzyNumber ) then
            A (it) = A (it) + 0.5*( B(it) - A(it) )
            B (it) = A (it)
            F (it) = A (it)
            FF (:) = HFBNumberofParticles(F, Delta, L2) - Particles(:)
            FA(it) = FF(it)
            FB(it) = FA(it)
            print '(" Fuzzy number guessed for it = ",i2," with e_F = ",  &
                    & f12.8," and N = ",f14.8)',it,A(it),FB(it)+Particles(it)
            NotFound(it) = .false.
          endif
          if ( eFermiFuzzy .lt. 9999.9_dp ) then
            A (it) = eFermiFuzzy
            B (it) = A (it)
            F  (:) = A (:)
            FF (:) = HFBNumberofParticles(F, Delta, L2) - Particles(:)
            FA(it) = FF(it)
            FB(it) = FA(it)
            print '(" Fuzzy number for it = ",i2," with e_F = ",f12.8,  &
                    & " and N = ",f14.8)',it,A(it),FB(it)+Particles(it)
            NotFound(it) = .false.
          endif

  end subroutine HFBFindFermiEnergyBisectionEmergency

  subroutine HFBLNFindFermiEnergyBisection(Fermi,L2,Delta,DeltaLN,Lipkin,DN2,  &
    &                                                  ConstrainDispersion,Prec)
  !-----------------------------------------------------------------------------
  ! MB 19/20/26: with the adaptation of HFBFindFermiEnergyBisection to handle
  ! the LN case, this subroutine is now obsolete and subject to future 
  ! suppression from the code.
  !-----------------------------------------------------------------------------
  ! The idea is simple. Note:
  !
  ! F  =  Fermi Level
  ! L2 =  LN parmeter
  ! N( F, L2) = deviation of particle number, as a function of F and L2
  !
  ! Note that also, for a given L2, the function N_{L2 fixed}(F) is monotonously
  ! ascending. Clearly, when F => -infinity, N = 0 and when F => infinity
  ! N = infinity (or in practice the number of particle states considered.)
  ! Thus, N_{L2 fixed}(F) has exactly ONE zero.
  !
  ! Now, we can guarantee finding the correct F for zeroing N_{L2 fixed}(F) by
  ! using a 1D bisection algorithm. (One can also use faster, less surefire
  ! methods). Denote this special F by F(L2).
  ! This in fact defines a 1D curve as a function of L2.
  !
  ! N(L2, F(L2)) = N(L2)
  !
  ! Thus now our problem is reduced to a one-dimensional one, which again can
  ! be bisected or treated more cleverly.
  !
  ! G(F, L2) = L2_{computed} - L2_{varied}
  !
  ! Here, we implement two methods for every 1D optimisation. For every single-one
  ! we have the secant method, which is usually enough. However, for every level
  ! we also implement a back-up bisection method in case the secant method fails.
  ! This should be enough to always find a solution to our problem.
  !-------------------------------------------------------------------------------
  ! This routine also tries to get out of problems by guessing KappaHFB for
  ! an isospin index for which pairing has disappeard. This function is only
  ! active when Lipkin-Nogami is active.
  ! TODO:
  !   decide if this guessing should also be done when no LN is present.
  !-------------------------------------------------------------------------------
  ! Notes:
  !
  ! 1) There is once case which merits special mention: if there are multiple
  !    zeros (as a function of L2). This method will not necessarily pick one
  !    of them consistently. (In some sense of `one' and `consistent', since
  !    the problem is different at every mean-field iteration).
  !    However, in practice I have never encountered this situation. M. Bender
  !    told me that B. Avez encountered this, but since there absolutely no notes
  !    or any other trace left(except a VERY WEIRD subroutine in CR8 that was
  !    never commented or, I suspect, finished), at the moment I ignore
  !    this possibility.
  !    In any case, I suspect that close to convergence, this routine will
  !    not `oscillate' between different L2, since this routine tries to stay
  !    close to the value of the previous mean-field iteration.
  !
  !-------------------------------------------------------------------------------

  real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
  complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
  complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
  logical, intent(in)                       :: Lipkin,COnstrainDispersion
  real(KIND=dp), intent(in)                 :: Prec, DN2(2)

  real(KIND=dp) :: InitialBracket(2,2), FA(2), FB(2)
  real(KIND=dp) :: L2Old(2)
  real(KIND=dp) :: tempL2(2), G(2), GOld(2), N(2)
  logical       :: Succes
  integer       :: iter, FailCount, flag(2), it

  Particles(1) = Neutrons
  Particles(2) = Protons

  if(Lipkin) then
    flag=0
    L2Old = L2 +0.001
    N = HFBNumberofParticles(Fermi, Delta, L2Old ) - Particles
    GOld  = LNCR8(Delta, DeltaLN,flag) - L2Old
  endif

  !Check if initial guess is good enough.
  N = HFBNumberofParticles(Fermi, Delta, L2 ) - Particles
  if(all(abs(N) .lt. Prec)) return
  do iter=1,HFBIter
    InitialBracket(:,1) = Fermi
    InitialBracket(:,2) = Fermi
    ! Try to find an initial interval that brackets the Fermi energy
    Succes = .false.
    FailCount = -1

    do while(.not. Succes)
      Succes = .true.
      FailCount = FailCount + 1
      InitialBracket(:,1) = InitialBracket(:,1) - 1_dp
      InitialBracket(:,2) = InitialBracket(:,2) + 1_dp

      FA = HFBNumberofParticles(InitialBracket(:,1), Delta, L2 ) - Particles
      FB = HFBNumberofParticles(InitialBracket(:,2), Delta, L2 ) - Particles

      if(any(FA*FB.gt.0.0_dp)) then
        Succes=.false.
      endif

      if(Failcount.gt.50) then
        print *, InitialBracket
        print *, Fa, FB
        print *, HFBNumberofParticles(InitialBracket(:,1), Delta, L2 ) - Particles
        call stp('Failed a lot.')
      endif
    enddo

    Fermi = FermiBisection(InitialBracket(:,1), InitialBracket(:,2),         &
    &                      Fa, FB, 1, Delta, L2, Prec)

    N = HFBNumberofParticles(Fermi, Delta, L2 ) - Particles

    !L2 is automatically correct if Lipkin-Nogami is inactive.
    if(.not.Lipkin) exit

    !Temporary save.
    tempL2= L2

    N = HFBNumberofParticles(Fermi, Delta, L2 ) - Particles

    !Evaluate G
    ! print *, 'L2 at start', L2
    flag=0
    G  = LNCr8(Delta, DeltaLN,flag) - L2

    !Check for convergence
    if(all(abs(G).lt.Prec)) then
      exit
    endif
    !Do a secant step
    do it=1,2
      if(flag(it).eq.0) then
        L2 = L2 - (L2 - L2Old)/(G - Gold)*G
      else
        call stp('L2 in bisection failed.')
      endif
    enddo

    !Replace the previous step
    Gold  = G
    L2Old = tempL2

    if(iter.eq.HFBIter) then
      call stp('No correct L2 found.')
    endif
  enddo

  end subroutine HFBLNFindFermiEnergyBisection

subroutine HFBFindFermiEnergyBroyden                                          &
    &         (Fermi,LnLambda,Delta,DeltaLN,Lipkin,DN2,ConstrainDispersion,Prec)
  !-----------------------------------------------------------------------------
  ! This routine varies the Fermi energy and LNLambda parameter to get the
  ! number of particles right, and the Lipkin-Nogami parameter equal to its
  ! calculated value.
  !
  ! This corresponds to the following system of equations
  !
  ! f(Fermi, LNLambda) = Number of Particles(Fermi, LNLambda)-Protons/Neutrons
  ! g(Fermi, LNLambda) = LNLambda(Fermi, LNLambda) - LNLambda guessed
  !
  ! We use a 2-dimensional Broyden method to solve this nonlinear system.
  !
  ! See
  !  http://en.wikipedia.org/wiki/Broyden%27s_method
  ! for some more info.
  !
  !                      ( Fermi_New - Fermi_Old      )
  ! DX =  column vector  (                            )
  !                      (LNLambda_new - LNLambda_old )
  !
  !                      ( f_new - f_old )
  !  F =  column vector  (               )
  !                      ( g_new - g_old )
  !
  ! In short, the algorithm is as follows:
  !
  ! 1) Choose an approximation for the Jacobian matrix J. Here constructed
  !    via some finite difference scheme.
  ! |--
  ! | 2) Find an update direction so that J * DX = - F(Fermi, Lambda)
  ! |   Since the problem is 2 dimensional, we simply invert J.
  ! | 3) Update Fermi energy and LNLambda
  ! | 4) Update the Jacobian:
  ! |     J = J + F( Fermi_New , LNLambda_New) * DX^T / |DX|**2
  ! | 5) Go back to 2 if not converged.
  ! |--
  !-----------------------------------------------------------------------------
  ! Nice to note is that if we do not do the updating of the Jacobian and
  ! instead always calculate a finite-difference Jacobian with a fixed
  ! difference step in Fermi and Lambda, this reduces to the
  ! CR8 method. (With a different damping scheme)
  !
  ! Also, if no Lipkin-Nogami is desired, this routine reduces to the
  ! secant method. ( See http://en.wikipedia.org/wiki/Secant_method )
  ! This method reduces also to the CR8 situation if we fix the discretisation
  ! in Fermi and Lambda.
  !
  !-----------------------------------------------------------------------------
  ! While it might seem like a good idea to let FermiHistory, LNHistory and
  ! the Jacobian persist across mean-field iterations, I have found out that
  ! in practice this is a bad idea. After a mean-field iteration, the Fermi
  ! energy is generally close to the correct one, but the approximation that
  ! was saved might actually be a very bad one.
  !-----------------------------------------------------------------------------
  1 format('Attention, unconverged Fermi solver.'/, &
  &        'Iterations : ', i5,                  /, &
  &        'Particles  : ', 2f15.8)
  2 format('Lipkin parameter deviation: ', 2f15.8)


  real(KIND=dp), intent(inout)              :: Fermi(2), LnLambda(2)
  complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
  complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
  logical, intent(in)                       :: Lipkin, ConstrainDispersion
  real(KIND=dp), intent(in)                 :: Prec, DN2(2)

  integer                                   :: it, iter, i, flag(2)

  ! Previous values of the Fermi energy and LNLambda parameter
  real(KIND=dp), save                   :: FermiHistory(2)=100.0_dp
  real(KIND=dp), save                   :: LNhistory(2)   =100.0_dp
  ! Checking if this is the first time this routine is called
  logical, save                         :: FirstTime=.true.
  !Jacobian matrix, or at least its current approximation and its inverse
  real(KIND=dp), save                   :: Jacobian(2,2,2), invJ(2,2,2)
  real(KIND=dp)                         :: det(2)
  ! Respectively different values of f and g
  real(KIND=dp)                         :: N(2) ,  NX(2), NY(2), NewN(2)
  real(KIND=dp)                         :: LN(2), LNX(2), LNY(2), NewLN(2)
  ! Update values for Fermi and LNLambda
  real(KIND=dp)                         :: FermiUpdate(2), LNUpdate(2),norm(2)
  !Flag for convergence for every isospin
  logical                               :: Converged(2)
  ! Step size
  real(KIND=dp)                         :: step(2) = 0.5_dp
  ! dispersion of N used to check nearness of HF solution
  real(KIND=dp)                         :: disper(2) = 0.0_dp

  !-----------------------------------------------------------------------------
  ! First time, take some guess for the histories
  FermiHistory = Fermi    + 0.01_dp
  LNHistory    = LNLambda + 0.01_dp
  Particles(1) = Neutrons
  Particles(2) = Protons
  Converged    = .false.

  !-----------------------------------------------------------------------------
  ! Block some quasiparticles
  ! don't block qp first time round, this can cause many problems. 
  ! - During the very first call in a program run the U and V cannot be trusted.
  ! - I'm suspicious that during the first call in a later mean-field iteration,
  !   the HFBColumns indexation is the blocked one from the previous iteration,
  !   but I have not checked in detail. In any event, it sometimes happens 
  !   that the blocked HF state cannot be found in neither the U nor V columns.
  !   Commenting the call out makes the blocking more stable.
  !-----------------------------------------------------------------------------
  ! if(allocated(qpexcitations)) then
  !   call      BlockQuasiParticles(Fermi)
  ! endif
  !-----------------------------------------------------------------------------
  ! Store old U and V for future reference
  ! MB 19/01/22 it seems that OldU and OldV are never used in the Broyden solver
  OldU        = U ; OldV = V
  OldColumns  = HFBColumns
  
  flag = 0
  !-----------------------------------------------------------------------------
  !Check were we find ourselves in the phasespace
  N     = HFBNumberofParticles(Fermi, Delta, LnLambda ) - Particles
 
  ! Store the sign of the Pfaffian of the HFB hamiltonian constructed during
  ! initialisation. This pfaffian will be used as future reference to follow the
  ! number parity in each isospin-parity block.
  if (all(oldpf .eq. 0)) then
    oldpf = Pfaffian_HFBHamil()   
  endif

  if(Lipkin) then
    LN   = LNLambda - LNCR8(Delta,DeltaLN, flag)
  else if(ConstrainDispersion) then
    LN   = Dispersion() - DN2
  else
    ! Make sure the algorithm is converged in LN when LN is not active.
    ! This is only to make sure no contamination occurs and convergence
    ! is not wrongly judged.
    LN   = 0.0_dp
  endif

  ! check if particle number is correct already here (as it can happen
  ! in HFB when pairing broke down). This avoids any update of the 
  ! Fermi energy in that case (which otherwise will change its value)
  ! disper = Dispersion() 
  do it=1,2
    if( abs(N(it)).lt.Prec  .and. abs(LN(it)).lt. Prec) Converged(it) = .true.
  enddo 

  ! Evaluating f and g for finite differences.
  NX = HFBNumberOfParticles(FermiHistory, Delta,LNLambda )-Particles
  if(Lipkin) then
    LNX  = LNLambda - LNCR8(Delta,DeltaLN, flag)
  else if(ConstrainDispersion)  then
    LNX  = Dispersion() - DN2
  endif
  if(Lipkin .or. ConstrainDispersion) then
    NY   = HFBNumberOfParticles(Fermi, Delta,LNHistory )-Particles
  endif
  if(Lipkin) then
    LNY  = LNHistory - LNCR8(Delta,DeltaLN, flag)
  else if(ConstrainDispersion) then
    LNY  = Dispersion() - DN2
  endif

  ! Construct a finite difference approximation of the Jacobian
  Jacobian(1,1,:)   = (N   -  NX )/(Fermi    - FermiHistory)
  if(Lipkin .or. ConstrainDispersion) then
    Jacobian(1,2,:) = (N   -  NY )/(LNLambda - LNHistory)
    Jacobian(2,1,:) = (LN  -  LNX)/(Fermi    - FermiHistory)
    Jacobian(2,2,:) = (LN  -  LNY)/(LNLambda - LNHistory)
  endif

  do iter=1,HFBIter

     !Invert the Jacobian
      if (Lipkin .or. ConstrainDispersion) then
        det = Jacobian(1,1,:)*Jacobian(2,2,:) - Jacobian(1,2,:)*Jacobian(2,1,:)
      else
        det = Jacobian(1,1,:)
      endif

      !Invert the jacobian
      if(Lipkin .or. ConstrainDispersion) then
        invJ(1,1,:) =   Jacobian(2,2,:)/det
        invJ(1,2,:) = - Jacobian(1,2,:)/det
        invJ(2,1,:) = - Jacobian(2,1,:)/det
        invJ(2,2,:) =   Jacobian(1,1,:)/det
      else
       ! set a ceiling to inverse Jacobian. This is 
       ! apparently often needed when HFB pairing breaks down
       do it=1,2
         if ( abs(det(it)) .lt. 0.001_dp ) then
           invJ(1,1,it) = sign(0.99_dp,det(it))
         else
           invJ(1,1,it) = 1.0_dp/det(it)
         endif
       enddo
      endif

      !Find a good direction to update in
      do it=1,2
        !Don't update if converged
        if(Converged(it)) cycle
        !Find Update directions
        FermiUpdate(it) = - invJ(1,1,it) * N(it)
        if(Lipkin .or. ConstrainDispersion) then
          FermiUpdate(it) = FermiUpdate(it)          - invJ(1,2,it) * LN(it)
          if(flag(it).ne.1) then
            LNupdate(it)    = - invJ(2,1,it) * N(it) - invJ(2,2,it) * LN(it)
          else
            ! Don't update LN if the calculation of Lambda_2 got into trouble
            !print *, 'Discarded step at iteration ', iter, ' for isospin ', it
            LNUpdate(it)    = 0.0_dp
          endif
          norm(it)        = FermiUpdate(it)**2 + LNUpdate(it)**2
        else
         norm(it)         = FermiUpdate(it)**2
        endif
      enddo

      FermiUpdate = step * FermiUpdate
      LNUpdate    = step * LNUpdate

      !Replace the history and update
      do it=1,2
        ! bugfix: this was not done here before, only elsewhere
        ! such that the Fermi energy continued to change even after
        ! convergence.
        ! Don't update anymore if converged
        if(Converged(it)) cycle
        ! Replace the history
        FermiHistory(it) = Fermi(it)
        Fermi(it)        = Fermi(it)    + FermiUpdate(it)
        if(Lipkin .or. ConstrainDispersion) then
          LNHistory(it) = LNLambda(it)
          LNLambda(it)  = LNLambda(it) + LNUpdate(it)
        endif
      enddo

      !Recalculate
      flag = 0
      N    = HFBNumberOfParticles(Fermi,Delta,LNLambda )  - Particles
      if(Lipkin) then
        LN   = LNLambda - LNCR8(Delta,DeltaLN, flag)
      else if(ConstrainDispersion) then
        LN   = Dispersion() - DN2
      endif

      !Convergence check
      do it=1,2
        if( abs(N(it)).lt.Prec ) Converged(it) = .true.
      enddo
  !   print '(" iter = ",i5," Converged = ",2l," numbers = ",2f16.8,6f12.6)', &
  !     & iter,Converged,(N + Particles),Fermi,invJ(1,1,1),invJ(1,1,2), &
  !     & invJ(1,1,1)*N(1),invJ(1,1,2)*N(2)
      if(all(converged)) then
        exit
      endif
      !Update the Jacobian
      if(iter.gt.1) then
        do it=1,2
          if(Converged(it)) cycle
          Jacobian(1,1,it) = Jacobian(1,1,it)  + N(it)*FermiUpdate(it)/norm(it)
          if(Lipkin .or. ConstrainDispersion) then
            Jacobian(1,2,it) = Jacobian(1,2,it)+ N(it)*LNUpdate(it)   /norm(it)
            Jacobian(2,1,it) = Jacobian(2,1,it)+LN(it)*FermiUpdate(it)/norm(it)
            Jacobian(2,2,it) = Jacobian(2,2,it)+LN(it)*LNUpdate(it)   /norm(it)
          endif
        enddo
      endif

      do it=1,2
        if( N(it).eq.N(it)+1) then
          print '(" it = ",i5)',it
          print '(" Particles = ",1es24.16)',Particles(it)
          print '(" N         = ",1es24.16)',N(it)
          call stp('Nan in the calculation of Fermi energies')
        endif
      enddo
      !Print a warning if not converged
      if(iter.eq.HFBIter) then
          print 1, HFBIter, N + Particles
          print 2, abs(LN)
      endif
  enddo
  ! check history of convergence
  ! because of the logic of the "all" command, we have to check for 
  ! failure, not success
  do i=1,6
    FermiNonConvergenceHistory(i) = FermiNonConvergenceHistory(i+1)
  enddo
  if ( any(abs(N) .gt. 0.9_dp) ) then
    FermiNonConvergenceHistory(7) = .true.
    ! go back to Bogoliubov matrix from the previous SCF iteration and hope
    ! for the best next time round when the s.p.e.s and gaps have changed.
   !do it=1,2
   !  if ( abs(N(it)) .gt. 0.9_dp ) then
   !    U(:,:,:,it)        = OldU(:,:,:,it)
   !    V(:,:,:,it)        = OldV(:,:,:,it)
   !    HFBColumns(:,:,it) = OldColumns(:,:,it)
   !  endif
   !enddo
  else
    FermiNonConvergenceHistory(7) = .false.
  endif
  ! print '(" FermiNonConvergenceHistory",7l)',FermiNonConvergenceHistory
  ! print *,all(FermiNonConvergenceHistory)
  if ( all(FermiNonConvergenceHistory) ) then
  ! print '(/," WARNING: would stop now!",/)'
    call stp(' HFBFindFermiEnergyBroyden: lost particle number for 7 iterations!')
  endif 

end subroutine HFBFindFermiEnergyBroyden

subroutine ConstructRHOHFBLimited(Vlim)
    !------------------------------------------------------------------
    ! Construct RhoHFB from only half of the V matrix
    ! without referencing columns.
    !
    ! Rho = V^* V^T
    !------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Vlim(2*nwt,2*nwt,2,2)
    integer                   :: i,j,k,P,it,N

    do it=1,Iindex
        do P=1,Pindex
            N = blocksizes(P,it)
            do j=1,N
                do i=1,N
                    RhoHFB(i,j,P,it) = 0.0_dp
                    do k=1,N
                        RhoHFB(i,j,P,it) = RhoHFB(i,j,P,it) + Vlim(i,k,P,it) * Vlim(j,k,P,it)
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine ConstructRhoHFBLimited

subroutine ConstructKappaHFBLimited(Ulim,Vlim)
    !------------------------------------------------------------------
    ! Construct KappaHFB from only half of the U and V matrices
    ! without referencing columns.
    !
    ! Kapap = V^* U^T
    !------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Vlim(2*nwt,2*nwt,2,2),Ulim(2*nwt,2*nwt,2,2)
    integer                   :: i,j,k,P,it,N

    do it=1,Iindex
        do P=1,Pindex
            N = blocksizes(P,it)
            do j=1,N
                do i=1,N
                    KappaHFB(i,j,P,it) = 0.0_dp
                    do k=1,N
                        KappaHFB(i,j,P,it) = KappaHFB(i,j,P,it) +       &
                        &                 Vlim(i,k,P,it) * Ulim(j,k,P,it)
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine ConstructKappaHFBLimited

subroutine InitializeUandV(Delta,DeltaLN,Fermi,L2)
    !-------------------------------------------------------------------------------------
    ! Take an initial guess at the U and V matrices, if they have not been read from file.
    ! This routine is mainly useful for retro-active getting some U and V matrices from
    ! older source files.
    !
    !--------------------------------------------------------------------------------------
    real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
    complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
    complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)

    integer :: it, P, i,j,N,ind(2,2)

    1 format ( 58('_'),/,                                             &
    &         'ATTENTION: U AND V HAVE NOT BEEN READ FROM FILE. ',/,  &
    &         'They have been reinitialized and may not be reliable.' &
    &         ,/,58('_'))

    print 1

    call ConstructHFBHamiltonian(Fermi, Delta, L2, HFBGauge)
    call DiagonaliseHFBHamiltonian

    !----------------------------------------------------------------------
    ! Do as in CR8: take only the first half of the matrix, and conjugate
    ! the first quarter.
    if(SC) then
        do it=1,Iindex
            do P=1,Pindex
                N = blocksizes(P,it)
                do i=1,N/2
                    HFBColumns(i,P,it) = 2*N - i + 1
                enddo
                do i=N/2+1,N
                    HFBColumns(i,P,it) = i
                enddo
            enddo
        enddo
    else
        !-----------------------------------------------------------------
        ! If signature is not conserved, just take all the positive
        ! qp energies and pray it will work.
        ind = 1
        do it=1,Iindex
            do P=1,Pindex
                N = blocksizes(P,it)
                do i=1,2*N
                    if(QuasiEnergies(i,P,it) .gt. 0.0_dp) then
                        HFBColumns(ind(P,it),P,it) = i
                        ind(P,it) = ind(P,it) + 1
                    endif
                enddo
            enddo
        enddo
    endif
    call ConstructRHOHFB(HFBColumns)
    call ConstructKappaHFB(HFBColumns)
  end subroutine InitializeUandV

  subroutine FermiBrent(X1,X2,FX1,FX2,Depth,Delta,L2,PrecIn,ZB,ZFB) 
  !-----------------------------------------------------------------------------
  ! For a given L2, this routine searches for the Fermi energy
  ! by Brent's methods https://en.wikipedia.org/wiki/Brent%27s_method
  ! which combines bisection, secant method and inverse quadratic interpolation.
  ! original source is probably
  ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
  ! for Finding a Zero of a Function", Algorithms for Minimization without
  ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall, 
  !-----------------------------------------------------------------------------
  ! This routine is written such that it can replace the call of the recursive
  ! function FermiBisection.
  !-----------------------------------------------------------------------------
  ! see pages 1188 - 1189 of http://apps.nrbook.com/fortran/index.html
  ! [William H. Press, Saul A. Teukolsky, William T. Vetterling and Brian P.
  ! Flannery, Numerical Recipes in Fortran in Fortran 90, Second Edition (1996)]
  !-----------------------------------------------------------------------------
  integer, intent(in)          :: Depth
  real(KIND=dp), intent(in)    :: PrecIn, L2(2)
  real(KIND=dp), intent(in)    :: X1(2) , X2(2) , FX1(2) , FX2(2) 
  real(KIND=dp), intent(inout) :: ZB(2) , ZFB(2)
  real(KIND=dp)                :: A(2) , B(2), C(2) , FA(2), FB(2) , FC(2)
  real(KIND=dp)                :: D(2) , E(2), S(2) , P(2)  , Q(2) , R(2) 
  real(KIND=dp)                :: Num(2) , Prec , Tol(2) , XM(2) 
  real(KIND=dp)                :: eps = 1.d-12
  integer                      :: it , FailCount , i
  logical                      :: Success , Found(2)
  complex(KIND=dp),intent(in), allocatable :: Delta(:,:,:,:)

  Prec = PrecIn  
  A (:) = X1 (:)
  B (:) = X2 (:)
  FA(:) = FX1(:)
  FB(:) = FX2(:)
  
  Found(:) = .false.

  do it=1,Iindex
    if ( A(it) .eq. B(it) ) then 
      !-----------------------------------------------------------------------
      ! This signals that FA(it) = FB(it) is zero within the tolerance.
      ! Either near-converged HFB or HF case of completely broken-down pairing
      ! which also satisfies FA(it) = FB(it) = 0 within an interval. The 
      ! latter case cnnot be handled by the algorithm below.
      !-----------------------------------------------------------------------
      Found(it) = .true. 
    else
      !-----------------------------------------------------------------------
      ! Sanity check: is there a zero of F in [A,B]?
      !-----------------------------------------------------------------------
      if ( FA(it)*FB(it) .gt. 0.0_dp ) then 
        print '(" it = ",i2," A = ",f13.8," FA = ",f14.8,   & 
                       &    " B = ",f13.8," FB = ",f14.8)', &
                       & it,A(it),FA(it),B(it),FB(it)
        call stp(' FermiBrent : function not bracketed!')
      endif
    endif
  enddo
  C (:) = B(:)
  FC(:) = FB(:)

  FailCount = -1  
  Success  = .false.
  do while(.not. Success)
    FailCount = FailCount + 1
    do it=1,Iindex
      !------------------------------------------------------------------------
      ! this isospin is converged already - don't update anything
      !------------------------------------------------------------------------
      if ( Found(it) ) cycle   

      if ( ( FB(it) .gt. 0.0_dp .and. FC(it) .gt. 0.0_dp ) .or. & 
         & ( FB(it) .lt. 0.0_dp .and. FC(it) .lt. 0.0_dp ) )  then
        C (it) = A (it)
        FC(it) = FA(it)
        D(it)  = B(it) - A(it)
        E(it)  = D(it)
      endif
      if ( abs(FC(it)) .lt. abs(FB(it)) ) then
        A (it) = B (it)
        FA(it) = FB(it)
        B (it) = C (it)
        FB(it) = FC(it)
        C (it) = A (it)
        FC(it) = FA(it)
      endif
      !----------------------------------------------------------------
      ! Convergence check
      !----------------------------------------------------------------
!      Tol(it)  = 2.0_dp * eps * abs(B(it)) + 0.5_dp * Prec
      Tol(it) = Prec ! * abs(B(it))
      XM (it)  = 0.5_dp * (C(it)-B(it))

      !----------------------------------------------------------------
      ! Note: the tolerance is on the precision of the Fermi energy,
      ! NOT the nearness of the particle number to the targeted value.
      !----------------------------------------------------------------
      if ( abs(XM(it)) .le. Tol(it) ) then
        ZB(it) =  B(it)
        Found(it) = .true. 
        cycle
      endif
      if ( abs(E(it)) .ge. Tol(it) .and. abs(FA(it)) .gt. abs(FB(it)) ) then
        S = FB(it)/FA(it)
        if ( A(it) .eq. C(it) ) then
          P(it) = 2.0_dp * XM(it) * S(it)
          Q(it) = 1.0_dp - S(it)
        else
          Q(it) = FA(it)/FC(it)
          R(it) = FB(it)/FC(it)
          P(it) = S(it) * (2.0_dp * XM(it) * Q(it) * (Q(it)-R(it)) & 
                  &    - (B(it)-A(it))*(R(it)-1.0_dp))
          Q(it) = (Q(it)-1.0_dp)*(R(it)-1.0_dp)*(S(it)-1.0_dp)
        endif
        if ( P(it) .gt. 0.0_dp ) Q(it) = -Q(it)
        P(it) = abs(P(it))
        if (2.0_dp * P(it) .lt. min(3.0_dp*XM(it)*Q(it) - abs(Tol(it)*Q(it)),abs(E(it)*Q(it)))) then
          E(it) = D(it)
          D(it) = P(it) / Q(it)
        else
          D(it) = XM(it)
          E(it) = D (it)
        endif
      else
        D(it) = XM(it)
        E(it) = D (it)
      endif
      A (it) = B (it)
      FA(it) = FB(it)
      B (it) = B(it) + merge(D(it),sign(Tol(it),XM(it)),abs(D(it)) .gt. Tol(it))
    enddo
    !--------------------------------------------------------------------------
    ! B(:) is present best guess for epsilon_F, FB(:) the corresponding 
    ! particle numbers.
    !--------------------------------------------------------------------------
    Num(:) = HFBNumberofParticles(B, Delta, L2)
    FB (:) = Num(:) - Particles

    !--------------------------------------------------------------------------
    ! check if both isospins are converged, if so, iteration will be left.
    !--------------------------------------------------------------------------
    if ( Found(1) .and. Found(2) ) Success = .true.

    !--------------------------------------------------------------------------
    ! diagnostic printing for convergence analysis (usually commented out)
    !--------------------------------------------------------------------------
    ! For a convergence analysis, uncomment this and the similar statement in
    ! the subroutine HFBFindFermiEnergyBisection() that calls this subroutine.
    !--------------------------------------------------------------------------
    ! NOTE: B is the the best guess for the zero of F. A has been the previous
    ! "closest" interval boundary that is not updated after Found(:)
    ! is set to .true. The actual zero might therefore be outside the 
    ! interval [A,B]. If so, the true zero is typically closer to B than 
    ! A is to B.
    !--------------------------------------------------------------------------
    ! Note further: as A and B are swapped from time to time, B might be 
    ! smaller than A when printed here
    !--------------------------------------------------------------------------
    ! print '(" FermiBrent ",i4,1l2,2(1l2,2(f13.8,es16.7),f14.8))',    &
    !     & FailCount,Success, &
    !     & Found(1),A(1),FA(1),B(1),FB(1),Num(1), & 
    !     & Found(2),A(2),FA(2),B(2),FB(2),Num(2)
    !--------------------------------------------------------------------------
    ! Iteration did not converge for at least one isospin - print warning
    !--------------------------------------------------------------------------
    if ( FailCount .gt. Depth ) then
      print '(/," Warning: FermiBrent not converged after ",i4," iterations")',FailCount
      print *, Tol, Prec
      do it=1,2
        if ( .not. Found(it) ) then
          print '("it = "i2,2(f13.8,es16.7),f14.8)',it,A(it),FA(it),B(it),FB(it),Num(it)
        endif
      enddo
      print '(" ")'
    endif

    !--------------------------------------------------------------------------
    ! The algorithm converges to the eps where N(eps)-N_0 changes sign within 
    ! some tolerance. This does not always mean that also the particle number 
    ! has approached the targeted value N(eps) = N_0, so better check for it.
    ! If it has not, there is possibly the problem that N(eps) - N_0 has a 
    ! discontinuity around the value N_0.
    ! Print warning, print some diagnostics of N(eps)-N_0 around its zero
    ! and cross fingers that the future update of h and Delta makes the 
    ! problem disappear.
    !--------------------------------------------------------------------------
    ! as the handling of blocked states seems to demand to stop iterating
    ! cases where particle number is decreasing with Fermi energy,
    ! this routine now has to handle also cases where FA(it) = FB(it) != 0.0
    !--------------------------------------------------------------------------
    do it=1,2
      if (  Found(it) .and. &
         &  FA(it) .ne. FB(it) .and. &
         &  abs(FB(it)) .gt. 1.d-3 ) then
        print '(/," WARNING: unusual behaviour of N(lambda) for it =",i2)',it
        print '("it = "i2,2(f13.8,es16.7),f14.8)', &
              &   it,A(it),FA(it),B(it),FB(it),Num(it)
        do i=1,21
          if ( A(it) .lt. B(it) ) then
            C(:) = A(:) + (B(:)-A(:)) * (i-5.0_dp) / 10.0_dp 
          else
            C(:) = B(:) + (A(:)-B(:)) * (i-5.0_dp) / 10.0_dp 
          endif
          Num(:) = HFBNumberofParticles(C, Delta, L2) 
          print '(f13.8,f14.8,es16.7)',C(it),Num(it),Num(it)-Particles(it)
        enddo
      endif
    enddo
  enddo
  !----------------------------------------------------------------------------
  ! communicate best guess for Fermi energy and corresponding N.
  !----------------------------------------------------------------------------
  ZB(:)  = B (:)
  ZFB(:) = FB(:) 
  end subroutine FermiBrent

  recursive function FermiBisection(A,B,FA,FB,Depth,Delta,L2,Prec) result(C)
  !-----------------------------------------------------------------------------
  ! For a given L2, this routine recursively searches for the Fermi energy
  ! by a simple bisection method. Very slow, but completely reliable.
  !-----------------------------------------------------------------------------
  integer, intent(in)       :: Depth
  real(KIND=dp), intent(in) :: Prec, L2(2)
  real(KIND=dp), intent(in) :: A(2), B(2), FA(2), FB(2)
  real(KIND=dp)             :: C(2), FC(2), NewA(2), NewB(2), newFA(2), newFB(2)
  integer                   :: it
  logical                   :: Converged(2)
  complex(KIND=dp),intent(in), allocatable :: Delta(:,:,:,:)

  Converged = .false.

  ! Check for maximum number of iterations.
  if(Depth .eq. HFBIter) then
    print *, 'A', A
    print *, 'B', B
    print *, 'FA', FA
    print *, 'FB', FB
    print *, Occupations(:, 1,1)
    print *, Occupations(:, 2,1)
    print *, Occupations(:, 1,2)
    print *, Occupations(:, 2,2)
    print *, 'Warning: Bisection failed.'
    return
  endif

  ! Assign correct values to each C
  do it=1,Iindex
    if(abs(FA(it)).lt.Prec ) then
      C(it) = A(it)
      Converged(it) = .true.
    elseif(abs(FB(it)).lt.Prec ) then
      C(it) = B(it)
      Converged(it) = .true.
    else
      C(it) = (A(it) + B(it))/2.0_dp
    endif
  enddo

  if(all(Converged)) return

  !If not converged, compute the value at the middle of the interval
  FC = HFBNumberofParticles(C, Delta, L2 )-Particles

  !Decide how to continue the recursion
  do it=1,Iindex
    if(FA(it)*FC(it).lt.0.0_dp) then
      newA (it) = A(it)
      newB (it) = C(it)
      newFA(it) = FA(it)
      newFB(it) = FC(it)
    else
      newA (it) = C(it)
      newB (it) = B(it)
      newFA(it) = FC(it)
      newFB(it) = FB(it)
    endif
  enddo
  C = FermiBisection(newA,newB,newFA,newFb,Depth+1,Delta,L2,Prec)

  end function FermiBisection

  subroutine ConstructHFBHamiltonian(lambda, Delta, LnLambda, Gauge)
  !-----------------------------------------------------------------------------
  ! Constructs the HFB hamiltonian.
  ! The argument lambda is the current guess for the Fermi energy.
  !-----------------------------------------------------------------------------
  ! Note that the Rho & Kappa from the past mean-field iteration are used. This
  ! to combat a hysteresis-effect for the Fermi-solver routines.
  !-----------------------------------------------------------------------------
  
  real(KIND=dp), intent(in)                :: lambda(2), LNLambda(2), Gauge(2)
  integer                                  :: i,j, it, ii,iii,P,N,k
  complex(KIND=dp), allocatable,intent(in) :: Delta(:,:,:,:)
   
  !-----------------------------------------------------------------------------
  HFBHamil = 0.0_dp
  do it=1,Iindex
    do P=1,Pindex
      N  = blocksizes(P,it)
      do i=1,N
        ii  = blockindices(i,P,it)
        iii = mod(ii-1,nwt)+1
        !-----------------------------------------------------------------------
        ! Diagonal contribution:
        ! epsilon - lambda  - 2 * Lambda_2 ( 1 - 2 * Rho)
        HFBHamil(i,i,P,it) = HFBasis(iii)%GetEnergy() -  Lambda(it)            &
        &                     - 2*(1.0_dp-LNFraction)* LNLambda(it)
        !------------------------------------------------------------------------
        ! Offdiagonal
        do j=1,N
          HFBHamil(i,j,P,it) = HFBHamil(i,j,P,it)                              &
          &       +  (1.0_dp-LNFraction)*4*LNLambda(it)* oldrhohfb(i,j,P,it)
        enddo
        do j=1,N
          !---------------------------------------------------------------------
          ! - h^*
          HFBHamil(i+N,j+N,P,it)= - conjg(HFBHamil(i,j,P,it))
        enddo
        !-----------------------------------------------------------------------
        do j=1,N
          ! Delta - 2 * Lambda_2 * Kappa
          HFBHamil(i,j + N,P,it) = Delta(i,j,P,it)                             &
          &                    - LNFraction*4*LNLambda(it)*oldkappahfb(i,j,P,it)
        enddo
      enddo
      !-------------------------------------------------------------------------
      !Add the generalized density matrix to the HFBHamiltonian
      !-------------------------------------------------------------------------
      ! ( Rho       Kappa   )
      ! ( -Kappa^*  1-Rho^* )
      !-------------------------------------------------------------------------
!      do i=1,N
!        ! The 1 of the 1-Rho^*
!        k = blockindices(i,P,it)
!        HFBHamil(i+N,i+N,P,it) = HFBHamil(i+N,i+N,P,it) + Gauge(it)
!        do j=1,N
!            ! Rho
!            HFBHamil(i,j,P,it)     = HFBHamil(i,j,P,it)     + Gauge(it)*OldRhoHFB(i,j,P,it)
!            ! Kappa
!            HFBHamil(i,j+N,P,it)   = HFBHamil(i,j+N,P,it)   + Gauge(it)*OldKappaHFB(i,j,P,it)
!            HFBHamil(i+N,j,P,it)   = HFBHamil(i+N,j,P,it)   - Gauge(it)*OldKappaHFB(i,j,P,it)
!            ! - Rho
!            HFBHamil(i+N,j+N,P,it) = HFBHamil(i+N,j+N,P,it) - Gauge(it)*Conjg(OldRhoHFB(i,j,P,it))
!        enddo
!      enddo
      !--------------------------------------------------------------------------
      ! Make sure everything is symmetric.
      ! This might strictly be a waste of CPU cycles, but anyone carelessly using
      ! the Hamiltonian will then not be confronted with nasty surprises.
      do i=1,2*N
        do j=i+1,2*N
            HFBHamil(j,i,P,it) = conjg(HFBHamil(i,j,P,it))
        enddo
      enddo
    enddo
  enddo

  if(all(HFBHamil.eq.0.0_dp)) call stp('HFBHamiltonian completely zero!')
  end subroutine ConstructHFBHamiltonian
  
  function Pfaffian_HFBHamil() result(pf)
    !---------------------------------------------------------------------------
    ! Calculate the pfaffian of the Hamiltonian in the Majorana 
    ! representation. 
    !---------------------------------------------------------------------------
    ! Currently only valid when time-simplex is conserved.
    !---------------------------------------------------------------------------

    use pfaff

    real(KIND=dp)                            :: pf(2,2)
    integer ,allocatable                     :: ipv(:,:)
    real(KIND=dp), allocatable               :: Hamcopy(:,:)
    integer                                  :: it, P, N

    !---------------------------------------------------------------------------
    if(.not.allocated(ipv)) then
      N = maxval(blocksizes)
      allocate(ipv(2*N,2))
    endif
    if(.not.allocated(Hamcopy)) then
      N = maxval(blocksizes)
      allocate(Hamcopy(2*N,2*N))
    endif
    !---------------------------------------------------------------------------
    pf = 0.0
    do it=1,Iindex
      do P=1,Pindex
      
        N = blocksizes(P,it)
      
        Hamcopy(1:2*N, 1:2*N) = & 
        &              MajoranaTransform(real(HFBHamil(1:2*N, 1:2*N, P,it))) 
      
        call PfaffianF(Hamcopy(1:2*N,1:2*N),2*N,2*N,ipv(1:2*N,:),pf(p,it)) 
      enddo
    enddo
  
  end function Pfaffian_HFBHamil
 
  function MajoranaTransform( Hin ) result(Maj)
    !---------------------------------------------------------------------------
    ! Transform a HFB hamiltonian into a Majorana representation.
    !-----------------------------------------------------------------------
    ! Currently only valid when time-simplex is conserved. 
    !---------------------------------------------------------------------------
    real(KIND=dp) :: Hin(:,:)
    real(KIND=dp) :: Maj(size(Hin,1),size(Hin,2))
    integer          :: N
    
    ! Trivial to do when S^T_y is conserved.
    if(.not. TSC) call stp('Majorana transform invalid when timesimplex is broken.')
  
    N = size(Hin,1)/2
    
    Maj(1:N    ,1:N    )     = 0
    Maj(N+1:2*N,N+1:2*N)     = 0
    Maj(1:N    ,N+1:2*N) =   Hin(1:N,1:N) - Hin(1:N,N+1:2*N) 
    Maj(N+1:2*N,1:N)     = - Maj(1:N    ,N+1:2*N)
  
  end function MajoranaTransform

  subroutine DiagonaliseHFBHamiltonian_Signature
    !-------------------------------------------------------------------------------
    ! Diagonalise the HFBHamiltonian to get the U and V eigenvectors when Signature
    ! is conserved. That signature is conserved implies that we can split the
    ! HFBhamiltonian into two parts. When Time-reversal is conserved we can even get
    ! away with diagonalising only one half.
    !-------------------------------------------------------------------------------
    integer                         :: i,j,k,it, P,jj,ii,jjj,S,Sig1,iii,sig2, Nmax,N
    integer                         :: ifail
    real(KIND=dp), allocatable,save ::  Eigenvectors(:,:),Eigenvalues(:)
    real(KIND=dp), allocatable,save :: WORK(:)
    real(KIND=dp), allocatable,save :: Temp(:,:)
    real(KIND=dp)                   :: test

    Nmax = maxval(blocksizes)
    !---------------------------------------------------------------------------
    ! Preliminary work.
    if(.not.allocated(WORK)) then
        allocate(Temp(Nmax,Nmax))         ; Temp = 0.0_dp
        allocate(Work(Nmax))              ; Work = 0.0_dp
        allocate(Eigenvectors(Nmax,Nmax)) ; Eigenvectors=0.0_dp
        allocate(Eigenvalues(Nmax))       ; Eigenvalues=0.0_dp
    endif

    U = 0.0_dp
    V = 0.0_dp
    do it=1,Iindex
        do P=1,Pindex
            N = blocksizes(P,it)

            Temp(1:N,1:N)        = 0.0_dp
            Eigenvalues(1:N)     = 0.0_dp
            Eigenvectors(1:N,1:N)= 0.0_dp
            Work(1:N)            = 0.0_dp
            !-----------------------------------------------------------------------
            ! Copy the time-reversed positive signature part to Temp
            do j=1,N/2
                do i=1,N/2
                    Temp(i,j)         = DBLE(HFBHamil(i,j,P,it))
                    Temp(j,i)         = Temp(i,j)
                    Temp(i+N/2,j+N/2) = DBLE(HFBHamil(i+3*N/2,j+3*N/2,P,it))
                    Temp(j+N/2,i+N/2) = Temp(i+N/2,j+N/2)
                    Temp(i,j+N/2)     = DBLE(HFBHamil(i,j+3*N/2,P,it))
                    Temp(j+N/2,i)     = Temp(i,j+N/2)
                enddo
            enddo
            !-------------------------------------------------------------------------
            ! Diagonalize
            call diagoncr8(temp,Nmax,N, Eigenvectors, Eigenvalues, Work, 'DiagHamil ',ifail)
            if(ifail.eq.1) then
              do i=1,N
                print *, Temp(i,1:N)
              enddo

              print *
              print *, 'Parity', P, 'Isospin', it
              do i=1,N
                print *, DBLE(HFBHamil(i,1:2*N,P,it))
              enddo
              call stp('Diagoncr8 failed')
            endif
            U(      1:N/2,    1:N  ,P,it) = Eigenvectors(       1:N/2   ,  1:N  )
            V(  N/2+1:N  ,    1:N  ,P,it) = Eigenvectors(   N/2+1:N     ,  1:N  )
            QuasiEnergies(1:N,P,it)       = EigenValues (1:N)

            if(TRC) then
                !------------------------------------------------------------
                ! We can get the other half of the eigenvectors by symmetry
                U(  N/2+1:N  ,N+1:2*N  ,P,it) =-Eigenvectors(       1:N/2   ,  1:N  )
                V(      1:N/2,N+1:2*N  ,P,it) = Eigenvectors(   N/2+1:N     ,  1:N  )
                QuasiEnergies(N+1:2*N,  P,it) = EigenValues(1:N)
            else
                Temp(1:N,1:N)        = 0.0_dp
                Eigenvalues(1:N)     = 0.0_dp
                Eigenvectors(1:N,1:N)= 0.0_dp
                Work(1:N)            = 0.0_dp
                do j=N+1,2*N
                    do i=N+1,2*N
                        Temp(i-N,j-N) = DBLE(HFBHamil(i-N/2,j-N/2,P,it))
                        Temp(j-N,i-N) = Temp(i-N,j-N)
                    enddo
                enddo
                call diagoncr8(temp,Nmax,N, Eigenvectors, Eigenvalues, Work, 'DiagHamil2',ifail)
                if(ifail.eq.1) call stp('2nd')
                U(  N/2+1:N  ,N+1:2*N  ,P,it) = Eigenvectors(    1:N/2 ,1:N)
                V(      1:N/2,N+1:2*N  ,P,it) = Eigenvectors(N/2+1:N   ,1:N)
                QuasiEnergies(N+1:2*N  ,P,it) = EigenValues(     1:N)
            endif
        enddo
    enddo
  
end subroutine DiagonaliseHFBHamiltonian_Signature

subroutine DiagonaliseHFBHamiltonian_NoSignature
    !-------------------------------------------------------------------------------
    ! Diagonalise the HFBHamiltonian to get the U and V eigenvectors when Signature
    ! is not conserved.
    !-------------------------------------------------------------------------------
    integer                         :: i,j,k,it, P,jj,ii,jjj,S,Sig1,iii,sig2, Nmax,N
    integer                         :: ifail
    real(KIND=dp), allocatable,save :: Eigenvectors(:,:),Eigenvalues(:)
    real(KIND=dp), allocatable,save :: WORK(:)
    real(KIND=dp), allocatable,save :: Temp(:,:)
    real(KIND=dp)                   :: test

    Nmax = 2*maxval(blocksizes)
    !---------------------------------------------------------------------------
    ! Preliminary work.
    if(.not.allocated(WORK)) then
        allocate(Temp(Nmax,Nmax))         ; Temp = 0.0_dp
        allocate(Work(Nmax))              ; Work = 0.0_dp
        allocate(Eigenvectors(Nmax,Nmax)) ; Eigenvectors=0.0_dp
        allocate(Eigenvalues(Nmax))       ; Eigenvalues=0.0_dp
    endif

    U = 0.0_dp
    V = 0.0_dp
    do it=1,Iindex
        do P=1,Pindex
            N = 2*blocksizes(P,it)
            Temp(1:N,1:N)        = 0.0_dp
            Eigenvalues(1:N)     = 0.0_dp
            Eigenvectors(1:N,1:N)= 0.0_dp
            Work(1:N)            = 0.0_dp

            !-----------------------------------------------------------------------
            ! Copy HFBHamiltonian to Temp
            do j=1,N
                do i=1,j
                    Temp(i,j) = DBLE(HFBHamil(i,j,P,it))
                    Temp(j,i) = Temp(i,j)
                enddo
            enddo
            !-------------------------------------------------------------------------
            ! Diagonalize
            call diagoncr8(temp,Nmax,N, Eigenvectors, Eigenvalues, Work, 'DiagHamil ',ifail)
            U(      1:N/2, 1:N ,P,it) = Eigenvectors(      1:N/2,  1:N  )
            V(      1:N/2, 1:N ,P,it) = Eigenvectors(  N/2+1:N  ,  1:N  )
            QuasiEnergies(1:N,P,it)   = EigenValues(1:N)
        enddo
    enddo
end subroutine DiagonaliseHFBHamiltonian_NoSignature

!  subroutine DiagonaliseHFBHamiltonian_ZHEEVR
!    !---------------------------------------------------------------------------
!    ! Alternative subroutine for diagonalising the HFB Hamiltonian.
!    !---------------------------------------------------------------------------
!    ! Subroutine that solves the HFB eigenvalue problem
!    !
!    !      ( U_k )        ( U_k )
!    !   H  (     )  = Ek  (     )
!    !      ( V_k )        ( V_k )
!    !
!    ! And stores the result in matrices U & V.
!    !---------------------------------------------------------------------------
!    ! NOTE
!    ! 1) Extra documentation for the diagonalisation routine
!    !    http://www.netlib.org/lapack/explore-html/d9/dd2/zheevr_8f.html
!    !---------------------------------------------------------------------------
!    integer, save                       :: size, m, lwork,lrwork,liwork
!    complex(KIND=dp), allocatable,save       :: Eigenvectors(:,:)
!    complex(KIND=dp), allocatable,save       :: Temp(:,:)
!    real(KIND=dp) , allocatable,save         :: Eigenvalues(:)
!    complex(KIND=dp), allocatable,save       :: Work(:), rwork(:)
!    integer, allocatable,save                :: iwork(:), isuppz(:)
!    integer                             :: Succes,i,j,jj,it,P, N, jjj, S
!    real(KIND=dp)                       :: Sig
!    real(KIND=dp)                       :: time0=0,time1=0,timetot=0

!    size = maxval(blocksizes)
!    !---------------------------------------------------------------------------
!    ! Preliminary work.
!    if(.not.allocated(WORK)) then
!        !Allocate several arrays for the LAPACK routine
!        allocate(Temp(2*size, 2*size))                   ; Temp = 0.0_dp
!        allocate(Eigenvectors(2*HFBSize, 2*HFBsize))     ; Eigenvectors= 0.0_dp
!        allocate(Eigenvalues(2*size))                    ; EigenValues=0.0_dp
!        allocate(Isuppz(4*size))                         ; ISUPPZ=0
!        ! Do a preliminary call to the LAPACK routine to find the optimum
!        ! WORK sizes.
!        LWORK = -1 ; LRWORK = -1 ; LIWORK = -1
!        allocate(work(1), rwork(1), iwork(1))
!        call ZHEEVR('V', 'A', 'U', 2*size, Temp,                            &
!        &            2*size, 0.0_dp , 0.0_dp,0,0,                           &
!        &            0.0_dp, 2*size, Eigenvalues, Eigenvectors, 2*HFBSize,  &
!        &            isuppz,work,lwork,rwork,                               &
!        &            lrwork, iwork, liwork, Succes)

!        lwork = ceiling(DBLE(WORK(1))) ; lrwork = ceiling(DBLE(RWORK(1)))
!        liwork =IWORK(1)
!        deallocate(WORK, RWORK, IWORK)
!        allocate(WORK(lwork), RWORK(lrwork), IWORK(liwork))
!    endif
!    !----------------------------------------------------------------------------
!    U = 0.0_dp ; V=0.0_dp
!    do it=1,Iindex
!      do P=1,Pindex
!        Temp = 0.0_dp
!        N = blocksizes(P,it)
!        do j=1,2*N
!          do i=1,2*N
!            Temp(i,j) = HFBHamil(i,j,P,it)
!          enddo
!        enddo
!        !------------------------------------------------------------------------
!        ! Temporary trick to separate signatures.

!        if(SC) then
!          do j=1,N
!            jj = blockindices(j,P,it)
!            jjj= mod(jj-1,nwt)+1
!            S = HFBasis(jjj)%GetSignature()
!            if(jjj .ne. jj) S = -S
!            Temp(j,j)     = Temp(j,j)     + 1000*S
!            Temp(j+N,j+N) = Temp(j+N,j+N) - 1000*S
!          enddo
!        endif
!        !------------------------------------------------------------------------

!        Succes = 0
!        m = N
!        call ZHEEVR('V', 'A', 'U', 2*N , Temp(1:2*N,1:2*N),2*N, 0.0_dp,          &
!          &          0.0_dp,0,0,0.0_dp, 2*N, Eigenvalues(1:2*N),                 &
!          &          Eigenvectors(1:2*N,1:2*N),2*N, isuppz, work, lwork, rwork,  &
!          &          lrwork,iwork, liwork,Succes)
!        if(Succes.ne.0) then
!          call stp("Error in diagonalising the HFB Hamiltonian.",                &
!          &        "ZHEEVR Errorcode", Succes)
!        endif
!        !-------------------------------------------------------------------------
!        ! We store all possible eigenvectors and later make the proper selection.
!        U(1:N,1:2*N,P,it) = Eigenvectors(  1:N    ,1:2*N)
!        V(1:N,1:2*N,P,it) = Eigenvectors(  N+1:2*N,1:2*N)
!        QuasiEnergies(1:2*N,P,it) = EigenValues(1:2*N)

!        !-------------------------------------------------------------------------
!        ! Calculate the quasiparticle signatures for further use
!        QuasiSignatures(:,P,it)=0.0_dp
!        do i=1,2*N
!          do j=1,N
!            jj  = blockindices(j,P,it)
!            jjj = mod(jj-1,nwt)+1
!            Sig = HFBasis(jjj)%GetSignatureR()
!            if(jj.ne.jjj) Sig = -Sig
!            !-----------------------------------------------------
!            ! The U components have the same signature
!            QuasiSignatures(i,P,it) = QuasiSignatures(i,P,it)    &
!            &                       + sig * abs(U(j,i,P,it))**2
!            !-----------------------------------------------------
!            ! The V components have opposite signature, thus
!            ! opposite sign.
!            QuasiSignatures(i,P,it) = QuasiSignatures(i,P,it)    &
!            &                       - sig * abs(V(j,i,P,it))**2
!          enddo
!        enddo
!        !------------------------------------------------------------------------
!        ! Extra for temporary trick
!        if(SC) then
!          QuasiEnergies(1:2*N, P, it) = QuasiEnergies(1:2*N, P, it) &
!          &                         - 1000*QuasiSignatures(1:2*N,P,it)
!        endif
!      enddo
!    enddo

!    if(SC) call InsertionSortQPEnergies()

!end subroutine DiagonaliseHFBHamiltonian_ZHEEVR

subroutine InsertionSortQPEnergies
    !-----------------------------------------------------------------------
    ! Sort the QPenergies, since the trick to separate in
    ! subroutine DiagonaliseHFBHamilaltonian does not conserve the order.
    ! We also sort the U and V matrices, as well as the Quasisignatures.
    !-----------------------------------------------------------------------
    integer :: i,j, it, P, N
    real(KIND=dp) :: TempE,TempSig
    complex(KIND=dp) :: Temp(HFBSize)

    do it=1,Iindex
        do P=1,Pindex
            N = blocksizes(P,it)
            do i=1,2*N
              j = i
              if(j.le.1) cycle
              !----------------------------------------------------------
              ! This condition might seem weird, but it ensures that the
              ! positive signature state ends up at a higher index in
              ! the array, since these states are shifted to higher
              ! qp-energies in the diagonalization routine.
              ! This of course fails when signature is not conserved.
              do while(QuasiEnergies(j-1,P,it)-QuasiEnergies(j,P,it) .gt. 1d-8)
                  TempE = QuasiEnergies(j,P,it)
                  QuasiEnergies(j,P,it) = QuasiEnergies(j-1,P,it)
                  QuasiEnergies(j-1,P,it) = TempE

                  ! Also change the U and V around
                  Temp(1:N)       = U(1:N,j,P,it)
                  U(1:N,j,P,it)   = U(1:N,j-1,P,it)
                  U(1:N,j-1,P,it) = Temp(1:N)

                  Temp(1:N)       = V(1:N,j,P,it)
                  V(1:N,j,P,it)   = V(1:N,j-1,P,it)
                  V(1:N,j-1,P,it) = Temp(1:N)

                  ! And don't forget the quasisignatures!
                  TempSig                   = QuasiSignatures(j,P,it)
                  Quasisignatures(j,P,it)   = QuasiSignatures(j-1,P,it)
                  Quasisignatures(j-1,P,it) = TempSig

                  j = j - 1
                  !Stop when we've reached the bottom of the list
                  if(j.eq.1) exit
              end do

          enddo
        enddo
    enddo

  end subroutine InsertionSortQPEnergies

  subroutine CheckUandVColumns(Columns)
    !---------------------------------------------------------------------------
    ! Routine that checks if the U's and V's obey relations 7.5 in R&S.
    ! Namely:
    !        U^{dagger} U + V^{\dagger}V  = 1  [Check(1)]
    !        U^{T}V       + V^{T} U       = 0  [Check(2)]
    !        U U^{\dagger}+ V^* V^T       = 1  [Check(3)]
    !        U V^{\dagger}+ V^* U^T       = 0  [Check(4)]
    !--------------------------------------------------------------------------
    integer             :: i, j, k, P, it, ii, jj, kk
    complex(KIND=dp)    :: Check(4), UdaggerU, VDaggerV,UTV,VTU, UUdagger
    complex(KIND=dp)    :: VstarVT, ref, UVDagger, VstarUT
    logical             :: Problem, Cont
    integer, intent(in) :: Columns(HFBsize,Pindex,Iindex)
    real(KIND=dp)       :: Prec = 1d-6

    if(any(Imag(U) .ne. 0.0_dp)) call stp('U imaginary')
    print *, 'Check'
    Problem = .false.
    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          jj = Columns(j,P,it)
          do i=1,blocksizes(P,it)
            ii = Columns(i,P,it)
            Check = 0.0_dp
            do k=1,blocksizes(P,it)
              kk = Columns(k,P,it)
              UdaggerU = conjg(U(k,ii,P,it))*U(k,jj,P,it)
              VdaggerV = conjg(V(k,ii,P,it))*V(k,jj,P,it)

              UTV      = U(k,ii,P,it)*V(k,jj,P,it)
              VTU      = V(k,ii,P,it)*U(k,jj,P,it)

              UUdagger = U(i,kk,P,it)*conjg(U(j,kk,P,it))
              VstarVT  = conjg(V(i,kk,P,it))*V(j,kk,P,it)

              UVdagger = U(i,kk,P,it)*conjg(V(j,kk,P,it))
              VstarUT  = conjg(V(i,kk,P,it))*U(j,kk,P,it)

              Check(1) = Check(1) + UdaggerU + VdaggerV
              Check(2) = Check(2) + UTV      + VTU
              Check(3) = Check(3) + UUdagger + VstarVT
              Check(4) = Check(4) + UVdagger + VstarUT
            enddo

            if(i.eq.j) then
              ref = 1.0_dp
            else
              ref = 0.0_dp
            endif

            if(abs(Check(1) - ref).gt.Prec) then
              print *, 'Check 1 failed at: ', ii,jj, Check(1),' P=',P,' it=', it
              Problem=.true.
            endif
            if(abs(Check(2)      ).gt.Prec) then
              print *, 'Check 2 failed at: ', ii,jj, Check(2),' P=',P,' it=', it
              Problem=.true.
            endif
            if(abs(Check(3) - ref).gt.Prec) then
              print *, 'Check 3 failed at: ', ii,jj, Check(3),' P=',P,' it=', it
              Problem=.true.
            endif
            if(abs(Check(4)      ).gt.Prec) then
              print *, 'Check 4 failed at: ', ii,jj, Check(4),' P=',P,' it=', it
              Problem=.true.
            endif
          enddo
        enddo
      enddo
    enddo
    if(Problem) then
        do it=1,Iindex
            do P=1,Pindex
                print *,' P = ', P, 'It = ', it
                print *
                do i=1,2*blocksizes(P,it)
                    write(*, "(i7)", ADVANCE='NO') i
                enddo
                print *
                do i=1,2*blocksizes(P,it)
                   Cont = .false.
                   do j=1,blocksizes(P,it)
                        if(HFBColumns(j,P,it) .eq. i) then
                            Cont=.true.
                            exit
                        endif
                   enddo
                   if(Cont) then
                    write(*, "(i7)", ADVANCE='NO') 1
                   else
                    write(*, "(i7)", ADVANCE='NO') 0
                   endif
                enddo
                print *

                do i=1,blocksizes(P,it)
                    write(*,"(99f7.3)"), DBLE(U(i,1:2*blocksizes(P,it),P,it))
                enddo
                print *
                 do i=1,blocksizes(P,it)
                    write(*,"(99f7.3)"), DBLE(V(i,1:2*blocksizes(P,it),P,it))
                enddo
                print *
                write(*, '(99f7.3)') QuasiSignatures(1:2*blocksizes(P,it),P,it)
                write(*, '(99f7.3)') QuasiEnergies(1:2*blocksizes(P,it),P,it)
            enddo
            enddo
        call stp('Checks failed for U and V.')
    endif
  end subroutine CheckUandVColumns

  subroutine constructRhoHFB(Columns)
  !-----------------------------------------------------------------------------
  ! This function constructs a density matrix from the input of columns from
  ! the U and V matrices.
  !
  ! rho   = V^* V^T
  !-----------------------------------------------------------------------------
    integer, intent(in) :: Columns(HFBSize,Pindex,Iindex)
    integer             :: i,j,k,it,P
    real :: s
    !---------------------------------------------------------------------------
    ! Actual computation
    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          do i=1,blocksizes(P,it)
            RhoHFB(i,j,P,it) = 0.0_dp
            do k=1,blocksizes(P,it)
              RhoHFB(i,j,P,it) = RhoHFB(i,j,P,it) +                            &
              &       Conjg(V(i,Columns(k,P,it),P,it))*V(j,Columns(k,P,it),P,it)
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine ConstructRHOHFB

  subroutine DiagonaliseRhoHFB!_diagoncr8
        !-----------------------------------------------------------------------
        ! Diagonalise the HFB density matrix using the diagon subroutine from cr8.
        !
        ! Notice that we do not bother to explicitly use the signature symmetry
        ! to simplify the diagonalization. This routine does not cost any time in
        ! practice, since the dimension of Rho is only half that of he HFB
        ! hamiltonian.
        !
        !-----------------------------------------------------------------------
        ! Note that we do not diagonalise the full HFB density. We first
        ! identify the subspace of the HFBasis where the pairing is active.
        ! 
        !-----------------------------------------------------------------------
        ! NOT SUITED FOR TIME SIMPLEX BREAKING CALCULATIONS!
        !-----------------------------------------------------------------------
        real(KIND=dp), allocatable :: Work(:), Eigenvalues(:)
        real(KIND=dp), allocatable :: Temp(:,:)
        real(KIND=dp), allocatable :: Eigenvectors(:,:)
        integer                    :: P, it, N, i, Nmax, ii, iii, ifail

        Nmax = maxval(Blocksizes)
         if(.not. allocated(Work)) then
            allocate(Eigenvectors(Nmax,Nmax)) ; Eigenvectors = 0.0_dp
            allocate(Work(Nmax*5))            ; Work = 0.0_dp
            allocate(Temp(Nmax,Nmax))         ; Temp = 0.0_dp
            allocate(Eigenvalues(Nmax))       ; Eigenvalues = 0.0_dp
        endif
        
        do it=1,Iindex

          do P=1,Pindex
!              !---------------------------------------------------------------
!              ! Identify the pairing subspace 
!              do j=1,N
!                if(any(abs(abs(RhoHFB(1:N,j,P,it))-1).lt.HFBNumcut)) then
!                    ! This level corresponds to a fully occupied HF level
!                    
!                
!                endif
!              enddo
              
              Temp          = 0.0_dp
              N             = blocksizes(P,it)
              Temp(1:N,1:N) = DBLE(RhoHFB(1:N,1:N,P,it))

              !---------------------------------------------------------------
              ! when the particle number of species "it" is zero, the canonical
              ! basis cannot be safely obtained by diagonalizing rho, as this
              ! matrix is formally zero and numerically numerical noise.
              ! Needed for trapped neutron droplets.
              !---------------------------------------------------------------
              if ( Particles(it) .lt. 0.000000001_dp ) then
                ! print '(" Particles ",i5,f16.8)',it,Particles(it)
                CanTransfo (1:N,1:N,P,it) = 0.0_dp
                do i=1,N
                  CanTransfo (i,i,P,it) = 1.0_dp
                enddo
                Occupations(1:N,P,it)     = 0.0_dp
                cycle
              endif

              !-----------------------------------------------------------------
              ! Since our matrices are not split into signature blocks, and in
              ! general the eigenvalues of the density matrix are degenerate,
              ! we need to ensure that states of different signature do not get
              ! mixed. We do this by artificially shifting the eigenvalues of
              ! the negative signatures.
              !
              ! Notice that, if V is some eigenvector of some matrix A with
              ! eigenvalue lambda, then V is also an eigenvalue of (A - b*I)
              ! where I is the identity matrix and b is a number.
              ! Its corresponding eigenvalue is then just shifted, but the
              ! the eigenvectors of A and (A - b*I) are the same.
              !
              ! In this way, we subtract two from the negative signature
              ! diagonal elements, thus putting the negative signature
              ! eigenvalues in the (-2,-1) range, well separated from the
              ! positive signature eigenvalues in the (0,1) range.
              !-----------------------------------------------------------------
              if(SC) then
                do i=1,N
                  ii  = blockindices(i,P,it)
                  iii = mod(ii-1,nwt)+1
                  if( ii .ne. iii .or. HFBasis(iii)%GetSignature().eq.-1) then
                    Temp(i,i) = Temp(i,i) - 2.0_dp
                  endif
                enddo
              endif
              call diagoncr8(temp,Nmax,N,Eigenvectors,Eigenvalues, Work,       &
              &                                              'DiagRho   ',ifail)

              CanTransfo (1:N,1:N,P,it) = Eigenvectors(1:N,1:N)
              Occupations(1:N,P,it)     = Eigenvalues(1:N)
              !-----------------------------------------------------------------
              ! Since we shifted the eigenvalues of the negative signature states
              ! by -2, we now need to find the actual occupation numbers.
              ! Notice the slight offset to make sure we are not accidentally
              ! adding two to positive signature eigenvalues
              if(SC) where( Occupations.lt.-0.1_dp) Occupations = Occupations + 2
              !-----------------------------------------------------------------
              ! Ensure that occupations are positive semi-definite
              !-----------------------------------------------------------------
              ! MB 19/02/26 in the weak-pairing limit, occupations might 
              ! become -eps because of the accumulation of eigenvalues around
              ! 0, which might have catastrophic consequences when this 
              ! concerns a continuum state that is peaked at large distances,
              ! as this might lead to negative densities rho(r) at large 
              ! distances that are negative with absolute values larger than
              ! the 1.d-20 shift of rho when calculating rho^alpha(r). For
              ! such cases the calculation of rho^alpha(r) produces a NaN
              ! somewhere on the mesh, which then makes the integral a NaN,
              ! and also potentials a NaN.I had many calculations where this 
              ! happened in magic nuclei in HFB sans LN.
          !   do i=1,N
          !     if ( Occupations(i,P,it) .lt. 1d-11 ) &
          !     &  print '(" DiagonaliseRhoHFB!_diagoncr8 1 : ",3i5,1es16.8)', &
          !     &  it,P,i,Occupations(i,P,it)
          !   enddo
              where (  Occupations .lt. 0.0_dp ) Occupations = 0.0_dp
          !   do i=1,N
          !     if ( Occupations(i,P,it) .lt. 1d-11 ) &
          !     &  print '(" DiagonaliseRhoHFB!_diagoncr8 2 : ",3i5,1es16.8)', &
          !     &  it,P,i,Occupations(i,P,it)
          !   enddo
          enddo
        enddo
        where (abs(CanTransfo) .lt. 1d-11) CanTransfo = 0.0_dp
  end subroutine DiagonaliseRhoHFB!_diagoncr8

  subroutine CheckRho(Rho)
    !--------------------------------------------------------------------------
    ! Subroutine to check a density matrix(doesnt have to be the current one)
    ! to see if it satisfies a number of things.
    ! 1) Check if it is Hermitian
    ! 2) Check that it doesn't couple different parities.
    ! 3) Check that it doesn't couple protons & neutrons.
    !--------------------------------------------------------------------------
    complex(KIND=dp), intent(in) :: Rho(HFBSize,HFBSize,Pindex,Iindex)
    integer                      :: i,ii,j,jj,it,P, S1, S2, iii, jjj

    !---------------------------------------------------------------------------
    !Hermeticity
    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          do i=1,j
            if(abs(Rho(i,j,P,it) - conjg(Rho(j,i,P,it))) .gt. 1d-5) then
              print *, 'Rho is not hermitian ',                                &
              &                          i,j, Rho(i,j,P,it),conjg(Rho(j,i,P,it))
            endif
          enddo
        enddo
      enddo
    enddo
    !---------------------------------------------------------------------------
    ! Rho does not couple signatures
    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          jj = blockindices(j,P,it)
          jjj= mod(jj-1,nwt)+1
          S1 = HFBasis(jjj)%GetSignature()
          if(jjj.ne.jj) S1 = -S1
          do i=1,blocksizes(P,it)
            ii = blockindices(i,P,it)
            iii= mod(ii-1,nwt)+1
            S2 = HFBasis(iii)%GetSignature()
            if(iii.ne.ii) S2 = -S2
            if(abs(Rho(i,j,P,it)) .gt. HFBNumCut .and. S1.ne.S2) then
              print *, i,j,S1,S2, Rho(i,j,P,it)
              call stp('Rho couples signatures')
            endif
          enddo
        enddo
      enddo
    enddo
  end subroutine CheckRho

  subroutine ConstructKappaHFB(Columns)
  !----------------------------------------------------------------------------
  ! This function constructs the anomalous density matrix from the input of
  ! columns from the U and V matrices.
  !
  ! Kappa   = V^* U^T
  !----------------------------------------------------------------------------
    integer, intent(in) :: Columns(HFBSize,Pindex,Iindex)
    integer             :: i,ii, j, jj, it, P, k, S1, S2, iii, jjj

    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          do i=j+1,blocksizes(P,it)
            KappaHFB(i,j,P,it) = 0.0_dp
            do k=1,blocksizes(P,it)
              KappaHFB(i,j,P,it) = KappaHFB(i,j,P,it) +                    &
              &   conjg(V(i,Columns(k,P,it),P,it))*U(j,Columns(k,P,it),P,it)
            enddo
            KappaHFB(j,i,P,it) = - KappaHFB(i,j,P,it)
          enddo
        enddo
      enddo
    enddo
  end subroutine ConstructKappaHFB

  subroutine HFBOccupations(Fermi, Delta,LNLambda,PairingDisp)
  !-----------------------------------------------------------------------------
  ! This routine diagonalises the HFB density matrix, thereby (hopefully) also
  ! breing Kappa into canonical form.
  ! Once diagonalised, the canonical basis is constructed as linear combination
  ! of the HF-basis states.
  !-----------------------------------------------------------------------------
  use wavefunctions

  real(KIND=dp), intent(inout)              :: PairingDisp(2)
  real(KIND=dp), intent(in)                 :: Fermi(2), LNLambda(2)
  complex(KIND=dp), intent(in), allocatable :: Delta(:,:,:,:)
  real(KIND=dp)                             :: mx
  real(KIND=dp)                             :: Energy,  RhoII, SR, overlaps(nwt,nwt)
  integer                                   :: it,P,i,j,S,ii,iii,loc(1),TS,jj,jjj,k
  integer                                   :: Columns(nwt,Pindex,Iindex)
  integer                                   :: P2, C, index, N

  PairingDisp = 0.0_dp

  !-----------------------------------------------------------------------------
  ! Actual diagonalisation
  call DiagonaliseRHOHFB
  !-----------------------------------------------------------------------------
  ! Check the diagonalisation
  if(all(Occupations.eq.0.0_dp)) then
    call stp('No occupations in the canonical basis!')
  endif

  if(any(Occupations - 1.0_dp.gt.1d-5)) then
    call stp('Some occupations are bigger than one in the canonical basis.')
  endif
  where(Occupations.gt.1.0_dp) Occupations=1.0_dp
  if(any(Occupations .lt. -HFBNumCut)) then
!     Notice that we allow some very small negative occupation numbers.
!     They are entirely due to numerical error, and such errors are present
!     too in CR8. However, they are masked by setting all negative elements of
!     RhoHFB density matrix to zero. But in MOCCa, RHOHFB can be complex...
    where(Occupations.lt.0.0_dp) Occupations = 0.0_dp
    print *, 'Some occupations are smaller than zero in the canonical basis.'
  endif

  if(TRC) then
    !In the case of Time-reversal invariance, double the occupations.
    Occupations = 2.0_dp * Occupations
  endif
  !-----------------------------------------------------------------------------
  ! Note that we only need to actively store nwt of the canonical basis
  Columns = FindCorrectColumns()

  index = 1
  do it=1,Iindex
    do P=1,Pindex
      do i=1,nwt
        if(Columns(i,P,it).eq.0) exit !No more columns to take here
        ! The correct eigenvector of RhoHFB is this one:
        C = Columns(i,P,it)
        !-----------------------------------------------------------------------
        !Find the quantum numbers of this wavefunction
        loc = maxloc(abs(CanTransfo(:,C,P,it)))
        
        ii  = Blockindices(loc(1),P,it)
        iii = mod(ii-1,nwt)+1
        !-----------------------------------------------------------------------
        !note that we DO need to filter non-stored spwfs, but only in the case
        !where time-reversal is conserved, but signature isn't.
        S   = HFBasis(iii)%GetSignature()
        if(ii.ne.iii) S = -S
        TS  = HFBasis(iii)%GetTimeSimplex()
        P2  = HFBasis(iii)%GetParity()
        !Set the value to zero.
        call ResetWf(Canbasis(index))
        call CanBasis(index)%SetParity(P2)
        call CanBasis(index)%SetSignature(S)
        call CanBasis(index)%SetTimeSimplex(TS)
        call CanBasis(index)%SetIsospin(2*it-3)
        Energy = 0.0_dp
        !-----------------------------------------------------------------------
        ! Detect if the canonical basis level is just a level in the HF basis
!        mx  = abs(CanTransfo(loc(1),C,P,it))
!        if(abs(mx - 1).lt.HFBNumcut) then
!            CanDeriv(index)           = mod(blockindices(loc(1),P,it)-1,nwt)+1
!        else
!            CanDeriv(index)           = 0
!        endif
        !-----------------------------------------------------------------------
        ! Make the transformation
        do j=1,blocksizes(P,it)
          jj = Blockindices(j,P,it)
          jjj= mod(jj-1,nwt)+1
          !Cutoff for numerical stability
          if(abs(CanTransfo(j,C,P,it)) .lt.HFBNumCut) cycle
          !---------------------------------------------------------------------
          ! Get the phase of the not-transformed levels correctly
!          if(CanDeriv(index).ne.0) then
!            CanTransfo(j,C,P,it) = abs(CanTransfo(j,C,P,it))
!          endif
          !---------------------------------------------------------------------
          !Transformation
          !---------------------------------------------------------------------
          ! Note that the time-reversing here is kinda wasting cpu cycles
          ! since we apply it on the entire wavefunction. However, I prefer
          ! this routine to act on wavefunctions, not the (more efficient)
          ! spinors, since in this way MOCCa checks the consistency of
          ! quantum numbers.
          if(jj.eq.jjj) then
            if(TSC) then
                do k=1,nx*ny*nz*4
                  CanBasis(index)%Value%Grid(k,1,1,1,1) = Canbasis(index)%Value%Grid(k,1,1,1,1) + &
                  &           DBLE(CanTransfo(j,C,P,it))* HFBasis(jjj)%Value%Grid(k,1,1,1,1)
                enddo
            else
                CanBasis(index)%Value = Canbasis(index)%Value +       CanTransfo(j,C,P,it)*&
                &                                                    HFBasis(jjj)%Value
            endif
          else
            ! Note that this should only happen when signature is not conserved,
            ! but TimeReversal is.
            if(TSC) then
                CanBasis(index)%Value = Canbasis(index)%Value + DBLE(CanTransfo(j,C,P,it))*&
                &                                          TimeReverse(HFBasis(jjj)%Value)
            else
                CanBasis(index)%Value = Canbasis(index)%Value +       CanTransfo(j,C,P,it)*&
                &                                          TimeReverse(HFBasis(jjj)%Value)
            endif
          endif
        enddo

        !---------------------------------------------------------------------
        ! Compute contribution to the dispersion in the particle number
        ! DN^2 = 2*Tr(Rho(1-Rho)) = 2*sum v^2 (1 - v^2) in the canonical basis
        ! since we just diagonalised rho.
        if(TRC) then
          !Attention to the factors 2.
          ! Occupation is doubled, but we sum over only half of the states
          PairingDisp(it)=PairingDisp(it)+                                   &
          &     Occupations(C,P,it)*(2.0_dp-Occupations(C,P,it))
        else
          PairingDisp(it)=PairingDisp(it)+                                   &
          &   2*Occupations(C,P,it)*(1.0_dp-Occupations(C,P,it))
        endif
        call Canbasis(index)%SetOcc(Occupations(C,P,it))
        call CanBasis(index)%SymmetryOperators()
        index = index + 1
      enddo
    enddo
  enddo
  if(index.ne.nwt+1) then
    call stp('Not enough canonical spwfs were constructed.')
  endif
  !stop
  !-----------------------------------------------------------------------------
  ! Another important observable is the single-particle energy of these
  ! wavefunctions, but this module has no access to the hPsi routine
  ! defined in the Imaginarytime module, so this is calculated in the Main
  ! part of the program, until I think of a better way.
  do it=1,Iindex
    do P=1,Pindex
      do i=1,blocksizes(P,it)
        RhoII = DBLE(RhoHFB(i,i,P,it))
        if(TRC) RhoII = 2* RhoII
        ii = blockindices(i,P,it)
        iii=mod(ii-1,nwt)+1
        call HFBasis(iii)%SetOcc(RhoII)
      enddo
    enddo
  enddo
  !-----------------------------------------------------------------------------
  ! Measure convergence
!  do it=1,2
!     do P=1,2   
!        N = blocksizes(P,it)
!        drho  (P,it) = sum( (RhoHFB(1:N,1:N,P,it)   - OldRhoHFB(1:N,1:N,P,it))**2)
!        dkappa(P,it) = sum( (KappaHFB(1:N,1:N,P,it) - OldKappaHFB(1:N, 1:N, P,it))**2)
!     enddo
!  enddo

  !-----------------------------------------------------------------------------
  !Mix the densities, if there is a saved density
  if( .not. all(OldRhoHFB.eq.0.0_dp) ) then
    RhoHFB   = HFBMix * RhoHFB   + (1.0_dp - HFBMix) * OldRhoHFB
  endif
  if( .not. all(KappaHFB.eq.0.0_dp) ) then
    KappaHFB = HFBMix * KappaHFB + (1.0_dp - HFBMix) * OldKappaHFB
  endif
  !-----------------------------------------------------------------------------
  ! save number parities of fully manipulated qp vacuum for later use
  call CalcNumberParities()

  !-----------------------------------------------------------------------------
  ! Save old density and anomalous density matrix.
  do it=1,Iindex
    do P=1,Pindex
      N = blocksizes(P,it)
      OldRhoHFB(1:N,1:N,P,it)   = RhoHFB(1:N,1:N,P,it)
      OldKappaHFB(1:N,1:N,P,it) = KappaHFB(1:N,1:N,P,it)
    enddo
  enddo

  end subroutine HFBOccupations

  function FindCorrectColumns() result(Columns)
  !-----------------------------------------------------------------------------
  ! A subroutine that finds the correct eigenvectors of the HFB density matrix
  ! RhoHFB when Time Reversal is conserved. Since we would like to store only
  ! half of the canonical basis, we need to figure out what eigenvectors are
  ! needed.
  ! The result of the function are the indices of the columns containing proper
  ! eigenvectors.
  !-----------------------------------------------------------------------------
  ! When signature is conserved, this is easy, we need only the states with
  ! positive signature.
  ! When signature is not conserved, we just take half of them, but we need
  ! to take care that we do not select a wavefunction and its timereverse.
  !-----------------------------------------------------------------------------
  ! Note that it is important that Transformation and Occupations are the result
  ! of a LAPACK diagonalisation routine: when signature is conserved we abuse
  ! the eigenvalue ordering.
  !-----------------------------------------------------------------------------
    integer          :: Columns(nwt,Pindex,Iindex), i,index,P,it,N,P2,j

    Columns=0
    !---------------------------------------------------------------------------
    if(.not. TRC) then
        do it=1,Iindex
          do P=1,Pindex
            !Keep all columns
            do i=1,blocksizes(P,it)
                Columns(i,P,it) = i
            enddo
          enddo
        enddo
    !---------------------------------------------------------------------------
    else
        P2 = 0
        if(SC) then
          !If signature is conserved, we look for the columns with positive
          ! signature and save them.
          do it=1,Iindex
            do P=1,Pindex
              index = 1
              do i=1,blocksizes(P,it)
                N = Blocksizes(P,it)/2
                !Check for any components with negative signature.
                if(any(abs(CanTransfo(N+1:2*N,i,P,it)).gt.HFBNumCut)) then
                    ! When a component with signature - is detected, this
                    ! wavefunction does not need to be stored.
                else
                    Columns(index,P,it) = i
                    index               = index + 1
                endif
              enddo
              P2 = P2 + index - 1
            enddo
          enddo
        else
        !-----------------------------------------------------------------------
        ! Since Time-reversal is still conserved, the eigenvalues of the density
        ! matrix come in time-reversed pairs. LAPACK routines (like ZHEEV) order
        ! them by increasing eigenvalue, and it is thus sufficient for us to take
        ! the 'first' one of each pair.
        ! This might seem risky in the case where additional degeneracies show up
        ! but ZHEEV manages to separate eigenvalues to incredible high accuracy
        ! (~10^{-16}, absolute, not relative!) and this program will probably
        ! never achieve additional degeneracies on that level of precision.
        !-----------------------------------------------------------------------
        ! Alternatively, if this ever gives rise to problems, this piece of code
        ! might be replaced with another that explicitly checks for if new
        ! columns are time-reversed partners of already saved ones. I (WR) tried
        ! this first, but was unable to get it working correctly, and out of
        ! frustration, I implemented this.
        !-----------------------------------------------------------------------
          do it=1,Iindex
            do P=1,Pindex
              do i=1,blocksizes(P,it)/2
                Columns(i,P,it) = 2*(i-1)+1
              enddo
            enddo
          enddo
        endif
    !---------------------------------------------------------------------------
    endif
  end function FindCorrectColumns

  subroutine GuessHFBMatrices(PairingType,Fermi)
  !-----------------------------------------------------------------------------
  ! In order to start calculations with bad initial guesses and to write
  ! something to file, this routine constructs Fileblocksizes (if needed)
  ! and takes a guess at KappaHFB.
  ! It also allocates and puts to zero RhoHFB, U, V and Cantransfo.
  !-----------------------------------------------------------------------------
  ! Up to August 2019, the algorithm used to guess kappa fails for broken
  ! signature. 
  !-----------------------------------------------------------------------------
  ! Note that this DESTROYS ALL INFORMATION on HFB that MOCCa currently has.
  !-----------------------------------------------------------------------------
  ! For HFB, the routine generates a local copy of overlapT(:,:,:,:). I (MB)
  ! did not manage to get access to the array in Pairing without generating
  ! circular dependences. The array allocated and destroyed afterwards.
  !-----------------------------------------------------------------------------

    integer, intent(in)        :: PairingType
    integer                    :: i, ii, it, P, sig1,sig2, j,jj,iii,jjj, S1, S2
    integer                    :: startit, stopit
    integer                    :: IndCon(HFBSize) 
    logical                    :: FoundCon(HFBSize) 
    real(KIND=dp), intent(in)  :: Fermi(2)
    real(KIND=dp)              :: Occ , maxO
    real(KIND=dp), allocatable :: overlapT(:,:,:,:)

    if(.not. allocated(RhoHFB)) then
      allocate(RhoHFB(HFBSize,HFBSize,2,2))    ; RhoHFB = 0.0_dp
    endif
    if(.not. allocated(U)) then
      allocate(U(HFBSize,2*HFBSize,2,2))         ; U      = 0.0_dp
    endif
    if(.not. allocated(V)) then
      allocate(V(HFBSize,2*HFBSize,2,2))         ; V      = 0.0_dp
    endif

    if(.not.allocated(CanTransfo)) then
      allocate(CanTransfo(HFBSize,HFBSize,2,2)) ; CanTransfo = 0.0_dp
    endif

    select case (PairingType)
    case(0)
      !------------------------------------------------------------------------
      ! HF Calculation: we can't readily guess zero, since that would
      ! not allow us to start HFB calculations from converged HF calculations.
      ! So we try something more clever.
      !------------------------------------------------------------------------
      if(.not.allocated(kappahfb)) allocate(KappaHFB(HFBSize,HFBSize,2,2))
      KappaHFB=0.0_dp
      call PrepareHFBModule()
      do it=1,Iindex
        do P=1,Pindex
          do j=1,Blocksizes(P,it)
           ! We guess some 'small' occupation
            Occ = 0.1_dp
            jj  = blockindices(j,P,it)
            jjj = mod(jj-1,nwt)+1
            sig1 = HFBasis(jjj)%GetSignature()
            if(jjj.ne.jj) sig1 = -sig1
            do i=1,Blocksizes(P,it)
              ii  = blockindices(i,P,it)
              iii = mod(ii-1,nwt)+1
              sig2 = HFBasis(iii)%GetSignature()
              if(ii.ne.iii) sig2 = -sig2
              !Don't pair same parity states
              if( sig1 * sig2 .gt. 0 ) cycle
              KappaHFB(i,j,P,it) =  Occ ! = sqrt( u * v )
              KappaHFB(j,i,P,it) = -Occ
            enddo
          enddo
        enddo
      enddo
      return

    case(1)
      !BCS Calculation
      if(.not.allocated(kappahfb)) allocate(KappaHFB(HFBSize,HFBSize,2,2))
      KappaHFB=0.0_dp
      call PrepareHFBModule()
      !--------------------------------------------------------------------------
      ! We now make an initial guess as in EVCR8
      ! The only elements that are non-zero are those between wavefunctions
      ! and their time-reversed partners.
      do it=1,Iindex
        do P=1,Pindex
          do i=1,Blocksizes(P,it)/2
            ii  = Blockindices(i,P,it)
            iii = mod(ii-1,nwt)+1
            !Factor of 2 because of time-reversal
            Occ= HFBasis(iii)%GetOcc()/2.0_dp
            KappaHFB(i,i+Blocksizes(P,it)/2,P,it) = Occ * (1.0_dp - Occ)

            !For numerical stability
            if(abs(KappaHFB(i,i+Blocksizes(P,it),P,it)).ge.0.0_dp) then
              KappaHFB(i,i+Blocksizes(P,it)/2,P,it) =                            &
              &                        sqrt(KappaHFB(i,i+Blocksizes(P,it)/2,P,it))
            else
              KappaHFB(i,i+Blocksizes(P,it),P,it) = 0.0_dp
            endif
            !Kappa is antisymmetric
            KappaHFB(i+Blocksizes(P,it)/2,i,P,it) =                              &
            &                              - KappaHFB(i,i+Blocksizes(P,it)/2,P,it)
          enddo
        enddo
      enddo

    case(2)
        !---------------------------------------------------------------------------
        ! HFB calculation: Kappa should already be ok, but when asked MOCCa should
        ! be able to make some guess.
        !---------------------------------------------------------------------------
        ! MB, 19/08/24: Note that the signs of kappa assume a phase convention,
        ! which might not be respected by the single-particle states. As a 
        ! consequence, pairing might be accidentally repulsive for some states.
        !---------------------------------------------------------------------------
        ! MB: This fails by construction when signature is broken in the 
        ! calculation and by the actual state.
        !---------------------------------------------------------------------------
   !    ! old pre 19/08/30
   !
   !    Occ = 0.1_dp
   !    do it=1,Iindex
   !      KappaHFB(:,:,:,it)=0.0_dp
   !      do P=1,Pindex
   !        do i=1,blocksizes(P,it)
   !          ii = Blockindices(i,P,it)
   !          iii = mod(ii-1,nwt)+1
   !          S1 = HFBasis(iii)%GetSignature()
   !          if(iii.ne.ii) S1 = - S1
   !          do j=1,i-1        ! MB, 19/08/30: should that not be j=i+1,blocksizes(P,it)
   !            jj = Blockindices(j,P,it)
   !            jjj = mod(jj-1,nwt)+1
   !            S2 = HFBasis(jjj)%GetSignature()
   !            if(jjj.ne.jj) S2 = - S2
   !            if(S2.eq.S1) cycle
   !            KappaHFB(i,j,P,it) = Occ
   !            KappaHFB(j,i,P,it) =-Occ
   !          enddo
   !        enddo
   !      enddo
   !    enddo

        !-----------------------------------------------------------------------
        ! Initialise BCS-like kappa by scanning for largest overlap between HF
        ! states |< psi_i | T psi_j >|^2 and setting kappa for these "pairs" 
        ! to 0.2 * cutoff_i * cutoff_j
        !-----------------------------------------------------------------------
        ! multiplication with cutoff factors limits gaps to the pairing window.
        !-----------------------------------------------------------------------
        ! When this routine is called, neither the cutoff factors nor the 
        ! overlaps of HF states have been calculated already, so we have to
        ! do this first.
        !-----------------------------------------------------------------------
        call HFBoverlapT2(overlapT)
        call ComputePairingCutoffs(Fermi)

        Occ = 0.25_dp
        do it=1,Iindex
          KappaHFB(:,:,:,it)=0.0_dp
          do P=1,Pindex
            IndCon  (:) = 0
            FoundCon(:) = .false.
            do i=1,blocksizes(P,it)
              maxO = -9999.999_dp
              do j=1,blocksizes(P,it)
                if ( FoundCon(j) ) cycle   ! this index is already paired 
                if ( abs(overlapT(i,j,P,it)) .gt. maxO ) then 
                  IndCon(i) = j
                  maxO = abs(overlapT(i,j,P,it))
                endif
              enddo
              FoundCon(IndCon(i)) = .true.
              ii = Blockindices(i,P,it)
              iii = mod(ii-1,nwt)+1
              jj = Blockindices(IndCon(i),P,it)
              jjj = mod(jj-1,nwt)+1
              KappaHFB(i,IndCon(i),P,it) =  Occ !* PCutoffs(iii) * PCutoffs(jjj)
              KappaHFB(IndCon(i),i,P,it) = -Occ !* PCutoffs(iii) * PCutoffs(jjj)
              ! print '(" GuessHFBMatrices ",2i3,2i5,8f7.3)',it,P,i,IndCon(i), &
              !                 & overlapT(i,IndCon(i),P,it),PCutoffs(iii),PCutoffs(jjj), &
              !                 & KappaHFB(i,IndCon(i),P,it),KappaHFB(IndCon(i),i,P,it)
            enddo
          enddo
        enddo
      deallocate(overlapT)

    case DEFAULT
        call stp('Unknow PairingType in WriteOutKappa.')
    end select

  end subroutine GuessHFBMatrices

  function LNCr8_nosig(Delta, DeltaLN, flag) result(LNLambda)
  !-----------------------------------------------------------------------------
  ! Function that calculates Lambda2 as in CR8, which is to mean in a way that
  ! is completely incomprehensible.
  !
  ! Flag is an input parameter that checks if the pairing does not break down.
  ! If it is different from zero on exit, this means that LNLambda at isospin
  ! it could not be calculated.
  !
  ! NOTE:
  ! - Don't use this routine when S^T_y is broken
  !-----------------------------------------------------------------------------
    integer, intent(inout) :: flag(2)
    complex(KIND=dp), intent(in), allocatable :: Delta(:,:,:,:)
    complex(KIND=dp), intent(in), allocatable :: DeltaLN(:,:,:,:)
    real(KIND=dp) :: LNLambda(2)
    integer       :: i,j,it, ii, k, P, iii
    real(KIND=dp) :: c2(2),c3(2),c4(2), trx(2), txd(2), ex(2), gkr(2), erx(2)
    real(KIND=dp) :: gkx(2),gky(2)
    real(KIND=dp) :: E
    complex(KIND=dp) :: x1, x2, x3,x4,x5
    real(KIND=dp)    :: hl1(2) , hl2(2), deno(2), xnum(2)

    complex(KIND=dp),allocatable :: Chi(:,:,:,:)
    complex(KIND=dp),allocatable :: Chika(:,:,:,:)
    complex(KIND=dp),allocatable :: gamka(:,:,:,:)


    c2 = 0.0_dp ; c3 =0.0_dp ; c4 = 0.0_dp
    trx= 0.0_dp ; txd=0.0_dp ; ex = 0.0_dp ; gkr = 0.0_dp ; erx = 0.0_dp
    gkx = 0.0_dp ; gky = 0.0_dp

    allocate(Chi(HFBSize,HFBSize,2,2))
    allocate(Chika(HFBSize,HFBSize,2,2))
    allocate(gamka(HFBSize,HFBSize,2,2))
    do i=1,HFBSize
        do j=1,HFBSize
                Chi(j,i,:,:) = 0
                Chika(j,i,:,:) = 0
                Gamka (j,i,:,:) = 0
        enddo
    enddo
    do it=1,Iindex
      do P=1,Pindex
        do i=1,blocksizes(P,it)
          c2(it) = c2(it) + DBLE(RhoHFB(i,i,P,it))
          do j=1,blocksizes(P,it)
            c2(it) = c2(it) - DBLE(RhoHFB(i,j,P,it)*RhoHFB(j,i,P,it))
          enddo
        enddo
      enddo
    enddo
    c2 = 2*c2
    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          do i=1,blocksizes(P,it)
            Chi(i,j,P,it) = RhoHFB(i,j,P,it)
            do k=1,blocksizes(P,it)
              Chi(i,j,P,it) = Chi(i,j,P,it)  -RhoHFB(i,k,P,it)*RhoHFB  (k,j,P,it)
            enddo
          enddo
        enddo
        do j=1,blocksizes(P,it)
          do i=1,blocksizes(P,it)
            do k=1,blocksizes(P,it)
              Chika(i,j,P,it)=Chika(i,j,P,it)+Chi(i,k,P,it)   *KappaHFB(k,j,P,it)
              Gamka(i,j,P,it)=Gamka(i,j,P,it)+RhoHFB(i,k,P,it)*KappaHFB(k,j,P,it)
            enddo
            Chika(i,j,P,it) = KappaHFB(i,j,P,it) - 8 * Chika(i,j,P,it)
            Gamka(i,j,P,it) = KappaHFB(i,j,P,it) - 2 * Gamka(i,j,P,it)
          enddo
        enddo
      enddo
    enddo

    do it=1,Iindex
      do P=1,Pindex
        do i=1,Blocksizes(P,it)
          x1 = 0.0_dp ; x2 = 0.0_dp ; x3 = 0.0_dp ; x4 = 0.0_dp ; x5 = 0.0_dp
          do j=1,Blocksizes(P,it)
            x1 = x1 + Chi(i,j,P,it)    * DBLE(Chi(j,i,P,it))
            x2 = x2 + RhoHFB(i,j,P,it) * Chi(j,i,P,it)
          enddo
          ii  = Blockindices(i,P,it)
          iii = mod(ii-1,nwt)+1
          E = HFBasis(iii)%GetEnergy()
          ex(it) = ex(it)  + dble(E*Chi(i,i,P,it))
          txd(it)= txd(it) + dble(x1)
          trx(it)= trx(it) + dble(x2)
          erx(it)= erx(it) + dble(E*x2)

          x3 = 0.0_dp ;  x4 = 0.0_dp ; x5 = 0.0_dp
          do j=1,Blocksizes(P,it)
            x3 = x3 + DeltaLN(i,j,P,it) * KappaHFB(i,j,P,it)
            x4 = x4 + Delta  (i,j,P,it) * Chika(i,j,P,it)
            x5 = x5 + DeltaLN(i,j,P,it) * Gamka(i,j,P,it)
          enddo
          gkr(it) = gkr(it) + dble(x3)
          gkx(it) = gkx(it) + dble(x4)
          gky(it) = gky(it) + dble(x5)
        enddo
      enddo
    enddo

    !print *, ex , txd, trx, erx, gkr, gkx, gky

    ! Note that gkr,gkx & gky are double that what they are in CR8. This is because we
    ! sum over all contributions, while they only sum over one signature.
    c3 = 2*c2 - 8  * trx
    c4 = 4*c2 - 48 * txd

    hl1= 2*(ex + gkr/2)
    hl2= 4*(ex - 2*erx) + 2*(gkx/2 + gky/2)
    deno= c4*c2+2*c2**3 -c3*c3
    xnum=hl2*c2-hl1*c3

    LNLambda = 0.0_dp
    do it=1,Iindex
      if(abs(deno(it)).gt.1d-8) then
        LNLambda(it) = xnum(it)/deno(it)
      else
        ! Signalling problem
        flag(it) = 1
        print *, 'c2,c3,c4', c2(it),c3(it),c4(it)
      endif
    enddo
  end function LNCr8_nosig

  function LNCR8_sig(Delta, DeltaLN, flag) result(LNLambda)
  !-----------------------------------------------------------------------------
  ! API parameters
  complex(KIND=dp), intent(in), allocatable :: Delta(:,:,:,:)
  complex(KIND=dp), intent(in), allocatable :: DeltaLN(:,:,:,:)
  integer, intent(inout)                    :: flag(2)
  real(KIND=dp)                             :: LNLambda(2)

  integer       :: i,j,it, ii, k, P, iii,N
  real(KIND=dp) :: c2(2),c3(2),c4(2), trx(2), txd(2), ex(2), gkr(2), erx(2)
  real(KIND=dp) :: gkx(2),gky(2)
  real(KIND=dp) :: E
  real(KIND=dp)    :: hl1(2) , hl2(2), deno(2), xnum(2)
  complex(KIND=dp) :: x1, x2, x3,x4,x5
  complex(KIND=dp),allocatable :: Chi(:,:,:,:)
  complex(KIND=dp),allocatable :: Chika(:,:,:,:)
  complex(KIND=dp),allocatable :: gamka(:,:,:,:)

  c2 = 0.0_dp ; c3 =0.0_dp ; c4 = 0.0_dp
  trx= 0.0_dp ; txd=0.0_dp ; ex = 0.0_dp ; gkr = 0.0_dp ; erx = 0.0_dp
  gkx = 0.0_dp ; gky = 0.0_dp ; 
  !Chi = 0.0_dp ; chika = 0.0_dp ; Gamka = 0.0_dp
  
  allocate(Chi(HFBSize,HFBSize,2,2))
  allocate(Chika(HFBSize,HFBSize,2,2))
  allocate(gamka(HFBSize,HFBSize,2,2))
  do i=1,HFBSize
        do j=1,HFBSize
                Chi(j,i,:,:) = 0
                Chika(j,i,:,:) = 0
                Gamka (j,i,:,:) = 0
        enddo
  enddo
  !return

  !-----------------------------------------------------------------------------
  ! Chi = Rho - Rho^2
  ! c2 = 2*Tr [ \rho ( 1 - \rho) ]
  do it=1,Iindex
    do P=1,Pindex
      N = blocksizes(P,it)
      Chi(1:N/2,1:N/2,P,it)     = RhoHFB(1:N/2,1:N/2,P,it)
      Chi(N/2+1:N,N/2+1:N,P,it) = RhoHFB(N/2+1:N,N/2+1:N,P,it)

      Chi(1:N/2,1:N/2,P,it) = Chi(1:N/2,1:N/2,P,it)         -                  &
      &                matmul(RhoHFB(1:N/2,1:N/2,P,it),RhoHFB(1:N/2,1:N/2,P,it))
      Chi(N/2+1:N,N/2+1:N,P,it) = Chi(N/2+1:N,N/2+1:N,P,it) -                  &
      &        matmul(RhoHFB(N/2+1:N,N/2+1:N,P,it),RhoHFB(N/2+1:N,N/2+1:N,P,it))

      do i=1,N
        c2(it) = c2(it) + DBLE(2*Chi(i,i,P,it))
      enddo
    enddo
  enddo
  !-----------------------------------------------------------------------------
  ! Chika = ( 1 - 8 * chi) kappa
  do it=1,Iindex
    do P=1,Pindex
      N = blocksizes(P,it)
      ChiKa(1:N/2,N/2+1:N,P,it) =                                              &
      &        matmul(Chi(1:N/2,1:N/2,P,it), DBLE(KappaHFB(1:N/2,N/2+1:N,P,it)))
      ChiKa(N/2+1:N,1:N/2,P,it) =                                              &
      &    matmul(Chi(N/2+1:N,N/2+1:N,P,it), DBLE(KappaHFB(N/2+1:N,1:N/2,P,it)))

      ChiKa(1:N/2,N/2+1:N,P,it) =                                              &
      &             KappaHFB(1:N/2,N/2+1:N,P,it) - 8 * ChiKa(1:N/2,N/2+1:N,P,it)
      ChiKa(N/2+1:N,1:N/2,P,it) =                                              &
      &             KappaHFB(N/2+1:N,1:N/2,P,it) - 8 * ChiKa(N/2+1:N,1:N/2,P,it)
    enddo
  enddo
  !-----------------------------------------------------------------------------
  ! GamKa =  (1 - 2 Rho)Kappa
  do it=1,Iindex
    do P=1,Pindex
      N = blocksizes(P,it)
      GamKa(1:N/2,N/2+1:N,P,it) =                                              &
      &     matmul(RhoHFB(1:N/2,1:N/2,P,it), DBLE(KappaHFB(1:N/2,N/2+1:N,P,it)))
      GamKa(N/2+1:N,1:N/2,P,it) =                                              &
      & matmul(RhoHFB(N/2+1:N,N/2+1:N,P,it), DBLE(KappaHFB(N/2+1:N,1:N/2,P,it)))

      GamKa(1:N/2,N/2+1:N,P,it) =                                              &
      &             KappaHFB(1:N/2,N/2+1:N,P,it) - 2 * GamKa(1:N/2,N/2+1:N,P,it)
      GamKa(N/2+1:N,1:N/2,P,it) =                                              &
      &             KappaHFB(N/2+1:N,1:N/2,P,it) - 2 * GamKa(N/2+1:N,1:N/2,P,it)
    enddo
  enddo
  do it=1,Iindex
    do P=1,Pindex
      N = blocksizes(P,it)
      do i=1,N
        x1 = 0.0_dp ; x2 = 0.0_dp ; x3 = 0.0_dp ; x4 = 0.0_dp ; x5 = 0.0_dp
        do j=1,N
          x1 = x1 + DBLE(Chi(i,j,P,it)) * DBLE(Chi(j,i,P,it))
          x2 = x2 + DBLE(RhoHFB(i,j,P,it)) * DBLE(Chi(j,i,P,it))
        enddo
        ii  = Blockindices(i,P,it)
        iii = mod(ii-1,nwt)+1
        E = HFBasis(iii)%GetEnergy()
        ex(it) = ex(it)  + dble(E*Chi(i,i,P,it))
        txd(it)= txd(it) + dble(x1)
        trx(it)= trx(it) + dble(x2)
        erx(it)= erx(it) + dble(E*x2)

        x3 = 0.0_dp ;  x4 = 0.0_dp ; x5 = 0.0_dp
        do j=1,Blocksizes(P,it)
          x3 = x3 + DeltaLN(i,j,P,it) * KappaHFB(i,j,P,it)
          x4 = x4 + Delta  (i,j,P,it) * Chika(i,j,P,it)
          x5 = x5 + DeltaLN(i,j,P,it) * Gamka(i,j,P,it)
        enddo
        gkr(it) = gkr(it) + dble(x3)
        gkx(it) = gkx(it) + dble(x4)
        gky(it) = gky(it) + dble(x5)
      enddo
    enddo
  enddo
  ! Note that gkr,gkx & gky are double that what they are in CR8. This is because we
  ! sum over all contributions, while they only sum over one signature.
  c3 = 2*c2 - 8  * trx
  c4 = 4*c2 - 48 * txd

  hl1= 2*(ex + gkr/2)
  hl2= 4*(ex - 2*erx) + 2*(gkx/2 + gky/2)
  deno= c4*c2+2*c2**3 -c3*c3
  xnum=hl2*c2-hl1*c3

  LNLambda = 0.0_dp
  do it=1,Iindex
    if(abs(deno(it)).gt.1d-8) then
      LNLambda(it) = xnum(it)/deno(it)
    else
      ! Signalling problem
      flag(it) = 1
      print *, 'c2,c3,c4', c2(it),c3(it),c4(it)
    endif
  enddo

  end function LNCR8_sig

  subroutine ReadBlockingInfo(Block)
    !---------------------------------------------------------------------------
    ! Subroutine that reads from input the parameters for blocking calculations.
    ! The integer Block decides the total number of quasiparticles.
    !
    !---------------------------------------------------------------------------
    ! Technical note: we cannot initialize the QPBlockind array here, since
    ! MOCCa is not guaranteed to have read the wavefunction file when this
    ! info is read.
    ! Instead BlockQP checks if this has been performed and if not calls
    ! the subroutine QPindices.
    !---------------------------------------------------------------------------

    integer, intent(in) :: Block
    integer             :: io , i, j
    integer, allocatable :: Blocked(:), Filled(:) , Partners(:)
    integer, allocatable :: Temp(:)

    Namelist /Blocking/ Blocked, Filled, Partners, aliyangle

    io = 0
    allocate(QPExcitations(Block))   ; QPExcitations   = 0
    allocate(QPFilled(Block))        ; QPFilled        = 1
    allocate(QPPartners(Block))      ; QPPartners      = 0
    allocate(QPBlockedBefore(Block)) ; QPBlockedBefore = 0

    allocate(Filled(Block))          ; Filled          = 1
    allocate(Partners(Block))        ; Partners        = 0
    allocate(QPisinDCP(Block))       ; QPisinDCP       = 0
    allocate(aliyangle(Block))       ; aliyangle       = 0

    ! Read the namelist
    allocate(Blocked(block))
    read(unit=*,NML=Blocking, iostat=io)
    if(io.ne.0) call stp('Error on reading blocked particle indices', 'Iostat', io)
 
    !---------------------------------------------------------------------------
    ! Change the aliyangle to radians
    aliyangle = aliyangle/180*pi
    
    if(any(Blocked.le.-5) .or. any(Blocked.gt.nwt)) then
       call stp('Invalid quasiparticle index.')
    endif
    if(any(Filled .lt. 0) .or. any(Filled .gt. 1)) then
       call stp('Invalid filling of blocked quasiparticle.')
    endif
    if(any(Partners.le.-1) .or. any(Partners.gt.nwt)) then
       call stp('Invalid quasiparticle partner index.')
    endif
    QPExcitations = Blocked
    QPFilled      = Filled 
    QPPartners    = Partners

    !---------------------------------------------------------------------------
    ! check if some qp's are in decorrelated pair, which the later on leads to
    ! the use of a dedicated (experimental) blocking procedure.
    ! There are several ways of specifying a decorrelated pair state. Asking for
    ! 
    !     &Blocking
    !     Blocked  = 069 , 161
    !     Filled   = 1 , 1
    !     Partners = 161 , 069
    ! or 
    !     &Blocking
    !     Blocked  = 069 , 161
    !     Filled   = 0 , 0
    !     Partners = 161 , 069
    !
    ! or permutations of the two columns, the blocking will be done by the 
    ! standard routine also handling broken pair blocking.
    !
    ! Asking for
    !
    !     &Blocking
    !     Blocked  = 069 , 069
    !     Filled   = 0 , 1
    !     Partners = 161 , 161
    !
    ! or permutations, the code will set the flag QPisinDCP(:) for both 
    ! blocked qps and use a different part of the routine BlockQuasiParticles
    ! to identify the columns to be exchanged. 
    ! Note: the latter only works if the blocked qp can be found in the pairing
    ! window where u^2 and v^2 are both non-zero.
    !---------------------------------------------------------------------------
    ! The piece of code below also checks if data asks to block the same qp
    ! twice, in which case the code stops. 
    ! The same qp cannot be blocked twice in a HFB state as after the first
    ! blocking it cannot be found anymore.
    !---------------------------------------------------------------------------
    do i=1,block
    do j=i+1,block
      if ( QPExcitations(i) .eq. QPEXcitations(j) ) then
        if ( QPFilled(i) .ne. QPFilled(j) ) then
          print '(" Data ask for decorrelated pair ",8i5)', &
          &   i,QPExcitations(i),QPPartners(i),QPFilled(i), &
          &   j,QPExcitations(j),QPPartners(j),QPfilled(j)
          QPisinDCP(i) = j
          QPisinDCP(j) = i
        else
          print '(" Cannot block same quasiparticle twice! ",8i5)', &
          &   i,QPExcitations(i),QPPartners(i),QPFilled(i), &
          &   j,QPExcitations(j),QPPartners(j),QPfilled(j)
          stop
        endif
      endif
      if ( QPExcitations(i) .eq. QPPartners(j) .and. & 
             & QPFilled(i) .ne. QPFilled(j) ) then
        print '(" Cannot block same quasiparticle twice! ",8i5)', &
        &   i,QPExcitations(i),QPPartners(i),QPFilled(i), &
        &   j,QPExcitations(j),QPPartners(j),QPfilled(j)
        stop
      endif
    enddo
    enddo
  end subroutine ReadBlockingInfo

  subroutine QPindices(Fermi)
    !---------------------------------------------------------------------------
    ! Subroutine that looks for the correct blockindices of the HFBasis
    ! wavefunctions.
    !
    !---------------------------------------------------------------------------

    integer :: i, N, P,it,j, index, indexc , fill , S, tempind, refIt, refP
    real(KIND=dp) :: diffe
    real(KIND=dp), intent(in) :: Fermi(2)

    if(.not. allocated(Blockindices)) then
        call stp('QPindices cannot work before the HFBmodule is properly'     &
            &  //' initialized.')
    endif

    N = size(QPexcitations)

    allocate(QPParities(N))        ; QPParities        = 0
    allocate(QPIsospins(N))        ; QPIsospins        = 0
    allocate(QPSignatures(N))      ; QPSignatures      = 0
    allocate(QPBlockInd(N))        ; QPBlockInd        = 0
    allocate(QPBlockPartnerInd(N)) ; QPBlockPartnerInd = 0
    allocate(QPBlockLog(N))        ; QPBlockLog        = 'x'
    do i=1,N
        index = QPExcitations(i)
 
        !-----------------------------------------------------------------------
        ! Checking for special values of the blocking
        !
        ! index = -1 : P = - 1, neutron
        ! index = -2 : P = + 1, neutron
        ! index = -3 : P = - 1, proton
        ! index = -4 : P = + 1, proton
        !
        !-----------------------------------------------------------------------
        select case(index)
        case(-1)
          ! Negative parity neutron
          refP = 1
          refIt= 1
        case(-2)
          ! Positive parity neutron
          refP = 2
          refIt= 1
        case(-3)
          ! Negative parity proton
          refP = 1
          refIt= 2
        case(-4)
          ! Positive parity proton
          refP = 2
          refIt= 2
        case DEFAULT
          refP = 0
          refIt = 0
        end select
        
        if(index.lt.0) then
            diffe = 1000000
            do j=1,nwt
                  It  = (HFBasis(j)%GetIsospin()+3)/2
                  P   = (HFBasis(j)%GetParity()+3)/2
                  
                  if(It .ne. refIt ) cycle
                  if(P  .ne. refP  ) cycle
                  
                  if(abs(HFBasis(j)%energy - Fermi(1)).lt.diffe) then
                        diffe = abs(HFBasis(j)%energy - Fermi(1))
                        tempind = j
                  endif
            enddo            
            index            = tempind
            QPExcitations(i) = index 
            QPFilled(i)      = 1
        endif

        ! Getting the Quantum Numbers
        P  = (HFBasis(index)%GetParity()+3)/2
        It = (HFBasis(index)%GetIsospin()+3)/2
        S  = HFBasis(index)%GetSignature()
        !Saving input
        QPParities(i)    = P
        QPIsospins(i)    = It
        QpSignatures(i)  = S
        !Change the HFBNumberParity accordingly
        HFBNumberparity(P,it) = -1 * HFBNumberparity(P,it)
    enddo

    index  = 0
    indexc = 0
    do i=1,N
        P = QPParities(i)
        it= QPIsospins(i)
        do j=1,blocksizes(P,it)
          if ( Blockindices(j,P,it) .eq. QPExcitations(i) ) index  = j
          if ( Blockindices(j,P,it) .eq. QPPartners(i)    ) indexc = j
        enddo
        QPBlockind(i)        = index
        QPBlockPartnerInd(i) = indexc
    enddo
  end subroutine QPindices

subroutine PrintBlocking
    !---------------------------------------------------------------------------
    ! Subroutine to print the info on the any blocked states.
    !
    !---------------------------------------------------------------------------

    integer :: N, i

    1 format('Blocking parameters')
    2 format('   Index  Filled  Partner  Parity Isospin Signature')
    3 format(2x,'-------------------------------------------------')
    4 format(i7, 2x, i6,2x,i6,2x,i6,2x,i6,2x,i6)
    5 format(2x,'-------------------------------------------------')
    6 format(2x,'Number parities')
    7 format(2x,'Parity',8x,' -',11x,' +')
    8 format(2x,'Neutron',7x,i2,11x,i2)
    9 format(2x,'Proton ',7x,i2,11x,i2)
   10 format(2x, 'Blocking non-self-consistent after gradient step.')
   11 format(2x, 'Blocking self-consistent in gradient step.')
   
   12 format(2x, ' Alirotation at iteration 0',/, &
   &         2x, '   Index   y-angle')
   13 format(i7,5x,f12.5 )
    print 1

    if(.not. allocated(QPExcitations)) return

    N = size(QPExcitations)
    print 2
    print 3
    do i=1,N
      print 4, QPExcitations(i), QPFilled(i), QPPartners(i), 2*QPParities(i)-3,2*QPIsospins(i)-3, QPSignatures(i)
    enddo
    print 5
    print *
    print 6
    print 7
    print 5
    print 8, HFBNumberparity(:,1)
    print 9, HFBNumberParity(:,2)
    
    if(trim(FermiSolver).eq.'GRADIENT') then
        print 5
        if(blockconsistent) then
            print 11
        else
            print 10
        endif
    endif
    print 5
    if(any(aliyangle.ne.0.0_dp)) then
        print 12
        do i=1,N
            ! Note that we already converted the ali-angle to radians
            print 13, QPexcitations(i),aliyangle(i)*180/pi
        enddo
        print 5  
    endif
  end subroutine PrintBlocking

 subroutine BlockQuasiParticles(Fermi)
    !---------------------------------------------------------------------------
    ! Subroutine that arranges the blocking of quasiparticles. It should be used
    ! after diagonalisation of the HFB harmiltonian, but before the calculation
    ! of the HFB density & anomalous density matrix.
    !
    ! In practice, this routine just changes columns in the V & U matrices.
    ! (See B. Ballys thesis.)
    !
    ! If we want to block quasiparticle 'a' then we change the U & V matrices as
    ! follows:
    !
    ! U^a_{ij} = U_{ij}   V^a_{ij} = V_{ij} if j!=a
    ! U^a_{ia} = V^*_{ia} V^a_{ij} = U^*_{ia}
    !---------------------------------------------------------------------------
    ! The identification of the quasiparticle state is rather intuitive, we just
    ! (de)excite the column in the HFB matrix which has the largest overlap
    ! with a given HF index.
    !---------------------------------------------------------------------------
    ! loc(1) is an array for historical reasons, but this is not needed anymore 
    ! in this subroutine.
    !---------------------------------------------------------------------------
    integer          :: N, i, index , j, P, it, C, K, loc(1), CT, B, D, DT
    integer          :: icc, indexc , ii , iii , iiii , icycle 
    integer          :: jj ,  numnegE , fill
    integer          :: AlreadyBlocked(2) = 0
    integer          :: minind(2,2,4)     = 0
    complex(KIND=dp) :: TempU(HFBSize), TempV(HFBSize)
    real(KIND=dp)    :: TempU2 , TempV2 , Overlap, theta , maxEq , minEq
    real(KIND=dp)    :: TempU2c, TempV2c, normU , normV, OverlapU, OverlapV
    real(KIND=dp), intent(in) :: Fermi(2)
    real(KIND=dp)    :: minE(2,2,4)       = 0.0

    N = size(QPExcitations)

    if(.not.allocated(QPBlockind)) call QPindices(Fermi)
  
    ! call checkUandVColumns(HFBcolumns)
    ! save non-manipulated U and V matrices before blocking for later diagnostic use
    if (.not.allocated(NonBlockedHFBColumns)) allocate( NonBlockedHFBColumns(2*HFBSize,2,2) )
    NonBlockedHFBColumns = HFBColumns

    do i=1,N
        index  = QPblockind(i)
        P      = QPParities(i)
        it     = QPIsospins(i)
        fill   = QPfilled(i)
        ! determine index of partner state in HF basis.
        ! blockindices(index,P,it) is the index of the blocked state in the HF
        ! basis. In principle, the index should remain the same throughout the
        ! iterations, but for the moment it is assumed that it is not (which 
        ! also has been the case for the determination of the partner level 
        ! based on Delta elsewhere in the code in the past)
        ii  = Blockindices(index,P,it)
        iii = mod(ii-1,nwt)+1

        indexc = 0
        if ( QPPartners(i) .eq. 0 ) then  ! partner level not defined in data
          icc  = HFBasis(iii)%GetPairPartner()
          do j=1,blocksizes(P,it)
            ! print '(" BlockQuasiParticles : ",10i4)', & 
            !    &         i,index,ii,iii,ic,j,Blockindices(j,P,it)
            if ( Blockindices(j,P,it) .eq. icc ) indexc = j
          enddo
          ! print '(" index indexc = ",8i5)',i,iii,ic,index,indexc
        else
          indexc = QPBlockPartnerInd(i)  
          ! print '(" index indexc = ",8i5)',i,index,indexc
        endif
        if (indexc .eq.0) print '(" WARNING: BlockQuasiParticles : indexc = 0")' 

        !----------------------------------------------------------------------
        ! It might happen for broken-down pairing that a QP is already blocked,
        ! meaning a HF state is filled and its partner state empty. Such 
        ! situation is tracked with this variable.
        !----------------------------------------------------------------------
        AlreadyBlocked(:) = 0

        !----------------------------------------------------------------------
        ! Identify the quasiparticle excitation
        !----------------------------------------------------------------------
        ! Note: During the first call of pairing in a MOCCa run, the 
        ! E_qp are zero and this routine should not be called as the U and V
        ! matrices and/or their indeation might not be realistic yet. The
        ! following test will (hopefully) identify cases where this is done.
        maxEq = maxval(abs(QuasiEnergies(:,P,it)))
        ! print '(" maxEq = ",1e12.4)',maxEq
        if ( maxEq < 0.000001_dp ) then
          print '(/," Warning: all E_qp are zero for it = ",i2," P = ",i2, & 
          &       " in BlockQuasiParticles",/)',it,P
        endif

        !----------------------------------------------------------------------
        ! Search among the currently occupied columns for which qp
        ! has the biggest overlap with the asked-for HF state. 
        !----------------------------------------------------------------------
        ! When dealing with broken-pair excitations:
        !
        ! Defining :
        !  b   HF state to be filled in the end with <jz> approx K
        ! -b   HF state with largest < b | T | -b > when scanning the matrix
        !      for constant b
        !  B   HF state to be empty in the end because it has same parity and 
        !      similar
        !       <j^2> as b and <jz> approx -K
        ! -B   HF state with largest < B | T | -B > when scanning the matrix
        !      for constant b
        !
        ! In many cases, one has B = -b and -B = b, but for strong time-reversal
        ! breaking and/or very different mixing of HF states of opposite
        ! signature, this is not always the case.
        ! 
        ! b and B have to be specified by the user on input, where defining the 
        ! partner state B is optional (but highly recommended). When B is not 
        ! read from input, the code uses -b instead for the partner state.
        ! 
        ! If the norm of the U component is basically zero, U cannot be scanned
        ! for b and V has to be scanned instead. 
        !
        ! When -b is different from B, scanning V_ij for -b will produce the 
        ! wrong blocked (i.e. occupied HF) state. 
        !
        ! Because of these two problems, the code has to distinguish four
        ! different possibilities to identify the qp state to be blocked.
        !
        ! To find the dominant HF state when the qp states are highly mixed,
        ! the overlap is renormalized with the norm of the component of the 
        ! qp state that is scanned.
        ! 
        ! 1) blocking the qp with the aim of having a filled blocked HF state.
        !
        ! If |U_mu^2| = int d^3r U_mu^*(r) U_mu(r) > 0.1, the code scans for 
        ! the largest absolute value of 
        !
        !        int d^3r psi^*_b(r) U_mu(r) / sqrt(|U_mu^2|)
        !
        ! by searching for the largest square of this expression. In the
        ! BCS-limit of HFB, where U_mu(r) = U_mumu psi_mu(r), the necessity for 
        ! the renormalization is evident. In this limit it is also evident
        ! why the qp state to be blocked cannot identified anymore in the 
        ! u_mu -> 0 limit. In this case, one has to scan the V_mu(r) component.
        ! Defining |V_mu^2| = int d^3r V_mu^*(r) V_mu(r) = 1 - |U_mu^2| (for
        ! normalized qp states), for |V_mu^2| = 1 - |U_mu^2| > 0.9, the 
        ! code scans the V_mu(r) for the largest absolute value of
        !
        !       int d^3r psi^*_B(r)  V_mu(r) / sqrt(|V_mu^2|)
        !
        !  if a partner level is specified on input, or
        ! 
        !       int d^3r psi^*_-b(r) V_mu(r) / sqrt(|V_mu^2|)
        !
        ! if there is not, by searching for the largest square of this 
        ! expression.
        ! In the BCS limit of HFB, where V_mu(r) = V_mu-mu psi_-mu(r), it is
        ! evident that one has to scan for the conjugate HF state -j. To this
        ! end, the label of the state to be identified with the conjugate one 
        ! should be specified on input. When it is not, the code will try to
        ! identify it by scanning the (square of complex) overlap(s) of the 
        ! HF state used as tag with the tme-reversed of all other HF states
        ! for the most likely candidate for the "conjugate state". 
        ! 
        ! 2) blocking the qp with the aim of having an empty blocked HF state.
        !
        ! If |U_mu^2| = int d^3r U_mu^*(r) U_mu(r) > 0.1, the code scans for 
        ! the largest absolute value of 
        !
        !        int d^3r psi^*_B(r) U_mu(r) / sqrt(|U_mu^2|)
        !
        !  if a partner level is specified on input, or
        !
        !        int d^3r psi^*_-b(r) U_mu(r) / sqrt(|U_mu^2|)
        !
        ! if there is not, by searching for the largest square of this 
        ! expression. 
        ! For |V_mu^2| = 1 - |U_mu^2| > 0.9, the code scans the V_mu(r) for 
        ! the largest absolute value of
        !
        !        int d^3r psi^*_b(r) V_mu(r) / sqrt(|V_mu^2|)
        !
        ! by searching for the largest square of this expression.
        !----------------------------------------------------------------------
        ! When dealing with decorrelated-pair excitations:
        !
        ! NOTE: There is no meaningful distinction between filled/empty states
        ! after blocking in this case.
        ! 
        ! NOTE: in the notation introduced above, for such state one always
        ! knows b and B as both have to be specified on input.
        !
        ! Decorrelated pairs can often be handled with the procedure described
        ! above for broken-pair blocking. 
        ! 
        ! For decorrelated-pair states identified through a specific format
        ! of the data (such that QPisinDCP is set, see comments in subroutine
        ! ReadBlockingInfo), the code branches out to an experimental
        ! handling of such states.
        !
        ! If |U_mu^2| = int d^3r U_mu^*(r) U_mu(r) > 0.1, the code scans for 
        ! the largest absolute values of 
        !
        ! .... to be described ...
        !----------------------------------------------------------------------
        loc(:)    = 0
        Overlap   = 0.0_dp 
        OverlapU  = 0.0_dp 
        OverlapV  = 0.0_dp 
        numnegE   = 0

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! scan for lowest positive qp energy in each block. These states might 
        ! have to be manipulated later on in order to arrive at a 0qp reference
        ! state with positive number parity (in each block).
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        if( PfSolver ) then
          minE  (P,it,:) = 2000000
          minind(P,it,:) = 0
          do j=1,blocksizes(P,it)
            if ( QuasiEnergies(HFBColumns(j,P,it),P,it) .lt. minE(P,it,1) ) then
              minE  (P,it,1) = QuasiEnergies(HFBColumns(j,P,it),P,it)
              minind(P,it,1) = HFBColumns(j,P,it)
            endif
          enddo

          do j=1,blocksizes(P,it)
            if (QuasiEnergies(HFBColumns(j,P,it),P,it) .lt. minE(P,it,2) .and. &
              & QuasiEnergies(HFBColumns(j,P,it),P,it) .gt. minE(P,it,1) ) then
              minE  (P,it,2) = QuasiEnergies(HFBColumns(j,P,it),P,it)
              minind(P,it,2) = HFBColumns(j,P,it)
            endif
          enddo

          do j=1,blocksizes(P,it)
            if (QuasiEnergies(HFBColumns(j,P,it),P,it) .lt. minE(P,it,3) .and. &
              & QuasiEnergies(HFBColumns(j,P,it),P,it) .gt. minE(P,it,2) .and. &
              & QuasiEnergies(HFBColumns(j,P,it),P,it) .gt. minE(P,it,1) ) then
              minE  (P,it,3) = QuasiEnergies(HFBColumns(j,P,it),P,it)
              minind(P,it,3) = HFBColumns(j,P,it)
            endif
          enddo

          do j=1,blocksizes(P,it)
            if (QuasiEnergies(HFBColumns(j,P,it),P,it) .lt. minE(P,it,4) .and. &
              & QuasiEnergies(HFBColumns(j,P,it),P,it) .gt. minE(P,it,3) .and. &
              & QuasiEnergies(HFBColumns(j,P,it),P,it) .gt. minE(P,it,2) .and. &
              & QuasiEnergies(HFBColumns(j,P,it),P,it) .gt. minE(P,it,1) ) then
               minE  (P,it,4) = QuasiEnergies(HFBColumns(j,P,it),P,it)
               minind(P,it,4) = HFBColumns(j,P,it)
            endif
          enddo
          ! do j=1,4
          !   print '(" QP ",i2," for it= ",i2," P= ",i2,i5,f12.4)', &
          !   &      j,it,P,minind(P,it,j),minE(P,it,j)
          ! enddo
        endif

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! check if a blocked QP in this isospin-parity block is already filled 
        ! and its partner non-filled, i.e. if time-reversal breaking itself 
        ! already broke the pair (in a non-paired or almost non-paired state).
        ! This case might require a different manipulation of the QP vacuum 
        ! than blocking in order to obtain the correct number parity.
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Note: The limit of 0.99 mostly targets cases at the breakdown of
        ! pairing.
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! question: should one not check both cases of "fill", independent on
        ! what is asked for in the data?
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        do j=1,blocksizes(P,it)
          normU = 0.0_dp 
          do jj=1,blocksizes(P,it)
            normU = normU + abs(U(jj,HFBColumns(j,P,it),P,it))**2
          enddo
          normV = abs(1.0_dp - normU) ! abs() to save the day in case it's -eps
          normU = sqrt(normU)
          normV = sqrt(normV)
          if ( fill .eq. 1 ) then
            TempU2 = DBLE(U(indexc,HFBColumns(j,P,it),P,it)**2)
            TempV2 = DBLE(V(index ,HFBColumns(j,P,it),P,it)**2)
          else
            TempU2 = DBLE(U(index ,HFBColumns(j,P,it),P,it)**2)
            TempV2 = DBLE(V(indexc,HFBColumns(j,P,it),P,it)**2)
          endif

          if ( TempV2 .gt. 0.99 ) AlreadyBlocked(1) = HFBColumns(j,P,it)
          if ( TempU2 .gt. 0.99 ) AlreadyBlocked(2) = HFBColumns(j,P,it)

          if ( TempV2 .gt. 0.7  .and. TempU2 .gt. 0.7  ) then
                print '(" Pair already broken? ",10i4,1f7.2,6e12.4)', &
              & i,it,P,index,indexc,                                &
              & Blockindices(index,P,it),Blockindices(indexc,P,it), &
              & j,HFBColumns(j,P,it),loc(1),                        &
              & QuasiEnergies(HFBColumns(j,P,it),P,it),             &
              & normU*normU,normV*normV,TempU2,TempV2
          endif

         ! now check the opposite "fill" case - diagnostics only
         if ( fill .eq. 0 ) then
            TempU2 = DBLE(U(indexc,HFBColumns(j,P,it),P,it)**2)
            TempV2 = DBLE(V(index ,HFBColumns(j,P,it),P,it)**2)
          else
            TempU2 = DBLE(U(index ,HFBColumns(j,P,it),P,it)**2)
            TempV2 = DBLE(V(indexc,HFBColumns(j,P,it),P,it)**2)
          endif

      !!  if ( TempV2 .gt. 0.99 ) AlreadyBlocked(1) = HFBColumns(j,P,it)
      !!  if ( TempU2 .gt. 0.99 ) AlreadyBlocked(2) = HFBColumns(j,P,it)

          if ( TempV2 .gt. 0.7  .and. TempU2 .gt. 0.7  ) then
                print '(" Pair wrongly broken? ",10i4,1f7.2,6e12.4)', &
              & i,it,P,index,indexc,                                &
              & Blockindices(index,P,it),Blockindices(indexc,P,it), &
              & j,HFBColumns(j,P,it),loc(1),                        &
              & QuasiEnergies(HFBColumns(j,P,it),P,it),             &
              & normU*normU,normV*normV,TempU2,TempV2
          endif
        enddo

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Indeed, the pair with the QP to be blocked is already broken by 
        ! time-reversal breaking. Check if further action is needed.
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if ( AlreadyBlocked(1) .gt. 0 .and. AlreadyBlocked(2) .gt. 0 ) then
      !   print '(" QP with ",2i5," already blocked. filled: ",i5," empty: ",i5)', &
      !     &   index,indexc,AlreadyBlocked(1),AlreadyBlocked(2)

          if( PfSolver ) then
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 1) number parity of the vacuum "NPVac(P,it)" equals the targeted
            !    number parity for this isospin-parity block. Don't do anything 
            !    concerning this blocked QP. As the pair with the targeted
            ! blocked state is broken already, there is nothing left to do.
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! NOTE: might need further refinement (it might generate additional 
            ! 2qp excitations), but I don't have example for that yet. There
            ! also might be complications when two QP from the same isospin-
            ! parity-block are to be blocked (decorrelated pairs!).
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if ( HFBNumberParity(P,it) .eq. (-1)**NPVac(P,it) ) then
       !      print '(" targeted number parity ",i2," is the one of the vacuum = ",i2, &
       !         &     " continue cautiously")',HFBNumberParity(P,it),(-1)**NPVac(P,it)
              cycle
            endif

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 2) number parity of the vacuum "NPVac(P,it)" differs from the 
            ! targeted one and there is a swapped QP, but the swapped QP is 
            ! not the "accidentally" blocked one. A configuration with the 
            ! targeted number parity and the targeted blocked state is 
            ! obtained by unswapping the previously swapped QP.
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if ( HFBNumberParity(P,it) .ne. (-1)**NPVac(P,it) .and. &
                 & QPswapped(P,it,1) .ne. 0 ) then
              if ( AlreadyBlocked(1) .ne. QPswapped(P,it,1) .and. &
                 & AlreadyBlocked(2) .ne. QPswapped(P,it,1) ) then
             !!  print '(" Unswap QP ",i5," to restore number parity")', &
             !!  & QPswapped(P,it,1)
                 loc(1) = QPswapped(P,it,1)
                 Overlap = 0.999
                 QPBlockLog(i) = 'S'
              endif
            endif
            if ( HFBNumberParity(P,it) .ne. (-1)**NPVac(P,it) .and. &
                 & QPswapped(P,it,1) .eq. 0 ) then
            !! print '(" Nothing swapped for it = ",i2," P = ",i2, & 
            !! & " number parity = ",i5)',it,P,NPVac(P,it)
               do j=1,4
                 if ( minind(P,it,j) .ne. AlreadyBlocked(1) .and. &
                    & minind(P,it,j) .ne. AlreadyBlocked(2) ) then
                    loc(1) = minind(P,it,j)
                    QPBlockLog(i) = 'M'
            !!      print '(" additionally block ",i5," with E_qp= ",f7.4)', &
            !!      &        minind(P,it,j),minE(P,it,j)
                    exit 
                 endif
               enddo
            endif

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 3) number parity of the nonblocked vacuum "NPVac(P,it)" differs 
            ! from the targeted one and the swapped QP is the "accidentally" 
            ! blocked one.
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! not sure if this can actually happen as the pair is broken, at
            ! least when HFBGauge = 0.
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if ( HFBNumberParity(P,it) .ne. NPVac(P,it) ) then
              do j=1,2
                if ( AlreadyBlocked(j) .eq. QPswapped(P,it,1) ) then
                  print '(" Already blocked QP ",i5," with it = ",i1, &
                  &       " and P = ",i2, "has been swapped before")',&
                  &       QPswapped(P,it,2),it,P
                endif
              enddo
            endif
          endif
        endif

        !======================================================================
        ! Above, we constructed a qp vacuum that has positive number parity
        ! which might have necessitated to swap the lowest QP of positive 
        ! energy in a given spin-isospin block with its partner state of 
        ! negative energy.
        ! It might happen that the swapped QP actually was meant to be blocked.
        ! As the active part of HFBColumns now contains its partner state of 
        ! opposite QP energy, the QP to be blocked cannot be found anymore in 
        ! the range of columns that build the QP vacuum. 
        ! To identify such cases, one has also to scan the columns that have 
        ! been swapped away and which are now outside of the scanning loop that 
        ! follows below. If the swapped-away QP has larger overlap with a 
        ! given HF state than all the non-swapped QPs, it will be unswapped 
        ! by the blocking.
        ! (Otherwise one would block the quasiparticle in the range of active 
        ! columns that had the second-largest overlap before swapping, which 
        ! is usually not the targeted configuration.)
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        if ( QPswapped(P,it,1) .ne. 0 ) then
          j=0
          normU = 0.0_dp
          do jj=1,blocksizes(P,it)
            normU = normU + abs(U(jj,QPswapped(P,it,2),P,it))**2
          enddo
          normV = abs(1.0_dp - normU) ! abs() to save the day in case it's -eps
          normU = sqrt(normU)
          normV = sqrt(normV)

          if ( fill .eq. 1 ) then
            TempU2 = DBLE(U(index ,QPswapped(P,it,2),P,it)**2)
            TempV2 = DBLE(V(indexc,QPswapped(P,it,2),P,it)**2)
          else
            TempU2 = DBLE(U(indexc,QPswapped(P,it,2),P,it)**2)
            TempV2 = DBLE(V(index ,QPswapped(P,it,2),P,it)**2)
          endif
         !!   print '(" Check swapped ",7i4,2i4,1f7.2,6e12.4)',     &
         !!   & i,it,P,index,indexc,                                &
         !!   & Blockindices(index,P,it),Blockindices(indexc,P,it), &
         !!   & QPswapped(P,it,2),loc(1),                           &
         !!   & QuasiEnergies(QPswapped(P,it,2),P,it),              &
         !!   & normU*normU,normV*normV,TempU2,TempV2

          if ( normU*normU .gt. 0.1_dp ) then
            if( TempU2/normU .gt. Overlap) then
               loc(1) = QPswapped(P,it,1)
               Overlap = TempU2/normU
               QPBlockLog(i) = 'U'
         !!   print '(" Block U ",7i4,"   S",2i4,1f7.2,6e12.4)',    &
         !!   & i,it,P,index,indexc,                                &
         !!   & Blockindices(index,P,it),Blockindices(indexc,P,it), &
         !!   & QPswapped(P,it,2),loc(1),                           &
         !!   & QuasiEnergies(QPswapped(P,it,2),P,it),              &
         !!   & normU*normU,TempU2,TempV2,TempU2/normU,             &
         !!   & Overlap 
            endif
          else
            if( TempV2/normV .gt. Overlap) then
               loc(1) = QPswapped(P,it,1)
               Overlap = TempV2/normV
               QPBlockLog(i) = 'V'
         !!   print '(" Block V ",7i4,"   S",2i4,1f7.2,6e12.4)',    &
         !!   & i,it,P,index,indexc,                                &
         !!   & Blockindices(index,P,it),Blockindices(indexc,P,it), &
         !!   & QPswapped(P,it,2),loc(1),                           &
         !!   & QuasiEnergies(QPswapped(P,it,2),P,it),              &
         !!   & normU*normU,TempU2,TempV2,TempV2/normV,             &
         !!   & Overlap 
            endif
          endif 
        endif

        !======================================================================
        ! Ok, the complications have been dealt with, now we scan the 
        ! active columns for the largest overlap with a given HF state.
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        do j=1,blocksizes(P,it)

          ! other squares of the overlap cannot be larger; stop looking further
          if (Overlap .gt. sqrt(0.51_dp) ) exit

          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! Check if this column has been exchanged by blocking before.
          ! If so, don't scan it again and go directly to the next column.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! In case of highly mixed qps this avoids situations where a column 
          ! is swapped by blocking because of its largest overlap with one of 
          ! the HF states used to tag the QP states to be blocked, and then 
          ! swapped back because the swapped qp has the larget overlap with 
          ! another HF state used to tag a second blocked QP. Such situation
          ! is formally impossible and only arises in the code because the 
          ! solumns are swapped sequentially instead of simultaneously.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          icycle = 0
          do iiii=1,i-1
            if (   QPParities(iiii) .eq. P  .and. &
                 & QPIsospins(iiii) .eq. it .and. &
                 & QPBlockedBefore(iiii) .eq. HFBColumns(j,P,it) ) then
          !!  print '(" already blocked, cycle",2i2,5i5)', &
          !!  &    it,P,iiii,HFBColumns(j,P,it)
              icycle = 1
            endif
          enddo
          if ( icycle .eq. 1 ) cycle

          normU = 0.0_dp 
          do jj=1,blocksizes(P,it)
            normU = normU + abs(U(jj,HFBColumns(j,P,it),P,it))**2
          enddo
          normV = abs(1.0_dp - normU) ! abs() to save the day in case it's -eps
          normU = sqrt(normU)
          normV = sqrt(normV)

          !====================================================================
          ! distinguish two cases: standard and decorrelated-pair
          ! excitations. We will handle the standard case first, which treats
          ! the blocked qps individually.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          if ( QPisinDCP(i) .eq. 0 ) then
            if ( fill .eq. 1 ) then
              TempU2 = DBLE(U(index ,HFBColumns(j,P,it),P,it)**2)
              TempV2 = DBLE(V(indexc,HFBColumns(j,P,it),P,it)**2)
            else
              TempU2 = DBLE(U(indexc,HFBColumns(j,P,it),P,it)**2)
              TempV2 = DBLE(V(index ,HFBColumns(j,P,it),P,it)**2)
            endif

          !!  print '(" Block ",10i4,1f7.2,6e12.4)',              &
          !!  & i,it,P,index,indexc,                                &
          !!  & Blockindices(index,P,it),Blockindices(indexc,P,it), &
          !!  & j,HFBColumns(j,P,it),loc(1),                        &
          !!  & QuasiEnergies(HFBColumns(j,P,it),P,it),             &
          !!  & normU*normU,normV*normV,TempU2,TempV2

            if ( normU*normU .gt. 0.1_dp ) then
              if( TempU2/normU .gt. Overlap ) then
                loc(1) = HFBColumns(j,P,it)
                Overlap = TempU2/normU
                QPBlockLog(i) = 'U'
          !!    print '(" Block U ",10i4,1f7.2,6e12.4)',              &
          !!    & i,it,P,index,indexc,                                &
          !!    & Blockindices(index,P,it),Blockindices(indexc,P,it), &
          !!    & j,HFBColumns(j,P,it),loc(1),                        &
          !!    & QuasiEnergies(HFBColumns(j,P,it),P,it),             &
          !!    & normU*normU,TempU2,TempV2,TempU2/normU,             &
          !!    & Overlap 
              endif

          !!  if ( TempU2/normU .gt. 0.02_dp ) then
          !!    print '(" Block U ",10i4,1f7.2,6e12.4)',              &
          !!    & i,it,P,index,indexc,                                &
          !!    & Blockindices(index,P,it),Blockindices(indexc,P,it), &
          !!    & j,HFBColumns(j,P,it),loc(1),                        &
          !!    & QuasiEnergies(HFBColumns(j,P,it),P,it),             &
          !!    & normU*normU,TempU2,TempV2,TempU2/normU,             &
          !!    & Overlap 
          !!  endif
            else
              if( TempV2/normV .gt. Overlap) then
                loc(1) = HFBColumns(j,P,it)
                Overlap = TempV2/normV
                QPBlockLog(i) = 'V'
          !!    print '(" Block V ",10i4,1f7.2,6e12.4)',              &
          !!    & i,it,P,index,indexc,                                &
          !!    & Blockindices(index,P,it),Blockindices(indexc,P,it), &
          !!    & j,HFBColumns(j,P,it),loc(1),                        &
          !!    & QuasiEnergies(HFBColumns(j,P,it),P,it),             &
          !!    & normU*normU,TempU2,TempV2,TempV2/normV,             &
          !!    & Overlap 
              endif
          !!  if ( TempV2/normV  .gt. 0.02_dp ) then
          !!    print '(" Block V ",10i4,1f7.2,6e12.4)',              &
          !!    & i,it,P,index,indexc,                                &
          !!    & Blockindices(index,P,it),Blockindices(indexc,P,it), &
          !!    & j,HFBColumns(j,P,it),loc(1),                        &
          !!    & QuasiEnergies(HFBColumns(j,P,it),P,it),             &
          !!    & normU*normU,TempU2,TempV2,TempV2/normV,             &
          !!    & Overlap 
          !!  endif
            endif 
          endif

          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! 
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! in case of decorrelated-pair states, one usually wants simulatenous
          ! large overlap of one tagged HF state with U and of the other tagged
          ! HF state with V of the same quasiparticle. 
          ! While for conserved time-reversal symetry this is always the case 
          ! for the largest overlaps, in case of broken time-reversal it is not.
          ! Particularly problematic seem to be K=1/2 and 3/2 levels that
          ! can become highly mixed.
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! 
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          if ( QPisinDCP(i) .gt. 0 ) then
            ! print '(" QPisinDCP ")'
            if ( normU*normU .gt. 0.01_dp .and. normV*normV .gt. 0.01_dp ) then
              if ( fill .eq. 1 ) then
                TempU2 = DBLE(U(indexc,HFBColumns(j,P,it),P,it)**2)
                TempV2 = DBLE(V(index ,HFBColumns(j,P,it),P,it)**2)
              else
                TempU2 = DBLE(U(index ,HFBColumns(j,P,it),P,it)**2)
                TempV2 = DBLE(V(indexc,HFBColumns(j,P,it),P,it)**2)
              endif
       !!     print '(" QPisinDCP : ",i5,2i2,5i5,1f7.2,8e12.4)',   &
       !!     & i,it,P,index,indexc,                               &
       !!     & j,HFBColumns(j,P,it),QPisinDCP(i),                 &
       !!     & QuasiEnergies(HFBColumns(j,P,it),P,it),            &
       !!     & normU*normU,normV*normV,                           &
       !!     & DBLE(U(index ,HFBColumns(j,P,it),P,it)**2)/normU,  &
       !!     & DBLE(V(indexc,HFBColumns(j,P,it),P,it)**2)/normV,  &
       !!     & DBLE(U(indexc,HFBColumns(j,P,it),P,it)**2)/normU,  &
       !!     & DBLE(V(index ,HFBColumns(j,P,it),P,it)**2)/normV

              if( TempV2/normV .gt. OverlapV ) then
                OverlapV = TempV2/normV
                loc(1) = HFBColumns(j,P,it)
                QPBlockLog(i) = 'D'
              endif
            endif
          endif

        enddo
   !!   if ( numnegE .gt. 1 ) then
   !!     print '(/," BlockQuasiParticles : it = ",i2," P = ",i2)',it,P
   !!     print '(  " there are ",i5," qp states with negative energy")',numnegE
   !!   endif
        if ( loc(1).lt.1 .or. loc(1).gt.2*blocksizes(P,it) )  then
          print '(/," Warning: for it = ",i2," P = ",i2," Blockindices ",2i5,  &
               &    " not found",/)',                                          &
               &     it,P,Blockindices(index,P,it),Blockindices(indexc,P,it)
        endif
        ! NOTE: loc(1) might be 0 if the qp state to be blocked has already
        ! negative qp energy (as in this case U has to be scanned for indexc
        ! and V for index). This does not cause harm as the check for the
        ! swap of indices made below does not find a column to exchange for C=0.
        C = loc(1)
 
        !-----------------------------------------------------------------------
        ! If an alirotation is requested, find the time-reversed partner.
        if(aliyangle(i).ne. 0.0_dp) then
            B = blocksizes(P,it)
            if(C.gt.B/2) then
                CT = C - B/2
            else
                CT = C + B/2
            endif           
            
            Theta = aliyangle(i)/2 
                        
            TempU(1:B) = U(1:B,C,P,it)
            TempV(1:B) = V(1:B,C,P,it)

            U(1:B,C,P,it)  =  cos(theta) * U(1:B,C,P,it)+ sin(theta) *  U(1:B,CT,P,it)
            V(1:B,C,P,it)  =  cos(theta) * V(1:B,C,P,it)+ sin(theta) *  V(1:B,CT,P,it)
            
            U(1:B,CT,P,it) = -sin(theta) * tempU(1:B)   + cos(theta) *  U(1:B,CT,P,it)
            V(1:B,CT,P,it) = -sin(theta) * tempV(1:B)   + cos(theta) *  V(1:B,CT,P,it)  
            
            D  = 2*B - C  + 1 
            DT = 2*B - CT + 1
            
            TempU(1:B) = U(1:B,D,P,it)
            TempV(1:B) = V(1:B,D,P,it)
            
            U(1:B,D,P,it)  =  cos(theta) * U(1:B,D,P,it)+ sin(theta) *  U(1:B,DT,P,it)
            V(1:B,D,P,it)  =  cos(theta) * V(1:B,D,P,it)+ sin(theta) *  V(1:B,DT,P,it)
            
            U(1:B,DT,P,it) = -sin(theta) * tempU(1:B)   + cos(theta) *  U(1:B,DT,P,it)
            V(1:B,DT,P,it) = -sin(theta) * tempV(1:B)   + cos(theta) *  V(1:B,DT,P,it)  
                    
        endif
        !-----------------------------------------------------------------------
        ! When found, exchange it with its conjugate partner.
        !-----------------------------------------------------------------------
        ! print '(" C = ",i4)',C
        do j=1,blocksizes(P,it)
            if(HFBColumns(j,P,it) .eq. C) then
              HFBColumns(j,P,it) = 2*blocksizes(P,it) - HFBColumns(j,P,it) + 1
        !!    print '(" QP exchange : ",2i2,2i5,f9.4)',it,P,j,  &
        !!    & HFBColumns(j,P,it),QuasiEnergies(HFBColumns(j,P,it),P,it)
              QPBlockedBefore(i) = HFBColumns(j,P,it)
            endif
        !!  print '(" QP spectrum after exchange : ",2i2,2i5,f9.4)', &
        !!  & it,P,j,HFBColumns(j,P,it),QuasiEnergies(HFBColumns(j,P,it),P,it)
        enddo
    enddo
  !  call checkUandVColumns(HFBcolumns)
  !  stop
  end subroutine BlockQuasiParticles

  subroutine DecomposeBlockedQuasiParticles()
  !-----------------------------------------------------------------------------
  ! print decomposition of blocked QP with its dominant HF components
  ! HFBColumns(i,P,it) translates block indexing to global indexing
  !-----------------------------------------------------------------------------
    integer        :: i, P, it, j , jHF , jHFc
    integer        :: ii , iii , jj , ifound , N , fill , icc
    integer        :: loc , index , indexc , domU(1), domV(1)
    integer        :: domUHF,domVHF,domUHFc,domVHFc
    real(KIND=dp)  :: normU , normV
    real(KIND=dp)  :: overlapU , overlapV , overlapUn , overlapVn
    real(KIND=dp)  :: TempU2 , TempV2 , overlap

  2 format(/,"   Index  Filled  Partner  Parity Isospin Signature LastScan")
  3 format(2x,"-----------------------------------------------------------")
  4 format(i7, 2x, i6,2x,i6,2x,i6,2x,i6,2x,i6,4x,a6)
  5 format(2x,"-----------------------------------------------------------")
  9 format(/,"  HFU : dominant HF state in U component", &
       &   /,"  HFV : dominant HF state in V component", &
       &   /,"  HFUc: HF state with largest overlap with time-reversed of HFU",  &
       &   /,"  HFVc: HF state with largest overlap with time-reversed of HFV",  &
       &   /,"  swapped            : swapped QP to ensure correct number parity of reference vacuum",   &
       &   /,"  swapped and blocked: blocking of a previously swapped QP, which changes number parity", &
       &   /,"              blocked: blocking, which changes number parity")
 10 format(/,4x,"it  P  iqp    E_qp   |U|^2    |V|^2   HFU  HFV HFUc HFVc", &
       &     "  |<HFU|>|^2 |<HF|V|>|^2",                                    &
       &   /,2x,82("-"))
 11 format(  " ")
 20 format(/,"  Decomposition of manipulated QP states ",                 &
       &   //,4x," it  P  iqp    E_qp  jHF jHFc   |<HF|U>|^2 |<HF|V>|^2"  &
       &    /,2x,56("-"))
 30 format(/,"  How the blocked qp states were selected from non-manipulated U and V",  &
       &  //,4x," comp it  P  HFb HFbc fill  iqp    E_qp    |<V|V>|^2 ! "               &
       &        "|<HFb|comp>|^2/|<comp|comp>|"                                          &
       &    /,3x,84("-"))

    if(.not. allocated(QPExcitations)) return

    N = size(QPExcitations)

    !--------------------------------------------------------------------------------
    ! repeat read data on blocking
    if ( N .gt. 0 ) then
      print 2
      print 3
      do i=1,N
        print 4,QPExcitations(i),QPFilled(i),QPPartners(i),           &
           &    2*QPParities(i)-3,2*QPIsospins(i)-3, QPSignatures(i), &
           &    QPBlockLog(i)
      enddo
      print 5
    endif
    !--------------------------------------------------------------------------------
    ! print diagnostic table with index information of dominant HF component in
    ! the U and V matrices before blocking
    ! This basically repeats the scan for qp states in BlockQuasiParticles()

    ifound = 0
    do i=1,N
        index  = QPblockind(i)
        P      = QPParities(i)
        it     = QPIsospins(i)
        fill   = QPfilled  (i)
        ii     = Blockindices(index,P,it)
        iii    = mod(ii-1,nwt)+1
        indexc = 0
        if ( QPPartners(i) .eq. 0 ) then  ! partner level not defined in data
                                          ! guess from < iiiHF | T | jjjHF >
          icc  = HFBasis(iii)%GetPairPartner()
          do j=1,blocksizes(P,it)
            if ( Blockindices(j,P,it) .eq. icc ) indexc = j
          enddo
          ! print '(" index indexc = ",8i5)',i,iii,ic,index,indexc
        else                               ! partner level read from data
          indexc = QPBlockPartnerInd(i)
          ! print '(" index indexc = ",8i5)',i,index,indexc
        endif
        if (indexc .eq.0) print '(" WARNING: indexc = 0")'

        ! Identify the quasiparticle excitation
        do j=1,blocksizes(P,it)
          normU = 0.0_dp
          do jj=1,blocksizes(P,it)
            normU = normU + abs(U(jj,NonBlockedHFBColumns(j,P,it),P,it))**2
          enddo
          normV = abs(1.0_dp - normU) ! abs() to save the day in case it's -eps
          normU = sqrt(normU)
          normV = sqrt(normV)
          if ( fill .eq. 1 ) then
            TempU2 = DBLE(U(index ,NonBlockedHFBColumns(j,P,it),P,it)**2)
            TempV2 = DBLE(V(indexc,NonBlockedHFBColumns(j,P,it),P,it)**2)
          else
            TempU2 = DBLE(U(indexc,NonBlockedHFBColumns(j,P,it),P,it)**2)
            TempV2 = DBLE(V(index ,NonBlockedHFBColumns(j,P,it),P,it)**2)
          endif
          if ( normU*normU .gt. 0.1_dp ) then
            if ( TempU2/normU .gt. 0.02_dp ) then
              if (ifound .eq. 0 ) print 30
              print '(2x,i3,"  U ",2i3,4i5,1f8.3,6f14.8)',          &
              & i,it,P,Blockindices(index,P,it),                    &
              & Blockindices(indexc,P,it),fill,                     &
              & NonBlockedHFBColumns(j,P,it),                       &
              & QuasiEnergies(NonBlockedHFBColumns(j,P,it),P,it),   &
              & normV*normV,TempU2,TempU2/normU
              ifound = ifound + 1
            endif
          else
            if ( TempV2/normV .gt. 0.02_dp ) then
              if (ifound .eq. 0 ) print 30
              print '(2x,i3,"  V ",2i3,4i5,1f8.3,6f14.8)',          &
              & i,it,P,Blockindices(index,P,it),                    &
              & Blockindices(indexc,P,it),fill,                     &
              & NonBlockedHFBColumns(j,P,it),                       &
              & QuasiEnergies(NonBlockedHFBColumns(j,P,it),P,it),   &
              & normV*normV,TempV2/normV
              ifound = ifound + 1
            endif
          endif
        enddo
    enddo
    if (ifound .gt. 0 ) print 11

    !--------------------------------------------------------------------------------
    ! print diagnostic table with index information of dominant HF component of 
    ! qp states with negative E_qp
    ifound = 0
    print 9
    if ( any(abs(HFBGauge) .gt. 0.001 )) then
      print '("  Attention: for HFBGauge(it) != 0, swapped and blocked QPs cannot be properly distinguished")'
    endif
    do it=1,Iindex
    do P =1,Pindex
      do i=1,blocksizes(P,it)
        ii = HFBColumns(i,P,it)
        if ( QuasiEnergies(ii,P,it)  < 0.0_dp .or. QPswapped(P,it,2) .eq. ii ) then
          if (ifound .eq. 0 ) print 10
          ifound = ifound + 1
          normU = 0.0_dp 
          do jj=1,blocksizes(P,it)
            normU = normU + abs(U(jj,HFBColumns(j,P,it),P,it))**2
          enddo
          normV = abs(1.0_dp - normU) ! abs() to save the day in case it's -eps
          normU = sqrt(normU) 
          normV = sqrt(normV) 
          ! Getting the dominant components of both U and V matrices in
          ! isospin-parity-block indexation
          domU    = maxloc(abs(U(:,ii,P,it)))
          domV    = maxloc(abs(V(:,ii,P,it)))
          ! corresponding indices in HF basis indexation
          domUHF  = blockindices(domU(1),P,it)
          domVHF  = blockindices(domV(1),P,it)
          ! indices of states with largest overlap <domUHF|T|domUHFc> and 
          ! <domVHF|T|domVHFc> in HF basis
          domUHFc = HFBasis(domUHF)%GetPairPartner()
          domVHFc = HFBasis(domVHF)%GetPairPartner()
          overlapU = abs(U(domU(1),ii,P,it))**2
          overlapV = abs(U(domV(1),ii,P,it))**2
          if ( QPswapped(P,it,1) .eq. ii )  then
            print '(3x,2i3,i5,1f8.3,2f9.6,4i5,2f12.8," swapped")',              &
              &      it,P,ii,QuasiEnergies(ii,P,it),                            &
              &      normU*normU,normV*normV,                                   &
              &      domUHF,domVHF,domUHFc,domVHFc,                             &
              &      overlapU,overlapV
            cycle
          endif
          if ( QPswapped(P,it,2) .eq. ii )  then
            print '(3x,2i3,i5,1f8.3,2f9.6,4i5,2f12.8," swapped and blocked")',  &
              &      it,P,ii,QuasiEnergies(ii,P,it),                            &
              &      normU*normU,normV*normV,                                   &
              &      domUHF,domVHF,domUHFc,domVHFc,                             &
              &      overlapU,overlapV
            cycle
          endif
          print '(3x,2i3,i5,1f8.3,2f9.6,4i5,2f12.8,"             blocked")',    &
            &      it,P,ii,QuasiEnergies(ii,P,it),                              &
            &      normU*normU,normV*normV,                                     &
            &      domUHF,domVHF,domUHFc,domVHFc,                               &
            &      overlapU,overlapV
        endif
      enddo
    enddo
    enddo
    print '(/,i3," manipulated quasiparticles found",/)',ifound

    !--------------------------------------------------------------------------------
    ! print diagnostic table with index information of dominant HF component in
    ! the U and V matrices after blocking
    ifound = 0
    do it=1,Iindex
    do P =1,Pindex
      do i=1,blocksizes(P,it)
        ii = HFBColumns(i,P,it)
        if ( QuasiEnergies(ii,P,it)  < 0.0_dp .or. QPswapped(P,it,2) .eq. ii ) then
          do j=1,blocksizes(P,it)
            if ( abs(U(j,ii,P,it))**2 > 0.01_dp ) then
              if (ifound .eq. 0 ) print 20
              ifound = ifound + 1
              ! corresponding indices in HF basis indexation
              jHF = blockindices(j,P,it)
              ! indices of states with largest overlap <jHF|T|jHFc> in HF basis
              jHFc = HFBasis(jHF)%GetPairPartner()
              print '("  U ",2i3,i5,1f8.3,2i5,2x,2f11.7)',           &
                &      it,P,ii,QuasiEnergies(ii,P,it),jHF,jHFc,     &
                &      abs(U(j,ii,P,it))**2,abs(V(j,ii,P,it))**2
            endif
            if ( abs(V(j,ii,P,it))**2 > 0.01_dp ) then
              if (ifound .eq. 0 ) print 20
              ifound = ifound + 1
              ! corresponding indices in HF basis indexation
              jHF = blockindices(j,P,it)
              ! indices of states with largest overlap <jHF|T|jHFc> in HF basis
              jHFc = HFBasis(jHF)%GetPairPartner()
              print '("  V ",2i3,i5,1f8.3,2i5,2x,2f11.7)',           &
                &      it,P,ii,QuasiEnergies(ii,P,it),jHF,jHFc,     &
                &      abs(U(j,ii,P,it))**2,abs(V(j,ii,P,it))**2
            endif
          enddo
        endif
      enddo
    enddo
    enddo
    if (ifound .gt. 0 ) print 11

  end subroutine DecomposeBlockedQuasiParticles

  subroutine PrintQP()
  !-----------------------------------------------------------------------------
  ! Print the quasiparticles obtained by the HFB proces.
  !-----------------------------------------------------------------------------
    integer             :: i, P, it, j, domU(1), domV(1)
    character(len=7)    :: Species(2)=(/ 'Neutron', 'Proton '/)
    logical             :: skip
    real(KIND=dp)       :: u2mv2, u2pv2, align(6,2*HFBSize,Pindex,IIndex)
    real(KIND=dp)       :: angquantum


    10  format (90 ('_'))
    20  format (90 ('-'))
     1  format ( a7, ' P=', i2 ' quasiparticles')
     2  format ('  n   <Rz>    E_qp      U    V   u^2-v^2 ',3x,      &
        &       '<Jx|T>',4x,'<Jy|T>',6x,'<Jz>', 5x, ' J ',  5x, 'Sx')
     3  format ('  n   <Rz>    E_qp      U    V   u^2-v^2 ',4x,      &
        &       '<Jx>',5x,'<Jy|T>',5x,'<Jz>', 7x, ' J ',  8x, 'Sx')
     4  format ('  n   <Rz>    E_qp      U    V   u^2-v^2 ',2x,      &
        &       '<Jx>',5x,'<Jy>',5x,'<Jz>', 5x, ' J ', 5x, 'Sx')
     !           n   <Rz>   E_qp
    99  format ( i3, f7.2 , f10.5,2x, i3,2x, i3, 6(3x, f7.2))

    if(.not.allocated(quasisimplex)) then
         allocate(QuasiSimplex(2*HFBSize,Pindex,Iindex)) ; Quasisimplex = 0.0_dp
    endif
    
    do it=1,Iindex
        do P=1,Pindex
            QuasiSimplex(:,P,it) = 0.0
            do j=1,2*blocksizes(P,it)
                do i=1,blocksizes(P,it)
                    Quasisimplex(j,P,it) = Quasisimplex(j,P,it)                &
                    & +HFBasis(blockindices(i,P,it))%xsimplexr*                &
                    & (dble(U(i,j,P,it))**2 - dble(V(i,j,P,it))**2) 
                enddo
            enddo
        enddo
    enddo

    ! un-commented the commented-out call of QPalignment(). I'm highly suspicious
    ! that the calculation of observables in that routine is either formally wrong 
    ! or, if formally correct, what is calculated is useless in practice. The most 
    ! interesting information would be how much a given qp contributes to many-body
    ! expectation values of one-body operators, which is entirely determined by the 
    ! V matrix only.
    align = QPalignment()
    do it=1,Iindex
        do P=1,Pindex
          print 10
          print 1, Species(it), 2*P-3

          if(SC.and.TSC) then
            print 2
          elseif(TSC) then
            print 3
          else
            print 4
          endif

          print 10
          do i=1,2*blocksizes(P,it)
            Skip =.true.
            do j=1,blocksizes(P,it)
                ! Only print the qp if it is 'selected' in the HFBColumns
                if( i.eq. HFBColumns(j,P,it) ) then
                    Skip = .false.
                    exit
                endif
            enddo
            if(abs(QuasiEnergies(i,P,it)).gt.QPPrintWindow ) Skip=.true.
            if(skip) cycle

            ! Getting the dominant components of both U and V matrices.
            domU = maxloc(abs(U(:,i,P,it)))
            domV = maxloc(abs(V(:,i,P,it)))

            ! Computing u^2 - v^2 and u^2 + v^2
            u2mv2 = 0.0_dp ; u2pv2 = 0.0_dp
            do j=1,blocksizes(P,it)
              u2mv2 = u2mv2 + abs(U(j,i,P,it))**2 - abs(V(j,i,P,it))**2
              u2pv2 = u2pv2 + abs(U(j,i,P,it))**2 + abs(V(j,i,P,it))**2
            enddo

            !----------------------------------------------------------
            ! Now find j so that <J^2> = j*(j+1)
            ! The solution is obviously:
            ! j = [- 1 + sqrt( 1 + 4 * <J^2>)]/2
            AngQuantum =                                                  &
            &-0.5_dp*(1.0_dp-sqrt( 1.0_dp + 4.0_dp*sum(align(4:6,i,P,it))))

            print 99, i, QuasiSignatures(i,P,it), QuasiEnergies(i,P,it)    &
            &      ,blockindices(domU(1),P,it), blockindices(domV(1),P,it) &
            &      ,u2mv2, 0.0d0, 0.0d0, align(3,i,P,it), angquantum,      &
            &      quasisimplex(i,P,it)
          enddo
        enddo
    enddo
    print 20
    call DecomposeBlockedQuasiParticles()
  end subroutine PrintQP

  function QPAlignment() result(align)
    !------------------------------------------------------------------------
    ! MB 19/01/18 This routine does not make sense to me. The j are one-body
    ! operators proportional to a^+ a, such that they can have only non-zero
    ! matrix elements with the part of the QP state proportional to a^+ | - >
    ! which is the V component. 
    !------------------------------------------------------------------------
    ! Function that computes the angular momentum of all quasiparticles.
    ! This should only be called at print-out iterations, since this
    ! is quite costly.
    !------------------------------------------------------------------------

    real(KIND=dp)    :: Jmatrix(HFBSize,HFBsize,6,2)
    real(KIND=dp)    :: align(6,2*HFBSize,Pindex,IIndex)
    integer          :: i,j,it,P,k,jj,kk, N
    complex(KIND=dp) :: Transfo(2*HFBSize, 2*HFBSize,Pindex,IIndex)
    logical          :: TRX=.false., TRY=.false., TRZ=.false.

    if(SC)          TRX=.true.
    if(TSC .or. SC) TRY=.true.
    TRZ = .false.

    !------------------------------------------------------------------------
    ! First, we get all relevant matrix elements of the form
    ! < i | J | j >
    ! in the HFBasis
    ! The matrix elements of the nwt wavefunctions that are always stored
    JMatrix=0.0_dp
    do j=1,nwt
        do i=1,j
            ! Care for the indices here; Angularmomentum takes !reversed!
            ! indices
            JMatrix(i,j,:,:) =   AngularMomentum(HFBasis(j),HFBasis(i),.true.,&
            &                                                       TRX,TRY,TRZ)
            JMatrix(j,i,:,1) =   JMatrix(i,j,:,1)
            JMatrix(j,i,:,2) = - JMatrix(i,j,:,2)
        enddo
    enddo

    if (TRC .and. .not.SC) then
       print '(/, "Calculation of J matrix elements not yet", &
        &      /," correctly implemented for signature",      &
        &      /," breaking and Time Reversal conserving",    &
        &      /," calculations.")'
      return
    endif

    ! Get the angular momentum of the time-reverse wavefunctions
    if(TRC) then
        do j=1,nwt
            do i=1,j
                JMatrix(i+nwt,j+nwt,1:3,:) = - JMatrix(i,j,1:3,:)
                JMatrix(i+nwt,j+nwt,4:6,:) =   JMatrix(i,j,4:6,:)
                JMatrix(j+nwt,i+nwt,:,1)   =   JMatrix(i+nwt,j+nwt,:,1)
                JMatrix(j+nwt,i+nwt,:,2)   =   JMatrix(i+nwt,j+nwt,:,2)
            enddo
        enddo
        if(.not.SC) call stp('Calculation of J matrix elements not yet' &
        &                  //' correctly implemented for signature'     &
        &                  //' breaking and Time Reversal conserving'   &
        &                  //' calculations.')
    endif

    ! Construct the transformation matrix, in order to not get confused with
    ! indices.
    do it=1,Iindex
        do P=1,Pindex
            N = blocksizes(P,it)
            Transfo(1:N    ,1:2*N,P,it) = U(1:N,1:2*N,P,it)
            Transfo(N+1:2*N,1:2*N,P,it) = V(1:N,1:2*N,P,it)
        enddo
    enddo

    align= 0.0_dp
    do it=1,iindex
        do P=1,Pindex
            N = blocksizes(P,it)
            do i=1,2*N
                do j=1,N
                    jj = blockindices(j,P,it)
                    do k=1,N
                        kk = blockindices(k,P,it)
                        Align(:,i,P,it) = Align(:,i,P,it) +   &
                        &                 Jmatrix(kk,jj,:,1)&
                        &                 *DBLE(Transfo(j,i,P,it)) &
                        &                 *DBLE(Transfo(k,i,P,it))
                    enddo
                enddo

               do j=N+1,2*N
                    jj = blockindices(j-N,P,it)
                    do k=N+1,2*N
                        kk = blockindices(k-N,P,it)

                        ! A minus sign for <Jx>, <Jy> and <Jz>
                        Align(1:3,i,P,it) = Align(1:3,i,P,it) -    &
                        &                 Jmatrix(kk,jj,1:3,1)     &
                        &                 *DBLE(Transfo(j,i,P,it)) &
                        &                 *DBLE(Transfo(k,i,P,it))
                        ! A plus sign for <Jx^2>, <Jy^2> and <Jz^2>
                        Align(4:6,i,P,it) = Align(4:6,i,P,it) +    &
                        &                 Jmatrix(kk,jj,4:6,1)     &
                        &                 *DBLE(Transfo(j,i,P,it)) &
                        &                 *DBLE(Transfo(k,i,P,it))
                    enddo
                enddo
            enddo
        enddo
    enddo
  end function QPAlignment

  subroutine PrintNumberParities
    !---------------------------------------------------------------------------
    ! Count and print the number parities of the HFB vacuum we are currently in.
    !
    ! The number parity is the multiplicity of the number of 1 eigenvalues of
    ! the RhoHFB matrix.
    !---------------------------------------------------------------------------
    integer :: P,it,j,i

    1 format(" Number Parities (from counting holes, don't trust this easily)")
   11 format(' Number Parities (from counting particles)') 
  111 format(' precision = 1d-5'   )
  112 format(' precision = 1d-6'   ) 
 1111 format(' Pfaffian of the HFB hamiltonian in the subblocks')
    3 format(' P =+1', 17x, i5,5x,i5)
    2 format(' P =-1', 17x, i5,5x,i5)
    4 format(' P = 0', 17x, i5,5x,i5)
    
    5 format(' Sx=+i', 17x, i5,5x,i5)
    6 format(' Sx=-i', 17x, i5,5x,i5)

    NP = 0
    do it=1,Iindex
      do P=1,Pindex
        NP(P,it) = blocksizes(P,it)/2
        do j=1,Blocksizes(P,it)    
          !--------------------------------------------------------------------
          ! For reasons I do not completely understand yet, the detection
          ! of number parity is less dependent on the method used when counting
          ! the total number of holes rather than the total number of particles.
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The cutoff precision also seems to be 'optimal' in the following 
          ! sense: pairs in the canonical basis have the same occupation up to
          ! about 10**-9 in some cases (but often better). Test case was 
          ! one neutron qp state in Th223.
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! 
          !--------------------------------------------------------------------
          if(abs(Occupations(j,P,it)).lt.1d-8) NP(P,it) = NP(P,it) - 1
          !if(abs(Occupations(j,P,it) - 1.0_dp).lt.1d-5) NP(P,it) = NP(P,it) + 1
          !if(abs(Occupations(j,P,it) - 2.0_dp).lt.1d-5) NP(P,it) = NP(P,it) + 1
        enddo
        if(TRC) NP(P,it) = NP(P,it)*2
        enddo
    enddo
    
    print *
    print 1
    if(PC) then
      print 2, NP(1,:)
      print 3, NP(2,:)
    else
      print 4, NP(1,:)
    endif

    NP = 0
    do it=1,Iindex
      do P=1,Pindex
        NP(P,it) = 0
        do j=1,Blocksizes(P,it) 
            !-------------------------------------------------------------------
            ! For completeness' sake, MOCCa also prints the number parities
            ! by counting eigenvalues close to 1
            !-------------------------------------------------------------------
            if(TRC) then
                if(abs(Occupations(j,P,it) - 2.0_dp).lt.1d-5) then
                  NP(P,it) = NP(P,it) + 1
                endif
            else
                if(abs(Occupations(j,P,it) - 1.0_dp).lt.1d-5) then
                  NP(P,it) = NP(P,it) + 1
                endif
            endif
        enddo
        if(TRC) NP(P,it) = NP(P,it)*2
      enddo
    enddo
    
    print *
    print 11
    print 111
    if(PC) then
      print 2, NP(1,:)
      print 3, NP(2,:)
    else
      print 4, NP(1,:)
    endif

    NP = 0
    do it=1,Iindex
      do P=1,Pindex
        NP(P,it) = 0
        do j=1,Blocksizes(P,it)    
          !---------------------------------------------------------------------
          ! For completeness' sake, MOCCa also prints the number parities
          ! by counting eigenvalues close to 1
          !---------------------------------------------------------------------
          if(TRC) then
              if(abs(Occupations(j,P,it) - 2.0_dp).lt.1d-6) then
                NP(P,it) = NP(P,it) + 1
              endif
          else
              if(abs(Occupations(j,P,it) - 1.0_dp).lt.1d-6) then
                NP(P,it) = NP(P,it) + 1
              endif
          endif
        enddo
        if(TRC) NP(P,it) = NP(P,it)*2
      enddo
    enddo
    
    print 112
    if(PC) then
      print 2, NP(1,:)
      print 3, NP(2,:)
    else
      print 4, NP(1,:)
    endif

    ! count number of levels with occuation printed as 1.000
    ! in the table of canonical basis states
    NP = 0
    do it=1,Iindex
      do P=1,Pindex
        NP(P,it) = 0
        do j=1,Blocksizes(P,it)    
          if(TRC) then
              if(abs(Occupations(j,P,it) - 2.0_dp).lt.1d-4) then
                NP(P,it) = NP(P,it) + 1
              endif
          else
              if(Occupations(j,P,it).gt.0.99994999999_dp) then
                NP(P,it) = NP(P,it) + 1
              endif
          endif
        enddo
        if(TRC) NP(P,it) = NP(P,it)*2
      enddo
    enddo
    
    if(PC) then
      print '(/," number parities at 5.e-5 ",6i4,/)', &
        & NP(1,1),NP(2,1),NP(1,1)+NP(2,1), &
        & NP(1,2),NP(2,2),NP(1,2)+NP(2,2)
    else
      print '(/," number parities at 5.e-5 ",6i4,/)', &
        & NP(1,1),NP(1,2)
    endif

  !-----------------------------------------------------------------------------
  print *
  print 1111
  if(PC) then
    print 2, int(pf(1,:)/abs(pf(1,:)))
    print 3, int(pf(2,:)/abs(pf(2,:)))
  else
    print 4, int(pf(1,:)/abs(pf(1,:)))
  endif

  !-----------------------------------------------------------------------------
  if(PC) then
    print '(/," Swapped qps : ",4i5,/)',QPswapped(1,1,1),QPswapped(2,1,1), &
         & QPswapped(1,2,1),QPswapped(2,2,1)
  else

  endif

  end subroutine PrintNumberParities

  subroutine CalcNumberParities
    !--------------------------------------------------------------------------
    ! count number of levels with occuation printed as 1.000
    ! in the table of canonical basis states
    !--------------------------------------------------------------------------
    integer :: P,it,j,i
 
    NP = 0
    do it=1,Iindex
      do P=1,Pindex
        NP(P,it) = 0
        do j=1,Blocksizes(P,it)
          if(TRC) then
              if(abs(Occupations(j,P,it) - 2.0_dp).lt.1d-4) then
                NP(P,it) = NP(P,it) + 1
              endif
          else
              if(Occupations(j,P,it).gt.0.99994999999_dp) then
                NP(P,it) = NP(P,it) + 1
              endif
          endif
        enddo
        if(TRC) NP(P,it) = NP(P,it)*2
      enddo
    enddo
 
  end subroutine CalcNumberParities

  subroutine EstimateNumberParitiesVac(P,it)
        !-----------------------------------------------------------------------
        ! Estimate number parity from V matrices in the QP basis.
        ! Stripped-down ConstructRHOHFB without storing it combined with
        ! and stripped-down DiagonaliseRhoHFB without storing eigenvalues etc.
        !-----------------------------------------------------------------------
        ! NOT SUITED FOR TIME SIMPLEX BREAKING CALCULATIONS!
        !-----------------------------------------------------------------------
        integer, intent(in)            :: P, it
        real(KIND=dp), allocatable     :: Work(:), Eigenvalues(:)
        real(KIND=dp), allocatable     :: Temp(:,:)
        real(KIND=dp), allocatable     :: Eigenvectors(:,:)
        complex(KIND=dp), allocatable  :: RhoTemp(:,:)
        integer                        :: N, i, j, k, Nmax, ii, iii, ifail

        Nmax = maxval(Blocksizes)
         if(.not. allocated(Work)) then
            allocate(Eigenvectors(Nmax,Nmax)) ; Eigenvectors = 0.0_dp
            allocate(Work(Nmax*5))            ; Work = 0.0_dp
            allocate(Temp(Nmax,Nmax))         ; Temp = 0.0_dp
            allocate(Eigenvalues(Nmax))       ; Eigenvalues = 0.0_dp
            allocate(RhoTemp(Nmax,Nmax))      ; RhoTemp = 0.0_dp
        endif
        
  !     do it=1,Iindex
  !     do P =1,Pindex
          RhoTemp = 0.0_dp
          Temp    = 0.0_dp
          N       = blocksizes(P,it)
          !---------------------------------------------------------------
          ! rho = V^* V^T (which in general is complex)
          !---------------------------------------------------------------
          do i=1,N
          do j=1,N
          do k=1,N
            RhoTemp(i,j) = RhoTemp(i,j)                              &
            &              +  Conjg(V(i,HFBColumns(k,P,it),P,it))    &
            &                     * V(j,HFBColumns(k,P,it),P,it)
          enddo
          enddo
          enddo
          Temp(1:N,1:N) = DBLE(RhoTemp(1:N,1:N))

          !---------------------------------------------------------------
          ! when the particle number of species "it" is zero, the canonical
          ! basis cannot be safely obtained by diagonalizing rho, as this
          ! matrix is formally zero and numerically numerical noise.
          ! Needed for trapped neutron droplets.
          !---------------------------------------------------------------
          if ( Particles(it) .lt. 0.000000001_dp ) then
            Eigenvalues(1:N) = 0.0_dp
            NPVac(P,it) = 0
            return
          endif

          if(SC) then
            do i=1,N
              ii  = blockindices(i,P,it)
              iii = mod(ii-1,nwt)+1
              if( ii .ne. iii .or. HFBasis(iii)%GetSignature().eq.-1) then
                Temp(i,i) = Temp(i,i) - 2.0_dp
              endif
            enddo
          endif
          call diagoncr8(temp,Nmax,N,Eigenvectors,Eigenvalues, Work,       &
          &                                              'DiagRho   ',ifail)
          if(SC) where( Eigenvalues .lt. -0.1_dp) Eigenvalues = Eigenvalues + 2
          !-----------------------------------------------------------------
          ! Ensure that occupations are positive semi-definite
          !-----------------------------------------------------------------
          where ( Eigenvalues .lt. 0.0_dp ) Eigenvalues = 0.0_dp

          !-----------------------------------------------------------------
          ! now calculate number parity
          !-----------------------------------------------------------------
          NPVac(P,it) = 0
          do i=1,N
            if(TRC) then
              if( abs(Eigenvalues(i) - 2.0_dp).lt.1d-4) then
                NPVac(P,it) = NPVac(P,it) + 1
              endif
            else
              if( Eigenvalues(i).gt.0.99994999999_dp) then
                NPVac(P,it) = NPVac(P,it) + 1
              endif
            endif
          enddo
          if(TRC) NPVac(P,it) = NPVac(P,it)*2

          !-----------------------------------------------------------------
          ! diagnostic printing (usually commented out)
          !-----------------------------------------------------------------
          ! print '(" EstimateNumberParitiesVac: ",2i2,5i5)',it,P, & 
          !         & HFBNumberParity(P,it),NPVac(P,it),(-1)**NPVac(P,it)

   !    enddo
   !    enddo

  end subroutine EstimateNumberParitiesVac

  function Dispersion()
    !---------------------------------------------------------------------------
    ! Calculates 2*Tr(Rho * (1 - Rho))
    real(KIND=dp) :: Dispersion(2)
    real(KIND=dp) :: Chi(HFBSize, HFBSize,2,2)
    integer       :: i,j,k,P,it,N

    Dispersion = 0.0_dp
    do it=1,Iindex
      do P=1,PIndex
        N = blocksizes(P,it)
        do j=1,N
          do i=1,N
            Chi(i,j,P,it) = - RhoHFB(i,j,P,it)
          enddo
          Chi(j,j,P,it) = Chi(j,j,P,it) + 1
        enddo

        do j=1,N
          do i = 1,N
            Dispersion(it)  = Dispersion(it) + RhoHFB(i,j,P,it) * Chi(j,i,P,it)
          enddo
        enddo
      enddo
    enddo
    Dispersion = 2*Dispersion
  end function Dispersion

  subroutine PrintHFBconvergence
    !----------------------------------------------------------------------------
    ! Prints out some convergence info on the HFB subproblem.
    !   a) size of the change in rho, kappa
    !   b) size of rho*rho - rho + kapppa * kappa^T
    !   c) size of rho*kappa - kappa * rho
    ! All of these should be small at convergence.
    !----------------------------------------------------------------------------

    real(KIND=dp) :: test1(2,2), test2(2,2)
    real(KIND=dp), allocatable :: A(:,:)
    integer       :: N,i,j, P,it
    

    1 format (' HFB convergenceo:   (N,-)   (N,+)    (P,-)    (P,+)')
    2 format ('  Change  dRho   = ',  4es9.2,/,  &
       &      '          dKappa = ',  4es9.2)
    3 format ('  r^2-r+k*k^T    = ',  4es9.2)
    4 format ('  r*k-k*r        = ',  4es9.2)

    print *
    print 1
    print 2, drho, dkappa

    do it=1,Iindex
        do P=1,Pindex
            N = blocksizes(P,it)
            allocate(A(N,N))
            A = matmul(rhoHFB(1:N,1:N,P,it), rhoHFB(1:N,1:N,P,it)) - rhoHFB(1:N,1:N,P,it)
            A = A + matmul(kappaHFB(1:N,1:N,P,it), transpose(kappaHFB(1:N,1:N,P,it)))

            test1(P,it) = sum(A**2)
            
            A = matmul(rhoHFB(1:N,1:N,P,it), kappaHFB(1:N,1:N,P,it)) 
            A = A - matmul(kappaHFB(1:N,1:N,P,it),rhoHFB(1:N,1:N,P,it))
            test2(P,it) = sum(A**2)
            deallocate(A)
        enddo
    enddo

    print 3, test1
    print 4, test2

  end subroutine PrintHFBconvergence

  subroutine diagoncr8 (a,ndim,n,v,d,wd, callrout,ifail)
  !.............................................................................
  !  diagonalization of a real symmetric matrix a(i,j)                         .
  !     input : a  n*n matrix           (with declared dimensions ndim*ndim)   .
  !             a(i,k) * v(k,j) = v(i,k) * d(k)                                .
  !     output: v block of eigenvectors (with declared dimensions ndim*ndim)   .
  !             d eigenvalues in ascending order                               .
  !             wd working array                                               .
  !             -------------------------------------------------              .
  !             the content of a is lost (actualy, a(i,i) = d(i))              .
  !                THIS IS A LIE, A IS NOT TOUCHED AND DEFINITELY              . 
  !                         a(i,i) ! = d(i)               .                    .
  !             -------------------------------------------------              .
  !             a(i,j) = v(i,k) * d(k) * v(j,k)                                .
  !                      v is an orthogonal matrix : v(i,k)*v(j,k) = delta_ij  .
  !..............................................................................
      integer,parameter        :: jstop=30
      real(KIND=dp), parameter :: eps=9.0d-12,epsd=1.0d-16,tol=1.0d-36
      real(KIND=dp), parameter :: zero=0.0d0,one=1.0d0,two=2.0d0
      character(len=10), intent(in)  :: callrout

      real(KIND=dp), intent(inout) :: a(ndim,ndim),v(ndim,ndim),d(ndim),wd(ndim)
      integer, intent(in)          :: ndim,n
      integer, intent(inout)       :: ifail

      integer       :: i,j,ml, m1,m,n1, k, ii, j1, l
      real(KIND=dp) :: p, scale, g,r, s, c, hh, b, f, h

      ifail = 0

      v(1,1) = one
      d(1)   = a(1,1)
      if (n.le.1) return

      do i=1,ndim*ndim
        v(i,1) = a(i,1)
      enddo

      do 9 ii=2,n
        i     = n+2-ii
        l     = i-1
        h     = zero
        scale = zero
        if (l.eq.1) go to 100
        do k=1,l
          scale = scale + abs(v(i,k))
        enddo
        if (scale.gt.tol) go to 3
  100   continue
        wd(i) = v(i,l)
        d(i) = h
        go to 9
    3   do k=1,l
          v(i,k) = v(i,k)/scale
          h      = h + v(i,k)**2
        enddo
        f      = v(i,l)
        g      =-sign(sqrt(h),f)
        wd(i)  = g*scale
        h      = h-f*g
        v(i,l) = f-g
        f      = zero
        do j=1,l
          v(j,i) = v(i,j)/(h*scale)
          g = zero
          do k=1,j
            g = g + v(j,k)*v(i,k)
          enddo
          j1 = j + 1
          if (j1.le.l) then
            do k=j1,l
              g = g + v(k,j)*v(i,k)
            enddo
          endif
          wd(j) = g/h
          f     = f + wd(j)*v(i,j)
        enddo
        hh = f/(h+h)
        do j=1,l
          f     = v(i,j)
          g     = wd(j) - hh*f
          wd(j) = g
          do k=1,j
            v(j,k) = v(j,k) - f*wd(k) - g*v(i,k)
          enddo
        enddo
        do k=1,l
          v(i,k) = scale * v(i,k)
        enddo
        d(i)=h
    9 continue

      d(1)  = zero
      wd(1) = zero
      do i=1,n
        l=i-1
        if ((abs(d(i)).ge.epsd).and.(l.ne.0)) then
          do j=1,l
            g = zero
            do k=1,l
              g = g + v(i,k)*v(k,j)
            enddo
            do k=1,l
              v(k,j) = v(k,j) - g*v(k,i)
            enddo
          enddo
        endif
        d(i)   = v(i,i)
        v(i,i) = one
        if (l.ne.0) then
          do j=1,l
            v(i,j) = zero
            v(j,i) = zero
          enddo
        endif
      enddo

      do i=2,n
        wd(i-1) = wd(i)
      enddo

      wd(n) = zero
      b     = zero
      f     = zero
      do 212 l=1,n
        j = 0
        h = eps * ( abs(d(l)) + abs(wd(l)) )
        if (b.lt.h) b = h
        m = l - 1
  202   m = m + 1
        if (m.gt.n) go to 203
        if (abs(wd(m))-b) 203,203,202
  203   continue
        if (m.eq.l) go to 211
  204   continue
        if (j.eq.jstop) then
          ifail = 1
          return
        endif
        j = j + 1
        p = (d(l+1)-d(l))/(two*wd(l))
        r = sqrt(p*p+one)
        h =  d(l)-wd(l)/(p+sign(r,p))
        do i=l,n
          d(i) = d(i) - h
        enddo
        f  = f+h
        p  = d(m)
        c  = one
        s  = zero
        m1 = m  - 1
        ml = m1 + l
        do ii=l,m1
          i = ml - ii
          g = c*wd(i)
          h = c*p
          if (abs(p).ge.abs(wd(i))) then
            c       = wd(i)/p
            r       = sqrt(c*c+one)
            wd(i+1) = s*p*r
            s       = c/r
            c       = one/r
          else
            c       = p/wd(i)
            r       = sqrt(c*c+one)
            wd(i+1) = s*wd(i)*r
            s       = one/r
            c       = c/r
          endif
          p      = c*d(i) - s*g
          d(i+1) = h + s*(c*g+s*d(i))
          do k=1,n
            h        = v(k,i+1)
            v(k,i+1) = s*v(k,i) + c*h
            v(k,i)   = c*v(k,i) - s*h
          enddo
        enddo
        wd(l) = s*p
        d(l)  = c*p
        if (abs(wd(l))-b) 211,211,204
  211   continue
        d(l) = d(l)+f
  212 continue

      n1 = n-1
      do i=1,n1
        k  = i
        p  = d(i)
        ii = i + 1
        do j=ii,n
          if (d(j).lt.p) then
            k = j
            p = d(j)
          endif
        enddo
        if (k.ne.i) then
          d(k) = d(i)
          d(i) = p
          do j=1,n
            p      = v(j,i)
            v(j,i) = v(j,k)
            v(j,k) = p
          enddo
        endif
      enddo

      return
      end subroutine diagoncr8

end module HFB
