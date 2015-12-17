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
  real(KIND=dp),parameter :: HFBNumCut=1d-10
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
   integer, allocatable :: HFBColumns(:,:,:)
  !-----------------------------------------------------------------------------
  ! Logical. Check all kinds of relations that should hold for correct HFB
  ! calculations. Only to be used for debugging purposes.
  logical,parameter :: HFBCheck=.false.
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
  real(KIND=dp) :: HFBGauge = 0.0_dp
  !----------------------------------------------------------------------------
  !Quasiparticle Energies: note that all of these get stored, even the occupied
  ! ones.
  ! Also the signatures of the quasiparticles are stored, since they are quite 
  ! a bore to get everytime they are needed.
  real(KIND=dp), allocatable :: QuasiEnergies(:,:,:), QuasiSignatures(:,:,:)
  !-----------------------------------------------------------------------------
  ! Maximum number of iterations that the Fermi solver can take
  integer :: HFBIter = 1
  !-----------------------------------------------------------------------------
  ! Logical that traces if quasiparticles excitations are defined in the HF or
  ! in the canonical basis.
  logical :: QPinHFBasis=.false.
  !-----------------------------------------------------------------------------
  ! Indices of the quasiparticle states that we want to block.
  ! NOTE THAT THESE ARE IN THE HFBasis, as this is truely the only
  ! sensible thing to do.
  ! The rest of the arrays take the quantum numbers. Only QPBlockind takes 
  ! the indices of the HF wavefunctions in their respective parity/isospin
  ! blocks.
  integer,allocatable :: QPExcitations(:)
  integer,allocatable :: QPParities(:)
  integer,allocatable :: QPIsospins(:)
  integer,allocatable :: QPSignatures(:)
  integer,allocatable :: QPblockind(:)
  !-----------------------------------------------------------------------------
  ! Integer that controls whether the program looks for even or odd HFB states.
  ! This needs to be stored for every parity-isospin block.
  integer :: HFBNumberParity(2,2)=0
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
  ! 'Distance in energy' for the printing of quasiparticle excitations
  real(Kind=dp) :: QPPrintWindow=100.0_dp
  !-----------------------------------------------------------------------------
  ! Whether we want to use Broyden (fast and not 100% guaranteed) or Bisection
  ! (slow and guaranteed) method to solve for the Fermi energy
  character(len=9) :: FermiSolver='Broyden'


contains
  subroutine PrepareHFBModule
    !---------------------------------------------------------------------------
    ! This subroutine nicely prepares the module at the start of the program.
    !---------------------------------------------------------------------------
    integer :: i, P, it, ii
    
    Pindex = 1 ; if(PC) Pindex = 2
    Iindex = 1 ; if(IC) Iindex = 2 

    if(.not.allocated(Blocksizes))   allocate(blocksizes(Pindex,Iindex))
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
      allocate(Occupations(HFbsize,Pindex,Iindex))        ; Occupations    = 0.0_dp
      allocate(QuasiEnergies(2*HFBSize,Pindex,Iindex))    ; QuasiEnergies  = 0.0_dp
      allocate(QuasiSignatures(2*HFBSize,Pindex,Iindex))  ; QuasiSignatures= 0.0_dp
      allocate(HFBColumns(HFBSize,Pindex,Iindex))         ; HFBColumns     = 0
    endif

  end subroutine PrepareHFBModule

  function HFBEnergy(Delta) result(Energy)
    !---------------------------------------------------------------------------
    ! Placeholder function for calculating the HFB energy.
    ! Note that this formula introduces an extra cutoff factor compared to R&S.
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
    !  Sum_{i,j} 0.5**fj *Gamma_{ij}*s*<r,s;r,-s|V|ij>
    !  where Gamma = (1 - 2 Rho) * Kappa
    !---------------------------------------------------------------------------
    use Spinors
    use PairingInteraction
    
    complex(KIND=dp), allocatable, intent(inout) :: Field(:,:,:,:)
    complex(KIND=dp), allocatable, intent(inout) :: FieldLN(:,:,:,:)
    complex(KIND=dp), allocatable, intent(inout) :: Delta(:,:,:,:)
    complex(KIND=dp)                             :: ActionOfPairing(nx,ny,nz)
    complex(KIND=dp),allocatable,save            :: Gamma(:,:,:,:)
    real(KIND=dp) :: maxK
    integer      :: i, it, j, ii, sig1, sig2, jj,k,P, iii, jjj,l
    real(KIND=dp):: Cutoff(2) 
    type(Spinor) :: Temp(2)

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
      do it=1,2
        do P=1,Pindex
          do j=1,HFBSize
            do i=1,HFBSize
              Gamma(i,j,P,it) = KappaHFB(i,j,P,it)
              do k=1,HFBsize
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
    maxK = 0.0_dp
    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          jj        = Blockindices(j,P,it)
          jjj       = mod(jj-1,nwt)+1
          sig1      = HFBasis(jjj)%GetSignature()
          Temp(1)   = HFBasis(jjj)%GetValue()
          Cutoff(1) = PCutoffs(jjj)
          if(Cutoff(1).lt.HFBNumCut) cycle
          if(TRC .and. jj .ne. jjj) then
            sig1 = - sig1
            !PairingInter now takes care of time-reversed problems
            !Temp(1) = TimeReverse(Temp(1))
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
            Temp(2)  = HFBasis(iii)%GetValue()

            !Note that this does automatically include a Time-reversal operator
            ! when appropriate
            ActionOfPairing = GetPairDensity(Temp(1),Temp(2))
            do l=1,nx*ny*nz
              Field(l,1,1,it) = Field(l,1,1,it)     -  Cutoff(1)*Cutoff(2)*      &
              &                   DBLE(KappaHFB(i,j,P,it))*ActionOfPairing(l,1,1)
            enddo
            if(allocated(FieldLN)) then
              do l=1,nx*ny*nz
                FieldLN(l,1,1,it) = FieldLN(l,1,1,it) -  Cutoff(1)*Cutoff(2)*    &
                &                 DBLE(Gamma(i,j,P,it))*ActionOfPairing(l,1,1)
              enddo
            endif          
          enddo
        enddo
      enddo
    enddo
    do it=1,2
      do i=1,nx*ny*nz
        Field(i,1,1,it) = Field(i,1,1,it) * DensityFactor(i,1,1,it)
      enddo
    enddo
    if(allocated(FieldLN)) then
      do i=1,nx*ny*nz
        FieldLN(i,1,1,it) = FieldLN(i,1,1,it) * DensityFactor(i,1,1,it)
      enddo
    endif
  end subroutine HFBPairingField

  subroutine HFBGaps(Delta,DeltaLN,PairingField,PairingFieldLN,Gaps,ConstantGap)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the Delta_{i j} for the HFB model.
    ! The formula (see Ring & Shuck, page 254 eq. (7.41))
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
    
    if(ConstantGap) call stp('Trying to do constant gap pairing in HFB!')
    
    Delta = 0.0_dp ; if(allocated(DeltaLN)) DeltaLN = 0.0_dp

    do it=1,Iindex
      do P=1,Pindex
        do i=1,Blocksizes(P,it)
          ii   = Blockindices(i,P,it)
          iii  = mod(ii-1,nwt)+1
          sig1 = HFBasis(iii)%GetSignature()
          Psi1 = HFBasis(iii)%GetValue()
          
          if(TRC .and. ii.ne.iii) then 
            Psi1 = TimeReverse(Psi1)
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
            if(sig1.ne.-sig2) cycle
            Cutoff(2) = PCutoffs(jjj)
            !if(Cutoff(2).lt.HFBNumCut) cycle

            if(Cutoff(1)*Cutoff(2) .lt. HFBNumCut) cycle

            Psi2 = HFBasis(jjj)%GetValue()
            if(TRC .and. jj.gt.Blocksizes(P,it)/2) then
              Psi2 = TimeReverse(Psi2)
            endif

            do k=1,nx*ny*nz
                TempReal(k,1,1) =                                                  &
                &             - Psi1%Grid(k,1,1,3,1) * Psi2%Grid(k,1,1,1,1)        &
                &             + Psi1%Grid(k,1,1,4,1) * Psi2%Grid(k,1,1,2,1)        &
                &             + Psi1%Grid(k,1,1,1,1) * Psi2%Grid(k,1,1,3,1)        &
                &             - Psi1%Grid(k,1,1,2,1) * Psi2%Grid(k,1,1,4,1)  
            
                ! Attention for the extra minus sign: the bra < i,j | indicates an 
                ! extra complex conjugation.
                TempIm(k,1,1)   =                                                  &
                &             + Psi1%Grid(k,1,1,4,1) * Psi2%Grid(k,1,1,1,1)        &
                &             + Psi1%Grid(k,1,1,3,1) * Psi2%Grid(k,1,1,2,1)        &
                &             - Psi1%Grid(k,1,1,2,1) * Psi2%Grid(k,1,1,3,1)        &
                &             - Psi1%Grid(k,1,1,1,1) * Psi2%Grid(k,1,1,4,1)
            enddo

            ! Do less computations if TSC is conserved...
            if(TSC) then
                Delta(i,j,P,it) =   dv*Cutoff(1)*Cutoff(2)*                    & 
                &            (sum(   TempReal * DBLE(PairingField(:,:,:,it)))  &
                &            -sum(   TempIm   * AIMAG(PairingField(:,:,:,it))))
            else
                Temp = DCMPLX(TempReal,TempIm)
                Delta(i,j,P,it) =                                              &
                &        dv*Cutoff(1)*Cutoff(2)*sum(Temp*PairingField(:,:,:,it))
            endif

            !Delta is antisymmetric
            Delta(j,i,P,it) = - Delta(i,j,P,it)
            if(allocated(DeltaLN)) then
                !Only activate when Lipkin-Nogami is active
                if(TSC) then
                    DeltaLN(i,j,P,it) =   dv*Cutoff(1)*Cutoff(2)*              & 
                    &   (sum(   TempReal * DBLE (PairingFieldLN(:,:,:,it)))    &
                    &   -sum(   TempIm   * aimag(PairingFieldLN(:,:,:,it))))
                else
                    DeltaLN(i,j,P,it) = Cutoff(1)*Cutoff(2)*dv*                &
                    &              sum(Temp * PairingFieldLN(:,:,:,it))
                endif
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
    ! The formula (see Ring & Shuck, page 254 eq. (7.41))
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
            &            -sum(   TempIm   * AIMAG(PairingField(:,:,:,it))))
            Delta(j,i,P,it) = - Delta(i,j,P,it)
   
            !Only valid when TimeSimplex is conserved
            DeltaLN(i,j,P,it) =   dv*Cutoff(1)*Cutoff(2)*              & 
            &   (sum(   TempReal * DBLE (PairingFieldLN(:,:,:,it)))    &
            &   -sum(   TempIm   * aimag(PairingFieldLN(:,:,:,it))))
            DeltaLN(j,i,P,it) =  - DeltaLN(i,j,P,it)
          enddo
        enddo
      enddo
    enddo
  end subroutine HFBGaps_TIMEREV

  function HFBNumberOfParticles(Lambda, Delta, LNLambda) result(N)
  !-----------------------------------------------------------------------------
  ! Diagonalise the HFB Hamiltonian as a function of Lambda to get the number
  ! of particles of species it.
  !-----------------------------------------------------------------------------
  ! This is NOT a pure function, the HFB hamiltonian gets constructed, as well
  ! as the HFB density and anomalous density matrix.
  !-----------------------------------------------------------------------------
    integer                   :: i, it, P, iter
    real(Kind=dp), intent(in) :: Lambda(2),LNLambda(2)
    real(KIND=dp)             :: N(2), N2(2,2)
    complex(KIND=dp), allocatable,intent(in) :: Delta(:,:,:,:)

    call ConstructHFBHamiltonian(Lambda, Delta, LNLambda,HFBGauge)
!     ! Diagonalisation of the HFBHamiltonian: computation of U & V matrices.
     call DiagonaliseHFBHamiltonian()
!     !LNFRACTION=0.001
!     print *, 'Before'
!     do it=1,1
!     do P=1,1
!     print *, QuasiEnergies(1:2*blocksizes(P,it),P,it)
!     enddo
!     enddo
!     print *
!     do it=1,1
!     do P=1,1
!     do i=1,blocksizes(P,it)
!     print *, DBLE(U(i,1:2*blocksizes(P,it),P,it))
!     enddo
!     enddo
!     enddo
!     print *
!     do it=1,1
!     do P=1,1
!     do i=1,blocksizes(P,it)
!     print *, DBLE(V(i,1:2*blocksizes(P,it),P,it))
!     enddo
!     enddo
!     enddo
   
!    call ConstructHFBHamiltonian(Lambda, Delta, LNLambda,HFBGauge + 1.0_dp)
!    ! Diagonalisation of the HFBHamiltonian: computation of U & V matrices.
!    call DiagonaliseHFBHamiltonian()

!    print *
!    print *, 'After'
!    do it=1,1
!        do P=1,1
!            print *, QuasiEnergies(1:2*blocksizes(P,it),P,it)
!        enddo
!    enddo
!    print *
!    do it=1,1
!      do P=1,1
!         do i=1,blocksizes(P,it)
!             print *, DBLE(U(i,1:2*blocksizes(P,it),P,it))
!         enddo
!      enddo
!    enddo
!    print *
!    do it=1,1
!      do P=1,1
!         do i=1,blocksizes(P,it)
!             print *, DBLE(V(i,1:2*blocksizes(P,it),P,it))
!         enddo
!      enddo
!    enddo
!    stop

    ! Construct the generalised density matrix and the anomalous one.
    call ConstructHFBstate()
              
    N = 0.0_dp ; N2 = 0.0_dp
    !Calculating total number of particles, by tracing the density matrix
    ! Sidenote: RHOHFB is hermitian, thus the imaginary parts of RHoHFB(i,i) are
    ! zero anyways.
    call DiagonaliseRhoHFB()
    do it=1,Iindex
      do P=1,Pindex
        do i=1,blocksizes(P,it)
          N(it) = N(it) + real(RhoHFB(i,i,P,it))
        enddo
      enddo
    enddo
  end function HFBNumberofParticles

  subroutine HFBFindFermiEnergyBisection(Fermi,L2,Delta,DeltaLN,Lipkin,Prec)
  !-------------------------------------------------------------------------------
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
  ! 2) This routine should be optimised in the sense that in 99% of cases the
  !    bisection method can be replaced by some Newton solver which would be 
  !    faster. TODO
  !
  !-------------------------------------------------------------------------------

  real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
  complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
  complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
  logical, intent(in)                       :: Lipkin
  real(KIND=dp), intent(in)                 :: Prec
  
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
      InitialBracket(:,1) = InitialBracket(:,1) - 0.1_dp
      InitialBracket(:,2) = InitialBracket(:,2) + 0.1_dp
      
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
    print *, 'L2 at start', L2
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
      ! Make some guess for kappa for this isospin, in order to kickstart
      ! the pairing.
        call writeOutKappa(2,it)
      endif
    enddo

    !Replace the previous step
    Gold  = G
    L2Old = tempL2
 
    if(iter.eq.HFBIter) then  
      call stp('No correct L2 found.')
    endif
  enddo

  end subroutine HFBFindFermiEnergyBisection

  subroutine HFBFindFermiEnergyBroyden(Fermi,LnLambda,Delta,DeltaLN,Lipkin,Prec)
  !---------------------------------------------------------------------------
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
  !---------------------------------------------------------------------------
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
  !---------------------------------------------------------------------------
  ! While it might seem like a good idea to let FermiHistory, LNHistory and
  ! the Jacobian persist across mean-field iterations, I have found out that 
  ! in practice this is a bad idea. After a mean-field iteration, the Fermi 
  ! energy is generally close to the correct one, but the approximation that 
  ! was saved might actually be a very bad one.
  !---------------------------------------------------------------------------
  1 format('Attention, unconverged Fermi solver.'/, &
  &        'Iterations : ', i5,                  /, &
  &        'Particles  : ', 2f12.8)
  2 format('Lipkin parameter deviation: ', 2f12.8)


  real(KIND=dp), intent(inout)              :: Fermi(2), LnLambda(2)
  complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
  complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
  logical, intent(in)                       :: Lipkin
  real(KIND=dp), intent(in)                 :: Prec
  integer                                   :: it, iter, flag(2)
  ! Previous values of the Fermi energy and LNLambda parameter
  real(KIND=dp), save                   :: FermiHistory(2)=100.0_dp
  real(KIND=dp), save                   :: LNhistory(2)=100.0_dp
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

  ! First time, take some guess for the histories
  FermiHistory = Fermi    + 0.01_dp
  LNHistory    = LNLambda - 0.01_dp    
  Particles(1) = Neutrons
  Particles(2) = Protons
  Converged    = .false.

  flag = 0
  !Check were we find ourselves in the phasespace
  N = HFBNumberofParticles(Fermi, Delta, LnLambda ) - Particles
  if(Lipkin) then
    LN   = LNLambda - LNCR8(Delta,DeltaLN, flag)
  else
    ! Make sure the algorithm is converged in LN when LN is not active.
    ! This is only to make sure no contamination occurs and convergence
    ! is not wrongly judged.
    LN   = 0.0_dp
  endif

  ! Evaluating f and g for finite differences.
  NX = HFBNumberOfParticles(FermiHistory, Delta,LNLambda )-Particles               
  if(Lipkin) LNX  = LNLambda - LNCR8(Delta,DeltaLN, flag)

  if(Lipkin) NY = HFBNumberOfParticles(Fermi, Delta,LNHistory )-Particles
  if(Lipkin) LNY  = LNHistory - LNCR8(Delta,DeltaLN, flag)

  ! Construct a finite difference approximation of the Jacobian
  Jacobian(1,1,:)   = (N   -  NX )/(Fermi    - FermiHistory)
  if(Lipkin) then
    Jacobian(1,2,:) = (N   -  NY )/(LNLambda - LNHistory)
    Jacobian(2,1,:) = (LN  -  LNX)/(Fermi    - FermiHistory)
    Jacobian(2,2,:) = (LN  -  LNY)/(LNLambda - LNHistory)
  endif

  do iter=1,HFBIter
     !Invert the Jacobian
      if (Lipkin) then
        det = Jacobian(1,1,:)*Jacobian(2,2,:) - Jacobian(1,2,:)*Jacobian(2,1,:) 
      else
        det = Jacobian(1,1,:)
      endif
    
      !Invert the jacobian
      if(Lipkin) then
       invJ(1,1,:) =   Jacobian(2,2,:)/det 
       invJ(1,2,:) = - Jacobian(1,2,:)/det
       invJ(2,1,:) = - Jacobian(2,1,:)/det
       invJ(2,2,:) =   Jacobian(1,1,:)/det 
      else
       invJ(1,1,:) = 1.0/det
      endif

      !Find a good direction to update in
      do it=1,2
        !Don't update if converged
        if(Converged(it)) cycle
        !Find Update directions
        FermiUpdate(it) = - invJ(1,1,it) * N(it)    
        if(Lipkin) then
          FermiUpdate(it) = FermiUpdate(it)        - invJ(1,2,it) * LN(it)
          if(flag(it).ne.1) then
            LNupdate(it)    = - invJ(2,1,it) * N(it) - invJ(2,2,it) * LN(it)
          else
            ! Don't update LN if the calculation of Lambda_2 got into trouble
            print *, 'Discarded step at iteration ', iter, ' for isospin ', it
            LNUpdate(it)    = 0.0_dp
          endif
          norm(it)        = FermiUpdate(it)**2 + LNUpdate(it)**2
        else
         norm(it)        = FermiUpdate(it)**2
        endif  
      enddo

      FermiUpdate = step * FermiUpdate
      LNUpdate    = step * LNUpdate

      !Replace the history and update
      do it=1,2
        ! Replace the history
        FermiHistory(it) = Fermi(it) 
        Fermi(it)        = Fermi(it)    + FermiUpdate(it)
        if(Lipkin) then
          LNHistory(it) = LNLambda(it)
          LNLambda(it)  = LNLambda(it) + LNUpdate(it) 
        endif
      enddo

      !Recalculate
      flag = 0
      N    = HFBNumberOfParticles(Fermi       ,Delta,LNLambda )  - Particles
      if(Lipkin) LN   = LNLambda - LNCR8(Delta,DeltaLN, flag)

      !Convergence check
      do it=1,2
        if( abs(N(it)).lt.Prec  .and. abs(LN(it)).lt. Prec) Converged(it)=.true.
      enddo
      if(all(converged)) then
        exit
      endif
      !Update the Jacobian
      if(iter.gt.1) then
        do it=1,2
          if(Converged(it)) cycle    
          Jacobian(1,1,it) = Jacobian(1,1,it)  +N(it)*FermiUpdate(it)/norm(it)
          if(Lipkin) then           
            Jacobian(1,2,it) = Jacobian(1,2,it)+N(it) *LNUpdate(it)   /norm(it)
            Jacobian(2,1,it) = Jacobian(2,1,it)+LN(it)*FermiUpdate(it)/norm(it)
            Jacobian(2,2,it) = Jacobian(2,2,it)+LN(it)*LNUpdate(it)   /norm(it)
          endif
        enddo
      endif

      !Print a warning if not converged
      if(iter.eq.HFBIter) then
          print 1, HFBIter, N + Particles
          print 2, abs(LN)
      endif
  enddo
  end subroutine HFBFindFermiEnergyBroyden

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
  real(KIND=dp), intent(in)                :: lambda(2), LNLambda(2), Gauge
  integer                                  :: i,j, it, ii,iii,P
  complex(KIND=dp), allocatable,intent(in) :: Delta(:,:,:,:)
  
  !----------------------------------------------------------------------------- 
  HFBHamil = 0.0_dp  
  do it=1,Iindex
    do P=1,Pindex
      do i=1,blocksizes(P,it)
        ii  = blockindices(i,P,it)
        iii = mod(ii-1,nwt)+1
        !-----------------------------------------------------------------------
        ! Diagonal contribution:
        ! epsilon - lambda  - 2 * Lambda_2 ( 1 - 2 * Rho)                    
        HFBHamil(i,i,P,it) = HFBasis(iii)%GetEnergy() -  Lambda(it)            &
        &                     - 2*(1.0_dp-LNFraction)* LNLambda(it)                            
        !------------------------------------------------------------------------
        ! Offdiagonal
        do j=1,blocksizes(P,it)
          HFBHamil(i,j,P,it) = HFBHamil(i,j,P,it)                              &
          &       +  (1.0_dp-LNFraction)*4*LNLambda(it)*   OldRhoHFB(i,j,P,it) & 
          &       +  Gauge                             *   OldRhoHFB(i,j,P,it)
        enddo
        do j=1,blocksizes(P,it)
          !---------------------------------------------------------------------
          ! - h^*
          HFBHamil(i+blocksizes(P,it),j+blocksizes(P,it),P,it)=                &
          &                                      -conjg(HFBHamil(i,j,P,it))
        enddo
        !-----------------------------------------------------------------------
        do j=1,blocksizes(P,it)
          ! Delta - 2 * Lambda_2 * Kappa
          HFBHamil(i,j + blocksizes(P,it),P,it) = Delta(i,j,P,it)              &
          &            - LNFraction*4*LNLambda(it)*      OldKappaHFB(i,j,P,it) &
          &            + Gauge                    *      OldKappaHFB(i,j,P,it) 
        enddo
      enddo
    enddo
  enddo

  if(all(HFBHamil.eq.0.0_dp)) call stp('HFBHamiltonian completely zero!')
  end subroutine ConstructHFBHamiltonian  

  subroutine DiagonaliseHFBHamiltonian!_DSYEVR
    !--------------------------------------------
    !
    !http://www.netlib.org/lapack/double/dsyevr.f
    !--------------------------------------------
    integer       :: i,j,k,it, P,jj,ii,jjj,S,Sig1,iii,sig2
    integer,save  :: INFO, LIWORK, LWORK, M, sizeB, N, count=0
    real(KIND=dp) :: ABSTOL, VL, VU
    real(KIND=dp), allocatable,save :: Temp(:,:), Eigenvectors(:,:),Eigenvalues(:)
    real(KIND=dp), allocatable,save :: WORK(:)
    integer, allocatable, save      :: ISUPPZ(:), IWORK(:)

    !---------------------------------------------------------------------------
    ! Preliminary work.
    if(.not.allocated(WORK)) then
        INFO = 0
        sizeB= maxval(blocksizes)
        allocate(Temp(2*sizeB,2*sizeB))         ; Temp         = 0.0_dp
        allocate(Eigenvectors(2*sizeB,2*sizeB)) ; Eigenvectors = 0.0_dp
        allocate(Eigenvalues(2*sizeB))          ; Eigenvalues  = 0.0_dp  
        allocate(ISUPPZ(4*sizeB))               ; ISUPPZ       = 0
        !Workspace query
        LWORK=-1 ; LIWORK =-1
        allocate(WORK(10)) ; allocate(IWORK(10))

        call DSYEVR( 'V', 'A', 'U', sizeB, Temp, sizeB, 0.0, 0.0, 0, 0,      &
        &                0.0, M, Eigenvalues, Eigenvectors, sizeB, ISUPPZ, &
        &                WORK, LWORK, IWORK, LIWORK, INFO )

        LWORK = int(WORK(1)) ; LIWORK = int(IWORK(1))
        deallocate(WORK, IWORK)
        allocate(WORK(LWORK), IWORK(LIWORK))
    endif

    U = 0.0_dp ; V=0.0_dp
    do it=1,Iindex
      do P=1,Pindex
        
        Temp = 0.0_dp
        Eigenvalues = 0.0_dp
        Eigenvectors= 0.0_dp
        N = blocksizes(P,it)
        do j=1,N/2
            do i=1,N/2
                Temp(i,j)         = DBLE(HFBHamil(i,j,P,it))
                Temp(i+N/2,j+N/2) = DBLE(HFBHamil(i+3*N/2,j+3*N/2,P,it))
                Temp(i,j+N/2)     = DBLE(HFBHamil(i,j+3*N/2,P,it))
            enddo
        enddo
        do j=N+1,2*N
            do i=N+1,2*N
                Temp(i,j)          = DBLE(HFBHamil(i-N/2,j-N/2,P,it))
            enddo
        enddo
        !------------------------------------------------------------------------
        ! We store all possible eigenvectors and later make the proper selection.
        
        !First signature block
        INFO = 0
        call DSYEVR( 'V', 'A', 'U', N, Temp(1:N,1:N), N, 0.0, 0.0, 0, 0,       &
        &      0.0, M, Eigenvalues(1:N), Eigenvectors(1:N,1:N), N, ISUPPZ,  &
        &      WORK, LWORK, IWORK, LIWORK, INFO )
        if(INFO .ne. 0) call stp('DSYERV failed.', 'INFO', INFO)
        
        U(      1:N/2,    1:N  ,P,it) = Eigenvectors(       1:N/2   ,  1:N  )
        V(  N/2+1:N  ,    1:N  ,P,it) = Eigenvectors(   N/2+1:N     ,  1:N  )

        QuasiEnergies(1:N,P,it) = EigenValues(1:N)

        !Second signature block
        if(.not.TRC) then
            INFO = 0
            call DSYEVR( 'V', 'A', 'U', N, Temp(N+1:2*N,N+1:2*N), N, 0.0, 0.0, 0, 0,&
            &      0.0, M, Eigenvalues(N+1:2*N), Eigenvectors(N+1:2*N,N+1:2*N), N,  &
            &      ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
            if(INFO .ne. 0) call stp('DSYERV failed.', 'INFO', INFO)
            U(  N/2+1:N  ,N+1:2*N  ,P,it) = Eigenvectors(     N+1:3*N/2 ,N+1:2*N)
            V(      1:N/2,N+1:2*N  ,P,it) = Eigenvectors(   3*N/2+1:2*N   ,N+1:2*N)
            QuasiEnergies(N+1:2*N,P,it)   = EigenValues(N+1:2*N)
        else
            U(  N/2+1:N  ,N+1:2*N  ,P,it) = Eigenvectors(       1:N/2   ,  1:N  )
            V(      1:N/2,N+1:2*N  ,P,it) = Eigenvectors(   N/2+1:N     ,  1:N  )
            QuasiEnergies(N+1:2*N,  P,it) = QuasiEnergies(1: N,P,it)
        endif
        
      enddo
    enddo

    call InsertionSortQPEnergies
  end subroutine DiagonaliseHFBHamiltonian!_DSYEVR
 
!   subroutine DiagonaliseHFBHamiltonian_ZHEEVR
!     !---------------------------------------------------------------------------
!     ! Alternative subroutine for diagonalising the HFB Hamiltonian.
!     !---------------------------------------------------------------------------
!     ! Subroutine that solves the HFB eigenvalue problem 
!     !
!     !      ( U_k )        ( U_k )
!     !   H  (     )  = Ek  (     )
!     !      ( V_k )        ( V_k )
!     !
!     ! And stores the result in matrices U & V.
!     !---------------------------------------------------------------------------
!     ! NOTE 
!     ! 1) Extra documentation for the diagonalisation routine 
!     !    http://www.netlib.org/lapack/explore-html/d9/dd2/zheevr_8f.html
!     !---------------------------------------------------------------------------
!     integer, save                       :: size, m, lwork,lrwork,liwork
!     complex(KIND=dp), allocatable,save       :: Eigenvectors(:,:)
!     complex(KIND=dp), allocatable,save       :: Temp(:,:)
!     real(KIND=dp) , allocatable,save         :: Eigenvalues(:)
!     complex(KIND=dp), allocatable,save       :: Work(:), rwork(:)
!     integer, allocatable,save                :: iwork(:), isuppz(:)
!     integer                             :: Succes,i,j,jj,it,P, N, jjj, S
!     real(KIND=dp)                       :: Sig
!     real(KIND=dp)                       :: time0=0,time1=0,timetot=0
    
!     size = maxval(blocksizes)
!     !---------------------------------------------------------------------------
!     ! Preliminary work.
!     if(.not.allocated(WORK)) then
!         !Allocate several arrays for the LAPACK routine
!         allocate(Temp(2*size, 2*size))                   ; Temp = 0.0_dp
!         allocate(Eigenvectors(2*HFBSize, 2*HFBsize))     ; Eigenvectors= 0.0_dp
!         allocate(Eigenvalues(2*size))                    ; EigenValues=0.0_dp
!         allocate(Isuppz(4*size))                         ; ISUPPZ=0  
!         ! Do a preliminary call to the LAPACK routine to find the optimum
!         ! WORK sizes.
!         LWORK = -1 ; LRWORK = -1 ; LIWORK = -1
!         allocate(work(1), rwork(1), iwork(1))
!         call ZHEEVR('V', 'A', 'U', 2*size, Temp,                            &
!         &            2*size, 0.0_dp , 0.0_dp,0,0,                           &
!         &            0.0_dp, 2*size, Eigenvalues, Eigenvectors, 2*HFBSize,  &
!         &            isuppz,work,lwork,rwork,                               &
!         &            lrwork, iwork, liwork, Succes)

!         lwork = ceiling(DBLE(WORK(1))) ; lrwork = ceiling(DBLE(RWORK(1))) 
!         liwork =IWORK(1)
!         deallocate(WORK, RWORK, IWORK)
!         allocate(WORK(lwork), RWORK(lrwork), IWORK(liwork))
!     endif
!     !----------------------------------------------------------------------------
!     U = 0.0_dp ; V=0.0_dp
!     do it=1,Iindex
!       do P=1,Pindex
!         Temp = 0.0_dp
!         N = blocksizes(P,it)
!         do j=1,2*N
!           do i=1,2*N
!             Temp(i,j) = HFBHamil(i,j,P,it)
!           enddo
!         enddo
!         !------------------------------------------------------------------------
!         ! Temporary trick to separate signatures.
        
!         if(SC) then
!           do j=1,N
!             jj = blockindices(j,P,it)
!             jjj= mod(jj-1,nwt)+1
!             S = HFBasis(jjj)%GetSignature()
!             if(jjj .ne. jj) S = -S
!             Temp(j,j)     = Temp(j,j)     + 1000*S
!             Temp(j+N,j+N) = Temp(j+N,j+N) - 1000*S
!           enddo
!         endif
!         !------------------------------------------------------------------------

!         Succes = 0
!         m = N
!         call ZHEEVR('V', 'A', 'U', 2*N , Temp(1:2*N,1:2*N),2*N, 0.0_dp,          &
!           &          0.0_dp,0,0,0.0_dp, 2*N, Eigenvalues(1:2*N),                 &
!           &          Eigenvectors(1:2*N,1:2*N),2*N, isuppz, work, lwork, rwork,  &
!           &          lrwork,iwork, liwork,Succes)
!         if(Succes.ne.0) then  
!           call stp("Error in diagonalising the HFB Hamiltonian.",                &
!           &        "ZHEEVR Errorcode", Succes)
!         endif
!         !-------------------------------------------------------------------------
!         ! We store all possible eigenvectors and later make the proper selection.       
!         U(1:N,1:2*N,P,it) = Eigenvectors(  1:N    ,1:2*N)
!         V(1:N,1:2*N,P,it) = Eigenvectors(  N+1:2*N,1:2*N)
!         QuasiEnergies(1:2*N,P,it) = EigenValues(1:2*N)
        
!         !-------------------------------------------------------------------------
!         ! Calculate the quasiparticle signatures for further use
!         QuasiSignatures(:,P,it)=0.0_dp
!         do i=1,2*N
!           do j=1,N
!             jj  = blockindices(j,P,it)
!             jjj = mod(jj-1,nwt)+1
!             Sig = HFBasis(jjj)%GetSignatureR()
!             if(jj.ne.jjj) Sig = -Sig
!             !-----------------------------------------------------
!             ! The U components have the same signature
!             QuasiSignatures(i,P,it) = QuasiSignatures(i,P,it)    & 
!             &                       + sig * abs(U(j,i,P,it))**2
!             !-----------------------------------------------------
!             ! The V components have opposite signature, thus 
!             ! opposite sign.
!             QuasiSignatures(i,P,it) = QuasiSignatures(i,P,it)    &
!             &                       - sig * abs(V(j,i,P,it))**2
!           enddo
!         enddo
!         !------------------------------------------------------------------------
!         ! Extra for temporary trick
!         if(SC) then
!           QuasiEnergies(1:2*N, P, it) = QuasiEnergies(1:2*N, P, it) &
!           &                         - 1000*QuasiSignatures(1:2*N,P,it)
!         endif       
!       enddo
!     enddo

!     if(SC) call InsertionSortQPEnergies()

!   end subroutine DiagonaliseHFBHamiltonian_ZHEEVR

subroutine InsertionSortQPEnergies
    !-----------------------------------------------------------------------
    ! Sort the QPenergies, since the trick to separate in 
    ! subroutine DiagonaliseHFBHamilaltonian does not conserve the order.
    ! We also sort the U and V matrices, as well as the Quasisignatures.
    !-----------------------------------------------------------------------
    integer :: i,j, it, P
    real(KIND=dp) :: TempE,TempSig
    complex(KIND=dp) :: Temp(HFBSize)

    do it=1,Iindex
        do P=1,Pindex
            do i=1,2*blocksizes(P,it)-1
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
                  Temp = U(:,j,P,it)
                  U(:,j,P,it)   = U(:,j-1,P,it)
                  U(:,j-1,P,it) = Temp

                  Temp = V(:,j,P,it)
                  V(:,j,P,it)   = V(:,j-1,P,it)
                  V(:,j-1,P,it) = Temp

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

  subroutine ConstructHFBState()
  !-----------------------------------------------------------------------------
  ! High-level routine that constructs the generalised HFB density matrix and 
  ! the anomalous density Kappa.
  !
  ! rho   = V^* V^T
  ! kappa = V^* U^T
  !
  ! However, this concerns the 'text' U and V matrices, and we still need to 
  ! select the qusiparticle states we want to occupy, i.e. which columns to take.
  ! This is ofcourse dependent on the number parity.
  ! To this end, we construct and diagonalise the HFB density multiple times,
  ! checking whether the number parity is the one we want. 
  !
  ! Since Rho = V^* V^T = 1 - U^* U^T it is easy to see that the EigenValues
  ! of Rho that are equal to 1 correspond to the zero eigenvalues of U. It is
  ! exactly the dimension of the null-space of U that determines the number
  ! parity of the HFB state. If this dimension is D then:
  !    \pi_n = (-1)^D
  !
  ! For more see:
  ! G. Bertsch et al., Phys. Rev. A 79, 043602 (2009)
  !-----------------------------------------------------------------------------
    integer             :: i,it,P, C,j
    ! This counts the null-space dimension of U for every signature,parity & 
    ! isospin block.
    integer             :: NullDimension(2,2)
    real(KIND=dp)       :: TotalSignature(2,2), prod
    complex(KIND=dp) :: Temp(HFBSize)
    !---------------------------------------------------------------------------
    HFBColumns  = 0
    do it=1,Iindex
      do P=1,Pindex
        !-----------------------------------------------------------------------
        ! Our first try is always just taking the positive energy quasiparticle
        ! columns. These are the last columns in U and V, since the Hamiltonian
        ! diagonalisation ordered the eigenvalues in ascending order.
        do i=1,blocksizes(P,it)
          HFBColumns(i,P,it) = i + blocksizes(P,it)
        enddo
!         print *, HFBColumns(1:blocksizes(P,it),P,it)
!         do i=1,blocksizes(P,it)
!             print *, QuasiEnergies(HFBColumns(i,P,it),P,it)
!         enddo
      enddo
    enddo
    ! Switch the columns we want to block
    if(allocated(QPExcitations)) call BlockQuasiParticles()

    !---------------------------------------------------------------------------
    ! Extra check if needed
    if(HFBCheck) call CheckUandVColumns(HFBColumns)    
    ! We now explicitly construct Rho first
    RhoHFB = constructRhoHFB(HFBColumns)
    !if(HFBCheck) call CheckRho(RhoHFB)
    !---------------------------------------------------------------------------
    !Note that the entire procedure to combat gapless superconductivity is super-
    !fluous when conserving time-reversal
    if(.not.TRC) then
        !-----------------------------------------------------------------------
        !Now we diagonalise Rho.
        call DiagonaliseRHOHFB
        NullDimension = 0 ; TotalSignature = 0.0_dp
        do it=1,Iindex
          do P=1,Pindex        
            !-------------------------------------------------------------------
            ! Now we look for eigenvalues of Rho that are VERY CLOSE to 1, and 
            ! count them for each isospin-parity block.
            ! Only one `problem' here: when is an eigenvalue close enough to 1? 
            do i=1,blocksizes(P,it)
              if(abs(Occupations(i,P,it) - 1) .lt. 1d-10) then
                NullDimension(P,it) = NullDimension(P,it) + 1
              endif
            enddo

            ! If the number-parity is even, no problem.
            if((-1)**NullDimension(P,it) .eq. HFBNumberParity(P,it)) cycle
            print *,' Gapless Superconductivity'
            !-------------------------------------------------------------------
            ! Now for the fixing part: if a block has wrong number parity, we 
            ! check for the lowest positive quasiparticle energy of that block 
            ! and 'excite' that quasiparticle by exchanging U and V for that qp.          
            ! Of course we need to check that we excite the qp with the correct 
            ! signature.
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            ! Now check for total signature
            do i=1,blocksizes(P,it)
              C = HFBColumns(i,P,it)
              TotalSignature(P,it)=TotalSignature(P,it)+QuasiSignatures(C,P,it)
            enddo
            
            do i=1,blocksizes(P,it)
              C = HFBColumns(i,P,it)
              if(abs(QuasiSignatures(C,P,it) - TotalSignature(P,it)/2 ).lt.1d-8) then      
                  HFBColumns(i,P,it) = 2*blocksizes(P,it) - C + 1
                  
                  print *, 'Excited', P, it, C, HFBColumns(i,P,it)
                  !print *, 'S, TS',  QuasiSignatures(C,P,it), TotalSignature(P,it)
                  !print *, 'NS', QuasiSignatures(HFBColumns(i,P,it),P,it)
                  print *, 'Nulldimension', Nulldimension(P,it), HFBNumberparity(P,it)
                  TotalSignature = 0.0_dp
                  do j=1,blocksizes(P,it)
                    C = HFBColumns(j,P,it)
                    TotalSignature(P,it)=TotalSignature(P,it)+QuasiSignatures(C,P,it)
                  enddo
                  print *, 'NewTS', TotalSignature

                  print *
                  exit
              endif
            enddo
           enddo
        enddo 
      if(HFBCheck) call CheckUandVColumns(HFBColumns)   
      !-------------------------------------------------------------------------
      ! Reconstruct the density
      RhoHFB = constructRhoHFB(HFBColumns)
    endif
    !---------------------------------------------------------------------------
    ! Construct the anomalous density matrix
    KappaHFB = constructKappaHFB(HFBColumns)

  end subroutine ConstructHFBState

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


  function constructRhoHFB (Columns) result(Rho)
  !-----------------------------------------------------------------------------
  ! This function constructs a density matrix from the input of columns from
  ! the U and V matrices.
  !
  ! rho   = V^* V^T
  !-----------------------------------------------------------------------------
    integer, intent(in) :: Columns(HFBSize,Pindex,Iindex)
    integer             :: i,j,k,it,P
    complex(KIND=dp)    :: Rho(HFBSize,HFBSize,Pindex,Iindex)

    Rho = 0.0_dp
    !---------------------------------------------------------------------------
    ! Actual computation
    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          do i=1,blocksizes(P,it)
            do k=1,blocksizes(P,it)
              Rho(i,j,P,it) = Rho(i,j,P,it) +                                  &
              &       Conjg(V(i,Columns(k,P,it),P,it))*V(j,Columns(k,P,it),P,it)      
            enddo
          enddo
        enddo
      enddo
    enddo 

  end function ConstructRHOHFB

  subroutine DiagonaliseRhoHFB
    !---------------------------------------------------------------------------
    ! Diagonalise RhoHFB and find both the occupation numbers and the 
    ! transformation matrix between HFbasis and the canonical basis.
    !
    !---------------------------------------------------------------------------
    complex(KIND=dp), allocatable, save :: Work(:)
    real(KIND=dp), allocatable,    save :: RWork(:)
    integer,                       save :: WorkSize
    complex(KIND=dp), allocatable, save :: Temp(:,:)
    integer                             :: Succes, P, it, N, i, ii, iii

    !---------------------------------------------------------------------------
    ! Some preparations for the diagonalisation routine further on.
    if(.not.allocated(WORK)) then
      ! This RWORK array has to have a specific size, see the ZHEEV docs at
      ! http://www.netlib.org/lapack/explore-html/d6/dee/zheev_8f.html
      allocate(RWORK(3*maxval(blocksizes)-2))
      allocate(Temp(HFBSize,HFBSize))
      Temp = 0.0_dp ;  RWORK = 0.0_dp
      WORKSIZE = -1 ; allocate(WORK(1))
      ! Preliminary call to determine WORKSIZE
      call ZHEEV('V', 'U', maxval(blocksizes), Temp,maxval(blocksizes),      &
      &           Occupations, WORK, WORKSIZE, RWORK, Succes)
      WORKSIZE = ceiling(real(WORK(1)))
      deallocate(WORK); allocate(WORK(WORKSIZE))
    endif
    
    do it=1,Iindex
      do P=1,Pindex          
          Temp          = 0.0_dp
          N             = blocksizes(P,it)
          Temp(1:N,1:N) = RhoHFB(1:N,1:N,P,it)
          
          !-------------------------------------------------------------------
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
          if(SC) then
            do i=1,blocksizes(P,it)       
              ii  = blockindices(i,P,it)
              iii = mod(ii-1,nwt)+1
              if( ii .ne. iii .or. HFBasis(iii)%GetSignature().eq.-1) then
                Temp(i,i) = Temp(i,i) - 2.0_dp
              endif
            enddo
          endif
          
          call ZHEEV('V', 'U', N, Temp(1:N,1:N), N, Occupations(1:N,P,it)    &
          &          , WORK, WORKSIZE, RWORK, Succes)            
          if(Succes.ne.0) then
            call stp('ZHEEV failed to diagonalise RHOHFB', 'Errorcode', Succes)
          endif
          CanTransfo(1:N,1:N,P,it) = Temp(1:N,1:N)
          !------------------------------------------------------------------
          ! Since we shifted the eigenvalues of the negative signature states
          ! by -2, we now need to find the actual occupation numbers.
          ! Notice the slight offset to make sure we are not accidentally
          ! adding two to positive signature eigenvalues
          if(SC) where( Occupations.lt.-0.1_dp) Occupations = Occupations + 2
          ! Numerical safety net
          !where( Occupations.gt.1.0_dp) Occupations=1.0_dp
          !where( Occupations.lt.0.0_dp) Occupations=0.0_dp
      enddo
    enddo
    where (abs(CanTransfo) .lt. 1d-11) CanTransfo = 0.0_dp
  end subroutine DiagonaliseRhoHFB

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

  function ConstructKappaHFB(Columns) result(Kappa)
  !----------------------------------------------------------------------------
  ! This function constructs the anomalous density matrix from the input of 
  ! columns from the U and V matrices.
  !
  ! Kappa   = V^* U^T
  !----------------------------------------------------------------------------
    integer, intent(in) :: Columns(HFBSize,Pindex,Iindex)
    integer             :: i,ii, j, jj, it, P, k, S1, S2, iii, jjj
    complex(KIND=dp)    :: Kappa(HFBSize,HFBSize,Pindex,Iindex)

    Kappa = 0.0_dp
    do it=1,Iindex
      do P=1,Pindex
        do j=1,blocksizes(P,it)
          do i=j+1,blocksizes(P,it)  
            do k=1,blocksizes(P,it)
              Kappa(i,j,P,it) = Kappa(i,j,P,it) +                            &
              &   conjg(V(i,Columns(k,P,it),P,it))*U(j,Columns(k,P,it),P,it)
            enddo
            Kappa(j,i,P,it) = - Kappa(i,j,P,it)
          enddo
        enddo 
      enddo
    enddo

    !---------------------------------------------------------------------------
    ! Check if Kappa couples correctly between signatures
    if(SC) then
      do it=1,Iindex
          do P=1,Pindex
              do i=1,blocksizes(P,it)
                  ii  = blockindices(i,P,it)
                  iii = mod(ii-1,nwt)+1
                  S1 = HFBasis(iii)%GetSignature()
                  if(ii.ne.iii) S1 = -S1
                  do j=1,blocksizes(P,it)
                      jj = blockindices(j,P,it)
                      jjj=mod(jj-1,nwt)+1
                      S2 = HFBasis(jjj)%GetSignature()
                      if(jj .ne. jjj) S2 = -S2
                      if(S1*S2 .ne. -1 .and. abs(Kappa(i,j,P,it)).gt.HFBNumCut) then
                          print *, i,j, S1, S2, Kappa(i,j,P,it)
                          print *, QuasiEnergies(:,P,it)
                          call stp('Kappa Signatures!')
                      endif
                  enddo
              enddo
          enddo
      enddo
    endif

  end function ConstructKappaHFB

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
  real(KIND=dp)                             :: Energy,  RhoII, SR
  integer                                   :: it,P,i,j,S,ii,iii,loc(1),TS,jj,jjj,k
  integer                                   :: Columns(nwt,Pindex,Iindex)
  integer                                   :: P2, C, index

  PairingDisp = 0.0_dp
  !-----------------------------------------------------------------------------
  !Mix the densities, if there is a saved density
  if( .not. all(OldRhoHFB.eq.0.0_dp) ) then
    RhoHFB   = HFBMix * RhoHFB   + (1.0_dp - HFBMix) * OldRhoHFB
  endif
  if( .not. all(KappaHFB.eq.0.0_dp) ) then
    KappaHFB = HFBMix * KappaHFB + (1.0_dp - HFBMix) * OldKappaHFB
  endif    
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
  if(any(Occupations .lt. -HFBNumCut)) then
    ! Notice that we allow some very small negative occupation numbers.
    ! They are entirely due to numerical error, and such errors are present
    ! too in CR8. However, they are masked by setting all negative elements of
    ! RhoHFB density matrix to zero. But in MOCCa, RHOHFB can be complex...
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
        !------------------------------------------------------------------------
        ! Make the transformation
        do j=1,blocksizes(P,it)
          jj = Blockindices(j,P,it)
          jjj= mod(jj-1,nwt)+1
          !Cutoff for numerical stability
          if(abs(CanTransfo(j,C,P,it)) .lt.HFBNumCut) cycle
          !----------------------------------------------------------------------
          !Transformation
          !----------------------------------------------------------------------
          ! Note that the time-reversing here is kinda wasting cpu cycles 
          ! since we apply it on the entire wavefunction. However, I prefer
          ! this routine to act on wavefunctions, not the (more efficient) 
          ! spinors, since in this way MOCCa checks the consistency of 
          ! quantum numbers.
          if(jj.eq.jjj) then
            if(TSC) then
                CanBasis(index) = Canbasis(index) + DBLE(CanTransfo(j,C,P,it))*&
                &                                                    HFBasis(jjj)
            else
                CanBasis(index) = Canbasis(index) +       CanTransfo(j,C,P,it)*&
                &                                                    HFBasis(jjj)
            endif
          else
            ! Note that this should only happen when signature is not conserved,
            ! but TimeReversal is.
            if(TSC) then
                CanBasis(index) = Canbasis(index) + DBLE(CanTransfo(j,C,P,it))*&
                &                                  TimeReverseSpwf(HFBasis(jjj))
            else
                CanBasis(index) = Canbasis(index) +       CanTransfo(j,C,P,it)*&
                &                                  TimeReverseSpwf(HFBasis(jjj))
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
        index = index + 1
      enddo
    enddo
  enddo
  if(index.ne.nwt+1) call stp('Not enough canonical spwfs were constructed.')
  !Do some stuff to get the canbasis in fighting condition
  do i=1,nwt !Here there has to be 'nwt' instead of HFBSize
    call CanBasis(i)%SymmetryOperators()
    call CanBasis(i)%CompDer()
  enddo
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
  ! Save old density and anomalous density matrix.
  OldRhoHFB = RhoHFB ; OldKappaHFB = KappaHFB
  
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
!     do i=1,blocksizes(1,2)
!         do j=1,blocksizes(1,2)
!             write (*,"(f8.3)",advance="no")  DBLE(CanTransfo(i,j,1,2)) 
!         enddo
!         print *
!     enddo
!     print *
!     print '(10f8.3)', QuasiEnergies(1:2*blocksizes(1,2),1,2)
!     print *
!     print *, Columns(1:blocksizes(1,2),1,2)
!     print *
  end function FindCorrectColumns

  subroutine WriteOutKappa(PairingType, OnlyIso)
  !-----------------------------------------------------------------------------
  ! Construct a certain Kappa matrix for writing to an output file. In this way
  ! new calculations need not start from a very bad initial guess.
  !-----------------------------------------------------------------------------
  ! In addition, this routine can supply a nonzero guess for kappa which can 
  ! also be very useful for not getting stuck in a trivial solution with HFB.
  !
  !-----------------------------------------------------------------------------
  ! OnlyIso is an optional integer argument that can make this guess only 
  ! for a specific isospin.
  !-----------------------------------------------------------------------------
    integer, intent(in)       :: PairingType
    integer                   :: i, ii, it, P, sig1,sig2, j,jj,iii,jjj, S1, S2
    integer                   :: startit, stopit
    integer, intent(in), optional :: OnlyIso
    real(KIND=dp)                 :: Occ
  
    select case (PairingType) 
    case(0)
      !------------------------------------------------------------------------
      ! HF Calculation: we can't readily guess zero, since that would
      ! not allow us to start HFB calculations from converged HF calculations.
      ! So we try something more clever. 
      !------------------------------------------------------------------------
      if(present(OnlyIso)) call stp('A guess for KappaHFB for only one isospin' &
      &  //'should not be made in a HF calculation.')
      allocate(KappaHFB(HFBSize,HFBSize,2,2)) ; KappaHFB=0.0_dp
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
      if(present(OnlyIso)) call stp('A guess for KappaHFB for only one isospin' &
      &  //'should not be made in a BCS calculation.')

      !BCS Calculation
      allocate(KappaHFB(HFBSize,HFBSize,2,2)); KappaHFB=0.0_dp
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
        Occ = 0.5_dp
        if(present(OnlyIso)) then
          startit=OnlyIso
          stopit =OnlyIso
        else
          startit=1
          stopit=Iindex
        endif
        do it=startit,stopit
          KappaHFB(:,:,:,it)=0.0_dp
          do P=1,Pindex
            do i=1,blocksizes(P,it)
              ii = Blockindices(i,P,it)
              iii = mod(ii-1,nwt)+1
              S1 = HFBasis(iii)%GetSignature()
              if(iii.ne.ii) S1 = - S1
              do j=1,i-1
                jj = Blockindices(j,P,it)
                jjj = mod(jj-1,nwt)+1
                S2 = HFBasis(jjj)%GetSignature()
                if(jjj.ne.jj) S2 = - S2
                if(S2.eq.S1) cycle
                KappaHFB(i,j,P,it) = Occ
                KappaHFB(j,i,P,it) =-Occ
              enddo
            enddo
          enddo
        enddo
        
    case DEFAULT
        call stp('Unknow PairingType in WriteOutKappa.')
    end select
  end subroutine WriteOutKappa
  
  function LNCr8(Delta, DeltaLN,flag) result(LNLambda)
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
  ! - Don't use this routine when R_z is broken
  !-----------------------------------------------------------------------------
    integer, intent(inout) :: flag(2)
    real(KIND=dp) :: LNLambda(2)
    integer       :: i,j,it, ii, k, P, iii   
    real(KIND=dp) :: c2(2),c3(2),c4(2), trx(2), txd(2), ex(2), gkr(2), erx(2) 
    real(KIND=dp) :: gkx(2),gky(2)
    real(KIND=dp) :: E
    complex(KIND=dp) :: x1, x2, x3,x4,x5
    complex(KIND=dp) :: Chi(HFBSize,HFBSize,2,2), Chika(HFBSize,HFBSize,2,2)
    complex(KIND=dp) :: Gamka(HFBSize,HFBSize,2,2)
    real(KIND=dp)    :: hl1(2) , hl2(2), deno(2), xnum(2)
    complex(KIND=dp), intent(in), allocatable :: Delta(:,:,:,:)
    complex(KIND=dp), intent(in), allocatable :: DeltaLN(:,:,:,:)

    c2 = 0.0_dp ; c3 =0.0_dp ; c4 = 0.0_dp
    trx= 0.0_dp ; txd=0.0_dp ; ex = 0.0_dp ; gkr = 0.0_dp ; erx = 0.0_dp
    gkx = 0.0_dp ; gky = 0.0_dp ; Chi = 0.0_dp ; chika = 0.0_dp ; Gamka = 0.0_dp
       
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
  end function LNCr8

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
    integer             :: io
    integer, allocatable :: Blocked(:)
    integer, allocatable :: Temp(:)
    
    Namelist /Blocking/ Blocked

    io = 0
    allocate(QPExcitations(Block)) ; QPExcitations=0
    
    ! Read the namelist
    allocate(Blocked(block))
    read(unit=*,NML=Blocking, iostat=io)
    if(io.ne.0) call stp('Error on reading blocked particle indices', 'Iostat', io)
    if(any(Blocked.le.0) .or. any(Blocked.gt.nwt)) then
        call stp('Invalid quasiparticle index.')
    endif

    QPexcitations=blocked
  end subroutine ReadBlockingInfo

  subroutine QPindices
    !---------------------------------------------------------------------------
    ! Subroutine that looks for the correct blockindices of the HFBasis 
    ! wavefunctions.
    !
    !---------------------------------------------------------------------------

    integer :: i, N, P,it,j, index, S

    if(.not. allocated(Blockindices)) then
        call stp('QPindices cannot work before the HFBmodule is properly'     &
            &  //' initialized.')
    endif

    N = size(QPexcitations)

    allocate(QPParities(N))    ; QPParities   =0
    allocate(QPIsospins(N))    ; QPIsospins   =0
    allocate(QPSignatures(N))  ; QPSignatures =0
    allocate(QPblockind(N))    ; QPBlockind = 0.0_dp

    do i=1,N
        index = QPExcitations(i)
        
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

    index = 0
    do i=1,N
        P = QPParities(i)
        it= QPIsospins(i)
        do j=1,blocksizes(P,it)
          if(blockindices(j,P,it) .eq. QPExcitations(i) ) index = j
        enddo
        QPBlockind(i) = index
    enddo
  end subroutine QPindices

subroutine PrintBlocking
    !---------------------------------------------------------------------------
    ! Subroutine to print the info on the any blocked states.
    !
    !---------------------------------------------------------------------------

    integer :: N, i

    1 format('Blocking parameters')
    2 format(2x,'  Index     Parity     Isospin      Signature')
    3 format(2x,'---------------------------------------------')
    4 format(i7, 5x, i6, 5x, i8,5x,i8)

    5 format(2x,'---------------------------------------------')
    6 format(2x,'Number parities')
    7 format(2x,'Parity',8x,' -',11x,' +')
    8 format(2x,'Neutron',7x,i2,11x,i2)
    9 format(2x,'Proton ',7x,i2,11x,i2)
    print 1

    N = size(QPExcitations)
    print 2
    print 3
    do i=1,N
      print 4, QPExcitations(i), 2*QPParities(i)-3,2*QPIsospins(i)-3, QPSignatures(i)
    enddo
    print 5
    print *
    print 6
    print 7
    print 5
    print 8, HFBNumberparity(:,1)
    print 9, HFBNumberParity(:,2)
    print 5
  end subroutine PrintBlocking

 subroutine BlockQuasiParticles()
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
    !
    !---------------------------------------------------------------------------
    ! The identification of the quasiparticle state is something else. It either
    ! proceeds directly in the canonical basis or in the HF Basis (not 
    ! implemented yet.)
    !---------------------------------------------------------------------------

    integer          :: N, i, index , j, P, it, C, K, loc(1)
    complex(KIND=dp) :: Temp(HFBSize)
    real(KIND=dp)    :: TempU2

    N = size(QPExcitations)

    if(.not.allocated(QPBlockind)) call QPindices()
  
    do i=1,N
        index = QPblockind(i)
        P     = QPParities(i)
        it    = QPIsospins(i)

        ! Identify the quasiparticle excitation
        loc = 0
        TempU2 = 0.0_dp
        ! Search among the currently occupied columns for which qp 
        ! has the biggest overlap with the asked-for HF state.
        do j=1,blocksizes(P,it)
            if(DBLE(U(index,HFBColumns(j,P,it),P,it)**2) .gt. TempU2) then
                loc = HFBColumns(j,P,it)
                TempU2 = DBLE(U(index,HFBColumns(j,P,it),P,it)**2)
            endif
        enddo
        C   = loc(1)
  
        do j=1,blocksizes(P,it)
            if(HFBColumns(j,P,it) .eq. C) then
                !print *, 'after', 2*blocksizes(P,it) - HFBColumns(j,P,it) +1
                HFBColumns(j,P,it) = 2*blocksizes(P,it) - HFBColumns(j,P,it) +1
                !print *, 'Blocked', HFBColumns(j,P,it) 
                !print *, 'Signature', QuasiSignatures(HFBColumns(j,P,it),P,it)
            endif
        enddo
    enddo
     
  end subroutine BlockQuasiParticles

  subroutine PrintQP()
  !-----------------------------------------------------------------------------
  ! Print the quasiparticles obtained by the HFB proces.
  !-----------------------------------------------------------------------------
    integer             :: i, P, it, j, domU(1), domV(1)
    character(len=7)    :: Species(2)=(/ 'Neutron', 'Proton '/)
    logical             :: skip
    real(KIND=dp)       :: u2mv2, u2pv2, align(6,2*HFBSize,Pindex,IIndex)
    real(KIND=dp)       :: angquantum


    10  format (80 ('_'))
    20  format (80 ('-'))
     1  format ( a7, ' P=', i2 ' quasiparticles')
     2  format ('  n   <Rz>   E_qp       U    V   u^2-v^2 ',3x,      &
        &       '<Jx|T>',4x,'<Jy|T>',6x,'<Jz>', 5x, ' J ')
     3  format ('  n   <Rz>   E_qp       U    V   u^2-v^2 ',2x,      &
        &       '<Jx>',5x,'<Jy|T>',3x,'<Jz>', 5x, ' J ')
     4  format ('  n   <Rz>   E_qp       U    V   u^2-v^2 ',2x,      &
        &       '<Jx>',5x,'<Jy>',5x,'<Jz>', 5x, ' J ')
     !           n   <Rz>   E_qp   
    99  format ( i3, f7.2 , f10.5,2x, i3,2x, i3, 5(3x, f7.2))

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
            &      ,blockindices(domU(1),P,it), blockindices(domV(1),P,it)&
            &      ,u2mv2, align(1:3,i,P,it), angquantum
          enddo
        enddo
    enddo
    print 20
  end subroutine PrintQP

  function QPAlignment() result(align)
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
            JMatrix(i,j,:,:) =   AngularMomentum(HFBasis(j),HFBasis(i), .true.,TRX,TRY,TRZ)
            JMatrix(j,i,:,1) =   JMatrix(i,j,:,1)
            JMatrix(j,i,:,2) = - JMatrix(i,j,:,2)
        enddo
    enddo

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

end module HFB
