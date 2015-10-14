module SpwfFactory
!-------------------------------------------------------------------------------
! Module that can alter the number of wavefunctions present on MOCCa file.
!
!-------------------------------------------------------------------------------

use GenInfo
use Wavefunctions
use SpwfStorage

  type(Spwf),allocatable :: TempStorage(:)
  integer                :: nwtcut=0

contains
  subroutine DecideToCut( final , initial )
    !---------------------------------------------------------------------------
    ! Subroutine that decides whether to throw away some wavefunctions and 
    ! what wavefunctions to keep. Final is the number of wavefunctions desired,
    ! initial is the number of wavefunctions that are currently in memory.
    !---------------------------------------------------------------------------
    integer, intent(in) :: final, initial
    integer             :: newnwtN, newnwtP
    integer, allocatable :: Protonwfs(:), Neutronwfs(:)
  
    !No transforming to be done.
    if(final.ge.initial) return 
    
    !Temporary copout: don't cut. These routines are not ready yet.
    if(final.lt.initial) call stp('Wrong number of wavefunctions.' & 
      & // "Don't try to cut wavefunctions with MOCCa yet.")

    !Sort the HFBasis, to find out what the maximum number of every particle
    !species is
    Protonwfs =OrderSpwfsISO( 1,.false.)
    Neutronwfs=OrderSpwfsISO(-1,.false.)
  
    ! Now look for good proton and neutron wavefunction numbers that sum to 
    ! final, but are sufficient to accomodate all protons and neutrons.
    newnwtN = neutrons  
    if(mod(int(neutrons),2).ne.0 .or. mod(int(protons),2).ne.0) then
      call stp(" DecideToCut can't handle odd number of particles.")
    endif
    if(TRC) then
      newnwtN = newnwtN/2
    endif
    NewnwtP = nwt - newnwtN
    do while(newnwtP.gt.0 .and. newnwtN.lt.size(Neutronwfs)) 
      newnwtN = newnwtN + 1
      newnwtP = nwt - newnwtN
    enddo
    if(newnwtP.gt.size(Protonwfs)) call stp('Too much proton wfs.')
       
    !Effectively cut wavefunctions
    call CutWaveFunctions(newnwtN, newnwtP)
  end subroutine DecideTocut

  subroutine CutWaveFunctions(newnwtN, newnwtP)
    !---------------------------------------------------------------------------
    ! Subroutine that takes the lowest newnwtN and newnwtP wavefunctions for 
    ! each isospin and deletes all superfluous ones. 
    !
    ! SideEffects
    ! - Canonical basis gets deleted, but properly reallocated
    !---------------------------------------------------------------------------
  
    integer, intent(in)  :: newnwtN, newnwtP
    integer, allocatable :: Protonwfs(:), Neutronwfs(:)
    integer              :: i 
    real(KIND=dp), parameter :: OccupationCut=0.01_dp
  
    1 format (' Attention wavefunction nr. ', i4, ' will be deleted.' /,       &
    &         ' However, its occupation number is :' f7.5)
  
  
    if(.not.IC) call stp('Routine CutWaveFunctions is not yet usable when'//   &
    &                    ' isospin is broken.')
    
    !Providing temporary storage
    allocate(TempStorage(newnwtN + newnwtP))
    
    !Sort the HFBasis
    Protonwfs =OrderSpwfsISO( 1,.false.)
    Neutronwfs=OrderSpwfsISO(-1,.false.)
    
    !Copying the wavefunctions to temporary storage
    do i=1,newnwtN
      TempStorage(i) = HFBasis(Neutronwfs(i))
    enddo
    
    do i=1,newnwtP
      TempStorage(newnwtN+i) = HFBasis(Protonwfs(i)) 
    enddo
    ! Print a warning if we delete wavefunctions that are occupied
    ! beyond a certain cutoff point
    do i=1,size(HFBasis)
      if(    any(Neutronwfs(1:newnwtN) .eq.i)                                  &
      & .or. any( Protonwfs(1:newnwtP) .eq.i)) cycle
      
      if(HFBasis(i)%GetOcc().ge.OccupationCut) print 1, i, HFBasis(i)%GetOcc()
    enddo
    
    deallocate(HFBasis)
    allocate(HFBasis(nwt))
    
    if(allocated(CanBasis)) then
      deallocate(CanBasis)
      allocate(CanBasis(nwt))
    endif
    
    !Copying the wavefunctions
    do i=1,nwt
      HFBasis(i) = TempStorage(i)     
    enddo
    
    DensityBasis => HFBasis
    deallocate(TempStorage)

  end subroutine CutWaveFunctions
  

end module
