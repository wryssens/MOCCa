module SpwfStorage
  use CompilationInfo
  use GenInfo
  use WaveFunctions

  implicit none

  save
  !-----------------------------------------------------------------------------
  !Total number of single-particle states, as well as the separate neutron and
  ! and proton count. (These last two are currently only used to print some info
  ! for the analysis scripts.)
  !-----------------------------------------------------------------------------
  integer :: nwt=0, nwn=0, nwp=0
  !-----------------------------------------------------------------------------
  !Fermi Energy and printing window around Fermi Energy
  !-----------------------------------------------------------------------------
  real*8  :: FermiEnergy, PrintingWindow=10.0_dp
  !-----------------------------------------------------------------------------
  ! Whether or not MOCCa is allowed to adjust the number of states taken
  ! into account.
  !-----------------------------------------------------------------------------
  logical :: AdjustNumber=.false.
  !-----------------------------------------------------------------------------
  ! All Single particle wavefunctions in the HF basis.
  !-----------------------------------------------------------------------------
  type (Spwf), public,allocatable,target :: HFBasis(:)
  !-----------------------------------------------------------------------------
  ! All Single particle wavefunctions in the Canonical Basis.
  ! Only allocated when HFB pairing is active.
  type (Spwf), public, allocatable,target :: CanBasis(:)
  integer, allocatable                    :: CanDeriv(:)
  !-----------------------------------------------------------------------------
  ! In addition, a pointer to either the HFBasis or the CanBasis, depending on
  ! the type of parity. It should point to the wavefunctions needed to construct
  ! the mean-field densities in the Densities module.
  !-----------------------------------------------------------------------------
  type (Spwf), public, pointer :: DensityBasis(:)
  !-----------------------------------------------------------------------------
  !Total angular momentum <J_x >, <J_y>, <J_z>, <J^2>, and their values at
  !the previous 7 iterations and the quadratic expectation value
  ! <J_x^2>, <J_y^2>, <J_z^2>.
  !
  ! AMIsoBlock separates the different contributions along as much blocks as
  ! possible.
  !-----------------------------------------------------------------------------
  real(KIND=dp), public :: TotalAngMom(3)= 0.0_dp, AngMomOld(3)         = 0.0_dp
  real(KIND=dp), public :: J2(3)         = 0.0_dp, AMIsoBlock(2,2,2,3)  = 0.0_dp
  real(KIND=dp), public :: JTR(3)        = 0.0_dp, JTI(3)               = 0.0_dp
  !-----------------------------------------------------------------------------
  !Memory for the total dispersion of the occupied(!) spwfs
  !-----------------------------------------------------------------------------
  real(KIND=dp), public :: TotalDispersion
  !-----------------------------------------------------------------------------
  ! Matrixelements of the nabla operator
  !-----------------------------------------------------------------------------
  real(KIND=dp), public,allocatable :: NablaMElements(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! J2Total
  !   Sum of the qudratic expected values <J_i>^2 and the value at the previous
  !   iteration.
  real(KIND=dp), public             ::     J2total(3)    = 0.0_dp
  real(KIND=dp), public             ::  OldJ2total(3)    = 0.0_dp
  !-----------------------------------------------------------------------------
  ! Add this number of spwfs into every parity-isospin block
  integer :: PlusSpwf=0
  !-----------------------------------------------------------------------------
  ! Ali-rotation to perform on every parity-isospin block on input of the spwfs.
  real(KIND=dp) :: aliy(2,2) = 0.0_dp

contains

  subroutine ReadSpwfStorageInfo()
    !---------------------------------------------------------------------------
    ! Subroutine that reads the necessary info on this module.
    !---------------------------------------------------------------------------
    NameList /SpwfStorage/ nwt, AdjustNumber, PrintingWindow, PlusSpwf, aliy
    read(unit=*, NML=SpwfStorage)

    if(nwt.le.0) call stp('Negative number of states.', 'Nwt', nwt)

    ! The admissible number of states to take into account is less strict
    ! when time Reversal is conserved
    if(TRC) then
      if((nwt).lt.(Neutrons+Protons)/2) then
        call stp('Not enough states.', 'nwt',nwt, 'Particles', Neutrons+Protons)
      endif
    else
      if((nwt).lt.(Neutrons+Protons)) then
        call stp('Not enough states.', 'nwt',nwt, 'Particles', Neutrons+Protons)
      endif
    endif
  end subroutine ReadSpwfStorageInfo

  subroutine ChangeNumberWaveFunctions( NewNumber )
    !---------------------------------------------------------------------------
    ! A subroutine that changes the total number of wavefunctions.
    !---------------------------------------------------------------------------
    ! ATTENTION!!
    ! Note that in the case of the number of wavefunctions diminishing,
    ! only the lowest occupied states for each isospin are copied!
    ! This means that this subroutine is not finished, I have not yet
    ! implemented the copying of states beyond the occupied ones when
    ! diminishing the space.
    !---------------------------------------------------------------------------
    integer, intent(in)     :: NewNumber
    integer                 :: i, N, M, NeutWf, ProtWf
    type(Spwf), allocatable :: TempStorage(:), TempCanStorage(:)
    integer, allocatable    :: OrderN(:), OrderP(:)

    NeutWf=1
    ProtWf=1
    if(allocated(HFBasis)) then
      !-------------------------------------------------------------------------
      !Ordering if neccessary the neutron & Proton wavefunction
      if(nwt .gt. NewNumber) then
        OrderP = OrderSpwfsISO( 1, .false.)
        OrderN = OrderSpwfsISO(-1, .false.)
      endif
      !-------------------------------------------------------------------------
      !Copying the previous functions
      N = size(HFBasis)
      allocate(TempStorage(N))
      do i=1,N
        TempStorage(i)= HFBasis(i)
      enddo
      if(allocated(CanBasis)) then
        allocate(TempCanStorage(N))
        do i=1,N
          TempCanStorage(i) = CanBasis(i)
        enddo
      endif
      !-------------------------------------------------------------------------
      !Changing the size of the basis arrays.
      deallocate(HFBasis)
      allocate(HFBasis(NewNumber))
      if(allocated(CanBasis)) then
        deallocate(CanBasis)          ; deallocate(CanDeriv)
        allocate(CanBasis(NewNumber)) ; allocate(CanDeriv(Newnumber))
      endif

      !-------------------------------------------------------------------------
      ! Recopy
      M = min(N,NewNumber)
      if (M .ge. N) then
        !No problem if the space gets enlarged
        do i=1,M
          HFBasis(i) = TempStorage(i)
        enddo
        if(allocated(CanBasis)) then
          do i=1,N
            CanBasis(i) = TempCanStorage(i)
          enddo
          DensityBasis => CanBasis
        else
          DensityBasis => HFBasis
        endif
      else
        call stp('No reliable way yet to make the size of HFBasis smaller!', 'N', N)
      endif
      nwt = NewNumber
    else
      !-------------------------------------------------------------------------
      !Else, just make space.
      allocate(HFBasis(NewNumber))
      !NOte that the Canbasis is not allocated here, but in the pairing module.
    endif
    return
  end subroutine ChangeNumberWaveFunctions

!  subroutine GramSchmidt_NEW
!     use GenInfo
!     use CompilationInfo

!     implicit none

!     integer      :: nw,mw, Signature, Parity, Isospin, i
!     real(KIND=dp):: Norm
!     type(Spinor) :: ValueOne,ValueTwo,Temp, Temp2
!     real(KIND=dp):: MatrixElement(2)

!     if(.not.allocated(Temp%Grid)) allocate(Temp%Grid(nx,ny,nz,4,1))

!      do nw=1,nwt
!        ! First normalise \Psi_{nw}
!        Norm=0.0_dp
!        do i=1,4*nx*ny*nz
!          Norm = Norm + HFBasis(nw)%Value%Grid(i,1,1,1,1)**2
!        enddo
!        Norm = Norm * dv
!        Norm = (1.0d0/sqrt(Norm))  
!        do i=1,4*nx*ny*nz
!          HFBasis(nw)%Value%Grid(i,1,1,1,1) = Norm* HFBasis(nw)%Value%Grid(i,1,1,1,1)
!        enddo
!      !call HFBasis(nw)%CompNorm()

!      Isospin   = HFBasis(nw)%GetIsospin()
!      Parity    = HFBasis(nw)%GetParity()
!      Signature = HFBasis(nw)%GetSignature()
!      !-------------------------------------------------------------------------
!      ! Then subtract the projection on \Psi_{nw} from all the following Spwf.
!      ! Re(\Psi(\sigma)_{mw}) =
!      ! Re(\Psi(\sigma)_{nw})-Re(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
!      !                      + Im(< \Psi_{nw}|\Psi_{mw} >)Im(\Psi(\sigma)_{nw})
!      ! Im(\Psi(\sigma)_{mw}) =
!      ! Im(\Psi(\sigma)_{mw})-Re(< \Psi_{nw}|\Psi_{mw} >) Im(\Psi(\sigma)_{nw})
!      !                      -Im(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
!      !-------------------------------------------------------------------------
!      ! Note that the imaginary part of the inproduct only needs to be taken
!      ! into account when Time Simplex is not conserved.
!      !-------------------------------------------------------------------------
!      do mw=nw+1,nwt
!        !-----------------------------------------------------------------------
!        ! Do not waste time in manipulating the grids if the matrixelement
!        ! is zero. This is fulfilled, for example, when the two SPWF do not
!        ! share all quantum numbers.
!        !-----------------------------------------------------------------------
!        if(HFBasis(mw)%GetIsospin()  .ne.Isospin)   cycle
!        if(HFBasis(mw)%GetParity()   .ne.Parity)    cycle
!        if(HFBasis(mw)%GetSignature().ne.Signature) cycle
!        
!        MatrixElement = 0.0_dp
!        do i=1,4*nx*ny*nz
!          MatrixElement(1)= MatrixElement(1) + &
!          & HFBasis(mw)%Value%Grid(i,1,1,1,1)*HFBasis(nw)%Value%Grid(i,1,1,1,1) 
!        enddo
!        
!        do i=1,4*nx*ny*nz
!          HFBasis(mw)%Value%Grid(i,1,1,1,1) = HFBasis(mw)%Value%Grid(i,1,1,1,1) - &
!          &                 MatrixElement(1) * HFBasis(nw)%Value%Grid(i,1,1,1,1)
!        enddo

!!        if(.not.TSC) then
!!         !Imaginary Part of MatrixElement (Zero when timesimplex is conserved).
!!         Temp2 = MultiplyI(HFBasis(nw)%Value)
!!         do i=1,4*nx*ny*nz
!!            HFBasis(mw)%Value%Grid(i,1,1,1,1) = HFBasis(mw)%Value%Grid(i,1,1,1,1) + &
!!            &                           MatrixElement(2) * Temp2%Grid(i,1,1,1,1)
!!         enddo
!!        endif

!!        ! If signature is broken, but time reversal is conserved, also
!!        ! orthogonalise against the time-reversed functions
!!        if(.not.SC .and. TRC) then
!!          Temp2 = TimeReverse(HFBasis(nw)%Value)
!!          MatrixElement(1)  = InproductSpinorReal(Temp2, HFBasis(mw)%Value)
!!          if(.not.TSC) then
!!             MatrixElement(2)  = InproductSpinorImaginary(Temp2, HFBasis(mw)%Value)
!!             Temp2 = MultiplyI(Temp2)
!!             Temp = Temp - MatrixElement(2)*Temp2
!!             Temp2 = -MultiplyI(Temp2)
!!          endif
!!          do i=1,4*nx*ny*nz
!!            HFBasis(mw)%Value%Grid(i,1,1,1,1) = HFBasis(mw)%Value%Grid(i,1,1,1,1) - &
!!            &                                   MatrixElement(1) * Temp2%Grid(i,1,1,1,1)
!!          enddo
!!          !Temp = Temp - MatrixElement(1)*Temp2
!!        endif
!        !Save the result to the corresponding wavefunction
!        !call HFBasis(mw)%SetGrid(Temp)
!      enddo
!    enddo

!  end subroutine GramSchmidt_NEW

  subroutine GramSchmidt_OLD
    !---------------------------------------------------------------------------
    ! This subroutine uses a Gram-Schmidt scheme to orthonormalise the Spwfs in
    ! the HF basis.
    !---------------------------------------------------------------------------
    use GenInfo
    use CompilationInfo

    implicit none

    integer      :: nw,mw, Signature, Parity, Isospin, i
    real(KIND=dp):: Norm
    type(Spinor) :: ValueOne,ValueTwo,Temp, Temp2
    real(KIND=dp):: MatrixElement(2)


    if(.not.allocated(Temp%Grid)) allocate(Temp%Grid(nx,ny,nz,4,1))

    do nw=1,nwt
      ! First normalise \Psi_{nw}
      call HFBasis(nw)%CompNorm()
      Norm=HFBasis(nw)%GetNorm()

      do i=1,4*nx*ny*nz
        HFBasis(nw)%Value%Grid(i,1,1,1,1) = (1.0d0/sqrt(Norm)) * HFBasis(nw)%Value%Grid(i,1,1,1,1)
      enddo
      !call HFBasis(nw)%CompNorm()

      Isospin   = HFBasis(nw)%GetIsospin()
      Parity    = HFBasis(nw)%GetParity()
      Signature = HFBasis(nw)%GetSignature()
      !-------------------------------------------------------------------------
      ! Then subtract the projection on \Psi_{nw} from all the following Spwf.
      ! Re(\Psi(\sigma)_{mw}) =
      ! Re(\Psi(\sigma)_{nw})-Re(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
      !                      + Im(< \Psi_{nw}|\Psi_{mw} >)Im(\Psi(\sigma)_{nw})
      ! Im(\Psi(\sigma)_{mw}) =
      ! Im(\Psi(\sigma)_{mw})-Re(< \Psi_{nw}|\Psi_{mw} >) Im(\Psi(\sigma)_{nw})
      !                      -Im(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
      !-------------------------------------------------------------------------
      ! Note that the imaginary part of the inproduct only needs to be taken
      ! into account when Time Simplex is not conserved.
      !-------------------------------------------------------------------------
      do mw=nw+1,nwt
        !-----------------------------------------------------------------------
        ! Do not waste time in manipulating the grids if the matrixelement
        ! is zero. This is fulfilled, for example, when the two SPWF do not
        ! share all quantum numbers.
        !-----------------------------------------------------------------------
        if(HFBasis(mw)%GetIsospin().ne.Isospin)     cycle
        if(HFBasis(mw)%GetParity().ne.Parity)       cycle
        if(HFBasis(mw)%GetSignature().ne.Signature) cycle

        MatrixElement=InProduct(HFBasis(mw),HFBasis(nw))
        !ValueTwo = HFBasis(mw)%GetValue()
        !________________________________________________________
        ! Working version with overloaded operators.
        ! However, this is slow due to repeated loading/storing.
        !Real Part of MatrixElement
        !Temp = ValueOne
        !Temp = (-MatrixElement(1))*Temp
        !Temp = Temp + ValueTwo
        !________________________________________________________

        do i=1,4*nx*ny*nz
          HFBasis(mw)%Value%Grid(i,1,1,1,1) = HFBasis(mw)%Value%Grid(i,1,1,1,1) - &
          &                 MatrixElement(1) * HFBasis(nw)%Value%Grid(i,1,1,1,1)
        enddo

        if(.not.TSC) then
         !Imaginary Part of MatrixElement (Zero when timesimplex is conserved).
         Temp2 = MultiplyI(HFBasis(nw)%Value)
         do i=1,4*nx*ny*nz
            HFBasis(mw)%Value%Grid(i,1,1,1,1) = HFBasis(mw)%Value%Grid(i,1,1,1,1) + &
            &                           MatrixElement(2) * Temp2%Grid(i,1,1,1,1)
         enddo
        endif

        ! If signature is broken, but time reversal is conserved, also
        ! orthogonalise against the time-reversed functions
        if(.not.SC .and. TRC) then
          Temp2 = TimeReverse(HFBasis(nw)%Value)
          MatrixElement(1)  = InproductSpinorReal(Temp2, HFBasis(mw)%Value)
          if(.not.TSC) then
             MatrixElement(2)  = InproductSpinorImaginary(Temp2, HFBasis(mw)%Value)
             Temp2 = MultiplyI(Temp2)
             Temp = Temp - MatrixElement(2)*Temp2
             Temp2 = -MultiplyI(Temp2)
          endif
          do i=1,4*nx*ny*nz
            HFBasis(mw)%Value%Grid(i,1,1,1,1) = HFBasis(mw)%Value%Grid(i,1,1,1,1) - &
            &                                   MatrixElement(1) * Temp2%Grid(i,1,1,1,1)
          enddo
          !Temp = Temp - MatrixElement(1)*Temp2
        endif
        !Save the result to the corresponding wavefunction
        !call HFBasis(mw)%SetGrid(Temp)
      enddo
    enddo

  end subroutine GramSchmidt_OLD

  subroutine GramSchmidt!_Energy
    !---------------------------------------------------------------------------
    ! This subroutine uses a Gram-Schmidt scheme to orthonormalise the Spwfs in
    ! the HF basis. Small point of interest: the orthonormalisation is done in 
    ! order of ascending energy within each isospin block, in order to avoid
    ! wasting CPU cycles reordering levels.
    !---------------------------------------------------------------------------
    use GenInfo
    use CompilationInfo

    implicit none

    integer              :: nw,mw, Signature, Parity, Isospin, i, it,k,l
    integer, allocatable :: indices(:)
    real(KIND=dp):: Norm
    type(Spinor) :: ValueOne,ValueTwo,Temp, Temp2
    real(KIND=dp):: MatrixElement(2)


    if(.not.allocated(Temp%Grid)) allocate(Temp%Grid(nx,ny,nz,4,1))

    do it=1,2
        ! Order spwfs of this isospin by their energy
        indices = OrderSpwfsISO(it*2 - 3, .false.) 
        
        do k=1,size(indices)
            nw = indices(k)
            
            ! First normalise \Psi_{nw}
            call HFBasis(nw)%CompNorm()
            Norm=HFBasis(nw)%GetNorm()
            do i=1,4*nx*ny*nz
                HFBasis(nw)%Value%Grid(i,1,1,1,1) = (1.0d0/sqrt(Norm)) * HFBasis(nw)%Value%Grid(i,1,1,1,1)
            enddo
            
            Isospin   = HFBasis(nw)%GetIsospin()
            Parity    = HFBasis(nw)%GetParity()
            Signature = HFBasis(nw)%GetSignature()
            
            !-------------------------------------------------------------------------
            ! Then subtract the projection on \Psi_{nw} from all the following Spwf.
            ! Re(\Psi(\sigma)_{mw}) =
            ! Re(\Psi(\sigma)_{nw})-Re(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
            !                      + Im(< \Psi_{nw}|\Psi_{mw} >)Im(\Psi(\sigma)_{nw})
            ! Im(\Psi(\sigma)_{mw}) =
            ! Im(\Psi(\sigma)_{mw})-Re(< \Psi_{nw}|\Psi_{mw} >) Im(\Psi(\sigma)_{nw})
            !                      -Im(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
            !-------------------------------------------------------------------------
            ! Note that the imaginary part of the inproduct only needs to be taken
            ! into account when Time Simplex is not conserved.
            !-------------------------------------------------------------------------
            do l=k+1,size(indices)
                mw = indices(l)
                
                !-----------------------------------------------------------------------
                ! Do not waste time in manipulating the grids if the matrixelement
                ! is zero. This is fulfilled, for example, when the two SPWF do not
                ! share all quantum numbers.
                !-----------------------------------------------------------------------
                if(HFBasis(mw)%GetIsospin().ne.Isospin)     cycle
                if(HFBasis(mw)%GetParity().ne.Parity)       cycle
                if(HFBasis(mw)%GetSignature().ne.Signature) cycle

                MatrixElement=InProduct(HFBasis(mw),HFBasis(nw))

                do i=1,4*nx*ny*nz
                    HFBasis(mw)%Value%Grid(i,1,1,1,1) = HFBasis(mw)%Value%Grid(i,1,1,1,1) - &
                    &                 MatrixElement(1) * HFBasis(nw)%Value%Grid(i,1,1,1,1)
                enddo

                if(.not.TSC) then
                    !Imaginary Part of MatrixElement (Zero when timesimplex is conserved).
                    Temp2 = MultiplyI(HFBasis(nw)%Value)
                    do i=1,4*nx*ny*nz
                        HFBasis(mw)%Value%Grid(i,1,1,1,1) = HFBasis(mw)%Value%Grid(i,1,1,1,1) + &
                        &                           MatrixElement(2) * Temp2%Grid(i,1,1,1,1)
                    enddo
                endif

                ! If signature is broken, but time reversal is conserved, also
                ! orthogonalise against the time-reversed functions
                if(.not.SC .and. TRC) then
                    Temp2 = TimeReverse(HFBasis(nw)%Value)
                    MatrixElement(1)  = InproductSpinorReal(Temp2, HFBasis(mw)%Value)

                    if(.not.TSC) then
                        MatrixElement(2)  = InproductSpinorImaginary(Temp2, HFBasis(mw)%Value)
                    endif   
                  
                    do i=1,4*nx*ny*nz
                        HFBasis(mw)%Value%Grid(i,1,1,1,1) =                    &
                        &                  HFBasis(mw)%Value%Grid(i,1,1,1,1) - &
                        &               MatrixElement(1) * Temp2%Grid(i,1,1,1,1)
                    enddo
                    if(.not. TSC)then
                        ! Also orthonormalize against the time-reversed partner, 
                        ! for the imaginary part of 
                        temp2 = multiplyI(temp2)
                        do i=1,4*nx*ny*nz
                            HFBasis(mw)%Value%Grid(i,1,1,1,1) =                &
                            &              HFBasis(mw)%Value%Grid(i,1,1,1,1) - &
                            &           MatrixElement(2) * Temp2%Grid(i,1,1,1,1)
                        enddo
                    endif
                endif
            enddo
        enddo
    enddo
  end subroutine GramSchmidt!_Energy
  subroutine DeriveAll()
    !---------------------------------------------------------------------------
    ! This subroutine derives all of the Spwf, using the subroutine CompDer.
    !---------------------------------------------------------------------------
    use Derivatives
    integer            :: i
    real(KIND=dp)      :: time(2), temp(nx), temp2(nx)

    do i=1,nwt
        call HFBasis(i)%CompDer()
    enddo

    if(allocated(CanBasis)) then
      do i=1,nwt
            call CanBasis(i)%CompDer()
      enddo
    endif
    return
  end subroutine DeriveAll
  
  subroutine N2LODerive()
    !---------------------------------------------------------------------------
    ! This subroutine derives all of the Spwf, using the subroutine CompDer.
    !---------------------------------------------------------------------------
    integer            :: i
    real(KIND=dp)      :: time(2)

    do i=1,nwt
        call HFBasis(i)%CompDer()
        call HFBasis(i)%CompSecondDer()
    enddo
    if(allocated(CanBasis)) then
      do i=1,nwt
        call CanBasis(i)%CompDer()
        call CanBasis(i)%CompSecondDer()
      enddo
    endif
    return
  end subroutine N2LODerive

  subroutine PrintSpwf(PairingType, Fermi, PrintAllSpwf)
    !---------------------------------------------------------------------------
    ! This subroutine prints out all info for all relevant wavefunctions, both
    ! canonical and hartree-fock basis.
    !---------------------------------------------------------------------------
    ! Note that this routine needs pairingtype as input, as it is not accessible
    ! to this module explicitly.
    !---------------------------------------------------------------------------
    integer, allocatable          :: Order(:)
    integer                       :: i, it, PrintType, ifFilled
    integer, intent(in)           :: PairingType
    character(len=160)            :: HFheader='', CanHeader=''
    real(KIND=dp), intent(in)     :: Fermi(2)
    real(KIND=dp)                 :: Distance    
    logical, intent(in)           :: PrintAllSpwf

    10 format (21 ('-'), ' Sp wavefunctions ', 41('-'))
    20 format (97 ('-'))
    30 format (97 ('_'),/,3x , 'Neutron wavefunctions')
    40 format (97 ('_'),/,3x , 'Proton  wavefunctions')
    50 format (97 ('_'),/,3x , 'HF Basis')
    60 format (97 ('_'),/,3x , 'Canonical Basis' )
    !---------------------------------------------------------------------------
    ! Different explanatory headers, who match with the formats in subroutine
    ! PrintHF and PrintCanonical.
    !---------------------------------------------------------------------------
    ! Headers for Hartree-fock basis
    !---------------------------------------------------------------------------
    ! Hartree-Fock calculations
    1  format (5x,'n',3x,'<P>',3x,'<Rz>',4x,'v^2',5x,'E_sp', 4x,'Var(h)',4x,   &
    &         'r^2',4x, 'Jx',4x,'JxT',3x,'Jy',4x,'JyT',3x,'Jz',4x,'JzT',3x,'J' &
    &          ,3x, 'Sx')
    !---------------------------------------------------------------------------
    ! BCS calculations
    2  format (6x,'n',4x,'<P>', 4x,'<Rz>',4x,'v^2',4x,'Delta',4x,'E_sp',4x,   &
    &          'Var(h)',3x,'r^2',1x, 'Jx',1x,'JxT',1x, 'Jy',1x,'JyT',1x,'Jz', &
    &          1x,'JzT',1x,'J', 3x, 'Sx' )
    !---------------------------------------------------------------------------
    ! HFB calculations
    ! MB 18/12/16 adding the last label 'SzT' to this produces a crash at the CC
    ! with the message:
    !
    ! forrtl: severe (66): output statement overflows record, unit -5, file
    ! Internal Formatted Write
    !
    3  format (5x,'n',3x,'<P>',3x,'<Rz>',3x,'Rhoii',2x,' Delta ', 1x, ' m ',3x,&
    &          'E_sp', 4x,'Var(h)',4x,                                         &
    &         'r^2',4x, 'Jx',5x,'JxT',4x,'Jy',5x,'JyT',4x,'Jz',5x,'JzT',5x,'J',&
    &         3x,'Simplex',2x,'Sx',5x,'SxT',4x,'Sy',5x,'SyT',                  &
    &         4X,'Sz' ,5x )
 
    ! Headers for Canonical basis
    !---------------------------------------------------------------------------
    4  format (5x,'n',3x,'<P>',3x,'<Rz>',4x,'v^2', 5x,'E_sp', 5x,  &
    &         'r^2',4x, 'Jx',4x,'JxT',5x,'Jy',4x,'JyT',5x,'Jz',4x,'JzT',4x,'J',&
    &         6x,'Sx')
    100 format (97('_'))

    print 10

    select case(PairingType)

    case(0)
       write(HFheader, fmt=1)
       printType = 1
    case(1)
       write(HFheader, fmt=2)
       PrintType = 2
    case(2)
       write(HFheader, fmt=3)
       write(CanHeader,fmt=4)
       PrintType = 3
    end select

    if(IC) then
        do it=1,2
            if(it.eq.1) print 30
            if(it.eq.2) print 40
            if(PairingType.eq.2) print 50
            !-------------------------------------------------------------------
            !1) Printing the HF basis
            !   a) Get the ordering of the Spwfs for every isospin
            Order = OrderSpwfsIso( (it*2) - 3, .false. )
            print *, HFheader
            print 100
            do i=1,size(Order)
                Distance = abs(HFBasis(Order(i))%GetEnergy() - Fermi(it))
                ! Print in any case if all SPWFS are asked for
                if(PrintAllSpwf) Distance = 0
                if(Distance.lt.PrintingWindow) then
                    ! We also print the sum of particles that would be
                    ! there if the orbitals were completely filled.
                    ifFilled = i
                    if(TRC) ifFilled=2*i
                    write(*, fmt='(i4,1x)', ADVANCE='no') ifFilled
                    call HFBasis(Order(i))%PrintHF(Order(i),PrintType)
                endif
            enddo
            !-------------------------------------------------------------------
            !2) Printing the Canonical basis
            !   a) Get the ordering of the Spwfs for every isospin
            if(PairingType.ne.2) cycle
            print 60
            Order = OrderSpwfsIso( (it*2) - 3, .true. )
            print *, CanHeader
            print 100
            do i=1,size(Order)
                Distance = abs(CanBasis(Order(i))%GetEnergy() - Fermi(it))
                ! Print in any case if all SPWFS are asked for
                if(PrintAllSpwf) Distance = 0
                if(Distance.lt.PrintingWindow) then
                    ! We also print the sum of particles that would be
                    ! there if the orbitals were completely filled.
                    ifFilled = i
                    if(TRC) ifFilled=2*i
                    write(*, fmt='(i4,1x)', ADVANCE='no') ifFilled
                    call CanBasis(Order(i))%PrintCanonical(Order(i), PrintType)
                endif
            enddo
        enddo
    else
        call stp('No Isospin breaking printing yet.')
    endif

    print 20

  end subroutine PrintSpwf

  subroutine UpdateAM(SaveOlderValues)
    !---------------------------------------------------------------------------
    ! This subroutine calculates the angular momentum in each direction for
    ! every Spwf, and sums it.
    !---------------------------------------------------------------------------
    integer :: wave, i,P,S,it
    logical, intent(in) :: SaveOlderValues
    !Save the values of the previous iteration

    if(SaveOlderValues) then
        AngMomOld   = TotalAngMom
        OldJ2Total  = J2Total
    endif

    TotalAngMom = 0.0_dp ; J2Total = 0.0_dp ; AMIsoblock = 0.0_dp
    JTR = 0.0_dp  ; JTI = 0.0_dp     !;  AMTIsoblock = 0.0_dp 
    
    do wave=1,nwt
      call HFBasis(wave)%DiagSpin()
      call HFBasis(wave)%DiagAng()
      ! When canonical basis is allocated we need to recalculate the angular
      ! momentum for that basis too.
      if(allocated(CanBasis)) then
        call  CanBasis(wave)%DiagAng()
      endif
    enddo

    !---------------------------------------------------------------------------
    !Only compute angular momentum when TimeReversal is not conserved
    if(.not. TRC) then
      do wave=1,nwt
        do i=1,3
          TotalAngMom(i) = TotalAngMom(i) + &
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%J(i)
          JTR(i)          = JTR(i)    +     & 
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%JTR(i)
          JTI(i)          = JTI(i)    +     & 
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%JTI(i)
          J2Total(i)     = J2Total(i) +     &
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%J2(i)

          S = (DensityBasis(wave)%GetSignature() + 3)/2
          P = (DensityBasis(wave)%GetParity()    + 3)/2
          it= (DensityBasis(wave)%GetIsospin()   + 3)/2

          AMIsoblock(S,P,it,i) = AMIsoblock(S,P,it,i) +                  &
          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%J(i)
!          AMTIsoblock(S,P,it,i)= AMTIsoblock(S,P,it,i) +                 &
!          & DensityBasis(wave)%GetOcc()*DensityBasis(wave)%JTR(i)
        enddo
      enddo
      !-------------------------------------------------------------------------
      ! Artificially set total J_x and J_y to zero when dictated by symmetries.
      ! In that case the above summation is not \sum v^2_i <Psi_i|J_x|\Psi_i>
      ! but rather \sum v^2_i <Psi_i|J_x T |\Psi_i>
      if(SC)          TotalAngMom(1) = 0.0
      if(SC .or. TSC) TotalAngMom(2) = 0.0
    endif
    !---------------------------------------------------------------------------
  end subroutine UpdateAM

  function OrderSpwfs(Canonical) result(Indices)
    !---------------------------------------------------------------------------
    ! Function that returns the order of the wavefunctions, that is: a list of
    ! integers denoting the indices of the wavefunctions ordered according to
    ! increasing energy.
    !---------------------------------------------------------------------------
    integer       :: i, HolePos, Indices(nwt), ToInsertIndex
    real(Kind=dp) :: Energies(nwt), ToInsert
    logical, intent(in) :: Canonical
    type(Spwf), pointer :: ToSort(:)

    if(Canonical) then
      ToSort => CanBasis
    else
      ToSort => HFBasis
    endif

    !Copying the energies
    do i=1,nwt
      Energies(i) = ToSort(i)%GetEnergy()
      Indices(i)  = i
    enddo
    !Sort the energies
    do i=2,nwt
      !Make a hole at index i
      ToInsert = Energies(i)
      HolePos  = i
      ToInsertIndex = Indices(i)
      do while(ToInsert.lt.Energies(HolePos-1))
        !Move the hole one place down
        Energies(HolePos) = Energies(HolePos-1)
        Indices(HolePos) = Indices(HolePos-1)
        HolePos = HolePos - 1
        if(HolePos.eq.1.0_dp) exit
      enddo
      !Insert the energy at the correct place
      Energies(HolePos) = ToInsert
      Indices(HolePos) = ToInsertIndex
    enddo
  end function OrderSpwfs

  function OrderSpwfsISO(Isospin, Can) result(Indices)
    !---------------------------------------------------------------------------
    ! Orders the wavefunctions, but keeps the neutrons and protons separate.
    ! Obviously needs to go when Isospin symmetry is broken.
    !---------------------------------------------------------------------------
    integer, allocatable :: Indices(:)
    integer              :: i, nwf, k, HolePos, ToInsertIndex
    integer, intent(in)  :: Isospin
    real(Kind=dp), allocatable :: Energies(:)
    real(Kind=dp)        :: ToInsert
    type(Spwf), pointer  :: ToSort(:)
    logical, intent(in)  :: Can

    if(Can) then
      ToSort => CanBasis
    else
      ToSort => HFBasis
    endif

    !Count the number of proton and neutron wavefunctions
    nwf = 0
    do i=1,nwt
      if(ToSort(i)%GetIsospin() .eq. Isospin) nwf = nwf +1
    enddo
    allocate(Indices(nwf), Energies(nwf))
    !Filling Energies & Indices
    k=1
    do i=1,nwt
      if(ToSort(i)%GetIsospin() .eq. Isospin) then
        Energies(k)= ToSort(i)%GetEnergy()
        Indices(k) = i
        k= k +1
      endif
    enddo

    !Sort the energies
    do i=2,nwf
      !Make a hole at index i
      ToInsert = Energies(i)
      HolePos  = i
      ToInsertIndex = Indices(i)
      do while(ToInsert.lt.Energies(HolePos-1))
        !Move the hole one place down
        Energies(HolePos) = Energies(HolePos-1)
        Indices(HolePos) = Indices(HolePos-1)
        HolePos = HolePos - 1
        if(HolePos.eq.1.0_dp) exit
      enddo
      !Insert the energy at the correct place
      Energies(HolePos) = ToInsert
      Indices(HolePos) = ToInsertIndex
    enddo
  end function OrderSpwfsISO

  subroutine CompNablaMElements()
    !---------------------------------------------------------------------------
    ! Computes the matrix elements of Nabla
    !
    !   < Psi_i | \nabla | \Psi_j >
    !
    ! Some notes:
    ! 1) They are, in general, complex!
    ! 2) Symmetries restrict them in weird ways:
    !    Same p             => <Px> = <Py> = <Pz> = 0
    !    Same signature     => <Px> = <Py> = 0
    !    Different signature=> <Pz> = 0
    !    Different isospin  => <Px> = <Py> = <Pz> = 0
    ! 3) We must not forget the time-reversed states when time reversal is
    !    conserved!
    !---------------------------------------------------------------------------
    ! Why do we need them?
    !  1) Calculation of 2-body COM correction in Energy module.
    !---------------------------------------------------------------------------
    ! This is one of the most complicated routines in MOCCa, due to the
    ! complications of the symmetries here. Any suggestions to make this
    ! less complicated are absolutely welcome!
    ! P.S. I nominate the subroutine pipj in CR8/EV8/EV4 as the place in the
    ! respective codes that needs clarification the most, as several
    ! transformation are totally implicit.
    !---------------------------------------------------------------------------
     use Spinors, only : InproductSpinor

     integer     :: i,j,m, mstart, mstop
     type(Spinor):: Value, Der

     !Allocating the correct number of NablaMElements
     if(.not.allocated(NablaMElements)) then
        allocate(NablaMElements(  nwt,  nwt,3,2))
     endif

     NablaMElements= 0.0_dp

     if(.not.TRC) then
        do i=1,nwt
          Value = DensityBasis(i)%GetValue()
          do j=1,nwt
            !-------------------------------------------------------------------
            ! Cycling if not the same isospin or same parity (when P is conserved)
            !-------------------------------------------------------------------
            if(DensityBasis(i)%GetIsospin().ne.DensityBasis(j)%GetIsospin()) cycle
            if(PC .and. (DensityBasis(i)%GetParity() .eq.DensityBasis(j)%GetParity()))  cycle
            !-------------------------------------------------------------------
            ! If signature is conserved we calculate either (Nx,Ny) or (Nz)
            if(SC) then            
                if(DensityBasis(i)%Signature .ne. DensityBasis(j)%Signature) then
                  Der = DensityBasis(j)%GetDer(1)
                  NablaMElements(i,j,1,1) =   InproductSpinorReal(Value,Der)
                  if(.not. TSC) NablaMElements(i,j,1,2) = InproductSpinorImaginary(Value,Der)
                  
                  Der = DensityBasis(j)%GetDer(2)
                  NablaMElements(i,j,2,2) =   InproductSpinorImaginary(Value,Der)
                  if(.not. TSC) NablaMElements(i,j,2,2) = InproductSpinorReal(Value,Der)
                else
                  Der = DensityBasis(j)%GetDer(3)
                  NablaMElements(i,j,3,1) =   InproductSpinorReal(Value,Der)
                  if(.not. TSC) NablaMElements(i,j,3,2) = InproductSpinorImaginary(Value,Der)
                endif
            else
                !---------------------------------------------------------------
                !  Calculate all matrix elements 
                Der = DensityBasis(j)%GetDer(1)
                NablaMElements(i,j,1,1) =   InproductSpinorReal(Value,Der)
                if(.not. TSC) NablaMElements(i,j,1,2) = InproductSpinorImaginary(Value,Der)               
                
                Der = DensityBasis(j)%GetDer(2)
                NablaMElements(i,j,2,2) =   InproductSpinorImaginary(Value,Der)
                if(.not. TSC) NablaMElements(i,j,2,1) = InproductSpinorReal(Value,Der)

                Der = DensityBasis(j)%GetDer(3)
                NablaMElements(i,j,3,1) =   InproductSpinorReal(Value,Der)
                if(.not. TSC) NablaMElements(i,j,3,2) = InproductSpinorImaginary(Value,Der)
                !---------------------------------------------------------------
            endif
          enddo
        enddo
     else
        do i=1,nwt
          Value = DensityBasis(i)%GetValue()
          do j=1,nwt
            !-------------------------------------------------------------------
            ! Cycling if not the same isospin or same parity
            !-------------------------------------------------------------------
            if(DensityBasis(i)%GetIsospin().ne.DensityBasis(j)%GetIsospin()) cycle
            if(PC .and. (DensityBasis(i)%GetParity() .eq.DensityBasis(j)%GetParity()))  cycle

            ! If signature is conserved we calculate either (x,y) or (z)
              
            Der = TimeReverse(DensityBasis(j)%GetDer(1))
            NablaMElements(i,j,1,1) =   InproductSpinorReal(Value,Der)
            Der = TimeReverse(DensityBasis(j)%GetDer(2))
            NablaMElements(i,j,2,2) =   InproductSpinorImaginary(Value,Der)
            Der =DensityBasis(j)%GetDer(3)
            NablaMElements(i,j,3,1) =   InproductSpinorReal(Value,Der)
              
            NablaMElements(i,j,1,2) = 0
            NablaMElements(i,j,2,1) = 0
            NablaMElements(i,j,3,2) = 0 
          enddo
        enddo     
     endif
  end subroutine CompNablaMElements
  
  subroutine AlirotateBasis(aliyangle)
    !---------------------------------------------------------------------------
    ! This subroutine alirotates all parity-isospin blocks available with a 
    ! given alirotation in the y-direction. 
    ! 
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: aliyangle(2,2)
    integer                   :: wave, Iindex, Pindex, it, P
    type(spinor)              :: Tphi, phi
    
    if(SC) then
        call stp('Alirotation breaks signature.')
    endif
    
    Iindex = 1 ; if(IC) Iindex = 2
    Pindex = 1 ; if(PC) Pindex = 2
    
    do it=1,Iindex
        do P=1,Pindex
            do wave=1,nwt
                !---------------------------------------------------------------
                ! Check if this particular wave-function is in this
                ! parity-isospin block
                if((HFBasis(wave)%isospin+3)/2 .ne. it) cycle
                if((HFBasis(wave)%parity+3)/2  .ne. P)  cycle
                !---------------------------------------------------------------
                ! Take the time-reverse                 
                Tphi = TimeReverse(HFBasis(wave)%Value)
                
                ! Transform
                phi = cos(aliyangle(P,it)/2) * HFBasis(wave)%Value             &
                &   + sin(aliyangle(P,it)/2) * Tphi
                
                HFBasis(wave)%Value = phi
            enddo        
        enddo
    enddo
    
  end subroutine Alirotatebasis
  
  subroutine NumberParityCounting()
    !---------------------------------------------------------------------------
    ! Correctly determine the quantum number of the many-body state under 
    ! consideration. 
    !  Step a) Collect the spwfs in the DensityBasis with occupation 1.
    !  Step b) Diagonalize the symmetry in this subspace.
    !  Step c) Count the eigenvalues.
    !---------------------------------------------------------------------------

    integer                       :: i,j,k, wave, N(2), it, nw, ref, ifail
    integer                       :: plus(2), minus(2), unclas(2)
    integer, allocatable          :: indices(:,:), isuppz(:), iwork(:)
    complex(KIND=dp), allocatable :: A(:,:,:), work(:)
    real(KIND=dp)                 :: prec = 1d-5, inprod(2)
    real(KIND=dp), allocatable    :: rwork(:), eigen(:,:)
    type(Spinor)                  :: left, right
    
    1 format ( ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ' )
    2 format ( ' Re-diagonalizing S_x among fully occupied levels.           ')
    3 format ( '              +i   -i unclassified     total ')
    4 format ( ' Neutrons: ', 2i5, i5,10x,i5)
    5 format ( ' Protons:  ', 2i5, i5,10x,i5)
    
    allocate(indices(nwt,2)) ; indices = 0
    
    ref = 1
    if(TRC) ref = 2
    
    N = 0
    do wave=1,nwt
        ! Find the spwfs with occupation = 1.
        if(abs(DensityBasis(wave)%Occupation - ref) .lt. prec) then
            it = (DensityBasis(wave)%Isospin + 3)/2
            N(it) = N(it) + 1
            ! Filling the indices array
            indices(N(it),it) = wave
        endif
    enddo
    
    
    ! Populating the array to be diagonalized
    nw = maxval(N)

    allocate(A(nw,nw,2))  ; A = 0.0_dp 
    allocate(eigen(nw,2)) ; eigen = 0
    allocate(isuppz(2*nw))
    allocate(work(2*nw))
    allocate(rwork(2*nw))
    allocate(iwork(2*nw))
    
    do it=1,2
        do j=1,N(it)
           do k=j,N(it)
               if(DensityBasis(indices(j,it))%Parity.ne.  &
               &  DensityBasis(indices(k,it))%Parity      ) cycle
               
               if(DensityBasis(indices(j,it))%Parity.ne.  &
               &  DensityBasis(indices(k,it))%Parity      ) cycle
               
               left      = DensityBasis(indices(j,it))%Value
               right     = ActionOfXSimplex(DensityBasis(indices(k,it))%Value)
               
               ! Note that 
               ! a) S_x is not hermitian, but iS_x is. 
               ! b) Re < j | S_x | k > is restricted by y-timesimplex, and thus
               !                       manually set to zero.
               !    This of course corresponds to i< j | S_x | k > being fully 
               !    real
               inprod(1) = InproductSpinorImaginary(left, right)
               A(j,k,it) = CMPLX(inprod(1), 0) 
           enddo
        enddo
    enddo 

    do it=1,2
       ifail = 0
       call zheev('N','U',N(it), A(1:N(it),1:N(it),it),N(it),eigen(1:N(it),it), work, 2*nw, rwork, 2*nw) 	
    enddo
    
    print 1
    print 2
    print 3
    ! Counting and printing the outcome
    do it=1,2
        plus(it) = 0 ; minus(it) = 0 ; unclas(it) = 0
        do j=1,N(it)
            if( abs(eigen(j,it) - 1)   .lt. 1d-2) then
                plus(it)  = plus(it) + 1
            elseif(abs(eigen(j,it) + 1).lt. 1d-2 ) then
                minus(it) = minus(it) + 1
            else
                unclas(it) = unclas(it) +1 
            endif
        enddo
    enddo
    print 4, plus(1), minus(1), unclas(1), plus(1) + minus(1) + unclas(1)
    print 5, plus(2), minus(2), unclas(2), plus(2) + minus(2) + unclas(2)
    print 1
  end subroutine NumberParityCounting

end module SpwfStorage
