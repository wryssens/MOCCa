module HartreeFock
!
!
!
!
!
!

    use CompilationInfo
    use GenInfo
    use SpwfStorage

    implicit none

    !----------------------------------------------------------
    ! Integers that describe the demanded HF-configuration.
    ! First index
    ! P (1,2) = (-,+)
    ! Second index
    ! Rz(1,2) = (-,+)
    ! Third index
    ! Isospin(1,2) = (N,P)
    integer :: HFConfiguration(2,2,2) = 0
    procedure(PickHFConfig),pointer :: HFFill !=> NaiveFill

    integer :: HFBlock=0
contains

    subroutine PickHFConfig(Configuration)
    !--------------------------
    ! Picks a HF config using the definition employed by HFODD.
    !
    !
    !---------------------------

        integer, intent(in) :: Configuration(2,2,2)
        integer             :: i,it,P,S,j
        integer             :: ConfCopy(2,2,2),order(nwt)

        !Make a copy
        ConfCopy = Configuration

        do i=1,nwt
            call HfBasis(i)%SetOcc(0.0_dp)
        enddo

        !Start counting down
        order = OrderSpwfs(.false.)

        !Loop over the Spwfs in ascending energy order and decrement ConfCopy
        do i=1,nwt
            P = (HFBasis(order(i))%GetParity()    + 3)/2
            it= (HFBasis(order(i))%GetIsospin()   + 3)/2
            S = (HFBasis(order(i))%GetSignature() + 3)/2

            if(ConfCopy(P,S,it) .gt. 0 ) then
                call HFBasis(order(i))%SetOcc(1.0_dp)
                ConfCopy(P,S,it) = ConfCopy(P,S,it) - 1
            endif
        enddo
        !Sanity check
        if(.not. all(ConfCopy .eq. 0)) then
          call stp('PickHFConfig did not find enough orbitals.')
        endif

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

    end subroutine PickHFConfig

    subroutine ReadHFConfig()
        !------------------------------------------------
        ! Read a specific HF configuration from input.
        ! and transforms it to something MOCCa can use.
        !
        !
        !-------------------------------------------------
        integer ::  Npp,Npm,Nmp,Nmm,Ppp,Ppm,Pmp,Pmm

        NameList /HFConfig/ Npp,Npm,Nmp,Nmm,Ppp,Ppm,Pmp,Pmm

        read(*,nml=HFConfig)

        HFConfiguration(1,1,1) = Nmm
        HFConfiguration(1,1,2) = Pmm

        HFConfiguration(2,1,1) = Npm
        HFConfiguration(2,1,2) = Ppm

        HFConfiguration(1,2,1) = Nmp
        HFConfiguration(1,2,2) = Pmp

        HFConfiguration(2,2,1) = Npp
        HFConfiguration(2,2,2) = Ppp

        !Sanity check
        if(.not. TRC) then
            if (sum(HFConfiguration(:,:,1)).ne.int(Neutrons)) then
                call stp('The hfconfiguration you wanted does not coincide with' &
                &      //' the neutron number you wanted.',                      &
                &        'HFConfig', sum(HFConfiguration(:,:,1)),                &
                &        'Neutrons', Neutrons)
            endif

            if (sum(HFConfiguration(:,:,2)).ne.int(Protons)) then
                call stp('The hfconfiguration you wanted does not coincide with' &
                &      //' the neutron number you wanted.',                      &
                &        'HFConfig', sum(HFConfiguration(:,:,2)),                &
                &        'Protons', Protons)
            endif   
        else
            if (2*sum(HFConfiguration(:,:,1)).ne.int(Neutrons)) then
                call stp('The hfconfiguration you wanted does not coincide with' &
                &      //' the neutron number you wanted.',                      &
                &        'HFConfig', 2*sum(HFConfiguration(:,:,1)),                &
                &        'Neutrons', Neutrons)
            endif

            if (2*sum(HFConfiguration(:,:,2)).ne.int(Protons)) then
                call stp('The hfconfiguration you wanted does not coincide with' &
                &      //' the neutron number you wanted.',                      &
                &        'HFConfig', 2*sum(HFConfiguration(:,:,2)),                &
                &        'Protons', Protons)
            endif   
        endif
    end subroutine ReadHFConfig

  function HFEnergy( Delta) result(E)
    !---------------------------------------------------------------------------
    ! Returns the pairing energy of a HF calculation, namely 0.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: E(2)
    complex(KIND=dp), allocatable,intent(in) :: Delta(:,:,:,:)

    E = 0.0_dp

    return
  end function HFEnergy

  subroutine NaiveFill(Configuration)
    !---------------------------------------------------------------------------
    ! This subroutine finds the orbitals with the lowest single particle
    ! energy and fills them, after sorting all the levels.
    ! Of course, nothing happens if the user has asked to freeze the occupation.
    !---------------------------------------------------------------------------
    integer :: i, n,p, ProtonUpperBound, NeutronUpperBound, order(nwt),j
    integer, intent(in) :: Configuration(2,2,2)

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
      if(HFBlock .ne. 0 ) then
        NeutronUpperBound = NeutronUpperBound -1
      endif
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
    ! Blocking a certain state
    if(HFBlock.ne.0) then
        call HFBasis(HFBlock)%SetOcc(1.0_dp)
    endif
    
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
  end subroutine NaiveFill

!   subroutine HFFermi
!   !-----------------------------------------------------------------------------
!   ! Find the fermi energy of a HF calculation. Trivial, since it is just the
!   ! energy of the highest occupied level.
!   !-----------------------------------------------------------------------------
!     integer       :: i, it
!     real(KIND=dp) :: E

!     Fermi = -1000000.0_dp
!     do i=1,nwt
!         it = (HFBasis(i)%GetIsospin() + 3)/2
!         E  = HFBasis(i)%GetEnergy()
!         if(HFBasis(i)%GetOcc().gt.0.0_dp.and.Fermi(it) .lt. E)then
!             Fermi(it) = E
!         endif
!     enddo

!   end subroutine HFFermi



end module HartreeFock
