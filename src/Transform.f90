module Transform
  !-----------------------------------------------------------------------------
  ! This module contains all the routines that transform wavefunctions in
  ! relevant ways. Breaking symmetries, translating wavefunctions etc...
  ! One single routine is also included that transforms densities. They are much
  ! easier to deal with.
  !-----------------------------------------------------------------------------
  ! Note that this module uses Coulombderivatives instead of ordinary ones,since
  ! they are more flexible. More precisely, the PointToLineExtension Routines in
  ! the Derivatives module work with fixed nx,ny,nz for efficiency reasons and
  ! are thus not reliable in the routines in this module since we are dealing
  ! with modified nx/ny/nz. The Coulomb equivalents for PointToLineExtension
  ! work with implicit array sizes and thus are better suited.
  ! This is absolutely not a textbook example of programming, but it is slightly
  ! more efficient than the textbook example :D
  !-----------------------------------------------------------------------------
  use CompilationInfo
  use GenInfo
  use Spinors
  use Densities
  use WaveFunctions
  use Spwfstorage
  use Spinors

  use CoulombDerivatives

  implicit none

  real(KIND=dp), allocatable :: interpolX(:,:,:), interpolY(:,:,:), interpolZ(:,:,:)

contains

  subroutine TransformInput(inTRC,inTSC,inIC,inPC,inSC,                        &
  &                         filenx,fileny,filenz,filenwt, filedx)
  !-----------------------------------------------------------------------------
  ! Transform the data read from the wavefunction file to a form appropriate
  ! for the different symmetry combinations. (This includes the densities!)
  !
  !-----------------------------------------------------------------------------
  logical, intent(in)       :: inTRC, inTSC, inIC, inPC, inSC
  integer, intent(in)       :: filenx,fileny,filenz,filenwt
  real(KIND=dp), intent(in) :: filedx
  type(Spwf), allocatable   :: HFBasisTransformed(:)
  integer                   :: wave, index
  logical                   :: interpolx, interpoly, interpolz

    !---------------------------------------------------------------------------
    ! Checking for Time Reversal breaking first!
    ! Note that we try to take the same ordering as CR8: first all the neutron
    ! wavefunctions, then all the proton wavefunctions.
    if(inTRC.eqv.TRC) then
      allocate(HFBasisTransformed(nwt))
      do wave=1,nwt
        HFBasisTransformed(wave) = CopyWaveFunction(HFBasis(wave))
      enddo
    else
      allocate(HFBasisTransformed(nwt))
      ! Count the neutron wavefunctions
      do wave=1,filenwt
        if(HFBasis(wave)%GetIsospin().gt.0) exit
      enddo
      index = wave - 1
      ! Create the neutron time-reversed pairs
      do wave=1,index
        call BreakTimeReversal(HFBasis(wave), HFBasisTransformed(wave),        &
        &                                     HFBasisTransformed(wave+index) )
      enddo

      ! Create the proton time-reversed pairs
      do wave=index+1,filenwt
        call BreakTimeReversal(HFBasis(wave),HFBasisTransformed(wave+index),   &
        &                                    HFBasisTransformed(wave+filenwt))
      enddo
    endif
    !---------------------------------------------------------------------------
    !Checking for Parity Breaking
    if(inPC.neqv.PC) then
      !Breaking Parity
      do wave=1,size(HFBasisTransformed)
        call BreakParity(HFBasisTransformed(wave) )
      enddo
    endif
    !---------------------------------------------------------------------------
    !Checking for Signature Breaking
    if(inSC.neqv.SC) then
      !Breaking Signature
      do wave=1,size(HFBasisTransformed)
        call BreakSignature(HFBasisTransformed(wave))
      enddo
    endif
    !---------------------------------------------------------------------------
    !Checking for Time Simplex Breaking
    if(inTSC.neqv.TSC) then
      !Breaking Parity
      do wave=1,size(HFBasisTransformed)
        call BreakTimeSimplex(HFBasisTransformed(wave))
      enddo
    endif
    !---------------------------------------------------------------------------
    !Applying the changes to the HFBasis
    call ChangeNumberWaveFunctions(nwt)
    do wave=1,nwt
      HFBasis(wave) = HFBasisTransformed(wave)
      call HFBasis(wave)%SymmetryOperators()
    enddo    
    !---------------------------------------------------------------------------
    ! Transforming the densities
    call TransformDensities(inSC, inTSC, inPC , filenx, fileny, filenz)
    !---------------------------------------------------------------------------
    ! Interpolation 
    interpolx = .false.
    interpoly = .false.
    interpolz = .false.

    if(( SC .eqv.  inSC  ).and.( nx.ne.  filenx)) interpolx = .true.
    if(( SC .neqv. inSC  ).and.( nx.ne.2*filenx)) interpolx = .true.

    if((TSC .eqv.  inTSC ).and.( ny.ne.  fileny)) interpoly = .true.
    if((TSC .neqv. inTSC ).and.( ny.ne.2*fileny)) interpoly = .true.

    if(( PC .eqv.  inPC  ).and.( nz.ne.  filenz)) interpolz = .true.
    if(( PC .neqv. inPC  ).and.( nz.ne.2*filenz)) interpolz = .true.

    if(interpolx .or. interpoly .or. interpolz) then
      call ConstructInterpolationFunctions(size(HFBasis(1)%Value%Grid,1)         &
      &                                   ,size(HFBasis(1)%Value%Grid,2)         &
      &                                   ,size(HFBasis(1)%Value%Grid,3)         &
      &                                   ,filedx)
      ! Interpolating the wavefunctions
      do wave=1,nwt
        HFBasis(wave) = InterpolateSpwf(HFBasis(wave))
      enddo
    endif
    !---------------------------------------------------------------------------
    ! Make sure the densities get recalculated
    Density = NewDensityVector()
    Recalc=.true.
  
  end subroutine TransformInput

  subroutine BreakParity(wf)
   !---------------------------------------------------------------------------
   ! This subroutine takes as argument a single particle wavefunction and
   ! returns the parity-broken version of that spwf.
   ! Note that this subroutine doubles the parameter nz!
    !---------------------------------------------------------------------------
    type(Spwf), intent(inout) :: wf
    type(Spinor)              :: Temp, Value
    integer                   :: p, TimeSimplex, Parity, Signature,TimeReversal
    integer                   :: i,j,k,l, s, oldnz
    real(KIND=dp)             :: Energy, Occ, Dispersion, AngMom(3)

    !Setting up the parameters of the wavefunction we are dealing with
    Parity      = wf%GetParity()
    TimeSimplex = wf%GetTimeSimplex()
    Signature   = wf%GetSignature()
    TimeReversal= wf%GetTimeReversal()
    Energy      = wf%GetEnergy()
    Occ         = wf%GetOcc()
    Dispersion  = wf%GetDispersion()

    do i=1,3
        AngMom(i) = wf%GetAngMoment(i)
    enddo

    Value = wf%GetValue()

    !Checking input
    if(Parity.eq.0) call stp("Parity is already broken!")

    !When breaking parity, we need to store the values for the entire z-axis,
    !instead of only half.
    !Note that the x and y dimensions need not be nx/ny when we are breaking
    !more than one symmetry.
    oldnz = size(Value%Grid,3)
    Temp  = NewSizeSpinor(size(Value%Grid,1), size(Value%Grid,2), 2*oldnz)
    ! Moving the old values to their correct place in the new spinor.
    Temp%Grid(:,:,oldnz+1:2*oldnz ,:,:) = Value%Grid(:,:,1:oldnz,:,:)

    !Using the symmetries of the wavefunction to construct the rest
    do l=1,4
      s = CompSignExtension(3,Parity,Signature,TimeSimplex,l)
      do j=1, size(Value%Grid,2)
        do i=1,size(Value%Grid,1)
          Temp%Grid(i,j,1:oldnz,l,1) =                                         &
          & s*Reverse(LineExtensionCoulombZ(Value%Grid(:,:,:,l,1), &
          & oldnz,Parity,TimeSimplex,Signature,i,j))
        enddo
      enddo
    enddo
    !Creating a new wavefunction, with the same quantum numbers,
    ! except for the Parity, which is now 0
    wf=NewWaveFunction(Temp,wf%GetIsoSpin(),TimeSimplex,0,Signature,TimeReversal)

    !Adding the other parameters
    call wf%SetEnergy(Energy)
    call wf%SetOcc(Occ)
    call wf%SetDispersion(Dispersion)
    call wf%SetAngMoment(AngMom)
  end subroutine BreakParity

  subroutine BreakSignature(wf)
    !---------------------------------------------------------------------------
    ! This subroutine takes as argument a single particle wavefunction and
    ! returns the signature-broken version of that spwf.
    !---------------------------------------------------------------------------
    use Spinors
    type(Spwf), intent(inout) :: wf
    type(Spinor)              :: Temp, Value
    integer                   :: p, TimeSimplex, Parity, Signature, TimeReversal
    integer                   :: i,j,k,l,s, oldnx
    real(KIND=dp)             :: AngMom(3), Energy, Occ, Dispersion

    !Setting up the parameters of the wavefunction we are dealing with
    Parity      = wf%GetParity()
    TimeSimplex = wf%GetTimeSimplex()
    Signature   = wf%GetSignature()
    TimeReversal= wf%GetTimeReversal()
    Energy      = wf%GetEnergy()
    Occ         = wf%GetOcc()
    Dispersion  = wf%GetDispersion()

    do i=1,3
      AngMom(i) = wf%GetAngMoment(i)
    enddo

    Value = wf%GetValue()

    !Checking input
    if(Signature.eq.0) call stp("Signature is already broken!")

    ! When breaking parity, we need to store the values for the entire x-axis,
    ! instead of only half.
    oldnx = size(Value%Grid,1)
    Temp  = NewSizeSpinor(2*oldnx,size(Value%Grid,2), size(Value%Grid,3))

    ! Moving the old values to their correct place in the new spinor.
    Temp%Grid(oldnx+1:2*oldnx,:,:,:,:) = Value%Grid(1:oldnx,:,:,:,:)

    !Using the symmetries of the wavefunction to construct the rest.
    do l=1,4
      s = CompSignExtension(1,Parity,Signature,TimeSimplex,l)
      do k=1,size(Value%Grid,3)
        do j=1,size(Value%Grid,2)
          Temp%Grid(1:oldnx,j,k,l,1) = real(s, KIND=dp)* &
          & Reverse(LineExtensionCoulombX(Value%Grid(:,:,:,l,1), &
          & oldnx,Parity,TimeSimplex,Signature,j,k))
        enddo
      enddo
    enddo
    !Creating a new wavefunction, with the same quantum numbers
    Wf=NewWaveFunction(Temp, wf%GetIsoSpin(),TimeSimplex,Parity,0,TimeReversal)

    !Adding the other parameters
    call wf%SetEnergy(Energy)
    call wf%SetOcc(Occ)
    call wf%SetDispersion(Dispersion)
    call wf%SetAngMoment(AngMom)
  end subroutine BreakSignature

  subroutine BreakTimeSimplex(wf)
    !---------------------------------------------------------------------------
    ! This subroutine takes as argument a single particle wavefunction and
    ! returns the timesimplex-broken version of that spwf.
    !---------------------------------------------------------------------------
    use Spinors
    type(Spwf), intent(inout) :: wf
    type(Spinor)              :: Temp, Value
    integer                   :: p, TimeSimplex, Parity, Signature, TimeReversal
    integer                   :: i,j,k,l,s, oldny
    real(KIND=dp)             :: Occ, Energy,Dispersion, AngMom(3)

    !Setting up the parameters of the wavefunction we are dealing with
    Parity      = wf%GetParity()
    TimeSimplex = wf%GetTimeSimplex()
    Signature   = wf%GetSignature()
    TimeReversal= wf%GetTimeReversal()
    Energy      = wf%GetEnergy()
    Occ         = wf%GetOcc()
    Dispersion  = wf%GetDispersion()

    do i=1,3
      AngMom(i) = wf%GetAngMoment(i)
    enddo

    Value = wf%GetValue()

    !Checking input
    if(TimeSimplex.eq.0) call stp("Time-simplex is already broken!")

    oldny = size(Value%Grid,2)
    !When breaking parity, we need to store the values for the entire y-axis.
    Temp=NewSizeSpinor(size(Value%Grid,1),2*oldny,size(Value%Grid,3))

    ! Moving the old values to their correct place in the new spinor.
    Temp%Grid(:, oldny+1:2*oldny,:,:,:) = Value%Grid(:,1:oldny,:,:,:)
   
    !Using the symmetries of the wavefunction to construct the rest
    do l=1,4
      s = CompSignExtension(2,Parity,Signature,TimeSimplex,l)
      do k=1,size(Value%Grid,3)
        do i=1,size(Value%Grid,1)
          Temp%Grid(i,1:oldny,k,l,1) = s*Reverse( &
          & LineExtensionCoulombY(Value%Grid(:,:,:,l,1), &
            & oldny,Parity,TimeSimplex,Signature,i,k))
        enddo
      enddo
    enddo
    ! Creating a new wavefunction, with the same quantum numbers,
    ! except for the Time Simplex, which is now0.
    wf = NewWaveFunction(Temp, wf%GetIsoSpin(),0,Parity,Signature,TimeReversal)

    !Adding the other parameters
    call wf%SetEnergy(Energy)
    call wf%SetOcc(Occ)
    call wf%SetDispersion(Dispersion)
    call wf%SetAngMoment(AngMom)
  end subroutine BreakTimeSimplex

  subroutine TransformDensities(inSC, inTSC, inPC, filenx,fileny,filenz)
  !-----------------------------------------------------------------------------
  ! Transform the densities as they were read on file to ones that are usable
  ! by the program.
  !
  ! WR: This is BY FAR the most ugly and cumbersome routine I've ever written...
  !-----------------------------------------------------------------------------
    use Densities, only : DensityVector, NewDensityVector, Density, FillOct
    use Derivatives

    type(DensityVector)       :: NewDen
    logical, intent(in)       :: inPC, inTSC, inSC
    integer, intent(in)       :: filenx,fileny,filenz
    integer                   :: startx,endx,starty,endy,startz,endz, i,j,k
    integer                   :: x,y,z, Oct

    ! Do nothing
    if((TSC.eqv.inTSC) .and. (PC.eqv.inPC) .and. (SC.eqv.TSC)) return

    !Create new densityvectortype with normal (nx,ny,nz) size
    NewDen = NewDensityVector()

    !---------------------------------------------------------------------------
    ! We divide the transformation by Octants.
    ! For the signs, best to see
    ! V.Hellemans, P.H. Heenen and M. Bender, Phys. Rev. C, 85, 014326 (2012)
    ! Appendix B, table 5
    !
    ! Octant(x >0 , y > 0, z > 0) (This is ofcourse just copying...)
    Oct = 1
    call FillOct(NewDen%Rho              , Density%Rho              , Oct,+1)
    call FillOct(NewDen%DerRho(:,:,:,1,:), Density%DerRho(:,:,:,1,:), Oct,+1)
    call FillOct(NewDen%DerRho(:,:,:,2,:), Density%DerRho(:,:,:,2,:), Oct,+1)
    call FillOct(NewDen%DerRho(:,:,:,3,:), Density%DerRho(:,:,:,3,:), Oct,+1)
    call FillOct(NewDen%LapRho           , Density%LapRho           , Oct,+1)
    call FillOct(NewDen%Tau              , Density%Tau              , Oct,+1)
    call FillOct(NewDen%NablaJ           , Density%NablaJ           , Oct,+1)

    if(allocated(NewDen%JMuNu)) then
      call FillOct(NewDen%JMunu(:,:,:,1,1,:), Density%JMuNu(:,:,:,1,1,:),Oct,1)
      call FillOct(NewDen%JMunu(:,:,:,1,2,:), Density%JMuNu(:,:,:,1,2,:),Oct,1)
      call FillOct(NewDen%JMunu(:,:,:,1,3,:), Density%JMuNu(:,:,:,1,3,:),Oct,1)
      call FillOct(NewDen%JMunu(:,:,:,2,1,:), Density%JMuNu(:,:,:,2,1,:),Oct,1)
      call FillOct(NewDen%JMunu(:,:,:,2,2,:), Density%JMuNu(:,:,:,2,2,:),Oct,1)
      call FillOct(NewDen%JMunu(:,:,:,2,3,:), Density%JMuNu(:,:,:,2,3,:),Oct,1)
      call FillOct(NewDen%JMunu(:,:,:,3,1,:), Density%JMuNu(:,:,:,3,1,:),Oct,1)
      call FillOct(NewDen%JMunu(:,:,:,3,2,:), Density%JMuNu(:,:,:,3,2,:),Oct,1)
      call FillOct(NewDen%JMunu(:,:,:,3,3,:), Density%JMuNu(:,:,:,3,3,:),Oct,1)
    endif

    if(.not.TRC) then
      call FillOct(NewDen%vecj(:,:,:,1,:), Density%vecj(:,:,:,1,:),Oct,1)
      call FillOct(NewDen%vecj(:,:,:,2,:), Density%vecj(:,:,:,2,:),Oct,1)
      call FillOct(NewDen%vecj(:,:,:,3,:), Density%vecj(:,:,:,3,:),Oct,1)

      call FillOct(NewDen%vecs(:,:,:,1,:), Density%vecs(:,:,:,1,:),Oct,1)
      call FillOct(NewDen%vecs(:,:,:,2,:), Density%vecs(:,:,:,2,:),Oct,1)
      call FillOct(NewDen%vecs(:,:,:,3,:), Density%vecs(:,:,:,3,:),Oct,1)

      call FillOct(NewDen%rots(:,:,:,1,:), Density%rots(:,:,:,1,:),Oct,1)
      call FillOct(NewDen%rots(:,:,:,2,:), Density%rots(:,:,:,2,:),Oct,1)
      call FillOct(NewDen%rots(:,:,:,3,:), Density%rots(:,:,:,3,:),Oct,1)

      call FillOct(NewDen%rotvecj(:,:,:,1,:), Density%rotvecj(:,:,:,1,:),Oct,1)
      call FillOct(NewDen%rotvecj(:,:,:,2,:), Density%rotvecj(:,:,:,2,:),Oct,1)
      call FillOct(NewDen%rotvecj(:,:,:,3,:), Density%rotvecj(:,:,:,3,:),Oct,1)

      call FillOct(NewDen%ders(:,:,:,1,1,:), Density%ders(:,:,:,1,1,:),Oct,1)
      call FillOct(NewDen%ders(:,:,:,1,2,:), Density%ders(:,:,:,1,2,:),Oct,1)
      call FillOct(NewDen%ders(:,:,:,1,3,:), Density%ders(:,:,:,1,3,:),Oct,1)
      call FillOct(NewDen%ders(:,:,:,2,1,:), Density%ders(:,:,:,2,1,:),Oct,1)
      call FillOct(NewDen%ders(:,:,:,2,2,:), Density%ders(:,:,:,2,2,:),Oct,1)
      call FillOct(NewDen%ders(:,:,:,2,3,:), Density%ders(:,:,:,2,3,:),Oct,1)
      call FillOct(NewDen%ders(:,:,:,3,1,:), Density%ders(:,:,:,3,1,:),Oct,1)
      call FillOct(NewDen%ders(:,:,:,3,2,:), Density%ders(:,:,:,3,2,:),Oct,1)
      call FillOct(NewDen%ders(:,:,:,3,3,:), Density%ders(:,:,:,3,3,:),Oct,1)

      if(allocated(NewDen%LapS)) then
        call FillOct(NewDen%LapS(:,:,:,1,:), Density%LapS(:,:,:,1,:),Oct,1)
        call FillOct(NewDen%LapS(:,:,:,2,:), Density%LapS(:,:,:,2,:),Oct,1)
        call FillOct(NewDen%LapS(:,:,:,3,:), Density%LapS(:,:,:,3,:),Oct,1)

        call FillOct(NewDen%graddivS(:,:,:,1,:),Density%gradDivS(:,:,:,1,:),Oct,1)
        call FillOct(NewDen%graddivS(:,:,:,2,:),Density%gradDivS(:,:,:,2,:),Oct,1)
        call FillOct(NewDen%graddivS(:,:,:,3,:),Density%graddivS(:,:,:,3,:),Oct,1)

        call FillOct(NewDen%divs(:,:,:,:), Density%divs(:,:,:,:),Oct,1)
      endif

      if(allocated(NewDen%VecF)) then
        call FillOct(NewDen%vecF(:,:,:,1,:), Density%vecF(:,:,:,1,:),Oct,1)
        call FillOct(NewDen%vecF(:,:,:,2,:), Density%vecF(:,:,:,2,:),Oct,1)
        call FillOct(NewDen%vecF(:,:,:,3,:), Density%vecF(:,:,:,3,:),Oct,1)
      endif

      if(allocated(NewDen%vecT)) then
        call FillOct(NewDen%vecT(:,:,:,1,:), Density%vecT(:,:,:,1,:),Oct,1)
        call FillOct(NewDen%vecT(:,:,:,2,:), Density%vecT(:,:,:,2,:),Oct,1)
        call FillOct(NewDen%vecT(:,:,:,3,:), Density%vecT(:,:,:,3,:),Oct,1)
      endif
    endif
    ! End of octant 1
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Octant ( x < 0 , y > 0, z > 0)
    ! Necessary when Signature is not conserved
    if(inSC.neqv.SC) then
      Oct = 2
      call FillOct(NewDen%Rho              , Density%Rho              , Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,1,:), Density%DerRho(:,:,:,1,:), Oct,-1)
      call FillOct(NewDen%DerRho(:,:,:,2,:), Density%DerRho(:,:,:,2,:), Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,3,:), Density%DerRho(:,:,:,3,:), Oct,+1)
      call FillOct(NewDen%LapRho           , Density%LapRho           , Oct,+1)
      call FillOct(NewDen%Tau              , Density%Tau              , Oct,+1)
      call FillOct(NewDen%NablaJ           , Density%NablaJ           , Oct,+1)

      if(allocated(NewDen%JMuNu)) then
        call FillOct(NewDen%JMunu(:,:,:,1,1,:),Density%JMuNu(:,:,:,1,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,1,2,:),Density%JMuNu(:,:,:,1,2,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,1,3,:),Density%JMuNu(:,:,:,1,3,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,1,:),Density%JMuNu(:,:,:,2,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,2,:),Density%JMuNu(:,:,:,2,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,3,:),Density%JMuNu(:,:,:,2,3,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,1,:),Density%JMuNu(:,:,:,3,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,3,2,:),Density%JMuNu(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,3,:),Density%JMuNu(:,:,:,3,3,:),Oct,-1)
      endif

      if(.not.TRC) then
        call FillOct(NewDen%vecj(:,:,:,1,:),Density%vecj(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%vecj(:,:,:,2,:),Density%vecj(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%vecj(:,:,:,3,:),Density%vecj(:,:,:,3,:),Oct,-1)

        call FillOct(NewDen%vecs(:,:,:,1,:), Density%vecs(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%vecs(:,:,:,2,:), Density%vecs(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%vecs(:,:,:,3,:), Density%vecs(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rots(:,:,:,1,:), Density%rots(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%rots(:,:,:,2,:), Density%rots(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%rots(:,:,:,3,:), Density%rots(:,:,:,3,:),Oct,-1)

        call FillOct(NewDen%rotvecj(:,:,:,1,:), Density%rotvecj(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%rotvecj(:,:,:,2,:), Density%rotvecj(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%rotvecj(:,:,:,3,:), Density%rotvecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%ders(:,:,:,1,1,:), Density%ders(:,:,:,1,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,1,2,:), Density%ders(:,:,:,1,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,1,3,:), Density%ders(:,:,:,1,3,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,1,:), Density%ders(:,:,:,2,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,2,:), Density%ders(:,:,:,2,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,3,:), Density%ders(:,:,:,2,3,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,1,:), Density%ders(:,:,:,3,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,2,:), Density%ders(:,:,:,3,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,3,:), Density%ders(:,:,:,3,3,:),Oct,+1)

        if(allocated(NewDen%LapS)) then
          call FillOct(NewDen%LapS(:,:,:,1,:), Density%LapS(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%LapS(:,:,:,2,:), Density%LapS(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%LapS(:,:,:,3,:), Density%LapS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%graddivS(:,:,:,1,:),Density%gradDivS(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%graddivS(:,:,:,2,:),Density%gradDivS(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%graddivS(:,:,:,3,:),Density%graddivS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%divs(:,:,:,:), Density%divs(:,:,:,:),Oct,+1)
        endif

        if(allocated(NewDen%VecF)) then
          call FillOct(NewDen%vecF(:,:,:,1,:), Density%vecF(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%vecF(:,:,:,2,:), Density%vecF(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%vecF(:,:,:,3,:), Density%vecF(:,:,:,3,:),Oct,+1)
        endif

        if(allocated(NewDen%vecT)) then
          call FillOct(NewDen%vecT(:,:,:,1,:), Density%vecT(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%vecT(:,:,:,2,:), Density%vecT(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%vecT(:,:,:,3,:), Density%vecT(:,:,:,3,:),Oct,+1)
        endif
      endif
    endif
    ! End of Octant 2
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Octant (x > 0, y < 0, z > 0)
    ! Necessary when TimeSimplex is not conserved
    if(inTSC.neqv.TSC) then
      Oct = 3
      call FillOct(NewDen%Rho              , Density%Rho              , Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,1,:), Density%DerRho(:,:,:,1,:), Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,2,:), Density%DerRho(:,:,:,2,:), Oct,-1)
      call FillOct(NewDen%DerRho(:,:,:,3,:), Density%DerRho(:,:,:,3,:), Oct,+1)
      call FillOct(NewDen%LapRho           , Density%LapRho           , Oct,+1)
      call FillOct(NewDen%Tau              , Density%Tau              , Oct,+1)
      call FillOct(NewDen%NablaJ           , Density%NablaJ           , Oct,+1)

      if(allocated(NewDen%JMuNu)) then
        call FillOct(NewDen%JMunu(:,:,:,1,1,:),Density%JMuNu(:,:,:,1,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,1,2,:),Density%JMuNu(:,:,:,1,2,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,1,3,:),Density%JMuNu(:,:,:,1,3,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,1,:),Density%JMuNu(:,:,:,2,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,2,:),Density%JMuNu(:,:,:,2,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,3,:),Density%JMuNu(:,:,:,2,3,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,3,1,:),Density%JMuNu(:,:,:,3,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,2,:),Density%JMuNu(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,3,:),Density%JMuNu(:,:,:,3,3,:),Oct,-1)
      endif

      if(.not.TRC) then
        call FillOct(NewDen%vecj(:,:,:,1,:),Density%vecj(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%vecj(:,:,:,2,:),Density%vecj(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%vecj(:,:,:,3,:),Density%vecj(:,:,:,3,:),Oct,-1)

        call FillOct(NewDen%vecs(:,:,:,1,:), Density%vecs(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%vecs(:,:,:,2,:), Density%vecs(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%vecs(:,:,:,3,:), Density%vecs(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rots(:,:,:,1,:), Density%rots(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%rots(:,:,:,2,:), Density%rots(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%rots(:,:,:,3,:), Density%rots(:,:,:,3,:),Oct,-1)

        call FillOct(NewDen%rotvecj(:,:,:,1,:), Density%rotvecj(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%rotvecj(:,:,:,2,:), Density%rotvecj(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%rotvecj(:,:,:,3,:), Density%rotvecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%ders(:,:,:,1,1,:), Density%ders(:,:,:,1,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,1,2,:), Density%ders(:,:,:,1,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,1,3,:), Density%ders(:,:,:,1,3,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,1,:), Density%ders(:,:,:,2,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,2,:), Density%ders(:,:,:,2,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,3,:), Density%ders(:,:,:,2,3,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,1,:), Density%ders(:,:,:,3,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,2,:), Density%ders(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,3,:), Density%ders(:,:,:,3,3,:),Oct,+1)

        if(allocated(NewDen%LapS)) then
          call FillOct(NewDen%LapS(:,:,:,1,:), Density%LapS(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%LapS(:,:,:,2,:), Density%LapS(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%LapS(:,:,:,3,:), Density%LapS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%graddivS(:,:,:,1,:),Density%gradDivS(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%graddivS(:,:,:,2,:),Density%gradDivS(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%graddivS(:,:,:,3,:),Density%graddivS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%divs(:,:,:,:), Density%divs(:,:,:,:),Oct,+1)
        endif

        if(allocated(NewDen%VecF)) then
          call FillOct(NewDen%vecF(:,:,:,1,:), Density%vecF(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%vecF(:,:,:,2,:), Density%vecF(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%vecF(:,:,:,3,:), Density%vecF(:,:,:,3,:),Oct,+1)
        endif

        if(allocated(NewDen%vecT)) then
          call FillOct(NewDen%vecT(:,:,:,1,:), Density%vecT(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%vecT(:,:,:,2,:), Density%vecT(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%vecT(:,:,:,3,:), Density%vecT(:,:,:,3,:),Oct,+1)
        endif
      endif
    endif
    ! End of Octant 3
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Octant (x > 0, y > 0, z < 0)
    ! Necessary when Parity is broken
    if(PC.neqv.InPC) then
      Oct = 4
      call FillOct(NewDen%Rho              , Density%Rho              , Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,1,:), Density%DerRho(:,:,:,1,:), Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,2,:), Density%DerRho(:,:,:,2,:), Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,3,:), Density%DerRho(:,:,:,3,:), Oct,-1)
      call FillOct(NewDen%LapRho           , Density%LapRho           , Oct,+1)
      call FillOct(NewDen%Tau              , Density%Tau              , Oct,+1)
      call FillOct(NewDen%NablaJ           , Density%NablaJ           , Oct,+1)

      if(allocated(NewDen%JMuNu)) then
        call FillOct(NewDen%JMunu(:,:,:,1,1,:), Density%JMuNu(:,:,:,1,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,1,2,:), Density%JMuNu(:,:,:,1,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,1,3,:), Density%JMuNu(:,:,:,1,3,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,1,:), Density%JMuNu(:,:,:,2,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,2,:), Density%JMuNu(:,:,:,2,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,3,:), Density%JMuNu(:,:,:,2,3,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,3,1,:), Density%JMuNu(:,:,:,3,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,3,2,:), Density%JMuNu(:,:,:,3,2,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,3,3,:), Density%JMuNu(:,:,:,3,3,:),Oct,-1)
      endif

      if(.not.TRC) then
        call FillOct(NewDen%vecj(:,:,:,1,:), Density%vecj(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%vecj(:,:,:,2,:), Density%vecj(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%vecj(:,:,:,3,:), Density%vecj(:,:,:,3,:),Oct,-1)

        call FillOct(NewDen%vecs(:,:,:,1,:), Density%vecs(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%vecs(:,:,:,2,:), Density%vecs(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%vecs(:,:,:,3,:), Density%vecs(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rots(:,:,:,1,:), Density%rots(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%rots(:,:,:,2,:), Density%rots(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%rots(:,:,:,3,:), Density%rots(:,:,:,3,:),Oct,-1)

        call FillOct(NewDen%rotvecj(:,:,:,1,:), Density%rotvecj(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%rotvecj(:,:,:,2,:), Density%rotvecj(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%rotvecj(:,:,:,3,:), Density%rotvecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%ders(:,:,:,1,1,:), Density%ders(:,:,:,1,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,1,2,:), Density%ders(:,:,:,1,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,1,3,:), Density%ders(:,:,:,1,3,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,1,:), Density%ders(:,:,:,2,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,2,:), Density%ders(:,:,:,2,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,3,:), Density%ders(:,:,:,2,3,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,1,:), Density%ders(:,:,:,3,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,2,:), Density%ders(:,:,:,3,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,3,:), Density%ders(:,:,:,3,3,:),Oct,-1)

        if(allocated(NewDen%LapS)) then
          call FillOct(NewDen%LapS(:,:,:,1,:), Density%LapS(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%LapS(:,:,:,2,:), Density%LapS(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%LapS(:,:,:,3,:), Density%LapS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%graddivS(:,:,:,1,:),Density%gradDivS(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%graddivS(:,:,:,2,:),Density%gradDivS(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%graddivS(:,:,:,3,:),Density%graddivS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%divs(:,:,:,:), Density%divs(:,:,:,:),Oct,-1)
        endif

        if(allocated(NewDen%VecF)) then
          call FillOct(NewDen%vecF(:,:,:,1,:), Density%vecF(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%vecF(:,:,:,2,:), Density%vecF(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%vecF(:,:,:,3,:), Density%vecF(:,:,:,3,:),Oct,+1)
        endif

        if(allocated(NewDen%vecT)) then
          call FillOct(NewDen%vecT(:,:,:,1,:), Density%vecT(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%vecT(:,:,:,2,:), Density%vecT(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%vecT(:,:,:,3,:), Density%vecT(:,:,:,3,:),Oct,+1)
        endif
      endif
    endif
    ! End of Octant 4
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Octant ( x < 0 , y < 0 , z > 0)
    ! Necessary when signature and time simplex are conserved on file, but both
    ! broken in the calculation
    if((inSC.neqv.SC) .and. (inTSC.neqv.inTSC)) then
      Oct = 5
      call FillOct(NewDen%Rho              , Density%Rho              , Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,1,:), Density%DerRho(:,:,:,1,:), Oct,-1)
      call FillOct(NewDen%DerRho(:,:,:,2,:), Density%DerRho(:,:,:,2,:), Oct,-1)
      call FillOct(NewDen%DerRho(:,:,:,3,:), Density%DerRho(:,:,:,3,:), Oct,+1)
      call FillOct(NewDen%LapRho           , Density%LapRho           , Oct,+1)
      call FillOct(NewDen%Tau              , Density%Tau              , Oct,+1)
      call FillOct(NewDen%NablaJ           , Density%NablaJ           , Oct,+1)

      if(allocated(NewDen%JMuNu)) then
        call FillOct(NewDen%JMunu(:,:,:,1,1,:), Density%JMuNu(:,:,:,1,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,1,2,:), Density%JMuNu(:,:,:,1,2,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,1,3,:), Density%JMuNu(:,:,:,1,3,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,1,:), Density%JMuNu(:,:,:,2,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,2,:), Density%JMuNu(:,:,:,2,2,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,3,:), Density%JMuNu(:,:,:,2,3,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,1,:), Density%JMuNu(:,:,:,3,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,2,:), Density%JMuNu(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,3,:), Density%JMuNu(:,:,:,3,3,:),Oct,+1)
      endif

      if(.not.TRC) then
        call FillOct(NewDen%vecj(:,:,:,1,:), Density%vecj(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%vecj(:,:,:,2,:), Density%vecj(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%vecj(:,:,:,3,:), Density%vecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%vecs(:,:,:,1,:), Density%vecs(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%vecs(:,:,:,2,:), Density%vecs(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%vecs(:,:,:,3,:), Density%vecs(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rots(:,:,:,1,:), Density%rots(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%rots(:,:,:,2,:), Density%rots(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%rots(:,:,:,3,:), Density%rots(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rotvecj(:,:,:,1,:), Density%rotvecj(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%rotvecj(:,:,:,2,:), Density%rotvecj(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%rotvecj(:,:,:,3,:), Density%rotvecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%ders(:,:,:,1,1,:), Density%ders(:,:,:,1,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,1,2,:), Density%ders(:,:,:,1,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,1,3,:), Density%ders(:,:,:,1,3,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,1,:), Density%ders(:,:,:,2,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,2,:), Density%ders(:,:,:,2,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,3,:), Density%ders(:,:,:,2,3,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,1,:), Density%ders(:,:,:,3,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,2,:), Density%ders(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,3,:), Density%ders(:,:,:,3,3,:),Oct,+1)

        if(allocated(NewDen%LapS)) then
          call FillOct(NewDen%LapS(:,:,:,1,:), Density%LapS(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%LapS(:,:,:,2,:), Density%LapS(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%LapS(:,:,:,3,:), Density%LapS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%graddivS(:,:,:,1,:),Density%gradDivS(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%graddivS(:,:,:,2,:),Density%gradDivS(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%graddivS(:,:,:,3,:),Density%graddivS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%divs(:,:,:,:), Density%divs(:,:,:,:),Oct,+1)
        endif

        if(allocated(NewDen%VecF)) then
          call FillOct(NewDen%vecF(:,:,:,1,:), Density%vecF(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%vecF(:,:,:,2,:), Density%vecF(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%vecF(:,:,:,3,:), Density%vecF(:,:,:,3,:),Oct,+1)
        endif

        if(allocated(NewDen%vecT)) then
          call FillOct(NewDen%vecT(:,:,:,1,:), Density%vecT(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%vecT(:,:,:,2,:), Density%vecT(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%vecT(:,:,:,3,:), Density%vecT(:,:,:,3,:),Oct,+1)
        endif
      endif
    endif
    ! End of Octant 5
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Octant ( x < 0 , y > 0, z < 0)
    ! only necessary when signature and parity are conserved on file, but both
    ! are broken in the calculation
    if((inPC.neqv.PC) .and. (inSC.neqv. SC)) then
      Oct = 6
      call FillOct(NewDen%Rho              , Density%Rho              , Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,1,:), Density%DerRho(:,:,:,1,:), Oct,-1)
      call FillOct(NewDen%DerRho(:,:,:,2,:), Density%DerRho(:,:,:,2,:), Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,3,:), Density%DerRho(:,:,:,3,:), Oct,-1)
      call FillOct(NewDen%LapRho           , Density%LapRho           , Oct,+1)
      call FillOct(NewDen%Tau              , Density%Tau              , Oct,+1)
      call FillOct(NewDen%NablaJ           , Density%NablaJ           , Oct,+1)

      if(allocated(NewDen%JMuNu)) then
        call FillOct(NewDen%JMunu(:,:,:,1,1,:), Density%JMuNu(:,:,:,1,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,1,2,:), Density%JMuNu(:,:,:,1,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,1,3,:), Density%JMuNu(:,:,:,1,3,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,1,:), Density%JMuNu(:,:,:,2,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,2,:), Density%JMuNu(:,:,:,2,2,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,3,:), Density%JMuNu(:,:,:,2,3,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,1,:), Density%JMuNu(:,:,:,3,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,3,2,:), Density%JMuNu(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,3,:), Density%JMuNu(:,:,:,3,3,:),Oct,+1)
      endif

      if(.not.TRC) then
        call FillOct(NewDen%vecj(:,:,:,1,:), Density%vecj(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%vecj(:,:,:,2,:), Density%vecj(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%vecj(:,:,:,3,:), Density%vecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%vecs(:,:,:,1,:), Density%vecs(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%vecs(:,:,:,2,:), Density%vecs(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%vecs(:,:,:,3,:), Density%vecs(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rots(:,:,:,1,:), Density%rots(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%rots(:,:,:,2,:), Density%rots(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%rots(:,:,:,3,:), Density%rots(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rotvecj(:,:,:,1,:), Density%rotvecj(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%rotvecj(:,:,:,2,:), Density%rotvecj(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%rotvecj(:,:,:,3,:), Density%rotvecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%ders(:,:,:,1,1,:), Density%ders(:,:,:,1,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,1,2,:), Density%ders(:,:,:,1,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,1,3,:), Density%ders(:,:,:,1,3,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,1,:), Density%ders(:,:,:,2,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,2,:), Density%ders(:,:,:,2,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,3,:), Density%ders(:,:,:,2,3,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,1,:), Density%ders(:,:,:,3,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,2,:), Density%ders(:,:,:,3,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,3,:), Density%ders(:,:,:,3,3,:),Oct,-1)

        if(allocated(NewDen%LapS)) then
          call FillOct(NewDen%LapS(:,:,:,1,:), Density%LapS(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%LapS(:,:,:,2,:), Density%LapS(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%LapS(:,:,:,3,:), Density%LapS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%graddivS(:,:,:,1,:),Density%gradDivS(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%graddivS(:,:,:,2,:),Density%gradDivS(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%graddivS(:,:,:,3,:),Density%graddivS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%divs(:,:,:,:), Density%divs(:,:,:,:),Oct,-1)
        endif

        if(allocated(NewDen%VecF)) then
          call FillOct(NewDen%vecF(:,:,:,1,:), Density%vecF(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%vecF(:,:,:,2,:), Density%vecF(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%vecF(:,:,:,3,:), Density%vecF(:,:,:,3,:),Oct,+1)
        endif

        if(allocated(NewDen%vecT)) then
          call FillOct(NewDen%vecT(:,:,:,1,:), Density%vecT(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%vecT(:,:,:,2,:), Density%vecT(:,:,:,2,:),Oct,-1)
          call FillOct(NewDen%vecT(:,:,:,3,:), Density%vecT(:,:,:,3,:),Oct,+1)
        endif
      endif
    endif
    !End of Octant 6
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Octant ( x > 0 , y < 0, z < 0)
    ! only necessary when time simplex and parity are conserved on file, but
    ! both are broken in the calculation
    if((inPC.neqv.PC) .and. (inTSC.neqv.TSC)) then
      Oct = 7
     call FillOct(NewDen%Rho              , Density%Rho              , Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,1,:), Density%DerRho(:,:,:,1,:), Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,2,:), Density%DerRho(:,:,:,2,:), Oct,-1)
      call FillOct(NewDen%DerRho(:,:,:,3,:), Density%DerRho(:,:,:,3,:), Oct,-1)
      call FillOct(NewDen%LapRho           , Density%LapRho           , Oct,+1)
      call FillOct(NewDen%Tau              , Density%Tau              , Oct,+1)
      call FillOct(NewDen%NablaJ           , Density%NablaJ           , Oct,+1)

      if(allocated(NewDen%JMuNu)) then
        call FillOct(NewDen%JMunu(:,:,:,1,1,:), Density%JMuNu(:,:,:,1,1,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,1,2,:), Density%JMuNu(:,:,:,1,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,1,3,:), Density%JMuNu(:,:,:,1,3,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,1,:), Density%JMuNu(:,:,:,2,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,2,:), Density%JMuNu(:,:,:,2,2,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,2,3,:), Density%JMuNu(:,:,:,2,3,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,3,1,:), Density%JMuNu(:,:,:,3,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,2,:), Density%JMuNu(:,:,:,3,2,:),Oct,+1)
        call FillOct(NewDen%JMunu(:,:,:,3,3,:), Density%JMuNu(:,:,:,3,3,:),Oct,-1)
      endif

      if(.not.TRC) then
        call FillOct(NewDen%vecj(:,:,:,1,:), Density%vecj(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%vecj(:,:,:,2,:), Density%vecj(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%vecj(:,:,:,3,:), Density%vecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%vecs(:,:,:,1,:), Density%vecs(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%vecs(:,:,:,2,:), Density%vecs(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%vecs(:,:,:,3,:), Density%vecs(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rots(:,:,:,1,:), Density%rots(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%rots(:,:,:,2,:), Density%rots(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%rots(:,:,:,3,:), Density%rots(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rotvecj(:,:,:,1,:), Density%rotvecj(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%rotvecj(:,:,:,2,:), Density%rotvecj(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%rotvecj(:,:,:,3,:), Density%rotvecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%ders(:,:,:,1,1,:), Density%ders(:,:,:,1,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,1,2,:), Density%ders(:,:,:,1,2,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,1,3,:), Density%ders(:,:,:,1,3,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,1,:), Density%ders(:,:,:,2,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,2,2,:), Density%ders(:,:,:,2,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,3,:), Density%ders(:,:,:,2,3,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,1,:), Density%ders(:,:,:,3,1,:),Oct,+1)
        call FillOct(NewDen%ders(:,:,:,3,2,:), Density%ders(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,3,:), Density%ders(:,:,:,3,3,:),Oct,-1)

        if(allocated(NewDen%LapS)) then
          call FillOct(NewDen%LapS(:,:,:,1,:), Density%LapS(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%LapS(:,:,:,2,:), Density%LapS(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%LapS(:,:,:,3,:), Density%LapS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%graddivS(:,:,:,1,:),Density%gradDivS(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%graddivS(:,:,:,2,:),Density%gradDivS(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%graddivS(:,:,:,3,:),Density%graddivS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%divs(:,:,:,:), Density%divs(:,:,:,:),Oct,-1)
        endif

        if(allocated(NewDen%VecF)) then
          call FillOct(NewDen%vecF(:,:,:,1,:), Density%vecF(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%vecF(:,:,:,2,:), Density%vecF(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%vecF(:,:,:,3,:), Density%vecF(:,:,:,3,:),Oct,+1)
        endif

        if(allocated(NewDen%vecT)) then
          call FillOct(NewDen%vecT(:,:,:,1,:), Density%vecT(:,:,:,1,:),Oct,-1)
          call FillOct(NewDen%vecT(:,:,:,2,:), Density%vecT(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%vecT(:,:,:,3,:), Density%vecT(:,:,:,3,:),Oct,+1)
        endif
      endif
    endif
    ! End of octant 7
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Octant ( x < 0 , y < 0, z < 0)
    ! only necessary when all spatial symmetries are conserved on file, but all
    ! of them are broken in the calculation
    if((inPC.neqv.PC) .and. (inSC.neqv. SC)) then
      Oct = 8
      call FillOct(NewDen%Rho              , Density%Rho              , Oct,+1)
      call FillOct(NewDen%DerRho(:,:,:,1,:), Density%DerRho(:,:,:,1,:), Oct,-1)
      call FillOct(NewDen%DerRho(:,:,:,2,:), Density%DerRho(:,:,:,2,:), Oct,-1)
      call FillOct(NewDen%DerRho(:,:,:,3,:), Density%DerRho(:,:,:,3,:), Oct,-1)
      call FillOct(NewDen%LapRho           , Density%LapRho           , Oct,+1)
      call FillOct(NewDen%Tau              , Density%Tau              , Oct,+1)
      call FillOct(NewDen%NablaJ           , Density%NablaJ           , Oct,+1)

      if(allocated(NewDen%JMuNu)) then
        call FillOct(NewDen%JMunu(:,:,:,1,1,:), Density%JMuNu(:,:,:,1,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,1,2,:), Density%JMuNu(:,:,:,1,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,1,3,:), Density%JMuNu(:,:,:,1,3,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,1,:), Density%JMuNu(:,:,:,2,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,2,:), Density%JMuNu(:,:,:,2,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,2,3,:), Density%JMuNu(:,:,:,2,3,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,1,:), Density%JMuNu(:,:,:,3,1,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,2,:), Density%JMuNu(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%JMunu(:,:,:,3,3,:), Density%JMuNu(:,:,:,3,3,:),Oct,-1)
      endif

      if(.not.TRC) then
        call FillOct(NewDen%vecj(:,:,:,1,:), Density%vecj(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%vecj(:,:,:,2,:), Density%vecj(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%vecj(:,:,:,3,:), Density%vecj(:,:,:,3,:),Oct,-1)

        call FillOct(NewDen%vecs(:,:,:,1,:), Density%vecs(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%vecs(:,:,:,2,:), Density%vecs(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%vecs(:,:,:,3,:), Density%vecs(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%rots(:,:,:,1,:), Density%rots(:,:,:,1,:),Oct,-1)
        call FillOct(NewDen%rots(:,:,:,2,:), Density%rots(:,:,:,2,:),Oct,-1)
        call FillOct(NewDen%rots(:,:,:,3,:), Density%rots(:,:,:,3,:),Oct,-1)

        call FillOct(NewDen%rotvecj(:,:,:,1,:), Density%rotvecj(:,:,:,1,:),Oct,+1)
        call FillOct(NewDen%rotvecj(:,:,:,2,:), Density%rotvecj(:,:,:,2,:),Oct,+1)
        call FillOct(NewDen%rotvecj(:,:,:,3,:), Density%rotvecj(:,:,:,3,:),Oct,+1)

        call FillOct(NewDen%ders(:,:,:,1,1,:), Density%ders(:,:,:,1,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,1,2,:), Density%ders(:,:,:,1,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,1,3,:), Density%ders(:,:,:,1,3,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,1,:), Density%ders(:,:,:,2,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,2,:), Density%ders(:,:,:,2,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,2,3,:), Density%ders(:,:,:,2,3,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,1,:), Density%ders(:,:,:,3,1,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,2,:), Density%ders(:,:,:,3,2,:),Oct,-1)
        call FillOct(NewDen%ders(:,:,:,3,3,:), Density%ders(:,:,:,3,3,:),Oct,-1)

        if(allocated(NewDen%LapS)) then
          call FillOct(NewDen%LapS(:,:,:,1,:), Density%LapS(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%LapS(:,:,:,2,:), Density%LapS(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%LapS(:,:,:,3,:), Density%LapS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%graddivS(:,:,:,1,:),Density%gradDivS(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%graddivS(:,:,:,2,:),Density%gradDivS(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%graddivS(:,:,:,3,:),Density%graddivS(:,:,:,3,:),Oct,+1)

          call FillOct(NewDen%divs(:,:,:,:), Density%divs(:,:,:,:),Oct,-1)
        endif

        if(allocated(NewDen%VecF)) then
          call FillOct(NewDen%vecF(:,:,:,1,:), Density%vecF(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%vecF(:,:,:,2,:), Density%vecF(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%vecF(:,:,:,3,:), Density%vecF(:,:,:,3,:),Oct,+1)
        endif

        if(allocated(NewDen%vecT)) then
          call FillOct(NewDen%vecT(:,:,:,1,:), Density%vecT(:,:,:,1,:),Oct,+1)
          call FillOct(NewDen%vecT(:,:,:,2,:), Density%vecT(:,:,:,2,:),Oct,+1)
          call FillOct(NewDen%vecT(:,:,:,3,:), Density%vecT(:,:,:,3,:),Oct,+1)
        endif
      endif
    endif
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Save the new density
    Density =  NewDen
  end subroutine TransformDensities

  subroutine TransformHFBMatrices(inU, inV, inRho, InKappa,inPC,inIC,          &
    &                                                InColumns, InputBlocksizes)
    !---------------------------------------------------------------------------
    ! Subroutine to transform all HFB matrices.
    !---------------------------------------------------------------------------

    use HFB

    complex(KIND=dp),intent(in)   :: InKappa(:,:,:,:), InRho(:,:,:,:)
    complex(KIND=dp),intent(in)   :: InU(:,:,:,:), InV(:,:,:,:)
    integer, intent(in)           :: InputBlockSizes(2,2), InColumns(:,:,:)
    logical, intent(in)           :: inPC, inIC
    integer                       :: sizes(2),it,i,j, check

    if(.not. inIC) then
        call stp('No rules to transform the HFB matrices yet when breaking Isospin.')
    endif

    if(PC .eqv. inPC) then
        ! File matches rundata ; no action necessary
        KappaHFB = InKappa
        RhoHFB   = InRho
        U        = inU
        V        = inV
        HFBColumns = inColumns
    else
        !---------------------------------------------------------------------
        ! On file, we have for Kappa:
        !
        !          (  0    K^+ )            (  0    K^- )
        !          (           )    and     (           )
        !          ( -K^+  0   )            ( -K^-  0   )
        !
        ! And for Rho:
        !
        !          (  R^+_+  0    )           (  R^-_+   0   )
        !          (              )   and     (              )
        !          (  0    R^+_-  )           (  0     R^-_- )
        !
        ! And for U
        !          ( U^+_+     0    )           (  U^-_+     0    )
        !          (                )   and     (                 )
        !          ( 0        U^+_- )           (  0         U^-_+)
        !
        ! And for V
        !          ( 0        V^+_-   )         (  0       V^-_+  )
        !          (                  )   and   (                 )
        !          ( V^+_+     0      )         (  V^-_+   U^-_+  )
        !
        !
        !---------------------------------------------------------------------
        !
        ! The correct forms are
        ! (using the order of wavefunctions as generated by the NIL8 code)
        !
        !                ( 0    0    K^+ 0  )
        !  Kappa =       ( 0    0    0   K^-)
        !                ( -K^+ 0    0   0  )
        !                ( 0   -K^-  0   0  )
        !
        !  Rho   =       ( R^+_+ 0      0     0    )
        !                ( 0     R^-_+  0     0    )
        !                ( 0     0      R^+_- 0    )
        !                ( 0     0      0     R^-_-)
        !
        !----------------------------------------------------------------------
        if(.not.allocated(HFBColumns)) then
          allocate(HFBColumns(HFBSize,1,2)) ; HFBColumns = 0
        endif

        if(.not. SC) then
          print *, 'WARNING: TRANSFORMATION RULES FOR U AND V MIGHT NOT BE '   &
          &      //'CORRECT WHEN BREAKING PARITY FROM A SIGNATURE BROKEN '     &
          &      //'CALCULATION.'
        endif

        do it=1,2
            sizes = InputBlocksizes(:,it)
            !-------------------------------------------------------------------
            ! Anomalous density matrix

            ! Parity plus, meaning P=2 for inputKappa and input Rho
            KappaHFB(sum(sizes)/2+1:sum(sizes)/2+sizes(2)/2,1:sizes(2)/2,1,it) &
            &                 = InKappa(sizes(2)/2+1:sizes(2),1:sizes(2)/2,2,it)
            ! Parity minus, meaning P=1 for inputKappa
            KappaHFB(sum(sizes)/2+sizes(2)/2+1:sum(sizes),                     &
            &                     sizes(2)/2+1:sum(sizes)/2,1,it )             &
            &                 = InKappa(sizes(1)/2+1:sizes(1),1:sizes(1)/2,1,it)

            ! Antisymmetrize Kappa
            do i=1,HFBSize
              do j=i+1, HFBSize
                KappaHFB(i,j,1,it) = - KappaHFB(j,i,1,it)
              enddo
            enddo

            !-------------------------------------------------------------------
            ! Density matrix
            ! Positive parity, positive signature
            RhoHFB(1:sizes(2)/2,1:sizes(2)/2,1,it)                             &
            &                            = InRho(1:sizes(2)/2,1:sizes(2)/2,2,it)
            ! Negative parity, positive signature
            RhoHFB(sizes(2)/2+1:sum(sizes)/2,sizes(2)/2+1:sum(sizes)/2,1,it)   &
            &                            = InRho(1:sizes(1)/2,1:sizes(1)/2,1,it)

            ! Positive parity, negative signature
            RhoHFB(sum(sizes)/2+1:sum(sizes)/2+sizes(1)/2,                     &
            &      sum(sizes)/2+1:sum(sizes)/2+sizes(1)/2,1,it)                &
            &          = InRho(sizes(2)/2+1:sizes(2),sizes(2)/2+1:sizes(2),2,it)
            ! Negative parity, negative signature
            RhoHFB(sum(sizes)/2+sizes(2)/2+1:sum(sizes),                       &
            &      sum(sizes)/2+sizes(2)/2+1:sum(sizes),1,it)                  &
            &          = InRho(sizes(1)/2+1:sizes(1),sizes(1)/2+1:sizes(1),1,it)

            !-------------------------------------------------------------------
            ! U and V matrices
            ! The order of these is not very important, as long as the
            ! HFBColumns matrix is correctly initialized.
            U(:,:,:,it) = 0.0d0 ; V(:,:,:,it) = 0.0d0

            ! Positive parity, positive signature
            U(1:sizes(2)/2, 1:sizes(2),1,it) =                                 &
            &                               inU(1:sizes(2)/2, 1:sizes(2),2,it)

            ! Negative parity, positive signature
            U(sizes(2)/2+1:sum(sizes)/2, sizes(2)+1:sum(sizes),1,it) =         &
            &                            inU( 1:sizes(1)/2,   1:sizes(1),1,it)

            ! Positive parity, negative signature
            U(sum(sizes)/2+1:sizes(2)+sizes(1)/2,                              &
            & sum(sizes)+1:sum(sizes) + sizes(2),1,it) =                       &
            &             inU(sizes(2)/2+1:sizes(2), sizes(2)+1:2*sizes(2),2,it)

            ! Negative parity, negative signature
            U(  sizes(2)+sizes(1)/2+1:  sum(sizes),                            &
            & 2*sizes(2)+sizes(1)  +1:2*sum(sizes),1,it)=                      &
            &             inU(sizes(1)/2+1:sizes(1), sizes(1)+1:2*sizes(1),1,it)
            !

            V(sum(sizes)/2 + 1:sum(sizes)/2 + sizes(2)/2, 1:sizes(2),1,it) =   &
            &                        inV(sizes(2)/2+1:sizes(2), 1:sizes(2),2,it)

            V(sizes(2)+sizes(1)/2+1:sum(sizes), sizes(2)+1:sum(sizes),1,it) =  &
            &                        inV(sizes(1)/2+1:sizes(1), 1:sizes(1),1,it)

            V(1:sizes(2)/2,                                                    &
            & sum(sizes) +1: sum(sizes) + sizes(2),1,it) =                     &
            &              inV(1:sizes(2),sizes(2)+1:2*sizes(2),2,it)

            V(sizes(2)/2+1:sum(sizes)/2,                                       &
            & sum(sizes) + sizes(2)+1: 2*sum(sizes),1,it) =                    &
            &              inV(1:sizes(1),sizes(1)+1:2*sizes(1),1,it)
            !

            do i=1,sizes(2)
              if(Incolumns(i,2,it) .gt. sizes(2)) then
                HFBColumns(i,1,it) = Incolumns(i,2,it) + sizes(1)
              else
                HFBColumns(i,1,it) = Incolumns(i,2,it)
              endif
            enddo

            do i=1,sizes(1)
              if(Incolumns(i,1,it) .gt. sizes(1)) then
                HFBColumns(i+sizes(2),1,it) = Incolumns(i,1,it) + 2*sizes(2)
              else
                HFBColumns(i+sizes(2),1,it) = Incolumns(i,1,it) + sizes(2)
              endif
            enddo

        enddo
    endif
end subroutine TransformHFBMatrices

!   function TransformKappa( InKappa, inPC, inIC, InputBlocksizes) result(OutKappa)
!     !----------------------------------------------------------------------------
!     ! Function that takes Kappa from a file and checks if it needs transformation
!     ! when either Parity or Isospin gets broken.
!     !
!     !----------------------------------------------------------------------------
!     use HFB
!     use Spwfstorage

!     complex(KIND=dp)              :: InKappa(:,:,:,:)
!     complex(KIND=dp), allocatable :: OutKappa(:,:,:,:), Temp(:)
!     logical, intent(in)           :: inPC, inIC
!     integer                       :: sizes(2), P, it, iplus,imin,i,j

!     integer, intent(in)           :: InputBlockSizes(2,2)

!     if(.not. inIC) then
!         call stp('No rules to transform Kappa yet when breaking Isospin.')
!     endif

!     if(PC .eqv. inPC) then
!         ! File matches rundata ; no action necessary
!         OutKappa = InKappa
!     else
!         !-----------------------------------------------------------------------
!         ! Parity is conserved on file, broken in run.
!         ! On file, we have two matrices:
!         !
!         !          (  0    K^+ )           (  0    K^- )
!         !          (           )   and     (           )
!         !          ( -K^+  0   )           ( -K^-  0   )
!         !
!         ! HOWEVER, taking as the following as NEW kappa is NOT ALLOWED!
!         !
!         !          (  0    K^+   0    0   )
!         !          ( -K^+   0    0    0   )
!         !          (  0     0    0    K^- )
!         !          (  0     0   -K^-  0   )
!         !
!         ! because the wavefunction indices in the HFB module need not
!         ! be in the correct order.
!         !
!         ! The correct form is (using the order of wavefunctions as generated
!         ! by the NIL8 code)
!         !         ( 0    0    K^+ 0  )
!         !         ( 0    0    0   K^-)
!         !         ( -K^+ 0    0   0  )
!         !         ( 0   -K^-  0   0  )
!         !
!         !-----------------------------------------------------------------------

!         ! We can take this size anyway, since P needs to be broken
!         ! and isospin conserved.
!         allocate(OutKappa(HFBSize,HFBSize,1,2)) ; OutKappa = 0.0_dp
!         do it=1,2
!             sizes = InputBlocksizes(:,it)

!             OutKappa(sum(sizes)/2+1:sum(sizes)/2+sizes(2)/2,1:sizes(2)/2,1,it)            &
!             & = InKappa(sizes(2)/2+1:sizes(2),1:sizes(2)/2,2,it)

!             OutKappa(sum(sizes)/2+sizes(2)/2+1:sum(sizes),sizes(2)/2+1:sum(sizes)/2,1,it )&
!             & = InKappa(sizes(1)/2+1:sizes(1),1:sizes(1)/2,1,it)

!             do i=1,HFBSize
!               do j=i+1, HFBSize
!                 OutKappa(i,j,1,it) = - OutKappa(j,i,1,it)
!               enddo
!             enddo
!         enddo
!     endif

!   end function TransformKappa

  subroutine Switch( A, index1, index2 )
    !--------------------------------------------------------------------------
    ! Small subroutine that switches indices on a matrix.
    !--------------------------------------------------------------------------
    complex(KIND=dp), intent(inout) :: A(:,:)
    complex(KIND=dp), allocatable   :: Temp(:)
    integer, intent(in)             :: index1,index2
    integer                         :: N

    N = size(A,1)
    allocate(Temp(N)) ; Temp = 0.0d0
    if(index1 .gt. N .or. index2 .gt. N) stop

    !First switch rows
    Temp = A(index1,:)
    A(index1,:) = A(index2,:)
    A(index2,:) = Temp

    !Then switch columns
    Temp = A(:,index1)
    A(:,index1) = A(:,index2)
    A(:,index2) = Temp

  end subroutine Switch

  subroutine Displace(Direction, wf, Coord)
    !---------------------------------------------------------------------------
    ! Wrapper routine for displacing the wavefunction in coordinate space.
    !
    !---------------------------------------------------------------------------
    type(Spwf)          :: wf
    integer, intent(in) :: Direction, Coord

    select case(Direction)

    case(1)
      call DisplaceX(wf,coord)
    case(2)
      call DisplaceY(wf,coord)
    case(3)
      call DisplaceZ(wf,coord)
    case DEFAULT
      call stp('Non-cartesian direction in Displace!')
    end select

  end subroutine Displace

  subroutine DisplaceX(wf,coord)
    !---------------------------------------------------------------------------
    ! Displace the wavefunction in the X-direction with 'coord' points.
    !---------------------------------------------------------------------------
    use Spinors

    type(Spwf)          :: wf
    integer, intent(in) :: coord
    integer             :: startx,endx,i
    type(Spinor)        :: Displaced, Temp

    if(Coord.eq.0) return
    !---------------------------------------------------------------------------
    ! Check for signature
    if(SC) then
      call stp('DisplaceX should not be called when signature is conserved!')
    endif
    !---------------------------------------------------------------------------
    ! Getting the loop indices
    if(Coord.lt.0) then
      startx =    - Coord + 1
      endx   = nx
    else
      startx = 1
      endx   = nx - Coord
    endif
    Temp = wf%GetValue()
    !---------------------------------------------------------------------------
    ! Displacing the complete wavefunction
    Displaced= NewSpinor()
    do i=startx,endx
      Displaced%Grid(i+Coord,:,:,:,:) =Temp%Grid(i,:,:,:,:)
    enddo
    call wf%SetGrid(Displaced)
  end subroutine DisplaceX

  subroutine DisplaceY(wf,coord)
    !---------------------------------------------------------------------------
    ! Displace the wavefunction in the y-direction with 'coord' points.
    !---------------------------------------------------------------------------
    use Spinors

    type(Spwf)          :: wf
    integer, intent(in) :: coord
    integer             :: starty,endy,j
    type(Spinor)        :: Displaced, Temp

    if(Coord.eq.0) return
    !---------------------------------------------------------------------------
    ! Check for TimeSimplex
    if(TSC) then
      call stp('DisplaceY should not be called when TimeSimplex is conserved!')
    endif
    !---------------------------------------------------------------------------
    ! Getting the loop indices
    if(Coord.lt.0) then
      starty =    - Coord + 1
      endy   = ny
    else
      starty = 1
      endy   = ny - Coord
    endif
    Temp = wf%GetValue()
    !---------------------------------------------------------------------------
    ! Displacing the complete wavefunction
    Displaced= NewSpinor()
    do j=starty,endy
      Displaced%Grid(:,j+Coord,:,:,:) =Temp%Grid(:,j,:,:,:)
    enddo
    call wf%SetGrid(Displaced)
  end subroutine DisplaceY

  subroutine DisplaceZ(wf,coord)
    !---------------------------------------------------------------------------
    ! Displace the wavefunction in the z-direction with 'coord' points.
    !---------------------------------------------------------------------------
    use SPinors
    type(Spwf)          :: wf
    integer, intent(in) :: coord
    integer             :: startz,endz,k
    type(Spinor)        :: Displaced, Temp

    if(Coord.eq.0) return
    !---------------------------------------------------------------------------
    ! Check for TimeSimplex
    if(PC) then
      call stp('DisplaceZ should not be called when Parity is conserved!')
    endif
    !---------------------------------------------------------------------------
    ! Getting the loop indices
    if(Coord.lt.0) then
      startz =    - Coord + 1
      endz   = nz
    else
      startz = 1
      endz   = nz - Coord
    endif
    Temp = wf%GetValue()
    !---------------------------------------------------------------------------
    ! Displacing the complete wavefunction
    Displaced= NewSpinor()
    do k=startz,endz
      Displaced%Grid(:,:,k+Coord,:,:) =Temp%Grid(:,:,k,:,:)
    enddo
    call wf%SetGrid(Displaced)
  end subroutine DisplaceZ

  subroutine SymDisplace(Direction,inwf, outwf1,outwf2,Coord)
    !---------------------------------------------------------------------------
    ! Wrapper subroutine for symmetrically displacing the input wavefunction.
    !
    ! The output is are the wavefunctions
    !
    ! sqrt(2) * outwf1 = [T_(x/y/z) (a) + T_(x/y/z) (-a)] inwf
    ! sqrt(2) * outwf2 = [T_(x/y/z) (a) - T_(x/y/z) (-a)] inwf
    !
    ! Note that outwf1 will in general share quantum numbers with inwf,
    ! but outwf2 will not!
    !---------------------------------------------------------------------------
    type(Spwf), intent(in) :: inwf
    type(Spwf), intent(out):: outwf1, outwf2
    integer, intent(in)    :: Coord, Direction

    select case(Direction)
    case(1)
      call SymDisplaceX(inwf,outwf1,outwf2,Coord)
    case(2)
      call SymDisplaceY(inwf,outwf1,outwf2,Coord)
    case(3)
      call SymDisplaceZ(inwf,outwf1,outwf2,Coord)
    end select

  end subroutine SymDisplace

  subroutine SymDisplaceX(inwf,outwf1,outwf2, Coord)
    !---------------------------------------------------------------------------
    ! Symmetrically displace inwf over Coord to outwf1 & outwf2.
    ! The output is
    !
    ! sqrt(2) * outwf1 = [T_(x) (a) + T_(x) (-a)] inwf
    ! sqrt(2) * outwf2 = [T_(x) (a) - T_(x) (-a)] inwf
    !
    ! Note that outwf will have opposite signature and parity with respect to
    ! outwf2.
    ! In addition, when time reversal is conserved, we don't want any negative
    ! signature states, so then outwf2 gets acted upon with the time-reversal
    ! operator.
    !---------------------------------------------------------------------------
    use Spinors

    type(Spwf), intent(in) :: inwf
    type(Spwf), intent(out):: outwf1, outwf2
    type(Spwf)             :: workwf
    integer, intent(in)    :: Coord
    integer                :: startx(2),endx(2), q, i,j,k
    type(Spinor)           :: Displaced(2), Symmetric(2), Temp, Final(2)

    if(Coord.eq.0) call stp("Don't displace by zero!")
    !---------------------------------------------------------------------------
    ! Make an initial copy
    outwf1 = inwf
    outwf2 = inwf
    !---------------------------------------------------------------------------
    ! Explicitly breaking signature if necessary
    if(SC) then
      nx = nx * 2
      workwf = inwf
      call BreakSignature(Workwf)
      temp = workwf%GetValue()
    else
      temp = inwf%GetValue()
    endif
    !---------------------------------------------------------------------------
    ! Getting the iteration indices
    startx(1)    = 1          ; startx(2)    = 1 + Coord
    endx(1)      = nx - Coord ; endx(2)      = nx
    Displaced(1) = NewSpinor()   ; Displaced(2) = NewSpinor()
    !---------------------------------------------------------------------------
    ! Applying T_x(-a) and T_x(+a)
    do i=startx(1),endx(1)
      Displaced(1)%Grid(i+Coord,:,:,:,:) =  Temp%Grid(i,:,:,:,:)
    enddo
    do i=startx(2), endx(2)
      Displaced(2)%Grid(i-Coord,:,:,:,:) =  Temp%Grid(i,:,:,:,:)
    enddo
    !---------------------------------------------------------------------------
    ! Taking symmetric combinations
    Symmetric(1) = 1.0/sqrt(2.0_dp)*(Displaced(1) + Displaced(2))
    Symmetric(2) = 1.0/sqrt(2.0_dp)*(Displaced(1) - Displaced(2))
    !---------------------------------------------------------------------------
    !Reverting the damage to nx
    if(SC) then
      nx = nx/2
      do q=1,2
        Final(q) = NewSpinor()
        do i=1, nx
          Final(q)%Grid(i,:,:,:,:) = Symmetric(q)%Grid(i+nx,:,:,:,:)
        enddo
      enddo
    else
      do q=1,2
        Final(q) = Symmetric(q)
      enddo
    endif
    !---------------------------------------------------------------------------
    ! Setting the result to the output wavefunctions, and change parity and
    ! signature of the second wavefunction.
    call outwf1%SetGrid(Final(1))
    if(TSC) then
      call outwf2%SetGrid(TimeReverse(Final(2)))
    else
      call outwf2%SetGrid(Final(2))
      call outwf2%SetSignature(-inwf%GetSignature())
    endif
    call outwf2%SetParity(- inwf%GetParity())
  end subroutine SymDisplaceX

  subroutine SymDisplaceY(inwf,outwf1,outwf2, Coord)
    !---------------------------------------------------------------------------
    ! Symmetrically displace inwf over Coord to outwf1 & outwf2.
    ! The output is
    !
    ! sqrt(2) * outwf1 = [T_(y) (a) + T_(y) (-a)] inwf
    ! sqrt(2) * outwf2 = [T_(y) (a) - T_(y) (-a)] inwf
    !
    ! Note that outwf2 will have opposite signature,time simplex and parity with
    ! respect to inwf.
    ! In addition, when time reversal is conserved, we don't want any negative
    ! signature states, so then outwf2 gets acted upon with the time-reversal
    ! operator.
    !---------------------------------------------------------------------------
    use Spinors
    type(Spwf), intent(in) :: inwf
    type(Spwf), intent(out):: outwf1, outwf2
    type(Spwf)             :: workwf
    integer, intent(in)    :: Coord
    integer                :: starty(2),endy(2), q, i,j,k
    type(Spinor)           :: Displaced(2), Symmetric(2), Temp, Final(2)

    if(Coord.eq.0) call stp("Don't displace by zero!")
    !---------------------------------------------------------------------------
    ! Make an initial copy
    outwf1 = inwf
    outwf2 = inwf
    !---------------------------------------------------------------------------
    ! Explicitly breaking signature if necessary
    if(TSC) then
      ny = ny * 2
      workwf = inwf
      call BreakTimeSimplex(Workwf)
      temp = workwf%GetValue()
    else
      temp = inwf%GetValue()
    endif
    !---------------------------------------------------------------------------
    ! Getting the iteration indices
    starty(1)    = 1          ; starty(2)    = 1 + Coord
    endy(1)      = ny - Coord ; endy(2)      = ny
    Displaced(1) = NewSpinor()   ; Displaced(2) = NewSpinor()
    !---------------------------------------------------------------------------
    ! Applying T_x(-a) and T_x(+a)
    do j=starty(1),endy(1)
      Displaced(1)%Grid(:,j+Coord,:,:,:) =  Temp%Grid(:,j,:,:,:)
    enddo
    do j=starty(2), endy(2)
      Displaced(2)%Grid(:,j-Coord,:,:,:) =  Temp%Grid(:,j,:,:,:)
    enddo
    !---------------------------------------------------------------------------
    ! Taking symmetric combinations
    Symmetric(1) = 1.0/sqrt(2.0_dp)*(Displaced(1) + Displaced(2))
    Symmetric(2) = 1.0/sqrt(2.0_dp)*(Displaced(1) - Displaced(2))
    !---------------------------------------------------------------------------
    !Reverting the damage to nx
    if(SC) then
      ny = ny/2
      do q=1,2
        Final(q) = NewSpinor()
        do j=1, ny
          Final(q)%Grid(:,j,:,:,:) = Symmetric(q)%Grid(:,j+ny,:,:,:)
        enddo
      enddo
    else
      do q=1,2
        Final(q) = Symmetric(q)
      enddo
    endif
    !---------------------------------------------------------------------------
    ! Setting the result to the output wavefunctions, and change parity,
    ! signature and time simplex of the second wavefunction.
    call outwf1%SetGrid(Final(1))
    if(TSC) then
      Final(2) = TimeReverse(Final(2))
    else
      call outwf2%SetSignature(-inwf%GetSignature())
    endif
    call outwf2%SetParity(- inwf%GetParity())
    if(TSC) then
      Final(2) = - Final(2)
    endif
    call outwf2%SetGrid(Final(2))
  end subroutine SymDisplaceY

  subroutine SymDisplaceZ(inwf,outwf1,outwf2, Coord)
    !---------------------------------------------------------------------------
    ! Symmetrically displace inwf over Coord to outwf1 & outwf2.
    ! The output is
    !
    ! sqrt(2) * outwf1 = [T_(z) (a) + T_(z) (-a)] inwf
    ! sqrt(2) * outwf2 = [T_(z) (a) - T_(z) (-a)] inwf
    !
    ! Note that outwf2 will have opposite parity with respect to inwf.
    !---------------------------------------------------------------------------
    use Spinors
    type(Spwf)             :: inwf
    type(Spwf), intent(out):: outwf1, outwf2
    type(Spwf)             :: workwf
    integer, intent(in)    :: Coord
    integer                :: startz(2),endz(2), q, i,j,k
    type(Spinor)           :: Displaced(2), Symmetric(2), Temp, Final(2)

    if(Coord.eq.0) call stp("Don't displace by zero!")
    !---------------------------------------------------------------------------
    ! Make an initial copy
    outwf1 = inwf
    outwf2 = inwf
    !---------------------------------------------------------------------------
    ! Explicitly breaking signature if necessary
    if(PC) then
      nz = nz * 2
      workwf = inwf
      call BreakParity(Workwf)
      temp = workwf%GetValue()
    else
      temp = inwf%GetValue()
    endif
    !---------------------------------------------------------------------------
    ! Getting the iteration indices
    startz(1)    = 1          ; startz(2)    = 1 + Coord
    endz(1)      = nz - Coord ; endz(2)      = nz
    Displaced(1) = NewSpinor()   ; Displaced(2) = NewSpinor()
    !---------------------------------------------------------------------------
    ! Applying T_x(-a) and T_x(+a)
    do k=startz(1),endz(1)
      Displaced(1)%Grid(:,:,k+Coord,:,:) =  Temp%Grid(:,:,k,:,:)
    enddo
    do k=startz(2), endz(2)
      Displaced(2)%Grid(:,:,k-Coord,:,:) =  Temp%Grid(:,:,k,:,:)
    enddo
    !---------------------------------------------------------------------------
    ! Taking symmetric combinations
    Symmetric(1) = 1.0/sqrt(2.0_dp)*(Displaced(1) + Displaced(2))
    Symmetric(2) = 1.0/sqrt(2.0_dp)*(Displaced(1) - Displaced(2))
    !---------------------------------------------------------------------------
    !Reverting the damage to nx
    if(PC) then
      nz = nz/2
      do q=1,2
        Final(q) = NewSpinor()
        do k=1, nz
          Final(q)%Grid(:,:,k,:,:) = Symmetric(q)%Grid(:,:,k+nz,:,:)
        enddo
      enddo
    else
      do q=1,2
        Final(q) = Symmetric(q)
      enddo
    endif
    !---------------------------------------------------------------------------
    ! Setting the result to the output wavefunctions, and change parity of the
    ! second wavefunction.
    call outwf2%SetParity(- inwf%GetParity())
    call outwf2%SetGrid(Final(2))
    call outwf1%SetGrid(Final(1))
  end subroutine SymDisplaceZ

  subroutine BreakTimeReversal(wf, NewWfOne,NewWfTwo)
  !-----------------------------------------------------------------------------
  ! This subroutine takes as input a wavefunction wf with conserved time
  ! reversal symmetry. Two new wavefunctions NewWfOne and newWFTwo get created
  ! with opposite signatures.
  !-----------------------------------------------------------------------------

    use Spinors

    type(Spwf), intent(in) :: wf
    type(Spwf), intent(out):: NewWFOne, NewWFTwo
    type(Spinor)           :: Temp

    !Checking input
    if(wf%GetTimeReversal().eq.0) call stp("Time Reversal is already broken!")

    !Copying the wavefunctions
    NewWfOne = CopyWaveFunction(wf)
    NewWfTwo = CopyWaveFunction(wf)
    !Changing time reversal quantum numbers
    call NewWfOne%SetTimeReversal(0)
    call NewWfTwo%SetTimeReversal(0)

    !Halving the occupation
    call NewWfOne%SetOcc(NewWfOne%GetOcc()/2.0_dp)
    call NewWfTwo%SetOcc(NewWfTwo%GetOcc()/2.0_dp)

    !Changing signature quantum number
    call NewWfOne%SetSignature( wf%signature)
    call NewWfTwo%SetSignature(-wf%signature)

    !Switching the components of nr.2 around
    ! Action of timereversal is -i *\sigma_y * K
    Temp = wf%GetValue()
    Temp = TimeReverse(Temp)

    call NewWfTwo%SetGrid(Temp)
  end subroutine BreakTimeReversal
  
  subroutine ConstructInterpolationFunctions(mx,my,mz,ex)
    !---------------------------------------------------------------------------
    ! Construct the values of the interpolation functions f_r associated with
    ! the old mesh at the points of the new mesh.
    ! (mx,my,mz, ex) => (nx,ny,nz,dx)
    !
    !
    ! Currently only valid when conserving all of the spatial symmetries
    ! 
    !---------------------------------------------------------------------------
    
    integer, intent(in)       :: mx,my,mz
    real(KIND=dp), intent(in) :: ex
    integer                   :: i,j,k
    real(KIND=dp)             :: fac, c,d,x1,x2, dh, ph
    real(KIND=dp), parameter  :: eps = 1d-3
    
    allocate(interpolX(nx,mx,2), interpolY(ny,my,2), interpolZ(nz,mz,2))
    InterpolX = 0.0d0
    InterpolY = 0.0d0
    InterpolZ = 0.0d0
    
    if(.not. SC) then
      call stp('Interpolation not yet supported for signature breaking calculations.')
    endif
    if(.not. TSC) then
      call stp('Interpolation not yet supported for timesimplex breaking calculations.')
    endif
    
    dh = dx/ex
    fac =   0.5d0 / mx
    x1  = - 0.5d0 * dh
    ph  = 0.5 * pi/mx
    do i=1,nx
      x1 = x1 + dh
      x2  = - 0.5d0
      do j=1,mx        
        x2 = x2 + 1
        if (abs(x1-x2).le.eps) then
          c = pi/ph
        else
          c = sin(pi * (x1 - x2))/sin(ph*(x1-x2))
        endif
        if (abs((x1+x2)/(2*mx)-1.0_dp).le.eps) then
          d= pi/ph
        else
          d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
        endif
        InterpolX(i,j,1) = fac * (c - d)
        InterpolX(i,j,2) = fac * (c + d)
      enddo
    enddo
    
    fac =   0.5d0 / my
    x1  = - 0.5d0 * dh
    ph  = 0.5 * pi/my
    do i=1,ny
      x1 = x1 + dh
      x2  = - 0.5d0
      do j=1,my        
        x2 = x2 + 1
        if (abs(x1-x2).le.eps) then
          c = pi/ph
        else
          c = sin(pi * (x1 - x2))/sin(ph*(x1-x2))
        endif
        if (abs((x1+x2)/(2*my)-1.0_dp).le.eps) then
          d= pi/ph
        else
          d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
        endif
        InterpolY(i,j,1) = fac * (c - d)
        InterpolY(i,j,2) = fac * (c + d)
      enddo
    enddo
    
    fac =   0.5d0 / mz
    x1  = - 0.5d0 * dh
    if(.not. PC) x1 = x1 - nz/2 * dh   
    ph  =   0.5 * pi/mz
    do i=1,nz
      x1 = x1 + dh
      x2 = - 0.5d0
      if(.not. PC) x2 = x2 - mz/2 
      do j=1,mz        
        x2 = x2 + 1
        if (abs(x1-x2).le.eps) then
          c = pi/ph
        else
          c = sin(pi*(x1-x2))/sin(ph*(x1-x2))
        endif        
        if(.not. PC) then 
          d = 0.0d0
        else
          if (abs((x1+x2)/(2*mz)-1.0_dp).le.eps) then
            d= pi/ph
          else
            d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
          endif
        endif
        
        InterpolZ(i,j,1) = fac * (c - d)
        InterpolZ(i,j,2) = fac * (c + d)
      enddo
    enddo
    
  end subroutine ConstructInterpolationFunctions
    
  function InterpolateSpwf(Phi) result(Psi)
    !---------------------------------------------------------------------------
    ! Interpolate a given wave-function from a certain mesh to 
    ! the new mesh (nx,ny,nz,dx), using the interpolation functions.
    !
    !---------------------------------------------------------------------------
    
    type(Spwf), intent(in) :: Phi
    type(Spwf)             :: Psi
    integer                :: oldnx,oldny,oldnz
    integer                :: i,j,k, ii,jj,kk, P,l, aP
    real(KIND=dp)          :: norm(2)
    real(KIND=dp),allocatable :: temp(:,:,:,:), alt(:,:,:,:)
    
    Psi = CopyWaveFunction(Phi)
    
    Psi%Value = NewSpinor()
    Psi%Der(1)= NewSpinor()
    Psi%Der(2)= NewSpinor()
    Psi%Der(3)= NewSpinor()
    Psi%Lap   = NewSpinor()
    
    oldnx = size(Phi%Value%Grid,1)
    oldny = size(Phi%Value%Grid,2)
    oldnz = size(Phi%Value%Grid,3)
    
    ! First interpolate the x-direction
    allocate(temp(nx,oldny,oldnz,4))
    allocate(alt (nx,ny,oldnz,4)) 
    do k=1,oldnz
      do j=1,oldny
        do i=1,nx 
          temp(i,j,k,:) = 0.0_dp
          do ii=1,oldnx
            temp(i,j,k,1) = temp(i,j,k,1) + InterpolX(i,ii,2) * Phi%Value%Grid(ii,j,k,1,1)
            temp(i,j,k,2) = temp(i,j,k,2) + InterpolX(i,ii,1) * Phi%Value%Grid(ii,j,k,2,1)
            temp(i,j,k,3) = temp(i,j,k,3) + InterpolX(i,ii,1) * Phi%Value%Grid(ii,j,k,3,1)
            temp(i,j,k,4) = temp(i,j,k,4) + InterpolX(i,ii,2) * Phi%Value%Grid(ii,j,k,4,1)
          enddo
        enddo
      enddo
    enddo
    ! Then interpolate in the y-direction
    do k=1,oldnz
      do i=1,nx
        do j=1,ny 
          alt(i,j,k,:) = 0.0_dp
          do jj=1,oldny
            alt(i,j,k,1) = alt(i,j,k,1) + InterpolY(j,jj,2) * temp(i,jj,k,1)
            alt(i,j,k,2) = alt(i,j,k,2) + InterpolY(j,jj,1) * temp(i,jj,k,2)
            alt(i,j,k,3) = alt(i,j,k,3) + InterpolY(j,jj,2) * temp(i,jj,k,3)
            alt(i,j,k,4) = alt(i,j,k,4) + InterpolY(j,jj,1) * temp(i,jj,k,4)
          enddo
        enddo
      enddo
    enddo
    
    !---------------------------------------------------------------------------
    ! Finally interpolate in the z-direction
    ! Depending on the parity of the wave-function, we need to use the 
    ! interpolation coefficients differently.
    P  = (Phi%Parity + 3)/2 
    if(PC) then
      aP = mod(P,2)+1
    else
      aP = P
    endif
    do j=1,ny
      do i=1,nx
        do k=1,nz
          Psi%Value%Grid(i,j,k,:,1) = 0.0_dp
          do kk=1,oldnz
            Psi%Value%Grid(i,j,k,1,1) = Psi%Value%Grid(i,j,k,1,1) + InterpolZ(k,kk,P ) * alt(i,j,kk,1)
            Psi%Value%Grid(i,j,k,2,1) = Psi%Value%Grid(i,j,k,2,1) + InterpolZ(k,kk,P ) * alt(i,j,kk,2)
            Psi%Value%Grid(i,j,k,3,1) = Psi%Value%Grid(i,j,k,3,1) + InterpolZ(k,kk,aP) * alt(i,j,kk,3)
            Psi%Value%Grid(i,j,k,4,1) = Psi%Value%Grid(i,j,k,4,1) + InterpolZ(k,kk,aP) * alt(i,j,kk,4)
          enddo
        enddo
      enddo
    enddo
    
!    do i=1,oldnx
!      print *, Phi%Value%Grid(i,1,1,1:4,1)
!    enddo
!    print *
!    do i=1,nx
!      print *, Psi%Value%Grid(i,1,1,1:4,1)
!    enddo
!    print *    
!    print *
!    norm = 0.0_dp
!    do l=1,4
!      do k=1,nz
!        do j=1,ny
!          do i=1,nx
!            norm(1)= norm(1) + Psi%Value%Grid(i,j,k,l,1)**2
!          enddo
!        enddo
!      enddo
!    enddo
!    
!    do l=1,4
!      do k=1,oldnz
!        do j=1,oldny
!          do i=1,oldnx
!            norm(2)= norm(2) + Phi%Value%Grid(i,j,k,l,1)**2
!          enddo
!        enddo
!      enddo
!    enddo
!    
!    print *, 'Norm', sqrt(norm(1) * dv), sqrt(norm(2) * dv*2)*2
!    print *

  end function InterpolateSpwf

end module Transform
