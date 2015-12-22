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
  use WaveFunctions
  use Spwfstorage
  use Spinors
  
  use CoulombDerivatives 

  implicit none

  !-----------------------------------------------------------------------------
  ! Interpolation coefficients for changing the mesh size, using the Lagrangian 
  ! formulas representation. Note that the coefficients have 2 components: 
  ! one for coordinates present in the Mesh module and one
  ! for coordinates that are accessible through symmetries.
  !-----------------------------------------------------------------------------
  real(KIND=dp), allocatable :: InterpolationCoefX(:,:,:)
  real(KIND=dp), allocatable :: InterpolationCoefY(:,:,:)
  real(KIND=dp), allocatable :: InterpolationCoefZ(:,:,:)

contains

  subroutine TransformInput(inTRC,inTSC,inIC,inPC,inSC,                        &
  &                         filenx,fileny,filenz,filenwt)
  !-----------------------------------------------------------------------------
  ! Transform the data read from the wavefunction file to a form appropriate
  ! for the different symmetry combinations. (This includes the densities!)
  !
  !-----------------------------------------------------------------------------
  logical, intent(in)       :: inTRC, inTSC, inIC, inPC, inSC
  integer, intent(in)       :: filenx,fileny,filenz,filenwt
  type(Spwf), allocatable   :: HFBasisTransformed(:)
  integer                   :: wave, index

  
  !-----------------------------------------------------------------------------
  !Checking for Time Reversal breaking first!
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
      call BreakTimeReversal(HFBasis(wave), HFBasisTransformed(wave),          &
      &                                     HFBasisTransformed(wave+index) )
    enddo

    ! Create the proton time-reversed pairs
    do wave=index+1,filenwt
      call BreakTimeReversal(HFBasis(wave),HFBasisTransformed(wave+index),     &
      &                                    HFBasisTransformed(wave+filenwt))
    enddo
  endif
  !-----------------------------------------------------------------------------
  !Checking for Parity Breaking
  if(inPC.neqv.PC) then
    !Breaking Parity
    do wave=1,size(HFBasisTransformed)
      call BreakParity(HFBasisTransformed(wave) )
    enddo
  endif
  !-----------------------------------------------------------------------------
  !Checking for Signature Breaking
  if(inSC.neqv.SC) then
    !Breaking Signature
    do wave=1,size(HFBasisTransformed)
      call BreakSignature(HFBasisTransformed(wave))
    enddo
  endif
  !-----------------------------------------------------------------------------
  !Checking for Time Simplex Breaking
  if(inTSC.neqv.TSC) then
    !Breaking Parity
    do wave=1,size(HFBasisTransformed)
      call BreakTimeSimplex(HFBasisTransformed(wave))
    enddo
  endif
  !-----------------------------------------------------------------------------
  !Applying the changes to the HFBasis
  call ChangeNumberWaveFunctions(nwt)
  do wave=1,nwt
    HFBasis(wave) = HFBasisTransformed(wave)
    call HFBasis(wave)%SymmetryOperators()
  enddo
  
  !-----------------------------------------------------------------------------
  ! Transforming the densities
  call TransformDensities(inSC, inTSC, inPC , filenx, fileny, filenz)

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
    integer                   ::  i,j,k,l, s
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
    Temp = NewSizeSpinor(size(Value%Grid,1), size(Value%Grid,2), nz)
    ! Moving the old values to their correct place in the new spinor.
    Temp%Grid(:,:,nz/2+1:nz ,:,:) = Value%Grid(:,:,1:nz/2,:,:)

    !Using the symmetries of the wavefunction to construct the rest
    do l=1,4
      s = CompSignExtension(3,Parity,Signature,TimeSimplex,l)
      do j=1, size(Value%Grid,2)
        do i=1,size(Value%Grid,1)
          Temp%Grid(i,j,1:nz/2,l,1) =                                          &
          & s*Reverse(LineExtensionCoulombZ(Value%Grid(:,:,:,l,1), &
          & nz/2,Parity,TimeSimplex,Signature,i,j))
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
    integer                   :: i,j,k,l,s
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
    Temp = NewSizeSpinor(nx,size(Value%Grid,2), size(Value%Grid,3))

    ! Moving the old values to their correct place in the new spinor.
    Temp%Grid(nx/2+1:nx,:,:,:,:) = Value%Grid(1:nx/2,:,:,:,:)

    !Using the symmetries of the wavefunction to construct the rest.
    do l=1,4
      s = CompSignExtension(1,Parity,Signature,TimeSimplex,l)
      do k=1,size(Value%Grid,3)
        do j=1,size(Value%Grid,2)
          Temp%Grid(1:nx/2,j,k,l,1) = real(s, KIND=dp)* & 
          & Reverse(LineExtensionCoulombX(Value%Grid(:,:,:,l,1), &
          & nx/2,Parity,TimeSimplex,Signature,j,k))
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
    integer                   :: i,j,k,l,s
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
    if(TimeSimplex.eq.0) call stp("Signature is already broken!")
          
    !When breaking parity, we need to store the values for the entire y-axis.
    Temp=NewSizeSpinor(size(Value%Grid,1),ny,size(Value%Grid,3))

    ! Moving the old values to their correct place in the new spinor.
    Temp%Grid(:, ny/2+1:ny,:,:,:) = Value%Grid(:,1:ny/2,:,:,:)

    !Using the symmetries of the wavefunction to construct the rest
    do l=1,4
      s = CompSignExtension(2,Parity,Signature,TimeSimplex,l)
      do k=1,size(Value%Grid,3)
        do i=1,size(Value%Grid,1)
          Temp%Grid(i,1:ny/2,k,l,1) = s*Reverse( &
          & LineExtensionCoulombY(Value%Grid(:,:,:,l,1), &
            & ny/2,Parity,TimeSimplex,Signature,i,k))
        enddo
      enddo
    enddo

    Temp = Conj(Temp)

    !Creating a new wavefunction, with the same quantum numbers, 
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

  function TransformKappa( InKappa, inPC, inIC, InputBlocksizes) result(OutKappa)
    !----------------------------------------------------------------------------
    ! Function that takes Kappa from a file and checks if it needs transformation
    ! when either Parity or Isospin gets broken.
    !
    !----------------------------------------------------------------------------
    use HFB
    use Spwfstorage

    complex(KIND=dp)              :: InKappa(:,:,:,:)
    complex(KIND=dp), allocatable :: OutKappa(:,:,:,:), Temp(:)
    logical, intent(in)           :: inPC, inIC
    integer                       :: sizes(2), P, it, iplus,imin,i,j

    integer, intent(in)           :: InputBlockSizes(2,2)

    if(.not. inIC) then
        call stp('No rules to transform Kappa yet when breaking Isospin.')
    endif

    if(PC .eqv. inPC) then
        ! File matches rundata ; no action necessary
        OutKappa = InKappa
    else
        !-----------------------------------------------------------------------
        ! Parity is conserved on file, broken in run.
        ! On file, we have two matrices:
        !
        !          (  0    K^+ )           (  0    K^- )
        !          (           )   and     (           )
        !          ( -K^+  0   )           ( -K^-  0   )
        ! 
        ! HOWEVER, taking as the following as NEW kappa is NOT ALLOWED!
        !
        !          (  0    K^+   0    0   )
        !          ( -K^+   0    0    0   )
        !          (  0     0    0    K^- )
        !          (  0     0   -K^-  0   )
        !
        ! because the wavefunction indices in the HFB module need not
        ! be in the correct order.
        !
        ! The correct form is (using the order of wavefunctions as generated
        ! by the NIL8 code)
        !         ( 0    0    K^+ 0  )
        !         ( 0    0    0   K^-)
        !         ( -K^+ 0    0   0  )
        !         ( 0   -K^-  0   0  )
        !
        !-----------------------------------------------------------------------

        ! We can take this size anyway, since P needs to be broken
        ! and isospin conserved.
        allocate(OutKappa(HFBSize,HFBSize,1,2)) ; OutKappa = 0.0_dp
        do it=1,2
            sizes = InputBlocksizes(:,it)
            
            OutKappa(sum(sizes)/2+1:sum(sizes)/2+sizes(2)/2,1:sizes(2)/2,1,it)            &
            & = InKappa(sizes(2)/2+1:sizes(2),1:sizes(2)/2,2,it)
            
            OutKappa(sum(sizes)/2+sizes(2)/2+1:sum(sizes),sizes(2)/2+1:sum(sizes)/2,1,it )&
            & = InKappa(sizes(1)/2+1:sizes(1),1:sizes(1)/2,1,it)
            
            do i=1,HFBSize
              do j=i+1, HFBSize
                OutKappa(i,j,1,it) = - OutKappa(j,i,1,it)
              enddo
            enddo
            
!             print *, 'it = ', it
!             print *, 'sizes', sizes
!             print *, 'In, negative parity'
!             do j=1,sizes(1)
!               print *, DBLE(Inkappa( j, 1:sizes(1), 1,it))
!             enddo
!             print * 
!             print *, 'In, positive parity'
!             do j=1,sizes(2)
!               print *, DBLE(Inkappa(j , 1:sizes(2), 2,it))
!             enddo
!             print *
!             print *, 'out, all parity'
!             do j=1,sizes(1) + sizes(2)
!               print *, DBLE(OutKappa(j, 1:sum(sizes), 1, it ))
!             enddo
!             print *
        enddo
    endif

  end function TransformKappa

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
    call NewWfOne%SetSignature(1)
    call NewWfTwo%SetSignature(-1)

    !Switching the components of nr.2 around
    ! Action of timereversal is -i *\sigma_y * K
    Temp = wf%GetValue()
    Temp = TimeReverse(Temp)

    call NewWfTwo%SetGrid(Temp)
  end subroutine BreakTimeReversal

end module Transform
