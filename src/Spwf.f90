!-------------------------------------------------------------------------------
! This module contains the type definition of the objects called Spwf.
!-------------------------------------------------------------------------------
module WaveFunctions

  use Compilationinfo
  use GenInfo
  use Mesh
  use Spinors
  implicit none  

  public

  type Spwf 
    private 

    !---------------------------------------------------------------------------
    ! Value
    !   Values of the associated wavefunction on the lattice
    ! Der
    !   Values of the partial derivatives of the WF.
    ! Lap 
    !   Value of the Laplacian of the WF.
    ! Occupation
    !   Occupation factor of the Spwf. In the case of HFB calculations, it is 
    !   the diagonal element of RhoHFB for HF basis and the occupation for the
    !   canonical basis.
    ! Dispersion
    !   Dispersion of the energy of the Spwf.
    ! Norm
    !   Norm of the wavefunction
    ! Angmoment
    !  Expectation value of the angular moment, in the three directions(x/y/z)
    ! J2
    !  Quadratic expectation values of the angular moment, <J_i^2>, in the three
    !  directions, in general not calculated.
    ! AngQuantum
    !  Quantum number of the angular momentum, meaning j so that <J^2> = j*(j+1)
    !
    !---------------------------------------------------------------------------   
    type(Spinor) :: Value           
    type(Spinor) :: Der(3)          
    type(Spinor) :: Lap             
    real(KIND=dp) :: Occupation     
    real(KIND=dp) :: Energy         
    real(KIND=dp) :: Dispersion     
    real(KIND=dp) :: Norm           
    real(KIND=dp) :: AngMoment(3)    
    real(KIND=dp) :: J2(3)
    real(KIND=dp) :: AngQuantum
    !---------------------------------------------------------------------------
    ! The following integers are the quantum numbers of the wavefunction.
    ! For each of them, 0 indicates that the symmetry is broken.
    ! They are defined by the following relations: 
    !
    ! [\check{S}^T_y \Psi] (x,y,z,\sigma) = 
    !                     \Psi^*(x,-y,z,\sigma) = \Psi(x,y,z,\sigma)
    !
    ! [\hat{P} \Psi] (x,y,z,\sigma)
    !                    = \Psi(-x,-y,-z,\sigma) = p \Psi(x,y,z,\sigma)
    !
    ! \hat{R}_z \Psi(x,y,z,\sigma) 
    !                    = -i \sigma \Psi(-x,-y,-z, \sigma)
    !                    = s \sigma \Psi(x,y,z,\sigma)
    !
    ! Note that signature in this code(s) and the usual quantum number \eta are 
    ! related: s = -i\eta. TimeSimplex can only be 1 or 0, since we can 
    ! arbitrarily choose the eigenvalue of that operator. TimeReversal can only 
    ! be 1 or 0, since there are no eigenvalues to be associated with this 
    ! operator.
    !---------------------------------------------------------------------------
    integer :: Isospin,TimeSimplex,Parity,Signature,TimeReversal
    !---------------------------------------------------------------------------
    !Pairing properties
    ! eqp   = quasi-particle energy
    !---------------------------------------------------------------------------
    real(KIND=dp) :: Eqp=1.0_dp
    !---------------------------------------------------------------------------
    ! Real version of the quantum numbers. They are the expected values of the 
    ! symmetry operators.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: IsospinR, TimeSimplexR, ParityR, SignatureR, TimeReversalR
    !---------------------------------------------------------------------------
    !sqrt(< R^2 >) for every wavefunction, very useful in debugging.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: RMSRadius=0.0_dp
    !---------------------------------------------------------------------------
    ! Delta : pairing gap for this wavefunction. Delta_{i, ibar} in BCS,
    ! maximum of the pairing gaps for this wavefunction in HFB. 
    !---------------------------------------------------------------------------
    complex(KIND=dp) :: Delta=0.0_dp
    !---------------------------------------------------------------------------
    ! PairPartner : integer indicating the most important pairing partner,
    ! of course only useful for HFB calculations.
    !---------------------------------------------------------------------------
    integer          :: PairPartner=0

  contains
    !---------------------------------------------------------------------------
    !Getter Procedures
    !---------------------------------------------------------------------------
    procedure, pass, public :: GetDer
    procedure, pass, public :: GetLap
    procedure, pass, public :: GetEnergy
    procedure, pass, public :: GetDispersion
    procedure, pass, public :: GetNorm
    procedure, pass, public :: GetOcc
    procedure, pass, public :: GetQNumbers
    procedure, pass, public :: GetIsospin
    procedure, pass, public :: GetSignature
    procedure, pass, public :: GetSignatureR
    procedure, pass, public :: GetParity
    procedure, pass, public :: GetParityR
    procedure, pass, public :: GetTimeSimplex
    procedure, pass, public :: GetTimeReversal
    procedure, pass, public :: GetValue
    procedure, pass, public :: GetAngMoment
    procedure, pass, public :: GetEqp
    procedure, pass, public :: GetAngQuantum
    procedure, pass, public :: GetAngMomentSquared
    !---------------------------------------------------------------------------
    !Setter Procedures
    !---------------------------------------------------------------------------
    procedure, pass, public :: SetEnergy
    procedure, pass, public :: SetDispersion
    procedure, pass, public :: SetIsospin        
    procedure, pass, public :: SetSignature        
    procedure, pass, public :: SetParity
    procedure, pass, public :: SetTimeSimplex
    procedure, pass, public :: SetTimeReversal
    procedure, pass, public :: SetGridMesh
    procedure, pass, public :: SetGridSpinor
    procedure, pass, public :: SetGridComponent
    procedure, pass, public :: SetAngMoment
    procedure, pass, public :: SetAngMomentSquared
    procedure, pass, public :: Seteqp
    procedure, pass, public :: SetDelta
    procedure, pass, public :: SetPairPartner
    procedure, pass, public :: SetSignatureR
    ! Overloading so creating wavefunctions is easier.
    generic :: SetGrid => SetGridMesh,SetGridSpinor,SetGridComponent 
    !---------------------------------------------------------------------------
    !Density Computing Routines
    !---------------------------------------------------------------------------
    procedure, pass, public :: GetTau
    procedure, pass, public :: GetDensity
    procedure, pass, public :: GetDerRho
    procedure, pass, public :: GetLapRho
    procedure, pass, public :: GetVecj
    procedure, pass, public :: GetVecs
    procedure, pass, public :: GetVecT
    procedure, pass, public :: GetVecF
    procedure, pass, public :: GetJMuNu
    procedure, pass, public :: GetNablaJ
    !---------------------------------------------------------------------------
    !Miscellaneous
    !---------------------------------------------------------------------------
    !procedure, pass, public :: PrintInfo
    procedure, pass, public :: PrintHF
    procedure, pass, public :: PrintCanonical
    procedure, pass, public :: CompNorm        
    procedure, pass, public :: SetOcc
    procedure, pass, public :: CompDer
    procedure, pass, public :: CompAngMoment
    procedure, pass, public :: SymmetryOperators
    procedure, pass, public :: CalcRMSRadius
    procedure, pass, public :: ResetWf
    !---------------------------------------------------------------------------
    !In/Output
    !---------------------------------------------------------------------------
    procedure, pass, public :: ReadSpwf
    procedure, pass, public :: WriteSpwf
    generic :: Read => ReadSpwf
    generic :: Write=> WriteSpwf
  end type Spwf
  !-----------------------------------------------------------------------------
  !Operators
  !-----------------------------------------------------------------------------
  interface operator (+)
    !Overloading "+" to be used on wavefunctions.
    module procedure AddSpwf
  end interface
  interface operator (-)
    !Overloading "-" to be used on wavefunctions.
    module procedure MinusSpwf
  end interface
  interface operator (*)
    !Overloading "*" to be used with numbers and wavefunctions.
    module procedure MultiplySpwf
  end interface
   
contains
  function AddSpwf (WF1, WF2) result (SumWF)
  !-----------------------------------------------------------------------------
  ! Sum operator for wavefunctions
  !-----------------------------------------------------------------------------
  ! Attention please: the derivatives are not added together! This is never 
  ! useful in MOCCa and thus was removed.
  !-----------------------------------------------------------------------------
    type(Spwf), intent(in) :: WF1, WF2
    type(Spwf)             :: SumWF
       
    !Some initial checks on symmetries.
    if(WF1%Parity.ne.WF2%Parity) then
      call stp("You can't add wavefunctions with different parity!")
    endif
    if(WF1%Signature.ne.WF2%Signature) then
      call stp("You can't add wavefunctions with different signature!")
    endif
    if(WF1%TimeSimplex.ne.WF2%TimeSimplex) then
      call stp("You can't add wavefunctions with different TimeSimplex!")
    endif
    if(WF1%Isospin.ne.WF2%Isospin) then
      call stp("You can't add wavefunctions with different isospin!")
    endif
 
    SumWF%Value    = WF1%Value  + WF2%Value
    
!    Completely unneccessary to add derivatives....
!    SumWF%Der(1)   = WF1%Der(1) + WF2%Der(1)
!    SumWF%Der(2)   = WF1%Der(2) + WF2%Der(2)
!    SumWF%Der(3)   = WF1%Der(3) + WF2%Der(3)
!    SumWF%Lap      = WF1%Lap    + WF2%Lap

    SumWF%Isospin      = WF1%Isospin
    SumWf%TimeSimplex  = WF1%TimeSimplex
    SumWF%Parity       = WF1%Parity
    SumWF%TimeReversal = WF1%TimeReversal
    SumWF%Signature    = WF1%Signature

    SumWF%Norm         = 0.0_dp
    SumWF%Energy       = 0.0_dp
    SumWF%Occupation   = 0.0_dp
    SumWF%Dispersion   = 0.0_dp
    SumWF%AngMoment    = 0.0_dp
    SumWF%J2           = 0.0_dp
  end function AddSpwf
  
  function MultiplySpwf(A, WF) result(AWF)
  
    type(Spwf), intent(in)       :: WF
    type(Spwf)                   :: AWF
    complex(KIND=dp), intent(in) :: A
    type(Spinor)                 :: Psi
  
    Psi = A * WF%Value
    AWF = NewWaveFunction(Psi, WF%Isospin, WF%TimeSimplex, WF%Parity, &
    &                     WF%Signature, WF%TimeReversal)
  
  end function MultiplySpwf
    
  function MinusSpwf (WF1, WF2) result (DiffWF)
  !-----------------------------------------------------------------------------
  ! Minus operator for wavefunctions
  !-----------------------------------------------------------------------------
    type(Spwf), intent(in) :: WF1, WF2
    type(Spwf)             :: DiffWF

    !Some initial checks on symmetries.
    if(WF1%Parity.ne.WF2%Parity) then
      call stp("You can't subtract wavefunctions with different parity!")
    endif
    if(WF1%Signature.ne.WF2%Signature) then
      call stp("You can't subtract wavefunctions with different signature!")
    endif
    if(WF1%TimeSimplex.ne.WF2%TimeSimplex) then
      call stp("You can't subtract wavefunctions with different TimeSimplex!")
    endif
    if(WF1%Isospin.ne.WF2%Isospin) then
      call stp("You can't subtract wavefunctions with different isospin!")
    endif

    DiffWF%Value    = WF1%Value - WF2%Value
    DiffWF%Der(1)   = WF1%Der(1) - WF2%Der(1)
    DiffWF%Der(2)   = WF1%Der(2) - WF2%Der(2)
    DiffWF%Der(3)   = WF1%Der(3) - WF2%Der(3)
    DiffWF%Lap      = WF1%Lap - WF2%Lap

    DiffWF%Isospin = WF1%Isospin
    DiffWf%TimeSimplex  = WF1%TimeSimplex
    DiffWF%Parity       = WF1%Parity
    DiffWF%TimeReversal = WF1%TimeReversal
    DiffWF%Signature    = WF1%Signature
      
    DiffWF%Norm=0.0_dp
    DiffWF%Energy=0.0_dp
    DiffWF%Occupation=0.0_dp
    DiffWF%Dispersion=0.0_dp
    DiffWF%AngMoment = 0.0_dp
    Diffwf%J2        = 0.0_dp
  end function MinusSpwf

  type(Spwf) function TimeReverseSpwf(Wf) result(TPsi)
    !--------------------------------------------------------------------------
    ! Function that returns the time-reversed partner of a wavefunction.
    !--------------------------------------------------------------------------
    ! Please note that this function only affects the actual value and quantum
    ! numbers, the derivatives and other quantities are unaffected! This is due
    ! to the fact that this function is only useful for HFB when time-reversal
    ! is conserved and Signature is broken. This function is only used in 
    ! the HFBoccupations subroutine. This should be modified in the future
    ! if more general use is needed.
    type(Spwf), intent(in) :: Wf
    
    TPsi%Value      = TimeReverse(WF%Value)
    TPsi%Isospin    = WF%Isospin
    TPsi%Parity     = WF%Parity
    TPsi%Signature  = -WF%Signature
    TPsi%TimeSimplex= WF%TimeSimplex
    TPsi%TimeReversal=WF%TimeReversal
    TPsi%Occupation = WF%Occupation
    TPsi%Energy     = WF%Energy
    TPsi%Norm       = WF%Norm
    TPsi%Dispersion = WF%Dispersion
    TPsi%AngMoment  = - WF%AngMoment
    TPsi%J2         = WF%J2
    TPsi%eqp        = WF%eqp
    TPsi%IsospinR   = WF%IsospinR
    TPsi%SignatureR = - WF%SignatureR
    TPsi%ParityR    = WF%ParityR
    TPsi%RMSRadius  = WF%RMSRadius
    TPsi%Delta      = WF%Delta
    TPsi%PairPartner= WF%Pairpartner
  end function TimeReverseSpwf

  subroutine ResetWF(WF)
    !--------------------------------------------------------------------------
    ! Zeroes the value of the wavefunction.
    !--------------------------------------------------------------------------
    class(Spwf), intent(inout) :: WF

    Wf%Value = NewSpinor()

  end subroutine ResetWf
  
  function NewWaveFunction                                                     &
  &(Value,Isospin,TimeSimplex,Parity,Signature,TimeReversal)result(WaveFunction)
    !---------------------------------------------------------------------------
    ! This function serves as a constructor for the SWPF class.
    ! Things that need to be computed are initialised to 0.
    !---------------------------------------------------------------------------
    type(Spinor), intent(in)   :: Value
    integer, intent(in)        :: Isospin,TimeSimplex,Parity
    integer, intent(in)        :: Signature,TimeReversal
    type(Spwf)                 :: WaveFunction

    WaveFunction%Value         = Value
    WaveFunction%Isospin       = Isospin
    WaveFunction%TimeSimplex   = TimeSimplex
    WaveFunction%Parity        = Parity
    WaveFunction%TimeReversal  = TimeReversal
    WaveFunction%Signature     = Signature
    WaveFunction%Norm          = 0.0_dp
    WaveFunction%Energy        = 0.0_dp
    WaveFunction%Occupation    = 0.0_dp
    WaveFunction%Dispersion    = 0.0_dp
    WaveFunction%AngMoment     = 0.0_dp
    WaveFunction%IsospinR      = 0.0_dp
    WaveFunction%SignatureR    = 0.0_dp
    WaveFunction%ParityR       = 0.0_dp
    WaveFunction%TimeSimplexR  = 0.0_dp
    WaveFunction%TimeReversalR = 0.0_dp
    WaveFunction%RMSRadius     = 0.0_dp
    WaveFunction%J2            = 0.0_dp
  end function NewWaveFunction
        
  function CopyWaveFunction(WF) result(Copy)
    !---------------------------------------------------------------------------
    ! This function makes a perfect copy of an Spwf. 
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: WF
    type (Spwf)             :: Copy

    Copy            = NewWaveFunction(WF%Value,WF%Isospin,WF%TimeSimplex,      &
    &                                 WF%Parity,WF%Signature,WF%TimeReversal)
    Copy%Der        = WF%Der
    Copy%Lap        = WF%Lap
    Copy%Occupation = WF%Occupation
    Copy%Energy     = WF%Energy
    Copy%Norm       = WF%Norm
    Copy%Dispersion = WF%Dispersion
    Copy%AngMoment  = WF%AngMoment
    Copy%J2         = WF%J2
    Copy%J2         = WF%J2
    Copy%eqp        = WF%eqp
    Copy%IsospinR   = WF%IsospinR
    Copy%SignatureR = WF%SignatureR
    Copy%ParityR    = WF%ParityR
    Copy%RMSRadius  = WF%RMSRadius
    Copy%Delta      = WF%Delta
    Copy%PairPartner=WF%Pairpartner

  end function CopyWaveFunction
        
  function GetValue(WF) result(Value)
    !---------------------------------------------------------------------------
    ! This function returns the value of the Spwf (which is a spinor)
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: WF
    type(Spinor)            :: Value

    Value = WF%Value
  end function GetValue
        
  subroutine SetGridSpinor(WaveFunction, Grid)
    !---------------------------------------------------------------------------
    !Function that assigns a spinor to Wavefunction Spinor
    ! Overloaded to "SetGrid" with SetGridMesh and SetGridComponent
    !---------------------------------------------------------------------------
    type(Spinor)  :: Grid
    class (Spwf)  :: WaveFunction
    
    WaveFunction%Value = Grid
  end subroutine SetGridSpinor
        
  subroutine SetGridMesh(WaveFunction, Grid)
    !---------------------------------------------------------------------------
    ! Function that assigns a Mesh to Wavefunction Spinor
    ! Overloaded to "SetGrid" with SetGridSpinor an SetGridComponent
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: Grid(:,:,:,:,:)
    class (Spwf)               :: WaveFunction

    WaveFunction%Value%Grid(:,:,:,:,:) = Grid
    return
  end subroutine SetGridMesh
        
  subroutine SetGridComponent(WaveFunction, Component,Grid,IsoComponent)
    !---------------------------------------------------------------------------
    ! Function that assigns values to one component of the Spinor of a spwf.
    ! Overloaded to "SetGrid" with SetGridSpinor and SetGridMesh
    !---------------------------------------------------------------------------                
    real(KIND=dp), intent(in)     :: Grid(:,:,:)
    class (Spwf)                  :: WaveFunction
    integer, intent(in)           :: Component
    integer, intent(in), optional :: IsoComponent
    if(.not.allocated(WaveFunction%Value%Grid)) WaveFunction%Value = NewSpinor()
    if(present(IsoComponent)) then
            call WaveFunction%Value%SetComponent(Component, Grid)
    else
            call WaveFunction%Value%SetComponent(Component, Grid, Isocomponent)
    endif
  end subroutine SetGridComponent

  function GetDensity(WaveFunction) result(Density)
    !---------------------------------------------------------------------------
    !Function returns the density \rho 
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: WaveFunction
    real(KIND=dp)           :: Density(nx,ny,nz)
    !This type of thing is exactly why the class "Spinor" exists.
    if(.not.allocated(Wavefunction%Value%Grid)) call stp('Den')
    Density = RealMultiplySpinor(WaveFunction%Value,WaveFunction%Value)
  
  end function GetDensity
        
  pure real(KIND=dp) function GetEnergy(WaveFunction) result(Energy)
    class (Spwf), intent (in) :: WaveFunction
    Energy=WaveFunction%Energy
  end function GetEnergy
        
  real(KIND=dp) function GetEqp(WaveFunction) result(Eqp)
    class (Spwf), intent (in) :: WaveFunction
    Eqp=WaveFunction%eqp
  end function Geteqp

  subroutine SetEqp(WaveFunction, eqp)
    class(Spwf), intent(inout) :: WaveFunction
    real(KIND=dp), intent(in)  :: eqp
    WaveFunction%eqp = eqp
  end subroutine Seteqp
  
  subroutine SetDelta(WaveFunction, Delta)
    class(Spwf), intent(inout) :: WaveFunction
    complex(KIND=dp), intent(in)  :: Delta
    WaveFunction%Delta = Delta
  end subroutine SetDelta
  
  subroutine SetPairPartner(WaveFunction, Partner)
    class(Spwf), intent(inout) :: WaveFunction
    integer, intent(in)        :: Partner
    WaveFunction%PairPartner = Partner
  end subroutine SetPairPartner
  
  real(KIND=dp) function GetDispersion(WaveFunction) result(Dispersion)
    class (Spwf), intent (in) :: WaveFunction
    Dispersion=WaveFunction%Dispersion
  end function GetDispersion

  real(KIND=dp) function GetAngMoment(WaveFunction, Direction) result(AngMoment)
    class (Spwf), intent (in) :: WaveFunction
    integer, intent(in)       :: Direction            
    AngMoment=WaveFunction%AngMoment(Direction)
    return
  end function GetAngMoment
  
  real(KIND=dp) function GetAngMomentSquared(Wavefunction, Direction) result(J2)
    class (Spwf), intent (in) :: WaveFunction
    integer, intent(in)       :: Direction            
    J2=WaveFunction%J2(Direction)
    return
  end function GetAngMomentSquared  
  
 real(KIND=dp) function GetAngQuantum(Wavefunction) result(J)
    class (Spwf), intent (in) :: WaveFunction       
    J =WaveFunction%AngQuantum
    return
  end function GetAngQuantum  
        
  subroutine SetAngMoment(WaveFunction, AngMoment)
    class(Spwf), intent(inout) :: WaveFunction
    real(KIND=dp), intent(in)  :: AngMoment(3)
    WaveFunction%AngMoment     = AngMoment
  end subroutine SetAngMoment
  
  subroutine SetAngMomentSquared(WaveFunction, AngMoment)
    class(Spwf), intent(inout) :: WaveFunction
    real(KIND=dp), intent(in)  :: AngMoment(3)
    WaveFunction%J2     = AngMoment
  end subroutine SetAngMomentSquared

  subroutine CompAngMoment(WaveFunction)
    !---------------------------------------------------------------------------
    ! Computes and set the correct values for the AngMoment.
    !---------------------------------------------------------------------------
    class(Spwf), intent(inout)    :: WaveFunction
    real(KIND=dp)                 :: Temp(6,2)
    Temp = 0.0_dp
    Temp = AngularMomentum(WaveFunction, WaveFunction, .true.)               

    !Angular Momentum in the three directions
    !Note: only the real parts!
    WaveFunction%AngMoment   = Temp(1:3,1)
    WaveFunction%J2          = Temp(4:6,1)
    !---------------------------------------------------------------------------
    ! Now find j so that <J^2> = j*(j+1)
    ! The solution is obviously:
    ! j = [- 1 + sqrt( 1 + 4 * <J^2>)]/2
    WaveFunction%AngQuantum=                                                   &
    &                  - 0.5_dp*(1.0_dp-sqrt( 1.0_dp + 4.0_dp*sum(Temp(4:6,1)))) 
  end subroutine CompAngMoment
 
  subroutine SetEnergy(WaveFunction,Energy) 
    !---------------------------------------------------------------------------
    ! Function sets the energy of a Wavefunction.
    !---------------------------------------------------------------------------
    real(KIND=dp) ,intent(in)     :: Energy
    class (Spwf), intent(inout) :: WaveFunction
    WaveFunction%Energy=Energy
  end subroutine SetEnergy
        
   subroutine SetDispersion(WaveFunction,Dispersion) 
    !---------------------------------------------------------------------------
    ! Function sets the energy of a Wavefunction.
    !---------------------------------------------------------------------------
    real(KIND=dp) ,intent(in)     :: Dispersion
    class (Spwf), intent(inout) :: WaveFunction
    WaveFunction%Dispersion=Dispersion
  end subroutine SetDispersion

  subroutine CompNorm(WaveFunction) 
    !---------------------------------------------------------------------------
    ! This subroutine computes the norm of the Spwf and assigns it to the
    ! wavefunction. 
    !---------------------------------------------------------------------------
    class (Spwf), intent(inout) :: WaveFunction
    ! int_{-\infty}^{\infty} |\Psi|**2 dv \approx 
    !                    dv*\Sum_{box} |\Psi(spin up)|**2 + |\Psi(spin down)|**2 
    !Using the spinor inproduct.
    WaveFunction%Norm = InproductSpinorReal(WaveFunction%Value,WaveFunction%Value) 
  end subroutine CompNorm

  real(KIND=dp) function GetNorm(WaveFunction) result(Norm)
    !---------------------------------------------------------------------------
    ! This function gets you the norm of the wavefunction.
    !---------------------------------------------------------------------------
    class (Spwf), intent(in) :: WaveFunction
    Norm=WaveFunction%Norm
  end function GetNorm
        
  subroutine SetOcc(WaveFunction, Occupation)
    !---------------------------------------------------------------------------
    ! This function assigns the occupation factor of a Spwf.
    !---------------------------------------------------------------------------
    class(Spwf) :: WaveFunction
    real(KIND=dp), intent(in) :: Occupation
    WaveFunction%Occupation=Occupation
  end subroutine SetOcc

  pure real(KIND=dp) function GetOcc(WaveFunction) result(Occupation)
    class(Spwf), intent(in) :: WaveFunction
    Occupation=WaveFunction%Occupation
    return
  end function GetOcc
        
  subroutine SetIsospin(WaveFunction, Isospin)
    !---------------------------------------------------------------------------
    ! Function that changes the isospin of a Spwf.
    !---------------------------------------------------------------------------
    integer :: Isospin
    class (Spwf) :: WaveFunction

    if((Isospin.ne.1).and.(Isospin.ne.-1)) then
        call stp('Invalid Isospin assignment')
    endif
    WaveFunction%Isospin=Isospin
    return
  end subroutine SetIsospin 

  pure integer function GetIsospin(WaveFunction) result(Isospin)
    !---------------------------------------------------------------------------
    ! Function that outputs the isospin of a Spwf.
    !---------------------------------------------------------------------------
    class (Spwf), intent(in) :: WaveFunction

    Isospin=Wavefunction%Isospin
  end function GetIsospin

  subroutine SetTimeReversal(WaveFunction, TimeReversal)
    !---------------------------------------------------------------------------
    ! Function that changes the isospin of a Spwf.
    !---------------------------------------------------------------------------
    integer      :: TimeReversal
    class (Spwf) :: WaveFunction
    if((TimeReversal.ne.1).and.(TimeReversal.ne.0)) then
        call stp('Invalid TimeReversal assignment')
    endif
    WaveFunction%TimeReversal=TimeReversal
  end subroutine SetTimeReversal 

  integer function GetTimeReversal(WaveFunction) result(TimeReversal)
    !---------------------------------------------------------------------------
    ! Function that outputs the isospin of a Spwf.
    !---------------------------------------------------------------------------
    class (Spwf), intent(in) :: WaveFunction

    TimeReversal=Wavefunction%TimeReversal
  end function GetTimeReversal

  subroutine SetTimeSimplex(WaveFunction, TimeSimplex)
    !---------------------------------------------------------------------------
    ! Function that changes the isospin of a Spwf.
    !---------------------------------------------------------------------------
    integer :: TimeSimplex
    class (Spwf) :: WaveFunction

    if((TimeSimplex.ne.1).and.(TimeSimplex.ne.0)) then
        call stp('Invalid TimeSimplex assignment')
    endif
    WaveFunction%TimeSimplex=TimeSimplex
    return
  end subroutine SetTimeSimplex 

  pure integer function GetTimeSimplex(WaveFunction) result(TimeSimplex)
    !---------------------------------------------------------------------------
    ! Function that outputs the isospin of a Spwf.
    !---------------------------------------------------------------------------
    class (Spwf), intent(in) :: WaveFunction

    TimeSimplex=Wavefunction%TimeSimplex
  end function GetTimeSimplex

  subroutine SetParity(WaveFunction, Parity)
    !---------------------------------------------------------------------------
    ! Function that changes the parity of a Spwf.
    !---------------------------------------------------------------------------
    integer :: Parity
    class (Spwf) :: WaveFunction

    if((Parity.ne.1).and.(Parity.ne.-1).and.(Parity.ne.0)) then
      call stp('Invalid Parity Assignment')
    endif
    WaveFunction%Parity=Parity
  end subroutine SetParity 

  pure integer function GetParity(WaveFunction) result(Parity)
    !---------------------------------------------------------------------------
    ! Function that outputs the parity of a Spwf.
    !---------------------------------------------------------------------------
    class (Spwf), intent(in) :: WaveFunction

    Parity=Wavefunction%Parity
  end function GetParity
  
  pure real(KIND=dp) function GetParityR(WaveFunction) result(ParityR)
    !---------------------------------------------------------------------------
    ! Function that outputs the parity of a Spwf.
    !---------------------------------------------------------------------------
    class (Spwf), intent(in) :: WaveFunction

    ParityR=Wavefunction%ParityR
  end function GetParityR

  subroutine SetSignature(WaveFunction, Signature)
    !---------------------------------------------------------------------------
    ! Function that changes the parity of a Spwf.
    !---------------------------------------------------------------------------
    integer, intent(in) :: Signature
    class (Spwf)        :: WaveFunction

    if((Signature.ne.1).and.(Signature.ne.0).and.(Signature.ne.-1)) then
        call stp('Invalid Signature Assignment')
    endif
    WaveFunction%Signature=Signature
    return
  end subroutine SetSignature 

  pure integer function GetSignature(WaveFunction) result(Signature)
    !---------------------------------------------------------------------------
    ! Function that outputs the parity of a Spwf.
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WaveFunction

    Signature=Wavefunction%Signature
  end function GetSignature

  pure real(KIND=dp) function GetSignatureR(WaveFunction) result(SignatureR)
    !---------------------------------------------------------------------------
    ! Function that outputs the parity of a Spwf.
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WaveFunction

    SignatureR=Wavefunction%SignatureR
  end function GetSignatureR

  subroutine SetSignatureR(WaveFunction, SR)
    !---------------------------------------------------------------------------
    ! Function that outputs the parity of a Spwf.
    !---------------------------------------------------------------------------
    class (Spwf)              :: WaveFunction
    real(KIND=dp), intent(in) :: SR
    Wavefunction%SignatureR = SR
  end subroutine SetSignatureR

  subroutine CompDer(WF)
    !---------------------------------------------------------------------------
    !This function computes all relevant derivatives of the wavefunction: first
    !order partial derivatives and the Laplacian.
    !---------------------------------------------------------------------------
    integer      :: P,S,TS
    class (Spwf) :: WF

    P=WF%Parity
    S=WF%Signature
    TS=WF%TimeSimplex
    
    WF%Der = DeriveSpinor(WF%Value,P,S,TS)
    WF%Lap = LapSpinor(WF%Value,P,S,TS)    
  end subroutine CompDer
        
 function GetDer(WF,Direction) result(Der)
    !---------------------------------------------------------------------------
    ! Function that returns the laplacian of the Spwf.
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: WF
    type(Spinor)            :: Der
    integer, intent(in)     :: Direction
    
    Der = WF%Der(Direction)
 end function GetDer
        
  function GetLap(WF) result(Lap)
    !---------------------------------------------------------------------------
    ! Function that returns the laplacian of the Spwf.
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: WF
    type(Spinor)            :: Lap
                  
    Lap = WF%Lap
  end function GetLap
             
  pure function InProduct(WFOne,WFTwo) result(IP)
    !---------------------------------------------------------------------------
    !This function computes the (possibly complex) inproduct of two Spwfs.
    !                   < WFOne|WFTwo >
    ! Note the ordering of the arguments!
    !---------------------------------------------------------------------------
    class (Spwf),intent(in)  :: WFOne,WFTwo
    real(KIND=dp)            :: IP(2) !Real and imaginary parts.

    IP =0.0_dp
    if(  (WFOne%Isospin.ne.WFTwo%Isospin).or.& 
    &     (WFOne%Parity.ne.WFTwo%Parity).or.&
    &  (WFOne%Signature.ne.WFTwo%Signature)) then
    !If the isospin, signature and the parity quantum numbers are not equal, 
    !the inproduct is exactly zero.
          return
    else
          IP(1) = InproductSpinorReal(WFOne%Value, WFTwo%Value)
          IP(2) = 0.0_dp
          if(.not.TSC) then
              !Only compute this if time simplex is not conserved.
              IP(2) = InproductSpinorImaginary(WFOne%Value, WFTwo%Value)                     
          endif
    endif
  end function InProduct
        
  pure function GetQNumbers(WF) result(Numbers)
    !---------------------------------------------------------------------------
    !This function returns the quantum numbers of the Wavefunction.
    !---------------------------------------------------------------------------
    class (Spwf), intent(in) :: WF
    integer     :: Numbers(5)

    Numbers=(/WF%Isospin,WF%TimeSimplex,WF%Signature,WF%Parity,WF%TimeReversal/)
    return
  end function GetQNumbers

  function AngMomOperator(WF, Direction) result(ActionOfJ)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the action of the angular momentum operator on a 
    ! wavefunction WF. 
    ! Input is the wavefunction, but the output is a spinor.
    !---------------------------------------------------------------------------
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !
    !       J_y = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 )
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: WF
    integer, intent(in)     :: Direction
    ! Temp1 contains the spin part, temp2 and temp3 the orbital parts.   
    type(Spinor)            :: ActionOfJ, Temp1, Temp2, Temp3, Temp4
    integer                 :: i
    
    Temp1=NewSpinor();Temp2=NewSpinor(); Temp3=NewSpinor();Temp4=NewSpinor()
    ActionOfJ = NewSpinor()
    
    select case(Direction)
    case (1)
      !The action of Jx
      Temp1 = Pauli(WF%Value,1)
      Temp1 = (1.0_dp/2.0_dp)* Temp1
      do i=1,4
        Temp2%Grid(:,:,:,i,1) = Mesh3D(2,:,:,:) * WF%Der(3)%Grid(:,:,:,i,1)
      enddo
      do i=1,4
        Temp3%Grid(:,:,:,i,1) = Mesh3D(3,:,:,:) * WF%Der(2)%Grid(:,:,:,i,1)
      enddo
    case (2)
      !The action of Jy
      Temp1 = Pauli(WF%Value,2)
      Temp1 = (1.0_dp/2.0_dp)* Temp1
      do i=1,4
        Temp2%Grid(:,:,:,i,1) = Mesh3D(3,:,:,:) * WF%Der(1)%Grid(:,:,:,i,1)
      enddo
      do i=1,4
        Temp3%Grid(:,:,:,i,1) = Mesh3D(1,:,:,:) * WF%Der(3)%Grid(:,:,:,i,1)
      enddo
    case (3)
      !The action of Jz
      Temp1 = Pauli(WF%Value,3)
      Temp1 = (1.0_dp/2.0_dp)*Temp1
      do i=1,4
        Temp2%Grid(:,:,:,i,1) = Mesh3D(1,:,:,:) * WF%Der(2)%Grid(:,:,:,i,1)
      enddo      
      do i=1,4
        Temp3%Grid(:,:,:,i,1) = Mesh3D(2,:,:,:) * WF%Der(1)%Grid(:,:,:,i,1)
      enddo
    end select
    
    Temp4 = Temp3 - Temp2
    Temp4 = MultiplyI(Temp4)
    ActionOfJ = Temp1 + Temp4
  end function AngMomOperator

  function AngularMomentum(WF2,WF1, quadratic) result(AngMom)
    !---------------------------------------------------------------------------
    ! This function returns the matrix element 
    !       < WF1 | \hat{J} | WF2 > 
    ! and the matrix elements of
    !       < WF1 | J_i**2 | WF2 >      
    ! if the quadratic is true.            
    ! Notice the ordering of the arguments in this expression!
    !---------------------------------------------------------------------------
    ! Note:
    ! - When dictated by time reversal symmetry, it returns zero. 
    !
    ! - The matrix elements < WF2 | \hat{J} | WF1 > are only guaranteed to be 
    !   real in the case of time simplex symmetry. This is the reason why J has 
    !   two columns: the real and imaginary part.
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: WF1,WF2
    type(Spinor)            :: Intermediate(3), Intermediate2(3)
    real(KIND=dp)           :: AngMom(6,2)
    integer                 :: i
    logical, intent(in)     :: Quadratic
    
    AngMom=0.0_dp
    if(WF1%Parity.ne.WF2%Parity .or. WF1%Isospin.ne.WF2%Isospin) then
        ! The matrix element is zero when opposite parities or isospins
        ! are involved.
        return
    endif

    do i=1,3
      Intermediate(i) = NewSpinor()
      Intermediate(i) = AngMomOperator(wf1,i)
    enddo
    !X Component
    if(SC                                                                      &
    &  .and. wf1%GetSignature().eq.wf2%GetSignature()) then
      !Good signature quantum numbers prohibit <Jx>
      AngMom(1,:) = 0.0_dp
    else    
      AngMom(1,1) = InproductSpinorReal(wf2%Value, Intermediate(1))
      AngMom(1,2) = InproductSpinorImaginary(wf2%Value, Intermediate(1))
    endif
    
    !Y Component
    if(SC                                                                      &
    &  .and. wf1%GetTimeSimplex().eq.wf2%GetTimeSimplex()) then
      !Good time simplex quantum numbers prohibit <Jy>
      AngMom(2,:) = 0.0_dp
    else
      AngMom(2,1) = InproductSpinorReal(wf2%Value, Intermediate(2))
      AngMom(2,2) = InproductSpinorImaginary(wf2%Value, Intermediate(2))
    endif 
    
    !Z Component
    if( SC  .and.                                                              &
        & wf1%GetSignature() .ne. wf2%Getsignature() ) then
        AngMom(3,:) = 0.0
    else

        AngMom(3,1) = InproductSpinorReal(wf2%Value, Intermediate(3))
        AngMom(3,2) = InproductSpinorImaginary(wf2%Value, Intermediate(3))
    endif
    if(Quadratic) then
      !J**2 components, only compute these if desired.
      do i=1,3
        Intermediate2(i) = AngMomOperator(wf2,i)
        AngMom(3+i,1)=InproductSpinorReal(Intermediate2(i),Intermediate(i))
        AngMom(3+i,2)=InproductSpinorImaginary(Intermediate2(i),Intermediate(i))                
      enddo
    endif  
  end function AngularMomentum

  subroutine PrintHF(WF,i,PrintType)
  !-----------------------------------------------------------------------------
  ! Prints the info of a normal Spwf function in the HFBasis, dependent on the
  ! value of PrintType.
  !-----------------------------------------------------------------------------
  ! PrintType
  ! = 1 => Minimal printout for HF calculations, no pairing.
  !        fmt =  n, <P> , (<S>) , v^2 , E , d^2h , (<Jx>), (<Jy>), <Jz>
  ! = 2 => Minimal printout for BCS calculations.
  !        fmt =  n, <P> , (<S>) , v^2, Delta , E , d^2h , (<Jx>), (<Jy>), <Jz>
  ! = 3 => Minimal printout for HFB calculations
  !        fmt =  n, <P> , (<S>) , Rho_{i,i}, Delta, Pairing Partner,
  !        E , d^2h , (<Jx>), (<Jy>), <Jz>
  !-----------------------------------------------------------------------------
  ! Note that: 
  ! * Values between brackets only get printed in the case of some 
  !   symmetry breakings.
  !  
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! Hartree-Fock calculations
  !           n        <P>        v^2       E_sp      d2H       <Jz>    < r^2 >
  1  format ( i3, 1x, f5.2, 1x, f7.4, 1x, f8.3, 1x, e9.2, 1x, f6.2,      f6.2 ) 
  !           n        <P>       <S>        v^2       E_sp      d2H       <Jz>
  11 format ( i3, 1x, f5.2, 1x, f5.2, 1x, f7.4 , 1x, f8.3, 1x, e9.2, 1x, f6.2, & 
  &           f6.2,1x, f6.2 )
  !           J        <r^2>
  !  
  !           n     <P> <S>        v^2       E_sp      d2H      <Jx> <Jz> <r^2>
  12 format ( i3, 2(1x, f5.2), 1x, f7.4 , 1x, f8.3, 1x, e9.2, 3(1x, f6.2) )
  !           n     <P> <S>        v^2       E_sp      d2H      <Jx> <Jy> <r^2>
  13 format ( i3, 2(1x, f5.2), 1x, f7.4 , 1x, f8.3, 1x, e9.2, 3(1x, f6.2) )
  
  !-----------------------------------------------------------------------------
  ! BCS calculations
  !           n        <P> <S>   v^2       Delta     E_sp    d2H <Jz>,J,<r^2>
  2  format (i3,1x, f5.2,  1x, f7.4, 1x, f7.2 ,1x, f8.3, 1x, e9.2, 3(1x,f6.2))
  !           n     <P> <S>    v^2  Delta     E_sp      d2H     <Jx> <Jz>J<r^2>
  21 format (i3,2(1x, f5.2), 1x,f7.4,1x,f7.2,1x, f8.3, 1x, e9.2, 4(1x, f6.2) ) 
  
  !-----------------------------------------------------------------------------
  ! HFB calculations
  !           n      <P>    RhoII  Delta , PairPartner E_sp d2H   <Jz> ,J,<r^2>
  3  format ( i3,1x,f5.2,1x,f7.4,1x,f7.2,1x, i3, 1x,f8.3,1x,e9.2,3(1x,f6.2)) 
  !           n      <P>     <S>    RhoII   Delta , PairPartner E_sp  d2H 
  31 format ( i3,1x,f5.2,1x,f5.2,1x,f7.4,1x,f7.2,1x,i3,1x,f8.3,1x,e9.2,        &
  &           3(1x,f6.2)) 
  !           <Jz>,J, <r^2>
  !           n      <P>     <S>    RhoII   Delta , PairPartner E_sp  d2H 
  32 format ( i3,1x,f5.2,1x,f5.2,1x,f7.4,1x,f7.2,1x,i3,1x,f8.3,1x,e9.2,        &
  &           4(1x,f6.2))
  !           <Jy>,<Jz>,J, <r^2>
  class(Spwf), intent(in) :: WF
  integer, intent(in)     :: i
  integer, intent(in)     :: PrintType
  real(KIND =dp)          :: PrintOcc
  
  PrintOcc = WF%Occupation
  ! Don't print double the occupation number when Time-reversal is conserved.
  if(TRC) PrintOcc = PrintOcc/2
  
  select case(PrintType)    
    case(1) !HF calculations
        !Time Reversal and signature conservation => <Rz> = 1
        if(TRC .and. SC .and. TSC) then
            print 1, i,WF%ParityR,PrintOcc,WF%Energy,WF%Dispersion,       &
            &          WF%AngMoment(3), WF%RMSRadius 
        elseif(SC.and.TSC) then
            print 11, i,WF%ParityR,WF%SignatureR,PrintOcc,WF%Energy,      &
            &           WF%Dispersion,WF%AngMoment(3),WF%AngQuantum,      &
            &           WF%RMSRadius
        elseif(TSC) then
        !No signature conservation
            print 12, i,WF%ParityR,WF%SignatureR,PrintOcc,WF%Energy,      &
            &           WF%Dispersion,WF%AngMoment(1),WF%AngMoment(3),    &
            &           WF%RMSRadius
        elseif(SC .and. .not. TSC) then
            print 13, i,WF%ParityR,WF%SignatureR,PrintOcc,WF%Energy,      &
            &           WF%Dispersion,WF%AngMoment(2),WF%AngMoment(3),    &
            &           WF%RMSRadius

        else
            call stp('no printout defined yet.')
        endif
        
    case(2) !BCS calculations
        if(TRC .and. SC) then
            print 2, i,WF%ParityR,PrintOcc,Real(WF%Delta), WF%Energy,     &
            &          WF%Dispersion,WF%AngMoment(3),WF%AngQuantum,       &
            &          WF%RMSRadius
        elseif(TRC) then
            print 21, i,WF%ParityR,WF%SignatureR,PrintOcc,Real(WF%Delta), &
            &          WF%Energy,WF%Dispersion,WF%AngMoment(1),           &
            &          WF%AngMoment(3),WF%AngQuantum, WF%RMSRadius
        else
            call stp('No Printout defined yet.')
        endif
        
    case(3) !HFB calculations
        if(TRC .and. SC .and. TSC) then
            print 3, i,WF%ParityR,PrintOcc,Real(WF%Delta),WF%PairPartner,     &
            &          WF%Energy,WF%Dispersion,WF%AngMoment(3),WF%AngQuantum, &
            &          WF%RMSRadius
        elseif(SC.and.TSC) then
            print 31,i,WF%ParityR,WF%SignatureR,PrintOcc,Real(WF%Delta),  &
            &          WF%PairPartner,WF%Energy,WF%Dispersion,WF%AngMoment(3), &
            &          WF%AngQuantum, WF%RMSRadius
        elseif(TSC) then
            print 32,i,WF%ParityR,WF%SignatureR,PrintOcc,Real(WF%Delta),       &
            &          WF%PairPartner,WF%Energy,WF%Dispersion,WF%AngMoment(1), &
            &          WF%AngMoment(3), WF%AngQuantum, WF%RMSRadius
        endif
    case DEFAULT
        call stp('Undefined printtype in subroutine PrintHF.')
  
  end select
  
  end subroutine PrintHF
  
  subroutine PrintCanonical(WF,i, PrintType)
  !-----------------------------------------------------------------------------
  ! Prints the info of a normal Spwf function in the Canonical basis, dependent 
  ! on the value of PrintType. The values and formats of PrintType should match
  ! with the ones from PrintHF.
  !-----------------------------------------------------------------------------
  ! PrintType
  ! = 3 => Minimal printout for HFB
  !        fmt =  n, <P> , (<S>) , v^2 , E , d^2h , (<Jx>), (<Jy>), <Jz>
  ! 
  !-----------------------------------------------------------------------------
  
  class(Spwf), intent(in) :: WF
  integer, intent(in)     :: i
  integer, intent(in)     :: PrintType
  
  !-----------------------------------------------------------------------------
  ! HFB calculations
  !           n      <P>    v^2      E_sp   <Jz> ,J,<r^2>
  3  format ( i3,1x,f5.2,1x,f7.4,1x,f8.3,3(1x,f6.2)) 
  !           n    <P>  <S>   v^2      E_sp   <Jz>,J, <r^2>
  31 format ( i3, 2(1x,f5.2) ,1x,f7.4,1x,f8.3,3(1x,f6.2)) 
  !           n    <P>  <S>   v^2      E_sp   <Jx> ,<Jz>,J, <r^2>
  32 format ( i3, 2(1x,f5.2) ,1x,f7.4,1x,f8.3,4(1x,f6.2) ) 
  real(KIND =dp)          :: PrintOcc
  
  PrintOcc = WF%Occupation
  ! Don't print double the occupation number when Time-reversal is conserved.
  if(TRC) PrintOcc = PrintOcc/2
  
  select case(PrintType)
  case(3)
    if(TRC.and.SC.and.TSC) then
        print 3, i , WF%ParityR, PrintOcc, WF%Energy, WF%AngMoment(3),    &
        &            WF%AngQuantum, WF%RMSRadius
    elseif(SC .and. TSC) then
        print 31, i , WF%ParityR, WF%SignatureR, PrintOcc, WF%Energy,     &
        &             WF%AngMoment(3), WF%AngQuantum, WF%RMSRadius
    elseif(TSC) then
        print 32, i , WF%ParityR, WF%SignatureR, PrintOcc, WF%Energy,     &
        &             WF%AngMoment(1),WF%AngMoment(3),  WF%AngQuantum,    &
        &             WF%RMSRadius
    endif
  case DEFAULT
    call stp('Undefined printtype in PrintCanonical.')
  end select
  
  end subroutine PrintCanonical
  
!-------------------------------------------------------------------------------

!---------------- Density-Computing Routines------------------------------------

  function GetTau(WaveFunction) result(tau)
    !---------------------------------------------------------------------------
    !Functions that returns the kinetic energy density of the spwf.
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WaveFunction
    real(KIND=dp)           :: tau(nx,ny,nz)
    integer                 :: i

    tau = 0.0_dp
    do i=1,3
      tau = tau + RealMultiplySpinor(WaveFunction%Der(i), WaveFunction%Der(i))
    enddo
  end function GetTau
        
  pure function GetDerRho(WF) result(DerRho)
    !---------------------------------------------------------------------------
    ! This function returns the contribution of this wavefunction to 
    ! the derivative of rho.
    !       \partial_{k} \rho = 2 \sum_{i} \Psi_i \partial_k \Psi_i
    !---------------------------------------------------------------------------
    !NOTE:
    ! This is a non-active subroutine, it is never called.
    !---------------------------------------------------------------------------
    class(Spwf),intent(in) :: WF
    real(KIND=dp)          :: DerRho(nx,ny,nz,3)
    integer                :: i

    do i=1,3
            DerRho(:,:,:,i) = RealMultiplySpinor(WF%Value , WF%Der(i))
    enddo
    DerRho = DerRho*2.0_dp
  end function GetDerRho
        
  function GetLapRho(WaveFunction) result(LapRho)
    !---------------------------------------------------------------------------
    ! Function that returns the value of the Laplacian of the density.
    ! \Delta \rho = 2*\sum_{\sigma} Re[
    !     \Delta \phi^*(\sigma) ]\phi(\sigma) 
    !   + \grad \phi^*(\sigma) \grad \phi(\sigma)
    !---------------------------------------------------------------------------
    !NOTE:
    ! This is a non-active subroutine, it is never called.
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: WaveFunction
    real(KIND=dp)           :: LapRho(nx,ny,nz)
    integer                 :: j

    LapRho = 2.0_dp*RealMultiplySpinor(WaveFunction%Lap,WaveFunction%Value)
    do j=1,3
     LapRho = LapRho +                                                         &
     &        2.0_dp*RealMultiplySpinor(WaveFunction%Der(j),WaveFunction%Der(j))             
    enddo
  end function GetLapRho

  function GetVecj(WaveFunction) result(Vecj)
    !---------------------------------------------------------------------------
    !Functions that returns the \vec{j} density of the spwf.
    ! \vec{j} = Im(\sum_{\sigma} \Psi^*(\sigma) \nabla \Psi(\sigma))
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WaveFunction
    real(KIND=dp)           :: Vecj(nx,ny,nz,3)
    integer                 :: m

    Vecj=0.0_dp
    ! vec(j)_{mu} = Im(\sum_{sigma} \Psi^*\partial_{mu}(\psi) )
    do m=1,3
      Vecj(:,:,:,m)= ImagMultiplySpinor(WaveFunction%Value, WaveFunction%Der(m))                     
    enddo
  end function GetVecj
        
  function GetRotVecJ(WF) result (RotVecJ)
    !---------------------------------------------------------------------------
    ! Function that returns the curl of the \vec{j} density.
    ! \nabla x \vec{j} = Im( \nabla x (\Psi^* \nabla \Psi))
    ! Note that no second derivatives figure since the curl of a gradient is 0.
    !---------------------------------------------------------------------------
    !NOTE:
    ! This is a non-active subroutine, it is never called.              
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WF
    real(KIND=dp)           :: RotVecj(nx,ny,nz,3)

    RotVecJ = 0.0_dp
    !X Component     
    RotVecJ(:,:,:,1) =  ImagMultiplySpinor(WF%Der(2) , WF%Der(3))
    !Y Component
    RotVecJ(:,:,:,2) =  ImagMultiplySpinor(WF%Der(3) , WF%Der(1))
    !Z Component
    RotVecJ(:,:,:,3) =  ImagMultiplySpinor(WF%Der(1) , WF%Der(2)) 
  end function GetRotVecJ
        
  function GetVecs(WaveFunction) result(Vecs)
    !---------------------------------------------------------------------------
    ! Function that returns the spin density \vec{s}
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WaveFunction
    real(KIND=dp)           :: Vecs(nx,ny,nz,3)
    integer                 :: l
    type(Spinor)            :: Psi

    Vecs=0.0_dp
    do l=1,3
            Psi = Pauli(WaveFunction%Value,l)
            Vecs(:,:,:,l) = RealMultiplySpinor(WaveFunction%Value,Psi)   
    enddo
  end function GetVecs
  
  function GetVecT(WF) result(VecT)
    !---------------------------------------------------------------------------
    ! Function that returns the density \vec{T}
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WF
    real(KIND=dp)           :: VecT(nx,ny,nz,3), Temp(nx,ny,nz)
    integer                 :: l
    type(Spinor)            :: Psi

    VecT(:,:,:,:)=0.0_dp

    if(WF%Signature.eq.0)  then
      !-------------------------------------------------------------------------
      !X component
      !T_x = 2*Re[\nabla \Psi^*(r,+)\cdot\nabla \Psi(r,-) ]
      ! This is constrained by signature symmetry.
      !-------------------------------------------------------------------------
      Temp =0.0_dp
      do l=1,3
        Psi  = Pauli(WF%Der(l),1)
        Temp = Temp + RealMultiplySPinor(WF%Der(l) , Psi)
      enddo
      VecT(:,:,:,1)  = Temp

    elseif((WF%Signature.eq.0).and.(WF%TimeSimplex.eq.0))then
      !-------------------------------------------------------------------------
      !Y component
      !T_y = 2*Im[\nabla \Psi^*(r,+)\cdot\nabla \Psi(r,-) ]
      !This is constrained by both time simplex and signature.
      !-------------------------------------------------------------------------
      Temp =0.0_dp
      do l=1,3
        Psi  = Pauli(WF%Der(l),2)
        Temp = Temp + RealMultiplySpinor(WF%Der(l) , Psi)
      enddo
      VecT(:,:,:,2)  = Temp
    endif
    !---------------------------------------------------------------------------
    !Z component
    ! T_z = \nabla \Psi^*(r,+)\cdot\nabla \Psi(r,+)
    !     - \nabla \Psi^*(r,-)\cdot\nabla \Psi(r,-)
    !---------------------------------------------------------------------------
    Temp =0.0_dp
    do l=1,3
      Psi  = Pauli(WF%Der(l),3)
      Temp = Temp + RealMultiplySpinor(WF%Der(l) , Psi)
    enddo
    VecT(:,:,:,3)  = Temp
  end function GetVecT
        
  function GetVecF(WaveFunction) result(VecF)
    !-------------------------------------------------------
    ! Function that returns the density \vec{F}
    !-------------------------------------------------------
    class (Spwf),intent(in) :: WaveFunction
    real(KIND=dp)           :: VecF(nx,ny,nz,3), Temp(nx,ny,nz)
    integer                 :: mu,nu
    type(Spinor)            :: Psi

    VecF=0.0_dp
    do mu=1,3  
      Temp = 0.0_dp
      do nu=1,3
        Psi  = Pauli(WaveFunction%Der(mu),nu) 
        Temp = Temp + RealMultiplySpinor(WaveFunction%Der(nu) , Psi)
      enddo
      VecF(:,:,:,mu) = Temp                            
    enddo
    return            
  end function GetVecF

  function GetJMuNu(WaveFunction) result(JMuNu)
    !---------------------------------------------------------------------------
    ! Function that returns the tensorial density J_{\mu\nu}
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WaveFunction
    real(KIND=dp)           :: JMuNu(nx,ny,nz,3,3)
    integer                 :: mu,nu
    type(Spinor)            :: P

    JMuNu=0.0_dp
    do nu=1,3
      do mu=1,3
        P = Pauli(WaveFunction%Der(mu),nu)                
        ! We take the imaginary part of the product of the spinors, 
        ! instead of multiplying by -I/2.
        JMuNu(:,:,:,mu,nu) =  ImagMultiplySpinor(WaveFunction%Value, P)
      enddo
    enddo
  end function GetJMuNu
        
  function GetNablaJ(WF) result(NablaJ)
    !---------------------------------------------------------------------------
    !Functions that returns \Nabla * J
    !---------------------------------------------------------------------------
    class (Spwf),intent(in) :: WF
    real(KIND=dp)           :: NablaJ(nx,ny,nz)
    integer                 :: i,j,k
    type(Spinor)            :: Psi

    NablaJ=0.0_dp
    
    ! This is very ugly, but I did not immediately find a way to write this more 
    ! elegantly without needing A LOT OF extra operations...
    i=1;j=2;k=3
    Psi  = Pauli(WF%Der(j),k)
    NablaJ = NablaJ + dble(LeviCivita(i,j,k))*ImagMultiplySpinor(WF%Der(i),Psi)

    i=1;j=3;k=2
    Psi  = Pauli(WF%Der(j),k)
    NablaJ = NablaJ + dble(LeviCivita(i,j,k))*ImagMultiplySpinor(WF%Der(i),Psi)

    i=3;j=2;k=1
    Psi  = Pauli(WF%Der(j),k)
    NablaJ = NablaJ + dble(LeviCivita(i,j,k))*ImagMultiplySpinor(WF%Der(i),Psi)

    i=3;j=1;k=2
    Psi  = Pauli(WF%Der(j),k)
    NablaJ = NablaJ + dble(LeviCivita(i,j,k))*ImagMultiplySpinor(WF%Der(i),Psi)

    i=2;j=1;k=3
    Psi  = Pauli(WF%Der(j),k)
    NablaJ = NablaJ + dble(LeviCivita(i,j,k))*ImagMultiplySpinor(WF%Der(i),Psi)

    i=2;j=3;k=1
    Psi  = Pauli(WF%Der(j),k)
    NablaJ = NablaJ + dble(LeviCivita(i,j,k))*ImagMultiplySpinor(WF%Der(i),Psi)
  end function GetNablaJ
        
!--------------- In- and Output ------------------------------------------------
  subroutine WriteSpwf(wf, unit)
    !---------------------------------------------------------------------------
    ! Custom subroutine for writing Spwfs to file. 
    ! This is necessary to keep the components of an Spwf private: the standard 
    ! write doesn't work with them.
    !---------------------------------------------------------------------------
    class(Spwf), intent(in) :: wf
    integer, intent(in)     :: unit

    call WriteSpinor(wf%Value,unit)
    !It is of no use to store the derivatives, unless asking for all the output

    write(unit) wf%Occupation, &
    &       wf%Energy, wf%Dispersion, wf%Norm, wf%AngMoment, &
    &       wf%Isospin,wf%TimeSimplex,wf%Parity,wf%Signature,wf%TimeReversal
       
  end subroutine WriteSpwf
  subroutine ReadSpwf(wf, unit, filenx,fileny,filenz)
    !---------------------------------------------------------------------------
    ! Custom subroutine for reading Spwfs from file. 
    ! This is necessary to keep the components of an Spwf private: the standard 
    ! read doesn't work with them.
    !---------------------------------------------------------------------------
    class(Spwf), intent(inout) :: wf
    integer, intent(in)        :: unit,filenx,fileny,filenz
    integer                    :: i

    call ReadSpinor(wf%Value,unit,filenx,fileny,filenz,1)
    !Allocating the derivatives
    do i=1,3
        wf%Der(i)=NewSpinor(); 
    enddo
    wf%Lap = NewSpinor()
    Read(unit) wf%Occupation, &
    &       wf%Energy, wf%Dispersion, wf%Norm, wf%AngMoment, &
    &       wf%Isospin,wf%TimeSimplex,wf%Parity,wf%Signature,wf%TimeReversal
    
    if(abs(wf%Isospin).ne.1) then
      call stp('Wrong isospin on file.')
    endif
    
  end subroutine ReadSpwf
        
  subroutine SymmetryOperators(wf)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the expected values of the symmetry 
    ! operators for the wavefunction.
    ! In addition, the RMS radius of the wavefunction is computed.
    !---------------------------------------------------------------------------
    class(Spwf), intent(inout) :: wf
    type(Spinor)               :: Temp

    if(wf%Isospin.ne. 0) wf%IsospinR = wf%GetIsospin()

    !Expectation value of Signature
    if(wf%Signature.ne. 0) then
      wf%SignatureR = wf%Signature
    else
      Temp = ActionOfSignature(wf%Value)
     !Note that <Psi | Rz | Psi> is no Hermitian
      wf%SignatureR = InproductSpinorImaginary(wf%Value, Temp)
    endif

    !Expectation value of Parity
    if(wf%Parity .ne. 0) then
      wf%ParityR = wf%Parity
    else
      Temp = ActionOfParity(wf%value, wf%Signature)
      wf%ParityR = InproductSpinorReal(wf%Value, Temp)
    endif

    !Expectation value of Time Simplex.
    ! Note that it will be strictly real...
    if(wf%TimeSimplex .ne. 0) then
      wf%TimeSimplexR = wf%TimeSimplex
    else
      Temp = ActionOfTimeSimplex(wf%value)
      wf%TimeSimplexR = InproductSpinorReal(wf%value, Temp)
    endif

    !Calculate the RMS radius
    call wf%CalcRMSRadius()

  end subroutine SymmetryOperators
  
  subroutine CalcRMSRadius(wf)
    !---------------------------------------------------------------------------
    ! Calculates and stores sqrt(<R^2>) for every wavefunction.
    !---------------------------------------------------------------------------
    use Mesh
    class(Spwf), intent(inout) :: wf

    wf%RMSradius = sqrt(sum(sum(Mesh3D**2,1) * wf%GetDensity())*dv)
  
  end subroutine CalcRMSRadius
end module WaveFunctions
