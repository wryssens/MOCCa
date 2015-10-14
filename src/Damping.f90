module Damping

  use Compilationinfo
  use Geninfo
  use Force
  use Spinors
  
  implicit none
  
  real(KIND=dp), allocatable :: ExpLapCoefX(:,:), ExpLapCoefY(:,:)
  real(KIND=dp), allocatable :: ExpLapCoefZ(:,:)
  real(KIND=dp)              :: ExpLapDiag(3,2)
 
  integer                    :: IterationCount=0
    
contains

  function InverseKinetic(Psi, E,P,S, TS, I) result(Phi)
  !-----------------------------------------------------------------------------
  ! Use a conjugate gradient algorithm to solve the following equation for phi:
  ! 
  ! (P^2 + E0) Phi = Psi 
  !
  ! But in practice the equation
  ! [\Delta + E0/(-hbm/2))] Phi = 1/(-hbm/2) Psi
  ! is better conditioned.
  !
  !-----------------------------------------------------------------------------
  use Spinors, only     : Lapspinor
  use Force
  use Derivatives, only : AssignFDCoefs,MaxFDOrder,                            &
  &                       MaxFDLapOrder, CoulombLapOrder
  !use MeanFields

  integer                   :: iter, it, m
  type(Spinor)              :: Psi, Phi
  real(KIND=dp), intent(in) :: E
  type(Spinor)              :: Residual, Direction, Update, Der(3)
  real(KIND=dp)             :: NewResNorm, OldResNorm, alpha, beta, Prec=1.0d-10
  integer, intent(in)       :: P,S, TS, I
  
  !-----------------------------------------------------------------------------
  ! Use first-order laplacians here for speed.
  call AssignFDCoefs(1, 1, CoulombLapOrder)
  
  Residual=NewSpinor() ; Direction=NewSpinor()
  Update=NewSpinor()   ; Phi =NewSpinor()
  it = (I + 3) /2
  !-----------------------------------------------------------------------------
  ! Step 0
  ! Note that we take as initial guess zero everywhere. One might think that 
  ! taking Psi itself as initial guess is better, but in for some unknown 
  ! reason it usually is slightly quicker to take zero as initial guess.
  ! Something like 18 iterations (for zero-guess) versus 25 (for Psi-guess).
  Residual   = (1.0_dp/(-hbm(it)/2.0_dp))*Psi
  Direction  = Residual
  NewResNorm = InproductSpinorReal(Direction, Direction)
  !-----------------------------------------------------------------------------
  !Starting the iterations
  do iter=1,100
    !---------------------------------------------------------------------------
    ! Calculating the update
    Update = LapSpinor(Direction,P,S,TS)                                       &
    &      +(E )/(-hbm(it)/(2.0_dp))*Direction
    !Update = (Bpot(:,:,:,it)/(hbm(it)/(2.0_dp))) * LapSpinor(Direction,P,S,TS) &
    !  &      -(E )/(hbm(it)/(2.0_dp))*Direction

!      Der = DeriveSpinor(Direction,P,S,TS)
!      do m=1,3
!        Update = Update + NablaBPot(:,:,:,m,it)/(hbm(it)/(2.0_dp))*Der(m)
!      enddo


    alpha  = NewResNorm/(InproductSpinorReal(Direction, Update))
    !---------------------------------------------------------------------------
    ! Applying the update
    Phi        = Phi      + alpha * Direction
    Residual   = Residual - alpha * Update
    OldResNorm = NewResNorm
    NewResNorm = InproductSpinorReal(Residual, Residual)
    !---------------------------------------------------------------------------
    !Finding the next search direction
    Beta       = NewResNorm/OldResNorm
    Direction  = Residual + beta * Direction 
    
    if(NewResNorm.lt. Prec) then
      exit
    endif
  enddo
  if(NewResNorm .gt. Prec) then
    call stp("CG in the damping step didn't converge!",                        &
    'Residual Norm', NewResNorm)
  endif
  IterationCount = IterationCount + iter
  !-----------------------------------------------------------------------------
  ! Resetting the derivation routines
  call AssignFDCoefs(MaxFDOrder, MaxFDLapOrder, CoulombLapOrder)
  
  end function InverseKinetic
  
  function AverageSpinor(Psi, P,S, TS, I) result(Phi)
    !---------------------------------------------------------------------------
    ! Averages the spinor Psi along every spatial direction, weighted with a 
    ! gaussian exp(-a Delta x^2).
    !---------------------------------------------------------------------------
    use Spinors
    use Derivatives
    integer, intent(in)       :: P,S, TS, I
    type(Spinor)              :: Psi
    type(Spinor)              :: Phi
    real(KIND=dp)             :: a
    integer                   :: j
    
    !---------------------------------------------------------------------------
    !Change the derivative coefficients
    deallocate(FDLap)
    allocate(FDLap(-2:2))
    !---------------------------------------------------------------------------
    ! Numerical parameter. Note that this is adjusted semi-empirically. The 
    ! lower, the faster convergence, but is is bounded from below. The 
    ! approximation to the exponential should be positive definite.
    !  See P.G. Reinhardt & Cusson, Nuc. Phys. A A378, 418-442
    !---------------------------------------------------------------------------
    a = 2.8_dp

    !Factor dx**2 du to laplacian routines dividing by dx**2
    do j=1,size(FDLap)/2
      FDLap(-j) = exp(-j**2*a*dx**2)*dx**2 *sqrt(a/pi)
      FDLap(j ) = FDLap(-j)
    enddo
    !Factor 3.0 due to 3D.
    FDLap(0)  = dx**2/3.0_dp  *sqrt(a/pi) 

    Phi = LapSpinor(Psi,P,S,TS)
    !-----------------------------------------------------------------------------
    ! Resetting the derivation routines
    call AssignFDCoefs(MaxFDOrder, MaxFDLapOrder, CoulombLapOrder)
  
  end function AverageSpinor
  
end module Damping
