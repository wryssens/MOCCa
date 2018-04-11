    module Damping

  use Compilationinfo
  use Geninfo
  use Force
  use Spinors
  
  implicit none
 
  !---------------------------------------------------------------------------
  ! Effective density mixing factor, for pringting purposes.
  real(KIND=dp)              :: preconfac = 0.0
  !---------------------------------------------------------------------------
  ! Auxiliary arrays. 
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
 
!===============================================================================
!===============================================================================

  !-----------------------------------------------------------------------------
  function PreconditionRho(drho, alpha, beta) result(invrho)
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    use Derivatives

    real*8 :: drho(nx,ny,nz,2), residual(nx,ny,nz), update(nx,ny,nz), aCG
    real*8 :: invrho(nx,ny,nz,2), direction(nx,ny,nz), bCG, amix
    real*8 :: newresnorm, oldresnorm, ar(nx,ny,nz), br(nx,ny,nz)
    
    real*8, intent(in) :: alpha(2), beta(2)
    
    integer:: it, iter, p, s, ts
   
    p = 0
    if(PC)  p  = 1
    s = 0
    if(SC)  s  = 1
    ts = 0
    if(TSC) ts = 1
   
    amix = preconfac

    !---------------------------------------------------------------------------
    invrho           = 0.0
    do it=1,2 
        Residual         = drho(:,:,:,it) 
        Direction        = Residual
        newresnorm       = sum(direction**2)*dv

        ar = alpha(it)
        br = beta(it)
        
        do iter=1,300
          update   = preconoperator(direction, ar, br,p,s,ts,1)

          aCG    = NewResNorm/(sum(Direction*update)*dv)
          
          invrho(:,:,:,it)   = invrho(:,:,:,it)  + aCG * Direction
          residual           = residual          - aCG * update
          
          oldresnorm = newresnorm
          newresnorm = sum(residual**2)*dv
          
          BCG      = NewResNorm/OldResNorm
          Direction  = Residual + bCG * Direction
          !print *, it, newresnorm
          if(newresnorm.lt.1d-8) exit
        enddo
    enddo
    !---------------------------------------------------------------------------
    print *, 'Inverted, iter = ', iter, newresnorm, sum(invrho(:,:,:,1))*dv,   &
    &                     sum(invrho(:,:,:,2))*dv
    
  end function PreconditionRho
  
    !-----------------------------------------------------------------------------
  function PreconditionPotential(pot,alpha, beta, p,s,ts,c) result(invpot)
    !---------------------------------------------------------------------------
    ! Precondition a potential instead of the densities. 
    !---------------------------------------------------------------------------
    use Derivatives

    real*8 :: pot(nx,ny,nz,2), residual(nx,ny,nz), update(nx,ny,nz), aCG
    real*8 :: invpot(nx,ny,nz,2), direction(nx,ny,nz), bCG
    real*8 :: newresnorm, oldresnorm
    real*8 :: alpha(nx,ny,nz,2), beta(nx,ny,nz,2)
    
    integer, intent(in) :: p,s,ts,c
    
    integer:: it, iter
    
    !---------------------------------------------------------------------------
    invpot           = 0.0
    do it=1,2 
        Residual         = pot(:,:,:,it) 
        Direction        = Residual
        newresnorm       = sum(direction**2)*dv

        !-----------------------------------------------------------------------
        do iter=1,300
          update  &
          &  = preconoperator(direction,alpha(:,:,:,it),beta(:,:,:,it),p,s,ts,c)

          aCG    = NewResNorm/(sum(Direction*update)*dv)
          
          invpot(:,:,:,it)   = invpot(:,:,:,it)  + aCG * Direction
          residual           = residual          - aCG * update
          
          oldresnorm = newresnorm
          newresnorm = sum(residual**2)*dv
          
          BCG      = NewResNorm/OldResNorm
          Direction  = Residual + bCG * Direction
          if(newresnorm.lt.1d-8) exit
        enddo
        !-----------------------------------------------------------------------
    enddo
    !---------------------------------------------------------------------------
!    print *, 'Inv. Pot., iter = ', iter, newresnorm, sum(invpot(:,:,:,1))*dv,  &
!    &                     sum(invpot(:,:,:,2))*dv
    
  end function Preconditionpotential
  
  function preconoperator(drho,alpha,beta,p,s,ts,c) result(Prho)
    !---------------------------------------------------------------------------
    ! Implement the preconditioning operator
    !
    ! Pf(r) = B(r)*f(r) + a(r) * Delta(f(r)) 
    ! 
    ! There is at the moment no routine that makes use of the (r) dependence
    ! of neither B, nor A, but future improvements might depend on it. 
    !---------------------------------------------------------------------------
    use Force
    use Derivatives
    
    real(KIND=dp), intent(in) :: drho(nx,ny,nz)
    real(KIND=dp)             :: Prho(nx,ny,nz)
    integer                   :: p, s, ts,c, it
    real(KIND=dp),intent(in)  :: alpha(nx,ny,nz), beta(nx,ny,nz)
    !---------------------------------------------------------------------------
    ! Symmetries of the problem
    !---------------------------------------------------------------------------

    Prho = Laplacian(drho, p,s,ts,c)
    Prho = beta*drho + alpha*Prho
    
  end function preconoperator
end module Damping

