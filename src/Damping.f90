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
  
  !-----------------------------------------------------------------------------
  function InverseRho(drho, rho) result(invrho)
    !---------------------------------------------------------------------------
    ! Precondition the difference in densities, in order to penalize the 
    ! highly oscillatory components.
    !
    !
    !---------------------------------------------------------------------------
    use Derivatives

    real*8 :: drho(nx,ny,nz,2), residual(nx,ny,nz), update(nx,ny,nz), alpha
    real*8 :: invrho(nx,ny,nz,2), direction(nx,ny,nz), beta, amix
    real*8 :: newresnorm, oldresnorm, rho(nx,ny,nz,2)
    real*8 :: bfactor(nx,ny,nz), meff
    integer:: it, iter, p, s, ts
   
    amix = preconfac
    !---------------------------------------------------------------------------
    ! Symmetries of the problem
    !---------------------------------------------------------------------------
    if(PC) then
        p = 1
    else
        p = 0
    endif
    if(TSC) then
        ts = 1
    else
        ts = 0
    endif
    if(SC) then
        s = 1
    else
        s = 0
    endif
    !---------------------------------------------------------------------------
    do it=1,2
          invrho(:,:,:,it) = 0.0
          Residual         = drho(:,:,:,it)
          Direction        = Residual
          newresnorm       = sum(direction**2)*dv
          
          do iter=1,500
              update   = (Direction -amix*Laplacian(Direction, p,s,ts,+1))
              

              alpha   = NewResNorm/(sum(Direction*update)*dv)
              
              invrho(:,:,:,it) = invrho(:,:,:,it) + alpha * Direction
              residual         = residual         - alpha * update
              
              oldresnorm = newresnorm
              newresnorm = sum(residual**2)*dv
              
              Beta       = NewResNorm/OldResNorm
              Direction  = Residual + beta * Direction

              if(newresnorm.lt.1d-9) exit
          enddo
    enddo
    !---------------------------------------------------------------------------
  end function InverseRho
  
!===============================================================================
!===============================================================================

  !-----------------------------------------------------------------------------
  function PreconditionRho(drho, alpha, beta) result(invrho)
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    use Derivatives

    real*8 :: drho(nx,ny,nz,2), residual(nx,ny,nz,2), update(nx,ny,nz,2), aCG
    real*8 :: invrho(nx,ny,nz,2), direction(nx,ny,nz,2), bCG, amix
    real*8 :: newresnorm, oldresnorm
    
    real*8, intent(in) :: alpha, beta
    integer:: it, iter
   
    amix = preconfac

    !---------------------------------------------------------------------------
    invrho           = 0.0
    Residual         = drho
    Direction        = Residual
    newresnorm       = sum(direction**2)*dv
          
    do iter=1,100
      update   = preconoperator(direction, alpha, beta)

      aCG    = NewResNorm/(sum(Direction*update)*dv)
      
      invrho   = invrho   + aCG * Direction
      residual = residual - aCG * update
      
      oldresnorm = newresnorm
      newresnorm = sum(residual**2)*dv
      
      BCG      = NewResNorm/OldResNorm
      Direction  = Residual + bCG * Direction
      print *, newresnorm
      if(newresnorm.lt.1d-8) exit
    enddo
    !---------------------------------------------------------------------------
    print *, 'Inverted, iter = ', iter, newresnorm, sum(invrho(:,:,:,1))*dv,   &
    &                     sum(invrho(:,:,:,2))*dv
    print *, 'Deviation', sum((preconoperator(invrho,alpha, beta) - drho)**2)*dv
    
  end function PreconditionRho
  
  function preconoperator(drho, alpha, beta) result(Prho)
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    use Force
    use Derivatives
    
    real(KIND=dp), intent(in) :: drho(nx,ny,nz,2)
    real(KIND=dp)             :: Prho(nx,ny,nz,2)
    integer                   :: p, s, ts, it
    real(KIND=dp)             :: alpha, beta
    !---------------------------------------------------------------------------
    ! Symmetries of the problem
    !---------------------------------------------------------------------------
    if(PC) then
        p = 1
    else
        p = 0
    endif
    if(TSC) then
        ts = 1
    else
        ts = 0
    endif
    if(SC) then
        s = 1
    else
        s = 0
    endif
    
    do it=1,2
        Prho(:,:,:,it) = Laplacian(drho(:,:,:,it), p,s,ts,+1)
    enddo
    
    Prho = beta*drho - alpha*Prho
  end function preconoperator
end module Damping

!===============================================================================
! Code ZOO
!===============================================================================
!  !-----------------------------------------------------------------------------
!  
!  function AverageSpinor(Psi, P,S, TS, I) result(Phi)
!    !---------------------------------------------------------------------------
!    ! Averages the spinor Psi along every spatial direction, weighted with a 
!    ! gaussian exp(-a Delta x^2).
!    !---------------------------------------------------------------------------
!    use Spinors
!    use Derivatives
!    integer, intent(in)       :: P,S, TS, I
!    type(Spinor)              :: Psi
!    type(Spinor)              :: Phi
!    real(KIND=dp)             :: a
!    integer                   :: j
!    
!    !---------------------------------------------------------------------------
!    !Change the derivative coefficients
!    deallocate(FDLap)
!    allocate(FDLap(-2:2))
!    !---------------------------------------------------------------------------
!    ! Numerical parameter. Note that this is adjusted semi-empirically. The 
!    ! lower, the faster convergence, but is is bounded from below. The 
!    ! approximation to the exponential should be positive definite.
!    !  See P.G. Reinhardt & Cusson, Nuc. Phys. A A378, 418-442
!    !---------------------------------------------------------------------------
!    a = 2.8_dp

!    !Factor dx**2 du to laplacian routines dividing by dx**2
!    do j=1,size(FDLap)/2
!      FDLap(-j) = exp(-j**2*a*dx**2)*dx**2 *sqrt(a/pi)
!      FDLap(j ) = FDLap(-j)
!    enddo
!    !Factor 3.0 due to 3D.
!    FDLap(0)  = dx**2/3.0_dp  *sqrt(a/pi) 

!    Phi = LapSpinor(Psi,P,S,TS)
!    !-----------------------------------------------------------------------------
!    ! Resetting the derivation routines
!    call AssignFDCoefs(MaxFDOrder, MaxFDLapOrder, CoulombLapOrder)
!  
!  end function AverageSpinor
!  
!  function AverageDensity(r) result(ar)
!    real(KIND=dp)             :: ar(nx,ny,nz,2), r(nx,ny,nz,2), a
!    integer                   :: it, i,j,k
!    
!    !---------------------------------------------------------------------------
!    ! Numerical parameter. Note that this is adjusted semi-empirically. The 
!    ! lower, the faster convergence, but is is bounded from below. The 
!    ! approximation to the exponential should be positive definite.
!    !  See P.G. Reinhardt & Cusson, Nuc. Phys. A A378, 418-442
!    !---------------------------------------------------------------------------
!    a = 0.25_dp
!    
!    do it=1,2
!       do k=1,nz
!        do j=1,ny 
!         do i=2,nx-1
!            ar(i,j,k,it) = (1-2*a)*r(i,j,k,it) + a*r(i+1,j,k,it)        &
!            &                                  + a*r(i-1,j,k,it)
!         enddo
!         ar(1,j,k,it)  = (1-2*a)*r(1,j,k,it) + a*r(2,j,k,it)             &
!         &                                  + a*r(1,j,k,it)
!         ar(nx,j,k,it) = (1-2*a)*r(nx,j,k,it) + a*r(nx-1,j,k,it) 
!        enddo
!       enddo
!       
!       !------------------------------------------------------------------------
!       do k=1,nz
!        do i=1,nx 
!         do j=2,ny-1
!            ar(i,j,k,it) = (1-2*a)*r(i,j,k,it) + a*r(i,j+1,k,it) &
!            &                                  + a*r(i,j-1,k,it)
!         enddo
!         ar(i,1,k,it)  = (1-2*a)*r(i,1,k,it)  + a*r(i,2,k,it)    &
!         &                                    + a*r(i,1,k,it)
!         ar(i,ny,k,it) = (1-2*a)*r(i,ny,k,it) + a*r(i,ny-1,k,it)
!        enddo
!       enddo
!       
!       !------------------------------------------------------------------------
!       do j=1,ny
!        do i=1,nx 
!         do k=2,nz-1
!            ar(i,j,k,it) = (1-2*a)*r(i,j,k,it) + a*r(i,j,k+1,it) &
!            &                                  + a*r(i,j,k-1,it)
!         enddo
!         ar(i,j,1,it)  = (1-2*a)*r(i,j,1,it)  + a*r(i,j,2,it)    &
!         &                                    + a*r(i,j,1,it)
!         ar(i,j,nz,it) = (1-2*a)*r(i,j,nz,it) + a*r(i,j,nz,it)
!        enddo
!       enddo
!       
!    enddo

!  end function AverageDensity
