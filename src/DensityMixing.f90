module DensityMixing
!-------------------------------------------------------------------------------
! Module that contains the routines for mixing the densities.
! Why a separate module: because the Energy module needs the densities
! and some strategies need the Energy module.
!-------------------------------------------------------------------------------
! General note concerning most (if not all) of these mixing strategies.
! The element 1 in the array DensityHistory is NEVER used to mix in advanced
! strategies. This is because it represents the Rho_in to the Rho_out currently
! in the Density densityvector. Since we should only mix Rho_out's, these
! never get included. This is always visible in the iteration indices.
!-------------------------------------------------------------------------------

  use CompilationInfo
  use Geninfo
  use Densities

  implicit none

contains

  subroutine MixDensities(Iteration)
    !---------------------------------------------------------------------------
    ! Mixes the densities according to different schemes.
    ! MixingScheme=
    ! 0) => Linear Mixing with DampingParam as parameter.
    ! 1) => Linear mixing of ONLY rho.
    ! 2) => Linear mixing of ONLY laplacian of rho and laplacian of S.
    ! 3) => 
    ! 4) => DIIS
    !---------------------------------------------------------------------------
    use Damping
    
    integer, intent(in)         :: Iteration
    real(KIND=dp), allocatable  :: temp(:,:,:,:), last(:,:,:,:)
    real(KIND=dp)               :: preconrho(nx,ny,nz,2)

    DensityChange = 1 - Density * DensityHistory(1)/(Density * Density)
    select case(MixingScheme)
      case(0)
        !-----------------------------------------------------------------------
        ! This one-liner is more complex than it seems: it mixes ALL of the
        ! densities. 
        Density = (1-DampingParam) * Density + DampingParam*DensityHistory(1)
      case(1)
        !-----------------------------------------------------------------------
        ! Only mix rho with estimated dampingparam.
        DampingParam = sum(abs(Cdrho))
        DampingParam = (6*dt/hbar/(dx**2)*DampingParam)*1.5
        DampingParam = 1 - 1/DampingParam

        print *, 'damping', dampingparam
        
        Density%rho = (1-DampingParam)*Density%rho + &
        &                DampingParam *DensityHistory(1)%rho
        
        ! Recalculate the derivatives of rho
        call RecalcRhoDerivatives()
      case(2)
        !-----------------------------------------------------------------------
        ! only mix laprho and laps
        Density%laprho = (1-DampingParam)*Density%laprho + &
        &                   DampingParam *DensityHistory(1)%laprho
      case(3)
        !-----------------------------------------------------------------------
        preconrho = Inverserho(Density%rho - DensityHistory(1)%rho, DensityHistory(1)%rho)
           
        Density%rho = DensityHistory(1)%rho + preconrho
    
        if(any(Density%rho .lt. 0)) then
            print *, minval(Density%rho)
            where(Density%rho.lt.0) Density%rho = 0
        endif   
        
        ! Recalculate the derivatives of rho
        call RecalcRhoDerivatives()        
      case(4)
        !-----------------------------------------------------------------------
        call DIIS(mod(Iteration,100))
      case DEFAULT
        call stp('Mixing scheme not supported!')
    end select

  end subroutine MixDensities

  subroutine RecalcRhoDerivatives
    !---------------------------------------------------------------------------
    ! Small Subroutine to recalculate the derivatives of the density rho, 
    ! for use after this density has been mixed.
    !
    !---------------------------------------------------------------------------
    use Force
    
    integer :: it, P, S, TS
    
    P=0         ; S=0          ; TS=0
    if(PC) P=1  ; if(SC) S=1  ; if(TSC) TS = 1
            
    Density%laprho(:,:,:,1) = Laplacian(Density%rho(:,:,:,1),P,S,TS,1)
    Density%laprho(:,:,:,2) = Laplacian(Density%rho(:,:,:,2),P,S,TS,1)
    
    do it=1,2
        Density%DerRho(:,:,:,1,it) = DeriveX(Density%Rho(:,:,:,it), P,S,TS,1)
        Density%DerRho(:,:,:,2,it) = DeriveY(Density%Rho(:,:,:,it), P,S,TS,1)
        Density%DerRho(:,:,:,3,it) = DeriveZ(Density%Rho(:,:,:,it), P,S,TS,1)
        Density%LapRho(:,:,:,it)   = Laplacian(Density%Rho(:,:,:,it),P,S,TS,1)
    enddo
    
    if(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp) then
        do it=1,2
            Density%LapLapRho(:,:,:,it)   = &
            & Laplacian(Density%LapRho(:,:,:,it),P,S,TS,1)
        
            Density%D2Rho(:,:,:,1,1,it) = &
            & DeriveX(Density%DerRho(:,:,:,1,it), -P,-S, TS,1)
            Density%D2Rho(:,:,:,2,1,it) = & 
            & DeriveY(Density%DerRho(:,:,:,1,it), -P,-S, TS,1)
            Density%D2Rho(:,:,:,3,1,it) = &
            & DeriveZ(Density%DerRho(:,:,:,1,it), -P,-S, TS,1)

            Density%D2Rho(:,:,:,1,2,it) = &
            & DeriveX(Density%DerRho(:,:,:,2,it), -P,-S, TS,2)
            
            Density%D2Rho(:,:,:,2,2,it) = &
            & DeriveY(Density%DerRho(:,:,:,2,it), -P,-S, TS,2)
            Density%D2Rho(:,:,:,3,2,it) = &
            & DeriveZ(Density%DerRho(:,:,:,2,it), -P,-S, TS,1)
            
            Density%D2Rho(:,:,:,1,3,it) = &
            & DeriveX(Density%DerRho(:,:,:,3,it), -P, S, TS,1)
            Density%D2Rho(:,:,:,2,3,it) = &
            & DeriveY(Density%DerRho(:,:,:,3,it), -P, S, TS,1)
            Density%D2Rho(:,:,:,3,3,it) = &
            & DeriveZ(Density%DerRho(:,:,:,3,it), -P, S, TS,1)
        enddo
    endif
  end subroutine RecalcRhoDerivatives

  subroutine DIIS(Iteration)
    !---------------------------------------------------------------------------
    ! Simple implementation of the DIIS density mixing scheme.
    !---------------------------------------------------------------------------
    type(DensityVector), save, allocatable :: Difference(:)
    real(KIND=dp), allocatable             :: Matrix(:,:), RHS(:)
    real(KIND=dp)                          :: Work(100)
    integer                                :: i,j, N, Succes
    integer, intent(in)                    :: Iteration
    integer, allocatable                   :: PivotInfo(:)
    type(DensityVector)                    :: TMP

    if(.not.allocated(Difference)) then
      allocate(Difference(PulayOrder))
      do i=1,PulayOrder
        Difference(i) = NewDensityVector()
      enddo
    endif

    !Saving the new difference
    do i=1,PulayOrder-1
      Difference(PulayOrder - i + 1 ) = Difference(PulayOrder - i)
    enddo
    Difference(1) = Density + (-1.0_dp)* DensityHistory(1)

    if(Iteration.le.1) then
      !Return linear damping when not enough info is present
      Density = (1-DampingParam)*Density+DampingParam*DensityHistory(1)
      return
    endif

    !N is the number of mixed iterations
    N = min(Iteration, PulayOrder)
    allocate(Matrix(N+1,N+1), RHS(N+1),PivotInfo(N+1))
    !Computing the matrix and RHS
    Matrix = 0.0_dp

    !DIIS matrix elements are constructed from the density differences.
    do i=1,N
      do j=1,i
        Matrix(i,j) = Difference(i) * Difference(j)
      enddo
    enddo
    RHS             =  0.0_dp
    Matrix(N+1,1:N) = -1.0_dp
    RHS(N+1)        = -1.0_dp

    call DSYSV('L',N+1,1,Matrix,N+1,pivotinfo,RHS,N+1,work,size(Work),Succes)
    if(Succes.ne.0.0_dp) then
      call stp('Error while solving the DIIS linear system.', 'Succes', Succes)
    endif

    !Construct the new density, but retain the density as the old history.
    !In this way we retain a better set of extrapolation vectors.
    TMP = Density
    Density = RHS(1) * Density
    do i=2,N
      Density = Density  + RHS(i) * DensityHistory(i)
    enddo
    DensityHistory(1) = TMP
    where(Density%Rho.lt.0.0_dp)   Density%Rho = 0.0_dp

    return
  end subroutine DIIS


end module DensityMixing
