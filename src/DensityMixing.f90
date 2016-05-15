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
    ! 1) => DIIS
    ! 2) => CDIIS (not yet)
    ! 3) => CEDIIS (only dim=2 for the moment).
    !---------------------------------------------------------------------------
    integer, intent(in) :: Iteration

    select case(MixingScheme)
      case(0)
        ! Do linear damping
        DensityChange =                                                        &
        &             sum(Density%Rho-DensityHistory(1)%Rho)*dv/(1-DampingParam)
        Density = (1-DampingParam) * Density + DampingParam*DensityHistory(1)
      case(1)
        call DIIS(mod(Iteration,100))
      case(2)
        call stp('CDIIS not properly implemented yet.')
      case(3)
        !Look for the energetically best mixing!
        !call NaiveCEDIIS(Iteration)
        call stp('NaiveCDIIS not properly implemented yet.')
      case DEFAULT
        call stp('Mixing scheme not supported!')
    end select

  end subroutine MixDensities

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

  subroutine CDIIS
    !---------------------------------------------------------------------------
    ! Naïve implementation of CDIIS.
    !
    !---------------------------------------------------------------------------

  end subroutine CDIIS

  subroutine CEDIIS
    !---------------------------------------------------------------------------
    ! Non-naïve version of CEDIIS, using a Downhill Simplex Method for finding
    ! the minimum for the density mixing.
    ! See Numerical Recipes by Press for more info or the wikipedia page
    ! http://en.wikipedia.org/wiki/Nelder-Mead_method
    ! for a good step-by-step explanation.
    !---------------------------------------------------------------------------
    !use Energy, only : EstimateEnergy

    integer                   :: i
    !Location of the N-dimensional minimum
    real(KIND=dp),allocatable :: MinSimplex(:)

  end subroutine CEDIIS

!  subroutine NaiveCEDIIS(Iteration)
!    !---------------------------------------------------------------------------
!    ! Very naïve implementation of 2D CEDIIS.
!    ! The strategy is to search for the minimum of the energy for different
!    ! mixing parameters. This is a quadratic programming problem, and thus
!    ! conceptually very hard to implement for D>2.
!    !---------------------------------------------------------------------------
!    use Energy, only : EstimateEnergy
!
!    integer             :: i, location(1)
!    integer, intent(in) :: Iteration
!    real(KIND=dp)       :: Samples(100)=0.0_dp, Values(100)
!    type(DensityVector) :: Mix

!    if(PulayOrder.gt.1) call stp("Can't do CEDIIS for bigger dimensions yet!")
!
!    !Initialise the Discretisation
!    if(Samples(1).eq.0.0_dp) then
!      do i=1,size(Samples)
!        Samples(i) = -1.0_dp/(2.0_dp * size(Samples)) + i*1.0_dp/(size(Samples))
!      enddo
!    endif
!
!    if(Iteration.le.1) then
!      !Return linear damping when not enough info is present
!      Density = 0.5_dp*Density+0.5_dp*DensityHistory(1)
!      return
!    endif

!    Mix = NewDensityVector()
!    !Estimate the mixing energy
!    do i=1,size(Samples)
!      Mix = Samples(i) * Density + (1.0_dp - Samples(i)) * DensityHistory(1)
!      Values(i) = EstimateEnergy(Mix)
!    enddo
!    location = minloc(Values)
!    !Temporary
!    print *, 'Found minimum estimated energy at ',                             &
!    &       Samples(location), 1 - Samples(Location), Values(Location)
!
!    !Mix in the optimal way!
!    Density=  Samples(Location(1))*Density +                                   &
!    &        (1.0_dp-Samples(Location(1)))*DensityHistory(1)
!  end subroutine NaiveCEDIIS
end module DensityMixing
