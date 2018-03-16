module Coulomb

  use CompilationInfo
  use GenInfo
  use CoulombDerivatives
  use Force, only : e2

  implicit none

  save
  !-----------------------------------------------------------------------------
  !Determines the Coulomb Solver used.
  !       0=> No Coulomb
  !       1=> Conjugate Gradients
  !       2=> Red-Black Gauss-Seidel
  !       3=> SOR
  !       4=> SOR with Cheby acceleration
  !       5=> Combination of CG and SOR
  !       6=> Naïve multigrid
  !-----------------------------------------------------------------------------
  integer, public :: CoulombSolver=1

  ! Maximum number of Coulomb iterations.
  integer, public :: CoulombMax=500

  !Number of extra points to take for Coulomb. Default=0
  integer, public :: CEX=0,CEY=0,CEZ=0

  !Maximum l of moments to consider for Coulomb Boundary conditions
  integer, public :: MaxLCoul = -100

  ! This number measures the total time spent solving the Coulomb Problem
  real(KIND=dp), public             :: CoulombTime=0._dp

  !-----------------------------------------------------------------------------
  ! Integers governing the size of the arrays, depending on symmetries.
  !-----------------------------------------------------------------------------
  ! mx/y/z
  !            represent the entire box, enlarged with both extra Coulomb
  !            points and Boundary points.
  ! potmx/y/z
  !            represent the enlarged box, meaning the regular box with
  !            added Coulomb Points, but without Boundary points
  !-----------------------------------------------------------------------------
  integer,public              :: mx,my,mz, potmx,potmy,potmz
  !-----------------------------------------------------------------------------
  ! Number of boundary conditions needed.This depends on the discretisation
  ! of the laplacian for the Coulomb Problem
  !-----------------------------------------------------------------------------
  integer                     :: BCNumber
  !-----------------------------------------------------------------------------
  ! These integers contain the Start (sx/y/z) and End (ex/y/z) indices of
  ! the "middle of the box". With the middle of the box
  ! I mean the Coulomb-extended box, with extra Coulomb points but without
  ! Boundary points. They are used for keeping track of the position of
  ! the nucleus: where the CoulombPotential belongs in PoissonExt.
  !-----------------------------------------------------------------------------
  integer                     :: sx,sy,sz, ex,ey,ez
  !-----------------------------------------------------------------------------
  ! These integers serve a similar purpose: they are the Start(potsx/y/z)
  ! and End (potex/y/z) indices of the normal box
  ! (without Coulomb points) within the Coulomb extended box(with Coulomb
  ! points, without Boundary points).
  ! They are heavily dependent on symmetries.
  !-----------------------------------------------------------------------------
  integer                     :: potsx,potsy,potsz,potex,potey,potez
  !-----------------------------------------------------------------------------
  ! The array containing the Coulomb Potential (in the original box, the
  ! box with extra points and the box with boundary conditions).
  ! Needs to be public, for the HFB module.
  !-----------------------------------------------------------------------------
  real(KIND=dp), allocatable, public :: CoulombPotential(:,:,:)
  real(KIND=dp), allocatable, public :: CPotComplete(:,:,:),Source(:,:,:)
  real(KIND=dp), allocatable, public :: CoulExchange(:,:,:)
  real(KIND=dp), allocatable, public :: CPotentialExt(:,:,:)
  !-----------------------------------------------------------------------------
  !The array containing the boundary conditions for the Coulomb Problem.
  ! Note that memory space is present for the entire box: those points
  ! should remain zero. Not the most efficient.
  !-----------------------------------------------------------------------------
  real(KIND=dp), allocatable,public  :: BoundaryConditions(:,:,:)

  !-----------------------------------------------------------------------------
  ! Arrays containing,
  !-----------------------------------------------------------------------------
  ! 1) the values of the spherical harmonics on the extended mesh and
  ! 2) the value of the radial coordinate r on the extended mesh.
  !-----------------------------------------------------------------------------
  ! Note:
  !  a) as in the module MESH, SpherHarmCoulomb contains not the imaginary
  !    and real parts of Y_{lm} but rather the imaginary and real parts of
  !               r^L \sqrt{4\pi/(2*l+1)} Y^{m}_{l}
  !  b) I know that this is not memory-efficient as a lot of the points
  !   in SpherHarmCoulomb already have been computed in the MOMENTS module.
  !   But since these need to computed in this module one-time only, I'm
  !   currently preferring an implementation that is a bit simpler
  !   over a more efficient implementation .
  !-----------------------------------------------------------------------------
  real(KIND=dp), allocatable, public :: SpherHarmCoulomb(:,:,:,:,:,:)
  real(KIND=dp), allocatable, public :: r(:,:,:)

  !Precision required of the Coulomb Solvers
  real(KIND=dp), public              :: Prec
  !-----------------------------------------------------------------------------
  ! OptimalValues for Omega for SOR, empirically determined for a charged
  ! sphere with 82 protons of radius 7fm/
  ! However, this is not guaranteed to always be the best value.
  !-----------------------------------------------------------------------------
  ! 1-point Lagrangian : 1.91
  ! 2-point Lagrangian : 1.5
  ! 3-point Lagrangian : 1.453
  ! 4-point Lagrangian : 1.39
  !-----------------------------------------------------------------------------
  real(KIND=dp),public,dimension(4) :: OmegaOptimal=                     &
  &                              (/ 1.91_dp, 1.5_dp, 1.453_dp , 1.39_dp /)

contains

  subroutine ReadCoulombInfo()
    !---------------------------------------------------------------------------
    ! Subroutine that reads user input regarding the treatment of the
    !  Coulomb problem.
    !---------------------------------------------------------------------------

    use Moments, only : MaxMoment

    NameList /Coulomb/     CoulombMax,CoulombSolver, CEX,CEY,CEZ, MaxLCoul

    read(unit=*, NML=Coulomb)

    if(CoulombMax.le.0) then
      call stp('Negative number of Coulomb Iterations', &
      &        'CoulombMax', CoulombMax )
    endif
    !Checking if MOCCa recognises the Coulomb Solver
    if(CoulombSolver.ge.7 .or. CoulombSolver.lt.0) then
      call stp('Values for CoulombSolver can only range in 1-7.', &
      &        'CoulombSolver', CoulombSolver)
    endif

    if(MaxLCoul.lt.0) then
      ! If not read, set equal to maxmoment used in the moments subroutine.
      MaxLCoul = MaxMoment
    endif
    if(MaxLCoul.gt.MaxMoment) then
      call stp('MaxLCoul cannot exceed MaxMoment.')
    endif

    !Allocate the memory and initialise various constants for Coulomb
    call SetupCoulomb
    return
  end subroutine ReadCoulombInfo

  subroutine PrintCoulombInfo
    use Force
    
    !---------------------------------------------------------------------------
    ! A subroutine for printing some info on the Coulomb Calculations.
    !
    !---------------------------------------------------------------------------
    1       format(21("-"), 'Coulomb Parameters' , 21('-'))
    2       format("Solver               = ", a25)
    3       format("Discretisation       = ", i2)
    4       format("Omega Parameter      = ", f8.3)
    5       format("Coulomb Exchange     = Slater Approximation")
    6       format("Error Tolerance      = ", e10.3)
    7       format("MaxLCoul             =", i3)
    8       format("Extra points (X/Y/Z) = (",i2,',',i2,',',i2,')' )

    9       format('Exchange included    =', a11)

    character(len=25) :: GS='Gauss-Seidel' , Sor='S. Overrelaxation',          &
    &   CG='Conjugate Gradients', Cheb='GS + Cheby Acceleration',              &
    &   MG='Multigrid', Com = 'Combination'

    print 1

    select case(CoulombSolver)
            case(0)
              print 2, 'No Coulomb Contribution'
            case(1)
                    print 2, CG
            case(2)
                    print 2, GS
            case(3)
                    print 2, Sor
                    print 4, OmegaOptimal(CoulombLapOrder)
            case(4)
                    print 2, Cheb
            case(5)
                    print 2, Com
            case(6)
                    print 2, MG
            case DEFAULT
                    print 2, "Unrecognised Coulomb Solver"
    end select

    print 3, CoulombLapOrder
    print 5
    print 6, Prec
    print 7, MaxLCoul
    print 8, CEX, CEY, CEZ
  
    if(CExchange .eq. 0) then
      print 9, 'None'
    elseif(Cexchange .eq. 1) then
      print 9, 'Perturbative'
    else
      print 9, 'Selfconsistent'
    endif
  
    return
  end subroutine PrintCoulombInfo

  function CompCoulombEnergy(Den) result(CEnergy)
    !---------------------------------------------------------------------------
    !Function that calculates the Coulomb Energy
    !---------------------------------------------------------------------------

    use Densities
    real(KIND=dp) :: CEnergy
    type(DensityVector), intent(in) :: Den

    select case(CoulombSolver)
    case(0)
      CEnergy =  0.0_dp
    case DEFAULT
      CEnergy =  sum(Den%Rho(:,:,:,2) * CoulombPotential)*dv*0.5_dp
    end select
    return
  end function CompCoulombEnergy

  function CompCoulombExchange(Den) result(CExEnergy)
    !---------------------------------------------------------------------------
    ! Function that calculate the exchange energy due to the
    ! electromagnetic interaction.
    ! For the moment, the Slater Approximation is used.
    !---------------------------------------------------------------------------
    use Densities

    type(DensityVector), intent(in) :: Den
    real(KIND=dp)                   :: CExEnergy, Factor

    select case(CoulombSolver)

    case(0)
      CExEnergy = 0.0_dp
    case DEFAULt
      Factor = -(3/pi)**(1/3._dp)*e2*0.75_dp*dv
      CExEnergy = Factor*sum(Den%Rho(:,:,:,2)*Den%Rho(:,:,:,2)**(1/3._dp))
    end select

  end function CompCoulombExchange

  subroutine SolveCoulomb()
    !---------------------------------------------------------------------------
    ! The driving routine of this module.
    ! It does the following
    !       1) If necessary, calls SetupCoulomb to allocate memory
    !          and initialize all variables
    !       2) Then CoulombBound is called to make an initial guess,
    !          with correct boundary conditions.
    !       3) Then passes this first guess to a solver routine,
    !          dependent on the value of Solver.
    !       4) Compute the Coulomb Exchange, for the moment in
    !          Slaters approximation
    !       5) Print some information
    !---------------------------------------------------------------------------
    use Densities, only : Density

    real(KIND=dp)       :: Time0,Time1, Factor
    !---------------------------------------------------------------------------
    ! If the variables are not allocated and initialised yet,
    ! call SetUpCoulomb.
    !---------------------------------------------------------------------------
    if(.not.allocated(CoulombPotential)) then
            call SetUpCoulomb
    endif
    !---------------------------------------------------------------------------
    ! This "first" guess consists of taking boundary conditions
    ! dependent on the multipole moments of the proton distribution.
    ! Then the interior of the box is filled with the Coulomb
    ! potential of the previous solution.
    !---------------------------------------------------------------------------
    call CoulombBound

    !For Timing the solver
    call cpu_time(Time0)

    !Selecting the appropriate coulomb Solver
    select case (CoulombSolver)

            case(0)
                !Don't include Coulomb
                CPotComplete=0.0_dp
                CoulExchange=0.0_dp
                CoulombPotential=0.0_dp
                return

            case(1)
                    !Solve the Coulomb problem using the Simplest
                    ! conjugate Gradients method.
                    ! This is the same as in CR8 and EV8.
                    call ConjugGrad ( CPotComplete,mx,my,mz,       &
                    &    sx,sy,sz,ex,ey,ez, Source, dx,CoulombMax, &
                    &    .false.,Prec)
            case(2)
                    !Use Red-Black Gauss-Seidel
                    call GaussSeidelSolve(.false.,.false.)
            case(3)
                    !Use red-black Succesive Over-relaxation without
                    !Chebyshev acceleration
                    call GaussSeidelSolve(.true.,.false.)
            case(4)
                    !Use red-black Succesive Over-relaxation with
                    ! Chebyshev acceleration
                    call GaussSeidelSolve(.true.,.true.)
            case(5)
                call Combination
            case(6)
                    !Use a multigrid method
                    call MGCoulomb
            case(7)
                    !Solve the Coulomb problem using the method
                    ! described by Reid.
                    call stp('This method is not implemented yet')
            case DEFAULT
                    call stp('CoulombSolver not recognized',       &
                    &        'Solver', CoulombSolver)
    end select

    !Adding the time spent to the total time
    call cpu_time(time1)
    CoulombTime = CoulombTime + Time1 - Time0
    !---------------------------------------------------------------------------
    !Putting the values in the correct spots.
    !---------------------------------------------------------------------------
    CPotentialExt    = CPotComplete(sx:ex,sy:ey,sz:ez)
    CoulombPotential = CPotentialExt(potsx:potex,potsy:potey,      &
                     & potsz:potez)

    ! Computing the Exchange Potential
    Factor = -(3/pi)**(1/3._dp)*e2
    CoulExchange = Factor*Density%Rho(:,:,:,2)**(1/3._dp)
    return
  end subroutine SolveCoulomb

  subroutine SetupCoulomb
  !-----------------------------------------------------------------------------
  ! This routine analyses the symmetries present and uses them to
  ! allocate the correct amount of space for the arrays in this
  ! module. In addition, the SpherHarmCoulomb and r arrays are
  ! calculated, as are the bookkeeping indices (pot)sx/y/z and
  ! (pot)ex/y/z.
  ! A global precision level is set.
  !-----------------------------------------------------------------------------

  use Moments, only : SphericalHarmonics, MaxMoment

  real(KIND=dp),allocatable :: Mesh3DExt(:,:,:,:)
  integer                    :: i,j,k

  !Computing the number of Boundary Conditions necessary
  if(.not.allocated(FDCoulomb)) then
    call stp('Coulomb FD coefficients not initialised!')
  endif
  BCNumber = (size(FDCoulomb)-1)/2
  !-----------------------------------------------------------------------------
  !First, we compute the sizes of the arrays the we will need.
  ! mx/y/z consists of the normal box (nx/y/z)
  !            + Extra Coulomb points(CEX/Y/Z)
  !            + Boundary Points (BCNumber)
  !-----------------------------------------------------------------------------
  mx    = nx + CEX + BCnumber
  my    = ny + CEY + BCNumber
  mz    = nz + CEZ + BCNumber
  !-----------------------------------------------------------------------------
  ! potmx/y/z consists of the normal box (nx/y/z)
  !              + Extra Coulomb points (CEX/Y/Z)
  !-----------------------------------------------------------------------------
  potmx = nx + CEX
  potmy = ny + CEY
  potmz = nz + CEZ
  !-----------------------------------------------------------------------------
  ! sx/y/z are the start points of the Coulomb box
  ! within the mx/y/z box
  !-----------------------------------------------------------------------------
  sx = 1
  sy = 1
  sz = 1
  !-----------------------------------------------------------------------------
  ! ex/y/z are the end points of the Coulomb box within
  ! the mx/y/z box
  !-----------------------------------------------------------------------------
  ex = nx + CEX
  ey = ny + CEY
  ez = nz + CEZ
  !-----------------------------------------------------------------------------
  !potsx/y/z are the start points of the normal box within
  !the Coulomb Box
  !-----------------------------------------------------------------------------
  potsx = 1
  potsy = 1
  potsz = 1
  !-----------------------------------------------------------------------------
  !potzx/y/z are the end points of the normal box within
  !the Coulomb Box
  !-----------------------------------------------------------------------------
  potex = nx
  potey = ny
  potez = nz

  if(.not.SC) then
          !We are dealing with the full x-axis here, so additional
          ! coulomb points plus extra boundary conditions.
          mx    = mx    + CEX + BCNumber
          potmx = potmx + CEX

          !The coulomb box within the mx/y/z box gets shifted in
          ! the x-direction:
          sx = sx       + BCNumber
          ex = ex + CEX + BCNumber

          ! The normal box within the Coulomb Box gets shifted in
          ! the x-direction
          potsx = potsx + CEX
          potex = potex + CEX
  endif
  if(.not.TSC) then
          !We are dealing with the full y-axis here, so additional
          ! coulomb points plus 4 extra boundary conditions.
          my    = my    + CEY + BCNumber
          potmy = potmy + CEY

          !The coulomb box within the mx/y/z box gets shifted in
          !the y-direction:
          sy = sy       + BCNumber
          ey = ey + CEY + BCNumber

          ! The normal box within the Coulomb Box gets shifted in
          !the y-direction
          potsy = potsy + CEY
          potey = potey + CEY
  endif
  if(.not.PC) then
          !We are dealing with the full z-axis here, so additional
          !coulomb points plus extra boundary conditions.
          mz    = mz    + CEZ + BCNumber
          potmz = potmz + CEZ

          !The coulomb box within the mx/y/z box gets shifted in
          !the z-direction:
          sz = sz       + BCNumber
          ez = ez + CEZ + BCNumber

          ! The normal box within the Coulomb Box gets shifted in
          ! the z-direction
          potsz = potsz + CEZ
          potez = potez + CEZ
  endif

  !-----------------------------------------------------------------------------
  !Computing the necessary precision. Note that this cannot be a
  !parameter, since it depends on the symmetries.
  !-----------------------------------------------------------------------------
  Prec = 1.d-9/(dx**3 * potmx * potmy * potmz)

  if(allocated(CoulombPotential)) then
    ! If the variables are already allocated and this
    ! routine is still called, it means that the dimensions
    ! of these things should be re-evaluated.
    deallocate(CoulombPotential)
    deallocate(BoundaryConditions)
    deallocate(r)
    deallocate(SpherHarmCoulomb)
    deallocate(CPotentialExt)
    deallocate(Source)
    deallocate(CPotComplete)
    deallocate(CoulExchange)
  endif

  !Giving each matrix the correct size.
  allocate(CoulombPotential(nx,ny,nz))
  allocate(BoundaryConditions(mx,my,mz))
  allocate(Mesh3DExt(3,mx,my,mz))
  allocate(r(mx,my,mz))
  allocate(SpherHarmCoulomb(mx,my,mz,0:MaxMoment,0:MaxMoment,2))
  allocate(CPotentialExt(potmx,potmy,potmz))
  allocate(Source(mx,my,mz))
  allocate(CPotComplete(mx,my,mz))
  allocate(CoulExchange(nx,ny,nz))

  !Initialising
  CoulombPotential   = 0.0_dp
  BoundaryConditions = 0.0_dp
  r                  = 0.0_dp
  Mesh3DExt          = 0.0_dp
  SpherHarmCoulomb   = 0.0_dp
  CPotentialExt      = 0.0_dp
  CPotComplete       = 0.0_dp
  CoulExchange       = 0.0_dp

  !Filling in all the coordinates for the extended mesh.
  do i=1,mx
          Mesh3DExt(1,i,:,:) = (1/2.0_dp +(i-1))*dx
  enddo

  do j=1,my
          Mesh3DExt(2,:,j,:) = (1/2.0_dp +(j-1))*dx
  enddo

  do k=1,mz
          Mesh3DExt(3,:,:,k) = (1/2.0_dp +(k-1))*dx
  enddo

  !Correcting the coordinates when symmetries are broken.
  if(.not. SC) then
          !If signature is broken, we are dealing with the full x-axis
          Mesh3DExt(1,:,:,:) = Mesh3DExt(1,:,:,:) - mx/2.0_dp*dx
  endif
  if(.not. TSC) then
          !If signature is broken, we are dealing with the full y-axis
          Mesh3DExt(2,:,:,:) = Mesh3DExt(2,:,:,:) - my/2.0_dp*dx
  endif
  if(.not. PC) then
          !If signature is broken, we are dealing with the full z-axis
          Mesh3DExt(3,:,:,:) = Mesh3DExt(3,:,:,:) - mz/2.0_dp*dx
  endif

  ! Computing the spherical harmonics at the extra points.
  call SphericalHarmonics(Mesh3DExt, mx,my,mz, SpherHarmCoulomb)

  !Computing the distances r, for use in Coulombbound.
  r = sqrt(Mesh3DExt(1,:,:,:)**2+ Mesh3DExt(2,:,:,:)**2+ Mesh3DExt(3,:,:,:)**2)

  !Make an initial guess for the complete potential
  CPotComplete     = InitialGuess()

  return
  end subroutine SetupCoulomb

  function InitialGuess()
    !---------------------------------------------------------------------------
    ! This function constructs a (hopefully) good initial guess for the Coulomb
    ! potential. It computes the rms radius of the proton density and then
    ! constructs the potential for a sphere with that radius with a uniform
    ! charge of the number of protons.
    ! It remains to be tested if this is a good initial guess for deformed
    ! systems, but since nuclear deformations are general not that big, I expect
    ! this will not do badly.
    !---------------------------------------------------------------------------
    use Mesh

    real(KIND=dp)              :: Radius, Charge, Factor
    real(KIND=dp)              :: InitialGuess(mx,my,mz)
    integer                    :: i,j,k

    Radius = 0.0_dp

    !The radius formula slightly outperforms the actual calculation,
    !for deformed nuclei.
    Radius = 1.25_dp * (Neutrons + Protons)**(1._dp/3._dp)

    Charge = Protons*e2
    Factor = Charge/(2*Radius)
    do k=sz,ez
      do j=sy,ey
        do i=sx,ex
          if(r(i,j,k).le.Radius) then
            !Inside the sphere V(r) = Q/R * (3  - r**2/R**2)
            InitialGuess(i,j,k) = Factor *(3 - r(i,j,k)**2/Radius**2)
          else
            !Outside the sphere V(r) = Q/R
            InitialGuess(i,j,k) = Charge/r(i,j,k)
          endif
        enddo
      enddo
    enddo

    return
  end function InitialGuess

  subroutine CoulombBound
    !---------------------------------------------------------------------------
    ! This routine computes the source term of the Poisson equation and makes an
    ! initial guess at the Potential. This guess is based on the previous
    ! solution (if available) and boundary conditions computed using the
    ! the multipole moments of the density.
    !---------------------------------------------------------------------------
    ! As a reference:
    !
    ! Poissons equation:
    !  \Delta \phi + e*rho/(\epsilon) = 0
    ! where \phi is the Coulomb Potential, Rho is the proton density, e is the
    ! proton charge and \epsilon is the permittivity of the vacuum.
    !
    ! In general the potential at a distance r, attributed to a multipole moment
    ! Q_lm of order l and magnetic quantum number m is:
    !   V_lm(r) = C_l 1/r^{l+1}  <Q_lm> Y_{lm}
    !
    ! where the expected values are taken with respect to the charge density,
    ! e.g. the proton density and the constant C is:
    !   C_l = e/(4\pi \epsilon_0) \sqrt{4\pi/(2*l+1)}
    !---------------------------------------------------------------------------

    use Moments
    use Densities, only : Density

    class(Moment),pointer :: Current
    real(KIND=dp)         :: factor, Value
    integer               :: l,m,i,j,k,Im
    logical               :: Cont

    !Initialising
    BoundaryConditions = 0.0_dp

    !Computing the Source term
    Source = 0.0_dp
    Source(sx+potsx-1:ex,sy+potsy-1:ey,sz+potsz-1:ez) = - 4*pi*e2*Density%Rho(:,:,:,2)

    !We will now compute the boundary values of the potential by looping over
    !all non-zero multipole moments. Starting down the linked list of moments
    nullify(Current)
    Current => Root

    Cont = .true.
    do while(Cont)

      !Attributes of the multipole moment
      l = Current%l

      m = Current%m
      if(Current%Impart) then
        Im = 2
      else
        Im = 1
      endif

      Value = Current%Value(2)*(4*pi/(2*l+1))

      factor = e2

      do k=1,mz
        do j=1,my
          do i=1,mx
            if( i.gt.ex .or. j.gt.ey .or. k.gt.ez .or.                &
            &   i.lt.sx .or. j.lt.sy .or. k.lt.sz) then
              BoundaryConditions(i,j,k) = BoundaryConditions(i,j,k) + &
              & Value*SpherHarmCoulomb(i,j,k,l,m,Im)*factor/(r(i,j,k)**(2*l+1))
            endif
          enddo
        enddo
      enddo

      !Transferring to the next moment in the list, until the r**2 is reached or
      ! the highest admissible L.
      if(Current%Next%l .ge. 0 .and. Current%Next%l .le. MaxLCoul) then
        Current => Current%Next
      else
      !Signalling that there is no more moment to use.
        Cont=.false.
      endif
    end do
    do k=1,mz
      do j=1,my
        do i=1,mx
          if(BoundaryConditions(i,j,k).ne.0.0_dp) then
            CPotComplete(i,j,k)=BoundaryConditions(i,j,k)
          endif
        enddo
      enddo
    enddo
    return
  end subroutine CoulombBound

  subroutine GaussSeidelSolve(OverRelax,Cheby)
    !---------------------------------------------------------------------------
    ! A driver routine for solving Poissons equation using a Gauss-Seidel
    ! Relaxation Scheme or a succesive overrelaxation scheme.
    !---------------------------------------------------------------------------
    logical, intent(in)       :: OverRelax,Cheby

    if(Cheby) then
      call stp('Chebychev Acceleration does not work yet!')
      !call ChebyAccel(CPotComplete,mx,my,mz,sx,sy,sz,ex,ey,ez,Source,dx,CoulombMax,Prec)
    else
      !Calling the Gauss-Seidel (overrelaxed) scheme without printouts.
      ! Set .false. to .true. for extra printouts.
      call GaussSeidel(CPotComplete,mx,my,mz,sx,sy,sz,ex,ey,ez,Source,dx,      &
      &                CoulombMax,OverRelax,.false.,Prec)
    endif
    return
  end subroutine GaussSeidelSolve

  subroutine GaussSeidel( Solution,gx,gy,gz,x1,y1,z1,x2,y2,z2, SourceTerm, dr, &
  &                        MaxIteration, OverRelax,iprint,Precis)
    !---------------------------------------------------------------------------
    ! An implementation of a red-black Gauss-Seidel algorithm.
    ! Note: dr is probably not correctly used yet.
    !---------------------------------------------------------------------------

    1 format ("Coulomb Calculation converged after", i4, " Iterations.")
    2 format ("Coulomb Calculation did not converge after", i4,                &
              " Iterations. Residual:", e15.8)
    integer,intent(in)                :: gx,gy,gz,MaxIteration,x1,x2,y1,y2,z1,z2
    real(KIND=dp), intent(in)         :: dr,SourceTerm(gx,gy,gz)
    real(KIND=dp), intent(inout)      :: Solution(gx,gy,gz)
    real(KIND=dp), intent(in),optional:: Precis
    logical, intent(in)               :: OverRelax, iprint
    real(KIND=dp)                     :: Residual(gx,gy,gz), Omega
    integer                           :: Iteration, parity=0,Signature=0,      &
    &                                    TimeSimplex=0

    !Checking symmetries
    if(SC) then
            Signature = OneI
    endif
    if(TSC) then
            TimeSimplex = OneI
    endif
    if(PC) then
            Parity = OneI
    endif

    Residual = 0.0_dp

    if(OverRelax) then
           Omega = OmegaOptimal((size(FDCoulomb)-1)/2) !Selecting the best Omega
    else
           Omega=1._dp
    endif

    do Iteration=1,MaxIteration
      !print*, sum(Residual**2), Iteration
      Residual = 0.0_dp

      Residual = Gauss_Seidel_Update &
      & (Solution,gx,gy,gz,SourceTerm,x1,x2,y1,y2,z1,z2,Omega,Signature,Parity,&
      & TimeSimplex,dr)

      !Checking for convergence if necessary
      if(present(Precis)) then
        !Juggling with indices, because the boundary conditions shouldn't
        ! be part of the residual
        if(sum(Residual**2).le.Precis) then
          if(iprint) print 1, Iteration
          exit
        elseif(Iteration.eq.MaxIteration) then
          print 2, Iteration, sum(Residual**2)
        endif
      endif
    enddo
  end subroutine GaussSeidel

  subroutine ConjugGrad &
  & ( Solution,gx,gy,gz,x1,y1,z1,x2,y2,z2, SourceTerm, dr,MaxIteration,iprint, &
  &   Precis)
    !---------------------------------------------------------------------------
    ! This subroutine solves the Coulomb problem. The technique is identical to
    ! the ones employed in EV8 and CR8 and is a straight-forward conjugate
    ! gradients algorithm.
    !
    ! A good explanation of the actual method that is used is given on Wikipedia.
    ! http://en.wikipedia.org/wiki/Conjugate_gradient_method
    ! (Note that in our case the operator A = \Delta and the vector b is the
    !  proton density (multiplied by a constant).)
    ! Probably the same explanation is given in this (ancient) book:
    ! M. Engeli et al., 1959, Birkhäuser Basel
    ! `Refined Iterative Methods for Computation of the Solution and the
    !  Eigenvalues of Self-Adjoint Boundary Value Problems'
    !---------------------------------------------------------------------------
    ! Actual Algorithm
    !
    !  Define the following quantities (indices k refer to iteration k):
    !            "Proton density "        \rho
    !            "Coulomb Potential"      \phi_k
    !            "Actual Residue"         r_k = \rho - \Delta \phi
    !            "Iterative Residue"      p_k
    !            "Conjugate coefficient"  a_k = (r_k^T r_k)/(p_k^T \Delta p_k^T)
    !            "Evolving Coefficient"   b_k = (r_(k+1)^T r_(k+1))/(r_k^T r_k)
    !  Note:
    !    *) Sorry if the names sound obscure, but there are no "accepted" names.
    !             I made them up to have something more tangible.
    !    *) The r_k and p_k are, of course, 3D Meshes, but can also be regarded
    !             as vectors.
    !       Notation like "r_k^T r_k" should thus be translated as
    !            "sum(r_k(i,j,k)* r_k(i,j,k) )".
    !
    !  The actual algorithm goes as follows:
    !    1) From a starting potential CoulombPotential (0 if this is the first
    !       time this algorithm is invoked),
    !       compute:
    !                r_0 = -e/(epsilon)*Rho - \Delta \Phi
    !                p_0 = r_0
    ! |--2)  Then compute
    ! |              a_k       = (r_k^T r_k)/(p_k^T \Delta p_k)
    ! |              phi_(k+1) = phi_k + a_k p_k
    ! |              r_(k+1)   = r_k - a_k \Delta p_k
    ! |  3) Check if sum(|r_k|) is smaller than the desired precision.
    ! |      If so, phi_(k+1) is the desired potential.
    ! |      If not compute:
    ! |              b_k       = (r_(k+1)^T r_(k+1))/(r_k^T r_k)
    ! |              p_(k+1)   = r_(k+1) + b_k p_k
    ! |
    ! |--4) Go to 2).
    !
    !---------------------------------------------------------------------------
    ! In case anyone wants to compare to CR8:
    !               \rho               = rho
    !               \phi               = w
    !               r_k                = z
    !               p_k                = p
    !               \Delta p_k         = t
    !               r_k^T r_k          = zz1
    !               r_(k+1)^T r_(k+1)  = zz2
    !               p_k^T \Delta p_k   = zt
    !               a_k                =  a  = zz1/zt
    !               b_k                =  c  = zz2/zz1
    !---------------------------------------------------------------------------
    1 format ("Coulomb Calculation converged after", i4, " Iterations.")
    2 format ("Coulomb Calculation did not converge after", i4,                &
    &         " Iterations. Residual:", e15.8)

    integer,intent(in)                :: gx,gy,gz,MaxIteration,x1,x2,y1,y2,z1,z2
    real(KIND=dp), intent(in)         :: dr,SourceTerm(gx,gy,gz)
    real(KIND=dp), intent(inout)      :: Solution(gx,gy,gz)
    real(KIND=dp), intent(in),optional:: Precis
    logical, intent(in)               :: iprint

    integer                    :: iteration,Signature=0,TimeSimplex=0,Parity=0
    real(KIND=dp), allocatable :: p_k(:,:,:), Temp(:,:,:),Temp2(:,:,:)
    real(KIND=dp), allocatable :: Residual(:,:,:)
    real(KIND=dp)              :: PoissonNorm, Integral, a_k, c_k
    real(KIND=dp)              :: NewPoissonNorm, Conversion
    !Checking the symmetries for use in the iterations.
    if(SC) then
      Signature = OneI
    endif
    if(TSC) then
      TimeSimplex = OneI
    endif
    if(PC) then
      Parity = OneI
    endif

    allocate(p_k(mx,my,mz))
    allocate(Temp(mx,my,mz))
    allocate(Temp2(mx,my,mz))
    allocate(Residual(mx,my,mz))

    !Conversion factor takes into account that the Laplacian subroutine works on
    ! a mesh with steps of dx
    Conversion = dx**2/dr**2

    !Computing the Residual
    Residual = - Conversion*CoulombLaplacian(Solution,mx,my,mz,Parity,         &
    &                                        Signature,TimeSimplex,OneI,       &
    &                                        .true., x1,y1,z1,x2,y2,z2)
    Residual = Residual + SourceTerm
    !The variable p_k is the conjugate direction. It starts out equal to our
    ! initial Residual.
    p_k=Residual
    !Poissonnorm is r_(k+1)^T r_(k+1)
    PoissonNorm = sum(Residual**2)

    do iteration = 1,MaxIteration
      !Applying lagrangian to p_k
      Temp = Conversion*CoulombLaplacian(p_k,mx,my,mz, Parity,Signature,       &
      &                                  TimeSimplex,OneI,&
      &                                  .true., x1,y1,z1,x2,y2,z2)

      !Integral is p_k^T \Delta p_k^T
      Integral = sum(Residual*Temp)
      !if(Integral.le.1d-3) call stp('Integral')
      !a = zz1/zt
      a_k = PoissonNorm/Integral

      !Increment The Coulomb Potential
      Solution = Solution + a_k*p_k

      !Increment the Poisson Equation
      Residual = Residual - a_k*Temp
      !NewPoissonNorm = zz2
      NewPoissonNorm = sum(Residual**2)
      if(sum(Residual**2).le.Precis) then
        !Check if convergence is reached.
        if(iprint) print 1, Iteration
        exit
      elseif(Iteration.eq.CoulombMax) then
        if(iprint) print 2, Iteration, PoissonNorm
      endif
      !c_k = zz2/zz1
      c_k = NewPoissonNorm/PoissonNorm
      PoissonNorm = NewPoissonNorm
      !Calculate the next conjugate gradient
      p_k = Residual    + c_k*p_k
    enddo

    return
  end subroutine ConjugGrad

  subroutine Combination()
    !---------------------------------------------------------------------------
    ! This subroutine uses Conjugate Gradients and SOR alternating, between
    ! the two. This is not very efficient.
    !---------------------------------------------------------------------------
    integer :: i

    do i=1,CoulombMax/100
      call ConjugGrad &
      & ( CPotComplete,mx,my,mz,sx,sy,sz,ex,ey,ez, Source, dx,50,              &
      &   .false.,Prec)
      call GaussSeidel &
      & ( CPotComplete,mx,my,mz,sx,sy,sz,ex,ey,ez, Source, dx,50,              &
      &   .true.,.false.,Prec)
    enddo
    return
  end subroutine Combination

  subroutine MGCoulomb()
    !---------------------------------------------------------------------------
    ! This subroutine uses the module Multigrid to solve the Coulomb equation.
    !---------------------------------------------------------------------------
    use MultiGrid, only : MGSolve

    call MGSolve (CpotComplete, Source, dx, CoulombMax, Prec)
  end subroutine MGCoulomb


  ! subroutine PrintPotentialRadial(FileName)
  !   !---------------------------------------------------------------------------
  !   ! This subroutine produces a graph of the Coulomb potential along
  !   ! the radial direction.
  !   !---------------------------------------------------------------------------
  !   use GnuFor
  !
  !   integer                         :: i, maxCoord, error
  !   real(KIND=dp), allocatable      :: ToPlotR(:), ToPlotCoulombPotential(:)
  !   character (len = *), intent(in) :: FileName
  !
  !   maxCoord = min(potmx,potmy,potmz)
  !
  !   allocate(ToPlotR(maxcoord))
  !   allocate(ToPlotCoulombPotential(maxcoord))
  !
  !   do i=1,maxCoord
  !          ToPlotR(i) = r(i,i,i)
  !          ToPlotCoulombPotential(i) = CoulombPotential(i,i,i)
  !   enddo
  !
  !   call write_xy_data ( "RadialPotential", maxcoord, ToPlotR,                 &
  !   &                    ToPlotCoulombPotential, error )
  !   if(error.ne.ZeroI) then
  !           call stp('Error in Write_xy_data', "Error", error)
  !   endif
  !
  !   call write_xy_plot ( "RadialPotentialCommand", "RadialPotential", FileName,&
  !   &                    "r", "Potential", error )
  !   if(error.ne.ZeroI) then
  !           call stp('Error in Write_xy_plot', "Error", error)
  !   endif
  !
  !   call Run_GnuPlot("RadialPotentialCommand")
  !
  !   deallocate(ToPlotR); deallocate(ToPlotCoulombPotential)
  !
  ! end subroutine PrintPotentialRadial
end module Coulomb
