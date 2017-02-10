module Derivatives

    use CompilationInfo
    use GenInfo
    use OptimizedDerivatives

    implicit none

    public

    logical :: OptDer=.true.

    !-----------------------------------------------------------------------------
    ! These abstract interface seem to serve no well-defined use, but they are
    ! necessary for compatibility with the PGI compilers.
    ! In addition, they force me to use allocatable types in the derivation
    ! routines, which makes me uncomfortable: in this way the compiler does not
    ! get all available info (namely that the arrays will always be of dimension
    ! (nx,ny,nz)), and I fear this impacts performance.
    abstract interface
      function Central(Grid, Parity, Signature, TimeSimplex,Component) result(Der)
        import :: dp
        integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
        real(KIND=dp), target, intent(in) :: Grid(:,:,:)
        real(KIND=dp),allocatable         :: Der(:,:,:)
      end function Central
    end interface

    abstract interface
      function LineExtension(Grid, N, Y, Z )  result(Ext)
        import :: dp
        integer, intent(in)                :: Y,Z,N
        real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
        real(KIND=dp)                      :: Ext(N)
      end function LineExtension
    end interface

    abstract interface
      function Lapinterface(Grid,Parity, Signature, TimeSimplex,Component) &
      &   result(Lap)
        import                            :: dp
        integer,intent(in)                :: Parity,Signature
        integer,intent(in)                ::TimeSimplex,Component
        real(KIND=dp), target, intent(in) :: Grid(:,:,:)
        real(KIND=dp), allocatable        :: Lap(:,:,:)
      end function Lapinterface
    end interface
    !---------------------------------------------------------------------------
    logical  :: BStack = .true.
    !---------------------------------------------------------------------------
    !Procedure pointers to laplacian and derivative routines.
    procedure(Central), pointer              :: DeriveX, DeriveY, DeriveZ
    procedure(Lapinterface), pointer         :: Laplacian
    procedure(LineExtension), pointer        :: ExtX,ExtY,ExtZ
    !---------------------------------------------------------------------------
    ! Coefficients of the finite order scheme employed.
    ! Note that with "1st-order" we mean a scheme with 2 points for the first
    ! derivative and 3 points for the second derivative in 1D.
    !---------------------------------------------------------------------------
    real(KIND=dp),save,allocatable :: FDCoef(:), FDLap(:),FDCoulomb(:)
    !---------------------------------------------------------------------------
    !Order of the discretisation of various derivatives.
    ! MaxFDOrder      = 1st order derivatives
    ! MaxFDLapOrder   = Laplacian
    ! CoulombLapOrder = Laplacian in the Coulomb Module
    !---------------------------------------------------------------------------
    integer, public :: MaxFDOrder=3, MaxFDLapOrder=4, CoulombLapOrder=2
    !---------------------------------------------------------------------------
    !Coefficients of the Lagrangian Derivatives
    !---------------------------------------------------------------------------
    real(KIND=dp),allocatable :: LagCoefsX(:), LagCoefsY(:), LagCoefsZ(:)
    real(KIND=dp),allocatable :: LagLapX(:), LagLapY(:),LagLapZ(:),LagLapDiag(:)
    !---------------------------------------------------------------------------
    ! Matrix for storing the signs for the derivatives.
    ! 1st index : Cartesian direction:  x, y, z
    ! 2nd index : Parity             : -1, 0, 1
    ! 3rd index : Signature          : -1, 0, 1
    ! 4th index : Time Simplex       : -1, 0, 1 (*)
    ! 5th index : Spinor component   :  1, 2, 3, 4
    ! (*) Note that Time Simplex in principle will never take negative values,
    !     but it is in general useful to also be able to calculate the correct
    !     signs for that case, as deriving some quantities with rather strange
    !     symmetry properties comes in handy for calculating the mean-fields.
    integer :: Signs(3,-1:1,-1:1,-1:1,4)=0
    !---------------------------------------------------------------------------
    ! Finite difference coefficients for the first and second derivative of a
    ! function, emplying 1st, 2nd, 3rd or 4th order schemes.
    ! Coefficients taken on 15/07/2013 from the tables in
    ! http://en.wikipedia.org/wiki/Finite_difference_coefficient
    !---------------------------------------------------------------------------
    real(KIND=dp),parameter,dimension(3) :: FirstOrderFirstDer =               &
    &                       (/ -0.5_dp, 0.0_dp, 0.5_dp/)
    real(KIND=dp),parameter,dimension(5) :: SecondOrderFirstDer=(/ &
    &  1.0_dp/12.0_dp, -2.0_dp/3.0_dp, 0.0_dp, 2.0_dp/3.0_dp, -1.0_dp/12.0_dp /)
    real(KIND=dp),parameter,dimension(7) :: ThirdOrderFirstDer =(/ &
    &  -1.0_dp/60.0_dp, 3.0_dp/20.0_dp,-3.0_dp/4.0_dp,0.0_dp,3.0_dp/4.0_dp,    &
    &  -3.0_dp/20.0_dp,1.0_dp/60.0_dp /)
    real(KIND=dp),parameter,dimension(9) :: FourthOrderFirstDer=(/ &
    &    1.0_dp/280.0_dp,-4.0_dp/105.0_dp,1.0_dp/5.0_dp,-4.0_dp/5.0_dp, &
    &    0.0_dp,4.0_dp/5.0_dp,-1.0_dp/5.0_dp,4.0_dp/105.0_dp,-1.0_dp/280.0_dp /)

    real(KIND=dp),parameter,dimension(3) :: FirstOrderSecondDer                &
    &                             =(/ 1.0_dp,-2.0_dp,1.0_dp /)
    real(KIND=dp),parameter,dimension(5) :: SecondOrderSecondDer=(/ &
    &   -1.0_dp/12.0_dp,4.0_dp/3.0_dp,-5.0_dp/2.0_dp,4.0_dp/3.0_dp,            &
    &   -1.0_dp/12.0_dp /)
    real(KIND=dp),parameter,dimension(7) :: ThirdOrderSecondDer =(/ &
    &    1.0_dp/90.0_dp,-3.0_dp/20.0_dp,3.0_dp/2.0_dp,-49.0_dp/18.0_dp,        &
    &    3.0_dp/2.0_dp,-3.0_dp/20.0_dp,1.0_dp/90.0_dp/)
    real(KIND=dp),parameter,dimension(9) :: FourthOrderSecondDer=(/ &
    & -9.0_dp/8064.0_dp, 128.0_dp/8064.0_dp, -1008.0_dp/8064.0_dp, 1.0_dp,     &
    & -14350.0_dp/8064.0_dp, 1.0_dp, -1008.0_dp/8064.0_dp, 128.0_dp/8064.0_dp, &
    & -9.0_dp/8064.0_dp /)
contains
    subroutine ReadDerivativesInfo()
      !-----------------------------------------------------------------------
      ! Subroutine that reads the user input regarding this module.
      !-----------------------------------------------------------------------
      NameList /Derivatives/ MaxFDOrder, MaxFDLapOrder, CoulombLapOrder, OptDer&
      &                     ,BStack

      !-----------------------------------------------------------------------
      ! Note that the ifort compiler segfaults here for the following:
      !
      ! nullify(DeriveX,DeriveY,DeriveZ,Laplacian)
      !-----------------------------------------------------------------------
      DeriveX => null() ; DeriveY => null()
      DeriveZ => null() ; Laplacian => null()

      !Specifications of the derivative routines
      read (unit=*, nml=Derivatives)

      !Assign the finite difference coefficients
      call AssignFDCoefs(MaxFDOrder, MaxFDLapOrder, CoulombLapOrder)
      !Assign the correct ways to calculate derivatives, as a function of the
      ! symmetries conserved.
      call CompSigns
      call AssignExtensions

    end subroutine ReadDerivativesInfo

    subroutine CompSigns
    !---------------------------------------------------------------------------
    ! Subroutine that initialises the Signs array for all the correct necessary
    ! signs for derivation.
    !---------------------------------------------------------------------------
        integer :: D, P , S, TS, C

        ! Down and dirty for loops for the old, working routine instead of
        ! manually entering the matrix.
        do C=1,4
            do TS=-1,1
                do S=-1,1
                    do P=-1,1
                        do D=1,3
                           Signs(D,P,S,TS,C)=CompSignExtension(D,P,S,TS,C)
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine CompSigns

    subroutine AssignExtensions
    !---------------------------------------------------------------------------
    ! Subroutine that assigns the correct functions to the procedure pointers
    ! ExtX, ExtY and ExtZ, dependent on the symmetries involved.
    !---------------------------------------------------------------------------

        if(SZC .and. (.not. SC)) then
                call stp('nonoptimized derivatives cannot be used (yet).')        
        endif
        if(SC) then
            if(TSC) then
                ExtX => LineExtensionX_SandTS
            else
                ExtX => LineExtensionX_S
            endif
        else
            ExtX => ZeroExtension
        endif

        if(TSC) then
            ExtY => LineExtensionY_TS
        else
            ExtY => ZeroExtension
        endif

        if(PC) then
            if(SC) then
                ExtZ => LineExtensionZ_PandS
            else
                if(TSC) then
                    ExtZ => LineExtensionZ_PandTS
                else
                    ExtZ => LineExtensionZ_P
                endif
            endif
        else
            ExtZ => ZeroExtension
        endif

    end subroutine AssignExtensions

    subroutine AssignFDCoefs(FDScheme, FDLapScheme, FDCoulombScheme)
      !-------------------------------------------------------------------------
      ! Assign the correct finite difference coefficients.
      !
      !-------------------------------------------------------------------------

      integer, intent(in) :: FDScheme, FDLapScheme, FDCoulombScheme
      integer             :: TrueNX, TrueNY, TrueNZ

      if(allocated(FDCoef   )) deallocate(FDCoef)
      if(allocated(FDLap    )) deallocate(FDLap)
      if(allocated(FDCoulomb)) deallocate(FDCoulomb)

      !Allocating Variables
      allocate(FDCoef(-FDScheme:FDScheme))
      allocate(FDLap(-FDLapScheme:FDLapScheme))
      allocate(FDCoulomb(-FDCoulombScheme:FDCoulombScheme))

      !Allocating procedure pointers
      if(FDScheme.eq.-1) then
        !Use Lagrangian Derivatives
        DeriveX   => DerLagX
        DeriveY   => DerlagY
        DeriveZ   => DerlagZ
        Laplacian => Laplacian_Lag
       else
        !Check for optimized derivatives
        if(OptDer) then
          if(SC .and. PC .and. TSC) then
              ! All spatial symmetries conserved
              DeriveX => Opt_X_EV8
              DeriveY => Opt_Y_EV8
              DeriveZ => Opt_Z_EV8
              Laplacian => Lapla_EV8
          elseif(SC .and. TSC) then
              ! Parity broken, signature and timesimplex conserved
              DeriveX => Opt_X_EV8
              DeriveY => Opt_Y_EV8
              DeriveZ => Opt_Z_EV4
              Laplacian => Lapla_EV4
          elseif(PC .and. TSC) then
              ! Signature broken, parity and timesimplex conserved
              DeriveX => Opt_X_NOSIG
              DeriveY => Opt_Y_EV8
              DeriveZ => CentralZ
              Laplacian => Laplacian_Central
          elseif(SZC .and. TSC) then          
              ! parity and signature broken, simplex conserved
              DeriveX => Opt_X_NOSIG
              DeriveY => Opt_X_EV8
              DeriveZ => Opt_Z_EV8
          elseif(TSC) then
              ! Parity broken, signature broken and timesimplex conserved
              DeriveX   => Opt_X_NOSIG
              DeriveY   => Opt_Y_EV8
              DeriveZ   => OPT_Z_EV4
              Laplacian => lapla_ev4_nosig         
          else
              DeriveX   => CentralX
              DeriveY   => CentralY
              DeriveZ   => CentralZ
              Laplacian => Laplacian_Central
          endif
        else
          DeriveX   => CentralX
          DeriveY   => CentralY
          DeriveZ   => CentralZ
          Laplacian => Laplacian_Central
        endif
      endif

      !Selecting the appropriate coefficients for the First Derivative
      select case (FDScheme)
        case(1)
                FDCoef = FirstOrderFirstDer
        case(2)
                FDCoef = SecondOrderFirstDer
        case(3)
                FDCoef = ThirdOrderFirstDer
        case(4)
                FDCoef = FourthOrderFirstDer
        case(-1)
                TrueNx = nx*2**SignatureInt
                TrueNy = ny*2**TimeSimplexInt
                TrueNz = nz*2**ParityInt
                call Inilag(TrueNx,TrueNy,TrueNz)

                deallocate(FDCoef, FDLap)
        case default
                call stp('This Finite Difference scheme is not supported.', &
               &  'FDScheme', FDScheme)
        end select
        !Selecting the appropriate coefficients for the Second Derivative
        select case (FDLapScheme)
          case(1)
                  FDLap  = FirstOrderSecondDer
          case(2)
                  FDLap  = SecondOrderSecondDer
          case(3)
                  FDLap  = ThirdOrderSecondDer
          case(4)
                  FDLap  = FourthOrderSecondDer
                  ! Note that the weird numbers for the second derivative
                  ! are "legacy" from cr8
                  FDLap  = FDLAP*(8064.0_dp)/(5040.0_dp)
          case(-1)
            if(FDScheme.ne.-1) &
          & call stp(          &
          & "Use lagrangian derivatives for the laplacian and the normal ones!")

          case Default
            call stp(          &
            &'This order for the laplacian is not supported.',                 &
            &'FDLapScheme', FDLapScheme)
        end select

        select case (FDCoulombScheme)
          case(1)
            FDCoulomb  = FirstOrderSecondDer
          case(2)
            FDCoulomb  = SecondOrderSecondDer
          case(3)
            FDCoulomb  = ThirdOrderSecondDer
          case(4)
            FDCoulomb  = FourthOrderSecondDer
          case Default
            call stp(&
            &'This Finite Difference scheme for Coulomb is not supported.',    &
            & 'FDCoulombScheme', FDCoulombScheme)
        end select

    end subroutine AssignFDCoefs

    pure function CompSignExtension(                                           &
    &                Direction,Parity,Signature,TimeSimplex,Component) result(S)
    !---------------------------------------------------------------------------
    ! This function computes the sign associated with the LineExtension,
    ! dependent on the symmetries of the component of the wavefunction.
    !---------------------------------------------------------------------------
    integer, intent(in) :: Direction, Parity,Signature,TimeSimplex,Component
    integer             :: S

    S=0

    !X Direction
    if(Direction.eq.1) then
        if(Signature.ne.0) then
            if(TimeSimplex.ne.0) then
                !Using Signature and Time Simplex
                if(Component.eq.1) then
                        !Spin up real
                        S= Signature*TimeSimplex
                elseif(Component.eq.2) then
                        !Spin up imaginary
                        S=-Signature*TimeSimplex
                elseif(Component.eq.3) then
                        !Spin down real
                        S=-Signature*TimeSimplex
                elseif(Component.eq.4) then
                        !Spin down Imaginary
                        S=Signature*TimeSimplex
                endif
            else
                !Using Signature
                if((Component.eq.1).or.(Component.eq.2)) then
                        !Spin up
                        S=Signature
                else
                        !Spin down
                        S=-Signature
                endif
            endif
        else
            S=0
        endif
    !Y Direction
    elseif(Direction.eq.2) then
        if(TimeSimplex.ne.0) then
            !Using Time Simplex
            if((Component.eq.1).or.(Component.eq.3)) then
                    !The Real Components
                    S=TimeSimplex
            else
                    !The Imaginary Components
                    S=-TimeSimplex
            endif
        else
            S=0
        endif
    ! Z Direction
    elseif(Direction.eq.3) then
        if(Parity.ne.0) then
            if(Signature.ne.0) then
                ! Using Signature and Parity
                if((Component.eq.1).or.(Component.eq.2)) then
                    !Spin up components
                    S= Signature*Parity
                elseif((Component.eq.3).or.(Component.eq.4)) then
                    !Spin down components
                    S=-Signature*Parity
                endif
            else
                if(TimeSimplex.ne.0) then
                    !Using Parity and TimeSimplex
                    if((Component.eq.1).or.(Component.eq.3)) then
                        !Real Components
                        S=Parity*TimeSimplex
                    else
                        !Imaginary Components
                        S=-Parity*TimeSimplex
                    endif
                else
                        !Using Parity
                        S=Parity
                endif
            endif
        else
                S=0
        endif
    else
        S=0
    endif

    return
  end function CompSignExtension

  function ZeroExtension(Grid, N, Y, Z) result(LineExtension)
  !---------------------------------------------------------------------------
  ! Placeholder function for returning 'no extension'.
  !---------------------------------------------------------------------------
      integer, intent(in)                :: Y,Z,N
      real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
      real(KIND=dp)                      :: LineExtension(N)

      LineExtension=0.0_dp

  end function ZeroExtension

  function LineExtensionX_SandTS(Grid, N, Y, Z )  result(LineExtensionX)
  !---------------------------------------------------------------------------
  ! Returns the lineextension at points Y and Z of the Grid, when both
  ! Signature and TimeSimplex are conserved.
  !---------------------------------------------------------------------------
      integer, intent(in)                :: Y,Z,N
      real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
      real(KIND=dp)                      :: LineExtensionX(N)

      LineExtensionX=Grid(1:N,Y,Z)

  end function LineExtensionX_SandTS

  function LineExtensionX_S(Grid, N, Y, Z )  result(LineExtensionX)
  !---------------------------------------------------------------------------
  ! Returns the lineextension in X direction at points Y and Z of the Grid,
  ! when Signature is conserved, but TimeSimplex is not.
  !---------------------------------------------------------------------------
      integer, intent(in)                :: Y,Z,N
      real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
      real(KIND=dp)                      :: LineExtensionX(N)

      LineExtensionX=Grid(1:N,ny - Y + 1,Z)

  end function LineExtensionX_S

  function LineExtensionY_TS(Grid, N, X, Z )  result(LineExtensionY)
  !---------------------------------------------------------------------------
  ! Returns the lineextension in the Y direction at points X and Z of the Grid
  ! when TimeSimplex is conserved.
  !---------------------------------------------------------------------------
      integer, intent(in)                :: X,Z,N
      real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
      real(KIND=dp)                      :: LineExtensionY(N)

      LineExtensionY=Grid(X,1:N,Z)

  end function LineExtensionY_TS

  function LineExtensionZ_PandS(Grid,N,X,Y) result(LineExtensionZ)
  !---------------------------------------------------------------------------
  ! Returns the lineextension in the Z direction at points X and Y of the Grid
  ! when Parity and Signature are conserved.
  !---------------------------------------------------------------------------
      integer, intent(in)                :: X,Y,N
      real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
      real(KIND=dp)                      :: LineExtensionZ(N)

      LineExtensionZ=Grid(X,Y,1:N)

  end function LineExtensionZ_PandS

  function LineExtensionZ_PandTS(Grid,N,X,Y) result(LineExtensionZ)
  !---------------------------------------------------------------------------
  ! Returns the lineextension in the Z direction at points X and Y of the Grid
  ! when Parity and TimeSimplex are conserved, but signature is not.
  !---------------------------------------------------------------------------
      integer, intent(in)                :: X,Y,N
      real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
      real(KIND=dp)                      :: LineExtensionZ(N)

      LineExtensionZ=Grid(nx-X+1,Y,1:N)

  end function LineExtensionZ_PandTS

  function LineExtensionZ_P(Grid,N,X,Y) result(LineExtensionZ)
  !---------------------------------------------------------------------------
  ! Returns the lineextension in the Z direction at points X and Y of the Grid
  ! when only Parity is conserved.
  !---------------------------------------------------------------------------
      integer, intent(in)                :: X,Y,N
      real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
      real(KIND=dp)                      :: LineExtensionZ(N)

      LineExtensionZ=Grid(nx-X+1,ny-Y+1,1:N)

  end function LineExtensionZ_P

  function CentralX(Grid, Parity, Signature, TimeSimplex,Component) result(Der)
    !---------------------------------------------------------------------------
    ! This function splits the derivation into 1D problems, using Central_1D.
    ! Take as example the X derivation. For each Y- and Z-coordinate the line
    ! of Grid with varying X-coordinate is taken.
    ! The program then looks for the correct extension of this line (by using
    ! symmetries if possible), and whether or not the line picks up a sign.
    ! The extension is found with PointToLineExtensionX, wich returns a pointer
    ! to the correct variables.
    ! All this is then put into finite difference formulas in Central_1D
    !---------------------------------------------------------------------------
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    real(KIND=dp)                     :: LineExtension((size(FDCoef)-1)/2)
    integer                           :: j,k, SignExtension,N

    allocate(Der(nx,ny,nz)); Der=0.0_dp
    !N is the number of points any of the LineExtension subroutines needs
    N = (size(FDCoef) - 1)/2
    SignExtension = Signs(1,Parity,Signature,TimeSimplex,Component)

    do k=1,nz
      do j=1,ny
        !LineExtension=PointToLineExtensionX(Grid, N, Parity, TimeSimplex,      &
        !&                                   Signature,j,k)
        LineExtension=ExtX(Grid,N,j,k)
        Der(:,j,k)=Central_1DX(Grid(:,j,k),LineExtension,N,SignExtension,1)
      enddo
    enddo
  end function CentralX

  function CentralY(Grid, Parity, Signature, TimeSimplex,Component) result(Der)
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    real(KIND=dp)                     :: LineExtension((size(FDCoef)-1)/2)
    integer                           :: i,k, SignExtension,N

    allocate(Der(nx,ny,nz)) ; Der=0.0_dp
    !N is the number of points any of the LineExtension subroutines needs
    N = (size(FDCoef) - 1)/2
    ! Y Derivatives
    SignExtension=Signs(2,Parity,Signature,TimeSimplex,Component)
    do k=1,nz
      do i=1,nx
        !LineExtension=PointToLineExtensionY(Grid, N,Parity, TimeSimplex,       &
        !&                                   Signature,i,k)
        LineExtension=ExtY(Grid,N,i,k)
        Der(i,:,k)=Central_1DY(Grid(i,:,k),LineExtension,N,SignExtension,1)
      enddo
    enddo
  end function CentralY

  function CentralZ(Grid, Parity, Signature, TimeSimplex,Component) result(Der)
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp), allocatable        :: Der(:,:,:)
    real(KIND=dp)                     :: LineExtension((size(FDCoef)-1)/2)
    integer                           :: i,j, SignExtension,N

    allocate(Der(nx,ny,nz)); Der=0.0_dp

    !N is the number of points any of the LineExtension subroutines needs
    N = (size(FDCoef) - 1)/2
    SignExtension=Signs(3,Parity,Signature,TimeSimplex,Component)
    do j=1,ny
        do i=1,nx
            !LineExtension=PointToLineExtensionZ(Grid,N, Parity, TimeSimplex,   &
            !                                    Signature,i,j)
            LineExtension=ExtZ(Grid,N,i,j)
            Der(i,j,:)=Central_1DZ(Grid(i,j,:),LineExtension,N,SignExtension,1)
        enddo
    enddo
  end function CentralZ

  function Central_1DX(Line,LineExtension,N,SignExtension,D) result(DerLine)
    !---------------------------------------------------------------------------
    ! This function computes the derivative of degree D of the 1D array Line.
    ! To do this, one needs the extension of Line: LineExtension.
    ! N is the size of this extension, while SignExtension is the sign of this
    ! extension. If LineExtension is completely zero, then there is no relevant
    ! symmetry.
    !---------------------------------------------------------------------------
    integer, intent(in)                   :: N,SignExtension,D
    !S and E are the first and last points that need to be derived.
    real(KIND=dp), intent(in)             :: Line(:)
    real(KIND=dp)                         :: DerLine(nx), Factor, p
    real(KIND=dp), intent(in)             :: LineExtension(N)
    integer                               :: i,j
    real(KIND=dp)                         :: Coefs(-N:N)

    !Deciding the correct coefficients
    select case(D)
      case(1)
              Coefs = FDCoef
      case(2)
              Coefs = FDLap
    end select

    DerLine=0.0_dp

    p =real(SignExtension, KIND=dp)
    !Things on the start of the line
    do i=1,N
        !Going forward is no problem, and going backward on the line itself
        do j=-i+1,N
          DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
        enddo
        !Going backward on the Lineextension
        do j=-N,-i
          Derline(i) = Derline(i) + p*coefs(j)*LineExtension(1-(i+j))
        enddo
    enddo

    !Things on the middle of the line
    do i=1+N, nx-N
      do j=-N,N
    DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
      enddo
    enddo

    !Things on the end of the line
    do i=nx-N+1, nx
        !Going backward is no problem, going forward till the end of the line
        do j=-N,nx-i
            DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
        enddo
    enddo

    factor = 1.0_dp/dx

    if(D.eq.2) factor = factor/dx
    DerLine=DerLine*factor
  end function Central_1DX

  function Central_1DY(Line,LineExtension,N,SignExtension,D) result(DerLine)
    !---------------------------------------------------------------------------
    ! This function computes the derivative of degree D of the 1D array Line.
    ! To do this, one needs the extension of Line: LineExtension.
    ! N is the size of this extension, while SignExtension is the sign of this
    ! extension. If LineExtension is completely zero, then there is no relevant
    ! symmetry.
    !---------------------------------------------------------------------------
    integer, intent(in)                   :: N,SignExtension,D
    !S and E are the first and last points that need to be derived.
    real(KIND=dp), intent(in)             :: Line(ny)
    real(KIND=dp)                         :: DerLine(ny), Factor, p
    real(KIND=dp), intent(in)             :: LineExtension(N)
    integer                               :: i,j
    real(KIND=dp)                         :: Coefs(-N:N)

    !Deciding the correct coefficients
    select case(D)
        case(1)
                Coefs = FDCoef
        case(2)
                Coefs = FDLap
    end select

    DerLine=0.0_dp

    p =real(SignExtension, KIND=dp)

    !Things on the start of the line
    do i=1,N
        !Going forward is no problem, and going backward on the line itself
        do j=-i+1,N
               DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
        enddo
        !Going backward on the Lineextension
        do j=-N,-i
              Derline(i) = Derline(i) + p*coefs(j)*LineExtension(1-(i+j))
        enddo
    enddo

    !Things on the middle of the line
    do i=1+N, ny-N
              do j=-N,N
                      DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
              enddo
    enddo

    !Things on the end of the line
    do i=ny-N+1, ny
        !Going backward is no problem, going forward till the end of the line
        do j=-N,ny-i
              DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
        enddo
    enddo

    factor = 1.0_dp/dx
    if(D.eq.2) factor = factor/dx
    DerLine=DerLine*factor

  end function Central_1DY

  function Central_1DZ(Line,LineExtension,N,SignExtension,D) result(DerLine)
    !---------------------------------------------------------------------------
    ! This function computes the derivative of degree D of the 1D array Line.
    ! To do this, one needs the extension of Line: LineExtension.
    ! N is the size of this extension, while SignExtension is the sign of this
    ! extension. If LineExtension is completely zero, then there is no relevant
    ! symmetry.
    !---------------------------------------------------------------------------
    integer, intent(in)                   :: N,SignExtension,D
    !S and E are the first and last points that need to be derived.
    real(KIND=dp), intent(in)             :: Line(nz)
    real(KIND=dp)                         :: DerLine(nz), Factor, p
    real(KIND=dp), intent(in)             :: LineExtension(N)
    integer                               :: i,j
    real(KIND=dp)                         :: Coefs(-N:N)

    !Deciding the correct coefficients
    select case(D)
      case(1)
          Coefs = FDCoef
      case(2)
          Coefs = FDLap
    end select

    DerLine=0.0_dp

    p =real(SignExtension, KIND=dp)

    !Things on the start of the line
    do i=1,N
        !Going forward is no problem, and going backward on the line itself
        do j=-i+1,N
               DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
        enddo

        !Going backward on the Lineextension
        do j=-N,-i
              Derline(i) = Derline(i) + p*coefs(j)*LineExtension(1-(i+j))
        enddo
    enddo

    !Things on the middle of the line
    do i=1+N, nz-N
              do j=-N,N
                      DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
              enddo
    enddo

    !Things on the end of the line
    do i=nz-N+1, nz
        !Going backward is no problem, going forward till the end of the line
        do j=-N,nz-i
              DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
        enddo
    enddo

    factor = OneR/dx
    if(D.eq.2) factor = factor/dx
    DerLine=DerLine*factor

  end function Central_1DZ

  function Laplacian_Central(Grid,Parity, Signature, TimeSimplex,Component)    &
  &        result(Lap)
    !---------------------------------------------------------------------------
    ! This function splits the calculation of the Laplacian into 1D problems,
    ! using Central_1D. The structure is completely analogous to the Central
    ! subroutine, only with an added option for conserving boundary conditions
    ! in a Coulomb calculation.
    !---------------------------------------------------------------------------
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp), allocatable        :: Lap(:,:,:)
    real(KIND=dp)                     :: DerX(nx,ny,nz)
    real(KIND=dp)                     :: DerZ(nx,ny,nz), DerY(nx,ny,nz)
    real(KIND=dp)                     :: LineExtension((size(FDLap)-1)/2)
    integer                           :: i,j,k, SignExtension,N

    allocate(Lap(nx,ny,nz)); Lap = 0.0_dp

    !N is the number of points any of the LineExtension subroutines needs
    N = (size(FDLap) - 1)/2

    Lap=0.0_dp
    DerX=0.0_dp
    DerY=0.0_dp
    DerZ=0.0_dp

    ! X Derivatives
    SignExtension=Signs(1,Parity,Signature,TimeSimplex,Component)
    do k=1,nz
      do j=1,ny
        LineExtension=ExtX(Grid,N,j,k)
        DerX(:,j,k)=Central_1DX(Grid(:,j,k),LineExtension,N,SignExtension,2)
      enddo
    enddo

    ! Y Derivatives
    SignExtension=Signs(2,Parity,Signature,TimeSimplex,Component)
    do k=1,nz
      do i=1,nx
        LineExtension=ExtY(Grid,N,i,k)
        DerY(i,:,k)=Central_1DY(Grid(i,:,k),LineExtension,N,SignExtension,2)
      enddo
    enddo

    !Z Derivatives
    SignExtension=Signs(3,Parity,Signature,TimeSimplex,Component)
    do j=1,ny
      do i=1,nx
        LineExtension=ExtZ(Grid,N,i,j)
        DerZ(i,j,:)=Central_1DZ(Grid(i,j,:), LineExtension,N,SignExtension,2)
      enddo
    enddo

    Lap = DerX+DerY+DerZ
    return
  end function Laplacian_Central

  function SecondDer_Central(Grid,Direction, Parity, Signature,                &
  &                          TimeSimplex,Component)                            &
  &        result(SecondDer)
    !---------------------------------------------------------------------------
    ! This function only calculates the second derivative in one dimension.
    !---------------------------------------------------------------------------
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    integer,intent(in)                :: Direction
    real(KIND=dp), target, intent(in) :: Grid(nx,ny,nz)
    real(KIND=dp)                     :: SecondDer(nx,ny,nz)
    real(KIND=dp)                     :: LineExtension((size(FDLap)-1)/2)
    integer                           :: i,j,k, SignExtension,N

    !N is the number of points any of the LineExtension subroutines needs
    N = (size(FDLap) - 1)/2

    SecondDer=0.0_dp

    select case(Direction)

    case(1)
      ! X Derivatives
      SignExtension=Signs(1,Parity,Signature,TimeSimplex,Component)
      do k=1,nz
        do j=1,ny
          LineExtension=ExtX(Grid,N,j,k)
          SecondDer(:,j,k)=Central_1DX(Grid(:,j,k),LineExtension,N,SignExtension,2)
        enddo
      enddo

    case(2)
      ! Y Derivatives
      SignExtension=Signs(2,Parity,Signature,TimeSimplex,Component)
      do k=1,nz
        do i=1,nx
          LineExtension=ExtY(Grid,N,i,k)
          SecondDer(i,:,k)=Central_1DY(Grid(i,:,k),LineExtension,N,SignExtension,2)
        enddo
      enddo

      !Z Derivatives
      SignExtension=Signs(3,Parity,Signature,TimeSimplex,Component)
      do j=1,ny
        do i=1,nx
          LineExtension=ExtZ(Grid,N,i,j)
          SecondDer(i,j,:)=Central_1DZ(Grid(i,j,:), LineExtension,N,SignExtension,2)
        enddo
      enddo
    end select
    return
  end function SecondDer_Central

  subroutine IniLag(XLineSize,YLineSize,ZLineSize)
    !---------------------------------------------------------------------------
    !This subroutine initialises the coefficients needed for performing Lagrange
    ! derivatives.  The integers X/Y/ZLineSize need to be the actual sizes of
    ! the entire box, not nx,ny and nz.
    ! According to
    ! D.Baye & P-H. Heenen, Generalised Meshes for Quantum Mechanical Problems
    !        J.Phys. A: Meth. Gen 19 (1986)
    !
    !    \Psi'(i)=\pi/N \sum_{j\not= i} (-1)^{i-j}[sin(pi(i-j)/N)]^{-1}
    !    \Psi''(i)= pi^2/3*(1-1/(N^2))\Psi(i) + 2 pi^2/(N^2)
    !              \sum_{j\not= i} cos (pi(i-j)/N) sin(pi(i-j)/N)^{-2}
    !     where N is the number of points on the line.
    !
    ! The arrays LagCoefX/Y/Z contain the values of the coefficients as function
    ! of point distance.
    ! So :
    !    LagCoefX(k) = \pi/(N*dx) (-1)^k[sin(pi(k)/N)]^{-1}
    !    LagLap(k)   = -1^{k+1}*2*pi^2/(N^2*dx) cos(pi(k)/N) sin(pi(k)/N)^{-2}
    !---------------------------------------------------------------------------
    integer, intent(in) :: XLineSize,YLineSize,ZLineSize
    integer             :: N,k
    real(KIND=dp)       :: Factor,dxFactor,Sine,Cosine

    if(allocated(LagCoefsX )) deallocate(LagCoefsX )
    if(allocated(LagCoefsY )) deallocate(LagCoefsY )
    if(allocated(LagCoefsZ )) deallocate(LagCoefsZ )
    if(allocated(LagLapX   )) deallocate(LagLapX   )
    if(allocated(LagLapY   )) deallocate(LagLapY   )
    if(allocated(LagLapZ   )) deallocate(LagLapZ   )
    if(allocated(LagLapDiag)) deallocate(LagLapDiag)

    allocate (LagCoefsX(XLineSize))
    allocate (LagCoefsY(YLineSize))
    allocate (LagCoefsZ(ZLineSize))
    allocate (LagLapX(XLineSize))
    allocate (LagLapY(YLineSize))
    allocate (LagLapZ(ZLineSize))
    allocate (LagLapDiag(3))

    !X Type
    N=XLineSize
    Factor=pi/N
    dxFactor=Factor/dx
    LagLapDiag(1)=-pi**2/(ThreeR*dx**2)*(OneR-OneR/(N)**2)
    do k=1,N-1
        Sine=sin(Factor*k)
        Cosine=cos(Factor*k)
        LagCoefsX(k)=dxFactor*(-1)**k/Sine
        LagLapX(k)=(-1)**(k+1)*TwoR*(dxFactor**2)*Cosine/(Sine**2)
    enddo
    !Y Type
    N=YLineSize
    Factor=pi/N
    dxFactor=Factor/dx
    LagLapDiag(2)=-pi**2/(ThreeR*dx**2)*(OneR-OneR/(N)**2)
    do k=1,N-1
        Sine=sin(Factor*k)
        Cosine=cos(Factor*k)
        LagCoefsY(k)=dxFactor*(-1)**k/Sine
        LagLapY(k)=(-1)**(k+1)*TwoR*(dxFactor**2)*Cosine/(Sine**2)
    enddo
    !Z Type
    N=ZLineSize
    Factor=pi/N
    dxFactor=Factor/dx
    LagLapDiag(3)=-pi**2/(ThreeR*dx**2)*(OneR-OneR/(N)**2)
    do k=1,N-1
        Sine=sin(Factor*k)
        Cosine=cos(Factor*k)
        LagCoefsZ(k)=dxFactor*(-1)**k/Sine
        LagLapZ(k)=(-1)**(k+1)*TwoR*(dxFactor**2)*Cosine/(Sine**2)
    enddo
  end subroutine IniLag

  function DerlagX(Grid, Parity, Signature, TimeSimplex,Component)result(Der)
    !---------------------------------------------------------------------------
    ! This function splits the derivation into 1D problems, using Central_1D.
    ! For example the X derivation. For each Y- and Z-coordinate the line
    ! with varying X-coordinate is taken. The program then looks for the correct
    ! extension of this line (by using symmetries) and if the line picks up a
    ! sign. The extension is found with PointToLineExtensionX, wich returns a
    ! pointer to the correct variables.
    ! All this is then put into eq. 5.3 in
    !      D. Baye & P. Heenen - J.Phys.A.: Math.Gen. 19 (1986)
    ! in Derlag_1d
    !---------------------------------------------------------------------------
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    real(KIND=dp)                     :: LineExtension(nx)
    integer                           :: j,k, SignExtension

    allocate(Der(nx,ny,nz)) ; Der = 0.0_dp
    ! X Derivatives
    SignExtension=Signs(1,Parity,Signature,TimeSimplex,Component)
    do k=1,nz
      do j=1,ny
        !LineExtension=PointToLineExtensionX(                                   &
        !&                           Grid, nx,Parity, TimeSimplex, Signature,j,k)
        LineExtension=ExtX(Grid,nx,j,k)
        Der(:,j,k)=DerLag_1D(                                                  &
        &                 Grid(:,j,k), nx,LineExtension,SignExtension,LagCoefsX)
      enddo
    enddo
    return
  end function DerlagX

  function DerLagY(Grid, Parity, Signature, TimeSimplex,Component)result(Der)
    !---------------------------------------------------------------------------
    ! This function splits the derivation into 1D problems, using Central_1D.
    ! For example the Y derivation. For each X- and Z-coordinate the line
    ! with varying Y-coordinate is taken. The program then looks for the correct
    ! extension of this line (by using symmetries) and if the line picks up a
    ! sign. The extension is found with PointToLineExtensionY, wich returns a
    ! pointer to the correct variables.
    ! All this is then put into eq. 5.3 in
    !      D. Baye & P. Heenen - J.Phys.A.: Math.Gen. 19 (1986)
    ! in Derlag_1d
    !---------------------------------------------------------------------------
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    real(KIND=dp)                     :: LineExtension(ny)
    integer                           :: i,k, SignExtension

    allocate(Der(nx,ny,nz)) ; Der = 0.0_dp
    ! Y Derivatives
    SignExtension=Signs(2,Parity,Signature,TimeSimplex,Component)
    do k=1,nz
      do i=1,nx
        !LineExtension=PointToLineExtensionY(                                   &
        !&                          Grid, ny, Parity, TimeSimplex, Signature,i,k)
        LineExtension=ExtY(Grid,ny,i,k)
        Der(i,:,k)=DerLag_1D(                                                  &
        &                 Grid(i,:,k), ny,LineExtension,SignExtension,LagCoefsY)
      enddo
    enddo
  end function DerlagY

  function DerlagZ(Grid, Parity, Signature, TimeSimplex,Component)result(Der)
    !---------------------------------------------------------------------------
    ! This function splits the derivation into 1D problems, using Central_1D.
    ! For example the Z derivation. For each X- and Y-coordinate the line
    ! with varying Z-coordinate is taken. The program then looks for the correct
    ! extension of this line (by using symmetries) and if the line picks up a
    ! sign. The extension is found with PointToLineExtensionZ, wich returns a
    ! pointer to the correct variables.
    ! All this is then put into eq. 5.3 in
    !      D. Baye & P. Heenen - J.Phys.A.: Math.Gen. 19 (1986)
    ! in Derlag_1d
    !---------------------------------------------------------------------------
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    real(KIND=dp)                     :: LineExtension(nz)
    integer                           :: i,j, SignExtension

    allocate(Der(nx,ny,nz)) ; Der = 0.0_dp
    ! Z Derivatives
    SignExtension=Signs(3,Parity,Signature,TimeSimplex,Component)
    do j=1,ny
      do i=1,nx
        !LineExtension=PointToLineExtensionZ(                                   &
        !&                           Grid, nz, Parity, TimeSimplex,Signature,i,j)
        LineExtension=ExtZ(Grid,nz,i,j)
        Der(i,j,:)=DerLag_1D(                                                  &
        &                  Grid(i,j,:),nz,LineExtension,SignExtension,LagCoefsZ)
      enddo
    enddo

  end function DerLagZ
  function DerLag_1D(Line,LineSize,LineExtension,SignExtension,Coefs)          &
  &   result(DerLine)
    !---------------------------------------------------------------------------
    ! This function calculates the 1-Dimensional derivative of a line of points
    ! using the Lagrange mesh derivatives.
    ! These derivatives are defined by eq. 5.3 in
    ! D. Baye & P. Heenen - J.Phys.A.: Math.Gen. 19 (1986)
    !---------------------------------------------------------------------------
   integer, intent(in)                    :: LineSize,SignExtension
    real(KIND=dp), intent(in)             :: Line(LineSize)
    real(KIND=dp), intent(in)             :: LineExtension(LineSize)
    real(KIND=dp), intent(in)             :: Coefs(:)
    real(KIND=dp)                         :: DerLine(LineSize)
    integer                               :: SignDistance,i,j
    do i=1, LineSize
      DerLine(i)=0.0_dp
      do j=1,LineSize
        if(i.ne.j) then
          SignDistance=(i-j)/abs(i-j)
          DerLine(i)=DerLine(i)+ SignDistance *Coefs(abs(i-j))*Line(j)
          if(LineExtension(j).ne.0.0_dp ) then
            !Very unelegant, but necessary to evade out-of-bounds errors
            ! on the coefs array when breaking symmetries.
            DerLine(i)=DerLine(i)+ SignExtension*Coefs(i+j-1)*LineExtension(j)
          endif
        else
          if(LineExtension(j).ne.0.0_dp ) then
            DerLine(i)=DerLine(i)+ SignExtension*Coefs(i+j-1)*LineExtension(j)
          endif
        endif
      enddo
    enddo
    return
  end function DerLag_1D

  function Laplacian_Lag(Grid, Parity, Signature, TimeSimplex,Component)       &
                       & result(Lap)
    !---------------------------------------------------------------------------
    ! This function splits the calculation of the Laplacian into 1D problems,
    ! using SecondDerLag.
    ! The structure is completely analogous to the Laplacian subroutine.
    !---------------------------------------------------------------------------
    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Lap(:,:,:)
    real(KIND=dp)                     :: DerX(nx,ny,nz)
    real(KIND=dp)                     :: DerY(nx,ny,nz),DerZ(nx,ny,nz)
    real(KIND=dp)                     :: LineExtensionX(nx), LineExtensionY(ny)
    real(KIND=dp)                     :: LineExtensionZ(nz)
    integer                           :: i,j,k, SignExtension

    allocate(Lap(nx,ny,nz))
    Lap=0.0_dp; DerX=0.0_dp; DerY=0.0_dp; DerZ=0.0_dp

    ! X Derivatives
    SignExtension=Signs(1,Parity,Signature,TimeSimplex,Component)
    do k=1,nz
      do j=1,ny
          LineExtensionX=ExtX(Grid,nx,j,k)
          DerX(:,j,k)=SecondDerLag(Grid(:,j,k), nx,LineExtensionX,SignExtension&
          &                       ,LagLapX,LagLapDiag(1))
      enddo
    enddo

    ! Y Derivatives
    SignExtension=Signs(2,Parity,Signature,TimeSimplex,Component)
    do k=1,nz
      do i=1,nx
          LineExtensionY=ExtY(Grid,ny,i,k)
          DerY(i,:,k)=SecondDerLag(Grid(i,:,k), ny,LineExtensionY,SignExtension&
          &                       ,LagLapY,LagLapDiag(2))
      enddo
    enddo

    !Z Derivatives
    SignExtension=Signs(3,Parity,Signature,TimeSimplex,Component)

    do j=1,ny
      do i=1,nx
          LineExtensionZ=ExtZ(Grid,nz,i,j)
          DerZ(i,j,:)=SecondDerLag(Grid(i,j,:), nz,LineExtensionZ,SignExtension&
                                  ,LagLapZ,LagLapDiag(3))
      enddo
    enddo

    Lap = DerX+DerY+DerZ

    return
  end function Laplacian_Lag

  function SecondDerLag(Line,LineSize,LineExtension, SignExtension             &
  &                    ,Coefficients, Diag)                                    &
  &                     result(DerLine)
    !---------------------------------------------------------------------------
    ! This function computes the second derivative of the Line variable, using
    ! Lagrange mesh derivatives.
    ! For the formulas, see
    ! D.Baye & P-H. Heenen, Generalised Meshes for Quantum Mechanical Problems
    !                J.Phys. A: Meth. Gen 19 (1986)
    !---------------------------------------------------------------------------
    integer, intent(in)        :: LineSize,SignExtension
    real(KIND=dp), intent(in)  :: Line(LineSize),Coefficients(:), Diag
    real(KIND=dp), intent(in)  :: LineExtension(LineSize)
    real(KIND=dp)              :: DerLine(LineSize)
    integer :: i,j,p

    p=SignExtension

    do i=1,LineSize
      DerLine(i)=0._dp
      do j=1,LineSize
        if(i.ne.j) then
          DerLine(i)=DerLine(i)+ Coefficients(abs(i-j))*Line(j)
          if(LineExtension(j).ne.0.0_dp ) then
            DerLine(i)=DerLine(i)+ p*Coefficients(i+j-1)*LineExtension(j)
          endif
        else
          DerLine(i)=DerLine(i)+ Diag*Line(j)
          if(LineExtension(j).ne.0.0_dp ) then
            DerLine(i)=DerLine(i) + p*Coefficients(i+j-1)*LineExtension(j)
          endif
        endif
      enddo
    enddo

  end function SecondDerLag
end module Derivatives
