module GenInfo

    use CompilationInfo

    implicit none

    save
    !---------------------------------------------------------------------------
    !Number of points in every direction and total number of points
    ! Default = 16
    integer, public :: nx=16,ny=16,nz=16, nv=16**3
    !---------------------------------------------------------------------------
    !Number of protons and neutrons in the nucleus
    real(KIND=dp), public :: Neutrons, Protons
    !---------------------------------------------------------------------------
    ! The Symmetries conserved in the code.
    ! TRC = Time Reversal Conservation
    ! TSC = Time Simplex Conservation
    ! PC = Parity Conservation
    ! SC = Signature Conservation
    ! IC = Isospin Conservation
    ! SZC = z-simplex conservation
    !---------------------------------------------------------------------------
    logical, public :: TRC,TSC,PC,SC,IC, SZC
    !---------------------------------------------------------------------------
    ! Integer codes for conservations of these symmetries. Pretty handy for
    ! calling derivative routines.
    integer, public :: SignatureInt, ParityInt, TimeSimplexInt
    !---------------------------------------------------------------------------
    ! Line element and volume element of the box. dx is in fm, dv in fm^3.
    real(KIND=dp), public  :: dx=0.8_dp
    real(KIND=dp), public  :: dv=(0.8_dp**3)
    !---------------------------------------------------------------------------
    !Time step in 10^-22 s
    real(KIND=dp), public  :: dt=0.012_dp, ReadjustTime=0.95_dp
    !---------------------------------------------------------------------------
    !Some useful numbers. Suffix determine whether we are talking about integers
    ! or reals.
    real(KIND=dp), parameter, public  :: ZeroR = 0.0_dp, OneR=1.0_dp,TwoR=2.0_dp
    real(KIND=dp),parameter, public   :: ThreeR=3.0_dp,FourR=4.0_dp,TenR=10.0_dp
    real(KIND=dp),parameter, public   :: prec4=0.1E-03_dp,prec8=0.1E-07_dp
    real(KIND=dp),parameter, public   :: prec13=0.1E-12_dp
    integer, parameter, public        :: ZeroI=0,OneI=1,TwoI=2,ThreeI=3,FourI=4
    real(KIND=dp), parameter, public  :: pi=4.0_dp*atan2(1.0_dp,1.0_dp)
    !---------------------------------------------------------------------------
    ! Determines if MOCCa runs in Testmodus. Default = 0
    integer, public :: TestRun=0
    !---------------------------------------------------------------------------
    !Maximum number of iterations and number of iterations to skip printing of
    ! the code in the evolve Subroutine
    integer, public :: MaxIter=0, PrintIter=0
    !---------------------------------------------------------------------------
    !Convergence parameters.
    real(KIND=dp),public  :: MomentPrec=1d-4, EnergyPrec=0.01E-08
    real(KIND=dp),public  :: PairingPrec=1d-4, CrankPrec=1d-4
    !---------------------------------------------------------------------------
    !Channel for the transmitting of errors.
    integer, public :: ErrorChannel
    character(len=256), public :: ErrorFileName='MOCCa.error'
    !---------------------------------------------------------------------------
    !Integer governing the Taylor expansion of exp(-dt/hbar * h)
    integer, public :: TaylorOrder=1
    !---------------------------------------------------------------------------
    !Integer governing the extent of Broyden Mixing
    integer, public :: BroydenOrder=1, BroydenSize
    !Mixing coefficient of Broyden
    real(Kind=dp)   :: BroydenAlpha=0.5_dp
    !---------------------------------------------------------------------------
    !Flags for the Damping procedures
    logical, public       :: InverseKineticDamping=.false.
    real(KIND=dp), public :: E0 =0.0_dp
    !---------------------------------------------------------------------------
    ! Determining the calculation of the propagationfactor.
    integer, public :: TimeStepType=0
    !---------------------------------------------------------------------------
    ! Logical, determining whether or not to make search directions conjugate
    logical, public :: CG=.false.
    !---------------------------------------------------------------------------
    ! real: determines whether or not to use the momentum
    real(KIND=dp), public :: Momentum=0.0_dp
    !---------------------------------------------------------------------------
    ! Type of iterative procedure to employ.
    ! IMTS = Imaginary time-step method
    ! Nest = Nesterov optimal descent
    character(len=20):: IterType='ImTS'
    !---------------------------------------------------------------------------
    ! Number of iterations after which to restart the Nesterov-memory
    integer :: Restart = 25
    !---------------------------------------------------------------------------
    ! Whether or not to recalculate the findings with the first implementation 
    ! of the N2LO functional of Becker, Meyer & Davesne.
    logical :: recalcN2LO = .false.
    
    interface Stp
        module procedure StpReal
        module procedure StpInteger
        module procedure StpIntegerReal
    end interface Stp

contains
  subroutine StpReal (Message, NameReal,PrintReal, NameReal2, PrintReal2)
    !---------------------------------------------------------------------------
    ! This subroutine can stop the program, writing a (hopefully) useful message
    ! to the output. The subroutine is also capable of printing the value of a
    ! real, so that the reason for program exit becomes somewhat clearer.
    ! (If the values are explained in
    ! the string "message", that is!)
    !
    ! It is overloaded with StpInteger to "Stp".
    !---------------------------------------------------------------------------
    character (len = *), intent(in)          :: Message
    character (len = *), intent(in),optional :: NameReal,NameReal2
    real(KIND=dp), intent(in), optional      :: PrintReal, PrintReal2

    write (ErrorChannel,*) ,' '
    write (ErrorChannel,*) ' =====  S T O P  ===== '
    write (ErrorChannel,*) ' '
    write (ErrorChannel,*) ' ',message
    write (ErrorChannel,*) ' '

    if(present(NameReal)) then
      write (ErrorChannel,*) , NameReal, PrintReal
      write (ErrorChannel,*) , ""
    endif
    if(present(NameReal2)) then
      write (ErrorChannel,*) , NameReal2, PrintReal2
      write (ErrorChannel,*) , ""
    endif

    close(ErrorChannel)

    stop
  end subroutine StpReal

  subroutine StpInteger(Message,NameInt1,PrintInt1,NameInt2,PrintInt2)
    !---------------------------------------------------------------------------
    ! This subroutine is the same as StpReal. The only difference is that this
    ! one prints up to two integer values, instead of a real.
    ! It is overloaded with StpReal to "Stp".
    !---------------------------------------------------------------------------
    character (len = *), intent(in)         :: Message, NameInt1
    integer, intent(in)                     :: PrintInt1
    character (len = *),intent(in),optional :: NameInt2
    integer, intent(in), optional           :: PrintInt2
    write (ErrorChannel,*),' '
    write (ErrorChannel,*),' =====  S T O P  ===== '
    write (ErrorChannel,*),' '
    write (ErrorChannel,*),' ',message
    write (ErrorChannel,*),' '
    write (ErrorChannel,*), NameInt1, PrintInt1
    write (ErrorChannel,*), ""

    if(present(PrintInt2)) then
      write (ErrorChannel,*), NameInt2, PrintInt2
      write (ErrorChannel,*), ""
    endif

    close (ErrorChannel)

    stop
  end subroutine StpInteger

  subroutine StpIntegerReal(Message,NameInt,PrintInt,NameReal,PrintReal)
    !---------------------------------------------------------------------------
    ! This subroutine is the same as StpReal. The only difference is that this
    ! one prints up to two integer values, instead of a real.
    ! It is overloaded with StpReal to "Stp".
    !---------------------------------------------------------------------------
    character (len = *), intent(in)         :: Message, NameInt
    integer, intent(in)                     :: PrintInt
    character (len = *),intent(in)          :: NameReal
    real(KIND=dp), intent(in)               :: PrintReal
    write (ErrorChannel,*),' '
    write (ErrorChannel,*),' =====  S T O P  ===== '
    write (ErrorChannel,*),' '
    write (ErrorChannel,*),' ',message
    write (ErrorChannel,*),' '
    write (ErrorChannel,*), NameInt, PrintInt
    write (ErrorChannel,*), ""

    write (ErrorChannel,*), NameReal, PrintReal
    write (ErrorChannel,*), ""

    close (ErrorChannel)

    stop
  end subroutine StpIntegerReal

  pure integer function LeviCivita(i,j,k)
    !---------------------------------------------------------------------------
    ! This function is a quick & dirty implementation of the LeviCivita symbol
    ! epsilon_{ijk}
    !---------------------------------------------------------------------------
    integer, intent(in) :: i,j,k

    if((i.eq.j).or.(j.eq.k).or.(k.eq.i)) then
        LeviCivita=0

    elseif(((i.eq.1).and.(j.eq.2).and.(k.eq.3)) &
     & .or.((i.eq.3).and.(j.eq.1).and.(k.eq.2)) &
     & .or.((i.eq.2).and.(j.eq.3).and.(k.eq.1))) then
        LeviCivita=1
    else
        LeviCivita=-1
    endif

    return
  end function LeviCivita

  subroutine SetSymmetryInteger
    !---------------------------------------------------------------------------
    ! Subroutine that sets ParityInt, SignatureInt and TimeSimplexInt.
    !---------------------------------------------------------------------------
    if(PC) then
            ParityInt      = 1
    else
            ParityInt      = 0
    endif
    if(SC) then
            SignatureInt   = 1
    else
            SignatureInt   = 0
    endif
    if(TSC) then
            TimeSimplexInt = 1
    else
            TimeSimplexInt = 0
    endif
  end subroutine SetSymmetryInteger

  subroutine AssignSymmetries(TimeReversal,Parity,Signature,TimeSimplex,Isospin,zsimplex)
    !---------------------------------------------------------------------------
    ! Subroutine that assigns the correct values to the following logicals:
    !       TRC,PC,SC,TSC and IC.
    ! The number of spatial Symmetries is also counted and dv is multiplied by
    !       2**(NumSpatSym)
    ! so that MOCCa uses the correct volme element.
    !
    ! The subroutine stops the program is the input is not understood.
    !---------------------------------------------------------------------------
    ! At the moment, the program also stops when T=0 pairing is indicated.
    !---------------------------------------------------------------------------

    integer,intent(in) :: TimeReversal,Parity,Signature,TimeSimplex,Isospin, zsimplex
    integer            :: NumSpatSym = 0

    ! Converting input to logicals, and stopping if input isn't clear
    if(TimeReversal.eq.1) then
        TRC = .true.
    elseif(TimeReversal.eq.0) then
        TRC = .false.
    else
      call stp(&
      & "MOCCa doesn't know if TimeReversal should be conserved.",             &
      "TimeReversal",TimeReversal)
    endif

    if(Parity.eq.1) then
            PC = .true.
            NumSpatSym = NumSpatSym + 1
    elseif(Parity.eq.0) then
            PC = .false.
    else
            call stp("MOCCa doesn't know if Parity should be conserved.",      &
            &        "Parity", Parity)
    endif

    if(Signature.eq.1) then
            SC = .true.
            NumSpatSym = NumSpatSym + 1
    elseif(Signature.eq.0) then
            SC = .false.
    else
            call stp("MOCCa doesn't know if Signature should be conserved.",   &
            &        "Signature", Signature)
    endif

!    if((zsimplex.eq.1) .and. ((parity.ne.1) .and. (signature.ne.1))) then
!            SZC = .true.
!            NumSpatSym = NumSpatSym + 1
!    elseif(zsimplex.eq.0) then
!            SZC = .false.
!    endif

    if(TimeSimplex.eq.1) then
            TSC = .true.
             NumSpatSym = NumSpatSym + 1
    elseif(TimeSimplex.eq.0) then
            TSC = .false.
    else
            call stp("MOCCa doesn't know if Time Simplex should be conserved.",&
            &        "Time Simplex", TimeSimplex)
    endif

    if(Isospin.eq.1) then
            IC = .true.
    elseif(Isospin.eq.0) then
            IC = .false.
            call stp("Breaking isospin symmetry is not yet option.",           &
            & "Isospin.", Isospin)
    else
            call stp("MOCCa doesn't know if Isospin should be conserved.",     &
            & "Isospin", Isospin)
    endif

    !Multiply dv by the appropriate power of 2
    dv = dv* 2.0_dp**(NumSpatSym)

    return
  end subroutine AssignSymmetries

  subroutine ReadGenInfo
    !---------------------------------------------------------------------------
    ! Subroutine that reads the user input for the variables in the GenInfo
    ! module.
    !---------------------------------------------------------------------------
    integer :: Parity=1, Isospin=1, Signature=1, TimeSimplex=1,TimeReversal=1
    integer :: zsimplex=1

    NameList /GenInfo/  TimeReversal, Parity, Signature, TimeSimplex, Isospin, &
                     & MaxIter, dt, TestRun, zsimplex,                         &
                     &  Neutrons, Protons, ErrorFileName, nx, ny,nz, dx,       &
                     &  PrintIter, MomentPrec, EnergyPrec, TaylorOrder,        &
                     &  InverseKineticDamping, E0, PairingPrec, CrankPrec,     &
                     &  IterType, Restart, recalcN2LO, momentum

    !Reading the symmetries that are conserved or broken.
    read (unit=*, nml=GenInfo)

    !Assigning the correct value to dv
    dv=dx**3

    !Reserve a channel for the output of errors.
    call get_unit(ErrorChannel)
    open(ErrorChannel,form="formatted", file=ErrorFileName)

    !Checking the number of points on the mesh
    if(nx.le.0) &
    & call stp('Number of points on the mesh is not valid.', 'nx', nx)
    if(ny.le.0) &
    & call stp('Number of points on the mesh is not valid.', 'ny', ny)
    if(nz.le.0) &
    & call stp('Number of points on the mesh is not valid.', 'nz', nz)

    !Checking the mesh distance
    if(dx.le.ZeroR) &
    & call stp('Mesh distance should not be zero or negative.', 'dx', dx)

    !Checking the number of protons and neutrons.
    if(Neutrons.lt.0) &
    & call stp('Number of neutrons is not valid.', 'Neutrons', Neutrons)
    if(Protons .lt.0) &
    & call stp('Number of protons is not valid.', 'Protons', Protons)

    !Checking the Symmetry parameters
    if(TimeReversal.lt.0 .or. TimeReversal.gt.1) &
    & call stp('Wrong value for TimeReversal', 'Time Reversal', TimeReversal)
    if(Parity.lt.0 .or. Parity.gt.1) &
    & call stp('Wrong value for Parity', 'Parity', Parity)
    if(Signature.lt.0 .or. Signature.gt.1) &
    & call stp('Wrong value for Signature', 'Signature', Signature)
    if(TimeSimplex.lt.0 .or. TimeSimplex.gt.1) &
    & call stp('Wrong value for TimeSimplex', 'TimeSimplex', TimeSimplex)
    if(Isospin.ne.1) &
    call stp('Wrong value for Isospin', 'Isospin', Isospin)

    !The number of iterations should be positive
    if(MaxIter.lt.0) &
    & call stp('Number of Iterations should not be negative',                  &
    &          'Iterations', MaxIter)

    ! There is an approximate relation between dt and dx.
    ! See the discussion in
    ! P. Bonche et al, Solution of the Skyrme HF + BCS equation on a 3D mesh ,
    ! Computer Physics Communications 171 (2005) 49–62

    ! Assigning the symmetries used.
    call AssignSymmetries(TimeReversal,Parity,Signature,TimeSimplex,Isospin, zsimplex)
    call SetSymmetryInteger

    !If printiter is not set, set it to MaxIter/10
    if(Printiter.eq.0) PrintIter = MaxIter/10
    PrintIter = max(PrintIter, 1)

    ! Make Itertype uppercase for comparison purposes.
    call to_upper(IterType, IterType)

  end subroutine ReadGenInfo

  subroutine get_unit ( iunit )
    !***************************************************************************
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 September 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !---------------------------------------------------------------------------
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ios
    integer ( kind = 4 ) iunit
    logical lopen

    iunit = 0

    do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
                iunit = i
                return
        end if
      end if
    end if

    end do

    return
  end subroutine get_unit

  function Reverse(Line)
    !---------------------------------------------------------------------------
    ! A small function that reverses the values on a line.
    !---------------------------------------------------------------------------

    real(Kind=dp), intent(in) :: Line(:)
    real(KIND=dp)             :: Reverse(size(Line))
    integer                   :: i,n

    n= size(Line)
    do i=1, n
            Reverse(i) = Line(n - i  + 1)
    enddo
  end function Reverse

  subroutine to_upper (str, string)
    !---------------------------------------------------------------------------
    ! Subroutine that changes a string to uppercase.
    !---------------------------------------------------------------------------
    character(*)        :: str
    character(len(str)) :: string

    Integer :: ic, i
    !Ugly but effective and independent of platform and implementation.
    character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    if(len(str) .ne. len(string)) then
    call stp('Strings of different length in to_upper')
    endif
    string = str
    do i = 1, len_trim(str)
    ic = INDEX(low, str(i:i)) !Note that ic = 0 when substring is not found
    if (ic > 0) then
    string(i:i) = cap(ic:ic)
    else
    string(i:i) = str(i:i)
    endif
    end do

  end subroutine to_upper
end module GenInfo
