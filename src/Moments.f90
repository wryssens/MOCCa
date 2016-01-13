module Moments
!-------------------------------------------------------------------------------
! This module contains all the multipole moments of the nuclear wavefunction and
! the subroutines to calculate them. ! The representation of the spherical 
! harmonics used is the one from Messiah, using real spherical harmonics.
!
! Thus the spherical harmonics are, in terms of spherical coordinates:
!
!        Y_{l}^{m} = (-1)^{m}[(2*l+1)/4\pi (l-m)!/(l+m)!]^{1/2} 
!                 *    P^{m}_{l}(cos(\theta)) e^{i m \phi}
!
!  where the P^{m}_{l} are the associated Legendre Polynomials and -l <= m <= l.
!  The multipole moments can thus be written as:
!      \langle \hat{Q}_{lm} \rangle = \langle r^l Y_{lm} \rangle
!
!  The code "talks" about the real and imaginary parts of these moments. 
!  The only relevant values of m are thus positive.
!
! Note that one is also able to constraint the total size of moments, thereby
! leaving all the `angles' free. 
! We define the total size as:
!
! Q_l = 4 \sqrt{pi}/\sqrt{2*l+1} *                                              
!        \sqrt{ \Sum_{m = 0,..,l} [ < Re Q_{lm} >**2 + < Im Q_{lm}**2 > ] }
!
! For the l = 2 multipole moments, this reduces to Q0 (see EV8 article), leaving
! gamma angle free. For higher l, there are obviously more angles.
! Internally, the code constrains Q_l^2 so that we don't have to deal with 
! square roots.
!
! In practice, the code constrains Q_l**2, since that means we don't have to 
! worry about the square root.
!
! Deformation parameters Beta_lm are also calculated and printed. Their 
! definition is as follows:
!  Beta_lm = 4 * pi /( 3 * R_0**l * A) < Q_{lm} >
!  where R0 = 1.2A**(1/3)
!
! For more conversion formulae see:
!   
!   W.Ryssens, V. Hellemans, M. Bender & P.-H. Heenen, CPC 187 (2015) 175-194
!
!
!-------------------------------------------------------------------------------
!  The contribution of constraints to the single-particle hamiltonian is
!  collected in  the ConstraintEnergy array.
!  Several different constraints are possible:
!  - Quadratic constraints (as in EV8 and the articles)
!  - K.Rutz Constraints as described in his PhD Thesis and in
!    R.Y. Cusson et al., J.Phys.A 320, 475-482
!
!  The Rutz constraints are recommended and default.
!
!-------------------------------------------------------------------------------
!    Symmetries that forbid certain multipole moments:
!        Time Simplex constrains the imaginary parts of every moment.
!        Signature constrains every moment with odd m.
!        Parity constrains every moment with odd l.
!-------------------------------------------------------------------------------

  use GenInfo
  use CompilationInfo    
  use Mesh

  implicit none

  save

  !-----------------------------------------------------------------------------  
  type Moment 
      !-------------------------------------------------------------------------
      !Principal and magnetic quantum number of the multipole moment
      integer :: l,m
      !Looking at the real or imaginary part
      logical :: ImPart
      !-------------------------------------------------------------------------
      !The SpherHarm variable contains the values of the associated Spherical 
      !Harmonic at each point of the mesh, 
      !multiplied by r^l and sqrt{4pi/{2l+1}}
      !-------------------------------------------------------------------------
      real(KIND=dp), allocatable  :: SpherHarm(:,:,:)
      !-------------------------------------------------------------------------
      ! Integer controlling the type of constraint on the multipole moment. 
      ! 0 - unconstrained
      ! 1 - Quadratic
      ! 2 - Rutz-type
      !-------------------------------------------------------------------------
      integer        :: ConstraintType 
      !-------------------------------------------------------------------------
      ! Calculated value of the moment with and without cutoff, 
      ! for neutrons and protons.
      !-------------------------------------------------------------------------
      real(KIND=dp)  :: Value(2), ValueCut(2)
      !-------------------------------------------------------------------------
      ! Calculated value of the moment SQUARED.
      ! So < Q_{lm}^2 >
      !------------------------------------------------------------------------- 
      real(KIND=dp)  :: Squared(2), SquaredCut(2) 
      !-------------------------------------------------------------------------
      !The program retains the previous 7 values of the moments
      !for checking of convergence.
      !-------------------------------------------------------------------------
      real(KIND=dp)  :: OldValue(2,7), OldValueCut(2,7) 
      !-------------------------------------------------------------------------
      ! Value of the constraint on the multipole moment and the intensity of 
      ! the constraint (this intensity is the Lagrange multiplier in case
      ! of normal and Rutz Constraints.)
      ! The value is the value for the `total' multipole moment
      real(KIND=dp), allocatable :: Intensity(:)
      real(KIND=dp), allocatable :: Constraint(:)
      !-------------------------------------------------------------------------
      ! True value of the constraint. This means that the real numbers
      ! Constraint are subject to readjustment, but TrueConstraint will always
      ! contain the actual constraint value demanded by the user.
      !-------------------------------------------------------------------------
      real(Kind=dp), allocatable :: TrueConstraint(:)             
      !-------------------------------------------------------------------------
      ! Logical whether the `size' of the total multipole moment for a given l 
      ! is constrained or not. In that case the values for constraint take
      ! the desired value for the total moment.
      logical                    :: Total = .false.
      !-------------------------------------------------------------------------
      ! Integer switch to determine whether we place the constraint on
      ! 
      ! 1) Total value of the multipole moment  
      ! 2) Values for neutrons & protons separately
      ! 3) Difference of proton & neutron values
      !-------------------------------------------------------------------------
      integer                    :: Isoswitch
      !-------------------------------------------------------------------------
      ! Pointers to the previous and next item in the linked list.
      !-------------------------------------------------------------------------
      type(Moment), pointer      :: Prev
      type(Moment), pointer      :: Next
      !-------------------------------------------------------------------------
      ! Deformation parameter Beta_lm associated with the moment.
      !-------------------------------------------------------------------------
      real(KIND=dp)         :: Beta(3)
    contains
      procedure, pass, public :: WriteMoment
      generic :: Write=> WriteMoment
      
  end type Moment
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !Maximum degree of the multipole components that are considered. Default = 2
  integer, public :: MaxMoment=6
  !-----------------------------------------------------------------------------
  ! Starting point for the linked list of moments.
  !-----------------------------------------------------------------------------
  type(Moment), public, pointer :: Root
  !-----------------------------------------------------------------------------
  !CutoffType
  ! 0 : Rutz-Type
  ! 1 : Spherical cut-off
  !-----------------------------------------------------------------------------
  integer, public :: CutoffType=0
  !-----------------------------------------------------------------------------
  !Parameters of the cutoff function.
  ! 1 + exp[(\DeltaR - radd)/acut]
  ! Default values are the ones from K. Rutz et al, Nucl. Phys. A590 (1995) 690
  !-----------------------------------------------------------------------------
  real(KIND=dp), public :: radd = 4._dp,acut=0.4_dp
  !-----------------------------------------------------------------------------
  !Cutoff function
  ! For 1) Neutrons, 2) Protons 
  real(KIND=dp), allocatable, public :: Cutoff(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! Contribution to the single-particle Hamiltonian associated with the 
  ! constraints is stored in this variable
  !-----------------------------------------------------------------------------
  real(Kind=dp), public, allocatable :: ConstraintEnergy(:,:,:,:)
  !-----------------------------------------------------------------------------
  !Damping associated with the calculation of <Q_lm>
  !-----------------------------------------------------------------------------
  real(Kind=dp)  :: Damping=0.90_dp
  !-----------------------------------------------------------------------------
  !Damping associated with the readjustment of quadratic constraints
  !-----------------------------------------------------------------------------
  real(KIND=dp)  :: ReadjustSlowDown=0.01_dp
  !-----------------------------------------------------------------------------
  !Pointer to the cutoff procedure chosen. 
  !-----------------------------------------------------------------------------
  procedure(RutzCutOff), pointer :: CompCutoff
  !-----------------------------------------------------------------------------
  ! Quadrupole moments in different representations:
  !
  ! Cartesian: QCartesian
  ! (Q,gamma): Q, G
  ! iq1, iq2 : iq1, iq2
  !-----------------------------------------------------------------------------
  real(kind=DP),public :: QCartesian(3,3,3), Q(3), G(3), iq1(3), iq2(3)
  !-----------------------------------------------------------------------------
  ! Logical determining whether or not to take momentconstraints from
  ! the wavefunction file.
  logical, public      :: ContinueMoment=.false.
  !-----------------------------------------------------------------------------
  ! Numerical parameters for the Rutz-constraint update scheme
  real(KIND=dp) :: c0=0.2_dp, d0=0.01_dp, epsilon=7.0_dp

  logical :: RutzToQuadratic = .false.

contains

  recursive function FindMoment(l,m,Impart, StartMoment) result(FoundMoment)
    !---------------------------------------------------------------------------
    ! Subroutine that finds a pointer to the moment, specified by multipole 
    ! values l,m and whether it is a real or imaginary part.
    ! For maximum efficiency, it is possible to specify a StartMoment, so that 
    ! the search starts at that specific moment. 
    ! Ex: MOCCa can find Im(Q_{2,2}) faster if it starts from Re(Q_{20}).
    !
    ! The result is null if the moment isn't found.
    !---------------------------------------------------------------------------
    integer, intent(in)                        :: l,m
    logical, intent(in)                        :: Impart
    type(Moment), intent(in), pointer, optional:: StartMoment        
    type(Moment), pointer                      :: Current, FoundMoment

    nullify(Current); nullify(FoundMoment)

    if(present(StartMoment)) then
      Current => StartMoment
    else
      Current => Root
    endif        

    do while(associated(Current%Next))
      Current => Current%Next           
      if((Current%l.eq.l) .and. (Current%m .eq. m) .and.                       &
      &               (Impart.eqv.Current%Impart)) then
        !Then we've found the correct moment.
        FoundMoment => Current
        return
      endif
    enddo
    if(present(StartMoment)) then         
      ! If this part of the code gets reached, we picked a bad starting point.
      ! So we start over again! 
      FoundMoment => FindMoment(l,m,Impart)   
    else
      FoundMoment => null()
    endif
    return
  end function FindMoment
  
  subroutine IniMoments()    
    !---------------------------------------------------------------------------
    ! This routine sets up the linked list of moments that MOCCa needs.
    ! 1) The first moment ( l=0, m=0 ) is created. It is called 'Root'.
    !    The expected value of this moment are obviously the number of particles
    ! 2) The values of all the desired spherical harmonics 
    !    (l=0,..., MaxMoment, m=0,...,MaxMoment) are computed.
    ! 3) The routine creates all the necessary moments using subroutine 
    !    NewMoment and links them in the following order:
    !    (l,m, R/I): (0,0), (1,0), (1,1,R), (1,1,I), (2,0), (2,1,R), (2,1,I),...
    !    Moments that are automatically zero (by symmetry) are not created by 
    !    NewMoment and are absent from the list.
    ! 4) MOCCa detects the moments that do not represent physical degrees of 
    !    freedom and initialises constraints on them.
    ! 5) After all this, a final addition to the linked list is made: the 
    !    radius squared. (This can be used to calculate the rms radius 
    !    obviously).
    !
    ! --------------------------------------------------------------------------
    integer :: l,m,ImPart
    ! Meshes containing the values of all associated Legendre Polynomials
    ! Some dimensions start at index 0, because l and m can both be 0
    real(KIND=dp)  :: SpherHarmMesh(nx,ny,nz,0:MaxMoment,0:MaxMoment,2)

    type(Moment),pointer :: NextMoment,Current
    ! This subroutine creates all the needed multipole moment objects.

    !Create the first Multipole moment: l=0, m=0
    nullify(Root) ;    allocate(Root)            
    Root%l=0      ;    Root%m=0
    allocate(Root%SpherHarm(nx,ny,nz))
    
    !This is the l=0,m=0 spherical harmonic
    Root%SpherHarm=1.0_dp/sqrt(4.0_dp*pi)
    Root%Impart=.false.
    Root%ConstraintType=0
    nullify(Root%Prev) ;  nullify(Root%Next)

    nullify(Current)   ;  allocate(Current)
    Current=>Root
    nullify(Current%Next) ;  nullify(Current%Prev) ; nullify(NextMoment)

    ! Finding the Spherical Harmonics       
    SpherHarmMesh=SphericalHarmonics(Mesh3D,nx,ny,nz)

    !Creating all the moments and assigning each moment the spherical harmonic
    do l=1, MaxMoment
      do m=0,l
        do ImPart=0,1                        
          NextMoment => NewMoment(l,m,ImPart)

          !Placing the moment in the list
          if(associated(NextMoment)) then
            if(Impart.eq.1.) then
              NextMoment%SpherHarm=SpherHarmMesh(:,:,:,l,m,2)
            else
              NextMoment%SpherHarm=SpherHarmMesh(:,:,:,l,m,1)
            endif

            NextMoment%Prev => Current
            Current%Next => NextMoment
            Current => NextMoment
            nullify(NextMoment)            
          endif
        enddo
      enddo
    enddo
    
    !Appending the radius squared...
    NextMoment   => NewMoment(-2,0,0)
    NextMoment%SpherHarm = sum(Mesh3D**2,1)
    
    Current%Next => NextMoment
    NextMoment%Prev => Current
    
    nullify(Current); nullify(NextMoment)
    return
  end subroutine IniMoments
  
  function SphericalHarmonics(Mesh,mx,my,mz) result(SpherHarmMesh)
    !---------------------------------------------------------------------------
    ! This function computes the values of all the spherical harmonics up to 
    ! l=Maxmoment for the coordinates in Mesh. The input Mesh is thus an array 
    ! of size (3, mx,my,mz) where the first column contains the coordinates of 
    ! the points and mx,my and mz describe the number of points where we want to
    ! calculate the spherical harmonics.
    !---------------------------------------------------------------------------
    ! The formula used is the one from Messiah:
    !    Y^{m}_{l} = (-1)^m [(2*l+1)/(4*\pi) (l-m)!/(l+m)!]^{1/2}
    !              *  P^m_{l}[cos(\theta)] e^{im\phi}
    ! where P^m_l{x} are the associated Legendre Polynomials.
    ! Note that the variable SpherHarm does not store the real and imaginary 
    ! parts of this expression, but rather
    !     r^L  Y^{m}_{l}
    !---------------------------------------------------------------------------
    
    integer, intent(in)        :: mx,my,mz
    real(KIND=dp), intent(in)  :: Mesh(3,mx,my,mz)
    
    real(KIND=dp)  :: LegendreMesh(mx,my,mz,0:MaxMoment,0:MaxMoment)
    real(KIND=dp)  :: SpherHarmMesh(mx,my,mz,0:MaxMoment,0:MaxMoment,2)
    !r, Sin(\theta),Cos(\theta),\phi, sin(m*\Phi) and cos(m*\Phi)
    real(KIND=dp)  :: r,cosTheta,sinTheta, phi,sinmPhi, cosmPhi
    real(KIND=dp)  :: fac, X,Y,Z
    integer        :: q,i,j,k,l,m
    real(KIND=dp)  :: factorialquotient  

    !LegendreMesh = 0.0_dp
    
    do k=1,mz
      do j=1,my
        do i=1,mx
          X = Mesh(1,i,j,k)
          Y = Mesh(2,i,j,k)
          Z = Mesh(3,i,j,k)
          
          !Finding the corresponding spherical moments
          r=sqrt(X**2+Y**2+Z**2)
          sinTheta=sqrt(X**2+Y**2)/r
          cosTheta=Z/r
          phi=atan2(Y,X)

          !Calculating the values of the spherical harmonics
          do m=0,MaxMoment
            ! Finding the associated Legendre Polynomials
            LegendreMesh(i,j,k,:,m)=legendre_pnm(MaxMoment, m, cosTheta)
            sinmPhi=sin(m*Phi)
            cosmPhi=cos(m*Phi)

            do l=m,MaxMoment ! Note that l>= m
              !Computing the factorials in the spherical harmonics using the 
              !product intrinsic.
              
              factorialquotient=1.0_dp
              do q=1,2*m
                factorialquotient = factorialquotient/( l + m - q + 1)
              enddo
                                      
              fac=(-1)**m*sqrt(factorialquotient)*r**l
              fac = fac * sqrt((2*l +1)/(4 * pi))             
              ! Real part of Y_{lm}
              SpherHarmMesh(i,j,k,l,m,1)=fac*LegendreMesh(i,j,k,l,m)*cosmPhi
              ! Imaginary part of Y_{lm}
              SpherHarmMesh(i,j,k,l,m,2)=fac*LegendreMesh(i,j,k,l,m)*sinmPhi           
            enddo
          enddo        
        enddo
      enddo
    enddo
    
    return
  end function SphericalHarmonics
  
  pure function NewMoment(l,m,ImPart)
    !---------------------------------------------------------------------------
    ! This subroutine allocates the necessary memory space for the multipole 
    ! moments. That is, if the symmetries of the problem allow such a moment to 
    ! be nonzero. If this is not the case, the moment is not created.
    !---------------------------------------------------------------------------
    type(Moment),pointer   :: NewMoment
    integer, intent(in)    :: l,m,ImPart
    
    nullify(NewMoment)

    if((m.eq.0).and.(Impart.eq.1)) then
      !Do not create a new moment for an imaginary part of moments with m=0
      return
    endif            

    if((mod(l,2).ne.0).and.PC) then
      !The multipole moments of odd order are exactly zero 
      !when Parity is conserved.
      return
    endif

    if((ImPart.eq.1.).and.TSC) then
      ! All Imaginary parts of the multipole moments are exactly zero when time 
      ! simplex is conserved. This can be seen by noting that phi = atan(X/Y) 
      ! and thus phi is odd under the action of timesimplex. Thus, everything 
      ! involving the sin of phi is odd under the action of time simplex.
      return
    endif

    if((mod(m,2).ne.0).and.SC) then
        ! Odd magnetic quantum numbers are prohibited by Signature Conservation
        ! This can be seen as follows:
        ! If x>0 and y>0, \Phi= atan(y/x)
        ! Acting with R_z gives x<0 and y<0, in which case phi=atan(y/x) - \pi
        ! Analogously,if x>0 and y<0, \Phi = atan(y/x)
        ! Then R_z gives \phi=atan(y/x) + \pi
        ! Thus under R_z \phi => \phi +- \pi
        ! and thus e^{im\phi} is only invariant when m is even.
        return
    endif        
    allocate(NewMoment)
    NewMoment%l=l
    NewMoment%m=m
    
    nullify(NewMoment%Prev, NewMoment%Next)
    allocate(NewMoment%SpherHarm(nx,ny,nz))
    
    if(ImPart.eq.1.) then
        NewMoment%ImPart=.true.
    else
        NewMoment%ImPart=.false.
    endif
    
    ! The parameters of the new multipole moment are set to zero by default.
    NewMoment%ConstraintType= 0
    NewMoment%Value         = 0.0_dp
    NewMoment%ValueCut      = 0.0_dp
    NewMoment%SpherHarm     = 0.0_dp
    NewMoment%OldValue      = 0.0_dp
    NewMoment%OldValueCut   = 0.0_dp
    NewMoment%Isoswitch     = 1
    NewMoment%Squared       = 0.0_dp
    NewMoment%SquaredCut    = 0.0_dp 
    
    return
  end function NewMoment
  
  subroutine PrintMomentInfo()
    !--------------------------------------------------------------------
    ! Subroutine that prints all the info on the moments that are used.
    !
    !--------------------------------------------------------------------

    1 format (21('-'), ' Multipole Moments', 21('-') )
    2 format ('Maximum l considered ' , i3 ) 
    3 format ('Constraints parameters: ')
    4 format ('  radd = ', e9.2 , '  acut = ', e9.2)
    5 format ('  Damping of constraint calculation  = ', f6.3)
    6 format ('  Damping of constraint readjustment = ', f6.3)
    7 format (/, 'Constraints obtained from file.')
    8 format (/, 'Constraints obtained from data.')
    9 format ('  Rutz Constraints Epsilon : ' f6.3)    
     
    print 1
    print 2 , MaxMoment
    print 3
    print 4, radd, acut
    print 5, Damping
    print 6, ReadjustSlowdown
    print 9, epsilon
    if(ContinueMoment) then
      print 7
    else
      print 8
    endif

    call PrintConstraints

  end subroutine PrintMomentInfo

  subroutine PrintConstraints
    !---------------------------------------------------------------------------
    ! A subroutine that prints info on the constraints of multipole moments.
    !---------------------------------------------------------------------------
    1    format (/,11("_"),'Constraints on Q_{lm} (neutron + proton)', 11("_"))
    11   format (/,11("_"),'Constraints on Q_{lm} (neutron,  proton)', 11("_"))
    111  format (/,11("_"),'Constraints on Q_{lm} (neutron - proton)', 11("_"))
    1111 format (/,11("_"),'Constraints on Q_{l } (neutron + proton)', 11("_"))
    
    2    format (17x,"Total      Type")
    21   format (17x,"Neutron    Proton      Type")    
    211  format (17x,"Difference Type")
    2111 format (17x,"Total      Type")
    
    4   format ('RMS Radius ', 1(1x,f10.2), 7x, A9)
    41  format ('RMS Radius ', 2(1x,f10.2), 7x, A9)
    411 format ('RMS Radius ', 1(1x,f10.2), 7x, A9)
    
    5    format (A2, " Q_{", 2i2,"}", 1(1x,f10.2), 6x, a9)
    51   format (A2, " Q_{", 2i2,"}", 2(1x,f10.2), 6x, a9)
    511  format (A2, " Q_{", 2i2,"}", 1(1x,f10.2), 6x, a9)
    5111 format (A2, " Q_{",  i2,"}", 2x, 1(1x,f10.2), 6x, a9)
    
    type(Moment), pointer :: Current
    character(len=2)      :: ReIm
    character(len=9)      :: ConType
    logical               :: PrintWorthy(4)
    
    !---------------------------------------------------------------------------
    ! First check what kind of constraints are present.
    Current => Root    
    PrintWorthy=.false.
    !First check if there is something worth printing, else return.
    do while(associated(Current%Next))
      Current => Current%Next  
      if(Current%ConstraintType.ne.0) then
        if(.not. Current%Total) then
          select case(Current%Isoswitch)
          case(1)
            PrintWorthy(1) = .true.
          case(2)
            PrintWorthy(2) = .true.
          case(3)
            PrintWorthy(3) = .true.
          end select
        else
          PrintWorthy(4) = .true.
        endif
      endif
    enddo 
    !---------------------------------------------------------------------------
    !Doing the constraints on proton+neutron multipole moments
    if (PrintWorthy(1)) then
      Current => Root   
      print 1
      print 2

      !Print the moment specific Values
      do while(associated(Current%Next))
        Current => Current%Next
        
        if(Current%Isoswitch.ne.1 .or. Current%ConstraintType.eq.0) cycle
        if(Current%Total) cycle
        
        select case (Current%ConstraintType)             
        case(1)
            ConType='Quadratic'
        case(2)
          ConType='K.Rutz   '             
        end select
        
        select case(Current%l)
        case (-2)
          !Printing the constraint on rms; some extra manipulation since
          ! the internal housekeeping is done with the true radius**2        
          print 4, sqrt(Current%TrueConstraint/(Protons+Neutrons)), ConType
        case DEFAULT
          if(.not.Current%Impart) ReIm = 'Re'
          if(     Current%Impart) ReIm = 'Im'         
          print 5,  ReIm, Current%l, Current%m,Current%Constraint,ConType
        end select
      enddo
      nullify(Current)
    endif
    !---------------------------------------------------------------------------
    !Doing the constraints on neutron & proton separately
    if (PrintWorthy(2)) then
      Current => Root   
      print 11
      print 21

      !Print the moment specific Values
      do while(associated(Current%Next))
        Current => Current%Next
        
        if(Current%Isoswitch.ne.2 .or. Current%ConstraintType.eq.0) cycle
        if(Current%Total) cycle
        
        select case (Current%ConstraintType)             
        case(1)
          ConType='Quadratic'
        case(2)
          ConType='K.Rutz   '             
        end select
        
        select case(Current%l)
        case (-2)
          !Printing the constraint on rms; some extra manipulation since
          ! the internal housekeeping is done with the true radius**2        
          print 41, sqrt(Current%TrueConstraint/(Protons+Neutrons)), ConType
        case DEFAULT
          if(.not.Current%Impart) ReIm = 'Re'
          if(     Current%Impart) ReIm = 'Im'         
          print 51,  ReIm, Current%l, Current%m,Current%Constraint,ConType
        end select
      enddo
      nullify(Current)
    endif
    
    !---------------------------------------------------------------------------
    !Doing the constraints on Q_n - Q_p
    if (PrintWorthy(3)) then
      Current => Root   
      print 111
      print 211

      !Print the moment specific Values
      do while(associated(Current%Next))
        Current => Current%Next
        
        if(Current%Isoswitch.ne.3 .or. Current%ConstraintType.eq.0) cycle
        if(Current%Total) cycle
        
        select case (Current%ConstraintType)             
        case(1)
          ConType='Quadratic'
        case(2)
          ConType='K.Rutz   '             
        end select
        
        select case(Current%l)
        case (-2)
          !Printing the constraint on rms; some extra manipulation since
          ! the internal housekeeping is done with the true radius**2        
          print 411, sqrt(Current%TrueConstraint/(Protons+Neutrons)), ConType
        case DEFAULT
          if(.not.Current%Impart) ReIm = 'Re'
          if(     Current%Impart) ReIm = 'Im'         
          print 511,  ReIm, Current%l, Current%m,Current%Constraint,ConType
        end select
      enddo
      nullify(Current)
    endif
    
    !---------------------------------------------------------------------------
    !Doing the constraints on Q_l
    if (PrintWorthy(4)) then
      Current => Root   
      print 1111
      print 2111

      !Print the moment specific Values
      do while(associated(Current%Next))
        Current => Current%Next
        
        if(Current%ConstraintType.eq.0) cycle
        if(.not. Current%Total        ) cycle
        
        select case (Current%ConstraintType)             
        case(1)
          ConType='Quadratic'
        case(2)
          ConType='K.Rutz   '             
        end select
        
        if(.not.Current%Impart) ReIm = 'Re'
        if(     Current%Impart) ReIm = 'Im'         
        ! Don't forget: constraint is the square.
        print 5111,  ReIm, Current%l, sqrt(Current%Constraint),ConType
      enddo
      nullify(Current)
    endif
    
    print *
    
  end subroutine PrintConstraints
  
  subroutine PrintMoment(ToPrint)
    !---------------------------------------------------------------------------
    ! This subroutine provides the printing of all relevant info of a Moment.
    !---------------------------------------------------------------------------
    ! Special printing rules are provided for 
    ! - l = 0 moment => Particle numbers
    ! - l =-2 moment => RMS radii
    !---------------------------------------------------------------------------
    
    type(Moment),pointer,intent(in) :: ToPrint
    character(len=2)                :: ReIm
      1 format (A2, ' Q_{', 2i2, '}', 3(1x,f15.4) ) 
      2 format ('Constrained',  2(1x,f15.4))
      3 format ('Constrained',  33x, f15.4)          
      4 format (' Particles ',  3(1x,f15.4))
      5 format (' RMS radius',  3(1x,f15.4))
      6 format ('Pulling to ',  33x, f15.4)
      7 format ('Constrained Difference: ', 3x, 2(1x,f15.4))
      8 format ('Parameter  ', 33x,  f15.4)
      9 format ('Parameter  ', 2(1x,f15.4))
      
    select case(ToPrint%l)
    
    case(0)
      ! Printing the total number of particles
      print 4, sqrt(4*pi)*ToPrint%Value, sqrt(4*pi)*sum(ToPrint%Value)
    
    case(-2)
      ! Printing RMS radii
      print 5, sqrt(ToPrint%Value(1)/Neutrons)                                &
      &      , sqrt(ToPrint%Value(2)/(Protons))                               &
      &      , sqrt(sum(ToPrint%Value)/(Neutrons+Protons))
    
      if(ToPrint%ConstraintType.ne.0) then
        if(ToPrint%Isoswitch.eq.1) then
          print 3 ,sqrt(abs(ToPrint%Constraint/(Neutrons + Protons)))
        elseif(ToPrint%Isoswitch.eq.2) then
          print 2, sqrt(abs(ToPrint%Constraint(1)/Neutrons)),                  &
          &        sqrt(abs(ToPrint%Constraint(2)/Protons))
        else
          print 7, ToPrint%Constraint(1),                                      &
          &        ToPrint%Value(1) - ToPrint%Value(2)               
        endif
      endif
    
    case DEFAULT
      !All other "normal" multipole moments
      if(.not.ToPrint%Impart) then
          ReIm = 'Re'
      else
          ReIm = 'Im'
      endif
        
      print 1, ReIm, ToPrint%l, ToPrint%m, ToPrint%Value(1), ToPrint%Value(2) &
      &        , Sum(ToPrint%Value)
      
      !if(.not.ToPrint%Total) then
        if(ToPrint%ConstraintType.ne.0) then
          select case(ToPrint%Isoswitch)
          case(1)
            print 3, ToPrint%TrueConstraint
            if(ToPrint%ConstraintType.eq.1) then
              print 6, ToPrint%Constraint
            endif
            print 8, ToPrint%Intensity(1)    
          case(2)
            print 2, ToPrint%TrueConstraint
            print 9, ToPrint%Intensity
          case(3)
            print 7, ToPrint%TrueConstraint, ToPrint%Value(1) - ToPrint%Value(2)
          end select
        endif
      !endif
    end select
    
  end subroutine PrintMoment
  
  subroutine PrintAllMoments()
      !-------------------------------------------------------------------------
      ! Subroutine that uses subroutine PrintMoment to print the information 
      ! about all the moments.
      !-------------------------------------------------------------------------
      type(Moment), pointer :: Current => null()
      integer               :: currentl
      real(KIND=dp)         :: ql(2)
    
    100 format (20('-'),' Multipole Moments ', 21('-'))
    101 format (60('-'))
      1 format (60('_'))
      2 format (16x,4x, 'Neutrons',8x, 'Protons',9x, 'Total')
      7 format ('Beta_{', 2i2 , '}', 3(1x,f15.4) )
      8 format ('Q_{',i2,'}',5x, 3(1x,f15.4))
      9 format ('Constrained', 33x, f15.4)
      Current => Root
      currentl = Current%l
      print 100
      print 2
      print 1
      
      !-------------------------------------------------------------------------
      ! Print all moments individually
      call printMoment(Current)
      do while(associated(Current%Next))
        Current => Current%Next
        !------------------------------------
        !Print a new line when getting new l.
        if(currentl .ne. Current%l) print *
        currentl = Current%l
        
        call PrintMoment(Current)    
      enddo
      print 1
      nullify(Current)
      
      !-------------------------------------------------------------------------
      ! Printing the total multipole moments Q_l
      print 2
      print 1
      do currentl=1, MaxMoment
        ql = CalculateTotalQl(currentl)
        if(all(ql.eq.0.0_dp)) cycle
        print 8, currentl,ql,sum(ql)
        
        Current => FindMoment(Currentl,0,.false.)  
        if(.not.associated(Current)) cycle       
        if(Current%ConstraintType.ne.0 .and. Current%Total) then
          print 9, sqrt(Current%Constraint)
        endif
      enddo
      print 1
      nullify(Current)
      !-------------------------------------------------------------------------   
      !Printing Beta deformation parameters
      print 2
      print 1
      Current => Root
      currentl = Current%l
      do while(associated(Current%Next))
        Current => Current%Next
        !------------------------------------
        !Print a new line when getting new l.
        if(currentl .ne. Current%l) print *
        currentl = Current%l
       
        if(Current%l.ne.1 .and. Current%l.ne.-2) then
            print 7, Current%l, Current%m,Current%Beta
        endif
      enddo
      nullify(Current)

      print 1
      
      call PrintQuadrupoleAlt
      print 101
      return
  end subroutine PrintAllMoments
  
  subroutine Calculate(ToCalculate,SaveOld)
    !---------------------------------------------------------------------------
    ! This subroutine calculates the multipole moment, both with and without 
    ! cutoff by integrating over the mesh, using the SpherHarm variable.
    ! The integration is then
    !    <Q_{lm}> = \int d^3x \rho(x,y,z)* SpherHarm(x,y,z)
    !---------------------------------------------------------------------------
    ! The integer saveold controls if previous results are saved or not to the
    ! Oldvalues array. 1 saves, 0 does not save.
    !---------------------------------------------------------------------------
    use Densities, only : Density

    type(Moment),pointer,intent(inout)  :: ToCalculate
    integer                             :: i,it
    integer                             :: SaveOld

    if(SaveOld.ne.0) then
      !Move the old values up in the list
      do i=1,6
          ToCalculate%OldValue(:,7-i+1)    = ToCalculate%OldValue(:,7-i)
          ToCalculate%OldValueCut(:,7-i+1) = ToCalculate%OldValueCut(:,7-i)
      enddo

      !Copy the present value to OldValue
      ToCalculate%OldValue(:,1)    = ToCalculate%Value
      ToCalculate%OldValueCut(:,1) = ToCalculate%ValueCut
    endif
      
    !Initialise
    ToCalculate%Value      = 0.0_dp
    ToCalculate%ValueCut   = 0.0_dp
    ToCalculate%Squared    = 0.0_dp
    ToCalculate%SquaredCut = 0.0_dp
    
    ! Calculate the new value for ordinary constraints
    do it=1,2
      ToCalculate%Value(it)    = ToCalculate%Value(it) + &
      &       sum(ToCalculate%SpherHarm(:,:,:)*Density%Rho(:,:,:,it))
      ToCalculate%ValueCut(it) = ToCalculate%ValueCut(it) + &
      &       sum(ToCalculate%SpherHarm(:,:,:)*Density%Rho(:,:,:,it)           &
      &       *CutOff(:,:,:,it))
    enddo
      
    !Do this anyway, since it does not cost any time at all
    !if(ToCalculate%ConstraintType.eq.2 .or. ToCalculate%Total) then    
      do it=1,2
        ToCalculate%Squared(it)    = ToCalculate%Squared(it)    + &
        &       sum(ToCalculate%SpherHarm(:,:,:)**2*Density%Rho(:,:,:,it))
        ToCalculate%SquaredCut(it) = ToCalculate%SquaredCut(it) + &
        &       sum(ToCalculate%SpherHarm(:,:,:)**2*Density%Rho(:,:,:,it)      &
        &       *CutOff(:,:,:,it))
      enddo
    !endif

    if(any(ToCalculate%Value.eq.ToCalculate%Value+1)) then
      call stp('Nan in the calculation of multipole moments.', 'l=',           &
      &         ToCalculate%l, 'm=', ToCalculate%m)
    endif

    ToCalculate%Value     =ToCalculate%Value*dv
    ToCalculate%ValueCut  =ToCalculate%ValueCut*dv
    ToCalculate%SquaredCut=ToCalculate%SquaredCut*dv
    ToCalculate%Squared   =ToCalculate%Squared*dv
    
    call CalcBeta(ToCalculate)
    
    return
  end subroutine Calculate
  
  subroutine CalculateAllMoments(SaveOld)
    !---------------------------------------------------------------------------
    ! Subroutine that 
    !   1) Calculates the values of all multipole moments
    !   2) Calculate the energy associated to the multipole constraints
    !   3) Calculate the quadrupole moments in different representations
    !---------------------------------------------------------------------------
    type(Moment), pointer :: Current
    integer, intent(in)   :: SaveOld
    
    nullify(Current)
    Current => Root
    
    !First, we need to calculate the cutoff function
    call CompCutoff
    call Calculate(Current,SaveOld)
    do while(associated(Current%Next))
        Current => Current%Next
        call Calculate(Current,SaveOld)    
    enddo

    call CalcConstraintEnergy()

    call CalcQuadrupoleAlt()
    
    return
  end subroutine CalculateAllMoments
  
  subroutine ReadjustAllMoments(NoRutz)
    !---------------------------------------------------------------------------
    ! Readjust all the multipole constraints. 
    ! If NoRutz is true, the non-Rutz constrained moments are updated. If NoRutz
    ! is false, the Rutz-type constraints are updated.
    !---------------------------------------------------------------------------
    type(Moment), pointer :: Current
    integer               :: NoRutz
    
    nullify(Current)
    Current => Root
    
    do while(associated(Current%Next))
        Current => Current%Next
        if(Current%ConstraintType.lt.1) cycle 
        if(NoRutz.eq.1 .and. Current%ConstraintType.eq.2) cycle
        if(NoRutz.eq.0 .and. Current%ConstraintType.ne.2) cycle
        call Readjust(Current)    
    enddo      
    nullify(Current)   

    if(RutzToQuadratic .and. NoRutz.eq.0) call ChangeRutzIntoQuadraticMoments
 
  end subroutine ReadjustAllMoments
  
  subroutine CalcConstraintEnergy()
    !---------------------------------------------------------------------------
    ! This function calculates the energy contribution associated with the 
    ! constraints on the multipole moments.
    !
    ! For quadratic constraints of the form the contribution is:
    !    U = 2*C*(<Q_{lm}> - A)* Q_{lm}
    ! For linear constraints (and Rutz-type Constraints) this becomes:
    !    U = C*Q_{lm}        
    !
    ! If the total moment Q_{l} is constrained, the contribution of every Q_{lm}
    ! (individually!) becomes:
    !    U = 2*C*(Sum_{q}[<Q_{lq}**2>] - A) *  Q_{lm}**2 
    ! or for Rutz-type Constraints:
    !    U = C * Q_{lm}**2
    !
    ! Note that there is damping of this energy in this routine with 
    ! the parameter damping:
    !       E_new = Damping* E_old + (1-Damping)*ConstraintContribution
    !---------------------------------------------------------------------------
    
    integer               :: it, Power
    real(KIND=dp)         :: Factor(2), Value(2), Desired(2)
    type(Moment), pointer :: Current
    
    Current => Root
    ConstraintEnergy = Damping*ConstraintEnergy
    
    !---------------------------------------------------------------------------
    ! Note that all of the above cases have the same structure:
    ! U = factor * ( Value - Desired) * Q_{lm}**(power)
    ! which is why we just calculate these parameters.
    
    do while(associated(Current%Next))
      Current => Current%Next
      Factor = 0.0_dp

      select case(Current%ConstraintType)
      case(0)
        !-----------------------------------------------------------------------
        ! Go to the next moment if this moment is not constrained
        cycle
      case(1)
        !-----------------------------------------------------------------------
        ! Quadratic constraints
        select case(Current%Isoswitch)
        case(1)
          ! Total value constrained
          Factor = 2.0_dp * Current%Intensity
          Desired= Current%Constraint(1)
          
          if(Current%Total) then
            Value  = sum(CalculateTotalQl(Current%l))**2
            power  = 2
          else
            Value  = sum(Current%Value)          
            power  = 1
          endif

        case(2)
          ! Proton & Neutron independently constrained
          Factor = 2.0_dp * Current%Intensity
          Desired= Current%Constraint
          
          if(Current%Total) then
            Value  = CalculateTotalQl(Current%l)
            power  = 2
          else
            Value  = Current%Value          
            power  = 1
          endif
           
        case(3)
          ! Difference of proton & neutron constrained
          Factor    = 2 * Current%Intensity
          Desired   = Current%Constraint(1)
          if(Current%Total) then
            Value  = CalculateTotalQl(Current%l)
            Value  = Value(1) - Value(2)
            power  = 2
          else
            Value  = Current%Value(1) - Current%Value(2)          
            power  = 1
          endif
        end select
      case(2)  
        !-----------------------------------------------------------------------
        ! Rutz self-correcting constraints
        Value   = 1.0_dp
        Desired = 0.0_dp          
        select case(Current%Isoswitch)
        case(1)
          Factor = Current%Intensity(1)
        case(2)
          Factor = Current%Intensity
        case(3)
          Factor(1) =   Current%Intensity(1)
          Factor(2) = - Current%Intensity(1)
        end select
        
        if(Current%Total) then
          power = 2
        else
          power = 1
        endif 
      case DEFAULT
            call stp('Not implemented constrainttype.')  
      end select
      !Put everything together
      !-------------------------------------------------------------------------   
      ! Multiplying by (1-Damp) if this is not the first time that
      ! ConstraintEnergy is calculated
      
      if(.not.all(ConstraintEnergy.eq.0.0_dp)) then
        Factor = Factor * (1.0_dp - Damping)
      endif                        
      do it=1,2
        ConstraintEnergy(:,:,:,it)=ConstraintEnergy(:,:,:,it)                  &
        & +                      Factor(it)*(Value(it) - Desired(it))*         &
        &                        Current%SpherHarm(:,:,:)**power*Cutoff(:,:,:,it)**power
      enddo   
    enddo
    nullify(Current)
  end subroutine CalcConstraintEnergy

  function CalculateTotalQl(l) result(ql)
    !---------------------------------------------------------------------------
    ! Calculates the current value of Q_{l} for a given value of l.
    ! Q_l = 4 \sqrt{pi}/\sqrt{2*l+1} *                                              
    !        \sqrt{ \Sum_{m = 0,..,l} [ < Re Q_{lm} >**2 + < Im Q_{lm}**2 > ] }
    !---------------------------------------------------------------------------
    
    integer, intent(in)   :: l
    type(Moment), pointer :: Current
    real(KIND=dp)         :: ql(2)
    
    ql = 0.0_dp
    Current => FindMoment(l,0,.false.)
    if(.not.associated(Current)) return
    
    ql = Current%Value**2
    do while(associated(Current%Next))
      Current  => Current%Next
      if(Current%l .ne.l) exit
      ql = ql + Current%Value**2
    enddo
    ql = sqrt(16 * pi/(2*l + 1) * ql)
    nullify(Current)
    
  end function CalculateTotalQl

  logical function ConverMultipole(ToCheck, Prec) result(Converged)
    !---------------------------------------------------------------------------
    ! This function checks if the associated moment has converged.
    ! If the moment is constrained to 0, this test employs absolute differences:
    !    abs(Q(i)) < Prec
    ! If not, this function uses relative differences
    !    [Q(i) - Q_{con}]/{Q_{con}} < Prec
    ! where Q(i) (i=1,8) stands for the last 8 values of the moment.
    !
    ! Only the total "constraint" is tested. 
    ! Isospin dependent test still need to implemented.
    !
    !---------------------------------------------------------------------------
    type(Moment),pointer,intent(in)     :: ToCheck
    real(KIND=dp), intent(in), optional :: Prec
    integer                             :: i
    real(KIND=dp)                       :: Difference(8), Eps
    
    if(present(Prec)) then
      ! Assigning Eps the demanded precision
      Eps=Prec
    else
      ! Assigning Eps the default precision of 10^{-4}
      Eps=Prec4
    endif

    Converged=.true.

    if(ToCheck%ConstraintType.eq.0) then
      ! Signal that the iterations have converged for unconstrained moments.
      return
    endif                

    if(sum(ToCheck%TrueConstraint).ne.0.0_dp) then
      !Constructing the relative differences
      Difference(1)=abs((sum(ToCheck%Constraint) - sum(ToCheck%Value))         &
      &                  /sum(ToCheck%Constraint))
      do i=2,8
        Difference(i)=abs((sum(ToCheck%Constraint) &
        & - sum(ToCheck%OldValue(:,i-1)))/sum(ToCheck%Constraint))
      enddo
    else
      ! Constructing the absolute differences
      Difference(1)=abs(sum(ToCheck%Constraint) - sum(ToCheck%Value))
      do i=2,8
        Difference(i)=abs(sum(ToCheck%Constraint)- sum(ToCheck%OldValue(:,i-1)))
      enddo
    endif    
    do i=1,8
      if(Difference(i).ge.Eps) then
        !Signal that the moment has not converged if the differences are too big.
        Converged=.false.
        return
      endif
    enddo

    return
  end function ConverMultipole
  
  logical function ConverMultipoleAll(MomentPrec) result(Converged)
    !---------------------------------------------------------------------------
    ! Subroutine that uses subroutine ConverMultipole to check if all the 
    ! multipole moments have converged to their constraint values.
    !---------------------------------------------------------------------------
    type(Moment), pointer     :: Current=>null()
    real(Kind=dp), intent(in) :: MomentPrec        

    Current => Root

    Converged=.true.

    do while(associated(Current%Next))
        Current => Current%Next
        Converged=Converged.and.ConverMultipole(Current, MomentPrec)    
    enddo
    return
  end function ConverMultipoleAll

  subroutine RutzCutOff()
    !---------------------------------------------------------------------------
    ! This function computes the density cut-off function for the computation of
    ! the multipole moments. Multipole constraints should be computed using a 
    ! cut-off function since the a constraint on a multipole operator of order l
    ! adds a potential term to the energy which goes as r^l. There always is a 
    ! direction on the deformation landscape which becomes unstable.     
    ! 
    ! The cutoff is defined as:
    ! 
    !                       1                                                     
    !          ------------------------------------                                       
    !          1 + exp[(DeltaR(\vec{r})-radd)/acut]  
    !
    ! Traditional values are
    !       radd = 4   fm                                            
    !       acut = 0.4 fm    
    !
    ! The definition of DeltaR is somewhat involved: it is the distance from 
    ! the meshpoint to a surface for constant density.
    !
    ! In practice the computation of DeltaR(i,j,k) is done as follows:
    !   if(Density(i,j,k).le.Treshold ) then
    !      DeltaR(i,j,k) =   min[|r(i,j,k) - r(i',j',k')|] 
    !                    with Density(i',j',k').ge.Treshold
    !   else then
    !      DeltaR(i,j,k) = - min[|r(i,j,k) - r(i',j',k')|] 
    !                    with Density(i',j',k').le.Treshold
    !    endif
    ! A good value for Treshold ought to be max(Density)/10.
    !
    ! Everything in this comment is taken from 
    !                    K. Rutz et al, Nucl. Phys. A590 (1995) 690. 
    !
    !---------------------------------------------------------------------------
      
    use Densities, only : Density
    
    real(KIND=dp)  :: Treshold(2), DeltaR(nx,ny,nz), Surface(3,7*nx*ny*nz),X,Y,Z
    real(KIND=dp)  :: InterX,InterY,InterZ, Distance
    integer        :: it,i,j,k,l, T, Sig(nx,ny,nz)
    
    if(.not.allocated(Cutoff)) allocate(Cutoff(nx,ny,nz,2))  
                                    
    do it=1,2
      Surface = 0.0_dp

      !Finding the treshold value. At the moment it is fixed to one tenth 
      !of the maximum density.
      Treshold(it) = maxval(Density%Rho(:,:,:,it))/10.0_dp
      
      !Taking a ridiculously large number as starting point
      DeltaR = 1.d12

      !T keeps count of the number of surface points the routine found
      T=0
      
      ! Sig is a sign that keeps track whether a point is on the inside or
      ! the outside of the equidensity surface.
      where(Density%Rho(:,:,:,it) .gt. Treshold(it))
              Sig = - 1
      elsewhere
              Sig =   1
      endwhere

      ! First, we construct the mesh coordinates of the equidensity surface.
      ! To do this we check for all points of the grid if the surface lies
      ! between them and their immediate neighbours.
      ! If this is the case, some linear extrapolation is done and the resulting
      ! X,Y and Z coordinates are saved to Surface. 
      do k=1,nz-1 
        Z = MeshZ(k)                 
        do j=1,ny-1
          Y = MeshY(j)
          do i=1,nx-1
           X = MeshX(i)
           
           !Initialising the interpolated values of the coordinates. 
           !This needs to be done in light of the last if in this 
           !loop-construction.                 
           InterX= 0.0_dp
           InterY= 0.0_dp
           InterZ= 0.0_dp
           
           if( ((Density%Rho(i  ,j,k,it).ge.Treshold(it)) .and.                &
           &     (Density%Rho(i+1,j,k,it).le.Treshold(it)) ) &
           &  .or. &
           &   ((Density%Rho(i  ,j,k,it).le.Treshold(it)) .and.                &
           &    (Density%Rho(i+1,j,k,it).ge.Treshold(it)))) then
            !In this case the surface lies somewhere between i and i+1                          
            InterX=X +                                                         &
            &  dx*(Treshold(it)-Density%Rho(i,j,k,it))/                        &
            & (Density%Rho(i+1,j,k,it)-Density%Rho(i,j,k,it))
            
            T = T + 1 
            
            Surface(1,T) = InterX
            Surface(2,T) = Y  
            Surface(3,T) = Z                  
           endif
           
           if( ((Density%Rho(i,j  ,k,it).ge.Treshold(it)) .and.                &
           &    (Density%Rho(i,j+1,k,it).le.Treshold(it)) ) &
           &  .or. &
           &   ((Density%Rho(i,j  ,k,it).le.Treshold(it)) .and.                &
           &    (Density%Rho(i,j+1,k,it).ge.Treshold(it)))) then
            !In this case the surface lies somewhere between j and j+1                 
            InterY = Y +                                                       &
            & dx*(Treshold(it)-Density%Rho(i,j,k,it))/                         &
            & (Density%Rho(i,j+1,k,it)-Density%Rho(i,j,k,it))
            
            T = T + 1
            
            Surface(1,T) = X
            Surface(2,T) = InterY  
            Surface(3,T) = Z  
           endif
           
           if( ((Density%Rho(i,j,k  ,it).ge.Treshold(it)) .and.                &
           &    (Density%Rho(i,j,k+1,it).le.Treshold(it)))                     &
           &  .or. &
           &   ((Density%Rho(i,j,k  ,it).le.Treshold(it)) .and.                &
           &   (Density%Rho(i,j,k+1,it).ge.Treshold(it)) ) ) then
            !In this case the surface lies somewhere between k and k+1                          
            InterZ = Z +                                                       &
            & dx*(Treshold(it)-Density%Rho(i,j,k,it))/                         &
            & (Density%Rho(i,j,k+1,it) - Density%Rho(i,j,k,it)) 
            T = T + 1
            Surface(1,T) = X
            Surface(2,T) = Y  
            Surface(3,T) = InterZ      
           endif                                                      
          enddo
        enddo
      enddo

      ! Now we can find the minimum distance from every point on the 
      ! mesh to the surface.
      do k=1,nz
        Z = MeshZ(k)
        do j=1,ny
          Y = MeshY(j)
          do i=1,nx
            X = MeshX(i)
            do l=1,T
              !Distance to the surface point
              Distance = (X - Surface(1,l))**2 + (Y -Surface(2,l))**2 +        &
              &          (Z - Surface(3,l))**2       
              DeltaR(i,j,k) = min(DeltaR(i,j,k) , Distance)  
            enddo
            DeltaR(i,j,k) = Sig(i,j,k)*sqrt(DeltaR(i,j,k) )
          enddo
        enddo
      enddo
      
      ! With this distance we can calculate the cutoff function.
      Cutoff(:,:,:,it) = OneR/(OneR + exp( (DeltaR - radd)/acut)  )                          
    enddo        
    return
  end subroutine RutzCutOff

  subroutine StandardCutoff
    !---------------------------------------------------------------------------
    ! This subroutine calculates a cutoff function in a density independent way.
    !        exp(-d)/(1+exp(-d))  if d > 0
    ! cutoff(r) = 
    !         1 /(1+exp(d))        if d < 0
    ! with d = (|r| - radd)/acut
    !
    !---------------------------------------------------------------------------

    real(kind=dp) :: X,Y,Z,d
    integer :: i,j,k

    if(.not.allocated(Cutoff)) allocate(Cutoff(nx,ny,nz,2)) 

    do k=1,nz
      Z = MeshZ(k)
      do j=1,ny
        Y = MeshY(j)
        do i=1,nx
          X = MeshX(i)
          d = (sqrt(X**2 + Y**2 + Z**2) - radd)/acut
          if( d .ge. 0.0d0) then
            Cutoff(i,j,k,:) = exp(-d)/(1 + exp(-d))
          else
            Cutoff(i,j,k,:) = 1.0d0/(1 + exp(d))
          endif
        enddo
      enddo
    enddo

  end subroutine StandardCutoff
  
  pure function LegacyQuad (iq1,iq2) result(Constraints)
    !---------------------------------------------------------------------------
    ! Function that converts the legacy iq1, iq2 input from EV8 to new style
    ! input on Q20 & Q22.
    ! See conversion coefficients in W.Ryssens et al., CPC (to be published)
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: iq1,iq2
    real(KIND=dp)             :: Constraints(2)
    real(KIND=dp)             :: Q20, Q22
    Constraints=0.0_dp
  
    ! Q20 = 1/2 *sqrt(5/(16*pi))* ( 2 * q1 + q2)
    Q20 = 0.5_dp *sqrt(5.0_dp/(16.0_dp * pi)) * ( 2 * iq1 + iq2 ) 
  
    ! Q22 = 3/2*q2 * sqrt(5/(4*pi)) * 1/(2 * sqrt(6)) 
    Q22 = 1.5_dp * iq2 * sqrt(5.0_dp/(16.0_dp * pi)) /(2 * sqrt(6.0_dp))
    
    Constraints(1) = Q20
    Constraints(2) = Q22
    
  end function LegacyQuad
  
  subroutine Readjust(ToReadjust)
    !---------------------------------------------------------------------------
    ! Subroutine that readjusts the constraint of a certain multipole moment.
    !---------------------------------------------------------------------------    
    type(Moment), pointer             :: ToReadjust
    real(KIND=dp)                     :: O2(2)
    integer                           :: it

    select case(ToReadjust%ConstraintType)
    
    case (0)
      call stp('Moment is not constrained, but still gets readjusted!')
    
    case (1) 
      !-------------------------------------------------------------------------
      !Quadratically constrained.     
      select case(ToReadjust%Isoswitch)
      case(1)
        ToReadjust%Constraint(1) = ToReadjust%Constraint(1) &
        & -ReadjustSlowDown*(sum(ToReadjust%Value)-ToReadjust%TrueConstraint(1))
      case(2)  
        ToReadjust%Constraint = ToReadjust%Constraint &
        & -ReadjustSlowDown * (   ToReadjust%Value - ToReadjust%TrueConstraint )
      case(3)
        ToReadjust%Constraint(1) = ToReadjust%Constraint(1) &
        & -ReadjustSlowDown * ( ToReadjust%Value(1) - ToReadjust%Value(2)      &
        &                   - ToReadjust%TrueConstraint(1))
      end select    
    
    case (2)
     !--------------------------------------------------------------------------
     ! Rutz-type constraints. See the article by Cusson or the PhD thesis
     ! of K. Rutz. 
     select case(ToReadjust%Isoswitch)
     
     case(1)    
      !Calculate < O^2 >
      O2 = sum(ToReadjust%Squared)

      if(ToReadjust%Total) O2 = O2 * 16 * pi/(2*ToReadjust%l + 1)
     
      ToReadjust%Intensity = ToReadjust%Intensity +                            &
      & epsilon *  (sum(ToReadjust%Value) - sum(ToReadjust%OldValue(:,1))) /   &
      & (O2(1) + d0)
    
     case(2)
      !Calculate < O^2 >
      do it=1,2
        O2(it) = ToReadjust%Squared(it)
        
        if(ToReadjust%Total) O2(it) = O2(it) * 16 * pi/(2*ToReadjust%l + 1)

        ToReadjust%Intensity(it) = ToReadjust%Intensity(it)+                   &
        & epsilon *  ( ToReadjust%Value(it) - ToReadjust%OldValue(it,1))/      &
        & (O2(it) + d0)     
      enddo
      
     case(3) 
      ! Calculate < O^2 >
      O2 = sum(ToReadjust%Squared)

      if(ToReadjust%Total) O2(it) = O2(it) * 16 * pi/(2*ToReadjust%l + 1)
      
      ToReadjust%Intensity = ToReadjust%Intensity +                            &
      & epsilon *  (ToReadjust%Value(1) - ToReadjust%Value(2)                  &
      &             - ToReadjust%OldValue(1,1) + ToReadjust%OldValue(2,1))   / &
      & (O2(1) + d0)
      
     case DEFAULT
      call stp('Non-supported value of constrainttype.')
      end select
    end select
    return
  end subroutine Readjust
  
  subroutine ConstrainNonPhysicalMoments
    !---------------------------------------------------------------------------
    ! A subroutine that automatically assigns the correct constraining 
    ! parameters to the multipole moments that do not represent physical degrees
    ! of freedom. Those multipole moments are:
    !
    ! Q10     => Z-coordinate of the center of mass
    ! Re(Q11) => X-coordinate of the center of mass
    ! Im(Q11) => Y-coordinate of the center of masss
    !
    ! Re(Q21), Im(Q21), Im(Q22) => orientation of the nucleus in the box
    !
    ! All these moments just take the default type: Rutz-type.
    !---------------------------------------------------------------------------
      type(Moment), pointer :: Current => null()

      !Find the center of mass coordinates
      Current => FindMoment(1,0,.false.)
      if(associated(Current)) then
          Current%ConstraintType=2          
          allocate(Current%Constraint(1), Current%TrueConstraint(1))
          allocate(Current%Intensity(1)); Current%Intensity=0.0_dp
          Current%Constraint=0.0_dp ; Current%TrueConstraint=0.0_dp
          Current%Isoswitch=1       ; Current%Total=.false.
      endif
      
      Current => FindMoment(1,1,.false.)
      if(associated(Current)) then
          Current%ConstraintType=2
          allocate(Current%Constraint(1), Current%TrueConstraint(1))
          allocate(Current%Intensity(1)); Current%Intensity=0.0_dp
          Current%Constraint=0.0_dp ; Current%TrueConstraint=0.0_dp
          Current%Isoswitch=1
      endif
      
      Current => FindMoment(1,1,.true.)
      if(associated(Current)) then
          Current%ConstraintType=2
          allocate(Current%Constraint(1), Current%TrueConstraint(1))
          allocate(Current%Intensity(1)); Current%Intensity=0.0_dp
          Current%Constraint=0.0_dp  ; Current%TrueConstraint=0.0_dp
          Current%Isoswitch=1
      endif
      !Find the rotational degrees of freedom
      Current => FindMoment(2,1,.false.)
      if(associated(Current)) then
         Current%ConstraintType=2
         allocate(Current%Constraint(1), Current%TrueConstraint(1))
         allocate(Current%Intensity(1)); Current%Intensity=0.0_dp
         Current%Constraint=0.0_dp ;  Current%TrueConstraint=0.0_dp
         Current%Isoswitch=1
      endif        

      Current => FindMoment(2,1,.true.)
      if(associated(Current)) then
          Current%ConstraintType=2
          allocate(Current%Constraint(1), Current%TrueConstraint(1))
          allocate(Current%Intensity(1)); Current%Intensity=0.0_dp
          Current%Constraint=0.0_dp ; Current%TrueConstraint=0.0_dp
          Current%Isoswitch=1
      endif        

      Current => FindMoment(2,2,.true.)
      if(associated(Current)) then
          Current%ConstraintType=2
          allocate(Current%Constraint(1), Current%TrueConstraint(1))
          allocate(Current%Intensity(1)); Current%Intensity=0.0_dp
          Current%Constraint=0.0_dp ; Current%TrueConstraint=0.0_dp
          Current%Isoswitch=1
      endif       
  end subroutine ConstrainNonPhysicalMoments
  
  subroutine CalcQuadrupoleAlt()
    !---------------------------------------------------------------------------
    ! Calculates the Quadrupole moments in their various representations.
    ! See the formulas in W. Ryssens et al. (to be published)
    !---------------------------------------------------------------------------
    ! The off-diagonal components of the quadrupole tensor are kept, but not 
    ! used at the moment....
    !---------------------------------------------------------------------------
    type(Moment),pointer :: Q20, RQ22 , IQ22, RQ21, IQ21
    integer              :: i
    
    ! Finding the moments. Note that moments that are not allowed will not be
    ! found and the pointers will be null.
    Q20 =>FindMoment(2,0,.false.     )
    RQ21=>FindMoment(2,1,.false., Q20)
    IQ21=>FindMoment(2,1,.true. , Q20)
    RQ22=>FindMoment(2,2,.false., Q20)
    IQ22=>FindMoment(2,2,.true. , Q20)  
    
    QCartesian = 0.0_dp
    
    !Cartesian Moments
    do i=1,2
      ! a) Diagonal elements
      QCartesian(1,1,i) = - (sqrt(5.0_dp/( 4.0_dp * pi)))**(-1)  * &
      &         ( Q20%Value(i) - sqrt(6.0_dp) * RQ22%Value(i) )
      QCartesian(2,2,i) = - (sqrt(5.0_dp/( 4.0_dp * pi)))**(-1)  * &
      &         ( Q20%Value(i) + sqrt(6.0_dp) * RQ22%Value(i) )
      QCartesian(3,3,i) =   (sqrt(5.0_dp/(16.0_dp * pi)))**(-1)  * &
      &           Q20%Value(i)
    
      ! b) Off-diagonal elements
      if(associated(RQ21)) then
        ! RQ21 = sqrt(15/(pi 4)) zx
        QCartesian(1,3,i) = RQ21%Value(i) * sqrt(4.0_dp * pi/ (15.0_dp)) 
        QCartesian(3,1,i) = QCartesian(1,3,i)
      endif
      if(associated(IQ21)) then
        !IQ21 = sqrt(15/(4 pi)) yz
        QCartesian(2,3,i) = IQ21%Value(i) * sqrt(4.0_dp * pi/ (15.0_dp))
        QCartesian(3,2,i) = QCartesian(2,3,i)
      endif
      if(associated(IQ22)) then
        !IQ22 = sqrt(15/(4*pi)) xy
        QCartesian(1,2,i) = IQ22%Value(i) * sqrt(4.0_dp * pi/ (15.0_dp))
        QCartesian(2,1,i) = QCartesian(1,2,i)
      endif
    enddo
    
    !Total moments
    QCartesian(:,:,3)=sum(QCartesian(:,:,1:2),3)
    
    ! Q & Gamma representation
    do i=1,3
      Q(i) = sqrt(2.0_dp/3.0_dp *                                              &
      &   ( QCartesian(1,1,i)**2 + QCartesian(2,2,i)**2 + QCartesian(3,3,i)**2))
      if(abs(QCartesian(3,3,i)).gt.1d-08 ) then
        G(i) = atan2((QCartesian(1,1,i) - QCartesian(2,2,i)),                  &
        &                                     (sqrt(3.0_dp)* QCartesian(3,3,i)))   
      else
        G(i) = 0.0_dp
      endif
    enddo  
    ! iq1 & iq2
    do i=1,3
      iq1(i) = Q(i) * cos(G(i)) - 1/sqrt(3.0_dp) * Q(i) * sin(G(i))
      iq2(i) = 2/sqrt(3.0_dp)*Q(i)*sin(G(i))
    enddo
  end subroutine CalcQuadrupoleAlt
    
  subroutine PrintQuadrupoleAlt()
    !---------------------------------------------------------------------------
    ! Prints out the Quadrupole moments in their various representations.
    ! This should be the numbers to compare with EV8/CR8 and EV4.
    !---------------------------------------------------------------------------
    
   10 format(' Cart. Moments', ' X' , 14x, 'Y', 14x, 'Z')  
    1 format( 'N', 7x ,3f15.7)
    2 format( 'P', 7x ,3f15.7)
    3 format( 'T', 7x ,3f15.7)
   11 format(' (Q, gamma)   ', ' Q', 14x, 'Gamma') 
   12 format(' (iq1, iq2)   ', ' iq1', 12x, 'iq2') 

    !real(kind=DP)        :: QAlt(3), GAlt(3)

    print 10
    !Printing only the diagonal components
    print 1, QCartesian(1,1,1), QCartesian(2,2,1), QCartesian(3,3,1)
    print 2, QCartesian(1,1,2), QCartesian(2,2,2), QCartesian(3,3,2)
    print 3, QCartesian(1,1,3), QCartesian(2,2,3), QCartesian(3,3,3)

  !-----------------------------------------------------------------------------
  !    !Double checking the results (they are not printed)
  !    do i=1,2
  !      QAlt(i) = sqrt(16.0_dp/5.0_dp*pi*(Q20%Value(i)**2 + Q22%Value(i)**2))
  !      GAlt(i) = atan2(sqrt(2.0_dp) * Q22%Value(i), Q20%Value(i))
  !    enddo
  !    QAlt(3) = QAlt(1) + QAlt(2)
  !    GAlt(3) = atan2(sqrt(2.0_dp) * sum(Q22%Value(:)), sum(Q20%Value(:)))
  !-----------------------------------------------------------------------------

    print 11
    !Note that Gamma is represented in radians, but is printed in degrees!
    print 1, Q(1), G(1)* 360.0_dp/(2*pi)!, QAlt(1), GAlt(1)* 360.0_dp/(2*pi)
    print 2, Q(2), G(2)* 360.0_dp/(2*pi)!, QAlt(2), GAlt(2)* 360.0_dp/(2*pi)
    print 3, Q(3), G(3)* 360.0_dp/(2*pi)!, QAlt(3), GAlt(3)* 360.0_dp/(2*pi)
    
    print 12
    print 1, iq1(1), iq2(1)
    print 2, iq1(2), iq2(2)
    print 3, iq1(3), iq2(3)
    
  end subroutine PrintQuadrupoleAlt
  
  subroutine WriteMoment(Mom,Ochan)
    !---------------------------------------------------------------------------
    ! Type-bound subroutine that writes all relevant values of multipole moment
    ! to file.
    !---------------------------------------------------------------------------
    class(Moment), intent(in) :: Mom
    integer, intent(in)       :: Ochan
    integer                   :: io
    
    io = 0
    if(Mom%ConstraintType.ne.0) then 
        write(OChan)               Mom%ConstraintType, Mom%Isoswitch           &
        &                        , Mom%l,Mom%m,Mom%Impart
        select case(Mom%Isoswitch) 
        case(1)
          ! To have a consistent length of the record
          write(OChan)              Mom%Intensity,     Mom%Intensity,          &
          &                         Mom%Constraint,    Mom%Constraint,         &
          &                         Mom%TrueConstraint,Mom%TrueConstraint,     &
          &                         Mom%Value         ,Mom%OldValue
        case(2)
          write(OChan)              Mom%Intensity, Mom%Constraint,             &
          &                         Mom%TrueConstraint, Mom%Value,             &
          &                         Mom%OldValue
        case(3)
          ! To have a consistent length of the record
          write(OChan)              Mom%Intensity,     Mom%Intensity,          &
          &                         Mom%Constraint,    Mom%Constraint,         &
          &                         Mom%TrueConstraint,Mom%TrueConstraint,     &
          &                         Mom%Value         ,Mom%OldValue  
        end select
    endif
    if(io.ne.0) call stp('Error while writing multipole moment to file!')
  end subroutine WriteMoment
  
  subroutine ReadMoment(IChan, ioerror)
    !---------------------------------------------------------------------------
    ! Non type-bound procedure that parses the next line for input on a moment.
    ! This subroutine finds the appropriate moment (if possible) and assigns it
    ! the correct values.
    !---------------------------------------------------------------------------
    integer, intent(in)    :: IChan
    integer, intent(inout) :: ioerror
    
    type (Moment), pointer :: Mom
    
    integer                :: l,m, ConStraintType, isoswitch
    real(KIND=dp)          :: Value(2), Intensity(2), Constraint(2),OldValue(2,7)
    real(KIND=dp)          :: TrueConstraint(2)
    logical                :: impart
    
    read(IChan, iostat=ioerror) ConstraintType, Isoswitch,l,m,impart
    if(ioerror.ne.0) return
    read(IChan, iostat=ioerror) Intensity,Constraint,TrueConstraint,           &
    &                           Value, OldValue
    if(ioerror.ne.0) return    
    !---------------------------------------------------------------------------
    ! Finding moment.
    Mom => FindMoment(l,m,Impart)
    if(.not.associated(Mom)) call stp("Tried to read moment from file that"    &
    &                             //  " isn't there.")
    !---------------------------------------------------------------------------
    ! Assigning correct values
    Mom%ConstraintType = ConstraintType
    Mom%Isoswitch      = Isoswitch
    
    Mom%Value    = Value
    Mom%OldValue = OldValue
    if(allocated(Mom%Constraint)) then
      deallocate(Mom%Constraint, Mom%TrueConstraint, Mom%Intensity)
    endif
    select case(Mom%isoswitch)
    case(1)
      allocate(Mom%Constraint(1), Mom%TrueConstraint(1),Mom%Intensity(1))
      Mom%Constraint     = Constraint(1)
      Mom%Intensity      = Intensity(1)  
      Mom%TrueConstraint = TrueConstraint(1)
    case(2)
      allocate(Mom%Constraint(2), Mom%TrueConstraint(2))
      Mom%Constraint     = Constraint
      Mom%Intensity      = Intensity
      Mom%TrueConstraint = TrueConstraint
    case(3)
      allocate(Mom%Constraint(1), Mom%TrueConstraint(1))
      Mom%Constraint     = Constraint(1)
      Mom%Intensity      = Intensity(1)
      Mom%TrueConstraint = TrueConstraint(1)    
    end select
  end subroutine ReadMoment
  
  function CheckForRutzMoments() result(Check)
  !-----------------------------------------------------------------------------
  ! Function that checks if any multipole moments are constrained according
  ! to the Rutz prescription. This implies a slightly different iteration 
  ! algorithm.
  !-----------------------------------------------------------------------------
    logical             :: Check
    type(Moment),pointer:: Current
    
    Check = .false.
    Current => Root
    
    do while(associated(Current%Next))
      Current => Current%Next
      if(Current%ConstraintType.eq.2) then
        Check = .true.
        exit
      endif
    enddo
  end function CheckForRutzMoments
  
  subroutine CalcBeta(Mom) 
  !-----------------------------------------------------------------------------
  ! Function that calculates the beta_lm deformation parameters associated with
  ! a multipole moment Q_lm.
  !-----------------------------------------------------------------------------
    type(Moment), intent(inout) :: Mom
    real(KIND=dp)               :: R, factor
  
    R = 1.2_dp  * (neutrons + protons)**(1.0_dp/3.0_dp)
    factor = 4.0_dp * pi /(3.0_dp * (neutrons+ protons) * R**(Mom%l))
  
    Mom%Beta(1:2) = factor*Mom%Value
    Mom%Beta(3)   = factor*sum(Mom%Value) 
  end subroutine CalcBeta
  
  subroutine ReadMomentData
    !---------------------------------------------------------------------------
    ! A separate subroutine that governs the input of data related to the 
    ! multipole moments.
    !
    ! Note that also legacy input is admitted with iq1 & iq2.
    !---------------------------------------------------------------------------           

    integer             :: iostat
    integer             :: l,m, ConstraintType, isoswitch
    logical             :: Impart
    real(KIND=dp)       :: Intensity, ConstraintNeutrons, ConstraintProtons
    real(KIND=dp)       :: Constraint, iq1, iq2
    real(KIND=dp)       :: iq1neutron, iq2neutron, iq1proton, iq2proton
    logical             :: MoreConstraints=.false., Total
    type(Moment),pointer:: Current, New
    
    real(KIND=dp), allocatable:: LegacyCon(:)
    
    NameList /MomentParam/ MaxMoment, radd, acut, MoreConstraints,Damping,     &
    &                      ReadjustSlowDown, CutoffType, ContinueMoment        &
    &                     ,c0,d0, epsilon
    
    NameList /MomentConstraint/ l,m,Impart, Isoswitch, Intensity,              &
    &                           ConstraintNeutrons,ConstraintProtons,          &
    &                           Constraint, MoreConstraints, ConstraintType,   &
    &                           iq1, iq2, Total

    nullify(Current)
    
    iostat = 0
    l = 1
  
    !Reading the parameters
    read(unit=*,NML=MomentParam )
    
    !Initialising the Linked List
    call IniMoments()

    allocate(ConstraintEnergy(nx,ny,nz,2))
    ConstraintEnergy = 0.0_dp

    !Choosing cutoff
    nullify(CompCutoff)
    select case (CutoffType)   
      case(0)
        !Rutz-Type cutoff
        CompCutoff => RutzCutoff
      case(1)
        !Spherical cutoff
        CompCutoff => StandardCutoff
      case default
        call stp('Invalid choice of cutofftype. It can either be 0 or 1.')
    end select

    !Constraining the non-physical degrees of freedom, if present
    call ConstrainNonPhysicalMoments()

    !Do-loop exits when MoreConstraints is indicated to be false.
    do while(MoreConstraints)
        !Initialisation
        l=0; m=0; Impart=.false.; Isoswitch = 1; Intensity=0.0_dp
        ConstraintNeutrons=0.0_dp; ConstraintProtons=0.0_dp;Constraint=0.0_dp
        MoreConstraints=.false.; ConstraintType=2; iq1=0.0_dp; iq2=0.0_dp
        iq1neutron=0.0_dp; iq2neutron=0.0_dp; iq1proton=0.0_dp
        iq2proton=0.0_dp; Total = .false.

        read(unit=*, NML=MomentConstraint)!, IOSTAT=iostat)
        if(iostat .ne. 0) then
          call stp("Input error for moment constraints."                       &
          &   //" Maybe you wanted constraints and didn't specify any? ")
        endif
        
        if(Total .and. Isoswitch.ne.1) then
          call stp('MOCCa is unsure what you want do to: constrain '    //&
          &        ' Q_l, but differently for protons in some way?')
        endif
        if(Total .and. l.eq.-2) then
          call stp('Constraining rms radii in the sense of Q_l is undefined.')
        endif
        
        !Finding the Correct Moment to constrain
        Current => FindMoment(l,m,Impart)
        if(.not.associated(Current)) then
          call stp("MOCCa can't find this moment!", 'l= ', l, ' m= ', m)
        endif
        
        if(            iq1.ne.0.0_dp .or. iq2.ne.0.0_dp                        &
        &   .or. iq1proton.ne.0.0_dp .or. iq2proton.ne.0.0_dp                  &
        &  .or. iq1neutron.ne.0.0_dp .or. iq2neutron.ne.0.0_dp) then
          if((l.ne.2) .or. (m.ne.2 .and. m .ne. 0)) then
            !Only allow iq1 & iq2 for Q20 and Q22
            
            call stp('You specified legacy iq1,iq2 input for moments that '    &
            &      //'are not Q20 or Q22.')            
          else
            !Converting input
            if(iq1proton.ne.0.0_dp .or. iq2proton.ne.0.0_dp .or.               &
               iq1neutron.ne.0.0_dp .or. iq1proton.ne.0.0_dp) then
               allocate(LegacyCon(4))
               LegacyCon(1:2) = LegacyQuad(iq1neutron, iq2neutron)
               LegacyCon(1:2) = LegacyQuad(iq1proton, iq2proton)
            else
               allocate(LegacyCon(2))
               LegacyCon = LegacyQuad(iq1,iq2)
            endif
          endif
        endif
   
        !Setting the parameters of the moment
        Current%ConstraintType = ConstraintType        
        !Reading the values for the constraints
        if(ConstraintType.ne.0) then       
            Current%Isoswitch = Isoswitch
            Current%Total     = Total
            if(Isoswitch .eq. 1 .or. Isoswitch.eq. 3) then                              
                ! To be able to override data from input
                if(allocated(Current%Constraint)) then
                    deallocate(Current%Constraint, Current%TrueConstraint,     &
                    &          Current%Intensity)
                endif
                allocate(Current%Constraint(1))
                allocate(Current%TrueConstraint(1))    
                allocate(Current%Intensity(1))                             
                !---------------------------------------------------------------
                ! Deal with legacy input
                if( allocated(LegacyCon)                                       &
                &  .and. l.eq. 2 .and. m.eq.2 .and. .not. Impart) then
                  Current%Constraint = LegacyCon(2)
                  Current%TrueConstraint = LegacyCon(2)
                elseif(allocated(LegacyCon)                                    &
                &  .and. l.eq. 2 .and. m.eq.0 .and. .not. Impart) then
                  Current%Constraint = LegacyCon(1)
                  Current%TrueConstraint = LegacyCon(1)               
                else     
                  if(Current%l .eq. -2) then
                    ! Practically we constrain the radius squared instead of the
                    ! rms radius.
                    Current%Constraint = Constraint**2 * (Protons + Neutrons)
                  else
                    if(Current%Total) then
                      !In practice, constrain the square
                      Current%Constraint = Constraint**2            
                    else
                      Current%Constraint = Constraint     
                    endif
                  endif
                  !Copying the true constraint
                  Current%TrueConstraint = Current%Constraint
                endif
                !---------------------------------------------------------------
            else
                if(allocated(Current%Constraint)) then
                    deallocate(Current%Constraint, Current%TrueConstraint,     &
                    &          Current%Intensity)
                endif
                allocate(Current%Constraint(2))
                allocate(Current%TrueConstraint(2))
                allocate(Current%Intensity(2))
                if(allocated(LegacyCon)                                        &
                &  .and. l.eq. 2 .and. m.eq.0 .and. .not. Impart) then
                  Current%Constraint(1) = LegacyCon(1)
                  Current%Constraint(2) = LegacyCon(2)
                elseif(allocated(LegacyCon)                                    &
                &  .and. l.eq. 2 .and. m.eq.2 .and. .not. Impart) then
                  Current%Constraint(1) = LegacyCon(3)
                  Current%Constraint(2) = LegacyCon(4)
                else
                  
                  if(Current%l .eq. -2) then
                    ! Practically we constrain the radius squared instead of the
                    ! rms radius.
                    Current%Constraint(1) = ConstraintNeutrons**2*Neutrons
                    Current%Constraint(2) = ConstraintProtons**2*Protons
                  else
                    if(Current%Total) then
                      ! In practice we constrain Ql**2
                      Current%Constraint(1) = ConstraintNeutrons**2
                      Current%Constraint(2) = ConstraintProtons**2
                    else
                      Current%Constraint(1) = ConstraintNeutrons
                      Current%Constraint(2) = ConstraintProtons           
                    endif
                  endif
                  if(Current%Total) then
                    Current%TrueConstraint(1) = Current%Constraint(1)
                    Current%TrueConstraint(2) = Current%Constraint(2)   
                  else
                    Current%TrueConstraint(1) = Current%Constraint(1)**2
                    Current%TrueConstraint(2) = Current%Constraint(2)**2
                  endif
                endif
            endif
            !-------------------------------------------------------------------
            ! If Total is demanded, make sure that all necessary moments have 
            ! the correct data.
            if(Current%Total) then
              !Find the starting moment
              New => Findmoment(Current%l, 0,.false.)
              do while(New%l .eq. Current%l)
                New%ConstraintType = Current%ConstraintType
                New%Constraint     = Current%Constraint
                New%TrueConstraint = Current%TrueConstraint
                New%Total          = Total
                New%Intensity      = Current%Intensity
                
                if(.not.associated(New%Next)) exit
                New => New%Next
              enddo              
            endif            
        endif       
        Current%Intensity = Intensity
     enddo
     nullify(Current)
     return
  end subroutine ReadMomentData

  subroutine ChangeRutzIntoQuadraticMoments
    !-----------------------------------------------------------------------------------
    ! In order to gain speed we try to replace Rutzconstraints by quadratic constraints
    ! once close enough to the solution.
    !
    ! The scheme is simple, close to convergence the Rutz Constraint will have the following
    ! form:
    ! 
    !  h - \lambda O
    !
    ! And we then seek a quadratic constraint that gives the same contribution to the
    ! single-particle Hamiltonian.
    ! Remember that if E - C(<O> - O_0)^2 then the contribution to the single-particle
    ! Hamiltonian is:
    !
    ! h - c( <O> - O_0 ) * O
    !
    ! Thus, we replace a Rutz-Constraint with a quadratic constraint at iteration i
    ! with
    !
    ! c = \lambda/( <O> - O_0 )
    !------------------------------------------------------------------------------------
    integer :: it
    real(KIND=dp) :: Treshold = 0.1, Deviation(2), O2(2)
    type(MOment), pointer :: Current

    Current => Root

    do while(associated(Current%Next))
      Current => Current%Next
      !Don't pay attention to non-Rutz moments
      if(Current%ConstraintType .ne. 2) cycle
      Deviation = abs(sum(Current%Value)  - sum(Current%Constraint))/(sum(abs(Current%Constraint)))
      !print *, Deviation, Current%Value, Current%Constraint
      !Check if is close enough
      
      if(all(Deviation.lt. Treshold)) then
        print *, 'Changed ConstraintType ', Current%l, Current%m, Deviation
        Current%ConstraintType = 1
        Current%Intensity      = abs(Current%Intensity/(sum(Current%Value) - sum(Current%Constraint)))
        print *, 'New Intensity', Current%Intensity
      endif
    enddo
  end subroutine ChangeRutzIntoQuadraticMoments

end module Moments
