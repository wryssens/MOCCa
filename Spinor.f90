module spinors

  use CompilationInfo
  use GenInfo
 
  type Spinor
  
    real(KIND=dp), allocatable :: Grid(:,:,:,:,:)
    
    contains
            
    procedure, public, pass :: SetComponent
    procedure, public, pass :: DeriveSpinor
    procedure, public, pass :: LapSpinor
    procedure, public, pass :: ActionOfTimeSimplex
    procedure, public, pass :: TimeReverse
    procedure, public, pass :: ActionOfParity  
    procedure, public, pass :: ActionOfSignature
  
  end type Spinor
  
  interface operator (*)        
    !Overloading "*" to be used on Spinors and various other object
    module procedure MultiplyScalar, MultiplyFunction,                         &
    &                MultiplySpinor, MultiplyComplex
    !---------------------------------------------------------------------------
    ! IMPORTANT:
    !   a)  This overloading only works for the following forms:
    !       1) Scalar * Spinor
    !       2) Function * Spinor
    !       This doesn't work for
    !       1) Spinor*Scalar
    !       2) Spinor*Function
    !
    !       Unless one decides to fix this by adding more subroutines
    !
    !   b) Notice that Spinor*Spinor actually calculates
    !       \psi^{\dagger} * \psi
    !      as that is the only meaningful way to define such a multiplication.
    !---------------------------------------------------------------------------
  end interface
    
  interface operator (+)
      !Overloading "+" to be used to add spinors.
      module procedure Add
  end interface
    
  interface operator (-)
    module procedure Subtract
    module procedure Negative
  end interface

  !Constructor interface
  interface spinor    
    module procedure newspinor 
  end interface

contains
  
  pure function NewSpinor( insizex, insizey, insizez) result (Psi)
    !---------------------------------------------------------------------------
    ! This function serves as constructor for the Spinor Class, taking 
    ! proper care of allocation and initialisation to zero.
    !---------------------------------------------------------------------------
    type(Spinor) :: Psi
    integer      :: sizex,sizey,sizez
    integer, intent(in), optional :: insizez, insizey, insizex

    sizex = nx
    sizey = ny
    sizez = nz
    
    if(present(insizex)) then
      sizex = insizex
      sizey = insizey
      sizez = insizez
    endif

    allocate(Psi%Grid(sizex,sizey,sizez,4,1))
    Psi%Grid=0.0_dp
  end function NewSpinor
  
  subroutine SetComponent(Psi, Component, Value, IsoComponent)
    !---------------------------------------------------------------------------
    ! This function changes the value of one component of the spinor Psi.
    !---------------------------------------------------------------------------
    class(Spinor)                :: Psi
    integer, intent(in)          :: Component
    integer, intent(in),optional :: IsoComponent
    real(KIND=dp), intent(in)    :: Value(nx,ny,nz)

    if(present(IsoComponent)) then
            Psi%Grid(:,:,:,Component,IsoComponent) = Value
    else
            Psi%Grid(:,:,:,Component,1) = Value
    endif
    
    return                
  end subroutine SetComponent
                   
  pure function Pauli(Psi,Direction) result(SigmaPsi)
    !---------------------------------------------------------------------------
    ! This function returns the action of the pauli matrices on Psi. 
    ! Direction: (1,2,3) = (x,y,z)
    !---------------------------------------------------------------------------
    class(Spinor),intent(in) :: Psi
    integer, intent(in)      :: Direction
    real(Kind=dp)            :: PauliSpinor(nx,ny,nz,4,1)
    type(Spinor)             :: SigmaPsi
                    
    if(Direction.eq.1) then
        !\sigma_x = ( 0  1 )
        !           ( 1  0 )      
        PauliSpinor(:,:,:,1,:) = Psi%Grid(:,:,:,3,:)
        PauliSpinor(:,:,:,2,:) = Psi%Grid(:,:,:,4,:)
        PauliSpinor(:,:,:,3,:) = Psi%Grid(:,:,:,1,:)
        PauliSpinor(:,:,:,4,:) = Psi%Grid(:,:,:,2,:)         
    elseif(Direction.eq.2) then                
        !\sigma_y = ( 0 -i )
        !           ( i  0 )          
        PauliSpinor(:,:,:,1,:) =   Psi%Grid(:,:,:,4,:)
        PauliSpinor(:,:,:,2,:) = - Psi%Grid(:,:,:,3,:)
        PauliSpinor(:,:,:,3,:) = - Psi%Grid(:,:,:,2,:)
        PauliSpinor(:,:,:,4,:) =   Psi%Grid(:,:,:,1,:)       
    elseif(Direction.eq.3) then
        !\sigma_z = ( 1  0 )
        !           ( 0 -1 )       
        PauliSpinor(:,:,:,1,:) =   Psi%Grid(:,:,:,1,:)
        PauliSpinor(:,:,:,2,:) =   Psi%Grid(:,:,:,2,:)
        PauliSpinor(:,:,:,3,:) = - Psi%Grid(:,:,:,3,:)
        PauliSpinor(:,:,:,4,:) = - Psi%Grid(:,:,:,4,:)
    endif
    SigmaPsi=NewSpinor()
    SigmaPsi%Grid = PauliSpinor
    return
  end function Pauli
  pure function MultiplyFunction(F,Psi) result(FPsi)
    !---------------------------------------------------------------------------
    ! This function multiplies the spinor Psi with a scalar function F.
    ! It is used for overloading the operator *.
    !---------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: FPsi
    real(KIND=dp),intent(in)  :: F(nx,ny,nz)
    integer                   :: i
    
    FPsi=NewSPinor()
    do i=1,4
      FPsi%Grid(:,:,:,i,1) = F*Psi%Grid(:,:,:,i,1)
    enddo
    return 
  end function MultiplyFunction
  pure function MultiplyScalar(S,Psi) result(SPsi)
    !---------------------------------------------------------------------------
    ! This function multiplies the spinor Psi with a scalar S.
    ! It is used for overloading the operator *.
    !---------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: SPsi
    real(KIND=dp),intent(in)  :: S
    
    SPsi=NewSPinor()
    SPsi%Grid = S*Psi%Grid
    
    return
  end function MultiplyScalar
  
  pure function MultiplyComplex(S,Psi) result(SPsi)
    !---------------------------------------------------------------------------
    ! This function multiplies the spinor Psi with a complex scalar S.
    ! It is used for overloading the operator *.
    !---------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: SPsi
    complex(KIND=dp),intent(in)  :: S

    SPsi = MultiplyI(Psi)
    SPsi = DReal(S) * Psi + Imag(S) * SPsi
    
    return
  end function MultiplyComplex
  
  pure function MultiplySpinor(Psi,Phi) result(PsiPhi)
    !---------------------------------------------------------------------------
    ! This function returns the scalar complex value that is the product of 
    ! two Spinors:
    !       PsiPhi = Psi^{\dagger}\Phi
    ! It is used for overloading "*".
    !---------------------------------------------------------------------------
    class(Spinor),intent(in) :: Psi, Phi
    real(KIND=dp)            :: PsiPhi(nx,ny,nz,2)
                            
    PsiPhi = 0.0_dp
    
    !Real part
    PsiPhi(:,:,:,1) =  sum( Psi%Grid(:,:,:,:,1)* Phi%Grid(:,:,:,:,1),4)
    !Imaginary Part
    PsiPhi(:,:,:,2) =  Psi%Grid(:,:,:,1,1)*Phi%Grid(:,:,:,2,1)                 &
                  & -  Psi%Grid(:,:,:,2,1)*Phi%Grid(:,:,:,1,1)                 &
                  & +  Psi%Grid(:,:,:,3,1)*Phi%Grid(:,:,:,4,1)                 &
                  & -  Psi%Grid(:,:,:,4,1)*Phi%Grid(:,:,:,3,1)

    return
  end function MultiplySpinor
  
  pure function RealMultiplySpinor(PSi, Phi) result(RePsiPhi)
    !-------------------------------------------------------------------------
    ! Computes the real part of 
    !       Psi^{dagger} Phi
    ! Not adapted yet for isospin symmetry breaking.
    !-------------------------------------------------------------------------
    class(Spinor),intent(in) :: Psi, Phi
    real(KIND=dp)            :: REPsiPhi(nx,ny,nz)
    real(KIND=dp)            :: GridPsi(nx,ny,nz,4,1), GridPhi(nx,ny,nz,4,1)
   
    GridPsi = Psi%Grid
    GridPhi = Phi%Grid

    RePsiPhi =   sum(GridPsi(:,:,:,:,1) * GridPhi(:,:,:,:,1),4)
    
  end function RealMultiplySpinor
  
  pure function ImagMultiplySpinor(Psi, Phi) result(ImPsiPhi)
    !-------------------------------------------------------------------------
    ! Computes the imaginary part of 
    !       Psi^{dagger} Phi
    ! Not adapted yet for isospin symmetry breaking.
    !-------------------------------------------------------------------------
    class(Spinor),intent(in) :: Psi, Phi
    real(Kind=dp)            :: GridPsi(nx,ny,nz,4,1), GridPhi(nx,ny,nz,4,1)
    real(KIND=dp)            :: ImPsiPhi(nx,ny,nz)

    GridPsi=Psi%Grid
    GridPhi=Phi%Grid

    ImPsiPhi =   GridPsi(:,:,:,1,1)*GridPhi(:,:,:,2,1)                       &
             & - GridPsi(:,:,:,2,1)*GridPhi(:,:,:,1,1)                       &
             & + GridPsi(:,:,:,3,1)*GridPhi(:,:,:,4,1)                       &
             & - GridPsi(:,:,:,4,1)*GridPhi(:,:,:,3,1)
  end function ImagMultiplySpinor
  
  pure function Add(psi, phi) result(PsiplusPhi)
    !-------------------------------------------------------------------------
    ! This function adds two spinors Psi and Phi.
    ! It is used for overloading + with spinors.
    !-------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi, Phi
    type(Spinor)              :: PsiplusPhi

    PsiPlusPhi=NewSPinor()
    PsiplusPhi%Grid = Psi%Grid + Phi%Grid
    
    return
  end function Add
  pure function Subtract(psi,phi) result(PsiMinPhi)
    !-------------------------------------------------------------------------
    ! This function subtracts two spinors Psi and Phi.
    ! It is used for overloading - with spinors.
    !-------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi, Phi
    type(Spinor)              :: PsiMinPhi

    PsiMinPhi = NewSPinor()
    PsiminPhi%Grid = Psi%Grid - Phi%Grid

    return
  end function Subtract
  
  pure function Negative(Psi) result(MinPsi)
    !-------------------------------------------------------------------------
    ! This function returns the negative of Psi.
    ! It is used for overloading - with spinors.
    !-------------------------------------------------------------------------                                              
    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: MinPsi
    
    MinPsi = NewSPinor()
    MinPsi%Grid = -Psi%Grid
    
  end function Negative
                
  pure function Conj(Psi) result(CPsi)
    !-------------------------------------------------------------------------
    ! This function calculates the conjugate of the spinor Psi.
    !-------------------------------------------------------------------------
          
    class(Spinor), intent(in) :: Psi
    real(Kind=dp)             :: ConjSpinor(nx,ny,nz,4,1), GridPsi(nx,ny,nz,4,1)
    type(Spinor)              :: CPsi
    
    CPsi=NewSPinor()
    GridPsi = Psi%Grid                        

    ConJSpinor(:,:,:,1,:) = GridPsi(:,:,:,1,:)
    ConJSpinor(:,:,:,2,:) =-GridPsi(:,:,:,2,:)
    ConJSpinor(:,:,:,3,:) = GridPsi(:,:,:,3,:)
    ConJSpinor(:,:,:,4,:) =-GridPsi(:,:,:,4,:)
    
    CPsi%Grid = ConjSpinor
    return
          
  end function Conj
  
  pure function MultiplyI(Psi) result(iPsi)
    !-------------------------------------------------------------------------
    ! This function calculates the action of the complex unit i on the spinor Psi.
    !-------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: iPsi

    IPsi = NewSPinor()  
    IPsi%Grid(:,:,:,1,:) = - Psi%Grid(:,:,:,2,:)
    IPsi%Grid(:,:,:,2,:) =   Psi%Grid(:,:,:,1,:)
    IPsi%Grid(:,:,:,3,:) = - Psi%Grid(:,:,:,4,:)
    IPsi%Grid(:,:,:,4,:) =   Psi%Grid(:,:,:,3,:)
    return
  end function MultiplyI
  
  pure function TimeReverse(Psi) result(Phi)
    !-------------------------------------------------------------------------
    ! Function that calculates the time reversal of the spinor Phi
    ! Phi = - i \sigma_y K \Psi
    !-------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: Phi
    real*8                    :: Grid1(nx,ny,nz,4,1), Grid2(nx,ny,nz,4,1)
    
    Phi = NewSPinor()
    Grid1 = Psi%Grid
    
    Grid2(:,:,:,1,:) =   Grid1(:,:,:,3,:)
    Grid2(:,:,:,2,:) = - Grid1(:,:,:,4,:)
    Grid2(:,:,:,3,:) = - Grid1(:,:,:,1,:)                
    Grid2(:,:,:,4,:) =   Grid1(:,:,:,2,:)

    Phi%Grid = Grid2

    return
  end function TimeReverse
  
  function ActionOfParity(Phi, Signature) result(Psi)
    !-------------------------------------------------------------------------
    ! Function that returns the spinor correspoding to the action of P on Phi.
    ! P Phi = Phi (-x, -y, -z, sigma)
    !
    ! Note that, in the case of signature conservation, the signature quantum 
    ! number is necessary to construct the action of Parity.
    !-------------------------------------------------------------------------
    class(SPinor), intent(in) :: Phi  
    type(Spinor)            :: Psi
    integer, intent(in)     :: Signature
    integer                 :: i,j,k            

    if(PC) call stp('ActionOfParity should not be used when parity is'         & 
    &               // ' conserved!')

    Psi = NewSPinor()
    if(SC) then
      !Only invert the z-coordinate
      do k=1,nz
          Psi%Grid(:,:,k,:,:) = Phi%Grid(:,:,nz-k+1,:,:)
      enddo
      !Invert the x- and y-coordinate by applying signature operator.
      Psi = dble(Signature)* Pauli(Psi,3)
    else
      if(TSC) then
        !Symmetry for the y-component
        do k=1,nz
          do i=1,nx
            Psi%Grid(i,:,k,:,:) = Phi%Grid(nx-i+1, :, nz-k+1, :,:)
          enddo
        enddo
        Psi = Conj(Psi)
      else
        !No symmetries
        do k=1,nz
          do j=1,ny
            do i=1,nx
              Psi%Grid(i,j,k,:,:) = Phi%Grid(nx-i+1, ny-i+1, nz-i+1, :,:)
            enddo
          enddo
        enddo
      endif
    endif
    return
  end function ActionOfParity
  
  function ActionOfSignature(Phi) result(Psi)
    !-------------------------------------------------------------------------
    ! R_z Phi = -i \hat{sigma}_z Phi(-x,-y,z)
    !
    ! Note that it is only safe to call this routine in the case of breaking of
    ! siganture; i.e. the storing of the entire x-axis.
    !
    !-------------------------------------------------------------------------
    class(Spinor), intent(in) :: Phi  
    type(Spinor)        :: Psi
    integer             :: i,j
    
    Psi = NewSPinor()

    if(SC) call stp('ActionOfSignature should not be used when signature is'   & 
    &               // ' conserved!')

    if(TSC) then
     !Don't invert the y-coordinate in the case of TimeSimplex conservation
     do i=1,nx
       Psi%Grid(i,:,:,:,:) = Phi%Grid(nx-i+1,:,:,:,:)
     enddo
     !Instead apply Timesimplex operator
     Psi = Pauli(Psi,3)
     Psi = MultiplyI(Psi)
     Psi = -Conj(Psi)
     
    else
      do j=1,ny
        do i=1,nx
          Psi%Grid(i,j,:,:,:) = Phi%Grid(nx-i+1,ny-j+1,:,:,:)
        enddo
      enddo
      Psi = Pauli(Psi,3)
      Psi = MultiplyI(Psi)
    endif
  end function ActionOfSignature
  
  function ActionOfTimeSimplex(Phi) result(Psi)
    !-------------------------------------------------------------------------
    ! This function computes the action of the Time Simplex operator on Phi.
    ! S^T_y Phi = Phi^*(x,-y,z,sigma)
    ! Note that it is only safe to call this routine in the case of breaking of
    ! signature; i.e. the storing of the entire x-axis.
    !-------------------------------------------------------------------------      
    class(Spinor), intent(in) :: Phi                    
    type(Spinor) :: Psi
    integer      :: j

    Psi = NewSPinor()                 
    do j=1,ny
        Psi%Grid(:,j,:,:,:) = Phi%Grid(:,ny-j+1,:,:,:)
    enddo            
    Psi = Conj(Psi)

  end function ActionOfTimeSimplex
  function DeriveSpinor(Psi,Parity,Signature,TimeSimplex) result(dPsi)
    !-------------------------------------------------------------------------
    ! This function uses the derivatives module to calculate the derivative of 
    ! the spinor Psi.
    ! Notice that the result is an array of three spinors.
    !       Parity, Signature and Timesimplex indicate the symmetries which the 
    !        spinor possesses.
    !-------------------------------------------------------------------------
    use Derivatives

    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: dPsi(3)
    integer, intent(in)       :: Parity,Signature,TimeSimplex
    real(KIND=dp)             :: DerivResult(nx,ny,nz,3), GridPsi(nx,ny,nz,4,1)
    integer                   :: i,l,k
    
    do i=1,3
         DPsi(i)  =NewSPinor()
    enddo
    GridPsi = Psi%Grid
    do k=1, size(Psi%Grid,5)
      do l=1,4
        DerivResult(:,:,:,1) = &
        & DeriveX(GridPsi(:,:,:,l,k),Parity,Signature,TimeSimplex,l)
        DerivResult(:,:,:,2) = &
        & DeriveY(GridPsi(:,:,:,l,k),Parity,Signature,TimeSimplex,l)
        DerivResult(:,:,:,3) = &
        & DeriveZ(GridPsi(:,:,:,l,k),Parity,Signature,TimeSimplex,l)
        do i=1,3
            dPsi(i)%Grid(:,:,:,l,k) = DerivResult(:,:,:,i)
        enddo 
      enddo
    enddo
    return    
  end function DeriveSpinor
  
  function SecondDerivativeSpinor(Psi,Direction,Parity,Signature,TimeSimplex)  &
  &                              result(dPsi)
  
    use Derivatives

    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: dPsi
    integer, intent(in)       :: Parity,Signature,TimeSimplex, Direction
    real(KIND=dp)             :: DerivResult(nx,ny,nz), GridPsi(nx,ny,nz,4,1)
    integer                   :: l,k
  
    dPsi = NewSPinor()
    
    GridPsi = Psi%Grid
    do k=1, size(Psi%Grid,5)
      do l=1,4
        DerivResult(:,:,:) = SecondDer_Central(GridPsi(:,:,:,l,k),             &
        & Direction,Parity,Signature,TimeSimplex,l)      
        dPsi%Grid(:,:,:,l,k) = DerivResult(:,:,:)
      enddo
    enddo
  end function SecondDerivativeSpinor
  
  function LapSpinor(Psi,Parity,Signature,TimeSimplex) result(LapPsi)
    !---------------------------------------------------------------------------
    ! This function computes the laplacian of the spinor Psi.
    ! Parity, Signature and Timesimplex indicate the symmetries which 
    ! the spinor possesses.
    !
    !---------------------------------------------------------------------------
    
    use Derivatives

    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: LapPsi
    real(Kind=dp)             :: LapGrid(nx,ny,nz,4,1),GridPsi(nx,ny,nz,4,1)
    integer, intent(in)       :: Parity,Signature,TimeSimplex
    integer                   :: i,k
    
    LapPsi = NewSPinor()

    GridPsi = Psi%Grid
    
    do k=1, size(Psi%Grid,5)
      do i=1,4
        LapGrid(:,:,:,i,k) = &
        & Laplacian(GridPsi(:,:,:,i,k),Parity,Signature,TimeSimplex,i)
      enddo
    enddo
    LapPsi%Grid = LapGrid
    return
  end function LapSpinor
  
  pure function InproductSpinorImaginary(Psi,Phi) result(Inproduct)
    !---------------------------------------------------------------------------
    ! This function computes the imaginary part of the inproduct of two spinors,
    ! meaning
    !       int d\vec{r} \Psi^\dagger \Phi
    !---------------------------------------------------------------------------
                    
    class(Spinor), intent(in) :: Psi,Phi
    real(KIND=dp)             :: Inproduct, Grid1(nx,ny,nz,4,1)
    real(KIND=dp)             :: Grid2(nx,ny,nz,4,1)

    Grid1 = Phi%Grid
    Grid2 = Psi%Grid

    Inproduct = dv*sum(Grid1(:,:,:,2,1)*Grid2(:,:,:,1,1) &
              &      - Grid1(:,:,:,1,1)*Grid2(:,:,:,2,1) &
              &      + Grid1(:,:,:,4,1)*Grid2(:,:,:,3,1) &
              &      - Grid1(:,:,:,3,1)*Grid2(:,:,:,4,1))
    
    return
  end function InproductSpinorImaginary
  
  pure function InproductSpinorReal(Psi,Phi) result(Inproduct)
    !---------------------------------------------------------------------------
    ! This function computes the real part of the inproduct of two spinors, 
    !       int d\vec{r} \Psi^\dagger \Phi
    !---------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi,Phi
    real(KIND=dp)             :: Inproduct
    real(KIND=dp)             :: Grid2(nx,ny,nz,4,1), Grid1(nx,ny,nz,4,1)
       
    Grid1 = Psi%Grid
    Grid2 = Phi%Grid                        

    Inproduct = sum(Grid1*Grid2)*dv
    
    return
  end function InproductSpinorReal
  
  pure function InproductSpinor(Psi,Phi) result(Inproduct)
    !---------------------------------------------------------------------------
    ! This function computes the inproduct of two spinors, 
    !       int d\vec{r} \Psi^\dagger \Phi
    !---------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi,Phi
    real(KIND=dp)             :: Inproduct(2)
    
    Inproduct(1) = InproductSpinorReal(Psi,Phi)
    Inproduct(2) = InproductSpinorImaginary(Psi,Phi)
  
  end function InproductSpinor
  
  pure function InproductSpinorComplex(Psi,Phi) result(Inproduct)
    !---------------------------------------------------------------------------
    ! This function computes the inproduct of two spinors, 
    !       int d\vec{r} \Psi^\dagger \Phi
    !---------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi,Phi
    complex(KIND=dp)          :: Inproduct
    
    Inproduct = dcmplx( InproductSpinorReal(Psi,Phi), &
                      & InproductSpinorImaginary(Psi,Phi))
  
  end function InproductSpinorComplex
  
  subroutine ReadSpinor(Psi,unit, SizeX,SizeY,SizeZ, IsoComp)
    !---------------------------------------------------------------------------
    ! Separate input routine for spinors. This is necessary, since the component
    ! of a spinor is allocatable.
    !---------------------------------------------------------------------------
    class(Spinor), intent(inout) :: Psi
    integer, intent(in)          :: unit, SizeX, SizeY, SizeZ, IsoComp
    real(KIND=dp)                :: Temp(SizeX,SizeY,SizeZ,4)
    
    read(unit) Temp
    
    if(.not.allocated(Psi%Grid)) then
            allocate(Psi%Grid(SizeX,SizeY,SizeZ,4,IsoComp))
    endif
    
    Psi%Grid(:,:,:,:,1) = Temp       
    
  end subroutine ReadSpinor
  
  subroutine WriteSpinor(Psi,unit)
    !---------------------------------------------------------------------------
    ! Separate output routine for spinors. This is necessary, since the 
    ! component of a spinor is allocatable.
    !---------------------------------------------------------------------------
    class(Spinor), intent(in)    :: Psi
    integer, intent(in)          :: unit

    if(.not.allocated(Psi%Grid)) then
      call stp("Can't write a nonallocated Spinor to file!")
    endif
    write (unit) Psi%Grid

  end subroutine WriteSpinor   
end module spinors
