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

  pure function NewSpinor() result(Psi)
    !---------------------------------------------------------------------------
    ! This function serves as constructor for the Spinor Class, taking 
    ! proper care of allocation and initialisation to zero.
    !---------------------------------------------------------------------------
    ! Strictly speaking this is inefficient, as the setting to zero takes some
    ! time. In time-critical places, a direct allocation is a better strategy.
    !---------------------------------------------------------------------------
    type(Spinor) :: Psi

    allocate(Psi%Grid(nx,ny,nz,4,1)); Psi%Grid=0.0_dp
  end function
  
  pure function NewSizeSpinor( insizex, insizey, insizez) result (Psi)
    !---------------------------------------------------------------------------
    ! COnstructor for a spinor with different mesh sizes, very practical in 
    ! the transform module.
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
  end function NewSizeSpinor
  
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
    integer                  :: i
                    
    if(Direction.eq.1) then
        !\sigma_x = ( 0  1 )
        !           ( 1  0 )
        do i=1,nx*ny*nz
          PauliSpinor(i,1,1,1,1) = Psi%Grid(i,1,1,3,1)
          PauliSpinor(i,1,1,2,1) = Psi%Grid(i,1,1,4,1)
          PauliSpinor(i,1,1,3,1) = Psi%Grid(i,1,1,1,1)
          PauliSpinor(i,1,1,4,1) = Psi%Grid(i,1,1,2,1)   
        enddo      
    elseif(Direction.eq.2) then                
        !\sigma_y = ( 0 -i )
        !           ( i  0 )
        do i=1,nx*ny*nz          
          PauliSpinor(i,1,1,1,1) =   Psi%Grid(i,1,1,4,1)
          PauliSpinor(i,1,1,2,1) = - Psi%Grid(i,1,1,3,1)
          PauliSpinor(i,1,1,3,1) = - Psi%Grid(i,1,1,2,1)
          PauliSpinor(i,1,1,4,1) =   Psi%Grid(i,1,1,1,1)
        enddo       
    elseif(Direction.eq.3) then
        !\sigma_z = ( 1  0 )
        !           ( 0 -1 )
        do i=1,nx*ny*nz       
          PauliSpinor(i,1,1,1,1) =   Psi%Grid(i,1,1,1,1)
          PauliSpinor(i,1,1,2,1) =   Psi%Grid(i,1,1,2,1)
          PauliSpinor(i,1,1,3,1) = - Psi%Grid(i,1,1,3,1)
          PauliSpinor(i,1,1,4,1) = - Psi%Grid(i,1,1,4,1)
        enddo
    endif
    allocate(SigmaPsi%Grid(nx,ny,nz,4,1))
    do i=1,nx*ny*nz*4
      SigmaPsi%Grid(i,1,1,1,1) = PauliSpinor(i,1,1,1,1)
    enddo
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
    integer                   :: i,k
    
    allocate(FPsi%Grid(nx,ny,nz,4,1))
    do k=1,4
      do i=1,nx*ny*nz
        FPsi%Grid(i,1,1,k,1) = F(i,1,1)*Psi%Grid(i,1,1,k,1)
      enddo
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
    integer                   :: i

    allocate(SPsi%Grid(nx,ny,nz,4,1)) 
    !No setting to zero as this is unneccessary and slow down...
    do i=1,4*nx*ny*nz
       SPsi%Grid(i,1,1,1,1) = S*Psi%Grid(i,1,1,1,1)
    enddo  
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
    SPsi = DBLE(S) * Psi + Imag(S) * SPsi
    
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
    integer                  :: i
                            
    PsiPhi = 0.0_dp
    
    !Real part
    do i=1,nx*ny*nz
      PsiPhi(i,1,1,1) = Psi%Grid(i,1,1,1,1) * Phi%Grid(i,1,1,1,1) &
      &               + Psi%Grid(i,1,1,2,1) * Phi%Grid(i,1,1,2,1) &
      &               + Psi%Grid(i,1,1,3,1) * Phi%Grid(i,1,1,3,1) &
      &               + Psi%Grid(i,1,1,4,1) * Phi%Grid(i,1,1,4,1)
    enddo
    !Imaginary Part
    do i=1,nx*ny*nz
      PsiPhi(i,1,1,1) = Psi%Grid(i,1,1,1,1) * Phi%Grid(i,1,1,2,1) &
      &               - Psi%Grid(i,1,1,2,1) * Phi%Grid(i,1,1,1,1) &
      &               + Psi%Grid(i,1,1,3,1) * Phi%Grid(i,1,1,4,1) &
      &               - Psi%Grid(i,1,1,4,1) * Phi%Grid(i,1,1,3,1)
    enddo


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
    integer                  :: i
  
    do i=1,nx*ny*nz
      RePsiPhi(i,1,1) = Psi%Grid(i,1,1,1,1) * Phi%Grid(i,1,1,1,1) &
      &               + Psi%Grid(i,1,1,2,1) * Phi%Grid(i,1,1,2,1) &
      &               + Psi%Grid(i,1,1,3,1) * Phi%Grid(i,1,1,3,1) &
      &               + Psi%Grid(i,1,1,4,1) * Phi%Grid(i,1,1,4,1)
    enddo

    !RePsiPhi =   sum(GridPsi(:,:,:,:,1) * GridPhi(:,:,:,:,1),4)
    
  end function RealMultiplySpinor
  
  pure function ImagMultiplySpinor(Psi, Phi) result(ImPsiPhi)
    !-------------------------------------------------------------------------
    ! Computes the imaginary part of 
    !       Psi^{dagger} Phi
    ! Not adapted yet for isospin symmetry breaking.
    !-------------------------------------------------------------------------
    class(Spinor),intent(in) :: Psi, Phi
    real(KIND=dp)            :: ImPsiPhi(nx,ny,nz)
    integer                  :: i

    do i=1,nx*ny*nz
      ImPsiPhi(i,1,1) =  Psi%Grid(i,1,1,1,1) * Phi%Grid(i,1,1,2,1) &
      &                - Psi%Grid(i,1,1,2,1) * Phi%Grid(i,1,1,1,1) &
      &                + Psi%Grid(i,1,1,3,1) * Phi%Grid(i,1,1,4,1) &
      &                - Psi%Grid(i,1,1,4,1) * Phi%Grid(i,1,1,3,1)
    enddo

  end function ImagMultiplySpinor
  
  pure function Add(psi, phi) result(PsiplusPhi)
    !-------------------------------------------------------------------------
    ! This function adds two spinors Psi and Phi.
    ! It is used for overloading + with spinors.
    !-------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi, Phi
    type(Spinor)              :: PsiplusPhi
    integer                   :: i
    real(KIND=dp)             :: temp(nx,ny,nz,4,1), tempPsi(nx,ny,nz,4,1)
    real(KIND=dp)             :: tempPhi(nx,ny,nz,4,1)

    allocate(PsiPlusPhi%Grid(nx,ny,nz,4,1))
    do i=1,4*nx*ny*nz
      PsiplusPhi%Grid(i,1,1,1,1) = Psi%Grid(i,1,1,1,1) + Phi%Grid(i,1,1,1,1)
    enddo
    
    return
  end function Add
  pure function Subtract(psi,phi) result(PsiMinPhi)
    !-------------------------------------------------------------------------
    ! This function subtracts two spinors Psi and Phi.
    ! It is used for overloading - with spinors.
    !-------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi, Phi
    type(Spinor)              :: PsiMinPhi
    integer                   :: i

    allocate(PsiMinPhi%Grid(nx,ny,nz,4,1))
    do i=1,nx*ny*nz*4
      PsiminPhi%Grid(i,1,1,1,1) = Psi%Grid(i,1,1,1,1) - Phi%Grid(i,1,1,1,1)
    enddo

    return
  end function Subtract
  
  pure function Negative(Psi) result(MinPsi)
    !-------------------------------------------------------------------------
    ! This function returns the negative of Psi.
    ! It is used for overloading - with spinors.
    !-------------------------------------------------------------------------                                              
    class(Spinor), intent(in) :: Psi
    type(Spinor)              :: MinPsi
    
    allocate(MinPsi%Grid(nx,ny,nz,4,1))
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
    integer                   :: i

    allocate(IPsi%Grid(nx,ny,nz,4,1))
    do i=1,nx*ny*nz
      IPsi%Grid(i,1,1,1,1) = - Psi%Grid(i,1,1,2,1)
      IPsi%Grid(i,1,1,2,1) =   Psi%Grid(i,1,1,1,1)
      IPsi%Grid(i,1,1,3,1) = - Psi%Grid(i,1,1,4,1)
      IPsi%Grid(i,1,1,4,1) =   Psi%Grid(i,1,1,3,1)
    enddo
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
    integer                   :: i
    
    allocate(Phi%Grid(nx,ny,nz,4,1))
    do i=1,nx*ny*nz
      Phi%Grid(i,1,1,1,1) =   Psi%Grid(i,1,1,3,1)
      Phi%Grid(i,1,1,2,1) = - Psi%Grid(i,1,1,4,1)
      Phi%Grid(i,1,1,3,1) = - Psi%Grid(i,1,1,1,1)                
      Phi%Grid(i,1,1,4,1) =   Psi%Grid(i,1,1,2,1)
    enddo
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
    integer                   :: i,l,k,n
    
    do i=1,3
         allocate(DPsi(i)%Grid(nx,ny,nz,4,1))
    enddo
    do k=1, size(Psi%Grid,5)
      do l=1,4
        dPsi(1)%Grid(:,:,:,l,k) = &
        & DeriveX(Psi%Grid(:,:,:,l,k),Parity,Signature,TimeSimplex,l)
        dPsi(2)%Grid(:,:,:,l,k) = &
        & DeriveY(Psi%Grid(:,:,:,l,k),Parity,Signature,TimeSimplex,l)
        dPsi(3)%Grid(:,:,:,l,k) = &
        & DeriveZ(Psi%Grid(:,:,:,l,k),Parity,Signature,TimeSimplex,l)
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
    
    do k=1, size(Psi%Grid,5)
      do l=1,4
        dPsi%Grid(:,:,:,l,k) = SecondDer_Central(Psi%Grid(:,:,:,l,k),             &
        & Direction,Parity,Signature,TimeSimplex,l)      
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
    
    allocate(LapPsi%Grid(nx,ny,nz,4,1))

    GridPsi = Psi%Grid
    !do k=1, size(Psi%Grid,5)
      do i=1,4
        LapPsi%Grid(:,:,:,i,1) = &
        & Laplacian(GridPsi(:,:,:,i,1),Parity,Signature,TimeSimplex,i)
      enddo
    !enddo
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
    integer                   :: i 

    Grid1 = Phi%Grid
    Grid2 = Psi%Grid

    Inproduct=0.0_dp
    do i=1,nx*ny*nz
        Inproduct = Inproduct + Grid1(i,1,1,2,1)*Grid2(i,1,1,1,1) &
                       &      - Grid1(i,1,1,1,1)*Grid2(i,1,1,2,1) &
                       &      + Grid1(i,1,1,4,1)*Grid2(i,1,1,3,1) &
                       &      - Grid1(i,1,1,3,1)*Grid2(i,1,1,4,1)
    enddo
    Inproduct = Inproduct * dv
    return
  end function InproductSpinorImaginary
  
  pure function InproductSpinorReal(Psi,Phi) result(Inproduct)
    !---------------------------------------------------------------------------
    ! This function computes the real part of the inproduct of two spinors, 
    !       int d\vec{r} \Psi^\dagger \Phi
    !---------------------------------------------------------------------------
    class(Spinor), intent(in) :: Psi,Phi
    real(KIND=dp)             :: Inproduct
    integer                   :: i

    Inproduct = 0.0_dp
    do i=1,nx*ny*nz*4
      Inproduct = Inproduct + Psi%Grid(i,1,1,1,1)*Phi%Grid(i,1,1,1,1)
    enddo
    Inproduct = Inproduct*dv
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
