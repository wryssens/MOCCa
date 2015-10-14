module CoulombDerivatives
  !-----------------------------------------------------------------------------
  ! This module contains all necessities for taking derivatives in the Coulomb 
  ! problem. Why are they separated from the normal derivatives?
  ! 1) To have routines that can operate on arrays that do are not of the size 
  !    (nx,ny,nz). In principle, this could also be done with the normal        
  !    derivation routines, but this interferes with their optimisation.
  ! 2) To have derivation routines that respect boundary conditions. 
  !    In principle, again, doable with the normal routines but interfering with
  !    optimisation of the code.
  !-----------------------------------------------------------------------------

  use GenInfo
  use Derivatives

  implicit none

contains

  function Gauss_Seidel_Update(Grid,SizeX,SizeY,SizeZ,SourceTerm,x1,x2,        &       
  &  y1,y2,z1,z2,Omega,Signature,Parity, TimeSimplex,dr)                       &
  &    result(Residual)
  !-----------------------------------------------------------------------------
  ! Function that calculates the laplacian for Gauss-Seidel methods in Coulomb. 
  ! It has its own function, since updates are 
  ! done "in place" and not sequential, as in the Laplacian_Central routine.
  !-----------------------------------------------------------------------------
  ! This routine could utilise a colored "ordering" of the grid, but this is not
  ! used at this moment. One should change the value of Colours to 2 to have a 
  ! red-black ordering.
  !-----------------------------------------------------------------------------
  ! W.H. Press et al., 
  !           Numerical Recipes, Cambridge University Press, Third Edition, 2007
  !-----------------------------------------------------------------------------
  integer, intent(in)                :: SizeX,SizeY,SizeZ,x1,x2,y1,y2,z1,z2
  integer, intent(in)                :: Signature,Parity, TimeSimplex
  real(KIND=dp), intent(in)          :: Omega, SourceTerm(SizeX,SizeY,SizeZ)
  real(KIND=dp), intent(inout)       :: Grid(SizeX,SizeY, SizeZ)
  real(KIND=dp), intent(in),optional :: dr
  integer                            :: N,i,j,k,c,Colours=1,m, SignExtension
  real(KIND=dp)                      :: Residual(sizeX, SizeY, SizeZ)
  real(KIND=dp), allocatable         :: LineExtension(:), Coefs(:)
  real(KIND=dp)                      :: FDX,FDY,FDZ, Coef1, Coef2, FDX3
  
  SignExtension=0
  
  !N is the number of points any of the LineExtensions needs to return       
  N = (size(FDCoulomb) - 1)/2
  allocate(LineExtension(N))
  LineExtension = 0.0_dp
  
  !Taking the finite difference coefficients. It seems wasteful to copy them, 
  ! but in this way MOCCa can keep the indices
  ! from -N to N
  allocate(Coefs(-N:N))
  Coefs = FDCoulomb
  
  if(x2.gt.(SizeX-N) .or. Y2.gt.(SizeY-N) .or. z2.gt.(SizeZ-N)) then
    call stp(                                                                  &
    &'Gauss_Seidel_Update cannot be used on grids without boundary conditions.')
  endif
  !Determining the coefficients
  coef1 = 1._dp/6._dp     
  if(present(dr)) then
      coef2 = Coef1*dr**2
  else
      coef2 = Coef1*dx**2
  endif

  !Looping over the different colors.
  do c=0, Colours-1                       
    do k=z1,z2
      do j=y1,y2
        do i=x1,x2
          !Check if this point has the correct "color"
          if(mod(i+j+k,Colours).ne.c) then
                  !If not, pass to the next one
                  cycle
          endif                                
          FDX  = 0.0_dp
          FDY  = 0.0_dp
          FDZ  = 0.0_dp
          FDX3 = 0.0_dp          
          if(i.le.N) then
            !Not necessary to call these routines in this case.
            SignExtension=CompSignExtension(1,Parity,Signature,TimeSimplex,OneI)
            LineExtension = LineExtensionCoulombX &
            & (Grid,N, Parity, TimeSimplex, Signature,j,k)
          endif
          !Doing the Finite-Difference process in the X Direction
          do m=-N,N
            if(i+m.le.ZeroI) then
                    FDX = FDX + SignExtension*Coefs(m)*LineExtension(1 - (m+i))    
            else
                    FDX = FDX + Coefs(m)*Grid(i+m,j,k)
            endif
          enddo          
          if(j.le.N) then
            !Not necessary to call these routines in this case.               
            SignExtension=CompSignExtension(2,Parity,Signature,TimeSimplex,OneI)
            LineExtension = LineExtensionCoulombY &
            & (Grid,N, Parity, TimeSimplex, Signature,i,k)
          endif
          !Doing the Finite Difference Process in the Y Direction
          do m=-N,N
            if(j+m.le.ZeroI) then
              FDY = FDY + SignExtension*Coefs(m)*LineExtension(1- (m+j))
            else
              FDY = FDY + Coefs(m)*Grid(i,j+m,k)
            endif
          enddo
          if(k.le.N) then
            !Not necessary to call these routines in this case.
            SignExtension=CompSignExtension(3,Parity,Signature,TimeSimplex,OneI)
            LineExtension = LineExtensionCoulombZ &
            & (Grid,N, Parity, TimeSimplex, Signature,i,j)
          endif
          !Doing the Finite Difference process in the Z Direction
          do m=-N,N
            if(k+m.le.ZeroI) then
                    FDZ = FDZ + SignExtension*Coefs(m)*LineExtension(1-(m+k))
            else
                    FDZ = FDZ + Coefs(m)*Grid(i,j,m+k)
            endif
          enddo                                        
          
          !Calculating the Residual in any case
          Residual(i,j,k) = &
          & Coef1*(FDX + FDY + FDZ) - Coef2*SourceTerm(i,j,k)
          
          !Applying the correction.
          Grid(i,j,k) = Grid(i,j,k) + Omega*Residual(i,j,k)
        enddo !x loop
      enddo !y loop
    enddo   !z loop        
  enddo   !Color loop
  
  deallocate(LineExtension)
  deallocate(Coefs)
  
  end function Gauss_Seidel_Update


  function CoulombLaplacian(Grid, SizeX,SizeY,SizeZ, Parity, Signature,        &
  &                         TimeSimplex,Component,BC,x1,y1,z1,x2,y2,z2)        &
  & result(Lap)
    !---------------------------------------------------------------------------
    ! This function splits the calculation of the Laplacian into 1D problems, 
    ! using Central_1D. The structure is completely analogous to the Central 
    ! subroutine, only with an added option for conserving boundary conditions 
    ! in a Coulomb calculation.
    !---------------------------------------------------------------------------    
    integer,intent(in)                :: SizeX,SizeY,SizeZ,Parity,Signature
    real(KIND=dp), target, intent(in) :: Grid(SizeX,SizeY,SizeZ)
    integer,intent(in)                :: TimeSimplex,Component
    real(KIND=dp)             :: Lap(SizeX,SizeY,SizeZ),DerX(SizeX,SizeY,SizeZ)
    real(KIND=dp)             :: DerZ(SizeX,SizeY,SizeZ),DerY(SizeX,SizeY,SizeZ)
    real(KIND=dp), allocatable:: LineExtension(:)
    integer       :: i,j,k, SignExtension,N, StartX,EndX,StartY,EndY,StartZ,EndZ
    
    ! Optional variables: if these are present, the laplacian subroutine does 
    ! not touch the points outside [x1-x2, y1-y2, z1-z2].
    ! This is useful for preserving boundary conditions in Coulomb calculations.
    logical, optional, intent(in) :: BC
    integer, optional, intent(in) :: x1,y1,z1,x2,y2,z2
    
    if(.not.allocated(FDCoulomb)) then
      call stp('Laplacian_Central gets called while FDCoef is not allocated.')
    endif
    if(present(BC)) then
      if(.not.present(z2)) then
        call stp("You can't call Laplacian_Central this way")
      endif
      StartX = x1
      EndX   = x2
      StartY = y1
      EndY   = y2
      StartZ = z1
      EndZ   = z2
    else
      !Taking the entire line if no boundary conditions are present
      StartX = 1
      EndX   = SizeX
      StartY = 1
      EndY   = SizeY
      StartZ = 1
      EndZ   = SizeZ 
    endif
    
    !N is the number of points any of the LineExtension needs to return       
    N = (size(FDCoulomb) - 1)/2
    
    allocate(LineExtension(N))
    
    Lap=0.0_dp
    DerX=0.0_dp
    DerY=0.0_dp
    DerZ=0.0_dp
    
    ! X Derivatives
    SignExtension=CompSignExtension(1,Parity,Signature,TimeSimplex,Component)
    do k=StartZ,EndZ
        do j=StartY,EndY
            LineExtension=LineExtensionCoulombX(Grid, N, Parity, TimeSimplex,  &
            &                                   Signature,j,k)
            DerX(:,j,k)=Coulomb_1D(Grid(:,j,k), SizeX,LineExtension,N,         &
            &                                   SignExtension,StartX,EndX)        
        enddo
    enddo
!    if(allocated(LineExtension)) then
!            deallocate(LineExtension) 
!    endif
!   
    ! Y Derivatives
    SignExtension=CompSignExtension(2,Parity,Signature,TimeSimplex,Component)
    do k=StartZ,EndZ
        do i=StartX,EndX
            LineExtension=LineExtensionCoulombY(Grid,N, Parity, TimeSimplex,   &
            &                                   Signature,i,k)
            DerY(i,:,k)=Coulomb_1D(Grid(i,:,k), SizeY,LineExtension,N,         & 
            &                      SignExtension,StartY,EndY)            
        enddo
    enddo
!    if(allocated(LineExtension)) then
!            deallocate(LineExtension)
!    endif

    !Z Derivatives
    SignExtension=CompSignExtension(3,Parity,Signature,TimeSimplex,Component)
    do j=StartY,EndY
        do i=StartX,EndX
            LineExtension=LineExtensionCoulombZ(Grid, N, Parity, TimeSimplex,  &
            &                                   Signature,i,j)
            DerZ(i,j,:)=Coulomb_1D(Grid(i,j,:), SizeZ,LineExtension,N,         &
            SignExtension,StartZ,EndZ)            
        enddo
    enddo

    Lap = DerX+DerY+DerZ    
    return
  end function CoulombLaplacian

  function Coulomb_1D(Line,LineSize,LineExtension,N,SignExtension,S,E)         &
  & result(DerLine)
    !---------------------------------------------------------------------------
    ! This function computes the derivative of second degree of Line.
    ! To do this, one needs the extension of Line: LineExtension.
    ! N is the size of this extension, while SignExtension is the sign of this 
    ! extension.
    ! If LineExtension is completely zero, then there is no relevant symmetry.
    !---------------------------------------------------------------------------
    integer, intent(in)                   :: LineSize,N,SignExtension,S,E 
    !S and E are the first and last points that need to be derived.
    real(KIND=dp), intent(in)             :: Line(LineSize) 
    real(KIND=dp)                         :: DerLine(LineSize), Factor, p
    real(KIND=dp), intent(in)             :: LineExtension(N)    
    integer                               :: i,j
    real(KIND=dp)                         :: Coefs(-N:N)
    Coefs = FDCoulomb
                  
    DerLine=0.0_dp        
    
    p =real(SignExtension, KIND=dp)
    !Using the Finite Difference formula for all the points that are not at the
    ! end-points of the line.
    do i=S,E
      do j=-N,N
          if((i+j).le.ZeroI) then 
            !These are the points at the start of the Line
            DerLine(i) = DerLine(i) + Coefs(j)*p*LineExtension(1 - (i+j))
          elseif((i+j).le.LineSize) then
              !These are the middle points of the line
              DerLine(i) = DerLine(i) + Coefs(j)*Line(i+j)
          !Note that there is no "if" for the endpoints of the line. 
          ! We assume that the wavefunctions are zero out there.
          endif        
      enddo
    enddo
    factor = OneR/(dx**2)
    DerLine=DerLine*factor             
  end function Coulomb_1D

  function LineExtensionCoulombX(Grid,N,Parity,TimeSimplex,Signature,Y,Z)      &
  & result(LineExtensionX)
    !---------------------------------------------------------------------------
    ! This subroutine finds the correct pointers for the line extensions.
    ! The pointers that are unnecessary (meaning not defined by symmetries) are 
    ! all set to zero.
    !---------------------------------------------------------------------------
    integer, intent(in)                :: Y,Z,N
    real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
    integer, intent(in)                :: Parity,Timesimplex,Signature
    real(KIND=dp)                      :: LineExtensionX(N)                  
    integer                            :: SizeY     

    SizeY=size(Grid,2)

    if(Signature.ne.ZeroI) then
      if(TimeSimplex.ne.ZeroI) then
        !Using Signature and Time Simplex
        LineExtensionX=Grid(1:N,Y,Z)
      else
        !Using Signature
        LineExtensionX=Grid(1:N,SizeY-Y+1,Z)
      endif
    else
      LineExtensionX = 0.0_dp
    endif
    return
  end function LineExtensionCoulombX
        
  function LineExtensionCoulombY(Grid,N,Parity,TimeSimplex,Signature,X,Z)      &
  & result(LineExtensionY)
    !---------------------------------------------------------------------------
    ! This subroutine finds the correct pointers for the line extensions.
    ! The pointers that are unnecessary (meaning not defined by symmetries) are 
    ! set to zero.
    !---------------------------------------------------------------------------
    integer, intent(in)                :: X,Z,N
    real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
    integer, intent(in)                :: Parity,Timesimplex,Signature
    real(KIND=dp)                      :: LineExtensionY(N)       
   
    if(TimeSimplex.ne.ZeroI) then
        !Using Time Simplex
        LineExtensionY = Grid(X,1:N,Z) 
    else
        LineExtensionY=0.0_dp
    endif
    return
  end function LineExtensionCoulombY
        
  function LineExtensionCoulombZ(Grid,N,Parity,TimeSimplex,Signature,X,Y)      &
  & result(LineExtensionZ)
    !---------------------------------------------------------------------------
    ! This subroutine finds the correct pointers for the line extensions.
    ! The pointers that are unnecessary (meaning not defined by symmetries) are 
    ! set to zero.
    !---------------------------------------------------------------------------
    integer, intent(in)                :: X,Y,N
    real(KIND=dp), target, intent(in)  :: Grid(:,:,:)
    integer, intent(in)                :: Parity,Timesimplex,Signature
    real(KIND=dp)                      :: LineExtensionZ(N)
    integer                            :: sizeX, SizeY

    SizeX= Size(Grid,1); SizeY=Size(Grid,2)

    if(Parity.ne.ZeroI) then
        if(Signature.ne.ZeroI) then
            ! Using Signature and Parity
            LineExtensionZ = Grid(X,Y,1:N)
        else
            if(TimeSimplex.ne.ZeroI) then
                !Using Parity and TimeSimplex
                LineExtensionZ = Grid(SizeX-X+1, Y, 1:N)
            else
                !Using Parity
                LineExtensionZ = Grid(SizeX-X+1, SizeY-Y+1,1:N)
            endif
        endif
    else
        LineExtensionZ=0.0_dp
    endif

    return
  end function LineExtensionCoulombZ
end module CoulombDerivatives
