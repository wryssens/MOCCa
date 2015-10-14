module MultiGrid

  !-----------------------------------------------------------------------------
  ! This module contains all the routines needed to solve a linear system using 
  ! a multigrid solver.
  !-----------------------------------------------------------------------------
  use CompilationInfo
  use Geninfo
  use CoulombDerivatives
  
  implicit none
          
contains
  subroutine MGSolve(Solution, Source, dr, MaxIteration, Precis)
    !---------------------------------------------------------------------------
    ! This subroutine is a driving routine for the Multigrid solver. Input is
    ! first guess at the Solution, the source term and the grid spacing dr.
    ! It also checks the input Solution: this needs to be a grid with dimension 
    ! [2**j+1, 2**k+1, 2**l+1] with j,k and l integers.
    !---------------------------------------------------------------------------
    
    1 format ("Multigrid Converged after ", i5," iterations.")
    2 format ("Multigrid didn't converge after ", i5,                          &
    &         " iterations. Norm of the Residual = ", e15.8)
    
    real(KIND=dp), intent(inout) :: Solution(:,:,:)
    real(KIND=dp), intent(in)    :: Source(:,:,:), dr, Precis
    integer, intent(in)   :: MaxIteration
    real(KIND=dp), allocatable   :: Residual(:,:,:)
    real(KIND=dp)                :: NormResidual,  Time0, Time1
    integer                      :: Iteration,fx,fy,fz
    
    !Fixing dimensions
    fx = Size(Solution,1)
    fy = Size(Solution,2)
    fz = Size(Solution,3)
    
    !Preparing memory
    allocate(Residual(fx,fy,fz))
    Residual     = ZeroR
    NormResidual = ZeroR
    
    !---------------------------------Checks on the input-----------------------
    if(Size(Source,1).ne.fx) then
            call stp('The first  size of Solution and Source are not the same.')
    endif
    if(Size(Source,2).ne.fy) then
            call stp('The second size of Solution and Source are not the same.')
    endif
    if(Size(Source,3).ne.fz) then
            call stp('The third  size of Solution and Source are not the same.')
    endif
    
    !Checking if fx,fy and fz are odd
    if(mod(fx,2).eq.ZeroI) &
    &  call stp('Multigrid solver cannot take even grid dimensions!', 'fx', fx)
    if(mod(fy,2).eq.ZeroI)&
    &  call stp('Multigrid solver cannot take even grid dimensions!', 'fy', fy)
    if(mod(fz,2).eq.ZeroI) &
    & call stp('Multigrid solver cannot take even grid dimensions!', 'fz', fz)

  !  !Checking if fx is some power of 2 + 1
  !  Test = 1
  !  do while (Test.lt.(fx-1)) 
  !          Test = Test*2
  !  enddo
  !  if(Test.ne.(fx-1)) &
  !  & call stp('The first dimension is not equal to "2^j + 1" for some j')
  !  !Checking if fy is some power of 2 + 1
  !  Test = 1
  !  do while (Test.lt.(fy-1)) 
  !          Test = Test*2
  !  enddo
  !  if(Test.ne.(fy-1)) &
  !  & call stp('The second dimension is not equal to "2^j + 1" for some j')
  !  !Checking if fz is some power of 2 + 1
  !  Test = 1
  !  do while (Test.lt.(fz-1)) 
  !          Test = Test*2
  !  enddo
  !  if(Test.ne.(fz-1)) &
  !  & call stp('The third dimension is not equal to "2^j + 1" for some j')
  !-----------------------------------------------------------------------------
      
      call Cpu_time(Time0)
      do Iteration=1,MaxIteration
        !Perform 1 Vcycle
        Residual = VCycle(Solution, Source,dr)
        
        !Compute the norm of the residual
        NormResidual = sum(Residual**2)
        print*, Iteration, NormResidual
        !Checking for convergence
        if(NormResidual.le.Precis) then
                print 1, Iteration
                exit
        elseif(Iteration.eq.MaxIteration) then
                print 2, Iteration, NormResidual
        endif
      enddo
      call CPu_Time(Time1)
      print*, "Actual Time", Time1-Time0
      return
  end subroutine MGSolve
              
  recursive function VCycle(Solution_f, Source_f, dr_f) result(Residual)
    real(KIND=dp), intent(inout) :: Solution_f(:,:,:)
    real(KIND=dp), intent(in)    :: Source_f(:,:,:), dr_f
    real(KIND=dp), allocatable   :: Residual(:,:,:), Source_c(:,:,:)
    real(KIND=dp), allocatable   :: Solution_c(:,:,:), Corr_c(:,:,:)
    real(KIND=dp)                :: dr_c, Omega=1.0_dp,NormResidual
    integer               :: fx,fy,fz,cx,cy,cz, PreIter,PostIter, Iteration
    logical               :: DimEven
    
    !Fixing the dimensions
    fx = Size(Solution_f,1)
    fy = Size(Solution_f,2)
    fz = Size(Solution_f,3)
    
    cx = (fx-1)/2 + 1
    cy = (fy-1)/2 + 1
    cz = (fz-1)/2 + 1
    !This is neccessary in every case
    allocate(Residual(fx,fy,fz))
    
    !Is one of the dimensions of the fine grid even?
    DimEven=(mod(fx,2).eq.ZeroI .or. mod(fy,2).eq.ZeroI .or. mod(fz,2).eq.ZeroI)
    
    
    !Checking if we are at the coarsest level or not. For the moment.
    ! Ofcourse, we cannot go deeper if one of the dimensions is even.
     if((min(cx,cy,cz).gt.3).and.(.not.DimEven )) then

      !Allocating the extra necessary memory                        
      allocate(Source_c(cx,cy,cz))
      allocate(Solution_c(cx,cy,cz))
      allocate(Corr_c(cx,cy,cz))
      
      !The coarse grid size will be twice the fine grid size
      dr_c = 2*dr_f
      
      PreIter = 2
      PostIter= 2
      
      !Presmoothing
      do Iteration = 1,PreIter
        Residual = Gauss_Seidel_Update(Solution_f,fx,fy,fz,Source_f,&
       & 1,fx-1,1,fy-1,1,fz-1,Omega,SignatureInt,ParityInt, TimeSimplexInt,dr_f)
        ! Note that this is only viable with one boundary condition. 
        ! Generalisation to laplacians with more points
        ! are not yet supported.
      enddo
      
      !Cycle up to a coarser level
      Source_c   = Restrict(Residual)
      Solution_c = ZeroR
      
      !Solve the coarser level
      Corr_c = VCycle(Solution_c, Source_c, dr_c) !Recursive call
      
      !Correct the fine level
      Solution_c = Solution_c - Prolong(Corr_c)
      Residual = ZeroR
      !PostSmoothing
      do Iteration = 1,PostIter
          Residual = Gauss_Seidel_Update(Solution_f,fx,fy,fz,Source_f,&
      & 1,fx-1,1,fy-1,1,fz-1,Omega,SignatureInt,ParityInt, TimeSimplexInt,dr_f)
              ! Note that this is only viable with one boundary condition. 
              ! Generalisation to laplacians with more points
              ! are not yet supported.
      enddo
      
      !deallocating everything
      deallocate(Source_c, Solution_c,Corr_c)
      
    else ! Here we are at the coarsest level.
     print*, "Deepest Level attained"
     print*, fx,fy,fz
     print*, ""
     !We solve the coarsest level using Gauss-Seidel
     Residual=ZeroR
     NormResidual = OneR
     do while(NormResidual .gt. 1e-5) 
          Residual = Gauss_Seidel_Update(Solution_f,fx,fy,fz,Source_f,&
      & 1,fx-1,1,fy-1,1,fz-1,Omega,SignatureInt,ParityInt, TimeSimplexInt,dr_f)
          NormResidual = sum(Residual**2)      
     enddo
    endif
    return
  end function VCycle
        
  function Restrict(Fine) result(Coarse)
    !---------------------------------------------------------------------------
    ! This function takes a grid Fine and restricts it with direct injection to 
    ! a grid Coarse. 
    ! Fine(1:a,1:b,1:c) => Coarse(1:(a-1)/2 + 1, 1:(b-1)/2 + 1,1:(c-1)/2 + 1)
    !
    ! This is the 3D Equivalent of the following process
    !
    !  x     x     x     x     x           x     o     x     o     x
    !  x     x     x     x     x           o     x     o     x     o
    !  x     x     x     x     x     =>    x     o     x     o     x
    !  x     x     x     x     x           o     x     o     x     o
    !  x     x     x     x     x           x     o     x     o     x
    ! where the o's are deleted.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: Fine(:,:,:)
    real(KIND=dp),allocatable :: Coarse(:,:,:)
    integer            :: fx,fy,fz,cx,cy,cz,i,j,k
    
    ! Fixing the dimensions of Fine and Coarse
    fx = Size(Fine,1)
    fy = Size(Fine,2)
    fz = Size(Fine,3)
    cx = (fx - 1)/2 + 1
    cy = (fy - 1)/2 + 1
    cz = (fz - 1)/2 + 1
    
    allocate(Coarse(cx,cy,cz))
    
    do k=1,cz
            do j=1,cy
                    do i=1,cx
                            Coarse(i,j,k) = Fine(2*i-1,2*j-1,2*k-1)
                    enddo
            enddo
    enddo
  end function Restrict

  function Prolong(Coarse) result(Fine)
    !---------------------------------------------------------------------------
    ! This function takes a grid Coarse and prolongs it with direct injection 
    ! and some interpolation to a grid Fine. 
    ! Coarse(1:a,1:b,1:c) => Fine(1:(a-1)*2 + 1, 1:(b-1)*2 + 1,1:(c-1)*2 + 1)
    !       x     i     x     i     x           x     o     x     o     x
    !       i     x     i     x     i           o     x     o     x     o
    !       x     i     x     i     x     <=    x     o     x     o     x
    !       i     x     i     x     i           o     x     o     x     o
    !       x     i     x     i     x           x     o     x     o     x
    !
    ! where the 'i's have been interpolated between neighbouring points.
    !---------------------------------------------------------------------------
                
    real(KIND=dp),allocatable :: Fine(:,:,:)
    real(KIND=dp),intent(in)  :: Coarse(:,:,:)
    integer            :: fx,fy,fz,cx,cy,cz,i,j,k

    !Fixing the dimensions
    cx = Size(Coarse,1)
    cy = Size(Coarse,2)
    cz = Size(Coarse,3)
    
    fx = cx*2 - 1
    fy = cy*2 - 1
    fz = cz*2 - 1
    
    allocate(Fine(fx,fy,fz))
          
    !Filling in the values for the non-boundary points
    do k=1,cz-1
      do j=1,cy-1
        do i=1,cx-1
          !Straight Injection
          Fine(2*i-1,2*j-1,2*k-1) = Coarse(i,j,k)
          !Interpolation for t
          Fine(2*i-1,2*j-1,2*k)   = 0.5_dp*(Coarse(i,j,k) + Coarse(i,j,k+1))
          Fine(2*i-1,2*j,2*k-1)   = 0.5_dp*(Coarse(i,j,k) + Coarse(i,j+1,k))
          Fine(2*i,2*j-1,2*k-1)   = 0.5_dp*(Coarse(i,j,k) + Coarse(i+1,j,k))
          Fine(2*i,2*j,2*k-1)     =1/3._dp*(Coarse(i,j,k) + Coarse(i+1,j,k)    &
          &                       + Coarse(i,j+1,k))
          Fine(2*i,2*j-1,2*k)     =1/3._dp*(Coarse(i,j,k) + Coarse(i+1,j,k)    &
          &                       + Coarse(i,j,k+1))
          Fine(2*i-1,2*j,2*k)     =1/3._dp*(Coarse(i,j,k) + Coarse(i,j+1,k)    &
          &                       + Coarse(i,j,k+1))
          Fine(2*i,2*j,2*k)       =0.25_dp* &
          & (Coarse(i,j,k) + Coarse(i,j+1,k) + Coarse(i,j,k+1)+Coarse(i+1,j,k))
        enddo
      enddo
    enddo
          
    !Filling in the values for the Boundary points
    !Z Boundary
    do j=1,cy-1
      do i=1,cx-1
        !Straight Injection
        Fine(2*i-1,2*j-1,fz) = Coarse(i,j,cz)
        !Interpolation
        Fine(2*i-1,2*j,fz)   = 0.5_dp*(Coarse(i,j,cz) + Coarse(i,j+1,cz))
        Fine(2*i,2*j-1,fz)   = 0.5_dp*(Coarse(i,j,cz) + Coarse(i+1,j,cz))
        Fine(2*i,2*j,  fz)   =1/3._dp*(Coarse(i,j,cz) + Coarse(i+1,j,cz)       &
        &                    + Coarse(i,j+1,cz))
      enddo
    enddo
    !Y Boundary
    do k=1,cz-1
      do i=1,cx-1
        !Straight Injection
        Fine(2*i-1,fy,2*k-1) = Coarse(i,cy,k)
        !Interpolation
        Fine(2*i-1,fy,2*k)   = 0.5_dp*(Coarse(i,cy,k) + Coarse(i,cy,k+1))
        Fine(2*i,fy,2*k-1)   = 0.5_dp*(Coarse(i,cy,k) + Coarse(i+1,cy,k))
        Fine(2*i,fy,2*k)     =1/3._dp*(Coarse(i,cy,k) + Coarse(i+1,cy,k)       &
        &                    + Coarse(i,cy,k+1))
      enddo
    enddo
    !X Boundary
    do k=1,cz-1
      do j=1,cy-1
        !Straight Injection
        Fine(fx,2*j-1,2*k-1) = Coarse(fx,j,k)
        !Interpolation
        Fine(fx,2*j-1,2*k)   = 0.5_dp*(Coarse(cx,j,k) + Coarse(cx,j,k+1))
        Fine(fx,2*j,2*k-1)   = 0.5_dp*(Coarse(cx,j,k) + Coarse(cx,j+1,k))
        Fine(fx,2*j,2*k)     =1/3._dp*(Coarse(cx,j,k) + Coarse(cx,j+1,k)       &
        & + Coarse(cx,j,k+1))
      enddo
    enddo
    
    !X & Y Boundary
    do k=1,cz-1
            !Straight Injection
            Fine(fx,fy,2*k-1) = Coarse(cx,cy,k)
            !Interpolation
            Fine(fx,fy,2*k)   = 0.5_dp*(Coarse(cx,cy,k) + Coarse(cx,cy,k+1))
    enddo
    
    !X & Z Boundary
    do j=1,cy-1
            !Straight Injection
            Fine(fx,2*j-1,fz) = Coarse(cx,j,cz)
            !Interpolation
            Fine(fx,2*j,fz)   = 0.5_dp*(Coarse(cx,j,cz) + Coarse(cx,j+1,cz))
    enddo
    
    !Y & Z Boundary
    do i=1,cx-1
            !Straight Injection
            Fine(2*i-1,fy,fz) = Coarse(i,cy,cz)
            !Interpolation
            Fine(2*i,fy,fz)   = 0.5_dp*(Coarse(i,cy,cz) + Coarse(I+1,cy,cz))
    enddo
    
    !X,Y,Z Boundary
    Fine(fx,fy,fz) = Coarse(cx,cy,cz)
          
    return
  end function Prolong

end module MultiGrid
