module Interpolation
!-------------------------------------------------------------------------------
! Module that contains the routines to interpolate any function on a Lagrange
! mesh defined by the set (mx,my,mz,ex) to a new function determined by
! (nx,ny,nz,dx).
! Note that the interpolation is thus always done from a different Lagrange mesh
! to the runtime parameters of MOCCa.
!-------------------------------------------------------------------------------

  use GenInfo

  implicit none

  real(KIND=dp), allocatable :: interpolX(:,:,:), interpolY(:,:,:), interpolZ(:,:,:)

contains

  subroutine ConstructInterpolationFunctions(mx,my,mz,ex)
    !---------------------------------------------------------------------------
    ! Construct the values of the interpolation functions f_r associated with
    ! the old mesh at the points of the new mesh.
    ! (mx,my,mz, ex) => (nx,ny,nz,dx)
    !
    !
    ! Currently only valid when conserving all of the spatial symmetries
    ! 
    !---------------------------------------------------------------------------
    
    integer, intent(in)       :: mx,my,mz
    real(KIND=dp), intent(in) :: ex
    integer                   :: i,j,k
    real(KIND=dp)             :: fac, c,d,x1,x2, dh, ph
    real(KIND=dp), parameter  :: eps = 1d-3
    
    allocate(interpolX(nx,mx,2), interpolY(ny,my,2), interpolZ(nz,mz,2))
    InterpolX = 0.0d0
    InterpolY = 0.0d0
    InterpolZ = 0.0d0
    
    dh = dx/ex
    fac =   0.5d0 / mx
    x1  = - 0.5d0 * dh
    ph  = 0.5 * pi/mx
    do i=1,nx
      x1 = x1 + dh
      x2  = - 0.5d0
      do j=1,mx        
        x2 = x2 + 1
        if (abs(x1-x2).le.eps) then
          c = pi/ph
        else
          c = sin(pi * (x1 - x2))/sin(ph*(x1-x2))
        endif
        if (abs((x1+x2)/(2*mx)-1.0_dp).le.eps) then
          d= pi/ph
        else
          d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
        endif
        InterpolX(i,j,1) = fac * (c - d)
        InterpolX(i,j,2) = fac * (c + d)
      enddo
    enddo
    
    fac =   0.5d0 / my
    x1  = - 0.5d0 * dh
    ph  = 0.5 * pi/my
    do i=1,ny
      x1 = x1 + dh
      x2  = - 0.5d0
      do j=1,my        
        x2 = x2 + 1
        if (abs(x1-x2).le.eps) then
          c = pi/ph
        else
          c = sin(pi * (x1 - x2))/sin(ph*(x1-x2))
        endif
        if (abs((x1+x2)/(2*my)-1.0_dp).le.eps) then
          d= pi/ph
        else
          d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
        endif
        InterpolY(i,j,1) = fac * (c - d)
        InterpolY(i,j,2) = fac * (c + d)
      enddo
    enddo
    
    fac =   0.5d0 / mz
    x1  = - 0.5d0 * dh
    ph  =   0.5 * pi/mz
    do i=1,nz
      x1 = x1 + dh
      x2 = - 0.5d0
      do j=1,mz        
        x2 = x2 + 1 
        if (abs(x1-x2).le.eps) then
          c = pi/ph
        else
          c = sin(pi*(x1-x2))/sin(ph*(x1-x2))
        endif
        if (abs((x1+x2)/(2*mz)-1.0_dp).le.eps) then
          d= pi/ph
        else
          d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
        endif
        InterpolZ(i,j,1) = fac * (c - d)
        InterpolZ(i,j,2) = fac * (c + d)
      enddo
    enddo
  
  end subroutine ConstructInterpolationFunctions  

end module Interpolation
