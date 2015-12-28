module OptimizedDerivatives
    !-----------------------------------------------------------------------
    ! Module containing derivative functions for specialized symmetry cases, 
    ! as I find the energy to make them.
    !
    !
    !-----------------------------------------------------------------------

    use CompilationInfo
    use Geninfo

    implicit none

contains
!----------------------------------------------------------------------------
! EV8 mode: all symmetries conserved and no negative signature.
!    *) 3rd order FD for the first order derivatives
!    *) 4rd order FD for the laplacian
!

  function Opt_X_EV8(Grid,Parity,Signature, TimeSimplex,Component) result(Der)
    
    integer,intent(in) :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    integer                           :: i,j,k, S

    allocate(Der(nx,ny,nz)) ; Der = 0.0_dp
    select case (Component)
    case(1)
        S =  Signature
    case(2)
        S = -Signature
    case(3)
        S = -Signature
    case(4)
        S =  Signature
    end select
    call derx(Grid, der, S)

  end function Opt_X_EV8

  subroutine derx (w,wx,ix)

      integer, intent(in):: ix
      integer            :: i,j,k
      real(KIND=dp), intent(in) ::w(nx,ny,nz)
      real(KIND=dp), intent(out)::wx(nx,ny,nz)
      real(KIND=dp)           :: ccdx  
      real(KIND=dp),parameter :: ca=45.0d0,cb=9.0d0,cc=60.0d0

      if (ix.gt.0) then
        do j=1,ny*nz
          wx(1,j,1)    = ca*(w(2,j,1)    - w(1,j,1)   ) &
     &                  -cb*(w(3,j,1)    - w(2,j,1)   ) &
     &                     +(w(4,j,1)    - w(3,j,1)   ) 
          wx(2,j,1)    = ca*(w(3,j,1)    - w(1,j,1)   ) &
     &                  -cb*(w(4,j,1)    - w(1,j,1)   ) &
     &                     +(w(5,j,1)    - w(2,j,1)   ) 
          wx(3,j,1)    = ca*(w(4,j,1)    - w(2,j,1)   ) &
     &                  -cb*(w(5,j,1)    - w(1,j,1)   ) &
     &                     +(w(6,j,1)    - w(1,j,1)   ) 
          wx(nx-2,j,1) = ca*(w(nx-1,j,1) - w(nx-3,j,1)) &
     &                  -cb*(w(nx,j,1)   - w(nx-4,j,1)) &
     &                     +(            - w(nx-5,j,1)) 
          wx(nx-1,j,1) = ca*(w(nx,j,1)   - w(nx-2,j,1)) &
     &                  -cb*(            - w(nx-3,j,1)) &
     &                     +(            - w(nx-4,j,1)) 
          wx(nx,j,1)   = ca*(            - w(nx-1,j,1)) &
     &                  -cb*(            - w(nx-2,j,1)) &
     &                     +(            - w(nx-3,j,1))
        enddo
      else
        do j=1,ny*nz
          wx(1,j,1)    = ca*(w(2,j,1)    + w(1,j,1)   ) &
     &                  -cb*(w(3,j,1)    + w(2,j,1)   ) &
     &                     +(w(4,j,1)    + w(3,j,1)   )
          wx(2,j,1)    = ca*(w(3,j,1)    - w(1,j,1)   ) &
     &                  -cb*(w(4,j,1)    + w(1,j,1)   ) &
     &                     +(w(5,j,1)    + w(2,j,1)   )
          wx(3,j,1)    = ca*(w(4,j,1)    - w(2,j,1)   ) &
     &                  -cb*(w(5,j,1)    - w(1,j,1)   ) &
     &                     +(w(6,j,1)    + w(1,j,1)   )
          wx(nx-2,j,1) = ca*(w(nx-1,j,1) - w(nx-3,j,1)) &
     &                  -cb*(w(nx,j,1)   - w(nx-4,j,1)) &
     &                     +(            - w(nx-5,j,1))
          wx(nx-1,j,1) = ca*(w(nx,j,1)   - w(nx-2,j,1)) &
     &                  -cb*(            - w(nx-3,j,1)) &
     &                     +(            - w(nx-4,j,1))
          wx(nx,j,1)   = ca*(            - w(nx-1,j,1)) &
     &                  -cb*(            - w(nx-2,j,1)) &
     &                     +(            - w(nx-3,j,1))
        enddo
      endif

      do i=4,nx-3
      do j=1,ny*nz
        wx(i,j,1)    = ca*(w(i+1,j,1) - w(i-1,j,1) )  &
     &                -cb*(w(i+2,j,1) - w(i-2,j,1) )  &
     &                   +(w(i+3,j,1) - w(i-3,j,1) )
      enddo
      enddo

      ccdx = 1.0_dp / (cc*dx)
      do i=1,nx*ny*nz
        wx(i,1,1) = wx(i,1,1) * ccdx
      enddo

      return
  end subroutine derx

  function Opt_Y_EV8(Grid,Parity,Signature, TimeSimplex,Component) result(Der)
    
    integer,intent(in) :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    integer                           :: i,j,k, S

    allocate(Der(nx,ny,nz)) ; Der = 0.0_dp
    select case (Component)
    case(1)
        S =  TimeSimplex
    case(2)
        S = -TimeSimplex
    case(3)
        S =  TimeSimplex
    case(4)
        S = -TimeSimplex
    end select
    call dery(Grid, der, S)

  end function Opt_Y_EV8

  subroutine dery (w,wy,iy)

      integer, intent(in):: iy
      integer            :: i,j,k
      real(KIND=dp), intent(in) ::w(nx,ny,nz)
      real(KIND=dp), intent(out)::wy(nx,ny,nz)
      real(KIND=dp)             :: ccdx

      real(KIND=dp),parameter :: ca=45.0d0,cb=9.0d0,cc=60.0d0

      if (iy.gt.0) then
        do k=1,nz
        do i=1,nx
          wy(i,1,k)    = ca*(w(i,2,k)    - w(i,1,k)   ) &
     &                  -cb*(w(i,3,k)    - w(i,2,k)   ) &
     &                     +(w(i,4,k)    - w(i,3,k)   )
          wy(i,2,k)    = ca*(w(i,3,k)    - w(i,1,k)   ) &
     &                  -cb*(w(i,4,k)    - w(i,1,k)   ) &
     &                     +(w(i,5,k)    - w(i,2,k)   )
          wy(i,3,k)    = ca*(w(i,4,k)    - w(i,2,k)   ) &
     &                  -cb*(w(i,5,k)    - w(i,1,k)   ) &
     &                     +(w(i,6,k)    - w(i,1,k)   )
          wy(i,ny-2,k) = ca*(w(i,ny-1,k) - w(i,ny-3,k)) &
     &                  -cb*(w(i,ny,k)   - w(i,ny-4,k)) &
     &                     +(            - w(i,ny-5,k)) 
          wy(i,ny-1,k) = ca*(w(i,ny,k)   - w(i,ny-2,k)) &
     &                  -cb*(            - w(i,ny-3,k)) &
     &                     +(            - w(i,ny-4,k))
          wy(i,ny,k)   = ca*(            - w(i,ny-1,k)) &
     &                  -cb*(            - w(i,ny-2,k)) &
     &                     +(            - w(i,ny-3,k))
        enddo
        enddo
      else
        do k=1,nz
        do i=1,nx
          wy(i,1,k)    = ca*(w(i,2,k)    + w(i,1,k)   ) &
     &                  -cb*(w(i,3,k)    + w(i,2,k)   ) &
     &                     +(w(i,4,k)    + w(i,3,k)   )
          wy(i,2,k)    = ca*(w(i,3,k)    - w(i,1,k)   ) &
     &                  -cb*(w(i,4,k)    + w(i,1,k)   ) &
     &                     +(w(i,5,k)    + w(i,2,k)   )
          wy(i,3,k)    = ca*(w(i,4,k)    - w(i,2,k)   ) &
     &                  -cb*(w(i,5,k)    - w(i,1,k)   ) &
     &                     +(w(i,6,k)    + w(i,1,k)   )
          wy(i,ny-2,k) = ca*(w(i,ny-1,k) - w(i,ny-3,k)) &
     &                  -cb*(w(i,ny,k)   - w(i,ny-4,k)) &
     &                     +(            - w(i,ny-5,k)) 
          wy(i,ny-1,k) = ca*(w(i,ny,k)   - w(i,ny-2,k)) &
     &                  -cb*(            - w(i,ny-3,k)) &
     &                     +(            - w(i,ny-4,k)) 
          wy(i,ny,k)   = ca*(            - w(i,ny-1,k)) &
     &                  -cb*(            - w(i,ny-2,k)) &
     &                     +(            - w(i,ny-3,k))
        enddo
        enddo
      endif

      do k=1,nz
      do j=4,ny-3
      do i=1,nx
        wy(i,j,k)    = ca*(w(i,j+1,k) - w(i,j-1,k) ) &
     &                -cb*(w(i,j+2,k) - w(i,j-2,k) ) &
     &                   +(w(i,j+3,k) - w(i,j-3,k) )
      enddo
      enddo
      enddo

      ccdx= 1/(cc*dx)
      do i=1,nx*ny*nz
        wy(i,1,1) = wy(i,1,1) *ccdx
      enddo

      return
  end subroutine dery

  function Opt_Z_EV8(Grid,Parity,Signature, TimeSimplex,Component) result(Der)
    
    integer,intent(in) :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    integer                           :: i,j,k, S

    allocate(Der(nx,ny,nz)) ; Der = 0.0_dp
    select case (Component)
    case(1)
        S =  Signature*Parity
    case(2)
        S =  Signature*Parity
    case(3)
        S = -Signature*Parity
    case(4)
        S = -Signature*Parity
    end select
    call derz(Grid, der, S)

  end function Opt_Z_EV8

  subroutine derz (w,wz,iz)

      integer, intent(in):: iz
      integer            :: i,j,k
      real(KIND=dp), intent(in) ::w(nx,ny,nz)
      real(KIND=dp), intent(out)::wz(nx,ny,nz)
      real(KIND=dp)             :: ccdx
      real(KIND=dp),parameter :: ca=45.0d0,cb=9.0d0,cc=60.0d0


      if (iz.gt.0) then
        do i=1,nx*ny
          wz(i,1,1)    = ca*(w(i,1,2)    - w(i,1,1)   ) &
     &                  -cb*(w(i,1,3)    - w(i,1,2)   ) &
     &                     +(w(i,1,4)    - w(i,1,3)   )
          wz(i,1,2)    = ca*(w(i,1,3)    - w(i,1,1)   ) &
     &                  -cb*(w(i,1,4)    - w(i,1,1)   ) &
     &                     +(w(i,1,5)    - w(i,1,2)   ) 
          wz(i,1,3)    = ca*(w(i,1,4)    - w(i,1,2)   ) &
     &                  -cb*(w(i,1,5)    - w(i,1,1)   ) &
     &                     +(w(i,1,6)    - w(i,1,1)   )
          wz(i,1,nz-2) = ca*(w(i,1,nz-1) - w(i,1,nz-3)) &
     &                  -cb*(w(i,1,nz)   - w(i,1,nz-4)) &
     &                     +(            - w(i,1,nz-5))
          wz(i,1,nz-1) = ca*(w(i,1,nz)   - w(i,1,nz-2)) &
     &                  -cb*(            - w(i,1,nz-3)) &
     &                     +(            - w(i,1,nz-4))
          wz(i,1,nz)   = ca*(            - w(i,1,nz-1)) &
     &                  -cb*(            - w(i,1,nz-2)) &
     &                     +(            - w(i,1,nz-3))
        enddo
      else
        do i=1,nx*ny
          wz(i,1,1)    = ca*(w(i,1,2)    + w(i,1,1)   ) &
     &                  -cb*(w(i,1,3)    + w(i,1,2)   ) &
     &                     +(w(i,1,4)    + w(i,1,3)   )
          wz(i,1,2)    = ca*(w(i,1,3)    - w(i,1,1)   ) &
     &                  -cb*(w(i,1,4)    + w(i,1,1)   ) &
     &                     +(w(i,1,5)    + w(i,1,2)   )
          wz(i,1,3)    = ca*(w(i,1,4)    - w(i,1,2)   ) &
     &                  -cb*(w(i,1,5)    - w(i,1,1)   ) &
     &                     +(w(i,1,6)    + w(i,1,1)   )
          wz(i,1,nz-2) = ca*(w(i,1,nz-1) - w(i,1,nz-3)) &
     &                  -cb*(w(i,1,nz)   - w(i,1,nz-4)) &
     &                     +(            - w(i,1,nz-5))
          wz(i,1,nz-1) = ca*(w(i,1,nz)   - w(i,1,nz-2)) &
     &                  -cb*(            - w(i,1,nz-3)) &
     &                     +(            - w(i,1,nz-4))
          wz(i,1,nz)   = ca*(            - w(i,1,nz-1)) &
     &                  -cb*(            - w(i,1,nz-2)) &
     &                     +(            - w(i,1,nz-3))
        enddo
      endif


      do k=4,nz-3
      do i=1,nx*ny
        wz(i,1,k)    = ca*(w(i,1,k+1) - w(i,1,k-1) ) &
     &                -cb*(w(i,1,k+2) - w(i,1,k-2) ) &
     &                   +(w(i,1,k+3) - w(i,1,k-3) )
      enddo
      enddo

      ccdx=1.0_dp/ (cc*dx)
      do i=1,nx*ny*nz
        wz(i,1,1) = wz(i,1,1) * ccdx
      enddo
      return
  end subroutine derz

  function Opt_Z_EV4(Grid,Parity,Signature, TimeSimplex,Component) result(Der)
    
    integer,intent(in) :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp),allocatable         :: Der(:,:,:)
    integer                           :: i,j,k, S
    allocate(Der(nx,ny,nz)) ; Der = 0.0_dp
    call derz_ev4(Grid, der)

  end function Opt_Z_EV4

  subroutine derz_ev4 (w,wz)

      integer            :: i,j,k
      real(KIND=dp), intent(in) ::w(nx,ny,nz)
      real(KIND=dp), intent(out)::wz(nx,ny,nz)
      real(KIND=dp)                     :: ccdx
      real(KIND=dp),parameter :: ca=45.0d0,cb=9.0d0,cc=60.0d0

      do i=1,nx*ny
        wz(i,1,1)    = ca*(w(i,1,2)               ) & 
     &                -cb*(w(i,1,3)               ) &
     &                   +(w(i,1,4)               )
        wz(i,1,2)    = ca*(w(i,1,3)   -w(i,1,1)   ) &
     &                -cb*(w(i,1,4)               ) &
     &                   +(w(i,1,5)               )
        wz(i,1,3)    = ca*(w(i,1,4)   -w(i,1,2)   ) &
     &                -cb*(w(i,1,5)   -w(i,1,1)   ) &
     &                   +(w(i,1,6)               )
        wz(i,1,nz-2) = ca*(w(i,1,nz-1)-w(i,1,nz-3)) &
     &                -cb*(w(i,1,nz)  -w(i,1,nz-4)) &
     &                   +(           -w(i,1,nz-5))
        wz(i,1,nz-1) = ca*(w(i,1,nz)  -w(i,1,nz-2)) &
     &                -cb*(           -w(i,1,nz-3)) &
     &                   +(           -w(i,1,nz-4))
        wz(i,1,nz)   = ca*(           -w(i,1,nz-1)) &
     &                -cb*(           -w(i,1,nz-2)) &
     &                   +(           -w(i,1,nz-3))
      enddo

      do k=4,nz-3
      do i=1,nx*ny
        wz(i,1,k)    = ca*(w(i,1,k+1) - w(i,1,k-1) ) &
     &                -cb*(w(i,1,k+2) - w(i,1,k-2) ) &
     &                   +(w(i,1,k+3) - w(i,1,k-3) )
      enddo
      enddo

      ccdx = 1.0_dp/(cc*dx)
      do i=1,nx*ny*nz
        wz(i,1,1) = wz(i,1,1) *ccdx
      enddo
      return
  end subroutine derz_EV4

  function Lapla_EV8(Grid,Parity, Signature, TimeSimplex,Component) result(Lap)

    integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
    real(KIND=dp), target, intent(in) :: Grid(:,:,:)
    real(KIND=dp), allocatable        :: Lap(:,:,:)
    real(KIND=dp)                     :: DerX(nx,ny,nz)
    real(KIND=dp)                     :: DerZ(nx,ny,nz)
    real(KIND=dp)                     :: DerY(nx,ny,nz)
    integer                           :: xp,yp,zp

    allocate(Lap(nx,ny,nz)) ; Lap = 0.0_dp
    select case(Component)
    case(1)
        xp = Signature ; yp = TimeSimplex ; zp = Signature*Parity
    case(2)
        xp =-Signature ; yp =-TimeSimplex ; zp = Signature*Parity
    case(3)
        xp =-Signature ; yp = TimeSimplex ; zp =-Signature*Parity
    case(4)
        xp = Signature ; yp =-TimeSimplex ; zp =-Signature*Parity
    end select 

    call lapla_X_ev8(Grid,DerX,xp)
    call lapla_Y_ev8(Grid,DerY,yp)
    call lapla_Z_ev8(Grid,DerZ,zp)

    Lap =  Derx + Dery + Derz
  end function Lapla_EV8
  
  subroutine lapla_X_ev8 (w,dw,xp)

      real(KIND=dp), parameter :: xxf=-8064.0d0,xxm=43050.0d0,xx2=1008.0d0 
      real(KIND=dp), parameter :: xx3=-128.0d0,xx4=9.0d0,xxd=5040.0d0
      integer, intent(in)      :: xp
      real(KIND=dp), intent(in):: w(nx,ny,nz)
      real(KIND=dp)            :: dw(nx,ny,nz)
      real(KIND=dp)            :: xf, xm, x2,x3,x4, xd
      integer                  :: i,j,k

      xf = xxf                  ! = -8064.0d0
      xm = xxm / (xf*3)         ! = 43050.0d0 / xf
      x2 = xx2 / xf             ! =  1008.0d0 / xf
      x3 = xx3 / xf             ! =  -128.0d0 / xf
      x4 = xx4 / xf             ! =     9.0d0 / xf
      xd = xxd / xf * (dx*dx)   ! =  5040.0d0 / xf * (dx*dx)
      xd = 1.0/xd

      do i=1,nx*ny*nz
        dw(i,1,1) = xm * w(i,1,1)
      enddo

      do j=1,ny*nz
        dw(1,j,1)    = dw(1,j,1)   +   w(1,j,1)*xp+x2*w(2,j,1)*xp &
     &                             +x3*w(3,j,1)*xp+x4*w(4,j,1)*xp
        dw(2,j,1)    = dw(2,j,1)   +   w(1,j,1)   +x2*w(1,j,1)*xp &
     &                             +x3*w(2,j,1)*xp+x4*w(3,j,1)*xp
        dw(3,j,1)    = dw(3,j,1)   +   w(2,j,1)   +x2*w(1,j,1)    &
     &                             +x3*w(1,j,1)*xp+x4*w(2,j,1)*xp 
        dw(4,j,1)    = dw(4,j,1)   +   w(3,j,1)   +x2*w(2,j,1)    &
     &                             +x3*w(1,j,1)   +x4*w(1,j,1)*xp
        dw(nx-1,j,1) = dw(nx-1,j,1)+   w(nx,j,1)
        dw(nx-2,j,1) = dw(nx-2,j,1)+   w(nx-1,j,1)+x2*w(nx,j,1)
        dw(nx-3,j,1) = dw(nx-3,j,1)+   w(nx-2,j,1)+x2*w(nx-1,j,1) &
     &                             +x3*w(nx,j,1)
      enddo
      do j=1,ny*nz
      do i=1,nx-4
         dw(i,j,1) = dw(i,j,1)   +   w(i+1,j,1) +x2*w(i+2,j,1) &
     &                           +x3*w(i+3,j,1) +x4*w(i+4,j,1)
      enddo
      enddo
      do j=1,ny*nz
      do i=5,nx
         dw(i,j,1) = dw(i,j,1)   +   w(i-1,j,1) +x2*w(i-2,j,1) &
     &                           +x3*w(i-3,j,1) +x4*w(i-4,j,1)
      enddo
      enddo

      do i=1,nx*ny*nz
        dw(i,1,1) = - dw(i,1,1) * xd
      enddo

      return
  end subroutine lapla_X_ev8

  subroutine lapla_Y_ev8 (w,dw,yp)

      real(KIND=dp), parameter :: xxf=-8064.0d0,xxm=43050.0d0,xx2=1008.0d0 
      real(KIND=dp), parameter :: xx3=-128.0d0,xx4=9.0d0,xxd=5040.0d0
      integer, intent(in)      :: yp
      real(KIND=dp), intent(in):: w(nx,ny,nz)
      real(KIND=dp)            :: dw(nx,ny,nz), w_t(ny,nx,nz), d_T(ny,nx,nz)
      real(KIND=dp)            :: xf, xm, x2,x3,x4, xd
      integer                  :: i,j,k
      
      xf = xxf                  ! = -8064.0d0
      xm = xxm / (xf*3)         ! = 43050.0d0 / xf
      x2 = xx2 / xf             ! =  1008.0d0 / xf
      x3 = xx3 / xf             ! =  -128.0d0 / xf
      x4 = xx4 / xf             ! =     9.0d0 / xf
      xd = xxd / xf * (dx*dx)   ! =  5040.0d0 / xf * (dx*dx)
      xd = 1.0_dp/xd

      do i=1,nx*ny*nz
        dw(i,1,1) = xm * w(i,1,1)
      enddo

      do i=1,nx
      do j=1,ny
      do k=1,nz
        w_t(j,i,k) =  w(i,j,k)
        d_t(j,i,k) = dw(i,j,k)
      enddo
      enddo
      enddo
      do k=1,nx*nz
        d_t(1,k,1)    = d_t(1,k,1)   +   w_t(1,k,1)*yp+x2*w_t(2,k,1)*yp &
     &                               +x3*w_t(3,k,1)*yp+x4*w_t(4,k,1)*yp
        d_t(2,k,1)    = d_t(2,k,1)   +   w_t(1,k,1)   +x2*w_t(1,k,1)*yp &
     &                               +x3*w_t(2,k,1)*yp+x4*w_t(3,k,1)*yp
        d_t(3,k,1)    = d_t(3,k,1)   +   w_t(2,k,1)   +x2*w_t(1,k,1)    &
     &                               +x3*w_t(1,k,1)*yp+x4*w_t(2,k,1)*yp
        d_t(4,k,1)    = d_t(4,k,1)   +   w_t(3,k,1)   +x2*w_t(2,k,1)    &
     &                               +x3*w_t(1,k,1)   +x4*w_t(1,k,1)*yp
        d_t(ny-1,k,1) = d_t(ny-1,k,1)+   w_t(ny,k,1)
        d_t(ny-2,k,1) = d_t(ny-2,k,1)+   w_t(ny-1,k,1)+x2*w_t(ny,k,1)
        d_t(ny-3,k,1) = d_t(ny-3,k,1)+   w_t(ny-2,k,1)+x2*w_t(ny-1,k,1) &
     &                               +x3*w_t(ny,k,1)
      enddo
      do k=1,nx*nz
      do j=1,ny-4
        d_t(j,k,1) = d_t(j,k,1)   +   w_t(j+1,k,1) +x2*w_t(j+2,k,1) &
     &                            +x3*w_t(j+3,k,1) +x4*w_t(j+4,k,1)
      enddo
      enddo
      do k=1,nx*nz
      do j=5,ny
        d_t(j,k,1) = d_t(j,k,1)   +   w_t(j-1,k,1) +x2*w_t(j-2,k,1) &
     &                            +x3*w_t(j-3,k,1) +x4*w_t(j-4,k,1)
      enddo
      enddo
      do i=1,nx
      do j=1,ny
      do k=1,nz
        dw(i,j,k) = d_t(j,i,k)
      enddo
      enddo
      enddo

      do i=1,nx*ny*nz
        dw(i,1,1) = - dw(i,1,1) * xd
      enddo

    end subroutine lapla_Y_ev8

    subroutine lapla_Z_ev8 (w,dw,zp)

        real(KIND=dp), parameter :: xxf=-8064.0d0,xxm=43050.0d0,xx2=1008.0d0 
        real(KIND=dp), parameter :: xx3=-128.0d0,xx4=9.0d0,xxd=5040.0d0
        integer, intent(in)      :: zp
        real(KIND=dp), intent(in):: w(nx,ny,nz)
        real(KIND=dp)            :: dw(nx,ny,nz)
        real(KIND=dp)            :: xf, xm, x2,x3,x4, xd
        integer                  :: i,j,k

        xf = xxf                  ! = -8064.0d0
        xm = xxm / (xf*3)         ! = 43050.0d0 / xf
        x2 = xx2 / xf             ! =  1008.0d0 / xf
        x3 = xx3 / xf             ! =  -128.0d0 / xf
        x4 = xx4 / xf             ! =     9.0d0 / xf
        xd = xxd / xf * (dx*dx)   ! =  5040.0d0 / xf * (dx*dx)
        xd = 1.0_dp/xd

        do i=1,nx*ny*nz
            dw(i,1,1) = xm * w(i,1,1)
        enddo

        do i=1,nx*ny 
            dw(i,1,1)    = dw(i,1,1)   +   w(i,1,1)*zp+x2*w(i,1,2)*zp     &
            &                             +x3*w(i,1,3)*zp+x4*w(i,1,4)*zp
            dw(i,1,2)    = dw(i,1,2)   +   w(i,1,1)   +x2*w(i,1,1)*zp     &
            &                             +x3*w(i,1,2)*zp+x4*w(i,1,3)*zp  
            dw(i,1,3)    = dw(i,1,3)   +   w(i,1,2)   +x2*w(i,1,1)        &
            &                             +x3*w(i,1,1)*zp+x4*w(i,1,2)*zp
            dw(i,1,4)    = dw(i,1,4)   +   w(i,1,3)   +x2*w(i,1,2)        &
            &                          +x3*w(i,1,1)   +x4*w(i,1,1)*zp
            dw(i,1,nz-1) = dw(i,1,nz-1)+   w(i,1,nz)                      
            dw(i,1,nz-2) = dw(i,1,nz-2)+   w(i,1,nz-1)+x2*w(i,1,nz)       
            dw(i,1,nz-3) = dw(i,1,nz-3)+   w(i,1,nz-2)+x2*w(i,1,nz-1)     &
            &                          +x3*w(i,1,nz)
        enddo
        do k=1,nz-4
            do i=1,nx*ny
                dw(i,1,k) = dw(i,1,k)   +   w(i,1,k+1) +x2*w(i,1,k+2)    &
                &                          +x3*w(i,1,k+3) +x4*w(i,1,k+4)
            enddo
        enddo
        do k=5,nz
            do i=1,nx*ny
                dw(i,1,k) = dw(i,1,k)   +   w(i,1,k-1) +x2*w(i,1,k-2)    &
                &                          +x3*w(i,1,k-3) +x4*w(i,1,k-4)
            enddo
        enddo

        do i=1,nx*ny*nz
            dw(i,1,1) = - dw(i,1,1) * xd
        enddo
    end subroutine lapla_Z_ev8

    function Lapla_EV4(Grid,Parity, Signature, TimeSimplex,Component) result(Lap)

        integer,intent(in)                :: Parity,Signature,TimeSimplex,Component
        real(KIND=dp), target, intent(in) :: Grid(:,:,:)
        real(KIND=dp), allocatable        :: Lap(:,:,:)
        real(KIND=dp)                     :: DerX(nx,ny,nz)
        real(KIND=dp)                     :: DerZ(nx,ny,nz)
        real(KIND=dp)                     :: DerY(nx,ny,nz)
        integer                           :: xp,yp,zp

        allocate(Lap(nx,ny,nz)) ; Lap = 0.0_dp
        select case(Component)
        case(1)
            xp = Signature ; yp = TimeSimplex ; zp = Signature*Parity
        case(2)
            xp =-Signature ; yp =-TimeSimplex ; zp = Signature*Parity
        case(3)
            xp =-Signature ; yp = TimeSimplex ; zp =-Signature*Parity
        case(4)
            xp = Signature ; yp =-TimeSimplex ; zp =-Signature*Parity
        end select 

        call lapla_X_ev8(Grid,DerX,xp)
        call lapla_Y_ev8(Grid,DerY,yp)
        call lapla_Z_ev4(Grid,DerZ)

        Lap =  Derx + Dery + Derz
      end function Lapla_EV4

      subroutine lapla_Z_EV4 (w,dw)

        real(KIND=dp), parameter :: xxf=-8064.0d0,xxm=43050.0d0,xx2=1008.0d0 
        real(KIND=dp), parameter :: xx3=-128.0d0,xx4=9.0d0,xxd=5040.0d0
        real(KIND=dp), intent(in):: w(nx,ny,nz)
        real(KIND=dp)            :: dw(nx,ny,nz)
        real(KIND=dp)            :: xf, xm, x2,x3,x4, xd
        integer                  :: i,j,k

        xf = xxf                  ! = -8064.0d0
        xm = xxm / xf             ! = 43050.0d0 / xf
        x2 = xx2 / xf             ! =  1008.0d0 / xf
        x3 = xx3 / xf             ! =  -128.0d0 / xf
        x4 = xx4 / xf             ! =     9.0d0 / xf
        xd = xxd / xf * (dx*dx)   ! =  5040.0d0 / xf * (dx*dx)

        do i=1,nx*ny*nz
            dw(i,1,1) = xm * w(i,1,1)/3
        enddo

        do i=1,nx*ny
            dw(i,1,2)    = dw(i,1,2)   +   w(i,1,1)
            dw(i,1,3)    = dw(i,1,3)   +   w(i,1,2)   +x2*w(i,1,1)
            dw(i,1,4)    = dw(i,1,4)   +   w(i,1,3)   +x2*w(i,1,2)+x3*w(i,1,1)
            dw(i,1,nz-1) = dw(i,1,nz-1)+   w(i,1,nz)
            dw(i,1,nz-2) = dw(i,1,nz-2)+   w(i,1,nz-1)+x2*w(i,1,nz)
            dw(i,1,nz-3) = dw(i,1,nz-3)+   w(i,1,nz-2)+x2*w(i,1,nz-1)+x3*w(i,1,nz)
        enddo
        do k=1,nz-4
            do i=1,nx*ny
                dw(i,1,k) = dw(i,1,k)   +   w(i,1,k+1) +x2*w(i,1,k+2)   &
                &                          +x3*w(i,1,k+3) +x4*w(i,1,k+4)
            enddo
        enddo
        do k=5,nz
            do i=1,nx*ny
                dw(i,1,k) = dw(i,1,k)   +   w(i,1,k-1) +x2*w(i,1,k-2)   &
                &                          +x3*w(i,1,k-3) +x4*w(i,1,k-4)
            enddo
        enddo

        do i=1,nx*ny*nz
            dw(i,1,1) = - dw(i,1,1) / xd
        enddo
    end subroutine lapla_Z_ev4

end module OptimizedDerivatives
