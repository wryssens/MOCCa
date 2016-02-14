module Mesh
!-------------------------------------------------------------------------------
! Note
!       For different combinations of symmetries we can make several choices on 
!       how to store the mesh, i.e. which points to store.
!       For example when breaking time simplex and conserving parity and 
!       signature (and thus breaking Time reversal) we can either
!             - Store the entire x-axis, storing half of both y- and z-axis; Or
!             - Store the entire y-axis, storing half of both x- and z-axis
!       In MOCCa, we make the following choice:
!             * Everything conserved : postive x/y/z-axis
!             * Breaking Signature: Entire x-axis
!             * Breaking Parity: Entire z-axis
!             * Breaking Time Simplex: Entire y-axis
!-------------------------------------------------------------------------------

  use CompilationInfo
  use GenInfo

  public
  
  save

  real(KIND=dp), public, allocatable :: MeshX(:), MeshY(:), MeshZ(:)
  ! These come in handy as representations of the \hat{x}_i operators
  real(KIND=dp), public, allocatable :: Mesh3DX(:,:,:)
  real(KIND=dp), public, allocatable :: Mesh3DY(:,:,:), Mesh3DZ(:,:,:)
  real(KIND=dp), public, allocatable :: Mesh3D(:,:,:,:) 

contains

  subroutine IniMesh
    !---------------------------------------------------------------------------
    ! This routine allocates the arrays and initialises the coordinates
    ! of the mesh.
    !---------------------------------------------------------------------------
    integer :: i,j,k

    if(allocated(MeshX)) then
        deallocate (MeshX,MeshY,MeshZ,Mesh3DX, Mesh3DY,Mesh3DZ, Mesh3D)
    endif

    !Allocating everything
    allocate(MeshX(nx)); allocate(MeshY(ny)); allocate(MeshZ(nz))
    allocate(Mesh3DX(nx,ny,nz)); allocate(Mesh3DY(ny,ny,nz)); 
    allocate(Mesh3DZ(nz,ny,nz))
    allocate(Mesh3D(3,nx,ny,nz))

    do i=1,nx
        MeshX(i) = dx/2.0_dp + dx*(i-1)
    enddo
    do j=1,ny
        MeshY(j) = dx/2.0_dp + dx*(j-1)
    enddo
    do k=1,nz
        MeshZ(k) = dx/2.0_dp + dx*(k-1)
    enddo

    !Shifting the mesh in function of symmetries.
    if(.not.TSC) then
        ! In our way 
        MeshY = MeshY - ny/2.0_dp*dx
    endif
    
    if((.not.SC)) then
         !If signature is not conserved, we need the full x-axis.
         MeshX = MeshX - nx/2.0_dp*dx
    endif
    
    if(.not.PC) then
        !If parity is not conserved, we need the full z-axis.
        MeshZ = MeshZ - nz/2.0_dp*dx
    endif
    
    ! Copying everything to the 3D Meshes
    do i=1,nx
        Mesh3DX(i,:,:) = MeshX(i)
        Mesh3D(1,i,:,:)= MeshX(i)
    enddo
    do j=1,ny
        Mesh3DY(:,j,:) = MeshY(j)
        Mesh3D(2,:,j,:)= MeshY(j)
    enddo    
    do k=1,nz
        Mesh3DZ(:,:,k) = MeshZ(k)
        Mesh3D(3,:,:,k)= MeshZ(k)
    enddo            

    return
  end subroutine IniMesh
        
 function legendre_pnm ( n, m, x ) result(cx)
      !-------------------------------------------------------------------------
      ! Taken on 14/05/13 from
      !    http://people.sc.fsu.edu/~jburkardt/f_src/polpak/polpak.html
      ! Under the GNU LGPL licence, which can be found at
      !    http://www.gnu.org/licenses/gpl.html
      !
      ! Minor changes from the original version:
      !    - double precision instead of single precision
      !    - changed from subroutine to function, for clarity and ease of use
      !    - changed to pure function
      !*************************************************************************
      !
      ! LEGENDRE_PNM evaluates the associated Legendre polynomial Pnm(x) of the 
      !first kind.
      !
      !
      !  Differential equation:
      !
      !    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
      !
      !  First terms:
      !
      !    M = 0  ( = Legendre polynomials of first kind P(N)(X) )
      !
      !    P00 =    1
      !    P10 =    1 X
      !    P20 = (  3 X**2 -   1)/2
      !    P30 = (  5 X**3 -   3 X)/2
      !    P40 = ( 35 X**4 -  30 X**2 +   3)/8
      !    P50 = ( 63 X**5 -  70 X**3 +  15 X)/8
      !    P60 = (231 X**6 - 315 X**4 + 105 X**2 -  5)/16
      !    P70 = (429 X**7 - 693 X**5 + 315 X**3 - 35 X)/16
      !
      !    M = 1
      !
      !    P01 =   0
      !    P11 =   1 * SQRT(1-X*X)
      !    P21 =   3 * SQRT(1-X*X) * X
      !    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
      !    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
      !
      !    M = 2
      !
      !    P02 =   0
      !    P12 =   0
      !    P22 =   3 * (1-X*X)
      !    P32 =  15 * (1-X*X) * X
      !    P42 = 7.5 * (1-X*X) * (7*X*X-1)
      !
      !    M = 3
      !
      !    P03 =   0
      !    P13 =   0
      !    P23 =   0
      !    P33 =  15 * (1-X*X)**1.5
      !    P43 = 105 * (1-X*X)**1.5 * X
      !
      !    M = 4
      !
      !    P04 =   0
      !    P14 =   0
      !    P24 =   0
      !    P34 =   0
      !    P44 = 105 * (1-X*X)**2
      !
      !  Recursion:
      !
      !    if N < M:
      !      P(N,M) = 0
      !    if N = M:
      !      P(N,M) = (2*M-1)!! * (1-X*X)**(M/2) where N!! means the product of
      !      all the odd integers less than or equal to N.
      !    if N = M+1:
      !      P(N,M) = X*(2*M+1)*P(M,M)
      !    if N > M+1:
      !      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
      !
      !  Restrictions:
      !
      !    -1 <= X <= 1
      !     0 <= M <= N
      !
      !  Special values:
      !
      !    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
      !    polynomial of the first kind equals the Legendre polynomial of the
      !    first kind.
      !
      !  Modified:
      !
      !    24 June 2001
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer N, the maximum first index of the Legendre
      !    polynomial, which must be at least 0.
      !
      !    Input, integer M, the second index of the Legendre polynomial,
      !    which must be at least 0, and no greater than N.
      !
      !    Input, real X, the point at which the polynomial is to be
      !    evaluated.  X must satisfy -1 <= X <= 1.
      !
      !    Output, real CX(0:N), the values of the first N+1 polynomials.
      !-------------------------------------------------------------------------
      implicit none

      real(KIND=dp)    :: cx(0:n)
      real(KIND=dp)     :: fact, somx2
      real(KIND=dp),  intent(in) :: x
      integer, intent(in) :: m,n
      integer :: i

      if ( m < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEGENDRE_PNM - Fatal error!'
      write ( *, '(a,i6)' ) '  Input value of M is ', m
      write ( *, '(a)' ) '  but M must be nonnegative.'
      stop
      end if

      if ( m > n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEGENDRE_PNM - Fatal error!'
      write ( *, '(a,i6)' ) '  Input value of M = ', m
      write ( *, '(a,i6)' ) '  Input value of N = ', n
      write ( *, '(a)' ) '  but M must be less than or equal to N.'
      stop
      end if

      if ( x < -1.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEGENDRE_PNM - Fatal error!'
      write ( *, '(a,g14.6)' ) '  Input value of X = ', x
      write ( *, '(a)' ) '  but X must be no less than -1.'
      stop
      end if

      if ( x > 1.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEGENDRE_PNM - Fatal error!'
      write ( *, '(a,g14.6)' ) '  Input value of X = ', x
      write ( *, '(a)' ) '  but X must be no more than 1.'
      stop
      end if

      cx(0:m-1) = 0.0E+00

      cx(m) = 1.0E+00
      somx2 = sqrt ( 1.0E+00 - x**2 )

      fact = 1.0E+00
      do i = 1, m
          cx(m) = - cx(m) * fact * somx2
          fact = fact + 2.0E+00
      end do

      if ( m == n ) then
          return
      end if

      cx(m+1) = x * real ( 2 * m + 1 ) * cx(m)

      do i = m+2, n
          cx(i) = ( real ( 2 * i - 1 ) * x * cx(i-1) &
          - real ( i + m - 1 ) * cx(i-2) ) / real ( i - m )
      end do

      return
  end function Legendre_pnm

end module Mesh
