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

        allocate(Der(nx,ny,nz)); Der = 0.0_dp
        
        ! The reflection sign around x = 0, in an EV8 setting is
        !     S
        ! 1   + 
        ! 2   - 
        ! 3   - 
        ! 4   +
        select case(Component)
        case(1)
            S =  1
        case(2)
            S = -1
        case(3)
            S = -1
        case(4)
            S =  1
        end select
        
        do k=1,nz
            do j=1,ny
                !Start of the line
                Der(1,j,k) =     45.0d0 * (Grid(2,j,k) - S*Grid(1,j,k)) &
                &          -      9.0d0 * (Grid(3,j,k) - S*Grid(2,j,k)) &
                &          +              (Grid(4,j,k) - S*Grid(3,j,k))
                Der(2,j,k) =     45.0d0 * (Grid(3,j,k) -   Grid(1,j,k)) &
                &          -      9.0d0 * (Grid(4,j,k) - S*Grid(1,j,k)) &
                &          +              (Grid(5,j,k) - S*Grid(2,j,k))
                Der(3,j,k) =     45.0d0 * (Grid(4,j,k) -   Grid(2,j,k)) &
                &          -      9.0d0 * (Grid(5,j,k) -   Grid(1,j,k)) &
                &          +              (Grid(6,j,k) - S*Grid(1,j,k)) 

                ! End of the line
                Der(nx-2,j,k) =  45.0d0 * (Grid(nx-1,j,k) - Grid(nx-3,j,k)) &
                &             -   9.0d0 * (Grid(nx,j,k)   - Grid(nx-4,j,k)) &
                &             +           (               - Grid(nx-5,j,k))
                Der(nx-1,j,k) =  45.0d0 * (Grid(nx  ,j,k) - Grid(nx-2,j,k)) &
                &             -   9.0d0 * (               - Grid(nx-3,j,k)) &
                &             +           (               - Grid(nx-4,j,k))
                Der(nx  ,j,k) =  45.0d0 * (               - Grid(nx-1,j,k)) &
                &             -   9.0d0 * (               - Grid(nx-2,j,k)) &
                &             +           (               - Grid(nx-3,j,k))
            enddo
        enddo

        do k=1,nz
            do j=1,ny
                do i=4,nx-3
                    Der(i,j,k) =  45.0d0 * (Grid(i+1,j,k) - Grid(i-1,j,k)) &
                    &          -  9.0d0  * (Grid(i+2,j,k) - Grid(i-2,j,k)) &
                    &          +           (Grid(i+3,j,k) - Grid(i-3,j,k))
                enddo
            enddo
        enddo

        Der = Der/(60.0 * dx)

    end function Opt_X_EV8

end module OptimizedDerivatives