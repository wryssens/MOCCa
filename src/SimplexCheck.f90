module SimplexCheck
!-------------------------------------------------------------------------------
! Transform all of the spwfs with an application of simplex.
! 
!
!
!
!-------------------------------------------------------------------------------

use spwfstorage
use spinors
use wavefunctions
use geninfo
use densities

implicit none

contains


    subroutine ApplySimplex
    
        integer :: wave,i,j,k
        type(Spinor) :: temp
        
!        do k=1,nz
!            do i=1,nx
!                Density%Rho(i,:,k,1) = Density%Rho(nx - i +1 ,:,nz - k +1 ,1)
!                Density%Rho(i,:,k,1) = Density%Rho(nx - i +1 ,:,nz - k +1 ,1)
!            enddo
!        enddo
!        ! Action of parity
        do wave=1,nwt
            HFBasis(wave)%Value = ActionOfParity(HFBasis(wave)%value, 0)
                     
            temp = NewSpinor() 
            do k=1,nz
                    temp%grid(1:nx,1:ny,k,1,1) = + HFBasis(wave)%value%grid(1:nx, 1:ny, nz - k + 1,3,1)
                    temp%grid(1:nx,1:ny,k,2,1) = - HFBasis(wave)%Value%grid(1:nx, 1:ny, nz - k + 1,4,1)
                    temp%grid(1:nx,1:ny,k,3,1) = + HFBasis(wave)%Value%grid(1:nx, 1:ny, nz - k + 1,1,1)
                    temp%grid(1:nx,1:ny,k,4,1) = - HFBasis(wave)%Value%grid(1:nx, 1:ny, nz - k + 1,2,1)
            enddo
            
            HFBasis(wave)%value = temp
        enddo
    end subroutine ApplySimplex

end module SimplexCheck
