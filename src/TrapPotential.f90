module TrapPotential
!-------------------------------------------------------------------------------
! Module for external potential that acts as a trap for nucleons.
! In the end, it's a linear constraint on r^2 without cutoff.
!-------------------------------------------------------------------------------

use Densities
use Mesh

implicit none

   real(KIND=dp) :: TrapConstant(2) = 0.0_dp
   real(KIND=dp) :: TrapHbarOmega = 0.0_dp
   real(KIND=dp),allocatable :: TrapUpot(:,:,:,:)

contains

  subroutine ReadTrapParameters
    !--------------------------------------------------------------------------
    ! Subroutine that reads the trap parameters
    ! and allocates the Trap potential
    !--------------------------------------------------------------------------

     integer :: i, j, k, it

     Namelist /TrapParameters/ TrapHbarOmega

     read(unit=*, NML=TrapParameters)

     print '(/,"************************************************")'
     print '(  "* Harmonic trap with hbar omega = ",1f12.6 ," *")',TrapHbarOmega
     print '(  "************************************************",/)'

     TrapConstant(:) = 0.5_dp * TrapHbarOmega * TrapHbarOmega / hbm(:)

     ! print '(" TrapConstant ",2f16.8)',TrapConstant

     !-------------------------------------------------------------------------
     ! construct trap potential
     !-------------------------------------------------------------------------
     allocate(TrapUPot(nx,ny,nz,2))

     do it=1,2
       do k=1,nz
       do j=1,ny
       do i=1,nx
          TrapUPot(i,j,k,it) = TrapConstant(it) &
                 & * ( MeshX(i)**2 + MeshY(j)**2 + MeshZ(k)**2 )
       enddo
       enddo
       enddo
     enddo
     
  end subroutine ReadTrapParameters

  subroutine PrintTrapEnergy(Etot)
    !-------------------------------------------------------------
    ! diagnostic printing 
    !-------------------------------------------------------------
    real(KIND=dp), intent(in) :: Etot
    real(KIND=dp) :: TrapEnergy
    real(KIND=dp) :: A

     print '(/,"***************************************************")'
     print '(  "* Harmonic trap with hbar omega = ",1f12.6 ,"    *")',TrapHbarOmega
     print '(  "***************************************************")'
     
     TrapEnergy = GetTrapEnergy()

     A = neutrons + protons 
     print '(" Energy without trap       ",1f24.12)',Etot-TrapEnergy
     print '(" Energy of trap            ",1f24.12)',TrapEnergy
     print '(" Etot                      ",1f24.12)',Etot
     print '(" Etot / A^4/3              ",1f24.12)',Etot/A**(4.0_dp/3.0_dp)
     print '(" Etot / (hbar omega A^4/3) ",1f24.12)',Etot/(TrapHbarOmega*A**(4.0_dp/3.0_dp))
 
  end subroutine PrintTrapEnergy


  function GetTrapUPot() result(TrapU)
    !-------------------------------------------------------------
    ! Return the contribution of the trap to the potential U(r)
    !-------------------------------------------------------------
    real(KIND=dp) :: TrapU(nx,ny,nz,2)

    TrapU(:,:,:,:) = TrapUPot(:,:,:,:)

  end function GetTrapUPot

  function GetTrapEnergy() result(TrapE)
    !-------------------------------------------------------------
    ! Return the contribution of the trap to the potential U(r)
    !-------------------------------------------------------------
    real(KIND=dp) :: TrapE
    integer       :: it
 
 !  TrapE = dv * sum(Density%Rho(:,:,:,1)*TrapUPot(:,:,:,1))
 !  print '("Summ. TrapE= ",1f14.8)',TrapE

    TrapE = dv * sum(Density%Rho(:,:,:,1)*TrapUPot(:,:,:,1)) &
      &    +dv * sum(Density%Rho(:,:,:,2)*TrapUPot(:,:,:,2))

  end function GetTrapEnergy


end module TrapPotential
