module PairingInteraction
!-------------------------------------------------------------------------------
! Module that contains everything related to the pairing interaction.
! Currently only a zero-range, density dependent interaction is offered.
! In addition, the pairing cutoff is taken to be 'part of' the pairing 
! interaction. In addition, this can calculate the effect of the pairing 
! interaction on a two-body state.
!-------------------------------------------------------------------------------

  use CompilationInfo
  use Geninfo
  use SpwfStorage
  use Densities
  
  implicit none

  public
  
  !-----------------------------------------------------------------------------
  ! Properties of the interaction.
  !-----------------------------------------------------------------------------
  ! Pairing strength for neutron & protons.
  real(KIND=dp) :: PairingStrength(2) = 0.0_dp
  ! Saturation Density. Put to 0.16 fm^{-3} for simplicity, for now.
  real(KIND=dp) :: RhoSat = 0.16_dp
  ! Density dependence of the pairing force.
  real(KIND=dp) :: alpha(2)=1.0_dp
  ! Determine the kind of cutoff used.
  ! CutType
  ! 1) Symmetric Fermi
  ! 2) Cosine
  integer       :: CutType=1

  ! Density-dependent factor for density-dependent pairing
  real(KIND=dp), allocatable :: DensityFactor(:,:,:,:)

  !-----------------------------------------------------------------------------
  !Again an abstract interface for the use of PGI compilers....
  abstract interface
    real(KIND=dp) function Cut(E, Lambda, it) result(Cutoff)
      import                    :: dp
      real(KIND=dp), intent(in) :: E, Lambda
      integer, intent(in)       :: it
    end function Cut
  end interface

  procedure(Cut), pointer :: PairingCutoff 
  !-----------------------------------------------------------------------------
  ! Cutoff for the pairing strength above and below the Fermi level.
  ! Traditionally set to 5Mev
  real(KIND=dp) :: PairingCut(2) = 5.0_dp, PairingMu(2) =0.5
  !-----------------------------------------------------------------------------
  ! Storage for all cutoffs. 
  real(KIND=dp), allocatable :: PCutoffs(:)
  
contains

  pure function GetPairDensity(wf1,wf2) result( ActionOfPairing)
    !---------------------------------------------------------------------------
    ! Calculate the action of the pairing interaction on a twobody state,
    ! represented by spinor 1 and 2.
    !
    ! So the result of this function is
    ! Sum_{s} s  < r , s ; r, -s | psi_1, psi_2 >
    ! 
    !---------------------------------------------------------------------------
    ! Note that this routine includes a TimeReversal operator if timereversal
    ! is conserved.
    !---------------------------------------------------------------------------
    ! Indices of the array ActionOfPairing.
    ! postion_x, position_y, position_z, real or imaginary
    !---------------------------------------------------------------------------

    type(Spinor), intent(in) :: wf1, wf2
    
    complex(KIND=dp)         :: ActionOfPairing(nx,ny,nz)
    real(KIND=dp)            :: Temp(2,nx,ny,nz)
    integer                  :: i

    !ActionOfPairing = 0.0_dp
    
    if(.not.TRC) then 
      !---------------------------------------------------------------------------
      ! Real part
      do i=1,nx*ny*nz
        Temp(1,i,1,1) =                                                      &
        &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,3,1)            &
        &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,4,1)            &
        &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,1,1)            &
        &             + wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,2,1)  
        !---------------------------------------------------------------------------
        ! Imaginary Part 
        Temp(2,i,1,1) =                                                      &
        &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,4,1)            &
        &             + wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,3,1)            &
        &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,2,1)            &
        &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,1,1)
      enddo
    else
      !---------------------------------------------------------------------------
      ! Real part
      do i=1,nx*ny*nz
        Temp(1,i,1,1) =                                                      &
        &             - wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,1,1)            &
        &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,2,1)            &
        &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,3,1)            &
        &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,4,1)  
        !---------------------------------------------------------------------------
        ! Imaginary Part 
        Temp(2,i,1,1) =                                                      &
        &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,2,1)            &
        &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,1,1)            &
        &             + wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,4,1)            &
        &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,3,1)
      enddo
    endif

    do i=1,nx*ny*nz
       ActionOfPairing(i,1,1) = dcmplx(Temp(1,i,1,1), Temp(2,i,1,1))
    enddo
    !if(all(ActionOfPairing .eq. 0.0 )) call stp('Action')
    return
  end function GetPairDensity
  
  subroutine ComputePairingCutoffs(Lambda)
  !-----------------------------------------------------------------------------
  ! Computes and stores all pairing cutoffs for further use.
  !
  !-----------------------------------------------------------------------------
    integer :: i, it
    real(KIND=dp), intent(in) :: Lambda(2)
    
    if(.not.allocated(PCutoffs)) then
        allocate(PCutoffs(nwt))
    endif
    
    do i=1,nwt
        it = (HFBasis(i)%GetIsospin() + 3)/2
        PCutoffs(i) = PairingCutoff(HFBasis(i)%GetEnergy(), Lambda(it), it)         
    enddo
  
  end subroutine ComputePairingCutoffs

  real(KIND=dp) function SymmetricFermi(E, Lambda, it) result(Cutoff)
    !---------------------------------------------------------------------------
    ! Calculates a cutoff that utilises two fermi functions, above and below the
    ! Fermi level.
    ! Cutoff is taken with a fermi function both above and below the fermi 
    !    energy.
    !        f^-2 = [1 + exp((epsilon - lambda - DeltaE)/mu)]
    !             * [1 + exp((- epsilon + lambda - DeltaE)/mu)]
    !    with mu and DeltaE being read from input.
    !    This is described in S.J. Krieger et al., Nucl. Phys.A517 (1990) 275
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: E, Lambda
    integer, intent(in)       :: it
    real(Kind=dp)             :: Up, Down
    
    Up   =   (E - Lambda - PairingCut(it))/PairingMu(it)
    Down = - 2*(E - Lambda + PairingCut(it))/PairingMu(it)
    Cutoff = sqrt(sqrt(1.0_dp/(1.0_dp + exp(Up))))
    Cutoff = Cutoff * sqrt(sqrt(1.0_dp/(1.0_dp + exp(Down))))

    return
  end function SymmetricFermi
  
  real(KIND=dp) function CosineCut(E, Lambda, it) result(Cutoff)
    !---------------------------------------------------------------------------
    ! Calculates a cutoff using a cosine function.
    ! 
    ! Define:
    ! xd = DeltaE - mu/2
    ! xu = DeltaE + mu/2
    !
    ! then
    ! f = 1                                     
    !   when abs(Epsilon - Lambda) <= xd
    ! f = 0.5 * cos ( (abs(Epsilon - lambda) -xd )*pi/2) + 0.5 
    !   when xd <= abs(Epsilon - Lambda) <= xu
    ! f = 0.0  
    !   when xu <= abs(Epsilon - Lambda)
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: E, Lambda
    integer, intent(in)       :: it
    real(Kind=dp)             :: Up, Down
    
    Down = PairingCut(it) - PairingMu(it)/2.0_dp
    Up   = PairingCut(it) + PairingMu(it)/2.0_dp
  
    if( abs(E - Lambda).le. Down) then
      Cutoff = 1.0_dp      
    elseif( abs(E-Lambda) .le. Up) then
      Cutoff = 0.5_dp * cos((abs(E-Lambda) - Down)*pi/PairingMu(it)) + 0.5
    else
      Cutoff = 0.0_dp
    endif
  
  end function Cosinecut

  subroutine CompDensityFactor
    !-----------------------------------------------------------------------------
    !Currently the only part of a general operator V that is implemented 
    ! is the Delta interaction with a density dependent part.
    ! It is described in
    !      S.J. Krieger et al., Nucl. Phys.A517 (1990) 275
    !
    ! This force is given by:
    !   v12 = 0.5 * strength * ( 1 - alpha * rho/rhosat) * (1-ps12)*delta(r12)
    !----------------------------------------------------------------------------
    integer :: i, it

    if(.not.allocated(DensityFactor)) allocate(DensityFactor(nx,ny,nz,2))
    do it=1,2
      do i=1,nx*ny*nz
         DensityFactor(i,1,1,it) =  - PairingStrength(it) *  ( 1.0_dp - alpha(it)/rhosat* &
         &                         (Density%Rho(i,1,1,1) + Density%Rho(i,1,1,2)))
      enddo
    enddo
  end subroutine CompDensityFactor
  
!   pure function PairingInter(wf1,wf2, it) result( ActionOfPairing)
!     !---------------------------------------------------------------------------
!     ! Calculate the action of the pairing interaction on a twobody state,
!     ! represented by spinor 1 and 2.
!     !
!     ! So the result of this function is
!     ! Sum_{s} s  < r , s ; r, -s | V |psi_1, psi_2 >
!     ! as a function of both spin and position.
!     ! Note that it is not antisymmetrised.
!     !---------------------------------------------------------------------------
!     ! Currently the only part of a general operator V that is implemented 
!     ! is the Delta interaction with a density dependent part.
!     ! It is described in
!     !      S.J. Krieger et al., Nucl. Phys.A517 (1990) 275
!     !
!     ! This force is given by:
!     !   v12 = 0.5 * strength * ( 1 - alpha * rho/rhosat) * (1-ps12)*delta(r12)
!     !---------------------------------------------------------------------------
!     ! Note that no cutoff factors are taken into account in this routine!
!     !---------------------------------------------------------------------------
!     ! Indices of the array ActionOfPairing.
!     ! postion_x, position_y, position_z, real or imaginary
!     !---------------------------------------------------------------------------
!     use Densities, only : Density
    
!     type(Spinor), intent(in) :: wf1, wf2
!     integer, intent(in)      :: it
    
!     complex(KIND=dp)         :: ActionOfPairing(nx,ny,nz)
!     real(KIND=dp)            :: factor(nx,ny,nz), Temp(nx,ny,nz,2)
!     integer                  :: i


!     ActionOfPairing = 0.0_dp ; factor = 0.0_dp
    
!     if(.not.TRC) then 
!       !---------------------------------------------------------------------------
!       ! Real part
!       do i=1,nx*ny*nz
!         Temp(i,1,1,1) =                                                      &
!         &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,3,1)            &
!         &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,4,1)            &
!         &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,1,1)            &
!         &             + wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,2,1)  
!         !---------------------------------------------------------------------------
!         ! Imaginary Part 
!         Temp(i,1,1,2) =                                                      &
!         &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,4,1)            &
!         &             + wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,3,1)            &
!         &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,2,1)            &
!         &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,1,1)
!       enddo
!     else
!       !---------------------------------------------------------------------------
!       ! Real part
!       do i=1,nx*ny*nz
!         Temp(i,1,1,1) =                                                      &
!         &             - wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,1,1)            &
!         &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,2,1)            &
!         &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,3,1)            &
!         &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,4,1)  
!         !---------------------------------------------------------------------------
!         ! Imaginary Part 
!         Temp(i,1,1,2) =                                                      &
!         &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,2,1)            &
!         &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,1,1)            &
!         &             + wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,4,1)            &
!         &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,3,1)
!       enddo
!     endif

!     !---------------------------------------------------------------------------
!     ! Density dependent interaction factor.
!     do i=1,nx*ny*nz
!       factor(i,1,1) =  - PairingStrength(it) *  ( 1.0_dp - alpha(it)/rhosat* &
!       &      (Density%Rho(i,1,1,1) + Density%Rho(i,1,1,2)))
!     enddo
    
!      do i=1,nx*ny*nz
!        ActionOfPairing(i,1,1) = factor(i,1,1) * dcmplx(Temp(i,1,1,1), Temp(i,1,1,2))
!      enddo

!   end function PairingInter

end module PairingInteraction
