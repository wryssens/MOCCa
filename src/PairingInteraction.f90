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
  real(KIND=dp) :: alpha(2) = 1.0_dp
  ! Exponent of the density dependence of the pairing force.
  real(KIND=dp) :: DDexp(2) = 1.0_dp
  ! relative strength of the Fayans-type gradient term
  real(KIND=dp) :: PairingGradientStrength(2) = 0.0_dp
  ! flag for considering the contribution of the pairing EDF to the 
  ! single-particle potential U(r)
  logical       :: PairingContributionUpot = .false.
  ! flag for (simple) old-school ULB pairing
  logical       :: PairingULB = .false.
  ! Energy associated with pairing obtained from direct integration of EDF
  real(KIND=dp) :: PairingEDF(2) = 0.0_dp
  ! invidual contribution to pairing EDF
  real(KIND=dp) :: PairingEDFDecomposition(3,2) = 0.0_dp
  ! averaged pairing gaps. The boolean flag signals if average pairing gaps
  ! can be defined (and therefore printed) or not.
  real(KIND=dp) :: AvGapuv  (2) = 0.0_dp
  real(KIND=dp) :: AvGapv2  (2) = 0.0_dp
  real(KIND=dp) :: ReSum_uv (2) = 0.0_dp
  real(KIND=dp) :: ImSum_uv (2) = 0.0_dp
  logical       :: PrintAvGaps = .false.

  ! Determine the kind of cutoff used.
  ! CutType
  ! 1) Symmetric Fermi
  ! 2) Cosine
  integer       :: CutType=1

  ! Density-dependent factor for density-dependent pairing
  real(KIND=dp), allocatable :: DensityFactor(:,:,:,:)
  ! Contribution from pair EDF  to single-particle potential Upot
  real(KIND=dp), allocatable :: UpotPairingContribution(:,:,:,:)
  real(KIND=dp)              :: UpotPairingRearrangementEnergy = 0.0_dp

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

  pure function GetPairDensity(wf1,wf2,TR) result( ActionOfPairing)
    !---------------------------------------------------------------------------
    ! Note: this is not returning the pair density, but the two-body wave
    ! function s entering the pair density.
    !---------------------------------------------------------------------------
    ! Calculate the action of the pairing interaction on a twobody state,
    ! represented by spinor 1 and 2.
    !
    ! So the result of this function is
    ! Sum_{s} s  < r , s ; r, -s | psi_1, psi_2 >
    ! 
    !---------------------------------------------------------------------------
    ! This routine includes a Time-reversal operator if asked for, by option
    ! TR.
    !---------------------------------------------------------------------------
    ! Indices of the array ActionOfPairing.
    ! postion_x, position_y, position_z, real or imaginary
    !---------------------------------------------------------------------------

    type(Spinor), intent(in) :: wf1, wf2
    
    complex(KIND=dp)         :: ActionOfPairing(nx,ny,nz)
    real(KIND=dp)            :: Temp(nx,ny,nz,2), l1,l2,l3,l4,r1,r2,r3,r4
    integer                  :: i
    logical, intent(in)      :: TR
    
    if(.not. TR) then 
      !---------------------------------------------------------------------------
      ! Real part
      do i=1,nx*ny*nz
        Temp(i,1,1,1) =                                                      &
        &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,3,1)            &
        &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,4,1)            &
        &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,1,1)            &
        &             + wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,2,1) 
      enddo

      do i=1,nx*ny*nz
        !---------------------------------------------------------------------------
        ! Imaginary Part 
        Temp(i,1,1,2) =                                                      &
        &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,4,1)            &
        &             + wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,3,1)            &
        &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,2,1)            &
        &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,1,1)
      enddo

      do i=1,nx*ny*nz
        ActionOfPairing(i,1,1) = dcmplx(Temp(i,1,1,1), Temp(i,1,1,2))
      enddo
    else
!      do i=1,nx*ny*nz
!        Temp(i,1,1,1) =                                                      &
!        &             - wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,1,1)            &
!        &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,2,1)            &
!        &             - wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,3,1)            &
!        &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,4,1) 
!        
!        Temp(i,1,1,2) =                                                      &
!        &               wf1%Grid(i,1,1,1,1) * wf2%Grid(i,1,1,2,1)            &
!        &             - wf1%Grid(i,1,1,2,1) * wf2%Grid(i,1,1,1,1)            &
!        &             + wf1%Grid(i,1,1,3,1) * wf2%Grid(i,1,1,4,1)            &
!        &             - wf1%Grid(i,1,1,4,1) * wf2%Grid(i,1,1,3,1)
!      enddo
!        do i=1,nx*ny*nz
!          ActionOfPairing(i,1,1) = dcmplx(Temp(i,1,1,1), Temp(i,1,1,2))
!        enddo
      do i=1,nx*ny*nz
        l1 = wf1%Grid(i,1,1,1,1) ; r1 = wf2%Grid(i,1,1,1,1)
        l2 = wf1%Grid(i,1,1,2,1) ; r2 = wf2%Grid(i,1,1,2,1)
        l3 = wf1%Grid(i,1,1,3,1) ; r3 = wf2%Grid(i,1,1,3,1)
        l4 = wf1%Grid(i,1,1,4,1) ; r4 = wf2%Grid(i,1,1,4,1)
        
        ActionofPairing(i,1,1) = dcmplx( & 
        & -l1*r1 - l2*r2 - l3*r3 - l4*r4, l1*r2 - l2*r1 + l3*r4 - l4*r3 )
      enddo
    endif

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
    !    This is described in S.J. Krieger et al., Nucl. Phys.A517 (1990) 27.
    !
    ! Note that the papers by Bonche et al, NPA443 and Krieger et al, NPA517,
    ! contain many misprints concerning the square of the cutoff:
    ! - the square of f is missing in Eq. (14) of Bonche et al NPA443
    ! - the denominator of Eqns. (6), (7), and (8) of Krieger et al NPA517
    !   should contain f^2_i Delta^2_i
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: E, Lambda
    integer, intent(in)       :: it
    real(Kind=dp)             :: Up, Down
    
    Up   =     (E - Lambda - PairingCut(it))/PairingMu(it)
    Down =   - (E - Lambda + PairingCut(it))/PairingMu(it)
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
    ! print '(" CosineCut ",i2,7f10.5)',it,E,Lambda,Up,Down, &
    !   &PairingCut(it),PairingMu(it),Cutoff
  
  end function Cosinecut

  subroutine CompDensityFactor
    !---------------------------------------------------------------------------
    ! Calculate 4 times the local pairing form factor 
    ! (see comments below for the factor 4)
    !---------------------------------------------------------------------------
    ! Currently the only part of a general operator V that is implemented 
    ! is the Delta interaction with a density dependent part.
    ! there are two options:
    !
    ! 1) "ULB pairing" as described in C. Rigollet, P. Bonche, H. Flocard, and 
    ! P.-H. Heenen, PRC 59, 3120 (1999). The underlying force is given by
    !   v12 = 0.5 * (1-ps12) * strength * ( 1 - alpha * rho/rhosat) * delta(r12)
    ! but the code actually allows for a charge-symmetry-breaking generalisation
    ! with different parameters "strength" and "alpha" for protons and neutrons.
    ! For alpha = 1 and rhosat = 0.16 fm^-3 this is a pure "surface" pairing EDF
    ! "Volume" and "mixed" pairing are obtained alpha = 0 and alpha = 1/2,
    ! respectively.
    !
    ! 2) Fayans pairing as proposed in S. A. Fayans and D. Zawischa, Phys. Lett.
    ! B383, 19 (1996) and further advocated in P.-G. Reinhard and W. Nazarewicz,
    ! Phys. Rev. C 95, 064328 (2017)
    !
    ! Note: the position-dependent pairing strength is DensityFactor(:,:,:) / 4,
    ! with the latter factor being absorbed by whatever DensityFactor(:,:,:)
    ! multiplies elsewhere in the code.
    !---------------------------------------------------------------------------
    integer       :: i, it , itp
    real(KIND=dp) :: Temp(nx,ny,nz) , Rho0(nx,ny,nz) , NablaRho02(nx,ny,nz)
    real(KIND=dp) :: fac

    if(.not.allocated(DensityFactor)) allocate(DensityFactor(nx,ny,nz,2))
    Rho0(:,:,:) = Density%Rho(:,:,:,1) + Density%Rho(:,:,:,2)

    if ( PairingULB ) then
      !=========================================================================
      ! old-school ULB pairing EDF
      do it=1,2
        DensityFactor(:,:,:,it)  & 
          & = - PairingStrength(it) * (1.0_dp - alpha(it)/rhosat * Rho0(:,:,:))
      enddo
    else
      !=========================================================================
      ! generalized pairing EDF

      ! square of the gradient of the isoscalar density needed later
      NablaRho02 = 0.0_dp
      if (  PairingGradientStrength(1) .ne. 0.0_dp .or.  &
          & PairingGradientStrength(2) .ne. 0.0_dp ) then 
        ! don't know what to do if nabla rho(r) is not allocated
        if ( .not.allocated(Density%DerRho) ) &
           &  call stp(" CompDensityFactor: DerRho not allocated!")
        do i =1,3
          NablaRho02(:,:,:) = NablaRho02(:,:,:) &
             & +          Density%DerRho(:,:,:,i,1)*Density%DerRho(:,:,:,i,1) &
             & +          Density%DerRho(:,:,:,i,2)*Density%DerRho(:,:,:,i,2) &
             & + 2.0_dp * Density%DerRho(:,:,:,i,1)*Density%DerRho(:,:,:,i,2)
        enddo
      endif
    ! print '(" rho D rho_t ",f12.5)', -90.806959 * dv*sum(NablaRho02(:,:,:))

      do it=1,2
        !-----------------------------------------------------------------------
        ! density-dependent part
        ! The possibility of alpha(1) != alpha(2) and ddexp(1) != ddexp(2) makes
        ! this slightly complicated.
        ! Note that the prefactor 1/4 of the EDF is taken care of by what this
        ! multiplies when calculating the Delta_ij
        if ( ddexp(it) .ne. 1.0_dp ) then
          fac  = alpha(it) / (rhosat**(ddexp(it)))
          Temp = Rho0**(ddexp(it)) 
        else
          fac  = alpha(it) / rhosat
          Temp = Rho0
        endif
        DensityFactor(:,:,:,it)  &
          & = - PairingStrength(it) *  ( 1.0_dp - fac * Temp(:,:,:) )

        !-----------------------------------------------------------------------
        ! gradient terms
        ! Note that the prefactor 1/4 of the EDF is taken care of by what this
        ! multiplies when calculating the Delta_ij
        if ( PairingGradientStrength(it) .ne. 0.0_dp ) then
          fac = - PairingGradientStrength(it) / (rhosat*rhosat)
          ! [nabla rho_n].[nabla rho_n] + [nabla rho_p].[nabla rho_p]
          DensityFactor(:,:,:,it) = DensityFactor(:,:,:,it) & 
            & + fac * NablaRho02(:,:,:)
        endif
      enddo
    endif
  end subroutine CompDensityFactor

  subroutine CompUpotPairingContribution(PairDensity)
    !---------------------------------------------------------------------------
    ! Calculate the contribution of the pairing EDF to the single-particle
    ! Hamiltonian
    !
    ! U_{pair,q}(r) = dE_{pair}/d rho_q(r)
    !
    ! See subroutine CompDensityFactor for further information about the EDFs.
    ! Note that for the pairing EDFs considered so far, which only depend on the
    ! isoscalar density rho_0(r), the contribution to U_{pair,q}(r) is the same 
    ! for protons and neutrons. By contrast, the contribution from protons and 
    ! neutrons might be different depending on choices made for coupling 
    ! constants of the pairing EDF.
    !---------------------------------------------------------------------------
    complex(KIND=dp), intent(in) :: PairDensity(nx,ny,nz,2)
    integer       :: i, it 
    real(KIND=dp) :: PairRho2(nx,ny,nz)
    real(KIND=dp) :: Rho0(nx,ny,nz)
    real(KIND=dp) :: NablaTerms(nx,ny,nz)
    real(KIND=dp) :: Temp(nx,ny,nz)
    real(KIND=dp) :: fac , eps

    ! nothing to do - what am I doing here ...?
    if ( .not.PairingContributionUpot ) return

    ! check if allocated - these should have been calculated
    if(.not.allocated(DensityFactor)) &
      & call stp('CompUpotPairingContribution: DensityFactor!') 

    ! check if allocated (which might be not the case during the first call)
    if ( .not.allocated(UpotPairingContribution)) then
      allocate(UpotPairingContribution(nx,ny,nz,2))
    endif

    !===========================================================================
    ! sum up the contributions from both isospins
    UpotPairingContribution        = 0.0_dp
    UpotPairingRearrangementEnergy = 0.0_dp

    if ( PairingULB ) then
      !=========================================================================
      ! old-school ULB pairing EDF
      do it=1,2   ! note: the loop is over contributions from it (and not to it)
        PairRho2(:,:,:) = DBLE (PairDensity(:,:,:,it))**2 &
          &              +DIMAG(PairDensity(:,:,:,it))**2
        ! note: there are two hidden signs in PairingStrength and alpha
        fac = 0.25 * PairingStrength(it) * alpha(it) / rhosat
        UpotPairingContribution(:,:,:,1)  &
          & = UpotPairingContribution(:,:,:,1) + fac * PairRho2(:,:,:)
      enddo
      UpotPairingContribution(:,:,:,2) = UpotPairingContribution(:,:,:,1)

    else
      !=========================================================================
      ! more general density dependence and gradient terms of the density
      Rho0(:,:,:) = Density%Rho(:,:,:,1) + Density%Rho(:,:,:,2)

      ! This parameter stabilizes the calculation of Rho**(ddexp - 1)
      ! Since usually ddexp < 1, this becomes unstable for small values of rho. 
      ! Thus we calculate it as Rho**ddexp/(Rho + eps) as is done in CalcUPot()  
      ! of MeanFields.f90
      eps = 1.d-20

      ! Note: the following loop is over the isospin of the term in the EDF,
      ! not the isospin of the potential it contributes to. For EDFs 
      ! considered so far, the contribution to the proton and neutron potential
      ! is equal.
      do it=1,2 
        PairRho2(:,:,:) = DBLE (PairDensity(:,:,:,it))**2 &
          &              +DIMAG(PairDensity(:,:,:,it))**2

        !-----------------------------------------------------------------------
        ! density-dependent contribution
        if ( ddexp(it) .ne. 1.0_dp ) then
          fac  = ddexp(it) / (rhosat**(ddexp(it)))
          Temp = PairRho2 * Rho0**(ddexp(it))/(Rho0+eps)
        else
          fac  = 1.0_dp / rhosat
          Temp = PairRho2
        endif
        ! there are two hidden "-" signe here, one from PairingStrength(it) ,
        ! the other from alpha(it) 
        fac = 0.25_dp * PairingStrength(it) * alpha(it) * fac
        UpotPairingContribution(:,:,:,1) &
        & = UpotPairingContribution(:,:,:,1) + fac * Temp(:,:,:)
     !  print '(" surface terms ",f16.4)', fac * dv * sum(Temp(:,:,:)*Rho0(:,:,:))

        !-----------------------------------------------------------------------
        ! the two contributions from gradient terms
        ! nabla.{ rho~*(r) rho~(r) [nabla rho_0(r)]}
        !  = {nabla[rho~*(r) rho~(r)]}.[nabla rho_0(r)]
        !       + rho~*(r) rho~(r)[Delta rho_0(r)]
        ! The factor -1/4 is by convention of the coupling constants
        ! where the "-" cancels the global sign of this contribution
        !-----------------------------------------------------------------------
        if ( PairingGradientStrength(it) .ne. 0.0_dp ) then
          fac = 2.0_dp * 0.25_dp * PairingGradientStrength(it) / (rhosat*rhosat)
          NablaTerms = 0.0_dp
          Temp(:,:,:) = DeriveX(PairRho2(:,:,:), ParityInt,SignatureInt,TimeSimplexInt,1)
          NablaTerms(:,:,:) =                     Temp(:,:,:) * (Density%DerRho(:,:,:,1,1) + Density%DerRho(:,:,:,1,2))
          Temp(:,:,:) = DeriveY(PairRho2(:,:,:), ParityInt,SignatureInt,TimeSimplexInt,1)
          NablaTerms(:,:,:) = NablaTerms(:,:,:) + Temp(:,:,:) * (Density%DerRho(:,:,:,2,1) + Density%DerRho(:,:,:,2,2))
          Temp(:,:,:) = DeriveZ(PairRho2(:,:,:), ParityInt,SignatureInt,TimeSimplexInt,1)
          NablaTerms(:,:,:) = NablaTerms(:,:,:) + Temp(:,:,:) * (Density%DerRho(:,:,:,3,1) + Density%DerRho(:,:,:,3,2))
          NablaTerms(:,:,:) = NablaTerms(:,:,:) + PairRho2(:,:,:) * (Density%LapRho(:,:,:,1) + Density%LapRho(:,:,:,2))
          UpotPairingContribution(:,:,:,1) = UpotPairingContribution(:,:,:,1) + fac * NablaTerms(:,:,:)
       !  print '(" NablaTerms ",f16.4)', fac * dv * sum(NablaTerms(:,:,:)*Rho0(:,:,:))
        endif 
      enddo
      ! as the form factor only depends on the isoscalar density, the contribution 
      ! to the proton potential equals the contribution to the neutron potential
      UpotPairingContribution(:,:,:,2) = UpotPairingContribution(:,:,:,1)
    endif

    !--------------------------------------------------------------------------
    ! calculate the corresponding rearrangement energy as
    !     1/2 int d^3 rho_0(r) sum_q U_{pair,q}(r)
    ! which has to be _s u b t r a c t e d_ from 1/2 sum_k epsilon_k
    ! As UpotPairingContribution(:,:,:,1) = UpotPairingContribution(:,:,:,2)
    ! this could be simplified, but the structure anticipates already form
    ! factors depending on isovector densities.
    fac = 0.5_dp * dv
    UpotPairingRearrangementEnergy &
      & = fac * sum(UpotPairingContribution(:,:,:,1)*Density%Rho(:,:,:,1)) &
      &  +fac * sum(UpotPairingContribution(:,:,:,2)*Density%Rho(:,:,:,2))

  end subroutine CompUpotPairingContribution

  subroutine CompPairingEDF(PairDensity)
    !---------------------------------------------------------------------------
    ! calculate pairing energy directly by integrating pairing EDF.
    !---------------------------------------------------------------------------
    ! At present only used as additional cross-check.
    !
    ! Also provides decomposition into the individual terms: 
    ! PairingEDFDecomposition(1,it): rho~*(r)rho~(r) 
    ! PairingEDFDecomposition(2,it): rho~*(r)rho~(r) rho_0^ddexp
    ! PairingEDFDecomposition(3,it): rho~*(r)rho~(r) [nabla rho_0].[nabla rho_0]
    !---------------------------------------------------------------------------
    complex(KIND=dp), intent(in) :: PairDensity(nx,ny,nz,2)
    real(KIND=dp)                :: PairRho2(nx,ny,nz)
    real(KIND=dp)                :: Rho0(nx,ny,nz)
    real(KIND=dp)                :: Temp(nx,ny,nz)
    real(KIND=dp)                :: NablaRho02(nx,ny,nz)
    real(KIND=dp)                :: fac
    integer                      :: i, it
    integer :: ix,iy,iz
    real(KIND=dp) :: summe

    if(.not.allocated(DensityFactor)) call stp('CompPairingEDF: DensityFactor!')

    !---------------------------------------------------------------------------
    ! square of the gradient of the isoscalar density
    NablaRho02 = 0.0_dp
    if (  PairingGradientStrength(1) .ne. 0.0_dp .or.  &
        & PairingGradientStrength(2) .ne. 0.0_dp ) then
      do i =1,3
        NablaRho02(:,:,:) = NablaRho02(:,:,:) &
              & +          Density%DerRho(:,:,:,i,1)*Density%DerRho(:,:,:,i,1) &
              & +          Density%DerRho(:,:,:,i,2)*Density%DerRho(:,:,:,i,2) &
              & + 2.0_dp * Density%DerRho(:,:,:,i,1)*Density%DerRho(:,:,:,i,2)
      enddo
    endif
  ! print '(" rho D rho_t EDF ",f12.5)', -90.806959 * dv*sum(NablaRho02(:,:,:))

    ! the isoscalar density
    Rho0(:,:,:) = Density%Rho(:,:,:,1) + Density%Rho(:,:,:,2)

    !---------------------------------------------------------------------------
    PairingEDFDecomposition = 0.0_dp
    do it=1,2 
      PairRho2(:,:,:) = DBLE (PairDensity(:,:,:,it))**2 &
        &              +DIMAG(PairDensity(:,:,:,it))**2

      !-------------------------------------------------------------------------
      ! calculate the total pairing energy directly by integrating pairing EDF 
      ! obtained as 0.25 * rho~*_q(r) rho~_q(r) * DensityFactor_q(r)  
      PairingEDF(it) = 0.25_dp * dv*sum(DensityFactor(:,:,:,it) * PairRho2(:,:,:))

      !------------------------------------------------------------------------
      ! prefactor -1/4 is by convention
      fac = -0.25_dp * PairingStrength(it)
      PairingEDFDecomposition(1,it) = fac * dv * sum(PairRho2(:,:,:))

      !------------------------------------------------------------------------
      ! prefactor +1/4 is by convention
      if ( ddexp(it) .ne. 1.0_dp ) then
        fac  = 1.0_dp / (rhosat**(ddexp(it)))
        Temp = Rho0**(ddexp(it))
      else
        fac  = 1.0_dp / rhosat
        Temp = Rho0
      endif
      fac = 0.25_dp * PairingStrength(it) * alpha(it) * fac
      PairingEDFDecomposition(2,it) = fac * dv * sum(PairRho2(:,:,:) * Temp(:,:,:))
      !-------------------------------------------------------------------------
      ! prefactor -1/4 is by convention
      if ( PairingGradientStrength(it) .ne. 0.0_dp ) then
        fac = - 0.25_dp * PairingGradientStrength(it) / (rhosat*rhosat)
        PairingEDFDecomposition(3,it) = fac * dv * sum(PairRho2(:,:,:) * NablaRho02(:,:,:))
      endif
    enddo
  end subroutine CompPairingEDF

  subroutine CompAverageGaps(PairDensity)
    !---------------------------------------------------------------------------
    ! Calculate average gaps in case of (dynamical) time-reversal symmetry
    ! (otherwise gaps are complex ...)
    !---------------------------------------------------------------------------
    ! PairPot is the Pairing Potential in russian representation
    ! AvGapuv : K. Yoshida M. Yamagami, and K. Matsuyanagi NPA 779 (2006) 99, 
    ! in particular Eq. (7). Earlier definitions (the earliest probably being
    ! J. Sauvage-Letessier, P. Quentin, H. Flocard, NPA 370 (1981) 231.)
    ! of this quantity that involve kappa instead of rho~ usually only work 
    ! in certain phase conventions where the u and v summed over are positive.
    !  <Delta uv > = - tr h~ rho~ / tr rho~ 
    !              = 2 E_pair / tr rho~ 
    ! (the latter assuming that all terms are bilinear in pair densities)
    ! AvGapv2 : J. Dobaczewski et al, NPA422 (1984) 103, Eq. (5.5)
    !  <Delta v2> = - tr h~ rho / tr rho
    ! This routine also calculates tr rho~ for later diagnostic printing. 
    !---------------------------------------------------------------------------
    complex(KIND=dp), intent(in) :: PairDensity(nx,ny,nz,2)
    real(KIND=dp)                :: PairPot(nx,ny,nz)
    real(KIND=dp)                :: Temp(nx,ny,nz)
    real(KIND=dp)                :: fac , sumuv , sumv2
    integer                      :: i, it

    PrintAvGaps = .false.
    AvGapuv  = 0.0_dp
    AvGapv2  = 0.0_dp
    ReSum_uv = 0.0_dp
    ImSum_uv = 0.0_dp

    do it=1,2 
      Temp(:,:,:)  = DBLE(PairDensity(:,:,:,it))
      ReSum_uv(it) = dv * sum(Temp(:,:,:))
      Temp(:,:,:)  = DIMAG(PairDensity(:,:,:,it))
      ImSum_uv(it) = dv * sum(Temp(:,:,:))
    enddo

    ! don't know how to average complex gaps, so check that they are not.
    ! test is only relevant when time-reversal is broken by the code, but 
    ! conserved as self-consistent symmetry in a program run.
    if( .not. TRC) then
      if  ( abs(ImSum_uv(1)) .gt. 1.d-6 .or. abs(ImSum_uv(2)) .gt. 1.d-6 ) &
        return
    endif

    ! passed the test, so there will be something to print
    PrintAvGaps = .true. 
    
    do it=1,2 
      Temp   (:,:,:) = DBLE(PairDensity(:,:,:,it))
      sumuv = dv * sum(Temp(:,:,:))
      sumv2 = dv * sum(Density%Rho(:,:,:,it))
      if ( abs(sumuv) .gt. 1.d-8 ) then
        PairPot(:,:,:) = 0.5_dp * DensityFactor(:,:,:,it)*Temp(:,:,:)
        AvGapuv(it) = - dv * sum(PairPot(:,:,:)*Temp(:,:,:))           / sumuv
        AvGapv2(it) = - dv * sum(PairPot(:,:,:)*Density%Rho(:,:,:,it)) / sumv2
        ! cross-check: this should be the pairing energy
        ! print '(" CompAverageGaps Epair",i2,f12.5)', &
        !  &  it,-0.5_dp * dv*sum(PairPot(:,:,:)*Temp (:,:,:))
      endif
    enddo

  end subroutine CompAverageGaps

  subroutine PrintPairingEnergies
    !---------------------------------------------------------------------------
    ! Print decomposition of pairing energies
    !---------------------------------------------------------------------------
    1 format (  ' Pairing Terms')
    2 format (   45x, ' N ',9x, ' P ',9x, ' T ')
    5 format (  ' rho~* rho*                             ',2x,f10.5,2x,f10.5,2x,f10.5)
    6 format (  ' rho~* rho* rho_0^g                     ',2x,f10.5,2x,f10.5,2x,f10.5)
    7 format (  ' rho~* rho* [nabla rho_0].[nabla rho_0] ',2x,f10.5,2x,f10.5,2x,f10.5)
    8 format (  ' sums                                   ',2x,f10.5,2x,f10.5,2x,f10.5)
   10 format (/,' Pairing EDF                            ',2x,f10.5,2x,f10.5,2x,f10.5)
   11 format (/,' Re tr rho~                             ',2x,f10.5,2x,f10.5)
   12 format (  ' Im tr rho~                             ',2x,f10.5,2x,f10.5)
   13 format (  ' <Delta uv>                             ',2x,f10.5,2x,f10.5)
   14 format (  ' <Delta v2>                             ',2x,f10.5,2x,f10.5)
   15 format (  ' Pairing Rearrangement Energy           ',2x,f10.5)
    !---------------------------------------------------------------------------
    print  1
    print  2
    print  5, PairingEDFDecomposition(1,1:2),sum(PairingEDFDecomposition(1,:))
    print  6, PairingEDFDecomposition(2,1:2),sum(PairingEDFDecomposition(2,:))
    if ( .not.PairingULB ) then
      print 7, PairingEDFDecomposition(3,1:2),sum(PairingEDFDecomposition(3,:))
    endif
 !  print  8, PairingEDFDecomposition(1,1)  & 
 !       &   +PairingEDFDecomposition(2,1)  &
 !       &   +PairingEDFDecomposition(3,1), &
 !       &    PairingEDFDecomposition(1,2)  &
 !       &   +PairingEDFDecomposition(2,2)  &
 !       &   +PairingEDFDecomposition(3,2)
    print 10, PairingEDF,sum(PairingEDF(1:2))

    print 11, ReSum_uv (1:2)
    if (.not.TRC) print 12, ImSum_uv (1:2)
    if ( PrintAvGaps ) then
      print 13, AvGapuv(1:2)
      print 14, AvGapv2(1:2)
    endif
    if ( PairingContributionUpot ) print 15, UpotPairingRearrangementEnergy

  end subroutine PrintPairingEnergies

  subroutine PrintPairingProfiles(PairDensity)
    !---------------------------------------------------------------------------
    ! diagnostic printing of various pairing-related densities and fields
    !---------------------------------------------------------------------------
    use Mesh, only   : Mesh3D

    complex(KIND=dp), intent(in) :: PairDensity(nx,ny,nz,2)
    real(KIND=dp)                :: UpotPC(nx,ny,nz,2)
    integer                      :: i , j , k , it

    UpotPC = 0.0_dp
    if ( allocated(UpotPairingContribution)) then
      UpotPC = UpotPairingContribution 
    endif

    i = 1 ; if (.not.SC ) i = nx/2+1
    j = 1 ; if (.not.TSC) j = ny/2+1

    print '(/," Pairing-related densities")'
    print '(  " cut at x = ",f8.3," y = ",f8.3)',Mesh3D(1,i,j,1),Mesh3D(2,i,j,1)
    print '(  "   k   z        rho_n           rho_p       ", & 
     &        " rho~*_n rho~_n  rho~*_ rho~_p       F_n             F_p    ", &
     &        "         U_n             U_p")'
    do k=1,nz
      print '(i4,f8.3,8es16.8)',                                             &
       & k,Mesh3D(3,i,j,k),                                                  &
       & Density%Rho(i,j,k,1),Density%Rho(i,j,k,2),                          &
       & DBLE (PairDensity(i,j,k,1))**2 + DIMAG(PairDensity(i,j,k,1))**2,    &
       & DBLE (PairDensity(i,j,k,2))**2 + DIMAG(PairDensity(i,j,k,2))**2,    &
       & 0.25_dp * DensityFactor(i,j,k,1), 0.25_dp * DensityFactor(i,j,k,2), &
       & UpotPC(i,j,k,1),UpotPC(i,j,k,2)
    enddo
    print '(" ")'

  end subroutine PrintPairingProfiles 


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
!     !       C. Rigollet et al, Phys. Rev. C 59, 3120 (1999).
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
