    module Energy
!-------------------------------------------------------------------------------
! This module contains all the variables to compute the value of the energy
! functional: basically everything related to the energy.
!
!Subroutines included in this module:
!
!       CompEnergy()
!       SpwfEnergy()
!       ConverEnergy()
!       CompKinetic()
!       SkyrmeTerms()
!       recoupleTensorNLO() 
!       PrintEnergy()
!       EnergyEV8Style
!       EnergyCR8Style
!
!-------------------------------------------------------------------------------
  use CompilationInfo
  use GenInfo
  use Densities
  use Force
  use SpwfStorage
  use Pairing
  use Meanfields

  implicit none

  save

  ! All the different Energies
  real(KIND=dp),public :: Kinetic(2),  CoulombEnergy
  real(KIND=dp),public :: CoulombExchange,CoMCorrection(2,2)
  real(KIND=dp),public :: TotalEnergy, OldEnergy(7), SpEnergy
  real(KIND=dp),public :: LNEnergy(2), Routhian, OldRouthian(7)

  ! Two different ways of calculating and treating Skyrme terms
  real(KIND=dp) :: Skyrmeterms(36),BTerm(21), N2LOterms(32), N3LOterms(12)
  real(KIND=dp) :: Skyrmeterms2dd(2)

  ! Recoupled tensor terms
  real(KIND=dp) :: Tensorterms(3,6)

  ! Signal to the Nesterov iteration whether or not the energy decreased
  integer       :: NesterovSignal=0

  abstract interface
    subroutine PrintEnergy_interface(Lagrange)
      logical, intent(in), optional :: Lagrange
    end subroutine
  end interface

  abstract interface
    subroutine CompEnergy_interface()
    end subroutine
  end interface
    
  procedure(PrintEnergy_interface),pointer :: PrintEnergy

contains

  subroutine CompEnergy()
    !---------------------------------------------------------------------------
    ! Subroutine that computes all the different energies of the main program
    ! state.
    !---------------------------------------------------------------------------
    ! Note that the TotalEnergy would give a slightly different result, due to
    ! the different computation of the kinetic energy
    !---------------------------------------------------------------------------
    use Coulomb
    use Pairing
    use Cranking
    use Moments
    
    integer :: i
    
    ! Making sure the PrintEnergy routine is associated
    if(.not.associated(PrintEnergy)) then
        if(trim(SkyrmeTreatment).eq.'BTERMS') then
            PrintEnergy => PrintEnergy_Bterms
        else
            PrintEnergy => PrintEnergy_termbyterm
        endif
    endif

    !Initialise the contributions to the energy
    CoulombEnergy=0.0_dp
    CoMCorrection=0.0_dp
    BTerm=0.0_dp
    Kinetic=0.0_dp

    !Shift the entries in OldEnergy by one place and put in
    !the previous value of the Energy.
    do i=0,5
        OldEnergy(7-i) = OldEnergy(6-i)
    enddo
    OldEnergy(1)=TotalEnergy

    !Calculate the Kinetic Energy
    call CompKinetic

    if(trim(SkyrmeTreatment).eq.'BTERMS') then
        !Calculate the Skyrme Terms by the B-s
        Bterm  = compSkyrmeTerms(Density)
        TotalEnergy = sum(BTerm)
    else
        !Calculate the Skyrme terms term-by-term
        SkyrmeTerms = compSkyrme(Density)
        TotalEnergy = sum(SkyrmeTerms)
    endif

    if (J2Terms) then
       Tensorterms = recoupleTensorNLO(Density)
    endif

    ! Calculate the N2LO terms
    N2LOterms = N2LO(Density)

    ! Calculate the N3LO term
    N3LOterms = N3LO(Density)

    !Pairing Energy
    PairingEnergy = CompPairingEnergy(Delta)
    !Lipkin-Nogami Energy
    if(Lipkin) then
      LNEnergy      = - LNLambda * PairingDisp
    endif

    !Calculating the CoulombEnergy
    CoulombEnergy   = CompCoulombEnergy(Density)
    
    if(Cexchange .ge. 1) then
      ! Calculate the coulomb exchange energy, both if pertubatively included 
      ! or selfconsistently included
      CoulombExchange = CompCoulombExchange(Density)
    else
      CoulombExchange = 0
    endif
    
    !COM Correction
    call CompCOMCorrection

    !Sum of all energies (Note that the Skyrme sum already was included)
    TotalEnergy = TotalEnergy        + sum(Kinetic)  + CoulombEnergy  +        &
    &             sum(CoMCorrection) + CoulombExchange + sum(PairingEnergy) +  &
    &             sum(LNEnergy) + sum(N2LOterms)

    !Calculate the Routhian too
    do i=0,5
        OldRouthian(7-i) = OldRouthian(6-i)
    enddo
    OldRouthian(1) = Routhian
    Routhian = TotalEnergy                                                     &
    !                         Contribution of multipole moment constraints
    &                      +    sum(ConstraintEnergy*Density%Rho)*dv           &
    !                         Contribution of the cranking constraints (Omega*J)
    &                      +    sum(CrankEnergy)

    !Energy due to the single particle states.
    SpEnergy = SpwfEnergy()

    ! Signal
    if(TotalEnergy .gt. OldEnergy(1)) then
      NesterovSignal = 1
    else 
      NesterovSignal = 0
    endif
    return
  end subroutine CompEnergy

  function N2LO(Den) result(terms)
    !---------------------------------------------------------------------------
    ! Calculates the N2LO contribution to the energy functional.
    !
    !---------------------------------------------------------------------------
    use force
    
    real(KIND=dp) :: terms(32)
    real(KIND=dp) :: rhotot (nx,ny,nz),     laprhotot(nx,ny,nz)
    real(KIND=dp) :: Rtauten(nx,ny,nz,3,3), Itauten(nx,ny,nz,3,3)
    real(KIND=dp) :: Tmu    (nx,ny,nz,3,2), dmuJmunu(nx,ny,nz,3,2)
    real(KIND=dp) :: D2rho  (nx,ny,nz,3,3)
    type(densityvector), intent(in) :: Den
    integer       :: it,i,j
    
    terms = 0.0_dp
    
    !---------------------------------------------------------------------------
    ! Don't calculate anything if all of the coupling coefficients are zero
    if(.not.allocated(Den%Rtaun2lo) ) return
    
    !---------------------------------------------------------------------------
    ! Some sums for ease of coding
    rhotot   = sum(Den%Rho,   4)
    laprhotot= sum(Den%LapRho,4)
    
    Rtauten = sum(den%RtauN2LO,6)
    Itauten = sum(den%ItauN2LO,6)
    D2rho   = sum(den%D2Rho,   6)
    DmuJmunu = 0.0
    do i=1,3
        do j=1,3
            DmuJmunu(:,:,:,i,:) = DmuJmunu(:,:,:,i,:)+Den%DJmunu(:,:,:,j,j,i,:)
        enddo
    enddo
    Tmu = 0.0
    do i = 1,3
        do j=1,3
            Tmu(:,:,:,i,:) = Tmu(:,:,:,i,:) + Den%ReKN2LO(:,:,:,j,j,i,:) 
        enddo
    enddo
    !---------------------------------------------------------------------------
    
    !---------------------------------------------------------------------------
    ! Delta rho Delta rho                              T-even
    terms(1) = sum(laprhotot**2)                                 * N2D2rho(1)
    do it=1,2
        terms(2) =  terms(2) + sum(Den%LapRho(:,:,:,it)**2)      * N2D2rho(2) 
    enddo
    !----------------------------------------------------------------------------
    !  rho Q                                           T-even 
    terms(3) = sum(sum(Den%rho,4) * sum(Den%QN2LO,4))            * N2rhoQ(1)
    do it=1,2
        terms(4) = terms(4) + sum(Den%rho(:,:,:,it) * Den%QN2LO(:,:,:,it)) 
    enddo
    terms(4) = terms(4)                                          * N2rhoQ(2)
    
    !----------------------------------------------------------------------------
    !  tau^2                                           T-even
    terms(5) = sum(sum(Den%tau,4)**2)                            * N2tau(1)
    do it=1,2
        terms(6) = terms(6) + sum(Den%tau(:,:,:,it)**2)          * N2tau(2)
    enddo
    !---------------------------------------------------------------------------
    ! Re Tau_munu * Re Tau_munu                        T_even
    terms(7) = 2*sum(Rtauten*Rtauten)                            * N2rtaumn(1)
    do it=1,2
        terms(8) = terms(8)+2*sum(Den%RTauN2LO(:,:,:,:,:,it)**2) * N2rtaumn(2)
    enddo
    !---------------------------------------------------------------------------
    ! - Re tau_munu * DmuDnuRho                        T_even
    terms(9) =   -  2*sum(Rtauten * D2rho)                       * N2tddr(1)
    do it=1,2
        terms(10) = terms(10) &
        &              + sum(Den%RTauN2LO(:,:,:,:,:,it)*Den%D2Rho(:,:,:,:,:,it))
    enddo
    terms(10) =  -  2*terms(10)                                  * N2tddr(2)

    !---------------------------------------------------------------------------
    ! - Im Tau_munu * Im Tau_munu                            T_odd
    if(.not.TRC) then
        terms(11) = - 2* sum(Itauten * Itauten)                  * N2itaumn(1)
        do it=1,2
            terms(12) = terms(12) - sum(Den%ITauN2LO(:,:,:,:,:,it)**2) 
        enddo
        terms(12) = 2*terms(12)                                  * N2itaumn(2)
    endif
    !---------------------------------------------------------------------------
    ! D_mu J_munu                                            T-even
    terms(13) = sum(sum(DmuJmunu,5)**2)                          * N2DJ(1)
    do it=1,2
        terms(14) = terms(14) + sum(DmuJmunu(:,:,:,:,it)**2) 
    enddo
    terms(14) = terms(14)                                        * N2DJ(2)
    !---------------------------------------------------------------------------
    ! J_munu V_munu                                          T-even
    terms(15) =-4*sum(sum(Den%JmuNu,6)*sum(Den%VN2LO,6))         * N2JV(1)
    do it=1,2
        terms(16) = terms(16) + &
        &               4 * sum(Den%JmuNu(:,:,:,:,:,it)*Den%VN2LO(:,:,:,:,:,it))
    enddo
    terms(16) =-terms(16)                                        * N2JV(2)
    !---------------------------------------------------------------------------
    ! Delta s^2                                             T-odd
    if(.not.TRC) then
        terms(17) = sum(sum(Den%laps,5)**2)                      * N2D2s(1)
        do it=1,2
            terms(18) = terms(18)+sum(Den%laps(:,:,:,:,it)**2)   * N2D2s(2)
        enddo
    endif
    !---------------------------------------------------------------------------
    ! D_mu j_mu s^2                                         T-odd
    ! This should be zero for local interactions by the way
    if(.not.TRC) then
        terms(19) = - sum(sum(Den%divvecj,4)**2)                 * N2Dvecj(1)
        do it=1,2
            terms(20) = terms(20)-sum(Den%divvecj(:,:,:,it)**2)  * N2Dvecj(2)
        enddo
    endif
    !---------------------------------------------------------------------------
    ! j_mu Pi_mu                                            T-odd
    if(.not.TRC) then
        terms(21) = -4*sum(sum(den%vecj,5) * sum(den%PiN2LO,5))  * N2jpi(1)
        do it=1,2
            terms(22) = terms(22) - &
            & 4*sum(den%vecj(:,:,:,:,it)*den%PiN2LO(:,:,:,:,it)) * N2jpi(2)
        enddo
    endif
    !---------------------------------------------------------------------------
    ! s_mu S_mu                                            T-odd
    if(.not.TRC) then
        terms(23) = sum(sum(Den%vecs,5) * sum(Den%SN2LO,5))      * N2sS(1)
        do it=1,2
            terms(24) = terms(24) + &
            & sum(Den%vecs(:,:,:,:,it) * Den%SN2LO(:,:,:,:,it))  * N2sS(2)
        enddo
    endif
    !---------------------------------------------------------------------------
    ! Tmu^2                                                T-odd
    if(.not.TRC) then
        terms(25) = sum(sum(Tmu,5)**2)                           * N2vecT(1)
        do it=1,2
            terms(26) = terms(26) + sum(Tmu(:,:,:,:,it)**2)      * N2vecT(2)
        enddo
    endif
    !---------------------------------------------------------------------------
    ! Re T_munuka^2                                        T-odd
    if(.not.TRC) then
        terms(27) =   2 * sum(sum(Den%ReKN2LO,7)**2)             * N2ReTmn(1)
        do it=1,2
            terms(28) = terms(28) + sum(Den%ReKN2LO(:,:,:,:,:,:,it)**2)
        enddo
        terms(28)= 2*terms(28)                                   * N2ReTmn(2)
    endif
    !---------------------------------------------------------------------------
    ! DmunuS * Re T_munuka
    if(.not.TRC) then
        terms(29) = -2*sum(sum(Den%ReKN2LO,7)*sum(Den%D2S,7))  * N2TmnD2s(1)
        do it=1,2
            terms(30) = terms(30) + &
            &sum(Den%ReKN2LO(:,:,:,:,:,:,it)*Den%D2S(:,:,:,:,:,:,it))
        enddo
        terms(30) = -2*terms(30)                               * N2TmnD2s(2)
    endif
    !---------------------------------------------------------------------------
    ! Im T_munuka^2                   T-even (!!!)
    terms(31)     = - 2 * sum(sum(Den%ImKN2LO,7)**2)           * N2ImTmn(1)
    do it=1,2
        terms(32) = terms(32) - &
        &              2*sum(Den%ImKN2LO(:,:,:,:,:,:,it)**2)   * N2ImTmn(2)
    enddo
    !---------------------------------------------------------------------------
    ! Multiply everyting by dv
    terms = terms * dv 
  end function N2LO
  
  function N3LO(Den) result(terms)
    !---------------------------------------------------------------------------
    ! Calculate the contributions of the N3LO densities to the functional.
    ! 
    !
    !---------------------------------------------------------------------------
    integer                         :: it
    real(KIND=dp)                   :: terms(12)
    type(DensityVector), intent(in) :: Den
    real(KIND=dp)                   :: laptau(nx,ny,nz,2), laplaplaprho
    
    terms = 0.0_dp

    if(t1n3.eq.0.0_dp .and. t2n3.eq.0.0_dp) return

    !---------------------------------------------------------------------------
    ! Intermediate stuff
    !---------------------------------------------------------------------------
    laptau = Den%DDRtaumn(:,:,:,1,1,:) + Den%DDRtaumn(:,:,:,2,2,:)             &
    &      + Den%DDRtaumn(:,:,:,3,3,:) 

    !---------------------------------------------------------------------------
    ! rho Delta Delta Delta rho                                        Time-even
    terms(1) = sum(sum(Den%rho,4) * sum(Den%D3rho,4))             *N3D3rho(1)
    do it=1,2
      terms(2)=terms(2)+sum(Den%rho(:,:,:,it)*Den%D3rho(:,:,:,it))*N3D3rho(2)
    enddo
    !---------------------------------------------------------------------------
    ! ReTau delta ReTau term                                           Time-even
    terms(3) = -3.0/2.0 * sum(sum(Den%tau,4) * sum(laptau,4))   *N3tauDtau(1)
    do it=1,2
        terms(4) = terms(4)   &
        &      -3.0/2.0 *sum(Den%tau(:,:,:,it)*laptau(:,:,:,it))*N3tauDtau(2)
    enddo
    !---------------------------------------------------------------------------
    ! Re Tau_munu Delta Re Tmunu
    terms(5) = -3.0/2.0 * sum(sum(Den%RTauN2LO,6) * sum(Den%DDRTaumn,6))       &
    &                                                            *N3tauDmntau(1)
    do it=1,2
        terms(6) = terms(6)   &
        & -3.0/2.0 *sum(Den%RtauN2LO(:,:,:,:,:,it)*Den%DDRTaumn(:,:,:,:,:,it)) &
        &                                                        *N3tauDmntau(2)
    enddo

    !---------------------------------------------------------------------------
    ! Delta rho Dmu Dnu Re Tau_munu                                  Time-even
    terms(7) = 3d0/2d0*sum(sum(Den%laprho,4) * sum(Den%D2RTau,4))* N3DrhoDtau(1)
    do it=1,2
        terms(8) = terms(8) + 3d0/2d0*&
        &      sum(Den%laprho(:,:,:,it) * Den%D2RTau(:,:,:,it))  * N3DrhoDtau(2)
    enddo
    !---------------------------------------------------------------------------
    ! Tau Dmu Dnu Re Tau_munu                                         Time-even
    terms(9) = -3 *sum(sum(Den%tau,4) * sum(Den%D2RTau,4))      * N3tauDmntau(1)
    do it=1,2
        terms(10) = terms(10) - 3* &
        &      sum(Den%tau(:,:,:,it) * Den%D2RTau(:,:,:,it))    * N3tauDmntau(2)
    enddo
    !---------------------------------------------------------------------------
    ! tau Delta Delta rho                                             Time-even
    terms(11) = +3d0/2d0 * sum(sum(Den%tau,4)*sum(Den%laplaprho,4))            &
    &                                                           *  N3tauDDrho(1)
    do it=1,2
        terms(12) = terms(12) + 3d0/2d0 *                                      &
        &          sum(Den%tau(:,:,:,it)*Den%laplaprho(:,:,:,it))              &
        &                                                       *  N3tauDDrho(2)
    enddo
    
    terms = terms*dv 
  end function

  function compSkyrmeTerms(Den) result(B)
    !---------------------------------------------------------------------------
    ! Calculates all the different Skyrme contributions.
    !    See P.Bonche, H.Flocard, P.H. Heenen; Nucl. Phys. A487(1987), 115-135
    !        V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
    !---------------------------------------------------------------------------
    type(DensityVector), intent(in) :: Den
    real(KIND=dp)                   :: B(21)
    integer                         :: it, l, m

    !Temporary storage for the total densities, as these figure quite often in
    !the expressions.
    real(KIND=dp) :: RhoT(nx,ny,nz), NablaJT(nx,ny,nz), TauT(nx,ny,nz)
    real(KIND=dp) :: vecJT(nx,ny,nz,3), LapRhoT(nx,ny,nz), RotST(nx,ny,nz,3)
    real(KIND=dp) :: JMuNuT(nx,ny,nz,3,3), VecTT(nx,ny,nz,3), LapST(nx,ny,nz,3)
    real(KIND=dp) :: DivST(nx,ny,nz), VecST(nx,ny,nz,3)
    real(KIND=dp) :: VecFT(nx,ny,nz,3), ref

    RhoT    = sum(Den%Rho,4)    ; TauT    = sum(Den%Tau,4)
    NablaJT = sum(Den%NablaJ,4) ; LapRhoT = sum(Den%LapRho,4)
    if(allocated(Den%JmuNu)) JMuNuT = sum(Den%JMuNu,6)

    if(.not.TRC) then
      VecJT=sum(Den%VecJ,5) ;  VecST=sum(Den%VecS,5) ; RotST=sum(Den%RotS,5)
      if(allocated(Den%VecT))  VecTT = sum(Den%VecT,5)
      if(allocated(Den%VecF))  VecFT = sum(Den%VecF,5)
      if(allocated(Den%LapS))  LapST = sum(Den%LapS,5)
      if(allocated(Den%DivS))  DivST = sum(Den%DivS,4)
    endif
    B =0.0_dp

    !B1 Terms
    B(1) = B1*sum(RhoT**2)*dv

    !B2 Term
    do it=1,2
      B(2) = B(2) + sum(Den%Rho(:,:,:,it)**2)
    enddo
    B(2) = B2*B(2)*dv

    !B3 Term
    B(3) = B3*(sum(RhoT*TauT))*dv
    if(.not.TRC) then
      B(3) = B(3) - B3*dv*sum(VecJT**2)
    endif
    !B4 Term
    do it=1,2
      B(4) = B(4) + sum(Den%Rho(:,:,:,it)*Den%Tau(:,:,:,it))
      if(.not.TRC) then
         B(4) = B(4) - sum(Den%vecj(:,:,:,:,it)**2)
       endif
    enddo
    B(4) = B4*dv*B(4)

    !B5 Terms
    B(5) = B5*sum(RhoT*LapRhoT)*dv

    !B6 Terms
    do it=1,2
      B(6) = B(6) + sum(Den%Rho(:,:,:,it)*Den%LapRho(:,:,:,it))
    enddo
    B(6) = B6*dv*B(6)

    !B7 Terms
    B(7) = B7a*sum(RhoT**(2+byt3a))*dv

    if(B7b.ne. 0.0_dp .or. B7b.ne. 0.0_dp) &
     & call stp(' compSkyrmeTerms not prepared for 2nd density dependence!')

    !B8 Terms
    do it=1,2
     B(8) = B(8) + sum(Den%Rho(:,:,:,it)**2*RhoT**byt3a)
    enddo
    B(8) = B8a*dv*B(8)

    !B9 Terms
    B(9) = sum(RhoT*NablaJT)
    if(.not.TRC) then
      B(9) = B(9) + sum(sum(vecJT*RotST,4))
    endif
    B(9) = B(9) * B9
    do it=1,2
      B(9) = B(9) + B9q*(sum(Den%Rho(:,:,:,it)*Den%NablaJ(:,:,:,it)))
      if(.not.TRC) then
        B(9) = B(9) + B9q* sum(sum(Den%vecJ(:,:,:,:,it)*Den%RotS(:,:,:,:,it),4))
      endif
    enddo
    B(9) = B(9) * dv

    if(.not.TRC) then
      B(10)= B10*sum(VecST**2)*dv
      do it=1,2
        B(11)= B(11) + B11*sum(Den%VecS(:,:,:,:,it)**2)*dv
      enddo
      !B12 Terms
      B(12)= B12a*sum( RhoT**byt3a*(vecsT(:,:,:,1)**2+vecsT(:,:,:,2)**2 +    &
                           VecsT(:,:,:,3)**2) )*dv

      !B13 Terms
      do it=1,2
       B(13)=  B(13) + sum(RhoT**byt3a*(Den%vecs(:,:,:,1,it)**2+             &
                               Den%vecs(:,:,:,2,it)**2+Den%Vecs(:,:,:,3,it)**2))
      enddo
      B(13) = B13a*dv*B(13)
    endif

    !B14 Terms
    if(B14.ne. 0.0_dp) then
      B(14)=(sum(JMuNuT**2))
      if(.not.TRC) then
        B(14) = B(14)  - sum(VecST(:,:,:,1)*VecTT(:,:,:,1) +                  &
          &                  VecST(:,:,:,2)*VecTT(:,:,:,2) +                  &
          &                  VecST(:,:,:,3)*VecTT(:,:,:,3))
      endif
      B(14) = B14*dv*B(14)
    endif

    !B15 Terms
    if(B15.ne. 0.0_dp) then
      B(15)= sum(Den%JMuNu(:,:,:,:,:,:)**2)
      ref = B15 * B(15) * dv
      if(.not.TRC) then
        B(15) = B(15) +(- sum(Den%VecS(:,:,:,1,:)*Den%VecT(:,:,:,1,:))    &
          &             - sum(Den%VecS(:,:,:,2,:)*Den%VecT(:,:,:,2,:))    &
          &             - sum(Den%VecS(:,:,:,3,:)*Den%VecT(:,:,:,3,:)) )
      endif
      B(15) = B15*dv*B(15)
    endif

    !B16 Terms
    if(B16.ne. 0.0_dp) then
      !Sum_{mu}J_{mumu}^2 terms
      B(16)= sum((JMuNuT(:,:,:,1,1)+JMuNuT(:,:,:,2,2)+JMuNuT(:,:,:,3,3))**2)
      do l=1,3
        do m=1,3
          ! Sum_{mu,nu} J_{mu,nu}J_{nu,mu}
          B(16) = B(16) + sum(JMuNuT(:,:,:,l,m)*JMuNuT(:,:,:,m,l))
        enddo
      enddo
      ! 2*s*F
      if(.not.TRC) then
        B(16) = B(16)  - 2.0*(sum(sum(VecST*VecFT,4)))
      endif
      B(16) = B(16)*B16*dv
    endif

    !B17 Terms
    if(B17.ne.0.0_dp) then
      B(17)= 0.0_dp
      do it=1,2
        !Sum_{mu}J_{mumu}^2 terms
        B(17) = B(17) + sum((Den%JMuNu(:,:,:,1,1,it)+Den%JMuNu(:,:,:,2,2,it)   &
        &             + Den%JMuNu(:,:,:,3,3,it))**2)
        do l=1,3
          do m=1,3
            ! Sum_{mu,nu} J_{mu,nu}J_{nu,mu}
            B(17) = B(17) + sum(Den%JMuNu(:,:,:,l,m,it)*Den%JMuNu(:,:,:,m,l,it))
          enddo
        enddo

        ! 2*s*F
        if(.not.TRC) then
          B(17) = B(17) - 2.0_dp*                                              &
          &                sum(sum(Den%VecS(:,:,:,:,it)*Den%VecF(:,:,:,:,it),4))
        endif
      enddo
      B(17)= B17*dv*B(17)
   endif

    !B18 Terms
    if(B18 .ne. 0.0_dp .and. .not.TRC) then
      B(18)= B18*dv*sum(sum(VecST*LapST,4))
    endif

    if(B19.ne. 0.0_dp .and..not.TRC) then
    !B19 Terms
      do it=1,2
        B(19)= B(19) &
        & + sum(sum(Den%VecS(:,:,:,:,it)*Den%LapS(:,:,:,:,it),4))
      enddo
      B(19) = B19*dv*B(19)
    endif
    if(B20 .ne. 0.0_dp .and. .not. TRC) then
      !B20 Terms
      B(20)= B20*dv*sum(DivST**2)
    endif

    if(B21 .ne. 0.0_dp .and. .not.TRC) then
      !B21 Terms
      B(21)= B21*dv*sum(Den%DivS(:,:,:,1)**2 + Den%DivS(:,:,:,2)**2)
    endif

    return
  end function compSkyrmeTerms

  function compSkyrme(Den) result(terms)

    type(DensityVector), intent(in) :: Den
    real(KIND=dp)                   :: terms(36)
    integer                         :: it, l, m

    !--------------------------------------------------------------------------
    !Temporary storage for the total densities, as these figure quite often in
    !the expressions.
    real(KIND=dp) :: RhoT(nx,ny,nz), NablaJT(nx,ny,nz), TauT(nx,ny,nz)
    real(KIND=dp) :: vecJT(nx,ny,nz,3), LapRhoT(nx,ny,nz), RotST(nx,ny,nz,3)
    real(KIND=dp) :: JMuNuT(nx,ny,nz,3,3), VecTT(nx,ny,nz,3), LapST(nx,ny,nz,3)
    real(KIND=dp) :: DivST(nx,ny,nz), VecST(nx,ny,nz,3)
    real(KIND=dp) :: VecFT(nx,ny,nz,3), ref

    RhoT    = sum(Den%Rho,4)    ; TauT    = sum(Den%Tau,4)
    NablaJT = sum(Den%NablaJ,4) ; LapRhoT = sum(Den%LapRho,4)
    if(allocated(Den%JmuNu)) JMuNuT = sum(Den%JMuNu,6)

    if(.not.TRC) then
      VecJT=sum(Den%VecJ,5) ;  VecST=sum(Den%VecS,5) ; RotST=sum(Den%RotS,5)
      if(allocated(Den%VecT))  VecTT = sum(Den%VecT,5)
      if(allocated(Den%VecF))  VecFT = sum(Den%VecF,5)
      if(allocated(Den%LapS))  LapST = sum(Den%LapS,5)
      if(allocated(Den%DivS))  DivST = sum(Den%DivS,4)
    endif

    Terms = 0

    ! B1 * \rho^2_{t}
    Terms(1) = B1*sum(RhoT**2)
    ! B2 * \rho^2_{q}
    do it=1,2
      Terms(2) = Terms(2) + sum(Den%Rho(:,:,:,it)**2)
    enddo
    Terms(2) = B2*Terms(2)

    ! B3 * \rho_t \tau_t
    Terms(3) = B3*(sum(RhoT*TauT))


    !B4 * \rho_t \tau_t
    do it=1,2
      Terms(4) = Terms(4) + B4* sum(Den%Rho(:,:,:,it)*Den%Tau(:,:,:,it))
    enddo

     ! B3 * \vecJ_{t}**2
    if(.not.TRC) Terms(5) = - B3*sum(VecJT**2)

    if(.not.TRC) then
        do it=1,2
            Terms(6) = Terms(6) - B4 * sum(Den%vecj(:,:,:,:,it)**2)
        enddo
    endif

    !B5 \rho_t \Delta \rho_t
    Terms(7) = B5*sum(RhoT*LapRhoT)

    !B6 * \rho_q \Delta \rho_q
    do it=1,2
      Terms(8) = Terms(8) + B6*sum(Den%Rho(:,:,:,it)*Den%LapRho(:,:,:,it))
    enddo

    !B7a \rho_t**(2+byt3a)
    Terms(9) = B7a*sum(RhoT**(2+byt3a))

    !B8a \rho_t**(byt3a) * \rho**2_q
    do it=1,2
     terms(10) = terms(10) + B8a*sum(Den%Rho(:,:,:,it)**2*RhoT**byt3a)
    enddo

    !B7b \rho_t**(2+byt3a) -- for evolutionary reasons added to end of array
    if(B7b .ne. 0.0_dp) then
      Terms(33) = B7b*sum(RhoT**(2+byt3b))
    endif

    !B8b \rho_t**(byt3b) * \rho**2_q -- for evolutionary reasons added to end of array
    if(B8b .ne. 0.0_dp) then
      do it=1,2
       terms(34) = terms(34) + B8b*sum(Den%Rho(:,:,:,it)**2*RhoT**byt3b)
      enddo
    endif

    !B9 * \rho_t \Nabla J_t
    Terms(11) = B9 * sum(RhoT*NablaJT)
    !B9q * \rho_q \Nabla J_q
    do it=1,2
      Terms(12) = Terms(12) + B9q*(sum(Den%Rho(:,:,:,it)*Den%NablaJ(:,:,:,it)))
    enddo

    !B9 * \vecJ_t * \Nabla \xtimes S_t
    if(.not.TRC) Terms(13) = B9*sum(sum(vecJT*RotST,4))

    ! B9q * \vecJ_q * \Nabla \xtimes \vecS_q
    if(.not.TRC) then
        do it=1,2
            Terms(14) = Terms(14) + B9q* sum(sum(Den%vecJ(:,:,:,:,it)*Den%RotS(:,:,:,:,it),4))
        enddo
    endif

    ! B10 \vecs^2_t
    if(.not.TRC) then
        Terms(15) = B10 *sum(VecST**2)
    endif

    ! B11 \vecs^2_q
    if(.not. TRC) then
        do it=1,2
            Terms(16)= Terms(16) + B11*sum(Den%VecS(:,:,:,:,it)**2)
        enddo
    endif

    !B12a \rho_t^{byt3a} * \vecs**2_t
    if(.not. TRC) then
        Terms(17)= B12a*sum( RhoT**byt3a*(sum(vecSt(:,:,:,:)**2,4)))
    endif

    !B13a \rho_t^{byt3a} * \vecs**2_q
    if(.not.TRC) then
        do it=1,2
            Terms(18)= Terms(18) + B13a* sum(RhoT**byt3a*(sum(Den%vecs(:,:,:,:,it)**2,4)))
        enddo
    endif

    !B12b \rho_t^{byt3b} * \vecs**2_t -- for evolutionary reasons added to end of array
    if(B12b .ne. 0.0_dp) then
      if(.not. TRC) then
        Terms(35)= B12b*sum( RhoT**byt3b*(sum(vecSt(:,:,:,:)**2,4)))
      endif
    endif

    !B13b \rho_t^{byt3b} * \vecs**2_q -- for evolutionary reasons added to end of array
    if(B13b .ne. 0.0_dp) then
      if(.not.TRC) then
        do it=1,2
          Terms(36)= Terms(36) + B13b* sum(RhoT**byt3b*(sum(Den%vecs(:,:,:,:,it)**2,4)))
        enddo
      endif
    endif

    !B14 \J_{mu mu}
    if(B14.ne. 0.0_dp) then
      Terms(19)= B14 * (sum(JMuNuT**2))
    endif

    ! B15 \sum J_{\mu \mu}_q
    if(B15.ne. 0.0_dp) then
      Terms(20)= B15 * sum(Den%JMuNu**2)
    endif

    ! - B14 \vecS * \vecT
    if(B14.ne. 0.0_dp) then
        if(.not. TRC) then
            Terms(21) = - B14 * sum(VecST*VecTT)
        endif
    endif

    ! B15 \vecs_q * \vecT_q
    if(B15.ne.0.0_dp .and. .not. TRC) then
        Terms(22) = - B15 * sum(Den%VecS*Den%VecT)
    endif

    ! B16 ( \sum J_{mu mu}_t ^2
    if(B16.ne.0.0_dp) then
        Terms(23) = B16 * sum((JMuNuT(:,:,:,1,1)+JMuNuT(:,:,:,2,2)+JMuNuT(:,:,:,3,3))**2)
    endif

     ! B17 ( \sum J_{mu mu}_q ^2
    if(B17.ne. 0.0_dp) then
        do it=1,2
            Terms(24) = Terms(24) + B17 * sum((Den%JMuNu(:,:,:,1,1,it)+ Den%JMuNu(:,:,:,2,2,it)    &
            &                                                         + Den%JMuNu(:,:,:,3,3,it))**2)
        enddo
    endif

    ! B16 \sum J_{\mu \nu}_q J_{\nu \mu}_q
    if(B16.ne.0.0_dp) then
        do l=1,3
            do m=1,3
                ! Sum_{mu,nu} J_{mu,nu}J_{nu,mu}
                Terms(25) = Terms(25) + B16 * sum(JMuNuT(:,:,:,l,m)*JMuNuT(:,:,:,m,l))
            enddo
        enddo
    endif

     ! B17 Sum_{mu,nu} J_{mu,nu}J_{nu,mu}
    if(B17.ne.0.0_dp) then
        do it=1,2
            do l=1,3
                do m=1,3
                    ! Sum_{mu,nu} J_{mu,nu}J_{nu,mu}
                    Terms(26) = Terms(26) + B17 * sum(Den%JMuNu(:,:,:,l,m,it)*Den%JMuNu(:,:,:,m,l,it))
                enddo
            enddo
        enddo
    endif

    ! B16 \vecs \vecF
    if(B16.ne. 0.0_dp .and. .not. TRC) then
        Terms(27) = - 2.0*B16*(sum(sum(VecST*VecFT,4)))
    endif


    ! B17 \vecs \vecF
    if(B17.ne.0.0_dp .and. .not. TRC) then
        do it=1,2
            Terms(28) = Terms(28) - 2.0_dp*B17*sum(Den%VecS(:,:,:,:,it)*Den%VecF(:,:,:,:,it))
        enddo
    endif

    !B18 s_t * \Delta s_t
    if(B18 .ne. 0.0_dp .and. .not.TRC) then
      Terms(29)= B18*sum(sum(VecST*LapST,4))
    endif

    !B19 s_q * \Delta s_q
    if(B19.ne. 0.0_dp .and..not.TRC) then
      do it=1,2
        Terms(30)= Terms(30)+B19*sum(Den%VecS(:,:,:,:,it)*Den%LapS(:,:,:,:,it))
      enddo
    endif

    ! B20 \nabla \cdot s_t
    if(B20 .ne. 0.0_dp .and. .not. TRC) then
      Terms(31)= B20*sum(DivST**2)
    endif

    ! B21 \nabla \cdot s_q
    if(B21 .ne. 0.0_dp .and. .not.TRC) then
      Terms(32)= B21*sum(Den%DivS(:,:,:,1)**2 + Den%DivS(:,:,:,2)**2)
    endif

    !-----------------
    ! Don't forget dv
    Terms = Terms * dv

  end function compSkyrme

  function recoupleTensorNLO(Den) result(Tterms)
    !--------------------------------------------------------------------------
    ! Recoupling of the cartesian NLO tensor terms into pseudoscalar, vector 
    ! and pseudotensor contributions. This function is an adaptation of what 
    ! has been implemented in EV4. Used for diagnostic purposes only. Only 
    ! called when true tensor terms are present in the EDF.
    ! For each of the three recoupled terms, there are 6 components for maximum 
    ! diagnostics: total (1), T=0 (2) T=1 (3) nn (4) pp (5) np (6). 
    ! These are separately calculated for the contributions from the central
    ! and the tensor force, but so far combined in the end for common printing.
    ! The explicit calculation of the (1) components is redundant, but kept
    ! for historical reasons.
    !--------------------------------------------------------------------------
    type(DensityVector), intent(in) :: Den
    real(KIND=dp)                   :: Tterms(3,6)
    integer                         :: it

    !--------------------------------------------------------------------------
    real(KIND=dp) :: J0(nx,ny,nz,2), J1(nx,ny,nz,3,2), J2(nx,ny,nz,3,3,2) 
    real(KIND=dp) :: ref
    real(KIND=dp) :: er(3)
    real(KIND=dp) :: eJ0c(6),eJ1c(6),eJ2c(6)
    real(KIND=dp) :: eJ0t(6),eJ1t(6),eJ2t(6)

    eJ0c(:) = 0
    eJ1c(:) = 0
    eJ2c(:) = 0
    eJ0t(:) = 0
    eJ1t(:) = 0
    eJ2t(:) = 0

    if(allocated(Den%JmuNu)) then

      !----------- recouple cartesian tensor to spherical tensors of rank 0,1,2
      ! note that the J2(:,:,:,2,1,it), J2(:,:,:,3,1,it), and J2(:,:,:,3,2,it)
      ! components of the pseudotensor are allocated but not calculated as the
      ! pseudotensor is symmetric by construction.
      !------------------------------------------------------------------------
      do it=1,2
        J0(:,:,:,it) = Den%JMuNu(:,:,:,1,1,it) &
        &             +Den%JMuNu(:,:,:,2,2,it) &
        &             +Den%JMuNu(:,:,:,3,3,it)
      enddo
      do it=1,2
        J1(:,:,:,1,it) = Den%JMuNu(:,:,:,2,3,it) - Den%JMuNu(:,:,:,3,2,it)
        J1(:,:,:,2,it) = Den%JMuNu(:,:,:,3,1,it) - Den%JMuNu(:,:,:,1,3,it)
        J1(:,:,:,3,it) = Den%JMuNu(:,:,:,1,2,it) - Den%JMuNu(:,:,:,2,1,it)
      enddo
      do it=1,2
        J2(:,:,:,1,1,it) = Den%JMuNu(:,:,:,1,1,it) - 1/3.0_dp * J0(:,:,:,it)
        J2(:,:,:,2,2,it) = Den%JMuNu(:,:,:,2,2,it) - 1/3.0_dp * J0(:,:,:,it)
        J2(:,:,:,3,3,it) = Den%JMuNu(:,:,:,3,3,it) - 1/3.0_dp * J0(:,:,:,it)
        J2(:,:,:,1,2,it) = 0.5_dp * ( Den%JMuNu(:,:,:,1,2,it) &
        &                            +Den%JMuNu(:,:,:,2,1,it) )
        J2(:,:,:,1,3,it) = 0.5_dp * ( Den%JMuNu(:,:,:,1,3,it) &
        &                            +Den%JMuNu(:,:,:,3,1,it) )
        J2(:,:,:,2,3,it) = 0.5_dp * ( Den%JMuNu(:,:,:,2,3,it) &
        &                            +Den%JMuNu(:,:,:,3,2,it) )
      enddo

      !------------------------------------ energy of the J^{(0)} J^{(0)} term 
      er(1) =          dv * sum(J0(:,:,:,1)*J0(:,:,:,1))  !   J^(0)_n J^(0)_n
      er(2) =          dv * sum(J0(:,:,:,2)*J0(:,:,:,2))  !   J^(0)_p J^(0)_p
      er(3) = 2.0_dp * dv * sum(J0(:,:,:,1)*J0(:,:,:,2))  ! 2 J^(0)_n J^(0)_p

      eJ0c(1) = 1/3.0_dp *(C14+C15       )*(er(1)+er(2))             & ! total
      &        +1/3.0_dp * C14            * er(3)                     
      eJ0c(2) = 1/3.0_dp *(C14+C15*0.5_dp)*(er(1)+er(2)+er(3))         ! T=0
      eJ0c(3) = 1/3.0_dp *(    C15*0.5_dp)*(er(1)+er(2)-er(3))         ! T=1
      eJ0c(4) = 1/3.0_dp *(C14+C15       )* er(1)                      ! nn
      eJ0c(5) = 1/3.0_dp *(C14+C15       )* er(2)                      ! pp
      eJ0c(6) = 1/3.0_dp * C14            * er(3)                      ! np
      eJ0t(1) = 1/3.0_dp *(T14+T15       )*(er(1)+er(2))             & ! total
      &        +1/3.0_dp * T14            * er(3)                    & 
      &        +4/3.0_dp *(B16+B17       )*(er(1)+er(2))             &
      &        +4/3.0_dp * B16            * er(3)
      eJ0t(2) = 1/3.0_dp *(T14+T15*0.5_dp)*(er(1)+er(2)+er(3))       & ! T=0
      &        +4/3.0_dp *(B16+B17*0.5_dp)*(er(1)+er(2)+er(3))
      eJ0t(3) = 1/3.0_dp *(    T15*0.5_dp)*(er(1)+er(2)-er(3))       & ! T=1
      &        +4/3.0_dp *(    B17*0.5_dp)*(er(1)+er(2)-er(3))
      eJ0t(4) = 1/3.0_dp *(T14+T15       )* er(1)                    & ! nn
      &        +4/3.0_dp *(B16+B17       )* er(1)
      eJ0t(5) = 1/3.0_dp *(T14+T15       )* er(2)                    & ! pp
      &        +4/3.0_dp *(B16+B17       )* er(2)
      eJ0t(6) = 1/3.0_dp * T14            * er(3)                    & ! np
      &        +4/3.0_dp * B16            * er(3)

      !------------------------------------ energy of the J^{(1)} J^{(1)} term
      er(1) =          dv * ( sum(J1(:,:,:,1,1)*J1(:,:,:,1,1)) &
      &                      +sum(J1(:,:,:,2,1)*J1(:,:,:,2,1)) &
      &                      +sum(J1(:,:,:,3,1)*J1(:,:,:,3,1)))
      er(2) =          dv * ( sum(J1(:,:,:,1,2)*J1(:,:,:,1,2)) &
      &                      +sum(J1(:,:,:,2,2)*J1(:,:,:,2,2)) &
      &                      +sum(J1(:,:,:,3,2)*J1(:,:,:,3,2)))
      er(3) = 2.0_dp * dv * ( sum(J1(:,:,:,1,1)*J1(:,:,:,1,2)) &
      &                      +sum(J1(:,:,:,2,1)*J1(:,:,:,2,2)) &
      &                      +sum(J1(:,:,:,3,1)*J1(:,:,:,3,2)))

      eJ1c(1) = 0.5_dp *(C14+C15       )*(er(1)+er(2)) &             ! total
      &        +0.5_dp * C14             *er(3)                     
      eJ1c(2) = 0.5_dp *(C14+C15*0.5_dp)*(er(1)+er(2)+er(3))         ! T=0
      eJ1c(3) = 0.5_dp *(    C15*0.5_dp)*(er(1)+er(2)-er(3))         ! T=1
      eJ1c(4) = 0.5_dp *(C14+C15       )* er(1)                      ! nn
      eJ1c(5) = 0.5_dp *(C14+C15       )* er(2)                      ! pp
      eJ1c(6) = 0.5_dp * C14            * er(3)                      ! np
      eJ1t(1) = 0.5_dp *(T14+T15       )*(er(1)+er(2)) &             ! total
      &        +0.5_dp * T14             *er(3)        &            
      &        -0.5_dp *(B16+B17       )*(er(1)+er(2)) &
      &        -0.5_dp * B16            * er(3)
      eJ1t(2) = 0.5_dp *(T14+T15*0.5_dp)*(er(1)+er(2)+er(3))       & ! T=0
      &        -0.5_dp *(B16+B17*0.5_dp)*(er(1)+er(2)+er(3))
      eJ1t(3) = 0.5_dp *(    T15*0.5_dp)*(er(1)+er(2)-er(3))       & ! T=1
      &        -0.5_dp *(    B17*0.5_dp)*(er(1)+er(2)-er(3))
      eJ1t(4) = 0.5_dp *(T14+T15       )* er(1)                    & ! nn
      &        -0.5_dp *(B16+B17       )* er(1)
      eJ1t(5) = 0.5_dp *(T14+T15       )* er(2)                    & ! pp
      &        -0.5_dp *(B16+B17       )* er(2)
      eJ1t(6) = 0.5_dp * T14            * er(3)                    & ! np
      &        -0.5_dp * B16            * er(3)

      !------------------------------------ energy of the J^{(2)} J^{(2)} term 
      er(1) = dv * (          sum(J2(:,:,:,1,1,1)*J2(:,:,:,1,1,1)) &
      &             +2.0_dp * sum(J2(:,:,:,1,2,1)*J2(:,:,:,1,2,1)) &
      &             +2.0_dp * sum(J2(:,:,:,1,3,1)*J2(:,:,:,1,3,1)) &
      &             +         sum(J2(:,:,:,2,2,1)*J2(:,:,:,2,2,1)) &
      &             +2.0_dp * sum(J2(:,:,:,2,3,1)*J2(:,:,:,2,3,1)) &
      &             +         sum(J2(:,:,:,3,3,1)*J2(:,:,:,3,3,1)))
      er(2) = dv * (          sum(J2(:,:,:,1,1,2)*J2(:,:,:,1,1,2)) &
      &             +2.0_dp * sum(J2(:,:,:,1,2,2)*J2(:,:,:,1,2,2)) &
      &             +2.0_dp * sum(J2(:,:,:,1,3,2)*J2(:,:,:,1,3,2)) &
      &             +         sum(J2(:,:,:,2,2,2)*J2(:,:,:,2,2,2)) &
      &             +2.0_dp * sum(J2(:,:,:,2,3,2)*J2(:,:,:,2,3,2)) &
      &             +         sum(J2(:,:,:,3,3,2)*J2(:,:,:,3,3,2)))
      er(3) = dv * ( 2.0_dp * sum(J2(:,:,:,1,1,1)*J2(:,:,:,1,1,2)) &
      &             +4.0_dp * sum(J2(:,:,:,1,2,1)*J2(:,:,:,1,2,2)) &
      &             +4.0_dp * sum(J2(:,:,:,1,3,1)*J2(:,:,:,1,3,2)) &
      &             +2.0_dp * sum(J2(:,:,:,2,2,1)*J2(:,:,:,2,2,2)) &
      &             +4.0_dp * sum(J2(:,:,:,2,3,1)*J2(:,:,:,2,3,2)) &
      &             +2.0_dp * sum(J2(:,:,:,3,3,1)*J2(:,:,:,3,3,2)))
      eJ2c(1) = (C14+C15       )*(er(1)+er(2))             & ! total
      &        + C14            * er(3)
      eJ2c(2) = (C14+C15*0.5_dp)*(er(1)+er(2)+er(3))         ! T=0
      eJ2c(3) = (    C15*0.5_dp)*(er(1)+er(2)-er(3))         ! T=1
      eJ2c(4) = (C14+C15       )* er(1)                      ! nn
      eJ2c(5) = (C14+C15       )* er(2)                      ! pp
      eJ2c(6) =  C14            * er(3)                      ! np
      eJ2t(1) = (T14+T15       )*(er(1)+er(2))             & ! total
      &        + T14            * er(3)                    &
      &        +(B16+B17       )*(er(1)+er(2))             &
      &        + B16            * er(3)
      eJ2t(2) = (T14+T15*0.5_dp)*(er(1)+er(2)+er(3))       & ! T=0
      &        +(B16+B17*0.5_dp)*(er(1)+er(2)+er(3))        
      eJ2t(3) = (    T15*0.5_dp)*(er(1)+er(2)-er(3))       & ! T=1
      &        +(    B17*0.5_dp)*(er(1)+er(2)-er(3))
      eJ2t(4) = (T14+T15       )* er(1)                    & ! nn
      &        +(B16+B17       )* er(1)
      eJ2t(5) = (T14+T15       )* er(2)                    & ! pp
      &        +(B16+B17       )* er(2)
      eJ2t(6) =  T14            * er(3)                    & ! np
      &        + B16            * er(3)

    endif

    !------------------------------------ map on variable to be returned
    Tterms(1,:) = eJ0c(:) + eJ0t(:) 
    Tterms(2,:) = eJ1c(:) + eJ1t(:)
    Tterms(3,:) = eJ2c(:) + eJ2t(:)

    return
  end function recoupleTensorNLO

  subroutine PrintEnergy_Bterms(Lagrange)
    !---------------------------------------------------------------------------
    ! This subroutine prints all relevant info that is contained in this
    ! module plus some extra diagnostics.
    !---------------------------------------------------------------------------
    ! Optional Argument Lagrange determines whether or not the energies have
    ! been calculated with the lagrange derivatives at the end.
    !---------------------------------------------------------------------------

    use Pairing, only : PairingType, Lipkin

    logical, intent(in), optional :: Lagrange

      1 format (22('-'),  ' Energies (MeV) ', 22('-'))
     11 format (18('-'),  ' Lagrange Energies (MeV) ', 17('-'))
      2 format (60('_'))
      3 format (60('-'))
    102 format (1x, 3('B', i2,2x,  f12.5, 2x))
    103 format (' Kinetic Energy ',/,                                          &
        &    2x,'Neutron', f12.5, ' Proton ', f12.5, ' Total ', f12.5)
    104 format (' Pairing Energy ',/,                                          &
        &    2x,'Neutron', f12.5, ' Proton ', f12.5, ' Total ', f12.5)
    105 format (' Lipkin-Nogami Energy ',/,                                    &
        &    2x,'Neutron', f12.5, ' Proton ', f12.5, ' Total ', f12.5)
    106 format (' Coulomb ',/,3x,'Direct', f12.5, ' Exch.  ', f12.5)
    107 format (' Total Energy ',/,                                            &
        & 3x,'from functional: ', f15.6, /,                                    &
        & 3x,'from Spwfs:      ', f15.6)
    108 format (' Differences:',/,                                             &
        & 3x,'Absolute:        ', f15.6,/                                      &
        & 3x,'Relative:        ', 4x,es15.6,/                                  &
        & 3x,'Abs. func vs sp  ', f15.6,/                                      &
        & 3x,'Rel. func vs sp  ', 4x,es15.6,/                                  &
        & 3x,'Density changed  ', 4x,es15.6      )
    109 format (' Total Energy ',/,                                            &
        & 3x,'Lagrange:        ', f15.6,/,                                     &
        & 3x,'Non-Lagrange:    ', f15.6,/,                                     &
        & 3x,'Abs. Difference: ', f15.6,/,                                     &
        & 3x,'Rel. Difference: ', 4x,es15.6)
    110 format(' 1-body COM Correction')
    111 format(' 2-body COM Correction')
    112 format(2x,'Neutron', f12.5, ' Proton ', f12.5, ' Total ', f12.5)
    113 format(' Routhian:', 10x, f15.6)

    !Header
    if(present(Lagrange)) then
      print 11
    else
      print 1
    endif

    print 102,  1, BTerm(1),   2, Bterm(2),   3, BTerm(3)
    print 102,  4, BTerm(4),   5, Bterm(5),   6, BTerm(6)
    print 102,  7, BTerm(7),   8, BTerm(8),   9, Bterm(9)
    print 102, 10, BTerm(10), 11, BTerm(11), 12, Bterm(12)
    print 102, 13, BTerm(13), 14, BTerm(14), 15, BTerm(15)
    print 102, 16, Bterm(16), 17, BTerm(17), 18, BTerm(18)
    print 102, 19, Bterm(19), 20, BTerm(20), 21, BTerm(21)

    print 103, Kinetic(1), Kinetic(2), sum(Kinetic)
    if(PairingType.ne.0) then
      print 104, PairingEnergy, sum(PairingEnergy)
      print *, 'wtf3', pairingenergy
      if(Lipkin) then
        print 105, LNEnergy, sum(LNEnergy)
      endif
    endif
    print 106, CoulombEnergy, CoulombExchange

    if(any(ComCorrection(1,:).ne.0.0_dp)) then
      print 110
      print 112, COMCorrection(1,:), sum(COMCorrection(1,:))
    endif
    if(any(COMCorrection(2,:) .ne. 0.0_dp)) then
      print 111
      print 112, COMCorrection(2,:), sum(COMCorrection(2,:))
    endif

    print 2

    !Don't print the energy differences when dealing with Lagrange reanalysis
    if(.not.present(Lagrange)) then

      print 107, TotalEnergy, SPEnergy
      print 113, Routhian
      if(OldEnergy(1).ne.0.0_dp) then
        print 108, abs(TotalEnergy - OldEnergy(1)),                            &
        &          abs((TotalEnergy - OldEnergy(1))/OldEnergy(1)),             &
        &          abs(TotalEnergy - SpEnergy),                                &
        &          abs((TotalEnergy - SpEnergy)/TotalEnergy),                  &
        &          DensityChange
      endif
    else
      !Note that oldenergy(1) now contains the last total energy that was
      !calculated with non-lagrange derivatives.
      print 109, TotalEnergy, OldEnergy(1), TotalEnergy - OldEnergy(1),        &
      &          abs((TotalEnergy - OldEnergy(1))/TotalEnergy)
      print 113, Routhian
    endif
    print 3

  end subroutine PrintEnergy_Bterms

  subroutine PrintEnergy_TermByTerm(Lagrange)
    !---------------------------------------------------------------------------
    ! This subroutine prints all relevant info that is contained in this
    ! module plus some extra diagnostics.
    !---------------------------------------------------------------------------
    ! Optional Argument Lagrange determines whether or not the energies have
    ! been calculated with the lagrange derivatives at the end.
    !---------------------------------------------------------------------------

    use Pairing, only : PairingType, Lipkin

    logical, intent(in), optional :: Lagrange
    real*8                        :: SKTodd, N2LOODD, dispersion
    integer                       :: i

      1 format (22('-'),  ' Energies (MeV) ', 22('-'))
     11 format (18('-'),  ' Lagrange Energies (MeV) ', 18('-'))
      2 format (60('_'))
      3 format (60('-'))

    100 format ('Skyrme Terms')
     12 format (2x,' rho^2_t     =', f12.5, 5x , ' rho^2_q     = ', f12.5, ' t', f12.5)
     13 format (2x,' rho*tau_t   =', f12.5, 5x , ' rho*tau_q   = ', f12.5, ' t', f12.5)
     14 format (2x,'       j^2_t =', f12.5, 5x , ' j^2_q       = ', f12.5, ' t', f12.5)
     15 format (2x,' rho D rho_t =', f12.5, 5x , ' rho D rho_q = ', f12.5, ' t', f12.5)
     16 format (2x,' rho^(2+a)_t =', f12.5, 5x , ' rho^(2+a)_q = ', f12.5, ' t', f12.5)
     17 format (2x,' rho^(2+b)_t =', f12.5, 5x , ' rho^(2+b)_q = ', f12.5, ' t', f12.5)
     18 format (2x,' rho d J_t   =', f12.5, 5x , ' rho d J_q   = ', f12.5, ' t', f12.5)
     19 format (2x,' j d x s_t   =', f12.5, 5x , ' j d x s_q   = ', f12.5, ' t', f12.5)
     20 format (2x,' s^2_t       =', f12.5, 5x , ' s^2_q       = ', f12.5, ' t', f12.5)
     21 format (2x,' rho^a s^2_t =', f12.5, 5x , ' rho^a s^2_q = ', f12.5, ' t', f12.5)
     22 format (2x,' rho^b s^2_t =', f12.5, 5x , ' rho^b s^2_q = ', f12.5, ' t', f12.5)
     23 format (2x,' JmnJmn_t    =', f12.5, 5x , ' JmnJmn_t    = ', f12.5, ' t', f12.5)
     24 format (2x,' s * T_t     =', f12.5, 5x , ' s * T_q     = ', f12.5, ' t', f12.5)
     25 format (2x,' JmmJmm_t    =', f12.5, 5x , ' JmmJmm_q    = ', f12.5, ' t', f12.5)
     26 format (2x,' JmnJnm_t    =', f12.5, 5x , ' JmnJnm_q    = ', f12.5, ' t', f12.5)
     27 format (2x,' s * F_t     =', f12.5, 5x , ' s * F_q     = ', f12.5, ' t', f12.5)
     28 format (2x,' s D s_t     =', f12.5, 5x , ' s D s_t     = ', f12.5, ' t', f12.5)
     29 format (2x,' (N s)^2_t   =', f12.5, 5x , ' (N s)^2_q   = ', f12.5, ' t', f12.5)

    200 format ('Recoupled Tensor Terms')
    201 format (2x,' J0 * J0_0   =', f12.5, 5x , ' J0 * J0_1   = ', f12.5, ' t', f12.5)
    202 format (2x,' J1 * J1_0   =', f12.5, 5x , ' J1 * J1_1   = ', f12.5, ' t', f12.5)
    203 format (2x,' J2 * J2_0   =', f12.5, 5x , ' J2 * J2_1   = ', f12.5, ' t', f12.5)
    204 format (2x,' J * J_0     =', f12.5, 5x , ' J * J_1     = ', f12.5, ' t', f12.5)
    205 format (2x,'              ',  12x , 5x , '               ',  12x,  ' t', f12.5)
    211 format (2x,' J0 * J0_0e  =',es12.4, 5x , ' J0 * J0_1e  = ',es12.4, ' t', es12.4)
    212 format (2x,' J1 * J1_0e  =',es12.4, 5x , ' J1 * J1_1e  = ',es12.4, ' t', es12.4)
    213 format (2x,' J2 * J2_0e  =',es12.4, 5x , ' J2 * J2_1e  = ',es12.4, ' t', es12.4)
    214 format (2x,' J * J_0e    =',es12.4, 5x , ' J * J_1e    = ',es12.4, ' t', es12.4)
    215 format (2x,'              ',  12x , 5x , '               ',  12x,  ' t', es12.4)

    300 format ('N2LO Terms')
     31 format (2x,' DrhoDrho_t  =', f12.5, 5x ,' DrhoDrho_q  = ', f12.5, ' t', f12.5)
     32 format (2x,' rhoQ_t      =', f12.5, 5x ,' rhoQ_q      = ', f12.5, ' t', f12.5)
     33 format (2x,' tau^2_t     =', f12.5, 5x ,' tau^2_t     = ', f12.5, ' t', f12.5)
     34 format (2x,' Re Tau_mn_t =', f12.5, 5x ,' Re Tau_mn_q = ', f12.5, ' t', f12.5)
     35 format (2x,' Tau_mn rho_t=', f12.5, 5x ,' Tau_mn rho_q= ', f12.5, ' t', f12.5)
     36 format (2x,' Im Tau_mn_t =', f12.5, 5x ,' Im Tau_mn_q = ', f12.5, ' t', f12.5)
     37 format (2x,' DmuJmunu_t  =', f12.5, 5x ,' DmuJmunu_q  = ', f12.5, ' t', f12.5)
     38 format (2x,' JVmunu_t    =', f12.5, 5x ,' JVmunu_q    = ', f12.5, ' t', f12.5)
     39 format (2x,' DsDs_t      =', f12.5, 5x ,' DsDs_q      = ', f12.5, ' t', f12.5)
     40 format (2x,' divj^2_t    =', f12.5, 5x ,' divj^2_q    = ', f12.5, ' t', f12.5)
     41 format (2x,' j Pi_t      =', f12.5, 5x ,' j Pi_q      = ', f12.5, ' t', f12.5)
     42 format (2x,' s S_t       =', f12.5, 5x ,' s S_q       = ', f12.5, ' t', f12.5)
     43 format (2x,' T^2_t       =', f12.5, 5x ,' T^2_q       = ', f12.5, ' t', f12.5)
     44 format (2x,' Re T_mn^2_t =', f12.5, 5x ,' Re T_mn^2_q = ', f12.5, ' t', f12.5)
     45 format (2x,' T_mn Ds_t   =', f12.5, 5x ,' T_mn Ds_q   = ', f12.5, ' t', f12.5)
     46 format (2x,' Im T_mn^2_t =', f12.5, 5x ,' Im T_mn^2_q = ', f12.5, ' t', f12.5)
    
    500 format ('N3LO terms')
     51 format (2x,' rhoDDDrho_t =', f12.5, 5x ,' rhoDDDrho_q = ', f12.5, ' t', f12.5)
     52 format (2x,' tauDtau_t   =', f12.5, 5x ,' tauDtau_q   = ', f12.5, ' t', f12.5)
     53 format (2x,' Re tmnDtmn_t=', f12.5, 5x ,' Re tmnDtmn_q= ', f12.5, ' t', f12.5)
     54 format (2x,' Drho DRtmn_t=', f12.5, 5x ,' Drho DRtmn_q= ', f12.5, ' t', f12.5)
     55 format (2x,' t DRtmn_t   =', f12.5, 5x ,' t DRtmn_q   = ', f12.5, ' t', f12.5)
     56 format (2x,' tauDDrho_t  =', f12.5, 5x ,' tauDDrho_q  = ', f12.5, ' t', f12.5)
     
    301 format (2x,' M^(Drho)    =', f12.5) 
    302 format (2x,' M^even(rho) =', f12.5) 
    303 format (2x,' M^even( s ) =', f12.5) 
    304 format (2x,' M^(Ds )     =', f12.5) 
    305 format (2x,' M^odd(rho)  =', f12.5) 
    306 format (2x,' M^odd( s )  =', f12.5) 
    
     58 format (2x,'Time-even    =', f12.5, 5x , ' Time-odd   = ', f12.5 )
     59 format (2x,'Skyrme Total =', f12.5)
    
    310 format (2x,'T-even N2LO  =', f12.5, 5x , 'T-odd N2LO  = ', f12.5 )
    311 format (2x,'N2LO total   =', f12.5)
     
    103 format (' Kinetic Energy ',/,                                          &
        &    2x,'Kin.  N', f12.5, ' Kin.  P', f12.5, ' Total ', f12.5)
    104 format (' Pairing Energy ',/,                                          &
        &    2x,'Pair. N', f12.5, ' Pair. P', f12.5, ' Total ', f12.5)
    105 format (' Lipkin-Nogami Energy ',/,                                    &
        &    2x,'  LN  N', f12.5, '    LN P', f12.5, ' Total ', f12.5)
    106 format (' Coulomb ',/,3x,'Direct', f12.5, ' Exch.  ', f12.5)
    107 format (' Total Energy ',/,                                            &
        & 3x,'from functional: ', f20.11, /,                                   &
        & 3x,'from Spwfs:      ', f20.11)
    108 format (' Differences:',/,                                             &
        & 3x,'Absolute:        ', f20.11,/                                     &
        & 3x,'Relative:        ', 4x,es20.11,/                                 &
        & 3x,'Abs. func vs sp  ', f20.11,/                                     &
        & 3x,'Rel. func vs sp  ', 4x,es20.11,/                                 &
        & 3x,'Density changed  ', 4x,es20.11      )
    109 format (' Total Energy ',/,                                            &
        & 3x,'Lagrange:        ', f20.11,/,                                    &
        & 3x,'Non-Lagrange:    ', f20.11,/,                                    &
        & 3x,'Spwf energy:     ', f20.11,/,                                    &
        & 3x,'Abs. Difference: ', f20.11,/,                                    &
        & 3x,'Rel. Difference: ', 4x,es20.11)
    110 format(' 1-body COM Correction')
    111 format(' 2-body COM Correction')
    112 format(2x,'COM1  N', f12.5, ' COM1  P', f12.5, ' Total ', f12.5)
    113 format(2x,'COM2  N', f12.5, ' COM2  P', f12.5, ' Total ', f12.5)
    114 format(' Routhian:', 10x, f15.6)
    115 format(' Disp. (weighted):', 6x, es20.11)
    !Header
    if(present(Lagrange)) then
      print 11
    else
      print 1
    endif

    print 100
    print 12, SkyrmeTerms(1),  SkyrmeTerms(2),  sum(SkyrmeTerms(1:2))
    print 13, SkyrmeTerms(3),  SkyrmeTerms(4),  sum(SkyrmeTerms(3:4))
    print 14, SkyrmeTerms(5),  SkyrmeTerms(6),  sum(SkyrmeTerms(5:6))
    print 15, SkyrmeTerms(7),  SkyrmeTerms(8),  sum(SkyrmeTerms(7:8))
    print 16, SkyrmeTerms(9),  SkyrmeTerms(10), sum(SkyrmeTerms(9:10))
    print 17, SkyrmeTerms(33), SkyrmeTerms(34), sum(SkyrmeTerms(33:34))
    print 18, SkyrmeTerms(11), SkyrmeTerms(12), sum(SkyrmeTerms(11:12))
    print 19, SkyrmeTerms(13), SkyrmeTerms(14), sum(SkyrmeTerms(13:14))
    print 20, SkyrmeTerms(15), SkyrmeTerms(16), sum(SkyrmeTerms(15:16))
    print 21, SkyrmeTerms(17), SkyrmeTerms(18), sum(SkyrmeTerms(17:18))
    print 22, SkyrmeTerms(35), SkyrmeTerms(36), sum(SkyrmeTerms(35:36))
    print 23, SkyrmeTerms(19), SkyrmeTerms(20), sum(SkyrmeTerms(19:20))
    print 24, SkyrmeTerms(21), SkyrmeTerms(22), sum(SkyrmeTerms(21:22))
    print 25, SkyrmeTerms(23), SkyrmeTerms(24), sum(SkyrmeTerms(23:24))
    print 26, SkyrmeTerms(25), SkyrmeTerms(26), sum(SkyrmeTerms(25:26))
    print 27, SkyrmeTerms(27), SkyrmeTerms(28), sum(SkyrmeTerms(27:28))
    print 28, SkyrmeTerms(29), SkyrmeTerms(30), sum(SkyrmeTerms(29:30))
    print 29, SkyrmeTerms(31), SkyrmeTerms(32), sum(SkyrmeTerms(31:32))
    print *
    
    if(any(Tensorterms.ne.0)) then
      print 200
      print 201,Tensorterms(1,2),Tensorterms(1,3),Tensorterms(1,1)
      print 202,Tensorterms(2,2),Tensorterms(2,3),Tensorterms(2,1)
      print 203,Tensorterms(3,2),Tensorterms(3,3),Tensorterms(3,1)
      print 204,Tensorterms(1,2)+Tensorterms(2,2)+Tensorterms(3,2), &
      &         Tensorterms(1,3)+Tensorterms(2,3)+Tensorterms(3,3), &
      &         Tensorterms(1,1)+Tensorterms(2,1)+Tensorterms(3,1)
   !   print 205,SkyrmeTerms(19)+SkyrmeTerms(20)+SkyrmeTerms(23) &
   !   &        +SkyrmeTerms(24)+SkyrmeTerms(25)+SkyrmeTerms(26)
      print *
      print 211,Tensorterms(1,2),Tensorterms(1,3),Tensorterms(1,1)
      print 212,Tensorterms(2,2),Tensorterms(2,3),Tensorterms(2,1)
      print 213,Tensorterms(3,2),Tensorterms(3,3),Tensorterms(3,1)
      print 214,Tensorterms(1,2)+Tensorterms(2,2)+Tensorterms(3,2), &
      &         Tensorterms(1,3)+Tensorterms(2,3)+Tensorterms(3,3), &
      &         Tensorterms(1,1)+Tensorterms(2,1)+Tensorterms(3,1)
   !   print 215,SkyrmeTerms(19)+SkyrmeTerms(20)+SkyrmeTerms(23) &
   !   &        +SkyrmeTerms(24)+SkyrmeTerms(25)+SkyrmeTerms(26)
       print *
    endif

    if(any(N2LOterms.ne.0)) then
        print 300
        print 31, N2LOterms(1:2)  , N2LOterms(1)+N2LOterms(2)  
        print 32, N2LOterms(3:4)  , N2LOterms(3)+N2LOterms(4)  
        print 33, N2LOterms(5:6)  , N2LOterms(5)+N2LOterms(6)  
        print 34, N2LOterms(7:8)  , N2LOterms(7)+N2LOterms(8)  
        print 35, N2LOterms(9:10) , N2LOterms(9)+N2LOterms(10) 
        print 36, N2LOterms(11:12), N2LOterms(11)+N2LOterms(12)
        print 37, N2LOterms(13:14), N2LOterms(13)+N2LOterms(14)
        print 38, N2LOterms(15:16), N2LOterms(15)+N2LOterms(16)
        print 39, N2LOterms(17:18), N2LOterms(17)+N2LOterms(18)
        print 40, N2LOterms(19:20), N2LOterms(19)+N2LOterms(20)
        print 41, N2LOterms(21:22), N2LOterms(21)+N2LOterms(22)
        print 42, N2LOterms(23:24), N2LOterms(23)+N2LOterms(24)
        print 43, N2LOterms(25:26), N2LOterms(25)+N2LOterms(26)
        print 44, N2LOterms(27:28), N2LOterms(27)+N2LOterms(28)
        print 45, N2LOterms(29:30), N2LOterms(29)+N2LOterms(30)
        print 46, N2LOterms(31:32), N2LOterms(31)+N2LOterms(32)
        print *
        print 301, sum(N2LOterms(1:2))
        print 302, sum(N2LOterms(3:12))
        print 303, sum(N2LOterms(13:16))
        print 304, sum(N2LOterms(17:18))
        print 305, sum(N2LOterms(19:22))
        print 306, sum(N2LOterms(23:32))
        print *
    endif
    
    if(any(N3LOterms.ne.0)) then
        print 500
        print 51, N3LOterms( 1: 2), sum(N3LOterms( 1: 2))
        print 52, N3LOterms( 3: 4), sum(N3LOterms( 3: 4))
        print 53, N3LOterms( 5: 6), sum(N3LOterms( 5: 6))
        print 54, N3LOterms( 7: 8), sum(N3LOterms( 7: 8))
        print 55, N3LOterms( 9:10), sum(N3LOterms( 9:10))
        print 56, N3LOterms(11:12), sum(N3LOterms(11:12))
        print *
    endif
    
    ! Time-even and odd parts
    SkTodd =  SkyrmeTerms(5)   + SkyrmeTerms(6)  + SkyrmeTerms(13) + &
    &         SkyrmeTerms(14)  + SkyrmeTerms(15) + SkyrmeTerms(16) + &
    &         SkyrmeTerms(17)  + SkyrmeTerms(18) + SkyrmeTerms(21) + &
    &         SkyrmeTerms(22)  + SkyrmeTerms(27) + SkyrmeTerms(28) + &
    &         SkyrmeTerms(29)  + SkyrmeTerms(30) + SkyrmeTerms(31) + &
    &         SkyrmeTerms(32)  + SkyrmeTerms(35) + SkyrmeTerms(36)
    
    N2LOODD = N2LOterms(11)    + N2LOterms(12)   + &
    &         N2LOterms(17)    + N2LOterms(18)   + N2LOterms(19)   + & 
    &         N2LOterms(20)    + N2LOterms(21)   + N2LOterms(22)   + &
    &         N2LOterms(23)    + N2LOterms(24)   + N2LOterms(25)   + &
    &         N2LOterms(26)    + N2LOterms(27)   + N2LOterms(28)   + &
    &         N2LOterms(29)    + N2LOterms(30)     

    print 58, sum(SkyrmeTerms) -  SkTodd , SkTodd
    print 59, sum(SkyrmeTerms)

    if(t1n2.ne.0.0) then
         print 310, sum(N2LOterms) - N2LOODD, N2LOODD
         print 311, sum(N2LOterms)
    endif
    
    print *

    print 103, Kinetic(1), Kinetic(2), sum(Kinetic)
    if(PairingType.ne.0) then
      print 104, PairingEnergy, sum(PairingEnergy)
      if(Lipkin) then
        print 105, LNEnergy, sum(LNEnergy)
      endif
    endif
    print 106, CoulombEnergy, CoulombExchange

    if(any(ComCorrection(1,:).ne.0.0_dp)) then
      print 110
      print 112, COMCorrection(1,:), sum(COMCorrection(1,:))
    endif
    !if(any(COMCorrection(2,:) .ne. 0.0_dp)) then
      print 111
      print 113, COMCorrection(2,:), sum(COMCorrection(2,:))
    !endif

    print 2

    !Don't print the energy differences when dealing with Lagrange reanalysis
    if(.not.present(Lagrange)) then

      print 107, TotalEnergy, SPEnergy
      print 114, Routhian
      if(OldEnergy(1).ne.0.0_dp) then
        print 108, abs(TotalEnergy - OldEnergy(1)),                            &
        &          abs((TotalEnergy - OldEnergy(1))/OldEnergy(1)),             &
        &          abs(TotalEnergy - SpEnergy),                                &
        &          abs((TotalEnergy - SpEnergy)/TotalEnergy),                  &
        &          DensityChange
      endif
      
      print *,'from potentials:', SpwfEnergy_func(Upot, Bpot, Wpot)
      
    else
      !Note that oldenergy(1) now contains the last total energy that was
      !calculated with non-lagrange derivatives.
      print 109, TotalEnergy, OldEnergy(1), SPEnergy,         &
      &          TotalEnergy - OldEnergy(1),                  &
      &          abs((TotalEnergy - OldEnergy(1))/TotalEnergy)
      print 114, Routhian
    endif
    !------------------------------------------------------------------------------
    ! Calculate sum of dispersions
    dispersion = 0.0_dp
    do i=1,nwt
        dispersion = dispersion + HFBasis(i)%Occupation*HFBasis(i)%Dispersion
    enddo
    print 115, Dispersion
    
    print 3

  end subroutine PrintEnergy_termbyterm

  subroutine CompKinetic()
    !---------------------------------------------------------------------------
    ! This subroutine computes the total kinetic energy,
    ! according to the following formula:
    !    E_k = -\hbar/2m \int d^3x \sum_{k} v_{k} \Psi_k^* \Delta \Psi_k
    !---------------------------------------------------------------------------
    ! Note that the 1-body c.o.m. correction is not taken into account here!
    ! This is markedly different from EV8, EV4 & CR8!
    !---------------------------------------------------------------------------
    integer          :: l, Isospin, it
    real(KIND=dp)    :: Occupation, Inproduct
    type(Spinor)     :: Value, Lap

    ! Kinetic Energy
    Kinetic = 0.0_dp
    do l=1,nwt
        Isospin    = DensityBasis(l)%GetIsospin()
        Occupation = DensityBasis(l)%GetOcc()

        it = (IsoSpin + 3)/2

        Value = DensityBasis(l)%GetValue()
        Lap   = DensityBasis(l)%GetLap()
        Inproduct = InproductSpinorReal(Value,Lap)

        Kinetic(it)= Kinetic(it) + Occupation*Inproduct
    enddo

    Kinetic=-Kinetic/2.0_dp * hbm
    return
  end subroutine CompKinetic

  function ConverEnergy(EnergyPrec) result(Converged)
    !---------------------------------------------------------------------------
    ! Function that checks if the energy has converged.
    ! For the moment, the criterion is:
    !
    ! |OldEnergy(i) - TotalEnergy|/|TotalEnergy| < EnergyPrec
    !  for the previous 7 iterations.
    !
    !---------------------------------------------------------------------------
    logical                   :: Converged
    integer                   :: i
    real(Kind=dp), intent(in) :: EnergyPrec

    Converged = .true.

    i = 1
    do while(Converged .and. i.le.7)
            Converged = Converged .and.                                        &
            &   ( abs((OldEnergy(i) - TotalEnergy)/TotalEnergy).le.EnergyPrec)
            i = i + 1
    enddo

    return
  end function ConverEnergy

  real(KIND=dp) function SpwfEnergy()
    !---------------------------------------------------------------------------
    ! Computes the energy as a function of the energies of the single particle
    ! states. Note the formula:
    !  E = 1/2 sum_i v2 epsilon_ii + 1/2 E_kin + E_pair
    !      + E_LN - E_constraint - E_rearr
    ! where
    !       - epsilon_i are the Spwf energies
    !       - v2 are the occupations
    !       - E_kin is the kinetic energy
    !       - E_pair is a contribution of the pairing.
    !       - E_LN is due to Lipkin-Nogami
    !       - E_constraint is the energy due to constraints.
    !---------------------------------------------------------------------------

    use Moments, only : ConstraintEnergy,  ConstraintEnergy_J0
    use Cranking,only : Omega, CrankEnergy
    use Spwfstorage, only: TotalAngMom
    integer           :: i, mu

    !Summing the energies of the single particle states.
    SpwfEnergy = 0.0_dp
    do i=1,nwt
      SpwfEnergy=SpwfEnergy+DensityBasis(i)%GetOcc()*DensityBasis(i)%GetEnergy()
    enddo
    ! Adding the kinetic energy and the rearrangement term due to the density
    ! dependence and Coulomb Exchange.
    if(trim(SkyrmeTreatment).eq.'TERMBYTERM') then
        SpwfEnergy = 0.5_dp*      (SpwfEnergy + sum(Kinetic))                       &
        &          - 0.5_dp*byt3a*(sum(SkyrmeTerms( 9:10)) + sum(SkyrmeTerms(17:18)))& 
        &          - 0.5_dp*byt3b*(sum(SkyrmeTerms(33:34)) + sum(SkyrmeTerms(35:36)))& 
        &          + CoulombExchange/3.d0
    elseif(trim(SkyrmeTreatment).eq.'BTERMS') then
        SpwfEnergy = 0.5_dp*(SpwfEnergy + sum(Kinetic))                   &
        &          - 0.5_dp*byt3a*(sum(Bterm(7:8)) + sum(Bterm(12:13)))   &
        &          + CoulombExchange/3.d0
    endif
    !Remove the contribution from the multipole constraints
    SpwfEnergy = SpwfEnergy - sum(ConstraintEnergy*Density%Rho)*dv/2.0_dp
    if(allocated(Density%Jmunu)) then
        do mu=1,3
            SpwfEnergy = SpwfEnergy -  &
            &   sum(ConstraintEnergy_J0*Density%Jmunu(:,:,:,mu,mu,:))*dv/2.0_dp
        enddo
    endif
     !Remove contribution from the cranking constraints
    SpwfEnergy = SpwfEnergy - sum(CrankEnergy)/2.0_dp
    !Add the pairing Energy and Lipkin-Nogami energy
    SpwfEnergy = SpwfEnergy + sum(PairingEnergy) + sum(LNEnergy)
    ! Always add the 1-body COMcorrection. In case it is used iteratively, it
    ! is double counted along with the kinetic energy!
    if(COM1body.gt.0) then
    	SpwfEnergy = SpwfEnergy  + sum(COMCorrection(1,:))/2.0_dp
    endif
    !Add the 2-COMcorrection energy when they are not used iteratively.
    if(COM2body.eq.1) then
    	SpwfEnergy = SpwfEnergy  + sum(COMCorrection(2,:))
    endif
    return
  end function SpwfEnergy
  
  !-----------------------------------------------------------------------------
  real(KIND=dp) function SpwfEnergy_func(U,B,W) result(SE)
    !---------------------------------------------------------------------------
    ! Calculates the total energy with the current densities for a given set
    ! of mean-field potentials.
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: U(nx,ny,nz,2), B(nx,ny,nz,2), W(nx,ny,nz,3,3,2)
    real(KIND=dp) :: exchange, compare, reducedmass(2)
    integer       :: it
   
    
    SE = 0.5* sum(Density%rho * U)*dv 
    SE = SE + 0.5*sum(Density%tau * B)*dv 
    
    Reducedmass= &
    &    (1.0_dp-nucleonmass/(neutrons*nucleonmass(1)+protons*nucleonmass(2)))    
     
    do it=1,2
      SE = SE + 0.25*hbm(it)*Reducedmass(it)*sum(density%tau(:,:,:,it))*dv
    enddo
     
    exchange = 0.5*byt3a * B7a*sum(sum(Density%rho,4)**(2+byt3a))*dv
    do it=1,2
     exchange = exchange + &
     &  0.5*byt3a*B8a*sum(Density%Rho(:,:,:,it)**2*sum(Density%rho,4)**byt3a)*dv
    enddo
    
    if( B7b.ne. 0.0_dp .or. B7b.ne. 0.0_dp ) then
      exchange = exchange + 0.5*byt3b * B7b*sum(sum(Density%rho,4)**(2+byt3b))*dv
      do it=1,2
        exchange = exchange + &
        &  0.5*byt3b*B8b*sum(Density%Rho(:,:,:,it)**2*sum(Density%rho,4)**byt3b)*dv
      enddo
    endif

    ! We are shamelessly taking the Coulomb exchange from last time, supposing
    ! it doesn't change.
    compare = Totalenergy  - CoulombExchange/3.d0 + exchange
    call CompKinetic()
    compare = compare - 0.5*sum(Kinetic) 

    if(COM1body.gt.0) then
      call CompComCorrection()
    	compare = compare  - sum(COMCorrection(1,:))/2.0
    endif
    
    SE = Compare - SE
    
  end function SpwfEnergy_func 
  !-----------------------------------------------------------------------------

  subroutine CompCOMCorrection()
    !---------------------------------------------------------------------------
    ! Subroutine that calculates the COMcorrection diagnostically, meaning
    ! NOT INCLUDED in the single-particle Hamiltonian.
    !---------------------------------------------------------------------------
    ! For the equations and some more discussion, see pg 455 and following in
    ! Ring & Schuck.
    ! A more clear explanation (and some relevant physics) can be found in
    ! M. Bender et al., Eur. Phys. J. A 7, 467-478 (2000)
    !
    !---------------------------------------------------------------------------
    use Force, only       : COM1Body, COM2Body, hbm
    use Spwfstorage, only : HFBasis, NablaMElements
    integer              :: i,j,it, m, actuali, actualj

    COMCorrection = 0.0_dp
    if(COM1Body .gt. 0) then
      !Deduce 1-body COM correction from the Kinetic Energy
      COMCorrection(1,:) = - Kinetic(:) * nucleonmass/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
    endif

    if(COM2Body .gt. 0) then
      !
      call compNablaMelements
      ! We sum carelessly over all single-particle wavefunctions, since the
      ! ones forbidden by symmetry are calculated as zero in the Spwfstorage
      ! module.
      do i=1,nwt
        it = (DensityBasis(i)%GetIsospin() + 3)/2
        do j=1,nwt
            if(it.ne.(DensityBasis(j)%GetIsospin()+3)/2) cycle
            COMCorrection(2,it) = COMCorrection(2,it) +                        &
            &             DensityBasis(i)%GetOcc() * DensityBasis(j)%GetOcc()  &
            &           *(sum(NablaMElements(i,j,:,:)**2))
        enddo
      enddo
      !Some more constants
      COMCorrection(2,:) = COMCorrection(2,:) * hbm * nucleonmass/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))

      if(TRC) then
        ! Take out the extra factor 4 due to the occupation factors being double
        ! what they should be.
        COMCorrection(2,:) = 0.25 * COMCorrection(2,:)
      else
        ! Take out the extra factor 2 due to double-counting everything.
        ! Note that this factor is completely different from the ones above,
        ! it is not half of the 4 in the T-conserved case.
        COMCorrection(2,:) = 0.5 * COMCorrection(2,:)
      endif
    endif
  end subroutine CompCOMCorrection

end module Energy
