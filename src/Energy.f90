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

  implicit none

  save

  ! All the different Energies
  real(KIND=dp),public :: Kinetic(2), BTerm(21), CoulombEnergy
  real(KIND=dp),public :: CoulombExchange,CoMCorrection(2,2)
  real(KIND=dp),public :: TotalEnergy, OldEnergy(7), SpEnergy, PairingEnergy(2)
  real(KIND=dp),public :: LNEnergy(2), Routhian
  
contains
    
  subroutine CompEnergy
    !---------------------------------------------------------------------------
    ! Subroutine that computes all the different energies of the main program 
    ! state.
    !---------------------------------------------------------------------------
    ! Note that the TotalEnergy would give a slightly different result, due to
    ! the different computation of the kinetic energy
    !---------------------------------------------------------------------------
    use Coulomb, only : CompCoulombEnergy,CompCoulombExchange
    use Pairing, only : CompPairingEnergy, Fermi, Delta,LNLambda,PairingDisp
    use Cranking, only: CrankEnergy
    use Moments, only : ConstraintEnergy
    integer :: i
    
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
    
    !Calculate the Skyrme Terms
    Bterm  = SkyrmeTerms(Density)
    
    !Pairing Energy  
    PairingEnergy = CompPairingEnergy(Delta)
    !Lipkin-Nogami Energy
    LNEnergy      = - LNLambda * PairingDisp 

    !Calculating the CoulombEnergy
    CoulombEnergy = CompCoulombEnergy(Density)
    CoulombExchange = CompCoulombExchange(Density)

    !COM Correction
    call CompCOMCorrection
    
    !Sum
    TotalEnergy = sum(Kinetic) + sum(Bterm) + CoulombEnergy +                  &
    &             sum(CoMCorrection) + CoulombExchange + sum(PairingEnergy) +  &
    &             sum(LNEnergy)
    
    !Calculate the Routhian too 
    Routhian = TotalEnergy                                                     &
    !                         Contribution of multipole moment constraints
    &                      +    sum(ConstraintEnergy*Density%Rho)*dv           &
    !                         Contribution of the cranking constraints (Omega*J)
    &                      +    sum(CrankEnergy)
    
    !Energy due to the single particle states.
    SpEnergy = SpwfEnergy()
    
    return 
  end subroutine CompEnergy
  
  function SkyrmeTerms(Den) result(B)
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
    
    RhoT = sum(Den%Rho,4); TauT = sum(Den%Tau,4)
    NablaJT = sum(Den%NablaJ,4) ; LapRhoT = sum(Den%LapRho,4)
    if(allocated(Den%JmuNu)) JMuNuT = sum(Den%JMuNu,6)
    
    if(.not.TRC) then
      VecJT=sum(Den%VecJ,5) ;  VecST=sum(Den%VecS,5) ; RotST=sum(Den%RotS,5) 
      if(allocated(Den%LapS)) then
        VecTT = sum(Den%VecT,5)
        LapST = sum(Den%LapS,5) ;  DivST = sum(Den%DivS,4)
      endif
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
        B(14) = B(14)  - sum(VecST(:,:,:,1)*VecTT(:,:,:,1)+                   &
        & VecST(:,:,:,2)*VecTT(:,:,:,2)+VecST(:,:,:,3)*VecTT(:,:,:,3))
      endif
      B(14) = B14*dv*B(14)
    endif

    !B15 Terms
    if(B15.ne. 0.0_dp) then
      B(15)= sum(Den%JMuNu(:,:,:,:,:,:)**2)
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

      print *, 'Pseudoscalar'
      print *, B(16) * B16 * dv
      ref = B(16) * B16 * dv

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

      print *, 'Tensor'
      print *, B(16) - ref
    endif

    !B17 Terms
    if(B17.ne.0.0_dp) then
      B(17)= 0.0_dp  
      ref =  0.0_dp    
      do it=1,2  
        !Sum_{mu}J_{mumu}^2 terms       
        B(17) = B(17) + sum((Den%JMuNu(:,:,:,1,1,it)+Den%JMuNu(:,:,:,2,2,it)   &
        &             + Den%JMuNu(:,:,:,3,3,it))**2)

        ref = ref +     sum((Den%JMuNu(:,:,:,1,1,it)+Den%JMuNu(:,:,:,2,2,it)   &
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
      ref = B17*dv*ref
      B(17)= B17*dv*B(17)

      print *, 'Pseudoscalar'
      print *, ref
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
  end function SkyrmeTerms
  
  subroutine PrintEnergy(Lagrange)
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
     11 format (18('-'),  ' Lagrange Energies (MeV) ', 18('-')) 
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

  end subroutine PrintEnergy
  
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
        !print *, l, occupation, Inproduct, -occupation*Inproduct/2.0_dp * hbm(1)
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

    use Moments, only : ConstraintEnergy
    use Cranking,only : Omega, CrankEnergy
    use Spwfstorage, only: TotalAngMom
    integer           :: i
    
    !Summing the energies of the single particle states.
    SpwfEnergy = 0.0_dp
    do i=1,nwt
      SpwfEnergy=SpwfEnergy+DensityBasis(i)%GetOcc()*DensityBasis(i)%GetEnergy()
    enddo
    ! Adding the kinetic energy and the rearrangement term due to the density
    ! dependence and Coulomb Exchange
    SpwfEnergy = 0.5_dp*(SpwfEnergy + sum(Kinetic)-byt3a*(Bterm(7) + BTerm(8)))&
    &            + CoulombExchange/3.d0
     !Remove the contribution from the multipole constraints
    SpwfEnergy = SpwfEnergy - sum(ConstraintEnergy*Density%Rho)*dv/2.0_dp
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
      ! We sum carelessly over all single-particle wavefunctions, since the 
      ! ones forbidden by symmetry are calculated as zero in the Spwfstorage
      ! module.
      do i=1,size(NablaMElements,1)
        actuali = mod(i,nwt) + i/nwt
        it = (HFBasis(actuali)%GetIsospin() + 3)/2
        do j=1,size(NablaMElements,2)
          if(i.eq.j) cycle
          actualj = mod(j,nwt) + j/nwt
          do m=1,3
            COMCorrection(2,it) = COMCorrection(2,it)                          &
            &            + sum(NablaMElements(i,j,m,:)*NablaMElements(j,i,m,:))&
            &            *HfBasis(actuali)%GetOcc() * HFBasis(actualj)%GetOcc()
          enddo
        enddo
      enddo
      !Some more constants
      COMCorrection(2,:) = COMCorrection(2,:) * nucleonmass/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
    endif
  end subroutine CompCOMCorrection

end module Energy
