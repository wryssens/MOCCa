module MeanFields
  !-----------------------------------------------------------------------------
  ! Module containing all the mean-fields generated by the Skyrme interaction,
  ! and the ways to compute them and their actions on the Spwfs.
  !-----------------------------------------------------------------------------
  ! In terms of a parametrised Skyrme force the single particle Hamiltonian
  !  looks as follows
  ! h_q = - \nabla B_q(r) \nabla + U_q(r) + \vec{S}_q(r) \dot \sigma 
  !       - i/2 \sum_{mu,nu = x}^{z} [W_{q;mu,nu}(r)\nabla_{mu}\sigma_{nu} 
  !        + \nabla_{mu}\sigma_{nu} W_{q;mu,nu}(r)]
  !       - i/2  [\vec{A}_q(r) \cdot \nabla +\nabla \cdot \vec{A}_q(r)]
  !       - nabla \cdot [\sigma \cdot C_q(r) ] \cdot \nabla 
  !       - \nabla \cdot \vec{D}_q(r) \sigma \cdot \nabla
  !
  ! where q is an isospin index and the formula is taken from 
  !       V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
  ! The formulas for the potentials can be found in this article, or in 
  !       P. Bonche et al., Nucl. Phys. A467 (1987), 115-135
  !-----------------------------------------------------------------------------
  ! Notes:
  ! - Apparently there is some reason (we haven't found yet) so that 
  !    (\div \vec{A}) and (Div W)
  !   Do not contribute to the single particle energies, meaning that
  !   \int d\vec{r} \psi_m^* (\div \vec{A}) \psi_n = 0 for all n,m
  !   \int d\vec{r} \psi_m^* \sum_{\mu}[\nabla_{\mu}W_{\mu,\nu}] \psi_n = 0 
  !           for all n,m
  !This does not mean that these potentials are zero everywhere, but only that
  !they do not contribute to the single particle energies. It is for this 
  !reason that they are not calculated.
  !-----------------------------------------------------------------------------
  
  use CompilationInfo
  use Geninfo
  use Densities, only : Density 
  use Wavefunctions
  use Force
  
  implicit none
  
  public
  
  !-----------------------------------------------------------------------------
  !Mean Field Potentials (and some of their relevant derivatives)
  real(KIND=dp),allocatable :: BPot(:,:,:,:), NablaBPot(:,:,:,:,:)
  real(KIND=dp),allocatable :: UPot(:,:,:,:),APot(:,:,:,:,:)
  real(KIND=dp),allocatable :: SPot(:,:,:,:,:)
  real(KIND=dp),allocatable :: Cpot(:,:,:,:,:),WPot(:,:,:,:,:,:)
  real(KIND=dp),allocatable :: DPot(:,:,:,:,:)
  real(KIND=dp),allocatable :: DerCPot(:,:,:,:,:,:),DivDPot(:,:,:,:)
  !-----------------------------------------------------------------------------
  !For experimentally calculating B in a more logical way.
  ! And including an astract interface for the PGI compilers (sigh....)
  abstract interface
    function B(Psi,NoKinetic)
      import
      type(Spwf), intent(in) :: Psi
      type(Spinor)           :: B
      logical, intent(in), optional :: NoKinetic
    end function B
  end interface

  procedure(B),pointer :: ActionOfB
  
contains
    
  !=============================================================================
  ! Computing all the potentials.
  !=============================================================================
  subroutine ConstructPotentials
    !---------------------------------------------------------------------------
    ! This subroutine constructs the following potentials for use in the 
    ! single-particle Hamiltonian.
    !  - B, U, S, W, A, C, D and Coulomb 
    !---------------------------------------------------------------------------                        
    
    use Coulomb, only: SolveCoulomb
    !Allocate the potentials if this has not already happened
    if(.not.allocated(UPot)) then
      allocate(BPot(nx,ny,nz,2), NablaBPot(nx,ny,nz,3,2),UPot(nx,ny,nz,2))
      allocate(SPot(nx,ny,nz,3,2),APot(nx,ny,nz,3,2))
      allocate(Cpot(nx,ny,nz,3,2),WPot(nx,ny,nz,3,3,2),DPot(nx,ny,nz,3,2))
      allocate(DerCPot(nx,ny,nz,3,3,2),DivDPot(nx,ny,nz,2))

      BPot     =0.0_dp
      NablaBPot=0.0_dp
      UPot     =0.0_dp
      APot     =0.0_dp
      SPot     =0.0_dp
      CPot     =0.0_dp
      Wpot     =0.0_dp
      DPot     =0.0_dp
      DerCPot  =0.0_dp
      DivDPot  =0.0_dp
    endif
    
    !Done here for the ifort compiler
    ActionOfB=> ActionOfBOld
    
    call CalcBPot()
    !Calculate the Coulomb Potential
    call SolveCoulomb()    
    
    !Construct the mean-field potentials
    call CalcUPot() !=> This includes the contribution of constraints & Coulomb
    call CalcWPot()

    if(.not.TRC) then
        call CalcSPot()                             
        call CalcAPot()
        call CalcCPot()
        call CalcDPot()
    endif             
  end subroutine ConstructPotentials
  
  subroutine DerivePotentials
  !-----------------------------------------------------------------------------
  ! Calculates the derivatives of several potentials we will need.
  ! 
  ! Gradient of B => NablaBPot
  ! All derivatives of C => DerCPot
  ! Divergence of D => DivDPot
  !-----------------------------------------------------------------------------
  ! Note: separate routine in case any one wants to implement mixing of 
  ! potentials.
  !-----------------------------------------------------------------------------
    use Derivatives
    
    integer :: it
    real(KIND=dp):: DerXDX(nx,ny,nz,2),DerYDY(nx,ny,nz,2), DerZDZ(nx,ny,nz,2) 
    
    do it=1,2
      !Gradient of B
      NablaBPot(:,:,:,1,it) = &
      & DeriveX(BPot(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt, 1)
      NablaBPot(:,:,:,2,it) = &
      & DeriveY(BPot(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt, 1)
      NablaBPot(:,:,:,3,it) = &
      & DeriveZ(BPot(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt, 1)
    enddo

    !Calculating the divergence of the Dpotential
    do it=1,2
        DerXDX(:,:,:,it) = &
        & DeriveX(DPot(:,:,:,1,it), ParityInt, SignatureInt,TimeSimplexInt, 1)
        DerYDY(:,:,:,it) = &
        & DeriveY(DPot(:,:,:,2,it), ParityInt, SignatureInt,TimeSimplexInt, 1)
        DerZDZ(:,:,:,it) = &
        & DeriveZ(DPot(:,:,:,3,it), ParityInt, SignatureInt,TimeSimplexInt, 1)
    enddo
    DivDPot = DerXDX + DerYDY + DerZDZ

    !The derivative of the cpotential
    do it=1,2                 
        DerCPot(:,:,:,1,1,it) = &
        & DeriveX(CPot(:,:,:,1,it), ParityInt,-SignatureInt,TimeSimplexInt,1)
        DerCPot(:,:,:,1,2,it) = &
        & DeriveY(CPot(:,:,:,1,it), ParityInt,-SignatureInt,TimeSimplexInt,1)
        DerCPot(:,:,:,1,3,it) = &
        & DeriveZ(CPot(:,:,:,1,it), ParityInt,-SignatureInt,TimeSimplexInt,1)
    
        DerCPot(:,:,:,2,1,it) = &
        & DeriveX(CPot(:,:,:,2,it), ParityInt,-SignatureInt,TimeSimplexInt,2)
        DerCPot(:,:,:,2,2,it) = &
        & DeriveY(CPot(:,:,:,2,it), ParityInt,-SignatureInt,TimeSimplexInt,2)
        DerCPot(:,:,:,2,3,it) = &
        & DeriveZ(CPot(:,:,:,2,it), ParityInt,-SignatureInt,TimeSimplexInt,2)
    
        DerCPot(:,:,:,3,1,it) = &
        & DeriveX(CPot(:,:,:,3,it), ParityInt, SignatureInt,TimeSimplexInt,1)
        DerCPot(:,:,:,3,2,it) = &
        & DeriveY(CPot(:,:,:,3,it), ParityInt, SignatureInt,TimeSimplexInt,1)
        DerCPot(:,:,:,3,3,it) = &
        & DeriveZ(CPot(:,:,:,3,it), ParityInt, SignatureInt,TimeSimplexInt,1)
    enddo
  end subroutine DerivePotentials
  !=============================================================================
  ! Subroutines for calculating all the different potentials.
  !============================================================================= 
  subroutine CalcBPot()
  !-----------------------------------------------------------------------------
  ! This subroutine computes the modified mass B_q(r).
  !       hbar^2/(2m_q^*) = hbar^2/(2m_q) + B_3 \rho + B_4 \rho_q
  !
  !-----------------------------------------------------------------------------
    integer :: it, at, i
    real(KIND=dp) :: Reducedmass(2)

    Reducedmass=(1.0_dp-nucleonmass/(neutrons*nucleonmass(1)+protons*nucleonmass(2)))

    do it=1,2
        at = 3 - it
        BPot(:,:,:,it) =B3*Density%Rho(:,:,:,at) + (B3+B4)*Density%Rho(:,:,:,it) 
        
        if(COM1body .eq. 2) then
          !Include 1-body C.O.M. correction
          Bpot(:,:,:,it) = Bpot(:,:,:,it) + hbm(it)/2.0_dp*Reducedmass(it)
        else
          !Don't include 1-body C.O.M. correction
          Bpot(:,:,:,it) = Bpot(:,:,:,it) + hbm(it)/2.0_dp
        endif

        do i=1,3
          NablaBPot(:,:,:,i,it) = B3 * Density%DerRho(:,:,:,i,at) +            &
          &                       (B3+B4)*Density%DerRho(:,:,:,i,it)
        enddo
    enddo
    return
  end subroutine CalcBPot
                
  subroutine CalcUPot()
  !-----------------------------------------------------------------------------
  ! This subroutine calculates the single particle potential U.
  ! U_q = 2 B1 rho + 2 B2 rho_q
  !       + B3(tau + i \nabla \cdot \vec{j} ) 
  !       + B4 (tau_q + i \nabla \cdot \vec{j}_q )
  !       + 2 B5 \Delta \rho + 2 B6 \Delta \rho_q
  !       + B7 (2 + \byt3a) \rho^{1+\byt3a}
  !       + B8 (\byt3a \rho^{byt3a - 1}(\rho_q^2 + \rho_q^2) 
  !             + 2\rho^{\byt3a} \rho_q)
  !       + B9 (\nabla \cdot \vec{J} + \nabla \cdot \vec{J}_q)
  !       + \byt3a \rho^{byt3a - 1} 
  !          * (B12a \vec{s}^2 + B_13(\vec{s}_q^2 + \vec{s}_{-q}^2 ))
  !       + Multipole Constraint Contributions
  !
  ! For protons, the coulomb potential is also added to this term:
  ! U_{coulomb} = V_{coulomb, direct} - V_{Coulomb, Exchange}
  ! For the moment, we use the Slater approximation for the exchange part:
  !                       V_{exchange} = e(3*\rho_p/\pi)^{1/3}
  !
  !-----------------------------------------------------------------------------
    use Moments, only : ConstraintEnergy         
    use Coulomb, only : CoulombPotential, CoulExchange
    integer        :: it, at
    
    !Temporary storage for total densities.
    real(KIND=dp)  :: RhoTot(nx,ny,nz),VecSTot(nx,ny,nz,3), NablaJTot(nx,ny,nz)
    real(KIND=dp)  :: eps
    
    RhoTot = sum(Density%Rho,4)
    if(.not.TRC) then
      VecsTot = sum(Density%VecS,5)
    endif
    NablaJTot = sum(Density%NablaJ, 4)

    ! This parameter is present to stabilize the calculation of Rho**(byt3a - 1)
    ! Since usually byt3a < 1, this becomes extremely unstable for small values 
    ! of rho. Thus we calculate it as Rho**byt3a/(Rho + eps)
    eps = 1.d-20
    
    !B1 & B2                     
    UPot(:,:,:,1)= (2.0_dp*B1 + 2.0_dp*B2)*Density%Rho(:,:,:,1)                &
    &            + 2.0_dp*B1*Density%Rho(:,:,:,2) 
    UPot(:,:,:,2)= (2.0_dp*B1 + 2.0_dp*B2)*Density%Rho(:,:,:,2)                &
    &            + 2.0_dp*B1*Density%Rho(:,:,:,1)

    do it=1,2
            at = 3 - it
            UPot(:,:,:,it) = UPot(:,:,:,it)                              &
            & + (B3+B4)*Density%Tau(:,:,:,it) + B3*Density%Tau(:,:,:,at) & 
            & + 2.0_dp*((B5 + B6)*Density%LapRho(:,:,:,it)               &
            &          + B5      *Density%LapRho(:,:,:,at))              & 
            & + (2+byt3a)*B7a*(RhoTot+eps)**(1+byt3a)                 & 
            & + B8a*(                                                 &
            & byt3a*(RhoTot**(byt3a))/(RhoTot+eps)*                   &
            & (Density%Rho(:,:,:,it)**2 + Density%Rho(:,:,:,at)**2)   &
            & + 2.0_dp*(RhoTot)**(byt3a)*Density%Rho(:,:,:,it)        &
            & )                                                       &
            & + B9*NablaJTot + B9q*Density%NablaJ(:,:,:,it)            

            if(.not.TRC) then
              Upot(:,:,:,it) = Upot(:,:,:,it)                         &
              & + byt3a*RhoTot**(byt3a)/(RhoTot + eps)*(                &
              &   B12a*sum(VecSTot**2,4)                                & 
              & + B13a*(sum(Density%VecS(:,:,:,:,it)**2,4)              &
              & +       sum(Density%VecS(:,:,:,:,at)**2,4)))
            endif   
    enddo
    ! Note that CoulExchange already carries a minus sign.
    UPot(:,:,:,2)= UPot(:,:,:,2) + CoulombPotential(:,:,:) + CoulExchange(:,:,:)
    
    !Add the contribution from Constraints on the Multipole Moment
    UPot = UPot + ConstraintEnergy
    return                
  end subroutine CalcUPot
                
  subroutine CalcAPot()
  !-----------------------------------------------------------------------------
  ! This subroutine calculates the A potential. 
  ! \vec{A}_q = - B3*2*vec{j} + B4*2*vec{j}_q
  !             + B9*(\nabla x \vec{s} + \nabla x \vec{s}_q
  !-----------------------------------------------------------------------------
  ! Notice that A in the tensor paper corresponds to - C in the Mg24 paper!
  ! This means that A corresponds to the pmx/y/z potential in the CR8 code, 
  ! without minus sign!
  !-----------------------------------------------------------------------------
    use Cranking, only : CrankApot
    
    integer :: it,at,i,j,k

    APot=0.0_dp

    do it=1,2
      at = 3 - it
      !B3 contribution
      APot(:,:,:,:,it) = - 2.0_dp*((B3+B4)*Density%Vecj(:,:,:,:,it) +          &
      &                                 B3*Density%VecJ(:,:,:,:,at))
      !B9 contribution
      APot(:,:,:,:,it) = APot(:,:,:,:,it)                                      &
      &                + ((B9+B9q)*Density%RotS(:,:,:,:,it)                    &
      &                +      B9  *Density%RotS(:,:,:,:,at))
    enddo
    !Add a cranking contribution
    APot = APot + CrankAPot()
    return
  end subroutine CalcAPot
                
  subroutine CalcWPot()
  !-----------------------------------------------------------------------------
  ! This subroutine calculates the W potential.
  !       \vec{W}_q = 
  !             -  \sum_{k=x}^z epsilon_{k,mu,nu} B9 \partial_k (\rho + \rho_q)
  !             + B14*2* J_{mu nu} + B15*2*J_{q,mu nu}
  !             + B16*2*(J_{nu mu} + \sum_{k=x}^z J_kk \delta_{mu,nu})
  !             + B17*2*(J_{q;nu mu} + \sum_{k=x}^z J_{q;kk} \delta_{mu,nu})
  !-----------------------------------------------------------------------------
    integer :: m,n,k,it

    WPot=0.0_dp
    do it=1,2                       
      do m=1,3
         do n=1,3
           !B9 Contribution
           do k=1,3 
             WPot(:,:,:,m,n,it) = WPot(:,:,:,m,n,it) - LeviCivita(k,m,n)*      &
             &        ((B9 + B9q)*Density%DerRho(:,:,:,k,it)                   &
             &        + B9*Density%DerRho(:,:,:,k,3-it))
           enddo        
         enddo
      enddo                               
    enddo 

    if(B14.ne.0.0_dp .or. B15.ne.0.0_dp) then
        do it=1,2                       
          !B14&15 Contribution
           WPot(:,:,:,:,:,it) =   WPot(:,:,:,:,:,it   )                &
           & + 2.0_dp *(B14+B15)*Density%Jmunu(:,:,:,:,:,it   )        &
           & + 2.0_dp * B14     *Density%Jmunu(:,:,:,:,:,3-it )
        enddo
   endif
   if(B16 .ne. 0.0_dp .or. B17 .ne. 0.0_dp) then
      do it=1,2                       
        do m=1,3
           do n=1,3
             !B16 & B17 Contribution 
             ! a) J_{nu mu}
             WPot(:,:,:,m,n,it) = WPot(:,:,:,m,n,it)                   &
             & + 2.0_dp*((B16+B17)*Density%Jmunu(:,:,:,n,m,it)         &
             & + B16*Density%Jmunu(:,:,:,n,m,3-it))
             !b) \delta_{mu,nu}
             if(m.eq.n) then
                do k=1,3
                  WPot(:,:,:,m,n,it) = WPot(:,:,:,m,n,it)              &
                  & + (B16+B17)*Density%Jmunu(:,:,:,k,k,it)            &
                  & + B16*Density%Jmunu(:,:,:,k,k,3-it)
                enddo
             endif
           enddo
        enddo                               
      enddo 
    endif        

    return
  end subroutine CalcWPot
                
  subroutine CalcSPot()
    !---------------------------------------------------------------------------
    ! This subroutine calculates the S potential (V potential in 24Mg paper).
    !
    ! \vec{S}_q =   B9 \nabla x \vec{j} + B9_q \nabla x \vec{j}_q 
    !             + B10*2*\vec{s} + B11*2*\vec{s}_q
    !             + B12a*2*rho^byt3a \vec{s} + B13a*2*rho^byt3a \vec{s}_q
    !             - B14*2*\vec{T} - B15*2*\vec{T}_q
    !             - B16*2*\vec{F} - B17*2*\vec{F}_q
    !             + B18*2*\Delta \vec{s} + B19*2*\Delta \vec{s}_q
    !             - B20*2*\nabla(\nabla \cdot \vec{s}) 
    !             - B21*2*\nabla(\nabla \cdot \vec{s})
    !
    ! Notice that I think the sign in Hellemans et al. is wrong for the B9 
    ! terms here.
    !---------------------------------------------------------------------------
    use Cranking, only : CrankSPot
   
    integer :: it,at,l,i,j,k
    real(KIND=dp)  :: RhoT(nx,ny,nz), RhoVecS(nx,ny,nz,3,2)
    
    Spot = 0.0_dp
    RhoT = Density%Rho(:,:,:,1) + Density%Rho(:,:,:,2)
    
    do it=1,2
      do l=1,3
        RhoVecS(:,:,:,l,it) = RhoT**(byt3a)*Density%VecS(:,:,:,l,it)
      enddo
    enddo
    do it=1,2
      at = 3 - it
      Spot(:,:,:,:,it) = &
       &          + (B9 + B9q)         *Density%RotVecJ(:,:,:,:,it)   &
       &          +  B9                *Density%RotVecJ(:,:,:,:,at)   &
       &          + 2.0_dp*((B10+B11)  *Density%Vecs(:,:,:,:,it)      &
       &          + B10                *Density%VecS(:,:,:,:,at))     &
       &          + 2.0_dp*((B12a+B13a)*        RhoVecS(:,:,:,:,it)   &
       &          + B12a               *        RhoVecS(:,:,:,:,at))
       
       if(B14.ne.0.0_dp .or. B15.ne.0.0_dp) then
        SPot(:,:,:,:,it) = SPot(:,:,:,:,it) &           
         &          -        ((B14+B15)  *Density%VecT(:,:,:,:,it)     &
         &          - B14                *Density%VecT(:,:,:,:,at))            
       endif
       if(B16.ne.0.0_dp .or. B17.ne.0.0_dp) then  
         SPot(:,:,:,:,it) = SPot(:,:,:,:,it) & 
         &          - 2.0_dp*((B16+B17)  *Density%VecF(:,:,:,:,it)             &
         &          - B16                *Density%VecF(:,:,:,:,at))
       endif
       if(B18.ne.0.0_dp.or.B19.ne.0.0_dp.or.B20.ne.0.0_dp .or. B21.ne.0.0_dp)then
         SPot(:,:,:,:,it) = SPot(:,:,:,:,it) &           
         &          + 2.0_dp*((B18+B19)  *Density%LapS(:,:,:,:,it)             &
         &          + B18                *Density%LapS(:,:,:,:,at))            &
         &          - 2.0_dp*((B20+B21)  *Density%GradDivS(:,:,:,:,it)         &
         &          - B20                *Density%GradDivS(:,:,:,:,at))
       endif
    enddo

    !Add the contribution from cranking
    SPot = SPot + CrankSPot()
    return
  end subroutine CalcSPot

  subroutine CalcDPot()
  !-----------------------------------------------------------------------------
  ! Subroutine that calculates the D potentials given by:
  !   D = - 2  b16  s_mu - 2  b17 s_qmu
  ! It also calculates the divergence of D.
  !-----------------------------------------------------------------------------
      integer      :: it
      
      DPot=0.0_dp
      DivDPot=0.0_dp
      if(B16.eq.0.0_dp .and. B17 .eq.0.0_dp) return
      
      !Calculating D itself
      do it=1,2
          DPot(:,:,:,:,it)=DPot(:,:,:,:,it)-2*B16*sum(Density%VecS(:,:,:,:,:),5)
      enddo

      do it=1,2
          DPot(:,:,:,:,it)=DPot(:,:,:,:,it)-2*B17*Density%VecS(:,:,:,:,it)
      enddo
      return
  end subroutine CalcDPot

  subroutine CalcCPot()
  !-----------------------------------------------------------------------------
  ! Subroutine that calculates the C potential, given by:
  !   C = - 2 B14 s_mu - 2 B15 s_qmu
  ! And also its derivatives.
  !-----------------------------------------------------------------------------
      integer :: it

      CPot=0.0_dp
      if(B14 .eq. 0.0_dp .and. B15 .eq. 0.0_dp) return
      do it=1,2
          CPot(:,:,:,:,it)=CPot(:,:,:,:,it)-2*B14*sum(Density%VecS(:,:,:,:,:),5)
      enddo
      do it=1,2
          CPot(:,:,:,:,it)=CPot(:,:,:,:,it)-2*B15*Density%VecS(:,:,:,:,it)
      enddo

      return
  end subroutine CalcCPot
  
  !=============================================================================
  !Subroutines for calculating the action of all potentials on the spwfs.
  !=============================================================================  
  function ActionOfU(Psi)
  !-----------------------------------------------------------------------------
  ! Subroutine that computes the action of the U potential on a wavefunction.
  !
  ! ActionOfU = U*Psi
  !-----------------------------------------------------------------------------

    type(Spwf), intent(in) :: Psi
    type(Spinor)           :: ActionOfU, Value
    
    integer :: it
    
    it = (Psi%GetIsoSpin() + 3)/2
    ActionOfU = Psi%GetValue()
    ActionOfU = UPot(:,:,:, it) * ActionOfU
          
    return                                       
  end function ActionOfU
  
  function ActionOfBNew(Psi,NoKinetic) result(ActionOfB)
  !-----------------------------------------------------------------------------
  ! Experimental way of calculating the action of B on a wavefunction.
  !-----------------------------------------------------------------------------
  ! Don't use this, as the ActionofBOld is more consistent for the spwf energies
  !-----------------------------------------------------------------------------

    type(Spwf), intent(in) :: Psi
    type(Spinor)           :: ActionOfB, Lap, Der(3), SumDer
    integer                :: it, m
    logical, intent(in), optional :: NoKinetic

    if(present(NoKinetic)) then
      call stp('NoKinetic not implemented for this subroutine!')
    endif

    it = (Psi%GetIsospin()+3)/2                        
    Lap       = Psi%GetLap()
    ActionOfB =   Lap
    ActionOfB = - BPot(:,:,:,it) * ActionOfB  
    do m=1,3
      Der(m) = Psi%GetDer(m)
      Der(m) = NablaBPot(:,:,:,m,it) * Der(m)
    enddo
    SumDer = Der(1) + Der(2) + Der(3)
    ActionOfB = ActionOfB - SumDer
    return
  end function ActionOfBNew

  function ActionOfBOld(Psi, NoKinetic) result(ActionOfB)
  !-----------------------------------------------------------------------------
  ! Subroutine that computes the action of the B potential on a wavefunction.
  !       
  !       ActionOfB = - \nabla B \nabla \Psi
  !-----------------------------------------------------------------------------
  ! This is currently computed in exactly the same way as in CR8. However, this 
  ! implies some stacked derivatives, and I'm not completely sure if they are 
  ! reliable.
  !-----------------------------------------------------------------------------
    use Derivatives
    
    type(Spwf), intent(in) :: Psi
    type(Spinor)           :: ActionOfB, Lap, Der(3), DerSec(3)

    integer :: it, m, at, Parity, Signature, TimeSimplex, P, S, TS
    logical, intent(in), optional :: NoKinetic
    real(KIND=dp) :: Reducedmass(2)

    Reducedmass = (1.0_dp-nucleonmass/(neutrons*nucleonmass(1)+protons*nucleonmass(2)))
    ! Getting the symmetry quantum Numbers
    P = Psi%GetParity()
    S = Psi%GetSignature()
    TS = Psi%GetTimeSimplex()
    it = (Psi%GetIsospin()+3)/2                        
    at = 3-it

    ActionOfB=NewSpinor()
    Lap      =NewSpinor()
    do m=1,3
      Der(m) = NewSpinor(); DerSec(m) = NewSpinor()
    enddo

    if(.not. present(NoKinetic)) then
      Lap       = Psi%GetLap()                                        
      ActionOfB = (- hbm(it)/2.0_dp ) * Lap
      if(COM1Body.ne.0) then
        ActionOfB =   Reducedmass(it)*ActionOfB
      endif
    endif                       
    !X component
    Der(1) = Psi%GetDer(1)
    Der(1) = (B3*Density%Rho(:,:,:,at) + (B3+B4)*Density%Rho(:,:,:,it)) * Der(1)

    !X component
    Parity = -P; Signature = -S; TimeSimplex =  TS

    do m=1,4
        DerSec(1)%Grid(:,:,:,m,1) = DeriveX(Der(1)%Grid(:,:,:,m,1), Parity,    &
        &                                   Signature,TimeSimplex,m)
    enddo
    
    ActionOfB = ActionOfB - DerSec(1)

    !Y Component

    Der(2) = Psi%GetDer(2)
    Der(2) = (B3*Density%Rho(:,:,:,at) + (B3+B4)*Density%Rho(:,:,:,it)) * Der(2)

    Parity = -P; Signature = -S; TimeSimplex = -TS

    do m=1,4
        DerSec(2)%Grid(:,:,:,m,1) = DeriveY(Der(2)%Grid(:,:,:,m,1), Parity,    &
        &                                   Signature,TimeSimplex,m)
    enddo
    
    ActionOfB = ActionOfB - DerSec(2)

    !Z Component

    Der(3) = Psi%GetDer(3)
    Der(3) = (B3*Density%Rho(:,:,:,at) + (B3+B4)*Density%Rho(:,:,:,it)) * Der(3)

    Parity = -P; Signature =  S; TimeSimplex =  TS

    do m=1,4
        DerSec(3)%Grid(:,:,:,m,1) = DeriveZ(Der(3)%Grid(:,:,:,m,1), Parity,    &
        &                                   Signature,TimeSimplex,m)
    enddo
    
    ActionOfB = ActionOfB - DerSec(3)

    return
  end function ActionOfBOld
                
  function ActionOfW(Psi)
    !---------------------------------------------------------------------------
    ! Function that computes the action of the W Potential on a wavefunction Psi
    !     ActionOfW = - i W_{mu,nu} \nabla_{mu} \sigma_{nu} (implied summation)
    !---------------------------------------------------------------------------
    !
    ! Note that the terms related to \nabla W are taken to be zero. This is 
    ! clearly justified when no tensor terms are present (i.e. b14,15,16,17 = 0)
    ! but unclear when those terms are present.
    !
    !---------------------------------------------------------------------------
    type(Spwf), intent(in) :: Psi
    type(Spinor)           :: ActionOfW, Der(3), P
    
    integer :: it, mu, nu
    
    ActionOfW=NewSpinor()
    it = (Psi%GetIsospin() + 3)/2
    do mu = 1,3
            Der(mu) = Psi%GetDer(mu)
            do nu=1,3
                    P = WPot(:,:,:,mu,nu,it)*Pauli(Der(mu),nu)
                    ActionOfW = ActionOfW + P
            enddo
    enddo
    ActionOfW = MultiplyI(ActionOfW)
    ActionOfW = (-1.0_dp)*ActionOfW
    return
  end function ActionOfW
                
  function ActionOfA(Psi)
    !---------------------------------------------------------------------------
    ! Function that computes the action of the A Potential
    !       - i \vec{A}_{q}(\vec{r}) \cdot \nabla \Psi
    !
    !---------------------------------------------------------------------------
    ! Note that 
    ! 1) A corresponds to -C from the Mg24 paper, meaning to the pmx/y/z 
    !    potential from cr8, without extra sign!
    ! 2) the terms involving the divergence of A are taken to be 0
    !---------------------------------------------------------------------------

    type(Spwf), intent(in) :: Psi
    integer                :: it, m
    type(Spinor)           :: ActionOfA, Der
    
    it = (Psi%getIsoSpin() + 3)/2
		ActionOfA=NewSpinor()
    do m=1,3
      Der = Psi%GetDer(m)
      ActionOfA = ActionOfA + APot(:,:,:,m,it)*Der
    enddo
    ActionOfA = MultiplyI(ActionOfA)
    ActionOfA = -ActionOfA
    return
  end function ActionOfA
                
  function ActionOfS(Psi)
    !---------------------------------------------------------------------------
    ! Function that computes the action of the S Potential.
    !       ActionOfS = \vec{S} \cdot \vec{\sigma} \Psi
    !---------------------------------------------------------------------------  
    integer                :: m, it
    type(Spwf), intent(in) :: Psi
    type(Spinor)           :: ActionOfS, Value
    
    it = (Psi%GetIsoSpin() + 3)/2
		ActionOfS = NewSpinor()
    do m=1,3
      Value = Psi%GetValue()
      Value = Pauli(Value, m)
      Value = SPot(:,:,:,m,it)*Value
      ActionOfS = ActionOfS + Value
    enddo
    return
  end function ActionOfS

  function ActionOfD(Psi)
    !---------------------------------------------------------------------------
    ! Function that computes the action of the D potential
    !   ActionOfD = - (\nabla \cdot D_q) * (\sigma \cdot \nabla) \Psi
    !---------------------------------------------------------------------------
    integer                 :: m, it
    type(Spwf) , intent(in) ::Psi
    type(Spinor)            ::ActionOfD, Der(3)

    it = (Psi%GetIsospin() + 3)/2

    do m=1,3
        Der(m) = Psi%GetDer(m)
        Der(m) = Pauli(Der(m),m)
    enddo
    ActionOfD = Der(1)
    ActionOfD = ActionOfD + Der(2)
    ActionOfD = ActionOfD + Der(3)

    ActionOfD = - DivDPot(:,:,:,it) * ActionOfD
  end function ActionOfD

  function ActionOfC(Psi)
    !---------------------------------------------------------------------------
    ! Function that computes the action of the C potential
    !  ActionOfC = - \nabla \cdot [\sigma \cdot C_q ] \nabla
    !---------------------------------------------------------------------------
    integer                 :: m, it, n
    type(Spwf) , intent(in) :: Psi
    type(Spinor)            :: ActionOfC, Value, Der(3), Lap, Temp

    Value = Psi%GetValue()
    it    = (Psi%GetIsospin() + 3)/2

    !Terms involving the laplacian
    Lap = Psi%GetLap()

    ActionOfC = NewSpinor()
    do m=1,3
        Temp = Pauli(Lap,m)
        ActionOfC = ActionOfC  - CPot(:,:,:,m,it) * Temp
    enddo

    !Terms involving the derivatives
    do m=1,3
        Der(m) = Psi%GetDer(m)
    enddo
    do m=1,3
        do n=1,3
            Temp = Pauli(Der(m), n)
            ActionOfC = ActionOfC - DerCPot(:,:,:,n,m,it)*Temp
        enddo
    enddo
  end function ActionOfC 
  
end module MeanFields
