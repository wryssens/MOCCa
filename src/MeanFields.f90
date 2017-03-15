module MeanFields
  !-----------------------------------------------------------------------------
  ! 
  ! Every mean-field density has a corresponding field and action that is 
  ! calculated in this module. Included at the moment are
  !
  !     Field    Density   Remarks                          T
  !     U        rho                                        even
  !     B_mn     tau_mn    Diagonal if functional is NLO    even (+odd for N2LO)
  !     W_mn     J_mn                                       even
  !     S        s                                          odd
  !     A        j                                          odd
  !     DN2LO    Q         N2LO only                        even
  !     Fs       S         N2LO only                        odd
  !     FT       Tmnk      N2LO only                        even + odd
  !     Fpi      Pi        N2LO only                        odd
  !     Fv       V_mn      N2LO only                        even
  !     C        Tk        Tensor only                      odd
  !     Dt       Fk        Tensor only                      odd
  !
  ! Where we use Dt to mean the tensor field D, since we are using D_m as 
  ! derivatives in this documentation.
  !
  ! The terms in the single-particle hamiltonian are
  !  
  !  h_q  = hb2m * D            (kinetic energy)
  !         +   U 
  !         -   Dm  B_mn Dn
  !         - i W_mn Dm sigma_n
  !         +   S_m sigma_m  
  !         -   A_m nabla_m
  !         +   Delta DN2LO Delta                          
  !         +   Delta Fs_m  Delta sigma_m 
  !         -   D_m Ft_mn D_n sigma_k
  !         + i D_n Fpi_m D_n D_m
  !         +   D_k F^V_mn D_k D_m sigma_n
  !         -   D_m sigma_n C_n D_m
  !         -   D_m Dt_m sigma_n D_n   
  !
  ! Note that this equation is only valid for local, i.e. gauge-invariant 
  ! functionals.
  !-----------------------------------------------------------------------------
  !
  ! Relevant papers: 
  ! P. Bonche et al., Nucl. Phys. A467 (1987), 115-135
  ! V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
  ! P. Becker et al., J. Phys. G. 40 (2015) 
  !
  ! Take care for various errors and different conventions in these papers.
  ! 
  ! Conventions in the Mg24 paper
  !  * MOCCa's A potential corresponds to - C, which is the pmx/y/z in CR8
  !  * MOCCa's S potential is the V potential 
  !
  ! Tensor paper
  !  * B9 contribution to S potential has a relative + sign, not a minus as in
  !    paper
  !  * The action of the D potential is not hermitian. See the routine ActionofD
  !
  ! Tools paper
  !  * The potential of the Vmn density is completely absent, while it should be
  !    included in the single-particle hamiltonian.
  !  * The imaginary part of the Tmnk density is time-reversal invariant, and 
  !    should thus be present in both the functional and the single-particle 
  !    hamiltonian.
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
  real(KIND=dp),allocatable :: BPot(:,:,:,:)          , NablaBPot(:,:,:,:,:)
  real(KIND=dp),allocatable :: UPot(:,:,:,:)          , APot(:,:,:,:,:)
  real(KIND=dp),allocatable :: SPot(:,:,:,:,:)
  real(KIND=dp),allocatable :: Cpot(:,:,:,:,:)        , WPot(:,:,:,:,:,:)
  real(KIND=dp),allocatable :: DPot(:,:,:,:,:)
  real(KIND=dp),allocatable :: DerCPot(:,:,:,:,:,:)   , DerDPot(:,:,:,:,:,:)
  real(KIND=dp),allocatable :: DivDpot(:,:,:,:)
  !------------------------------------------------------------------------------
  ! N2LO potentials
  real(KIND=dp),allocatable :: Bmunu(:,:,:,:,:,:)     , DN2LO(:,:,:,:)
  real(KIND=dp),allocatable :: Xpot(:,:,:,:,:,:)      , DXpot(:,:,:,:,:,:,:)
  real(KIND=dp),allocatable :: ReTfield(:,:,:,:,:,:,:), ImTfield(:,:,:,:,:,:)
  real(KIND=dp),allocatable :: ReDTfield(:,:,:,:,:,:) , PiField(:,:,:,:,:)
  real(KIND=dp),allocatable :: DmuBmunu(:,:,:,:,:,:)  , SN2LOField(:,:,:,:,:)
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
    use Derivatives
    !Allocate the potentials if this has not already happened
    if(.not.allocated(UPot)) then
      !-------------------------------------------------------------------------
      ! Time-even potentials
      allocate(BPot(nx,ny,nz,2), NablaBPot(nx,ny,nz,3,2), UPot(nx,ny,nz,2))
      allocate(Wpot(nx,ny,nz,3,3,2))
      BPot     =0.0_dp ; NablaBPot=0.0_dp ; Upot = 0.0_dp
      Wpot     =0.0_dp 
      if(.not. TRC) then
        !-----------------------------------------------------------------------
        ! Time-odd potentials
        allocate(SPot(nx,ny,nz,3,2),APot(nx,ny,nz,3,2))
        SPot     =0.0_dp ; Apot = 0.0_dp
      endif
     
      !-------------------------------------------------------------------------
      ! N2LO potentials, time-even 
      if(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp) then
        allocate(Bmunu(nx,ny,nz,3,3,2))     ; allocate(DmuBmunu(nx,ny,nz,3,3,2))
        allocate(DN2LO(nx,ny,nz,2))         ; allocate(Xpot(nx,ny,nz,3,3,2))
        allocate(DXpot(nx,ny,nz,3,3,3,2))   ; allocate(ImTfield(nx,ny,nz,3,3,2))
        
        Bmunu = 0.0_dp    ; DmuBmunu  = 0.0_dp
        DN2LO = 0.0_dp    ; Xpot    = 0.0_dp
        DXpot = 0.0_dp    ; ImTField= 0.0_dp
        
        if(.not. TRC) then
            !-------------------------------------------------------------------
            ! Time-odd N2LO potentials
            allocate(SN2LOField(nx,ny,nz,3,2))   ; SN2LOField = 0.0_dp
            allocate(ReTfield(nx,ny,nz,3,3,3,2)) ; ReTfield   = 0.0_dp 
            allocate(ReDTfield(nx,ny,nz,3,3,2))  ; ReDTfield  = 0.0_dp
            allocate(PiField(nx,ny,nz,3,2))      ; PiField    = 0.0_dp
        endif
      endif
      
      !-------------------------------------------------------------------------
      ! Tensor potentials, all time-odd
      if(B15.ne.0.0_dp .or. B16.ne.0.0_dp .or. &
      &  B17.ne.0.0_dp .or. B18.ne.0.0_dp) then
        allocate(Cpot(nx,ny,nz,3,2),DPot(nx,ny,nz,3,2))
        allocate(DerCPot(nx,ny,nz,3,3,2),DerDPot(nx,ny,nz,3,3,2))
        allocate(DivDpot(nx,ny,nz,2))
        
        CPot     =0.0_dp ; DPot     =0.0_dp
        DerCPot  =0.0_dp ; DerDPot  =0.0_dp
        DivDpot  =0.0_dp
      endif    
    endif

    !---------------------------------------------------------------------------
    ! Decide on the way to calculate the action of the kinetic field B_mn.
    ! Done here for the ifort compiler
    if (BStack) then
      if(t1n2 .ne. 0.0_dp .or. t2n2 .ne. 0.0_dp) then
        ! Derivatives are stacked multiple times for the N2LO functional, 
        ! thus it is only safe to calculate it with Lagrange derivatives.
        call stp('Cannot stack derivatives for the B-field with N2LO.')
      else
        ActionOfB=> ActionOfBOld
      endif
    else
      if(t1n2 .ne. 0.0_dp .or. t2n2 .ne. 0.0_dp) then
        ActionOfB=> ActionOfBN2LO
      else
        ActionOfB=> ActionOfBNew
      endif
    endif
    !---------------------------------------------------------------------------    
    !Calculate the Coulomb Potential
    call SolveCoulomb()
    !Construct the mean-field potentials
    call CalcUPot() !=> This includes the contribution of constraints & Coulomb
    call CalcWPot()
    call CalcBPot()  
        
    if(.not.TRC) then
        call CalcSPot()
        call CalcAPot()
        call CalcCPot()
        call CalcDPot()
    endif
    !---------------------------------------------------------------------------
    if(t1n2 .ne. 0.0_dp .or. t2n2 .ne. 0.0_dp) then
        ! All of the N2LO potentials
        call CalcBMunu()
        call CalcXpot()
        call calcDN2LO()
        call calcTfield()
        if(.not. TRC) then
            call calcPifield()
            call calcSN2LOField()
        endif
    endif
  end subroutine ConstructPotentials

  !=============================================================================
  ! Subroutines for calculating all the different potentials.
  !=============================================================================
  subroutine CalcBPot()
  !-----------------------------------------------------------------------------
  ! This subroutine computes the modified mass B_q(r).
  !       B  =  B_3 \rho + B_4 \rho_q
  !
  ! Note that this is only the diagonal component of the B_munu effective mass
  ! when N2LO terms are present.
  !-----------------------------------------------------------------------------
    use Derivatives
    
    integer :: it, at, i
    real(KIND=dp) :: Reducedmass(2)

    Reducedmass=(1.0_dp-nucleonmass/(neutrons*nucleonmass(1)+protons*nucleonmass(2)))

    do it=1,2
        at = 3 - it
        BPot(:,:,:,it) =B3*Density%Rho(:,:,:,at) + (B3+B4)*Density%Rho(:,:,:,it)
    enddo
    
    do it=1,2
      !Gradient of B
      NablaBPot(:,:,:,1,it) = &
      & DeriveX(BPot(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt, 1)
      NablaBPot(:,:,:,2,it) = &
      & DeriveY(BPot(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt, 1)
      NablaBPot(:,:,:,3,it) = &
      & DeriveZ(BPot(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt, 1)
    enddo
    return
  end subroutine CalcBPot
  
  subroutine calcBmunu()
    !---------------------------------------------------------------------------
    ! Calculate the Bmn tensor for the N2LO functional. Note that it collects 
    ! the NLO part from the calcBpot routine.
    !
    ! B_mn = B(r) delta_mn + 2 C^(4 Drho) tau delta_mn
    !                      + 4 C^(4 Mrho) tau_mn
    !                      - 2 C^(4 Mrho) D_m D_n rho 
    ! 
    ! This routine also calculates the derivatives of B_mn
    !
    ! DmuBmunu = sum_m D_m B_mn                     
    !---------------------------------------------------------------------------

    use Derivatives
    
    integer :: mu,nu, it, at

    Bmunu = 0.0_dp
    !---------------------------------------------------------------------------
    ! N1LO diagonal elements 
    do mu=1,3
        Bmunu(:,:,:,mu,mu,:) = Bpot
    enddo
    
    !---------------------------------------------------------------------------
    ! tau delta_munu terms
    do it=1,2
        do mu=1,3
            Bmunu(:,:,:,mu,mu,it) = Bmunu(:,:,:,mu,mu,it)                      &
            &                      +2*N2tau(1)*sum(Density%tau,4)              &
            &                      +2*N2tau(2)*    Density%tau(:,:,:,it) 
        enddo
    enddo   
        
    !---------------------------------------------------------------------------
    ! Off-diagonal terms
    do it=1,2
        at = 3 - it
        do mu=1,3
            do nu=1,3
                Bmunu(:,:,:,mu,nu,it) = Bmunu(:,:,:,mu,nu,it)                     &
                &    +4*(N2rtaumn(1) + N2rtaumn(2))* Density%RTauN2LO(:,:,:,mu,nu,it) &
                &    -2*(N2tddr(1) + N2tddr(2))    * Density%D2Rho   (:,:,:,mu,nu,it) & 
                &    +4*(N2rtaumn(1))              * Density%RTauN2LO(:,:,:,mu,nu,at) &
                &    -2*(N2tddr(1))                * Density%D2Rho   (:,:,:,mu,nu,at)       
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    ! Derivatives of B_munu
    ! Note that B_munu shares all of the symmetries of t_munu
    ! (at least its real part)
    do it=1,2
        DmuBmunu(:,:,:,1,1,it) =                                          &
        &  DeriveX(Bmunu(:,:,:,1,1,it), ParityInt, SignatureInt, TimeSimplexInt,1)
        DmuBmunu(:,:,:,2,2,it) =                                          &
        &  DeriveY(Bmunu(:,:,:,2,2,it), ParityInt, SignatureInt, TimeSimplexInt,1)
        DmuBmunu(:,:,:,3,3,it) =                                          &
        &  DeriveZ(Bmunu(:,:,:,3,3,it), ParityInt, SignatureInt, TimeSimplexInt,1)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        DmuBmunu(:,:,:,1,2,it) =                                          &
        &  DeriveX(Bmunu(:,:,:,1,2,it), ParityInt, SignatureInt, TimeSimplexInt,2)
        DmuBmunu(:,:,:,2,1,it) =                                          &
        &  DeriveY(Bmunu(:,:,:,2,1,it), ParityInt, SignatureInt, TimeSimplexInt,2)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        DmuBmunu(:,:,:,3,1,it) =                                          &
        &  DeriveZ(Bmunu(:,:,:,3,1,it), ParityInt,-SignatureInt, TimeSimplexInt,1)
        DmuBmunu(:,:,:,1,3,it) =                                          &
        &  DeriveX(Bmunu(:,:,:,1,3,it), ParityInt,-SignatureInt, TimeSimplexInt,1)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        DmuBmunu(:,:,:,2,3,it) =                                          &
        &  DeriveY(Bmunu(:,:,:,2,3,it), ParityInt,-SignatureInt, TimeSimplexInt,2)
        DmuBmunu(:,:,:,3,2,it) =                                          &
        &  DeriveZ(Bmunu(:,:,:,3,2,it), ParityInt,-SignatureInt, TimeSimplexInt,2)
    enddo
  end subroutine calcBmunu

  subroutine calcDN2LO()
    !---------------------------------------------------------------------------
    ! Calculates the N2LO potential D, associated with the N2LO density Q
    !
    ! DN2LO = C^(4 Mrho) rho
    !---------------------------------------------------------------------------
    
    use Derivatives
    use Force
    
    integer :: it
    DN2LO = 0.0_dp
    
    do it=1,2
        DN2LO(:,:,:,it)    = N2rhoQ(1)*sum(Density%Rho,4)                      &
        &                  + N2rhoQ(2)*    Density%Rho(:,:,:,it) 
    enddo
    
  end subroutine calcDN2LO
  
  subroutine calcSN2LOField()
    !---------------------------------------------------------------------------
    ! Calculates the N2LO potential Fs, associated with the N2LO density (big) S
    !
    ! Fs_m = C^(4, Ms) s_m
    !
    !---------------------------------------------------------------------------
    use Derivatives
    use Force
    
    integer :: it
    SN2LOField = 0.0_dp
    
    do it=1,2
        SN2LOField(:,:,:,:,it) = N2sS(1)*sum(Density%vecs,5)                   &
        &                      + N2sS(2)*    Density%vecs(:,:,:,:,it) 
    enddo
  end subroutine calcSN2LOField

  subroutine CalcUPot()
      !-------------------------------------------------------------------------
      ! This subroutine calculates the single particle potential U.
      !
      !        Time-even part
      !        --------------
      ! U   =   2 B1  rho      + 2 B2 rho_q
      !       +   B3  tau      +   B4 tau_q
      !       + 2 B5 Delta rho + 2 B6 Delta rho_q
      !       + (2+\byt3a) B7 rho^(alpha+1) 
      !       + (2+\byt3a) B8 rho^(alpha-1) [rho_n^2 + rho_p^2]  
      !       +   B9  D_k Jmn  + B9 D_k Jmn_q
      ! 
      !        Time-odd part
      !        -------------
      ! U   =    B12 alpha*rho^alpha s^2 + B13 alpha*rho^alpha [s_n^2 + s_p^2]
      !
      !
      !        N2LO part
      !        -------------
      ! U   =  2 C^(4Dr) Delta Delta rho + C^(4 Mrho) Q 
      !      - 2 C^(4 Mrho) D_m D_n tau_mn
      !-------------------------------------------------------------------------
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
            UPot(:,:,:,it) = UPot(:,:,:,it)                                    &
            & + (B3+B4)*Density%Tau(:,:,:,it) + B3*Density%Tau(:,:,:,at)       &
            & + 2.0_dp*((B5 + B6)*Density%LapRho(:,:,:,it)                     &
            &          + B5      *Density%LapRho(:,:,:,at))                    &
            & + (2+byt3a)*B7a*(RhoTot+eps)**(1+byt3a)                          &
            & + B8a*(                                                          &
            & byt3a*(RhoTot**(byt3a))/(RhoTot+eps)*                            &
            & (Density%Rho(:,:,:,it)**2 + Density%Rho(:,:,:,at)**2)            &
            & + 2.0_dp*(RhoTot)**(byt3a)*Density%Rho(:,:,:,it)                 &
            & )                                                                &
            & + B9*NablaJTot + B9q*Density%NablaJ(:,:,:,it)

            if(t1n2 .ne. 0.0_dp .or. t2n2 .ne. 0.0_dp) then
                !---------------------------------------------------------------
                ! n2lo contribution
                !---------------------------------------------------------------
                upot(:,:,:,it) = upot(:,:,:,it)                                &
                &              + 2*N2D2rho(1)*sum(density%laplaprho,4)         &
                &              + 2*N2D2rho(2)*    density%laplaprho(:,:,:,it)  &
                &              +   N2rhoQ(1) *sum(density%qn2lo,4)             &
                &              - 2*N2tddr(1) *sum(density%d2rtau,4)            &
                &              +   N2rhoQ(2)*    density%qn2lo(:,:,:,it)       &
                &              - 2*N2tddr(2)*    density%d2rtau(:,:,:,it)      
            endif
    enddo

    if(.not. TRC) then
       do it=1,2
          at = 3- it
          Upot(:,:,:,it) = Upot(:,:,:,it) + byt3a*RhoTot**(byt3a)/(RhoTot + eps) &
          & * ( B12a*sum(VecSTot**2,4) + B13a*sum(Density%VecS(:,:,:,:,it)**2,4) &
          &                            + B13a*sum(Density%VecS(:,:,:,:,at)**2,4))
       enddo
    endif

    ! Note that CoulExchange from the Coulomb module already carries a minus.
    UPot(:,:,:,2)= UPot(:,:,:,2) + CoulombPotential(:,:,:) + CoulExchange(:,:,:)

    !Add the contribution from Constraints on the Multipole Moment
    UPot = UPot + ConstraintEnergy
    return
  end subroutine CalcUPot

  subroutine CalcAPot()
  !-----------------------------------------------------------------------------
  ! This subroutine calculates the A potential.
  !
  !               Time-odd
  !               ---------
  ! A_q = - 2 B3 j     + 2 B4  j_q
  !       +   B9 D x s +   B9q D x s_q)
  !
  !               N2LO
  !               ---------
  ! A_q =   4 C^(4 Mrho)_q Pi
  !
  ! Also added is the contribution from the cranking constraints.
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
    
    if(t1n2 .ne. 0.0_dp .or. t2n2.ne.0.0_dp) then
        !N2LO contribution
        do it=1,2
            at = 3 - it
            Apot(:,:,:,:,it) = Apot(:,:,:,:,it) &
            &                - (N2jpi(1) + N2jpi(2))*Density%PiN2LO(:,:,:,:,it) &
            &                -  N2jpi(1)            *Density%PiN2LO(:,:,:,:,at) 
        enddo
    endif
    
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
             &        + B9       *Density%DerRho(:,:,:,k,3-it))
           enddo
         enddo
      enddo
    enddo

    if(B14.ne.0.0_dp .or. B15.ne.0.0_dp) then
        do it=1,2
           !B14&15 Contribution
           WPot(:,:,:,:,:,it) =   WPot(:,:,:,:,:,it)                   &
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
             & + 2.0_dp * (B16+B17)*Density%Jmunu(:,:,:,n,m,it)        &
             & + 2.0_dp *  B16     *Density%Jmunu(:,:,:,n,m,3-it)

           enddo
           !b) \delta_{mu,nu}
           do k=1,3
                WPot(:,:,:,m,m,it) = WPot(:,:,:,m,m,it)                &
                & + 2 * (B16+B17) * Density%Jmunu(:,:,:,k,k,it)         &
                & + 2 *  B16      * Density%Jmunu(:,:,:,k,k,3-it)
           enddo
        enddo
      enddo
    endif
    
    if(BN2LO(7).ne.0.0_dp .or. BN2LO(8).ne.0.0_dp) then
        !-----------------------------------------------------------------------
        ! N2LO contribution
        do it=1,2
            Wpot(:,:,:,:,:,it) = Wpot(:,:,:,:,:,it)                    &
            &                  - 4 * N2JV(1) * sum(Density%VN2LO,6)    & 
            &                  - 4 * N2JV(2) *     Density%VN2LO(:,:,:,:,:,it) 
        enddo
        !-----------------------------------------------------------------------
        ! Note that the N2LO term proportional to B^(4)_14 and 15 is not 
        ! implemented, since it should be zero for local interactions anyway.
        !-----------------------------------------------------------------------
    endif

    return
  end subroutine CalcWPot
  
  subroutine CalcXPot()
    !-------------------------------------------------------------------------
    ! Calculates the X_munu N2LO potential. 
    ! X_munu \sim J_munu
    !-------------------------------------------------------------------------
    integer :: it, mu,nu 

    Xpot= 0.0_dp

    do it=1,2
        ! Note that the final derivative is first in the array
        Xpot (:,:,:,:,:,it)  = - 4 * N2JV(1) * sum(Density%Jmunu,6)          &
        &                      - 4 * N2JV(2) * Density%Jmunu(:,:,:,:,:,it)
        DXpot(:,:,:,:,:,:,it)= - 4 * N2JV(1) * sum(Density%DJmunu,7)         &
        &                      - 4 * N2JV(2) * Density%DJmunu(:,:,:,:,:,:,it)
    enddo
      
  end subroutine CalcXpot

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
    !---------------------------------------------------------------------------
    use Cranking, only : CrankSPot

    integer        :: it,at,l,i,j,k
    real(KIND=dp)  :: RhoT(nx,ny,nz), RhoVecS(nx,ny,nz,3,2)

    Spot = 0.0_dp
    RhoT = Density%Rho(:,:,:,1) + Density%Rho(:,:,:,2)

    do it=1,2
      do l=1,3
        RhoVecS(:,:,:,l,it) = (RhoT**(byt3a))*Density%VecS(:,:,:,l,it)
      enddo
    enddo
    do it=1,2
      at = 3 - it
      Spot(:,:,:,:,it) = &
       &          + (B9 + B9q)         *Density%RotVecJ(:,:,:,:,it)   &
       &          +  B9                *Density%RotVecJ(:,:,:,:,at)   &
       &          + 2.0_dp*((B10+B11)  *Density%Vecs(:,:,:,:,it)      &
       &          +          B10       *Density%VecS(:,:,:,:,at))     &
       &          + 2.0_dp*((B12a+B13a)*        RhoVecS(:,:,:,:,it)   &
       &          +          B12a      *        RhoVecS(:,:,:,:,at))

       if(B14.ne.0.0_dp .or. B15.ne.0.0_dp) then
        SPot(:,:,:,:,it) = SPot(:,:,:,:,it) &
        &          -         (B14+B15)  *Density%VecT(:,:,:,:,it)     &
        &          -          B14       *Density%VecT(:,:,:,:,at)
       endif
       if(B16.ne.0.0_dp .or. B17.ne.0.0_dp) then
         SPot(:,:,:,:,it) = SPot(:,:,:,:,it) &
         &          - 2.0_dp* (B16+B17)  *Density%VecF(:,:,:,:,it)             &
         &          - 2.0_dp*  B16       *Density%VecF(:,:,:,:,at)
       endif
       if(B18.ne.0.0_dp.or.B19.ne.0.0_dp.or.B20.ne.0.0_dp .or. B21.ne.0.0_dp)then
         SPot(:,:,:,:,it) = SPot(:,:,:,:,it) &
         &          + 2.0_dp*((B18+B19)  *Density%LapS(:,:,:,:,it)             &
         &          + B18                *Density%LapS(:,:,:,:,at))            &
         &          - 2.0_dp*((B20+B21)  *Density%GradDivS(:,:,:,:,it)         &
         &          - B20                *Density%GradDivS(:,:,:,:,at))
       endif
    enddo

    if( BN2LO(5).ne.0.0_dp .or. BN2LO(6) .ne. 0.0_dp) then
        !-----------------------------------------------------------------------
        ! N2LO contribution
        !-----------------------------------------------------------------------
        do it=1,2
            at = 3 - it
            Spot(:,:,:,:,it)= Spot(:,:,:,:,it)   &
            &             +  2*(N2D2S(1)+N2D2S(2))*Density%LapLapS(:,:,:,:,it)&
            &             +  2*N2D2S(1)          *Density%LapLapS(:,:,:,:,at)&
            &             +   (N2sS(1)+N2sS(2))*Density%SN2LO(:,:,:,:,it)    &
            &             +    N2sS(1)         *Density%SN2LO(:,:,:,:,at)    &
            &   + 2*(N2TmnD2S(1) + N2TmnD2S(2))*Density%ReD2TN2LO(:,:,:,:,it)&
            &   + 2* N2TmnD2S(1)               *Density%ReD2TN2LO(:,:,:,:,at)
        enddo
    endif

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
      use Derivatives
      
      integer      :: it

      if(B16.eq.0.0_dp .and. B17 .eq.0.0_dp) return

      DPot   =0.0_dp
      DerDpot=0.0_dp



      !Calculating D itself
      do it=1,2
          DPot(:,:,:,:,it)= -2 * B16 * sum(Density%VecS(:,:,:,:,:),5) &
          &                 -2 * B17 *     Density%VecS(:,:,:,:,it)
      enddo

      !Calculating the derivatives of the Dpotential
        do it=1,2
            DerDPot(:,:,:,1,1,it) = &
            & DeriveX(DPot(:,:,:,1,it), ParityInt,-SignatureInt,TimeSimplexInt, 1)
            DerDPot(:,:,:,2,1,it) = &
            & DeriveY(DPot(:,:,:,1,it), ParityInt,-SignatureInt,TimeSimplexInt, 1)
            DerDPot(:,:,:,3,1,it) = &
            & DeriveZ(DPot(:,:,:,1,it), ParityInt,-SignatureInt,TimeSimplexInt, 1)

            DerDPot(:,:,:,1,2,it) = &
            & DeriveX(DPot(:,:,:,2,it), ParityInt,-SignatureInt,TimeSimplexInt, 2)
            DerDPot(:,:,:,2,2,it) = &
            & DeriveY(DPot(:,:,:,2,it), ParityInt,-SignatureInt,TimeSimplexInt, 2)
            DerDPot(:,:,:,3,2,it) = &
            & DeriveZ(DPot(:,:,:,2,it), ParityInt,-SignatureInt,TimeSimplexInt, 2)

            DerDPot(:,:,:,1,3,it) = &
            & DeriveX(DPot(:,:,:,3,it), ParityInt, SignatureInt,TimeSimplexInt, 1)
            DerDPot(:,:,:,2,3,it) = &
            & DeriveY(DPot(:,:,:,3,it), ParityInt, SignatureInt,TimeSimplexInt, 1)
            DerDPot(:,:,:,3,3,it) = &
            & DeriveZ(DPot(:,:,:,3,it), ParityInt, SignatureInt,TimeSimplexInt, 1)
        enddo
        do it=1,2
            DivDPot(:,:,:,it) = DerDpot(:,:,:,1,1,it) + DerDpot(:,:,:,2,2,it)    &
            &                 + DerDpot(:,:,:,3,3,it)
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

      if(B14 .eq. 0.0_dp .and. B15 .eq. 0.0_dp) return
      CPot=0.0_dp
      do it=1,2
          CPot(:,:,:,:,it) = CPot(:,:,:,:,it)                                  &
          &                - B14*(Density%VecS(:,:,:,:,1)+Density%VecS(:,:,:,:,2))
      enddo
      do it=1,2
          CPot(:,:,:,:,it) = CPot(:,:,:,:,it)                                   &
          &                - B15*Density%VecS(:,:,:,:,it)
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

    type(Spwf), intent(in) :: Psi
    type(Spinor)           :: ActionOfB, Lap, Der(3), SumDer
    integer                :: it, m
    real(KIND=dp)          :: Reducedmass(2)
    logical, intent(in), optional :: NoKinetic

    if(present(NoKinetic)) then
      call stp('NoKinetic not implemented for this subroutine!')
    endif

    if(COM1body .eq. 2) then
        Reducedmass = (1.0_dp-nucleonmass/(neutrons*nucleonmass(1)+protons*nucleonmass(2)))
    else
        Reducedmass = 1.0_dp
    endif
    
    it        = (Psi%GetIsospin()+3)/2
    Lap       = Psi%GetLap()
    ActionOfB = Lap
    ActionOfB = - BPot(:,:,:,it) * ActionOfB
    
    !---------------------------------------------------------------------------
    ! Note that we explicitly include the action of the constant here, so that
    ! it does not end up in the derivative of B, nablapot. Lagrangian 
    ! derivatives cannot handle that, as the constant is not in the basis.
    ActionOfB = ActionOfB - hbm(it)/2.0_dp*Reducedmass(it) * Lap
    
    do m=1,3
        Der(m)  = Psi%GetDer(m)
        Der(m)  = NablaBPot(:,:,:,m,it) * Der(m)
    enddo
    SumDer    = Der(1) + Der(2) + Der(3)
    ActionOfB = ActionOfB - SumDer
    
    return
  end function ActionOfBNew
  
  function ActionOfBN2LO(Psi,NoKinetic) result(ActionOfB)
  !-----------------------------------------------------------------------------
  ! Calculate the action of the effective mass potential B_munu, which is 
  ! now a tensor in the case of N2LO terms being active.
  !-----------------------------------------------------------------------------
    type(Spwf), intent(in)        :: Psi
    type(Spinor)                  :: ActionOfB, temp, ImAction
    integer                       :: it, mu,nu, at
    real(KIND=dp)                 :: Reducedmass(2)
    logical, intent(in), optional :: NoKinetic

    if(present(NoKinetic)) then
      call stp('NoKinetic not implemented for this subroutine!')
    endif

    if(COM1body .eq. 2) then
        Reducedmass = (1.0_dp-nucleonmass/(neutrons*nucleonmass(1)+protons*nucleonmass(2)))
    else
        Reducedmass = 1.0_dp
    endif

    it = (Psi%GetIsospin()+3)/2
    ActionOfB = - (hbm(it)/2.0_dp*Reducedmass(it)) * Psi%Lap
    
    do mu=1,3
        do nu=1,3  
            ActionOfB = ActionOfB -    Bmunu(:,:,:,mu,nu,it)*Psi%SecondDer(mu,nu)
            ActionOfB = ActionOfB - DmuBmunu(:,:,:,mu,nu,it)*Psi%Der      (nu)
        enddo
    enddo
    
    if(TRC) return
    
    !---------------------------------------------------------------------------
    ! time-odd contribution to the action. Note that we have not calculated a 
    ! separate potential for this, as this field is essentially a scaled density
    !
    ! Note that Im Tau_munu is antisymmetric so,
    !
    ! nabla_mu tau_munu nabla_nu = [ nabla_mu tau_munu ] nabla_nu
    ImAction = NewSpinor()
    at = 3 - it
    do nu=1,3
        ImAction = ImAction + (4*(N2itaumn(1) + N2imtmn(2))*Density%DmuItau(:,:,:,nu,it)  &
        &                   +  4* N2itaumn(1) *             Density%DmuItau(:,:,:,nu,at)) &
        &                   *Psi%Der(nu)
    enddo
    ActionofB = ActionOfB - MultiplyI(ImAction)
    
  end function ActionOfBN2LO
  
  function ActionOfDN2LO(Psi) result(ActionOfD)
  !-----------------------------------------------------------------------------
  ! Calculate the action of the N2LO D potential 
  !-----------------------------------------------------------------------------
    type(Spwf), intent(in)        :: Psi
    type(Spinor)                  :: ActionOfD
    integer                       :: it, mu,nu

    it = (Psi%GetIsospin()+3)/2
    
    !---------------------------------------------------------------------------
    !Note that I have chosen to 'stack' derivatives here, since I've spent a 
    ! day on this, only to realize that 
    !  
    !  Lap ( D Lap ) != Lap (D) Lap + D Lap Lap
    !
    ! Only Lagrange derivatives are responsible anyway, and this is no problem
    ! for them. 
    !---------------------------------------------------------------------------
    ActionOfD =  LapSpinor(DN2LO(:,:,:,it)* Psi%Lap, Psi%parity, psi%signature, psi%timesimplex)
    
    return
  end function ActionOfDN2LO
  
  function ActionOfSN2LO(Psi) result(ActionOfS)
    type(Spwf), intent(in)        :: Psi
    type(Spinor)                  :: ActionOfS
    integer                       :: it, mu

    it = (Psi%GetIsospin()+3)/2
    ActionofS = NewSpinor()
    do mu=1,3
        ActionOfS = ActionOfS +                                                &
        & Pauli(LapSpinor(SN2LOField(:,:,:,mu,it)*Psi%Lap,                     &
        &       Psi%parity, psi%signature, psi%timesimplex), mu)    
    enddo
    
    return
  end function ActionOfSN2LO

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
  
  function ActionOfX(Psi) 
    !---------------------------------------------------------------------------
    ! Part of the single-particle hamiltonian, field X in the Meyer-Becker notes
    ! and Fv in my notation.
    !
    ! D_ka F^V_{mn} D_k D_m sigma_nu = [ D_ka F^V_{mn} ] D_ka D_m sigma_nu  + 
    !                                         F^V_{mn}   Lap  D_m sigma_nu
    !
    ! and in the end we add an extra factor i!
    !---------------------------------------------------------------------------
    type(Spwf), intent(in) :: Psi
    type(Spinor)           :: ActionOfX, temp(3)
  
    integer :: it, mu,nu, ka
  
    ActionOfX = NewSpinor()
    it = (Psi%GetIsospin() + 3)/2
    
    do mu=1,3
        Temp(mu) = newspinor()
    enddo
    ! D_m Lap Psi
    call DeriveSpinor(Psi%Lap, temp, Psi%Parity, Psi%Signature, Psi%TimeSimplex)
        
    
    do nu=1,3
        do mu=1,3
            do ka=1,3
                ActionOfX = ActionOfX + &
                &      DXpot(:,:,:,ka,mu,nu,it) * Pauli(Psi%SecondDer(ka,mu),nu)
            enddo
        enddo
        do mu=1,3
            ActionOfX = ActionOfX + &
            &          Xpot(:,:,:,mu,nu,it) * Pauli(temp(mu),nu)
        enddo
    enddo
    
    ActionOfX = MultiplyI(ActionofX)
  end function ActionOfX

  function ActionOfA(Psi)
    !---------------------------------------------------------------------------
    ! Function that computes the action of the A Potential
    !       - i \vec{A}_{q}(\vec{r}) \cdot \nabla \Psi
    !
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
    ActionOfA = -MultiplyI(ActionOfA)

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
      ActionOfS = ActionOfS + SPot(:,:,:,m,it)*Pauli(Value, m)
    enddo
    return
  end function ActionOfS

  function ActionOfD(Psi)
    !---------------------------------------------------------------------------
    ! Function that computes the action of the D potential
    ! Original from Veerles paper:
    !
    ! ActionOfD = - (\nabla \cdot D_q) * (\sigma \cdot \nabla) \Psi
    !
    ! But the corrected one is fully:
    !
    ! ActionOfD = -1/2 * [  Term_1 + Term_2 + 2*Term_3 ]
    !
    !
    ! Where
    !
    ! Term_1 = \nabla \cdot D \sigma \cdot \nabla \Psi
    ! Term_2 = \partial_\nu D_\mu \sigma_\nu \partial \mu \Psi
    !          (with sum on \mu and \nu)
    ! Term_3 = D_\mu \sigma_\nu \partial_\mu \partial_\nu \Psi
    !
    !---------------------------------------------------------------------------
    integer                 :: m, it,n,P,S,TS
    type(Spwf) , intent(in) :: Psi
    type(Spinor)            :: ActionOfD, Der(3), Term1, Term2,Term3
    type(Spinor)            :: DPsi(3,3)
    real(KIND=dp)           :: DivD(nx,ny,nz,2)

    it = (Psi%GetIsospin() + 3)/2
    do m=1,3
        Der(m) = Psi%GetDer(m)
    enddo

    Term1=NewSpinor()
    do m=1,3
        Term1 = Term1 + Pauli(Der(m),m)
    enddo
    Term1 = DivDpot(:,:,:,it) * Term1

    Term2=NewSpinor()
    do m=1,3
        do n=1,3
            Term2 = Term2 + DerDPot(:,:,:,n,m,it)*Pauli(Der(m),n)
        enddo
    enddo

    Term3=NewSpinor()

    call DeriveSpinor(Der(1), DPsi(1:3,1),-Psi%Parity,-Psi%Signature, Psi%TimeSimplex)
    call DeriveSpinor(Der(2), DPsi(1:3,2),-Psi%Parity,-Psi%Signature,-Psi%TimeSimplex)
    call DeriveSpinor(Der(3), DPsi(1:3,3),-Psi%Parity, Psi%Signature, Psi%TimeSimplex)

    ! Slightly wasteful usage of CPU time here. We overwrite the previously
    ! calculated second derivate d_mu d_mu since it is pretty inaccurate
    ! when not represented correctly. If CPU time is important, this can be commented
    ! out.
!     DPsi(1,1) = SecondDerivativeSpinor(Psi%Value,1,-Psi%Parity,-Psi%Signature, Psi%TimeSimplex)
!     DPsi(2,2) = SecondDerivativeSpinor(Psi%Value,2,-Psi%Parity,-Psi%Signature,-Psi%TimeSimplex)
!     DPsi(3,3) = SecondDerivativeSpinor(Psi%Value,3,-Psi%Parity, Psi%Signature, Psi%TimeSimplex)

    do m=1,3
        do n=1,3
            Term3 = Term3 + Dpot(:,:,:,m,it) * Pauli(DPsi(m,n) ,n)
        enddo
    enddo

    ActionOfD = -0.5_dp*(Term1 + Term2 + 2.0_dp*Term3)

  end function ActionOfD

  function ActionOfC(Psi)
    !---------------------------------------------------------------------------
    ! Function that computes the action of the C potential
    !  ActionOfC = - \nabla \cdot [\sigma \cdot C_q ] \nabla
    !---------------------------------------------------------------------------
    ! Notice that there are two ways of calculating this, which numerically
    ! are NOT EQUIVALENT.
    !
    ! 1) Calculate
    !       \sigma \cdot C \nabla \Psi.
    !    Then apply \nabla on that.
    ! 2) Calculate
    !       \sigma \cdot C \Delta \Psi + (\sigma \cdot \nabla C) \nabla Psi
    !
    ! We use the first formulation, as this is the one CR8 uses.
    !
    ! Why are they not equivalent?
    ! The first reason is that the \nabla and \Delta operators usually are not
    ! represented to the same order in finite difference formulas (respectively
    ! order 3 and order 4). This is likely the biggest source of the difference.
    !
    !---------------------------------------------------------------------------
    use Derivatives

    integer                 :: m, it, n,l
    type(Spwf) , intent(in) :: Psi
    type(Spinor)            :: ActionOfC, Der(3), Lap, Temp(3)

    it    = (Psi%GetIsospin() + 3)/2
    ActionOfC = NewSPinor()

    do n=1,3
        Der(n)  = Psi%GetDer(n)
        Temp(n) = NewSpinor()
        do m=1,3
            Temp(n) = Temp(n) - CPot(:,:,:,m,it)*Pauli(Der(n),m)
        enddo
    enddo

    !---------------------------------------------------------------------------
    ! Small explanation on the symmetries here, as I spent quite some time
    ! getting this right.
    !
    ! 1) Notice that \sum_m C_m \sigma_m is a scalar, having +1 for all quantum
    !    numbers. (Because \sigma_m is the generator of C_m ofcourse.)
    ! 2) There are derivatives already present here.
    !    \nabla_x , \nabla_y , \nabla_z all anticommute with parity
    !    \nabla_x and \nabla_y anticommute with signature
    !    \nabla_z commutes with signature
    !    \nabla_y anticommutes with Time-Simplex
    !    \nabla_x and \nabla_z commute with TimeSimplex
    !---------------------------------------------------------------------------

    do l=1,4
        Temp(1)%Grid(:,:,:,l,1) =                                                     &
    &   DeriveX(Temp(1)%Grid(:,:,:,l,1), -Psi%Parity,-Psi%Signature, Psi%TimeSimplex,l)
        Temp(2)%Grid(:,:,:,l,1) =                                                     &
    &   DeriveY(Temp(2)%Grid(:,:,:,l,1), -Psi%Parity,-Psi%Signature,-Psi%TimeSimplex,l)
        Temp(3)%Grid(:,:,:,l,1) =                                                     &
    &   DeriveZ(Temp(3)%Grid(:,:,:,l,1), -Psi%Parity, Psi%Signature, Psi%TimeSimplex,l)
    enddo

    ActionOfC = ActionOfC + Temp(1) + Temp(2) + Temp(3)

   end function ActionOfC
   
   subroutine calcTfield()
        !-----------------------------------------------------------------------
        ! Calculates the field associated with real and imaginary parts of the 
        ! Tmnk N2LO density.
        !
        ! Re Ft_mnk= 2 C^(4 Ms)_t T_k delta_mn + 4 C^(4 Ms) Re T_mnk 
        !          - 2 C^(4 Ms) D_m D_n s_k
        !
        ! Im Ft_mnk= 4 C^(4 Ms) Im T_mnk 
        !
        ! This routine also calculates the derivative of Re Ft_mnk for use in 
        ! the actionofT routine.
        ! 
        ! ReDTField = sum_mu D_mu F_munuka
        !-----------------------------------------------------------------------
        use Derivatives
        integer :: it, mu
        
        ! The imaginary field is time-even
        do it=1,2
            ImTfield(:,:,:,:,:,it)=  4*N2ImTmn(1)*sum(Density%ImDTN2LO,6)      &
            &                       +4*N2ImTmn(2)*Density%ImDTN2LO(:,:,:,:,:,it)
        enddo

        if(TRC) return
        
        ! The real field is time-odd
        do it=1,2
            ReTfield(:,:,:,:,:,:,it) = &
            &                      4*N2ReTmn(1)*sum(Density%ReKN2LO,7)         &
            &                    + 4*N2ReTmn(2)*Density%ReKN2LO(:,:,:,:,:,:,it)&
            &                    - 2*N2TmnD2S(1)*sum(Density%D2S,7)            &
            &                    - 2*N2TmnD2S(2)*Density%D2S(:,:,:,:,:,:,it)     
            do mu=1,3
                ReTField(:,:,:,mu,mu,:,it)= 2*N2vecT(1)*sum(Density%VecT,5)    &
                &                          +2*N2vecT(2)*Density%VecT(:,:,:,:,it)
            enddo
            !-------------------------------------------------------------------
            ! Deriving the T-field
            ! ReDTField = sum_mu D_mu F_munuka
            !-------------------------------------------------------------------
            ! kappa = 1, nu = 1, sum over mu
            ReDTField(:,:,:,1,1,it) =                                          &
            & DeriveX(ReTfield(:,:,:,1,1,1,it), parityint, -signatureint, timesimplexint,1)
            ReDTField(:,:,:,1,1,it) = ReDTField(:,:,:,1,1,it) +                &
            & DeriveY(ReTfield(:,:,:,2,1,1,it), ParityInt, -SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,1,1,it) = ReDTField(:,:,:,1,1,it) +                &
            & DeriveX(ReTfield(:,:,:,3,1,1,it),-ParityInt,  SignatureInt, TimeSimplexInt,1)
            !-------------------------------------------------------------------
            ! kappa = 2, nu = 1, sum over mu
            ReDTField(:,:,:,1,2,it) =                                     &
            & DeriveX(ReTField(:,:,:,1,1,2,it), ParityInt, -SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,1,2,it) = ReDTField(:,:,:,1,2,it) +                                                      &
            & DeriveY(ReTField(:,:,:,2,1,2,it), ParityInt, -SignatureInt, TimeSimplexInt,1)
            ReDTField(:,:,:,1,2,it) = ReDTField(:,:,:,1,2,it) +                                                      &
            & DeriveZ(ReTField(:,:,:,3,1,2,it), ParityInt,  SignatureInt, TimeSimplexInt,2)
            !-------------------------------------------------------------------
            ! kappa = 3, nu = 1, sum over mu
            ReDTField(:,:,:,1,3,it) =                                     &
            & DeriveX(ReTfield(:,:,:,1,1,3,it), ParityInt,  SignatureInt, TimeSimplexInt,1)
            ReDTField(:,:,:,1,3,it) = ReDTField(:,:,:,1,3,it) +                                                      &
            & DeriveY(ReTfield(:,:,:,2,1,3,it), ParityInt,  SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,1,3,it) = ReDTField(:,:,:,1,3,it) +                                                      &
            & DeriveZ(ReTfield(:,:,:,3,1,3,it), ParityInt, -SignatureInt, TimeSimplexInt,1)
            !-------------------------------------------------------------------
            ! kappa = 1, nu = 2, sum over mu
            ReDTField(:,:,:,2,1,it)  =                                         &
            & DeriveX(ReTfield(:,:,:,1,2,1,it), ParityInt, -SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,2,1,it)  = ReDTField(:,:,:,2,1,it) +               &
            & DeriveY(ReTfield(:,:,:,2,2,1,it), ParityInt, -SignatureInt, TimeSimplexInt,1)
            ReDTField(:,:,:,2,1,it) = ReDTField(:,:,:,2,1,it) +                &
            & DeriveZ(ReTfield(:,:,:,3,2,1,it), ParityInt,  SignatureInt, TimeSimplexInt,2)
            !-------------------------------------------------------------------
            ! kappa = 2, nu = 2, sum over mu
            ReDTField(:,:,:,2,2,it)  =                                         &
            & DeriveX(ReTfield(:,:,:,1,2,2,it), ParityInt, -SignatureInt, TimeSimplexInt,1)
            ReDTField(:,:,:,2,2,it) = ReDTField(:,:,:,2,2,it) +                &
            & DeriveY(ReTfield(:,:,:,2,2,2,it), ParityInt, -SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,2,2,it) = ReDTField(:,:,:,2,2,it) +                &
            & DeriveZ(ReTfield(:,:,:,3,2,2,it), ParityInt,  SignatureInt, TimeSimplexInt,1)
            !-------------------------------------------------------------------
            ! kappa = 3, nu = 2, sum over mu 
            ReDTField(:,:,:,2,3,it) =                                          &
            & DeriveX(ReTfield(:,:,:,1,2,3,it), ParityInt, -SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,2,3,it) = ReDTField(:,:,:,2,3,it) +                &
            & DeriveY(ReTfield(:,:,:,2,2,3,it), ParityInt, -SignatureInt, TimeSimplexInt,1)
            ReDTField(:,:,:,2,3,it) = ReDTField(:,:,:,2,3,it) +                &
            & DeriveZ(ReTfield(:,:,:,3,2,3,it), ParityInt,  SignatureInt, TimeSimplexInt,2)
            !-------------------------------------------------------------------
            ! kappa= 1, nu = 3, sum over mu
            ReDTField(:,:,:,3,1,it) =                                          &
            & DeriveX(ReTfield(:,:,:,1,3,1,it), ParityInt,  SignatureInt, TimeSimplexInt,1)
            ReDTField(:,:,:,3,1,it) = ReDTField(:,:,:,3,1,it) +                &
            & DeriveY(ReTfield(:,:,:,2,3,1,it), ParityInt,  SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,3,1,it) = ReDTField(:,:,:,3,1,it) +                &
            & DeriveZ(ReTfield(:,:,:,3,3,1,it), ParityInt, -SignatureInt, TimeSimplexInt,1)
            !-------------------------------------------------------------------
            ! kappa = 2, nu = 3
            ReDTField(:,:,:,3,2,it) =                                          &
            & DeriveX(ReTfield(:,:,:,1,3,2,it), ParityInt,  SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,3,2,it) = ReDTField(:,:,:,3,2,it) +                &
            & DeriveY(ReTfield(:,:,:,2,3,2,it), ParityInt,  SignatureInt, TimeSimplexInt,1)
            ReDTField(:,:,:,3,2,it) = ReDTField(:,:,:,3,2,it) +                &
            & DeriveZ(ReTfield(:,:,:,3,3,2,it), ParityInt, -SignatureInt, TimeSimplexInt,2)
            !-------------------------------------------------------------------
            ! (nu,ka) = (3,3)
            ReDTField(:,:,:,3,3,it) =                                          &
            & DeriveX(ReTfield(:,:,:,1,3,3,it), ParityInt, -SignatureInt, TimeSimplexInt,1)
            ReDTField(:,:,:,3,3,it) = ReDTField(:,:,:,3,3,it) +                &
            & DeriveY(ReTfield(:,:,:,2,3,3,it), ParityInt, -SignatureInt, TimeSimplexInt,2)
            ReDTField(:,:,:,3,3,it) = ReDTField(:,:,:,3,3,it) +                &
            & DeriveZ(ReTfield(:,:,:,3,3,3,it), ParityInt,  SignatureInt, TimeSimplexInt,1)
        enddo
   end subroutine calcTfield
   
   function ActionOfImTField(Psi) result(ActionOfT)
        !-----------------------------------------------------------------------
        ! Action of the imaginary T-field.
        !
        ! - i [ D_m Im Ft_mnk ] D_n sigma_k
        ! 
        !-----------------------------------------------------------------------
        type(Spwf), intent(in) :: Psi
        type(Spinor)           :: ActionOfT
        integer :: nu, kappa, it
        
        it    = (Psi%GetIsospin() + 3)/2
        ActionOfT = newspinor()
        
        do kappa=1,3
            do nu=1,3
                ActionOfT = ActionofT -                                        &
                &          ImTfield(:,:,:,nu,kappa,it)*Pauli(Psi%Der(nu), kappa)
            enddo
        enddo
        ActionOfT = MultiplyI(ActionOfT) 
   end function ActionOfImTField

   function ActionOfReTField(Psi) result(ActionOfT)
        !-----------------------------------------------------------------------
        ! Action of the real part of the T-field.
        !
        ! - D_m Re Ft_mnk D_n sigma_k = - [ D_m Re Ft_mnk] D_n     sigma_k
        !                               -       Re Ft_mnk  D_m D_n sigma_k
        !-----------------------------------------------------------------------
        type(Spwf), intent(in) :: Psi
        type(Spinor)           :: ActionOfT
        integer :: nu, kappa, it, mu
        
        it    = (Psi%GetIsospin() + 3)/2
        ActionOfT = newspinor()
        
        do kappa=1,3
            do nu=1,3
                do mu=1,3
                    ActionOfT = ActionOfT - ReTfield(:,:,:,mu,nu,kappa,it) *   &
                    &                       Pauli(Psi%SecondDer(mu,nu), kappa)
                enddo
                ActionOfT = ActionOfT -                                        &
                &       ReDTfield(:,:,:,nu,kappa,it) * Pauli(Psi%Der(nu), kappa) 
            enddo
        enddo
        
   end function ActionOfReTField

   subroutine CalcPiField()
    !---------------------------------------------------------------------------
    ! Calculate the field associated with the N2LO density Pi
    !
    ! Fpi_m = 4 C^(4 Mrho) j_m 
    !---------------------------------------------------------------------------
    
    integer :: mu, it, at
    
    do it=1,2
        at = 3 - it
        PiField(:,:,:,:,it) = 4*(N2jpi(1)+N2jpi(2))*Density%vecj(:,:,:,:,it)   &
        &                   + 4*N2jpi(1)           *Density%vecj(:,:,:,:,at) 
    enddo
   end subroutine CalcPiField

   function ActionOfPi(Psi)
     !-----------------------------------------------------------------------
     ! Action of the Pi field
     !
     ! + i D_n Fpi_m D_n D_m
     !-----------------------------------------------------------------------
     type(Spwf), intent(in) :: Psi
     type(Spinor)           :: ActionOfPi, temp, ActionofPiX, ActionOfPiY, ActionOfPiZ
     integer                :: nu, kappa, it, mu
   
     it = (Psi%GetIsospin() + 3)/2
     
     ActionOfPi = NewSpinor()
     
     temp=NewSpinor()
     do mu=1,3
        temp = temp + PiField(:,:,:,mu,it) * Psi%SecondDer(mu, 1)
     enddo  
     call DeriveSpinor_X(temp, ActionofPiX,-Psi%Parity,-Psi%Signature, Psi%TimeSimplex)
     
     temp=NewSpinor()
     do mu=1,3
        temp = temp + PiField(:,:,:,mu,it) * Psi%SecondDer(mu, 2)
     enddo  
     call DeriveSpinor_Y(temp, ActionofPiY,-Psi%Parity,-Psi%Signature,-Psi%TimeSimplex)
     
     temp=NewSpinor()
     do mu=1,3
        temp = temp + PiField(:,:,:,mu,it) * Psi%SecondDer(mu, 3)
     enddo  
     call DeriveSpinor_Z(temp, ActionOfPiZ,-Psi%Parity, Psi%Signature, Psi%TimeSimplex)
         
     ActionOfPi = ActionOfPiX + ActionOfPiY + ActionOfPiZ
   end function ActionofPi
end module MeanFields
