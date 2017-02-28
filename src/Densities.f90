module Densities
!------------------------------------------------------------------------------
! Conventions are in line (as much as possible)
!                        with P. Bonche et al., Nucl. Phys. A467 (1987), 115-135
!                        and V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
!-------------------------------------------------------------------------------

  use CompilationInfo
  use Geninfo
  use Derivatives

  implicit none

  !---------------------------------------------------------------------------
  ! Damping parameter: only used when not doing Pulay or Broyden
  real(KIND=dp), save :: DampingParam=0.75_dp
  !---------------------------------------------------------------------------
  ! Change of the density at each iteration, integrated and divided by
  ! 1 - Damping
  real(KIND=dp),  save :: DensityChange=0.0_dp
  !---------------------------------------------------------------------------
  ! RhoIn & RhoOut: DIIS mixing densities...
  real(KIND=dp),  allocatable :: RhoIn(:), RhoOut(:)
  !-----------------------------------------------------------------------------
  ! Density mixing scheme info: order etc.
  integer                          :: PulayOrder=0, MixingScheme=0
  !-----------------------------------------------------------------------------
  ! Flag to determine whether MOCCa calculates the densities at iteration 0
  ! or not.
  logical :: Recalc=.false.

  type DensityVector
    !---------------------------------------------------------------------------
    ! Custom type for containing all densities.
    ! Very handy for mixing densities.
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Normal density and kinetic density and their derivatives
    real(KIND=dp), allocatable :: Rho(:,:,:,:), Tau(:,:,:,:)
    real(KIND=dp), allocatable :: DerRho(:,:,:,:,:),LapRho(:,:,:,:)
    !---------------------------------------------------------------------------
    ! Current density
    real(KIND=dp), allocatable :: Vecj(:,:,:,:,:), divvecj(:,:,:,:)
    real(KIND=dp), allocatable :: Rotvecj(:,:,:,:,:)
    !---------------------------------------------------------------------------
    ! Spin density and its rotor, laplacian and divergence.
    real(KIND=dp), allocatable :: Vecs(:,:,:,:,:), RotS(:,:,:,:,:)
    real(KIND=dp), allocatable :: LapS(:,:,:,:,:), DivS(:,:,:,:)
    real(KIND=dp), allocatable :: GradDivS(:,:,:,:,:)
    real(KIND=dp), allocatable :: DerS(:,:,:,:,:,:)
    !---------------------------------------------------------------------------
    ! \vec{T} density, Kinetic Vector Density
    real(KIND=dp), allocatable :: VecT(:,:,:,:,:)
    !---------------------------------------------------------------------------
    ! \vec{F} density
    real(KIND=dp), allocatable :: VecF(:,:,:,:,:)
    !---------------------------------------------------------------------------
    ! J_{\mu,\nu} Tensor density and the divergence of its vector component
    real(KIND=dp), allocatable :: JMuNu(:,:,:,:,:,:),NablaJ(:,:,:,:)
    !---------------------------------------------------------------------------
    ! N2LO densities needed.
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable ::  RtauN2LO(:,:,:,:,:,:)
    real(KIND=dp), allocatable ::  ItauN2LO(:,:,:,:,:,:)
    real(KIND=dp), allocatable ::   ReKN2LO(:,:,:,:,:,:,:)
    real(KIND=dp), allocatable ::   ImKN2LO(:,:,:,:,:,:,:)
    real(KIND=dp), allocatable :: LapLapRho(:,:,:,:)
    real(KIND=dp), allocatable ::    D2RTau(:,:,:,:)
    real(KIND=dp), allocatable ::    PiN2LO(:,:,:,:,:)
    real(KIND=dp), allocatable ::     VN2LO(:,:,:,:,:,:)
    real(KIND=dp), allocatable ::     QN2LO(:,:,:,:)
    real(KIND=dp), allocatable ::     SN2LO(:,:,:,:,:)
    real(KIND=dp), allocatable ::     D2Rho(:,:,:,:,:,:)
    real(KIND=dp), allocatable ::       D2S(:,:,:,:,:,:,:)
    real(KIND=dp), allocatable ::    DJmunu(:,:,:,:,:,:,:)
    
  end type DensityVector

  !-----------------------------------------------------------------------------
  ! Keeping the history of the densities, useful for mixing!
  type(DensityVector),  allocatable :: DensityHistory(:)

  !-----------------------------------------------------------------------------
  ! Current density vector.
  type(DensityVector),public,  save               :: Density

  ! Constructor interface
  interface DensityVector
    module procedure NewDensityVector
  end interface

  interface operator (+)
    !Overloading "+" to be used to add densities.
    module procedure Add
  end interface

  interface operator (*)
    module procedure MultiplyScalar
    module procedure Inproduct
  end interface
contains

  pure function NewDensityVector(newnx,newny,newnz) result(R)
    !---------------------------------------------------------------------------
    ! Constructor for the DensityVector class.
    ! The sizes are optional, if not supplied nx,ny,nz will be used.
    !---------------------------------------------------------------------------
    use Force

    integer, intent(in),optional :: newnx,newny,newnz
    integer                      :: sizex, sizey, sizez
    
    type(DensityVector) :: R

    if(present(newnx)) then
       sizex = newnx ; sizey = newny ; sizez = newnz
    else
      sizex = nx ; sizey = ny ; sizez = nz
    endif
   
    ! Time-even densities!
    allocate(R%Rho(sizex,sizey,sizez,2), R%DerRho(sizex,sizey,sizez,3,2))
    allocate(R%LapRho(sizex,sizey,sizez,2))
    allocate(R%Tau(sizex,sizey,sizez,2),R%NablaJ(sizex,sizey,sizez,2))
    R%Rho    =0.0_dp ; R%LapRho=0.0_dp ; R%DerRho=0.0_dp
    R%Tau    =0.0_dp ; R%NablaJ =0.0_dp

    if((B14.ne.0.0_dp .or. B15.ne.0.0_dp .or. B16.ne.0.0_dp.or.B17.ne.0.0_dp)  & 
    & .or. t1n2 .ne. 0.0_dp) then
      !Time-even density with Jmunu!
      allocate(R%JMuNu(sizex,sizey,sizez,3,3,2))
      R%Jmunu = 0.0_dp
    endif

    !Time-odd densities
    if(.not.TRC) then
      allocate(R%Vecj(sizex,sizey,sizez,3,2),R%RotS(sizex,sizey,sizez,3,2))
      allocate(R%Vecs(sizex,sizey,sizez,3,2),R%Rotvecj(sizex,sizey,sizez,3,2))
      allocate(R%DerS(sizex,sizey,sizez,3,3,2))
      R%Vecj=0.0_dp ; R%RotS = 0.0_dp
      R%Vecs=0.0_dp ; R%RotVecj =0.0_dp ; R%DerS=0.0_dp ;

      !Derivatives of s
      if((B18 .ne. 0.0_dp .or. B19 .ne. 0.0_dp .or. &
      &   B20.ne.0.0_dp .or. B21.ne.0.0_dp)    .or. &
      &   t1n2.ne.0.0_dp) then
        allocate(R%LapS(sizex,sizey,sizez,3,2))    ;   R%LapS=0.0_dp ;
        allocate(R%GradDivS(sizex,sizey,sizez,3,2));   R%GradDivS=0.0_dp
        allocate(R%DivS(sizex,sizey,sizez,2))      ;   R%DivS = 0.0_dp
      endif

      !F & T densities.
      if((B14.ne.0.0_dp .or. B15.ne.0.0_dp).or.(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp)) then
        allocate(R%VecT(sizex,sizey,sizez,3,2))
        R%VecT=0.0_dp ;
      endif
      if((B16.ne.0.0_dp .or. B17.ne.0.0_dp).or.(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp)) then
        allocate(R%VecF(sizex,sizey,sizez,3,2))
        R%VecF=0.0_dp ;
      endif
    endif

    if(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp) then
        allocate(R%RtauN2LO(sizex,sizey,sizez,3,3,2)) ; R%RtauN2LO = 0.0_dp
        allocate(R%ItauN2LO(sizex,sizey,sizez,3,3,2)) ; R%ItauN2LO = 0.0_dp
        allocate(R%ReKN2LO(sizex,sizey,sizez,3,3,3,2)); R%ReKN2LO   = 0.0_dp
        allocate(R%ImKN2LO(sizex,sizey,sizez,3,3,3,2)); R%ImKN2LO   = 0.0_dp
        allocate(R%PiN2LO(sizex,sizey,sizez,3,2))     ; R%PiN2LO  = 0.0_dp
        allocate(R%VN2LO(sizex,sizey,sizez,3,3,2))    ; R%VN2LO  = 0.0_dp
        allocate(R%QN2LO(sizex,sizey,sizez,2))        ; R%QN2LO  = 0.0_dp
        allocate(R%SN2LO(sizex,sizey,sizez,3,2))      ; R%SN2LO  = 0.0_dp
        allocate(R%D2Rho   (sizex,sizey,sizez,3,3,2)) ; R%D2Rho     = 0.0_dp
        allocate(R%D2S   (sizex,sizey,sizez,3,3,3,2)) ; R%D2S       = 0.0_dp
        allocate(R%DJmunu(sizex,sizey,sizez,3,3,3,2)) ; R%DJmunu  = 0.0_dp
        allocate(R%LapLapRho(sizex,sizey,sizez,2))    ; R%LapLapRho = 0.0_dp
        allocate(R%D2Rtau(sizex,sizey,sizez,2))       ; R%D2Rtau     = 0.0_dp
        allocate(R%divvecj (sizex,sizey,sizez,2))     ; R%divvecj   = 0.0_dp
    endif

  end function NewDensityVector

  function Add(Den1, Den2) result(Sum)
    !---------------------------------------------------------------------------
    ! Function for overloading the + operator for density vectors.
    !
    !---------------------------------------------------------------------------

    type(DensityVector), intent(in) :: Den1, Den2
    type(DensityVector)             :: Sum

    Sum = NewDensityVector()

    Sum%Rho   = Den1%Rho    + Den2%Rho
    Sum%Tau   = Den1%Tau    + Den2%Tau
    Sum%NablaJ= Den1%NablaJ + Den2%NablaJ
    Sum%DerRho= Den1%DerRho + Den2%DerRho
    Sum%LapRho= Den1%LapRho + Den2%LapRho

    if(allocated(Sum%Jmunu)) then
      Sum%JMuNu = Den1%Jmunu + Den2%JMunu
    endif
    
    if(allocated(Sum%RtauN2LO)) then
        Sum%RtauN2LO   = Den1%RTauN2LO     + Den2%RTauN2LO
        Sum%ItauN2LO   = Den1%ITauN2LO     + Den2%ITauN2LO
        Sum%D2RTau     = Den1%D2RTau       + Den2%D2RTau
        Sum%D2Rho      = Den1%D2Rho        + Den2%D2Rho
        Sum%LapLapRho  = Den1%LapLapRho    + Den2%LapLapRho
        Sum%DJmunu     = Den1%DJmunu       + Den2%DJmunu
        Sum%ReKN2LO    = Den1%ReKN2LO      + Den2%ReKN2LO
        Sum%ImKN2LO    = Den1%ImKN2LO      + Den2%ImKN2LO
        Sum%PiN2LO     = Den1%PiN2LO       + Den2%PiN2LO
        Sum%VN2LO      = Den1%VN2LO        + Den2%VN2LO
        Sum%QN2LO      = Den1%QN2LO        + Den2%QN2LO
        Sum%SN2LO      = Den1%SN2LO        + Den2%SN2LO
        Sum%D2S        = Den1%D2S          + Den2%D2S
    endif
    
    if(.not.TRC) then
      Sum%Vecj   = Den1%VecJ     + Den2%VecJ
      Sum%RotVecj= Den1%RotVecJ  + Den2%RotVecJ
      Sum%Vecs   = Den1%vecs     + Den2%vecs
      Sum%RotS   = Den1%RotS     + Den2%RotS
      if(allocated(Sum%LapS)) then
        Sum%LapS = Den1%LapS + Den2%LapS
        Sum%DerS = Den1%DerS + Den2%DerS
        Sum%GradDivS = Den1%GradDivS + Den2%GradDivS
        Sum%DivS = Den1%DivS + Den2%DivS
      endif
    endif

    if(allocated(Sum%VecT)) then
      Sum%VecT = Den1%VecT + Den2%VecT
    endif

    if(allocated(Sum%VecF)) then
      Sum%VecF = Den1%VecF + Den2%VecF
    endif

  end function Add

  function MultiplyScalar(A,Den) result(Prod)
    !---------------------------------------------------------------------------
    ! Function for overloading the * operator for density vectors.
    !
    !---------------------------------------------------------------------------

    type(DensityVector), intent(in) :: Den
    type(DensityVector)             :: Prod
    real(KIND=dp), intent(in)       :: A

    Prod = NewDensityVector()

    Prod%Rho   = A*Den%Rho
    Prod%Tau   = A*Den%Tau
    Prod%NablaJ= A*Den%NablaJ
    Prod%DerRho= A*Den%DerRho
    Prod%LapRho= A*Den%LapRho
    
    if(allocated(Prod%RtauN2LO)) then
        Prod%RtauN2LO   = A*Den%RTauN2LO
        Prod%ItauN2LO   = A*Den%ITauN2LO
        Prod%D2RTau     = A*Den%D2RTau
        Prod%LapLapRho  = A*Den%LapLapRho
        Prod%D2Rho      = A*Den%D2Rho
        Prod%ReKN2LO    = A*Den%ReKN2LO 
        Prod%ImKN2LO    = A*Den%ImKN2LO 
        Prod%PiN2LO     = A*Den%PiN2LO  
        Prod%VN2LO      = A*Den%VN2LO   
        Prod%QN2LO      = A*Den%QN2LO   
        Prod%SN2LO      = A*Den%SN2LO
        Prod%D2S        = A*Den%D2S 
        Prod%DJmunu     = A*Den%DJmunu 
    endif

    if(allocated(Prod%Jmunu)) then
      Prod%JMuNu = A*Den%Jmunu
    endif
    if(.not.TRC) then
      Prod%Vecj   = A*Den%VecJ
      Prod%RotVecj= A*Den%RotVecJ
      Prod%Vecs   = A*Den%vecs
      Prod%RotS   = A*Den%RotS

      if(allocated(Prod%LapS)) then
        Prod%LapS     = A*Den%LapS
        Prod%GradDivS = A*Den%GradDivS
        Prod%DivS     = A*Den%DivS
        Prod%DerS     = A*Den%DerS
      endif
    endif

    if(allocated(Prod%VecT)) then
      Prod%VecT = A*Den%VecT
    endif

    if(allocated(Prod%VecF)) then
      Prod%VecF = A*Den%VecF
    endif

  end function MultiplyScalar

  real(KIND=dp) function Inproduct(Den1, Den2)
    !---------------------------------------------------------------------------
    ! Inproduct for overloading * with densities.
    !
    !---------------------------------------------------------------------------
    type(DensityVector), intent(in) :: Den1, Den2

    Inproduct = 0.0_dp
    Inproduct = sum(Den1%Rho * Den2%Rho) + sum(Den1%Tau * Den2%Tau)
    if(allocated(Den1%JMuNu)) then
      Inproduct = Inproduct + sum(Den1%JMuNu * Den2%JMuNu)
    else
      Inproduct = Inproduct + sum(Den1%NablaJ* Den2%NablaJ)
    endif
    if(allocated(Den1%Vecj)) then
      Inproduct = Inproduct + sum(Den1%VecJ * Den2%VecJ)
    endif
    if(allocated(Den1%Vecs)) then
      Inproduct = Inproduct + sum(Den1%Vecs * Den2%Vecs)
    endif
    if(allocated(Den1%VecF)) then
      Inproduct = Inproduct + sum(Den1%VecF * Den2%VecF)
    endif
    if(allocated(Den1%VecT)) then
      Inproduct = Inproduct + sum(Den1%VecT * Den2%VecT)
    endif
  end function Inproduct

  subroutine ReadDensitInfo
    !---------------------------------------------------------------------------
    !Reading the relevant info from input.
    !---------------------------------------------------------------------------

    NameList /densit/ DampingParam, PulayOrder, MixingScheme, Recalc

    read(unit=*, NML=densit)
    if(MixingScheme.ne.0 .and. PulayOrder.lt.1) then
      call stp('PulayOrder should be strictly positive!')
    endif

  end subroutine ReadDensitInfo

  subroutine IniDen
    !---------------------------------------------------------------------------
    ! This routine allocates the necessary memory space and puts all the
    ! densities equal to zero.
    !---------------------------------------------------------------------------

    integer :: i

    !Always keep at least one history.
    allocate(DensityHistory(max(2,PulayOrder)))
    do i=1, max(2,PulayOrder)
      DensityHistory(1) = NewDensityVector()
    enddo

    !Initialise the standard densities
    Density = NewDensityVector()

  end subroutine IniDen

  subroutine UpdateDensities(NoSave,OnlyRho)
    !---------------------------------------------------------------------------
    ! Master routine for computing and mixing the densities at a new iteration.
    !
    !---------------------------------------------------------------------------
    integer             :: i
    logical, intent(in), optional :: OnlyRho
    integer, intent(in) :: NoSave
    !---------------------------------------------------------------------------
    !Move the history correctly
    if(NoSave.ne.1) then
      do i=1,max(2,PulayOrder)-1
        DensityHistory(max(2,PulayOrder)-i+1)=DensityHistory(max(2,PulayOrder)-i)
      enddo
      DensityHistory(1) = Density
    endif
    !---------------------------------------------------------------------------
    !Compute the new density.
    if(present(OnlyRho)) then
      call ComputeDensity(Density,.true.)
    else
      call ComputeDensity(Density,.false.)
    endif

  end subroutine UpdateDensities

  subroutine ComputeDensity(DenIn, OnlyRho, Response)
    !---------------------------------------------------------------------------
    ! Compute all the different densities from the curent wavefunctions.
    ! If Response is set and true, this routine computes the Lipkin-Nogami
    ! response densities instead of normal ones.
    !---------------------------------------------------------------------------

    use SpwfStorage, only : DensityBasis, nwt
    use Force
    use Derivatives, only : Laplacian, DeriveX,DeriveY,DeriveZ

    type(Densityvector), intent(inout) :: Denin
    integer                            :: it, i,k
    real(KIND=dp)                      :: Occupation
    logical, intent(in)                :: OnlyRho
    logical, intent(in), optional      :: Response
    logical                            :: moreders
    real(KIND=dp) :: DerVecj(nx,ny,nz,3,3), DivS(nx,ny,nz,2), temp(nx,ny,nz)
    
    !Reset the values of the densities.
    if(OnlyRho) then
      DenIn%Rho = 0.0_dp
    else
        DenIn = NewDensityVector()
    endif

    if(.not.associated(DensityBasis)) then
      call stp('Density basis not correctly associated at the start of '       &
      &        // 'subroutine ComputeDensity.')
    endif

    !---------------------------------------------------------------------------
    !Compute all the densities to sum.
    do i=1,nwt
      Occupation = DensityBasis(i)%GetOcc()

      ! For the Lipkin-Nogami response densities, the expressions are the
      ! same, only the weights are different.
      ! Now the weights are 2*u^2 v^2 = 2*(1 - v^2 ) v^2
      if(present(Response)) then
        if(Response) Occupation = 2*Occupation * (1 - Occupation)
      endif
      it = (DensityBasis(i)%GetIsospin() + 3 )/2

      DenIn%Rho(:,:,:,it) = DenIn%Rho(:,:,:,it) +                              &
      &                     Occupation * DensityBasis(i)%GetDensity()
      if(OnlyRho) cycle

      DenIn%Tau(:,:,:,it) = DenIn%Tau(:,:,:,it) +                              &
      &                     Occupation * DensityBasis(i)%GetTau()

      !Decide wether we need Jmunu completely, or only NablaJ
      if((B14.ne.0.0_dp.or.B15.ne.0.0_dp.or.B16.ne.0.0_dp .or.B17.ne.0.0_dp)   &
      &  .or. (t1n2 .ne. 0.0_dp)) then
        DenIn%Jmunu(:,:,:,:,:,it) = DenIn%Jmunu(:,:,:,:,:,it) +                &
        &                           Occupation * DensityBasis(i)%GetJmunu()
        !Temporary
        DenIn%NablaJ(:,:,:,it)    = DenIn%NablaJ(:,:,:,it) +                   &
        &                           Occupation * DensityBasis(i)%GetNablaJ()
      else
        DenIn%NablaJ(:,:,:,it)    = DenIn%NablaJ(:,:,:,it) +                   &
        &                           Occupation * DensityBasis(i)%GetNablaJ()
      endif

      !Only compute time-odd densities when Time Reversal is broken
      if(.not.TRC) then
        DenIn%Vecj(:,:,:,:,it)= DenIn%Vecj(:,:,:,:,it) +                       &
        &                       Occupation*DensityBasis(i)%GetVecj()
        DenIn%Vecs(:,:,:,:,it)= DenIN%Vecs(:,:,:,:,it) +                       &
        &                       Occupation*DensityBasis(i)%GetVecs()
        if((B14.ne.0.0_dp .or. B15.ne.0.0_dp).or.allocated(DenIn%RtauN2LO)) then
          DenIn%VecT(:,:,:,:,it)= DenIn%VecT(:,:,:,:,it) +                     &
          &                       Occupation*DensityBasis(i)%GetVecT()
        endif
        if((B16 .ne. 0.0_dp .or. B17.ne.0.0_dp).or.allocated(DenIn%RtauN2LO)) then
          DenIn%VecF(:,:,:,:,it)= DenIn%VecF(:,:,:,:,it) +                     &
          &                       Occupation*DensityBasis(i)%GetVecF()
        endif
      endif
      
      if (allocated(DenIn%RtauN2LO)) then
        ! calculate N2LO densities
        DenIn%RTauN2LO(:,:,:,:,:,it)  = DenIn%RTauN2LO(:,:,:,:,:,it) +         &
        &                              Occupation * DensityBasis(i)%GetRTauN2LO()
        DenIn%ImKN2LO(:,:,:,:,:,:,it) = DenIn%ImKN2LO(:,:,:,:,:,:,it) +        &
        &                              Occupation * DensityBasis(i)%GetImKN2LO()
        DenIn%QN2LO(:,:,:,it)         = DenIn%QN2LO(:,:,:,it) +                &
        &                              Occupation * DensityBasis(i)%GetQN2LO()
        DenIn%VN2LO(:,:,:,:,:,it)     = DenIn%VN2LO(:,:,:,:,:,it) +            &
        &                              Occupation * DensityBasis(i)%GetVN2LO()
        
        if(.not. TRC) then
            DenIn%ITauN2LO(:,:,:,:,:,it)  = DenIn%ITauN2LO(:,:,:,:,:,it) +     &
            &                         Occupation * DensityBasis(i)%GetITauN2LO()
            DenIn%ReKN2LO(:,:,:,:,:,:,it) = DenIn%ReKN2LO(:,:,:,:,:,:,it) +    &
            &                         Occupation * DensityBasis(i)%GetREKN2LO()
            DenIn%PiN2LO(:,:,:,:,it)      = DenIn%PiN2LO(:,:,:,:,it) +         &
            &                         Occupation * DensityBasis(i)%GetPiN2LO()
            DenIn%SN2LO(:,:,:,:,it)       = DenIn%SN2LO(:,:,:,:,it) +          &
            &                         Occupation * DensityBasis(i)%GetSN2LO()     
        endif
      endif
  
    enddo
    if(all(DenIn%Rho.eq.0.0_dp)) call stp('Rho is zero')
    if(onlyRho) return

    !---------------------------------------------------------------------------
    ! Compute all necessary derivatives of these densities.
    !Derivative and laplacian of Rho.
    do it=1,2
        DenIn%DerRho(:,:,:,1,it) = &
        & DeriveX(DenIn%Rho(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
        DenIn%DerRho(:,:,:,2,it) = &
        & DeriveY(DenIn%Rho(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
        DenIn%DerRho(:,:,:,3,it) = &
        & DeriveZ(DenIn%Rho(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
        DenIn%LapRho(:,:,:,it)   = &
        & Laplacian(DenIn%Rho(:,:,:,it),ParityInt,SignatureInt,TimeSimplexInt,1)
    enddo
    
    if (allocated(DenIn%RtauN2LO)) then
        do it=1,2
            DenIn%LapLapRho(:,:,:,it)   = &
            & Laplacian(DenIn%LapRho(:,:,:,it),ParityInt,SignatureInt,TimeSimplexInt,1)
        
            ! Full tensor of second derivatives of rho        
            DenIn%D2Rho(:,:,:,1,1,it) = &
            & DeriveX(DenIn%DerRho(:,:,:,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
            DenIn%D2Rho(:,:,:,2,1,it) = & 
            & DeriveY(DenIn%DerRho(:,:,:,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
            DenIn%D2Rho(:,:,:,3,1,it) = &
            & DeriveZ(DenIn%DerRho(:,:,:,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)

            DenIn%D2Rho(:,:,:,1,2,it) = &
            & DeriveX(DenIn%DerRho(:,:,:,2,it), -ParityInt,-SignatureInt,-TimeSimplexInt,1)
            
            DenIn%D2Rho(:,:,:,2,2,it) = &
            & DeriveY(DenIn%DerRho(:,:,:,2,it), -ParityInt,-SignatureInt,-TimeSimplexInt,1)
            DenIn%D2Rho(:,:,:,3,2,it) = &
            & DeriveZ(DenIn%DerRho(:,:,:,2,it), -ParityInt,-SignatureInt,-TimeSimplexInt,1)
            
            DenIn%D2Rho(:,:,:,1,3,it) = &
            & DeriveX(DenIn%DerRho(:,:,:,3,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
            DenIn%D2Rho(:,:,:,2,3,it) = &
            & DeriveY(DenIn%DerRho(:,:,:,3,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
            DenIn%D2Rho(:,:,:,3,3,it) = &
            & DeriveZ(DenIn%DerRho(:,:,:,3,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
            !-------------------------------------------------------------------
            ! Sum of derivatives of Tau_munu
            !  D2Rtau = sum_{mu nu} D_mu D_nu Re tau_mu_nu
            !  Note that the similar quantity of the imaginary part of tau is 
            !  never needed.
            DenIn%D2RTau(:,:,:,it) = 0.0_dp
            !-------------------------------------------------------------------
            ! Real part of tau, which is the only part that contributes to the
            ! fields
            temp =          &
            &      DeriveX(DenIn%RTauN2LO(:,:,:,1,1,it), ParityInt, SignatureInt, TimeSimplexInt,1)
            temp = temp   + &
            &      DeriveY(DenIn%RTauN2LO(:,:,:,1,2,it), ParityInt, SignatureInt, TimeSimplexInt,2)
            temp = temp   + &
            &      DeriveZ(DenIn%RTauN2LO(:,:,:,1,3,it), ParityInt,-SignatureInt, TimeSimplexInt,1)
            
            DenIn%D2RTau(:,:,:,it) =                                           &
            &           DeriveX(temp,-ParityInt,-SignatureInt, TimeSimplexInt,1)
            
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            temp =          &
            &      DeriveX(DenIn%RTauN2LO(:,:,:,2,1,it), ParityInt, SignatureInt, TimeSimplexInt,2)
            temp = temp   + &
            &      DeriveY(DenIn%RTauN2LO(:,:,:,2,2,it), ParityInt, SignatureInt, TimeSimplexInt,1)
            temp = temp   + &
            &      DeriveZ(DenIn%RTauN2LO(:,:,:,2,3,it), ParityInt,-SignatureInt, TimeSimplexInt,2)
            
            DenIn%D2RTau(:,:,:,it) =  DenIn%D2RTau(:,:,:,it) +                 &
            &           DeriveY(temp,-ParityInt,-SignatureInt, TimeSimplexInt,2)
            
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            temp =                                                             &
            &      DeriveX(DenIn%RTauN2LO(:,:,:,3,1,it), ParityInt,-SignatureInt, TimeSimplexInt,1)
            temp = temp +                                                      &
            &      DeriveY(DenIn%RTauN2LO(:,:,:,3,2,it), ParityInt,-SignatureInt, TimeSimplexInt,2)
            temp = temp +                                                      &
            &      DeriveZ(DenIn%RTauN2LO(:,:,:,3,3,it), ParityInt, SignatureInt, TimeSimplexInt,1)
            
            DenIn%D2RTau(:,:,:,it) =  DenIn%D2RTau(:,:,:,it) +                 &
            &           DeriveZ(temp,-ParityInt, SignatureInt, TimeSimplexInt,1)
            !-------------------------------------------------------------------
            ! Derivatives of nabla_mu J_munu
            !-------------------------------------------------------------------
            ! Diagonal component D_ka J_mumu
            !-------------------------------------------------------------------
            DenIn%DJmunu(:,:,:,1,1,1,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,1,1,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,2,1,1,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,1,1,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,3,1,1,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,1,1,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            
            DenIn%DJmunu(:,:,:,1,2,2,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,2,2,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,2,2,2,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,2,2,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,3,2,2,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,2,2,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            
            DenIn%DJmunu(:,:,:,1,3,3,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,3,3,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,2,3,3,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,3,3,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,3,3,3,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,3,3,it),-ParityInt,+SignatureInt,TimeSimplexInt,2)
            
            !-------------------------------------------------------------------
            ! Components of J_x^(1) => J_yz, J_zy
            !-------------------------------------------------------------------
            DenIn%DJmunu(:,:,:,1,2,3,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,2,3,it),-ParityInt,-SignatureInt,+TimeSimplexInt,1)
            DenIn%DJmunu(:,:,:,2,2,3,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,2,3,it),-ParityInt,-SignatureInt,+TimeSimplexInt,1)
            DenIn%DJmunu(:,:,:,3,2,3,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,2,3,it),-ParityInt,-SignatureInt,+TimeSimplexInt,1)
            
            DenIn%DJmunu(:,:,:,1,3,2,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,3,2,it),-ParityInt,-SignatureInt,+TimeSimplexInt,1)
            DenIn%DJmunu(:,:,:,2,3,2,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,3,2,it),-ParityInt,-SignatureInt,+TimeSimplexInt,1)
            DenIn%DJmunu(:,:,:,3,3,2,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,3,2,it),-ParityInt,-SignatureInt,+TimeSimplexInt,1)
            !-------------------------------------------------------------------
            ! Components of J_y^(1) => J_xz, J_zx
            !-------------------------------------------------------------------
            DenIn%DJmunu(:,:,:,1,3,1,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,3,1,it),-ParityInt,-SignatureInt, TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,2,3,1,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,3,1,it),-ParityInt,-SignatureInt, TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,3,3,1,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,3,1,it),-ParityInt,-SignatureInt, TimeSimplexInt,2)
            
            DenIn%DJmunu(:,:,:,1,1,3,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,1,3,it),-ParityInt,-SignatureInt, TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,2,1,3,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,1,3,it),-ParityInt,-SignatureInt, TimeSimplexInt,2)
            DenIn%DJmunu(:,:,:,3,1,3,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,1,3,it),-ParityInt,-SignatureInt, TimeSimplexInt,2)
            !-------------------------------------------------------------------
            ! Components of J_z^(1) => J_xy, J_yx
            !-------------------------------------------------------------------            
            DenIn%DJmunu(:,:,:,1,2,1,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,2,1,it),-ParityInt,+SignatureInt,+TimeSimplexInt,1)
            DenIn%DJmunu(:,:,:,2,2,1,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,2,1,it),-ParityInt,+SignatureInt,+TimeSimplexInt,1)
            DenIn%DJmunu(:,:,:,3,2,1,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,2,1,it),-ParityInt,+SignatureInt,+TimeSimplexInt,1)
            
            DenIn%DJmunu(:,:,:,1,1,2,it) =                                     &
            & DeriveX(DenIn%JMuNu(:,:,:,1,2,it),-ParityInt,+SignatureInt,+TimeSimplexInt,1)
            DenIn%DJmunu(:,:,:,2,1,2,it) =                                     &
            & DeriveY(DenIn%JMuNu(:,:,:,1,2,it),-ParityInt,+SignatureInt,+TimeSimplexInt,1)
            DenIn%DJmunu(:,:,:,3,1,2,it) =                                     &
            & DeriveZ(DenIn%JMuNu(:,:,:,1,2,it),-ParityInt,+SignatureInt,+TimeSimplexInt,1)
        enddo
    endif
    !Computing NablaJ by derivatives in the case of tensor interactions
    if(B14.ne.0.0_dp.or.B15.ne.0.0_dp.or.B17.ne.0.0_dp.or.B16.ne.0.0_dp) then
      !Temporary
      !DenIn%NablaJ = DeriveJMunu(Denin%JMuNu)
    endif

    !Computing the derivatives of vecj & vecs
    if(.not.TRC) then
        ! Derive vecj
        do it=1,2
          !See Table V in V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
          DerVecJ(:,:,:,2,1)  = &
          & DeriveY(DenIn%Vecj(:,:,:,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)
          DerVecJ(:,:,:,3,1)  = &
          & DeriveZ(DenIn%Vecj(:,:,:,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)
          DerVecJ(:,:,:,1,2)  = &
          & DeriveX(DenIn%Vecj(:,:,:,2,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
          DerVecJ(:,:,:,3,2)  = &
          & DeriveZ(DenIn%Vecj(:,:,:,2,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)
          DerVecJ(:,:,:,1,3)  = &
          & DeriveX(DenIn%Vecj(:,:,:,3,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
          DerVecJ(:,:,:,2,3)  = &
          & DeriveY(DenIn%Vecj(:,:,:,3,it), -ParityInt, SignatureInt, TimeSimplexInt,2)

          DenIn%RotVecJ(:,:,:,1,it) =  DerVecJ(:,:,:,2,3) - DerVecJ(:,:,:,3,2)
          DenIn%RotVecJ(:,:,:,2,it) =  DerVecJ(:,:,:,3,1) - DerVecJ(:,:,:,1,3)
          DenIn%RotVecJ(:,:,:,3,it) =  DerVecJ(:,:,:,1,2) - DerVecJ(:,:,:,2,1)
          
          if(allocated(DenIn%RtauN2LO)) then
            DerVecJ(:,:,:,1,1)  = &
            & DeriveX(DenIn%Vecj(:,:,:,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)
            DerVecJ(:,:,:,2,2)  = &
            & DeriveY(DenIn%Vecj(:,:,:,2,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
            DerVecJ(:,:,:,3,3)  = &
            & DeriveY(DenIn%Vecj(:,:,:,3,it), -ParityInt, SignatureInt, TimeSimplexInt,2)
            
            DenIn%divvecj(:,:,:,it) = DerVecJ(:,:,:,1,1) + &
                                    & DerVecJ(:,:,:,2,2) + &
                                    & DerVecJ(:,:,:,3,3)
          endif
        enddo

        MoreDers=.true.
        if(B18.eq.0.0_dp.and.B19.eq.0.0_dp.and.B20.eq.0.0_dp .and. B21.eq.0.0_dp) then
          MoreDers = .false.
        endif
        if(t1n2.ne.0.0_dp) then
          MoreDers = .true.
        endif

        do it=1,2
          !-------------------------------------------------------------------------
          ! See for instance Table V in
          !                      V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
          !-------------------------------------------------------------------------
          ! Double checked on 02/01/2016
          !--------------------- X Components---------------------------------------
          if(MoreDers) then
            DenIn%LapS(:,:,:,1,it) = Laplacian(DenIn%VecS(:,:,:,1,it),             &
            &                              ParityInt,-SignatureInt,TimeSimplexInt,1)
          endif
          DenIn%DerS(:,:,:,1,1,it) = &
          & DeriveX(DenIn%VecS(:,:,:,1,it),ParityInt,-Signatureint,TimeSimplexInt,1)
          DenIn%DerS(:,:,:,2,1,it) = &
          & DeriveY(DenIn%VecS(:,:,:,1,it),ParityInt,-Signatureint,TimeSimplexInt,1)
          DenIn%DerS(:,:,:,3,1,it) = &
          & DeriveZ(DenIn%VecS(:,:,:,1,it),ParityInt,-Signatureint,TimeSimplexInt,1)

          !--------------------- Y Components--------------------------------------
          if(MoreDers) then
            DenIn%LapS(:,:,:,2,it) = Laplacian(DenIn%VecS(:,:,:,2,it),             &
            &                             ParityInt,-SignatureInt,TimeSimplexInt,2)
          endif
          DenIn%DerS(:,:,:,1,2,it) = &
          &DeriveX(DenIn%VecS(:,:,:,2,it),ParityInt,-Signatureint,TimeSimplexInt,2)
          DenIn%DerS(:,:,:,2,2,it) = &
          &DeriveY(DenIn%VecS(:,:,:,2,it),ParityInt,-Signatureint,TimeSimplexInt,2)
          DenIn%DerS(:,:,:,3,2,it) = &
          &DeriveZ(DenIn%VecS(:,:,:,2,it),ParityInt,-Signatureint,TimeSimplexInt,2)
          !--------------------- Z Components--------------------------------------
          if(MoreDers) then
            DenIn%LapS(:,:,:,3,it) = Laplacian(DenIn%VecS(:,:,:,3,it),             &
            &                              ParityInt,SignatureInt,TimeSimplexInt,1)
          endif
          DenIn%DerS(:,:,:,1,3,it) = &
          & DeriveX(DenIn%VecS(:,:,:,3,it),ParityInt,Signatureint,TimeSimplexInt,1)
          DenIn%DerS(:,:,:,2,3,it) = &
          & DeriveY(DenIn%VecS(:,:,:,3,it),ParityInt,Signatureint,TimeSimplexInt,1)
          DenIn%DerS(:,:,:,3,3,it) = &
          & DeriveZ(DenIn%VecS(:,:,:,3,it),ParityInt,Signatureint,TimeSimplexInt,1)
          !------------------------------------------------------------------------
          ! Rotation of S
          DenIn%RotS(:,:,:,1,it) = DenIn%DerS(:,:,:,2,3,it)-DenIn%DerS(:,:,:,3,2,it)
          DenIn%RotS(:,:,:,2,it) = DenIn%DerS(:,:,:,3,1,it)-DenIn%DerS(:,:,:,1,3,it)
          DenIn%RotS(:,:,:,3,it) = DenIn%DerS(:,:,:,1,2,it)-DenIn%DerS(:,:,:,2,1,it)

          if(MoreDers) then
            DivS(:,:,:,it)=0.0_dp
            do i=1,3
              DivS(:,:,:,it) = DenIn%DivS(:,:,:,it) + DenIn%DerS(:,:,:,i,i,it)
            enddo
          endif
          !------------------------------------------------------------------------
          ! Second derivatives of S
          if(allocated(DenIn%RtauN2LO)) then
            DenIn%D2S(:,:,:,1,1,1,it) = &
            & DeriveX(DenIn%DerS(:,:,:,1,1,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
            DenIn%D2S(:,:,:,2,1,1,it) = &
            & DeriveY(DenIn%DerS(:,:,:,1,1,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
            DenIn%D2S(:,:,:,3,1,1,it) = &
            & DeriveZ(DenIn%DerS(:,:,:,1,1,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
            
            DenIn%D2S(:,:,:,1,1,2,it) = &      ! ParityInt,-Signatureint, TimeSimplexInt,2
            & DeriveX(DenIn%DerS(:,:,:,1,2,it), -ParityInt, SignatureInt, TimeSimplexInt,2)
            DenIn%D2S(:,:,:,2,1,2,it) = &
            & DeriveY(DenIn%DerS(:,:,:,1,2,it), -ParityInt, SignatureInt, TimeSimplexInt,2)
            DenIn%D2S(:,:,:,3,1,2,it) = &
            & DeriveZ(DenIn%DerS(:,:,:,1,2,it), -ParityInt, SignatureInt, TimeSimplexInt,2)
            
            DenIn%D2S(:,:,:,1,1,3,it) = &      ! ParityInt, Signatureint, TimeSimplexInt,1
            & DeriveX(DenIn%DerS(:,:,:,1,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
            DenIn%D2S(:,:,:,2,1,3,it) = &
            & DeriveY(DenIn%DerS(:,:,:,1,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
            DenIn%D2S(:,:,:,3,1,3,it) = &
            & DeriveZ(DenIn%DerS(:,:,:,1,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
            
            DenIn%D2S(:,:,:,2,2,1,it) = &      ! ParityInt,-Signatureint, TimeSimplexInt,1
            & DeriveY(DenIn%DerS(:,:,:,2,1,it), -ParityInt, SignatureInt, TimeSimplexInt,2)
            DenIn%D2S(:,:,:,3,2,1,it) = &
            & DeriveZ(DenIn%DerS(:,:,:,2,1,it), -ParityInt, SignatureInt, TimeSimplexInt,2)
            
            DenIn%D2S(:,:,:,2,2,2,it) = &      ! ParityInt,-Signatureint, TimeSimplexInt,2
            & DeriveY(DenIn%DerS(:,:,:,2,2,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
            DenIn%D2S(:,:,:,3,2,2,it) = &
            & DeriveZ(DenIn%DerS(:,:,:,2,2,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
          
            DenIn%D2S(:,:,:,2,2,3,it) = &      ! ParityInt, Signatureint, TimeSimplexInt,1
            & DeriveY(DenIn%DerS(:,:,:,2,3,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)
            DenIn%D2S(:,:,:,3,2,3,it) = &
            & DeriveZ(DenIn%DerS(:,:,:,2,3,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)
          
            DenIn%D2S(:,:,:,3,3,1,it) = &      ! ParityInt,-Signatureint, TimeSimplexInt,1
            & DeriveZ(DenIn%DerS(:,:,:,3,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)
            
            DenIn%D2S(:,:,:,3,3,2,it) = &      ! ParityInt,-Signatureint, TimeSimplexInt,2
            & DeriveZ(DenIn%DerS(:,:,:,3,2,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)
            
            DenIn%D2S(:,:,:,3,3,3,it) = &      ! ParityInt, Signatureint, TimeSimplexInt,1)
            & DeriveZ(DenIn%DerS(:,:,:,3,3,it), -ParityInt, SignatureInt, TimeSimplexInt,1)
          
            DenIn%D2S(:,:,:,1,2,1,it) = DenIn%D2S(:,:,:,2,1,1,it)
            DenIn%D2S(:,:,:,1,3,1,it) = DenIn%D2S(:,:,:,3,1,1,it)
            DenIn%D2S(:,:,:,2,3,1,it) = DenIn%D2S(:,:,:,3,2,1,it)
            
            DenIn%D2S(:,:,:,1,2,2,it) = DenIn%D2S(:,:,:,2,2,2,it)
            DenIn%D2S(:,:,:,1,3,2,it) = DenIn%D2S(:,:,:,3,2,2,it)
            DenIn%D2S(:,:,:,2,3,2,it) = DenIn%D2S(:,:,:,3,2,2,it)
            
            DenIn%D2S(:,:,:,1,2,3,it) = DenIn%D2S(:,:,:,2,3,3,it)
            DenIn%D2S(:,:,:,1,3,3,it) = DenIn%D2S(:,:,:,3,3,3,it)
            DenIn%D2S(:,:,:,2,3,3,it) = DenIn%D2S(:,:,:,3,2,3,it)
            
          endif
        enddo

        if(.not. MoreDers) return
         !This one is dubious: several numerical derivatives are "stacked".
            do it=1,2
              ! See for instance Table V in
              ! V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
              DenIn%GradDivS(:,:,:,1,it) = DeriveX( &
              &          DivS(:,:,:,it), -ParityInt, SignatureInt, TimeSimplexInt, 1)
              DenIn%GradDivS(:,:,:,2,it) = DeriveY( &
              &          DivS(:,:,:,it), -ParityInt, SignatureInt, TimeSimplexInt, 1)
              DenIn%GradDivS(:,:,:,3,it) = DeriveZ( &
              &          DivS(:,:,:,it), -ParityInt, SignatureInt, TimeSimplexInt, 1)
            enddo
        endif
  end subroutine ComputeDensity

  function DeriveJmuNu(JmuNu) result(NablaJ)
    !---------------------------------------------------------------------------
    ! Compute nabla J by deriving JmuNu
    !---------------------------------------------------------------------------
    use Derivatives, only : DeriveX,DeriveY,DeriveZ

    integer       :: i, it
    real(KIND=dp) :: ContractionJ(nx,ny,nz,3,2), DerContraction(nx,ny,nz,3)
    real(KIND=dp), intent(in) :: Jmunu(nx,ny,nz,3,3,2)
    real(KIND=dp)             :: NablaJ(nx,ny,nz,2)
    
    call stp('This is bugged!')

    ! Summing the vector part
    ContractionJ(:,:,:,1,:)= JMuNu(:,:,:,2,3,:) - JMuNu(:,:,:,3,2,:)
    ContractionJ(:,:,:,2,:)= JMuNu(:,:,:,3,1,:) - JMuNu(:,:,:,1,3,:)
    ContractionJ(:,:,:,3,:)= JMuNu(:,:,:,1,2,:) - JMuNu(:,:,:,2,1,:)

    !Deriving every component of the vector
    NablaJ = 0.0_dp
    do it=1,2
      do i=1,3
        DerContraction(:,:,:,i) = DeriveX( &!!Maybe the bug is the non-changing X?
        &   ContractionJ(:,:,:,i,it),-ParityInt,-SignatureInt, TimeSimplexInt,1)
      enddo
      NablaJ(:,:,:,it) = sum(DerContraction,4)
    enddo
  end function DeriveJmuNu

  function DeriveVecJ(VecJ) result(RotVecJ)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the rotation of the \vec{j}
    !---------------------------------------------------------------------------
    use Derivatives, only : DeriveX,DeriveY,DeriveZ

    integer       :: it
    real(KIND=dp) :: DerVecj(nx,ny,nz,3,3), RotVecJ(nx,ny,nz,3,2)
    real(KIND=dp), intent(in) :: VecJ(nx,ny,nz,3,2)

    do it=1,2
      !See Table V in V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
      DerVecJ(:,:,:,2,1)  = &
      & DeriveY(Vecj(:,:,:,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)

      DerVecJ(:,:,:,3,1)  = &
      & DeriveZ(Vecj(:,:,:,1,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)

      DerVecJ(:,:,:,1,2)  = &
      & DeriveX(Vecj(:,:,:,2,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)

      DerVecJ(:,:,:,3,2)  = &
      & DeriveZ(Vecj(:,:,:,2,it), -ParityInt,-SignatureInt, TimeSimplexInt,2)

      DerVecJ(:,:,:,1,3)  = &
      & DeriveX(Vecj(:,:,:,3,it), -ParityInt,-SignatureInt, TimeSimplexInt,1)

      DerVecJ(:,:,:,2,3)  = &
      & DeriveY(Vecj(:,:,:,3,it), -ParityInt, SignatureInt, TimeSimplexInt,2)

      RotVecJ(:,:,:,1,it) =  DerVecJ(:,:,:,2,3) - DerVecJ(:,:,:,3,2)
      RotVecJ(:,:,:,2,it) =  DerVecJ(:,:,:,3,1) - DerVecJ(:,:,:,1,3)
      RotVecJ(:,:,:,3,it) =  DerVecJ(:,:,:,1,2) - DerVecJ(:,:,:,2,1)
    enddo
  end function DeriveVecJ

  subroutine DeriveVecs(DenIn)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the laplacian, the rotator, the divergence and
    ! the gradient of the divergence of \vec{s}.
    !---------------------------------------------------------------------------
    use Derivatives, only : Laplacian, DeriveX,DeriveY,DeriveZ
    use Force
    integer :: i,it
    type(DensityVector), intent(inout) :: DenIn
    real(KIND=dp)                      :: DivS(nx,ny,nz,2)
    logical                            :: MoreDers

    MoreDers=.true.
    if(B18.eq.0.0_dp.and.B19.eq.0.0_dp.and.B20.eq.0.0_dp .and. B21.eq.0.0_dp) then
      MoreDers = .false.
    endif

    do it=1,2
      !-------------------------------------------------------------------------
      ! See for instance Table V in
      !                      V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
      !-------------------------------------------------------------------------
      ! Double checked on 02/01/2016
      !--------------------- X Components---------------------------------------
      if(MoreDers) then
        DenIn%LapS(:,:,:,1,it) = Laplacian(DenIn%VecS(:,:,:,1,it),             &
        &                              ParityInt,-SignatureInt,TimeSimplexInt,1)
      endif

      DenIn%DerS(:,:,:,1,1,it) = &
      & DeriveX(DenIn%VecS(:,:,:,1,it),ParityInt,-Signatureint,TimeSimplexInt,1)
      DenIn%DerS(:,:,:,2,1,it) = &
      & DeriveY(DenIn%VecS(:,:,:,1,it),ParityInt,-Signatureint,TimeSimplexInt,1)
      DenIn%DerS(:,:,:,3,1,it) = &
      & DeriveZ(DenIn%VecS(:,:,:,1,it),ParityInt,-Signatureint,TimeSimplexInt,1)

      !--------------------- Y Components--------------------------------------
      if(MoreDers) then
        DenIn%LapS(:,:,:,2,it) = Laplacian(DenIn%VecS(:,:,:,2,it),             &
        &                             ParityInt,-SignatureInt,TimeSimplexInt,2)
      endif
      DenIn%DerS(:,:,:,1,2,it) = &
      &DeriveX(DenIn%VecS(:,:,:,2,it),ParityInt,-Signatureint,TimeSimplexInt,2)
      DenIn%DerS(:,:,:,2,2,it) = &
      &DeriveY(DenIn%VecS(:,:,:,2,it),ParityInt,-Signatureint,TimeSimplexInt,2)
      DenIn%DerS(:,:,:,3,2,it) = &
      &DeriveZ(DenIn%VecS(:,:,:,2,it),ParityInt,-Signatureint,TimeSimplexInt,2)

      !--------------------- Z Components--------------------------------------
      if(MoreDers) then
        DenIn%LapS(:,:,:,3,it) = Laplacian(DenIn%VecS(:,:,:,3,it),             &
        &                              ParityInt,SignatureInt,TimeSimplexInt,1)
      endif
      DenIn%DerS(:,:,:,1,3,it) = &
      & DeriveX(DenIn%VecS(:,:,:,3,it),ParityInt,Signatureint,TimeSimplexInt,1)
      DenIn%DerS(:,:,:,2,3,it) = &
      & DeriveY(DenIn%VecS(:,:,:,3,it),ParityInt,Signatureint,TimeSimplexInt,1)
      DenIn%DerS(:,:,:,3,3,it) = &
      & DeriveZ(DenIn%VecS(:,:,:,3,it),ParityInt,Signatureint,TimeSimplexInt,1)
      !------------------------------------------------------------------------

      DenIn%RotS(:,:,:,1,it) = DenIn%DerS(:,:,:,2,3,it)-DenIn%DerS(:,:,:,3,2,it)
      DenIn%RotS(:,:,:,2,it) = DenIn%DerS(:,:,:,3,1,it)-DenIn%DerS(:,:,:,1,3,it)
      DenIn%RotS(:,:,:,3,it) = DenIn%DerS(:,:,:,1,2,it)-DenIn%DerS(:,:,:,2,1,it)

      if(MoreDers) then
        DivS(:,:,:,it)=0.0_dp
        do i=1,3
          DivS(:,:,:,it) = DenIn%DivS(:,:,:,it) + DenIn%DerS(:,:,:,i,i,it)
        enddo
      endif
    enddo

    if(.not. MoreDers) return
     !This one is dubious: several numerical derivatives are "stacked".
    do it=1,2
      ! See for instance Table V in
      ! V. Hellemans et al., Phys. Rev. C 85 (2012), 014326
      DenIn%GradDivS(:,:,:,1,it) = DeriveX( &
      &          DivS(:,:,:,it), -ParityInt, SignatureInt, TimeSimplexInt, 1)
      DenIn%GradDivS(:,:,:,2,it) = DeriveY( &
      &          DivS(:,:,:,it), -ParityInt, SignatureInt, TimeSimplexInt, 1)
      DenIn%GradDivS(:,:,:,3,it) = DeriveZ( &
      &          DivS(:,:,:,it), -ParityInt, SignatureInt, TimeSimplexInt, 1)
    enddo
  end subroutine DeriveVecs

  subroutine WriteDensity(Den, unit)
  !-----------------------------------------------------------------------------
  ! Writes the density vector den to a file opened on unit.
  !-----------------------------------------------------------------------------
  ! Note that we don't have to store all the different densities on file as a
  ! lot of them are derivatives of others. However, for simplicity, we do that
  ! anyway: the added storage cost is on the order of a couple of Spwfs at most.
  !-----------------------------------------------------------------------------
    type(DensityVector), intent(in) :: Den
    integer, intent(in)             :: unit
    integer                         :: i
    !---------------------------------------------------------------------------
    ! Things that are certainly allocated
    write(unit) Den%Rho, Den%DerRho, Den%LapRho
    write(unit) Den%Tau
    write(unit) Den%NablaJ

    if(allocated(Den%Jmunu)) then
      write(unit) Den%JMuNu
    else
      !Write an empty line
      write(unit)
    endif
    !---------------------------------------------------------------------------
    ! Things that are not necessarily allocated
    if(.not.TRC) then
      write(unit) Den%vecj, Den%RotVecj
      write(unit) Den%Vecs, Den%Rots, Den%DerS

      if(allocated(Den%LapS)) then
        write(unit) Den%LapS, Den%DivS, Den%GradDivS
      else
        write(unit)
      endif

      if(allocated(Den%VecT)) then
        write(unit) Den%VecT
      else
        write(unit)
      endif
      if(allocated(Den%VecF)) then
        write(unit) Den%VecF
      else
        write(unit)
      endif
    else
      do i=1,5
        write(unit)
      enddo
    endif

  end subroutine WriteDensity

  subroutine ReadDensity(Den, unit, filenx,fileny,filenz, flag)
  !-----------------------------------------------------------------------------
  ! Assigns a density vector den values read from file unit.
  ! Everything gets newly allocated.
  !
  ! Flag signals if problems have been encountered.
  !
  !-----------------------------------------------------------------------------
    type(DensityVector), intent(inout) :: Den
    integer, intent(in)                :: unit
    integer                            :: io, i, SizeOfDen
    integer, intent(in)                :: filenx, fileny, filenz
    integer, intent(inout)             :: flag

    Den = NewDensityVector(filenx,fileny,filenz)

    Flag = 0

    read(unit,iostat=flag) Den%Rho, Den%DerRho, Den%LapRho
    read(unit,iostat=flag) Den%Tau
    read(unit,iostat=flag) Den%NablaJ

    if(allocated(Den%Jmunu)) then
      io = 0
      SizeOfDen = size(Den%JMuNu)
      ! Read the density as a vector. My clever scheme of trying to read
      ! the densities fails if this is read as a matrix.
      read(unit,iostat=io) Den%JMuNu(1:SizeOfDen,1,1,1,1,1)
      if(io.ne.0) then
        Den%JMuNu = 0.0_dp
      endif
    else
      read(unit)
    endif

    !Only attempt to read these densities when Time-reversal is not conserved.
    ! The densities however are not guaranteed to be on file, which is why we
    ! utilize some iostat tricks.
    ! Note that the IOSTAT end-of-record code on FORTRAN is -2.
    if( TRC) then
      io = 0
      read(unit, iostat=io) Den%vecj,Den%RotVecj
      if(io.ne.0) then
        !End-of-record problem, no densities on file
        Den%Vecj    = 0.0_dp
        Den%Rotvecj = 0.0_dp
      endif

      io = 0
      read(unit, iostat=io) Den%vecS,Den%RotS, Den%DerS
      if(io.ne.0) then
        !End-of-record problem, no densities on file
        Den%Vecs    = 0.0_dp
        Den%RotS    = 0.0_dp
        Den%DerS    = 0.0_dp
      endif

      if(allocated(Den%LapS)) then
        io = 0
        read(unit, iostat=io) Den%LapS,Den%DivS, Den%GradDivS
        if(io.ne.0) then
          !End-of-record problem, no densities on file
          Den%Laps    = 0.0_dp
          Den%DivS    = 0.0_dp
          Den%GradDivS= 0.0_dp
        endif
      else
        read(unit)
      endif

      if(allocated(Den%VecT)) then
        io = 0
        read(unit, iostat=io) Den%VecT
        if(io.ne.0) then
          !End-of-record problem, no densities on file
          Den%VecT    = 0.0_dp
        endif
      else
        read(unit)
      endif

      if(allocated(Den%VecF)) then
        io = 0
        read(unit, iostat=io) Den%VecF
        if(io.ne.0) then
          !End-of-record problem, no densities on file
          Den%VecF    = 0.0_dp
        endif
      else
        read(unit)
      endif
    else
      do i=1,5
        read(unit)
      enddo
    endif
  end subroutine ReadDensity

  subroutine FillOct( NewDen, OldDen , Octant, S)
  !-----------------------------------------------------------------------------
  ! Subroutine that copies the content of OldDen to NewDen into the correct
  ! quadrant with a sign S. This is useful for breaking the symmetries of
  ! densities read on file.
  !
  ! The octant that needs to be filled are numbered as follows:
  !
  ! 1) x > 0, y > 0, z > 0
  ! 2) x < 0, y > 0, z > 0
  ! 3) x > 0, y < 0, z > 0
  ! 4) x > 0, y > 0, z < 0
  ! 5) x < 0, y < 0, z > 0
  ! 6) x < 0, y > 0, z < 0
  ! 7) x > 0, y < 0, z < 0
  ! 8) x < 0, y < 0, z < 0
  !
  ! Note that the final dimension in the arrays is ofcourse an isospin index.
  !-----------------------------------------------------------------------------
  ! To anyone debugging this, I apologise in advance for the abstraction of this
  ! routine. Some hints:
  !
  ! To fill in an octant, we first select the indices of the 'filling' loop and
  ! whether or not we have to 'invert' the values from input.
  ! 1) 'Inverting' means here that we have to traverse the input in opposite
  !    direction: i.e. for parity-breaking transformations we need to
  !
  !     NewDen(1,1,1) = (+/-) OldDen(1,1,nz  )
  !     NewDen(1,1,2) = (+/-) OldDen(1,1,nz-1)
  !     ......
  !
  !    which explains the integers sx,sy,sz. They control 'inversion' for each
  !    direction and are ofcourse dependent on the octant.
  ! 2) Indices of the filling loop:
  !    these are controlled by three integers: the offsets and the start/end.
  !    The start/end indices relate to the new density and thus are the indices
  !    to fill for the octant. The offset is then the difference between new
  !    and old indices.
  !-----------------------------------------------------------------------------
    integer, intent(in) :: Octant, S
    real(KIND=dp), intent(inout)   :: NewDen(nx,ny,nz,2)
    real(KIND=dp), intent(in)      :: OldDen(:,:,:,:)
    integer                        :: sizex,sizey,sizez, i,j,k
    integer                        :: offsetx,offsety,offsetz, sx, sy, sz
    integer                        :: startx,endx,starty,endy,startz,endz,x,y,z

    sizex = size(OldDen,1) ; sizey = size(OldDen,2) ; sizez = size(OldDen,3)

    !Checking which octants are there
    if(sizex.eq.nx) offsetx = 0
    if(sizex.ne.nx) offsetx = nx/2
    if(sizey.eq.ny) offsety = 0
    if(sizey.ne.ny) offsety = ny/2
    if(sizez.eq.nz) offsetz = 0
    if(sizez.ne.nz) offsetz = nz/2

    select case(Octant)
    case(1)
      ! x > 0, y > 0, z > 0
      startx = offsetx+1 ; endx = nx    ; sx =  1
      starty = offsety+1 ; endy = ny    ; sy =  1
      startz = offsetz+1 ; endz = nz    ; sz =  1
    case(2)
      ! x < 0, y > 0, z > 0
      startx = 1         ; endx = sizex ; sx = -1
      starty = offsety+1 ; endy = ny    ; sy =  1
      startz = offsetz+1 ; endz = nz    ; sz =  1
    case(3)
      ! x > 0, y < 0, z > 0
      startx = offsetx+1 ; endx = nx    ; sx =  1
      starty = 1         ; endy = sizey ; sy = -1
      startz = offsetz+1 ; endz = nz    ; sz =  1
    case(4)
      ! x > 0, y > 0, z < 0
      startx = offsetx+1 ; endx = nx    ; sx =  1
      starty = offsety+1 ; endy = ny    ; sy =  1
      startz = 1         ; endz = sizez ; sz = -1
    case(5)
      ! x < 0, y < 0, z > 0
      startx = 1         ; endx = sizex ; sx = -1
      starty = 1         ; endy = sizey ; sy = -1
      startz = offsetz+1 ; endz = nz    ; sz =  1
    case(6)
      ! x < 0, y > 0, z < 0
      startx = 1         ; endx = sizex ; sx = -1
      starty = offsety+1 ; endy = ny    ; sy =  1
      startz = 1         ; endz = sizez ; sz = -1
    case(7)
      ! x > 0, y < 0, z < 0
      startx = offsetx+1 ; endx = nx    ; sx =  1
      starty = 1         ; endy = sizey ; sy = -1
      startz = 1         ; endz = sizez ; sz = -1
    case(8)
      ! x < , y < 0 , z < 0
      startx = 1         ; endx = sizex ; sx = -1
      starty = 1         ; endy = sizey ; sy = -1
      startz = 1         ; endz = sizez ; sz = -1
    case DEFAULT
      call stp('FillInOctant does not recognise the Octant.')
    end select
    !---------------------------------------------------------------------------
    !Do the filling in
    do k=startz,endz
      z = sz * ( k - offsetz) + (1 - sz)/2
      do j=starty,endy
        y = sy * ( j - offsety) + (1 - sy)/2
        do i=startx,endx
          x = sx * ( i - offsetx) + (1 - sx)/2
          if(y.eq.0) then
              print *,j,sy,y,starty,endy,offsety
          endif
          NewDen(i,j,k,:) = S * OldDen(x,y,z,:)
        enddo
      enddo
    enddo
  end subroutine FillOct

end module Densities
