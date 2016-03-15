module BCS
!-------------------------------------------------------------------------------
! Module that provides the solution of the BCS equations to MOCCa.
!
!-------------------------------------------------------------------------------

  use CompilationInfo
  use Geninfo
  use SpwfStorage
  use PairingInteraction

  implicit none

contains

    subroutine BCSGaps(Delta,DeltaLN,PairingField,PairingFieldLN,Gaps,ConstantGap)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the Delta_{i ibar} for the BCS model.
    ! The formula (see Ring & Shuck, page 231 and 232, eq. 6.50 and eq.6.54)
    ! is the following:
    !
    !   Delta_{i ibar} =
    !         - 0.5 * sum_j \bar{v}_{i,ibar,j,jbar} * Delta_{j,jbar}
    !         / sqrt( epsilon_j**2 + Delta_{j,jbar}**2)
    !
    ! where we replace the Delta_{j,jbar} by the values of the gaps of
    ! the previous mean-field iteration.
    !---------------------------------------------------------------------------
    ! Currently just using the straight BCS formula
    !
    ! Delta_{i, ibar} = \int dr <i,ibar | Delta(r) | r >
    !---------------------------------------------------------------------------
    ! Of course, when using constant gap pairing, the gaps are not calculated
    ! but set to the constant values in Gaps (with added cutoffs)
    !
    !---------------------------------------------------------------------------
    complex(KIND=dp), allocatable,intent(inout) :: Delta(:,:,:,:)
    complex(KIND=dp), allocatable,intent(inout) :: DeltaLN(:,:,:,:)
    complex(KIND=dp), allocatable,intent(in)    :: Pairingfield(:,:,:,:)
    complex(KIND=dp), allocatable,intent(in)    :: PairingfieldLN(:,:,:,:)
    real(KIND=dp), intent(in)                   :: Gaps(2)
    logical, intent(in)                         :: ConstantGap
    integer                                     :: i, it
    real(KIND=dp)                               :: Cutoff

    Delta = 0.0_dp
    if(ConstantGap) then
        !Assign constant Deltas
        do i=1,nwt
            Cutoff = PCutoffs(i)
            it = (HFBasis(i)%GetIsospin()+3)/2
            Delta(i,1,1,1)=Cutoff**2*Gaps(it)
        enddo
    else
        !Calculate the Deltas
        do i=1,nwt
            Cutoff = PCutoffs(i)
            it = (HFBasis(i)%GetIsospin()+3)/2
            Delta(i,1,1,1) = Cutoff**2*dv*                                     &
            &            sum(HFBasis(i)%GetDensity()*PairingField(:,:,:,it))
        enddo
    endif
    end subroutine BCSGaps

    subroutine BCSPairingField(Field, FieldLN, Delta)
    !---------------------------------------------------------------------------
    ! Computes the pairing field due to a BCS model pairing.
    !  Delta(x,y,z) =
    !  Sum_{i,ibar,s}  0.5*fi*fibar*u_i*v_i * s * <r,s;r,-s|V|ij>
    !---------------------------------------------------------------------------
    ! FieldLN and DeltaLN are dummies for now.
    !
    !---------------------------------------------------------------------------
      use Spinors

      complex(KIND=dp), allocatable, intent(inout) :: Field(:,:,:,:)
      complex(KIND=dp), allocatable, intent(inout) :: FieldLN(:,:,:,:)
      complex(KIND=dp), allocatable,intent(inout)  :: Delta(:,:,:,:)
      integer      :: i, it,j
      real(KIND=dp):: Cutoff, factor
      type(Spinor) :: Psi
      complex(KIND=dp)                             :: Temp(nx,ny,nz)

      Field  = 0.0_dp
      do i=1,nwt
          it = (HFBasis(i)%GetIsospin() +3)/2
          Cutoff = PCutoffs(i)
          Psi  = HFBasis(i)%GetValue()
          Temp = GetPairDensity(Psi,Psi,.true.)
          factor = 0.5_dp*Cutoff**2*Delta(i,1,1,1)/(HFBasis(i)%GetEqp())
          do j=1,nx*ny*nz
            Field(j,1,1,it) =  Field(j,1,1,it) + factor*Temp(j,1,1)
          enddo
      enddo
      do it=1,2
        do j=1,nx*ny*nz
          Field(j,1,1,it) = DensityFactor(j,1,1,it) * Field(j,1,1,it)
        enddo
      enddo
    end subroutine BCSPairingField

   subroutine BCSFindFermiEnergy (Fermi,LNLambda,Delta,DeltaLN,Lipkin,         &
     &                                      DN2, ConstrainDispersion, FermiPrec)
   !----------------------------------------------------------------------------
   ! Function that finds the a) correct Fermi energy and b) the correct second
   ! Kamlah moment (when desired) for the BCS system.
   ! The number of particles of species it is also calculated and returned.
   ! Do this for one isospin index only, for consistency.
   !----------------------------------------------------------------------------
   ! At the moment, this takes a lot of dummy arguments....
   !----------------------------------------------------------------------------

    real(KIND=dp),intent(inout) :: Fermi(2), LNLambda(2)
    real(KIND=dp), intent(in)   :: FermiPrec, DN2(2)
    complex(KIND=dp), allocatable,intent(in)  :: Delta(:,:,:,:),DeltaLN(:,:,:,:)
    logical, intent(in)         :: Lipkin, ConstrainDispersion

    real(KIND=dp)               :: LambdaSums(2,2), eqp, nom, Particles(2)
    integer                     :: i
    integer                     :: it

    Particles(1) = Neutrons; Particles(2) = Protons

    call BCSQPEnergies(Fermi,Delta,LNLambda)

    LambdaSums=0.0_dp
    !Sum the quasiparticle energies to get the correct lambda.
    do i=1,nwt
     it = (HFBasis(i)%GetIsospin() + 3)/2
     eqp = HFBasis(i)%GetEqp()

     !nom is epsilon' from W.Ryssens et al.
     nom = HFBasis(i)%GetEnergy() + 4*LNLambda(it)*(HFBasis(i)%GetOcc()-0.5_dp)

     Lambdasums(1,it)= lambdasums(1,it) +(1.0_dp - nom/eqp)
     lambdasums(2,it)= lambdasums(2,it) + 1.0_dp/eqp

     ! The number of particles is found by summing occupation numbers.
     ! v^2  = 0.5 * ( 1 - [epsilon' - lambda]/Eqp )
     ! and a factor of two is present because of timereversal.
     !N(it) = N(it) + (1.0_dp - (nom - Fermi(it))/eqp)
    enddo

    do it=1,2
        Fermi(it) =  (Particles(it) - lambdasums(1,it))/lambdasums(2,it)
    enddo

   end subroutine BCSFindFermiEnergy

   subroutine BCSQPEnergies(Fermi, Delta, LNLambda)
   !----------------------------------------------------------------------------
   ! This subroutine calculates the BCS quasiparticle enrgies as a function of
   ! the Fermi energies and the pairing gaps.
   !
   ! eqp_k = sqrt[(epsilon_k - \lambda)**2 + Delta_{k, kbar}**2]
   !
   ! which is formula (6.72) on page 235 in Ring & Shuck.
   ! This gets modified when the second Kamlah moment is nonzero.
   ! The expression becomes:
   !
   ! eqp_k = sqrt[(epsilon_k - \lambda + 4*Lambda_2*(v^2 - 0.5))**2
   !                         + Delta_{k, kbar}**2]
   !
   ! as in W.Ryssens et al., CPC 187 (2015) 175–194
   !
   !----------------------------------------------------------------------------
   ! Note that I'm completely unsure what this becomes when Time Simplex is
   ! not conserved....
   !----------------------------------------------------------------------------
    real(KIND=dp), intent(in)                 :: Fermi(2), LNLambda(2)
    complex(KIND=dp), intent(in), allocatable :: Delta(:,:,:,:)
    integer                   :: i
    integer                   :: it
    real(KIND=dp)             :: epsilon, lambda
    real(KIND=dp)             :: eqp(nwt), LNTerm

    eqp = 0.0_dp
    do i=1,nwt
        it = (HFBasis(i)%GetIsospin()+3)/2

        epsilon = HFBasis(i)%GetEnergy()
        lambda  = Fermi(it)
        LNTerm  = 4 * LNLambda(it) * (HFBasis(i)%GetOcc()**2 - 0.5_dp)
        eqp(i)  = sqrt((epsilon - lambda + LNTerm)**2 + abs(Delta(i,1,1,1))**2)
        call HFBasis(i)%SetEqp(eqp(i))
    enddo
   end subroutine BCSQPEnergies

   subroutine BCSOccupations(Fermi, Delta,LNLambda, PairingDisp)
   !----------------------------------------------------------------------------
   ! Find the occupation numbers of the HFBasis from the quasiparticle energies
   ! and the Fermi energies.
   !
   ! v^2_k = 0.5 * 1 - (Epsilon - Lambda)/(E_{qp})
   !
   ! or formula (6.51) on page 231 in Ring & Schuck.
   !----------------------------------------------------------------------------
   ! This routine also calculates the dispersion in the particle number:
   ! \Delta N^2 = Tr( rho (1-rho)) = sum v^2_i ( 1 - v^2_i)
   !
   !----------------------------------------------------------------------------
    real(KIND=dp), intent(in)   :: Fermi(2),LNLambda(2)
    real(KIND=dp), intent(inout):: PairingDisp(2)
    integer                     :: i,it
    real(KIND=dp)               :: Occupation, eqp
    complex(KIND=dp), intent(in), allocatable :: Delta(:,:,:,:)

    call BCSQPEnergies(Fermi,Delta,LNLambda)

    PairingDisp=0.0_dp
    do i=1,nwt
        it  = (HFBasis(i)%GetIsospin()+3)/2
        eqp = HFBasis(i)%GetEqp()
        !There is already an intrinsic factor two here due to Time-reversal
        Occupation = (1 - (HFBasis(i)%GetEnergy() - Fermi(it))/eqp)
        if(Occupation.lt.0.0_dp) Occupation=0.0_dp
        call HFBasis(i)%SetOcc(Occupation)

        !We remove a factor two from the occupations, but need a global factor
        !two since we only sum over the positive signature states.
        ! DN^2 = 2*Tr(rho(1-rho))
        PairingDisp(it)=PairingDisp(it)+2*Occupation*(1.0_dp-Occupation/2.0_dp)
    enddo

   end subroutine BCSOccupations

   function BCSEnergy(Delta) result(Energy)
    !---------------------------------------------------------------------------
    ! This function computes the pairing energy for the BCS state.
    !
    !     Energy =- 1/2 * sum_i [ (f_i^2 Delta_i^2)/eqp(i)]
    !---------------------------------------------------------------------------
    complex(KIND=dp), allocatable,intent(in) :: Delta(:,:,:,:)
    integer                                  :: i, it
    real(KIND=dp)                            :: Energy(2)

    Energy = 0.0_dp

    !Testing original
    do i=1,nwt
      it = (HFBasis(i)%GetIsospin() + 3)/2
      Energy(it) = Energy(it)+ abs(Delta(i,1,1,1))**2 / HFBasis(i)%GetEqp()
    enddo
    Energy = - 0.5 * Energy
    return
   end function BCSEnergy

end module BCS
