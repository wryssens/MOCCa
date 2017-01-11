module GradientHFB
  !-----------------------------------------------------------------------------
  ! Module that offers a new fermi-solver for HFB calculations.
  ! See the following papers for more info:
  ! *) J.L. Egido et al. , Nucl. Phys. A 594 (1995) 70 - 86
  ! *) L.M. Robledo et al., Phys. Rev. C 84, 014312 (2011)
  !
  ! The main advantage is that this solver controls the possible qp-excitations
  ! in a more transparent way: it doesn't admit them at all.
  !
  ! Another advantage is that this method is truely variational.
  !
  !-----------------------------------------------------------------------------
  ! Technical notes:
  !
  ! *) This module liberally uses the matmul intrinsic. For typical matrix sizes
  !    that are relevant to this problem ( 10 ~ 100) this is the fastest way
  !    I've tested to multiply matrices. It is possible that LAPACK routines
  !    (GEMM routines) might do better for larger matrices, but I highly doubt
  !    this.
  !-----------------------------------------------------------------------------
  ! TODO:
  ! * Optimize the signature breaking routines.
  ! *
  !-----------------------------------------------------------------------------
  use HFB

  implicit none

  !-----------------------------------------------------------------------------
  ! Local storage for the U and V
  real(KIND=dp),allocatable  :: GradU(:,:,:,:), GradV(:,:,:,:)

  !-----------------------------------------------------------------------------
  ! Store single-particle hamiltonian and delta in block-form
  real(KIND=dp), allocatable :: hblock(:,:,:,:)   ,dblock(:,:,:,:)
  !-----------------------------------------------------------------------------
  ! Number of positive signature states in every isospin/parity block
  ! The number of negative states then is blocksizes - sigblocks
  integer :: SigBlocks(2,2) = 0

  !-----------------------------------------------------------------------------
  ! Matrices that I don't want to reallocate everytime.
  real(KIND=dp),allocatable  :: Gradient(:,:,:,:) ,OldGrad(:,:,:,:)
  real(KIND=dp),allocatable  :: aN20(:,:,:,:)
  real(KIND=dp),allocatable  :: oldgradU(:,:,:,:) ,oldgradV(:,:,:,:)
  real(KIND=dp),allocatable  :: Direction(:,:,:,:), OldDir(:,:,:,:)

  real(KIND=dp) :: over

  abstract interface
    function H20Interface(Ulim,Vlim,hlim,Dlim,S) result(H20)
      import :: dp
      integer, intent(in)        :: S
      real(KIND=dp), intent(in)  :: Ulim(:,:), Vlim(:,:), dlim(:,:), hlim(:,:)
      real(KIND=dp)              :: H20(size(Ulim,1),size(Ulim,1))
    end function
  end interface
  abstract interface
    function ConstructRhoInterface(Vlim,S) result(Rho)
      import :: dp
      real(KIND=dp), intent(in) :: Vlim(:,:)
      real(KIND=dp)             :: Rho(size(Vlim,1),size(Vlim,1))
      integer, intent(in)       :: S
    end function
  end interface
  abstract interface
    function ConstructKappaInterface(Ulim,Vlim,S) result(Kappa)
      import :: dp
      real(KIND=dp), intent(in) :: Ulim(:,:), Vlim(:,:)
      real(KIND=dp),allocatable :: Kappa(:,:)
      integer, intent(in)       :: S
    end function
  end interface
  abstract interface
    subroutine GradUpdateInterface(step, U1,V1,Grad,U2,V2, S)
      import :: dp
      real(KIND=dp), intent(in) :: U1(:,:), V1(:,:), step
      real(KIND=dp), intent(in) :: Grad(:,:)
      real(KIND=dp),intent(out) :: U2(:,:), V2(:,:)
      integer, intent(in)       :: S
    end subroutine
  end interface

  abstract interface
    subroutine orthointerface(U,V,S)
        import :: dp
        real(KIND=dp), intent(inout) :: U(:,:), V(:,:)
        integer, intent(in)          :: S
    end subroutine
  end interface
  procedure(H20Interface), pointer            :: H20
  procedure(constructrhoInterface), pointer   :: constructrho
  procedure(constructkappaInterface), pointer :: constructkappa
  procedure(gradupdateInterface), pointer     :: gradupdate
  procedure(orthoInterface), pointer          :: ortho

contains

  subroutine GetUandV
    !---------------------------------------------------------------------------
    ! Get U and V values from the HFB module.
    ! Also very important: counts the number of of (U,V) vectors in every block.
    !
    !
    !---------------------------------------------------------------------------

    integer :: it,P,i,j,N, ind(2,2), Rzindex,k, S
    logical :: incolumns

    !----------------------------------------
    ! Initialize the matrices.
    !----------------------------------------
    N = maxval(blocksizes)

    allocate(GradU(N,N,Pindex,Iindex)) ; GradU = 0.0_dp
    allocate(GradV(N,N,Pindex,Iindex)) ; GradV = 0.0_dp

    if(all(HFBColumns.eq.0)) then
      ! Safeguard for when we read legacyinput
      if(SC) then
          do it=1,Iindex
              do P=1,Pindex
                  S = blocksizes(P,it)
                  do i=1,S/2
                      HFBColumns(i,P,it) = 2*S - i + 1
                  enddo
                  do i=S/2+1,S
                      HFBColumns(i,P,it) = i
                  enddo
              enddo
          enddo
      else
        call stp('No legacy for signature broken.')
      endif
    endif
     !print *, HFBColumns(1:blocksizes(1,2),1,2)
     !print *, HFBColumns(1:blocksizes(2,2),2,2)

    if(allocated(qpexcitations)) then
      call BlockQuasiParticles
    endif

    do it=1,Iindex
        do P=1,Pindex
            N = blocksizes(P,it)
            ind(P,it) = 1
            do j=1,2*N
                !----------------------------------------------
                ! This loop makes sure that we go through the
                ! big matrix in order. Otherwise the signature
                ! block might end up somewhere we don't expect
                ! them.
                !----------------------------------------------
                incolumns = .false.
                do k=1,N
                    if (j .eq. HFBColumns(k,P,it)) then
                        incolumns=.true.
                        exit
                    endif
                enddo
                if(incolumns) then
                  GradV(1:N,ind(P,it),P,it) = DBLE(V(1:N,j,P,it))
                  GradU(1:N,ind(P,it),P,it) = DBLE(U(1:N,j,P,it))

                  !if(any(abs(GradU(1:N/2,ind(P,it),P,it)).gt.1.0d-12) .or. .not. SC) then
                  !  sigblocks(P,it) = sigblocks(P,it) +1
                  !endif
                  ind(P,it) = ind(P,it)+ 1
                endif
            enddo
        enddo
    enddo
    sigblocks = 0
    do i=1,nwt
        P  =( HFBasis(i)%GetParity()    + 3)/2
        it =( HFBasis(i)%GetIsospin()   + 3)/2
        if( ( HFBasis(i)%GetSignature() .eq. -1) .or. (.not. SC)) then
                sigblocks(P,it) = sigblocks(P,it) +1 
        endif
    enddo
    
    !print *, 'Negative parity'
    !print *, HFBColumns(1:blocksizes(1,2),1,2)

    !do i=1,blocksizes(1,2)
    ! print *, DBLE(GradU(i,1:blocksizes(1,2),1,2))
    !enddo
    !print *

    !do i=1,blocksizes(1,2)
     !print *, DBLE(GradV(i,1:blocksizes(1,2),1,2))
    !enddo
    !print *
    
    !print *, 'Positive parity'
    !print *, HFBColumns(1:blocksizes(2,2),2,2)

    !do i=1,blocksizes(2,2)
    ! print *, DBLE(GradU(i,1:blocksizes(2,2),2,2))
    !enddo
    !print *

    !do i=1,blocksizes(2,2)
    ! print *, DBLE(GradV(i,1:blocksizes(2,2),2,2))
    !enddo
    !print *
    
    !print *
    !i = 1
    !j = 2
    !call checkUandV(DBLE(GRADU(1:blocksizes(j,i), 1:blocksizes(j,i),j,i)), &
    !&               DBLE(GradV(1:blocksizes(j,i), 1:blocksizes(j,i),j,i)))

    !stop
    
  end subroutine GetUandV

  subroutine OutHFBModule
      !-------------------------------------------------------------------------
      ! Put our solution into the HFB module correctly.
      !-------------------------------------------------------------------------

      integer :: i,j,P,it,N, S

      do it=1,Iindex
          do P=1,Pindex
              N = blocksizes(P,it)
              S = sigblocks(P,it)
              do i=1,N
                  do j=1,S
                      HFBColumns(j,P,it) = j
                      U(i,j,P,it) = GradU(i,j,P,it)
                      V(i,j,P,it) = GradV(i,j,P,it)

                      ! The conjugate states
                      U(i,2*N-j+1,P,it) = GradV(i,j,P,it)
                      V(i,2*N-j+1,P,it) = GradU(i,j,P,it)
                  enddo
                  do j=1,N-S
                      HFBColumns(j+S,P,it) = j + N
                      U(i,j+N,P,it) = GradU(i,j+S,P,it)
                      V(i,j+N,P,it) = GradV(i,j+S,P,it)

                      ! The conjugate states
                      U(i,N-j+1,P,it) = GradV(i,j+S,P,it)
                      V(i,N-j+1,P,it) = GradU(i,j+S,P,it)
                  enddo
              enddo
              RhoHFB(1:N,1:N,P,it)   = ConstructRho(GradV(1:N,1:N,P,it),       &
              &                                     sigblocks(P,it))

              KappaHFB(1:N,1:N,P,it) = ConstructKappa(GradU(1:N,1:N,P,it),     &
              &                                       GradV(1:N,1:N,P,it),     &
              &                                       sigblocks(P,it))
          enddo
      enddo
      

!    print *
!    do j =1, blocksizes(1,2)
!      print *, DBLE(U(j,1:2*blocksizes(1,2),1,2))
!    enddo
!    print *
!    do j =1, blocksizes(1,2)
!      print *, DBLE(V(j,1:2*blocksizes(1,2),1,2))
!    enddo
!    print *
!      ! stop
  end subroutine OutHFBModule

  subroutine BlockHFBHamil(Delta, Fermi, L2)
      !-------------------------------------------------------------------------
      ! Construct the HFBHamiltonian and save it into practical blocks for
      ! our solver.
      !-------------------------------------------------------------------------

      complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
      real(KIND=dp), intent(in)                 :: Fermi(2), L2(2)
      real(KIND=dp)                             :: z(2) = 0.0_dp
      integer                                   :: N,i,j,P,it, Rzindex, Rz

      ! Construct
      call constructHFBHamiltonian(Fermi,Delta,L2,HFBGauge)

      N = maxval(blocksizes)

      if(.not.allocated(hblock)) then
          allocate(hblock(N,N,Pindex,Iindex))
          allocate(dblock(N,N,Pindex,Iindex))
      endif
      ! Get the single-particle hamiltonian and pairing matrix Delta
      ! (note that this already includes any LN contribution)
      do it=1,Iindex
          do P=1,Pindex
              N = blocksizes(P,it)
              hblock(1:N,1:N,P,it) = HFBHamil(1:N,1:N     ,P,it)
              dblock(1:N,1:N,P,it) = HFBHamil(1:N,N+1:2*N ,P,it)
          enddo
      enddo

  end subroutine BlockHFBHamil

  subroutine HFBFermiGradient(Fermi,L2,Delta,DeltaLN,Lipkin,DN2,               &
    &                                                  ConstrainDispersion,Prec)
    !---------------------------------------------------------------------------
    ! Solves the HFB problem using (conjugate) gradient technique.
    ! The idea is to vary the HFB state | Phi \rangle in the space spanned by
    ! the Thouless transformation of | \Phi.
    !
    ! Computationally, this amounts to minimising the HFB energy within this
    ! space with constraints on the average particle number ( and on the
    ! dispersion when LN is active.)
    !---------------------------------------------------------------------------
    1 format ('----------------------------------------------')
    2 format (' WARNING: Gradient HFBsolver did not converge.')
    3 format (' Particles: '      , 2f15.5)
    4 format (' Norm of gradient:', e15.7)

    implicit none

    !---------------------------------------------------------------------------
    ! Variables for the FindFermi API.
    real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
    complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
    complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
    logical, intent(in)                       :: Lipkin, ConstrainDispersion
    real(KIND=dp), intent(in)                 :: Prec, DN2(2)

    !---------------------------------------------------------------------------
    real(KIND=dp) :: gamma(2,2), maxstep, step, L20Norm(2), LN(2), Disp(2)
    real(KIND=dp) :: N20norm(2), par(2), slope,  OldFermi(2)
    real(KIND=dp) :: gradientnorm(2,2), oldnorm(2,2), PR(2,2), z(2)=0.0
    integer       :: i,j,P,it,N, Rzindex,iter, inneriter, first=1, succes(2)
    logical       :: converged(2)

    ! Make sure we search for the right number of particles
    Particles(1) = neutrons
    Particles(2) = protons

    !---------------------------------------------------------------------------
    ! Initialize the module
    oldnorm = 0.0d0 ; gradientnorm = 0.0d0
    ! Step size works typically. 0.025 would be too big and smaller is slower.
    step    = 0.020_dp ; OldFermi = 0.0d0
    converged = .false.
    if(.not. allocated(GradU)) then
        N = maxval(blocksizes)
        !-----------------------------------------------------------------------
        ! Get U and V from the HFB module
        call GetUandV()
        OldgradU = GradU
        OldgradV = GradV

        allocate(Direction(N,N,pindex,Iindex)); Direction= 0.0_dp
        allocate(Gradient(N,N,Pindex,Iindex)) ; Gradient = 0.0_dp
        allocate(OldGrad(N,N,Pindex,Iindex))  ; Gradient = 0.0_dp
        allocate(aN20(N,N,Pindex,Iindex))     ; aN20     = 0.0_dp
        allocate(OldDir(N,N,Pindex,Iindex))   ; OldDir   = 0.0_dp

        !-----------------------------------------------------------------------
        ! Point to the right routines
        if(SC) then
          H20            => H20_sig
          ConstructRho   => ConstructRho_sig
          ConstructKappa => ConstructKappa_sig
          GradUpdate     => GradUpdate_sig
          ortho          => ortho_sig
        else
          H20            => H20_nosig
          ConstructRho   => ConstructRho_nosig
          ConstructKappa => ConstructKappa_nosig
          GradUpdate     => GradUpdate_nosig
          ortho          => ortho_nosig
        endif
    endif

    !---------------------------------------------------------------------------
    ! Start of the iterative solver
    do iter=1,HFBIter
      !if(any(isnan(gradU)) .or. any(isnan(gradV))) then
      !  print *, iter
      !  call stp('iter')
      !endif
      !do it=1,2      
      !  do P=1,2
      !    print *, 'P, it', P, it
   ! call checkUandV(DBLE(GRADU(1:blocksizes(P,it), 1:blocksizes(P,it),P,it)), &
   ! &               DBLE(GradV(1:blocksizes(P,it), 1:blocksizes(P,it),P,it)))
   !     enddo
   !   enddo


      !-------------------------------------------------------------------------
      ! Construct the density and anomalous density matrix.
      ! Only necessary when LN is active, since they then contribute to the
      ! HFBHamiltonian
      do it=1,Iindex
        do P=1,Pindex
          N = blocksizes(P,it)
          RhoHFB(1:N,1:N,P,it)  = ConstructRho(  GradV(1:N,1:N,P,it),          &
          &                                      sigblocks(P,it))
          KappaHFB(1:N,1:N,P,it)= ConstructKappa(GradU(1:N,1:N,P,it),          &
          &                                      GradV(1:N,1:N,P,it),          &
          &                                      sigblocks(P,it))
        enddo
      enddo
      !-------------------------------------------------------------------------
      ! Construct the matrices in block form.
      call BlockHFBHamil(Delta, z, L2)
      !-------------------------------------------------------------------------
      ! Calculate the gradients H20, N20 and LN20
      do it=1,Iindex
        if(converged(it)) cycle ! Don't wast CPU cycles

        n20norm(it) = 0.0d0
        do P=1,Pindex
            N = blocksizes(P,it)
            ! Save the old gradient
            Oldgrad(1:N,1:N,p,it)  = Gradient(1:N,1:N,P,it)

            Gradient(1:N,1:N,P,it) = H20(GradU(1:N,1:N,P,it),                  &
            &                            GradV(1:N,1:N,P,it),                  &
            &                            hblock(1:N,1:N,P,it),                 &
            &                            Dblock(1:N,1:N,P,it),                 &
            &                            sigblocks(P,it))

            aN20(1:N,1:N,P,it) = N20(GradU(1:N,1:N,P,it),GradV(1:N,1:N,P,it),  &
            &                                           sigblocks(P,it))
            do j=1,N
              do i=1,N
                N20norm(it) = N20norm(it) + aN20(i,j,P,it) * aN20(i,j,P,it)
              enddo
            enddo
        enddo
      enddo
      !-------------------------------------------------------------------------
      ! Readjust Fermi and L2
      succes = 0
      if(Lipkin .and. mod(iter,1).eq.1) then
         LN = Lncr8(Delta,DeltaLN,succes)
      elseif(ConstrainDispersion) then 
         Disp = Dispersion()      
      endif
      if(any(succes.ne.0)) call stp("lncr8 failed")
      do it=1,Iindex
          par(it) = 0.0
          do P=1,Pindex
            N = blocksizes(P,it)
            do i=1,N
              do j=1,N
                par(it) = par(it) + GradV(i,j,P,it) * GradV(i,j,P,it)
              enddo
            enddo
          enddo
          ! This is why we don't include the Fermi contribution directly in the
          ! HFBHamiltonian. In this way, we get more control over the conjugate
          ! descent in the Fermi direction and (hopefully) faster convergence.
          !print *, 'N20Norm', N20Norm
          if(abs(N20Norm(it)).gt.1d-10) then
            Fermi(it) = Fermi(it) + (Particles(it) - par(it))/N20norm(it)
          else
            Fermi(it) = Fermi(it) + 0.1*(Particles(it) - par(it))
          endif
          ! We don't take this luxury for the L2 variable, as it is kind of a
          ! nonsense to add it anyway.
          if(Lipkin) then
            L2(it) = L2(it) + 0.1*(LN(it)   - L2(it))
          elseif(ConstrainDispersion) then
            L2(it) = L2(it) - 0.01*(Disp(it) - DN2(it))
!            if(it.eq.1) print *, 'iter, it',it, L2(it), Disp(it), DN2(it)
          endif
      enddo
      !-------------------------------------------------------------------------
      ! Construct the total gradient
      do it=1,Iindex
          if(converged(it)) cycle
          do P=1,Pindex
              N = blocksizes(P,it)
              Gradient(1:N,1:N,P,it) = Gradient(1:N,1:N,P,it)                  &
              &                                 - Fermi(it) * aN20(1:N,1:N,P,it)
              ! Calculate the norm of the new gradient to
              ! 1) check for convergence and
              ! 2) calculate the conjugate direction
              ! PR is the scalar product between this gradient and the old one,
              ! it is necessary to form the Polak-Ribière formula for nonlinear
              ! conjugate gradient descent.
              gradientnorm(P,it) = 0.0 ; Pr(P,it) = 0.0_dp
              do i=1,N
                do j=1,N
                  gradientnorm(p,it) = gradientnorm(P,it) +                    &
                  &                      Gradient(i,j,P,it) * Gradient(j,i,P,it)
                  Pr(P,it) = Pr(P,it) +                                        &
                  &                      Gradient(i,j,P,it) * OldGrad(j,i,P,it)
                enddo
              enddo
          enddo

          !---------------------------------------------------------------------
          ! Check if this isospin is converged or not. We check on
          ! 1) the gradient is small enough
          ! 2) the number of particles is close enough to the desired value
          ! 3) the Lambda_2 value is close enough to the calculated one
          !    (if Lipkin is active)
          if(all(sqrt(abs(GradientNorm(:,it))).lt.prec) .and.                      &
          &     abs(Par(it) - Particles(it)).lt.Prec  ) then
           converged(it)=.true.
           if(Lipkin) then
             if(abs(L2(it) - LN(it)).gt.Prec ) then
               converged(it) = .false.
             endif
           endif
         endif
      enddo
      !-------------------------------------------------------------------------
      ! Construct the conjugate direction and update the U and V matrices.
      do it=1,Iindex
          if(converged(it)) cycle
          do P=1,Pindex
              N = blocksizes(P,it)
              Direction(1:N,1:N,P,it) = Gradient(1:N,1:N,P,it)
              if(oldnorm(P,it).ne.0.0d0) then
                !---------------------------------------------------------------
                ! Polak-Ribière formula for the conjugate gradient. This
                ! outperforms the standard formula in my tests for heavy nuclei
                ! (226Ra) but is worse for light nuclei (24Mg). However, since
                ! the computational burden for these light nuclei is so small
                ! already, I don't care for them.
                gamma(P,it) = 0.0_dp != gradientnorm(p,it)/oldnorm(p,it)
                !gamma(P,it) = gamma(P,it) - PR(p,it)/oldnorm(p,it)
                !---------------------------------------------------------------
                ! Update
                !Direction(1:N,1:N,P,it) = Direction(1:N,1:N,P,it)  +           &
                !&                             gamma(P,it) * OldDir(1:N,1:N,P,it)
              endif
              !-----------------------------------------------------------------
              ! Check if this direction is indeed a descent direction.
              ! In a linear problem, this is guaranteed, but I'm not too sure
              ! about this case. So, to be safe we check.
              ! If it is not a descent direction, we take the gradient.
              slope = 0.0d0
              do i=1,N
                do j=1,N
                  slope = slope + Direction(i,j,P,it) * Gradient(j,i,P,it)
                enddo
              enddo

              if(slope.gt.0) then
                Direction(1:N,1:N,P,it) = Gradient(1:N,1:N,P,it)
                !call stp('Not a descent direction', 'slope', slope)
              endif
              !-----------------------------------------------------------------
              oldgradU(1:N,1:N,P,it) = GradU(1:N,1:N,P,it)
              oldgradV(1:N,1:N,P,it) = GradV(1:N,1:N,P,it)
              if(any(isnan(gradient))) call stp('grad')
              !-----------------------------------------------------------------
              ! Line search for a good step.
              ! In practice, this does nothing but slow the process in my
              ! experience. Fixed-step for some reason works fine.
              call Linesearch(oldgradU(1:N,1:N,P,it), oldgradV(1:N,1:N,P,it),  &
              &               Gradient(1:N,1:N,P,it), Direction(1:N,1:N,P,it), &
              &               step,   gradU(1:N,1:N,P,it), gradV(1:N,1:N,P,it),&
              &             hblock(1:N,1:N,P,it),Dblock(1:N,1:N,P,it),Fermi(it)&
              &            ,Sigblocks(P,it))

              OldDir(1:N,1:N,P,it)  = Direction(1:N,1:N,P,it)
              oldnorm(P,it)         = gradientnorm(P,it)
          enddo
      enddo
      !-------------------------------------------------------------------------
      ! Detect convergence or divergence?
      !print *, 'gradientnorm', gradientnorm
      if(any(.not. converged) .and. iter.eq.HFBIter) then
        print 1
        print 2
        print 3, Par
        print 4, sum(abs(GradientNorm))
        print 1
      endif
      if(all(converged)) exit
   enddo
   ! Make the HFB module happy again and hand control back.
   call OutHFBModule
   call CalcQPEnergies(GradU,gradV,Fermi,L2,Delta)
 end subroutine HFBFermiGradient

  subroutine Linesearch(OldU, OldV, Grad, Direction, maxstep, Ulim, Vlim, hlim,&
    &                   Dlim, Fermi,S)
    !---------------------------------------------------------------------------
    ! This linesearch algorithm is deprecated: it failed to do better than a
    ! fixed step algorithm.
    ! So, if anyone wants to try again, feel free, the API is there.
    !---------------------------------------------------------------------------
    ! What does not work:
    !   *) backtracking line-search with Wolfe Conditions on the energy
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: OldU(:,:), OldV(:,:), hlim(:,:), Dlim(:,:)
    real(KIND=dp), intent(in) :: Grad(:,:), Fermi
    real(KIND=dp), intent(out):: Ulim(:,:), Vlim(:,:)
    real(KIND=dp), intent(inout) :: maxstep,Direction(:,:)
    integer, intent(in)       :: S

    call GradUpdate(maxstep, oldU,oldV,Direction,Ulim,Vlim,S)

  end subroutine Linesearch

  function ConstructRho_nosig(Vlim,S) result(Rho)
    !---------------------------------------------------------------------------
    ! Construct the part of the density matrix when signature is not conserved.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) ::Vlim(:,:)
    integer, intent(in)       :: S
    real(KIND=dp)             :: Rho(size(Vlim,1),size(Vlim,1))

    integer :: N,i,j,k

    N = size(Vlim,1)
    rho = 0.0d0
    do j=1,N
      do i=j,N
        do k=1,N
          Rho(i,j) = Rho(i,j) + Vlim(i,k) * Vlim(j,k)
        enddo
        Rho(j,i) = Rho(i,j)
      enddo
    enddo

  end function ConstructRho_nosig

  function ConstructRho_sig(Vlim,S) result(Rho)
    !---------------------------------------------------------------------------
    ! Construct the part of the density matrix when signature is conserved.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Vlim(:,:)
    real(KIND=dp)             :: Rho(size(Vlim,1),size(Vlim,1))
    integer, intent(in)       :: S
    integer :: N,i,j,k

    N = size(Vlim,1)

    rho = 0.0d0
    do j=1,N/2
      do i=j,N/2
        do k=S+1,N
          Rho(i,j)  = Rho(i,j)+ Vlim(i,k) * Vlim(j,k)
        enddo
        Rho(j,i)         = Rho(i,j)
      enddo
    enddo

    do j=N/2+1,N
      do i=j,N
        do k=1,S
          Rho(i,j)  = Rho(i,j)+ Vlim(i,k) * Vlim(j,k)
        enddo
        Rho(j,i)         = Rho(i,j)
      enddo
    enddo

  end function ConstructRho_sig

  function ConstructKappa_nosig(Ulim,Vlim,S) result(Kappa)
    !---------------------------------------------------------------------------
    ! Construct the part of kappa when signature is broken.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Ulim(:,:), Vlim(:,:)
    real(KIND=dp),allocatable :: Kappa(:,:)
    integer, intent(in)       :: S

    integer :: N,i,j,k
    if(.not.allocated(Kappa)) then
      N = maxval(blocksizes)
      allocate(Kappa(N,N))
    endif

    N = size(Vlim,1)
    Kappa = 0.0d0
    do j=1,N
      do i=j+1,N
        do k=1,N
          Kappa(i,j) = Kappa(i,j) + Vlim(i,k) * Ulim(j,k)
        enddo
        Kappa(j,i) = - Kappa(i,j)
      enddo
    enddo

  end function ConstructKappa_nosig

  function ConstructKappa_sig(Ulim,Vlim,S) result(Kappa)
    !---------------------------------------------------------------------------
    ! Construct the part of kappa when signature is conserved.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Ulim(:,:), Vlim(:,:)
    real(KIND=dp),allocatable :: Kappa(:,:)
    integer, intent(in)       :: S

    integer :: N,i,j,k
    if(.not.allocated(Kappa)) then
      N = maxval(blocksizes)
      allocate(Kappa(N,N))
    endif

    N = size(Vlim,1)
    Kappa = 0.0d0
    do j=1,N/2
      do i=N/2+1,N
        do k=1,S
          Kappa(i,j)         = Kappa(i,j) + Vlim(i,k) * Ulim(j,k)
        enddo
        Kappa(j,i) = - Kappa(i,j)
      enddo
    enddo

  end function ConstructKappa_sig

  subroutine GradUpdate_nosig(step, U1,V1,Grad,U2,V2, S)
    !---------------------------------------------------------------------------
    ! Update U and V with a gradient step of size 'step'.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: U1(:,:), V1(:,:), step
    real(KIND=dp), intent(in) :: Grad(:,:)
    integer, intent(in)       :: S
    real(KIND=dp),intent(out):: U2(:,:), V2(:,:)
    integer                   :: i,j,P,it,N,k

    N = size(U1,1)

    do j=1,N
        do i=1,N
            U2(i,j) = U1(i,j) ; V2(i,j) = V1(i,j)
            do k=1,N
                U2(i,j) = U2(i,j) - step * V1(i,k)*Grad(k,j)
                V2(i,j) = V2(i,j) - step * U1(i,k)*Grad(k,j)
            enddo
        enddo
    enddo
    !Don't forget to orthonormalise
    call ortho(U2,V2,S)
   
  end subroutine GradUpdate_nosig

  subroutine GradUpdate_sig(step, U1,V1,Grad,U2,V2, S)
    !---------------------------------------------------------------------------
    ! Update U and V with a gradient step of size 'step'.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: U1(:,:), V1(:,:), step
    real(KIND=dp), intent(in) :: Grad(:,:)
    real(KIND=dp),intent(out) :: U2(:,:), V2(:,:)
    integer, intent(in)       :: S
    integer                   :: i,j,P,it,N,k

    N = size(U1,1)

    U2(1:N/2  ,1:S)   = U1(1:N/2  ,1:S)  - step *                               &
    &                     matmul(V1(1:N/2, S+1:N), Grad(S+1:N,1:S))
    U2(1+N/2:N,S+1:N) = U1(1+N/2:N,S+1:N) - step *                             &
    &                     matmul(V1(N/2+1:N, 1:S), Grad(1:S,S+1:N))
    V2(1:N/2  ,S+1:N) = V1(1:N/2  ,S+1:N) - step *                             &
    &                      matmul(U1(1:N/2, 1:S), Grad(1:S,S+1:N))
    V2(N/2+1:N ,1:S)  = V1(N/2+1:N ,1:S) - step *                              &
    &                      matmul(U1(N/2+1:N, S+1:N), Grad(S+1:N,1:S))

    !Don't forget to orthonormalise
    call ortho(U2,V2,S)

  end subroutine GradUpdate_sig

function H20_nosig(Ulim,Vlim,hlim,Dlim,S) result(H20)
    !-----------------------------------------------------------
    ! Get the gradient of the HFBHamiltonian H20.
    !
    !
    ! In abusive notation
    !
    ! H20 =
    !       U_1^{dagger} h_+ V_2 + U_1^{\dagger} \Delta U2
    !     - V_1^{dagger} h_+ U_2 - V_1^{\dagger} \Delta V2
    !-----------------------------------------------------------
    real(KIND=dp), intent(in)  :: Ulim(:,:), Vlim(:,:), dlim(:,:), hlim(:,:)
    real(KIND=dp)              :: H20(size(Ulim,1),size(Ulim,1))
    integer, intent(in)        :: S

    !-----------------------------------------------------
    ! Auxiliary arrays
    real(KIND=dp):: DsV(size(Ulim,1),size(Ulim,1))
    real(KIND=dp):: DsU(size(Ulim,1),size(Ulim,1))
    real(KIND=dp):: hV(size(Ulim,1),size(Ulim,1))
    real(KIND=dp):: hU(size(Ulim,1),size(Ulim,1))

    integer                    :: N, i,j,k

    H20 = 0.0_dp
    DsU = 0.0_dp
    DsV = 0.0_dp
    hU = 0.0_dp
    hV = 0.0_dp

    N = size(Ulim,1)
    !-----------------------------------------------------
    ! Construct auxiliary arrays
    do j=1,N
        do i=1,N
            hV(i,j) = 0.0_dp ; hU(i,j) = 0.0_dp
            do k=1,N
                hV(i,j) = hV(i,j) + hlim(i,k) * Vlim(k,j)
                hU(i,j) = hU(i,j) + hlim(i,k) * Ulim(k,j)
            enddo
            DsV(i,j) = 0.0_dp ; DsU(i,j) = 0.0
            do k=1,N
                DsV(i,j) = DsV(i,j) + Dlim(i,k) * Vlim(k,j)
                DsU(i,j) = DsU(i,j) + Dlim(i,k) * Ulim(k,j)
            enddo
        enddo
    enddo
    !-------------------------------------------------------
    ! Construct H20
     do j=1,N
         do i=j+1,N
             H20(i,j) = 0.0_dp
             do k=1,N
                 H20(i,j) = H20(i,j)  +  Ulim(k,i) * hV (k,j) &
                 &                    -  Vlim(k,i) * hU (k,j) &
                 &                    -  Vlim(k,i) * DsV(k,j) &
                 &                    +  Ulim(k,i) * DsU(k,j)
             enddo
             ! H20 is antisymmetric
             H20(j,i) = -H20(i,j)
         enddo
     enddo
  end function H20_nosig

  function H20_sig(Ulim,Vlim,hlim,Dlim,S) result(H20)
    !-----------------------------------------------------------
    ! Get the gradient of the HFBHamiltonian H20.
    !
    !
    ! In abusive notation
    !
    ! H20 =
    !       U_1^{dagger} h_+ V_2 + U_1^{\dagger} \Delta U2
    !     - V_1^{dagger} h_+ U_2 - V_1^{\dagger} \Delta V2
    !-----------------------------------------------------------
    integer, intent(in)        :: S
    real(KIND=dp), intent(in)  :: Ulim(:,:), Vlim(:,:), dlim(:,:), hlim(:,:)
    real(KIND=dp)              :: H20(size(Ulim,1),size(Ulim,1))
    !-----------------------------------------------------
    ! Auxiliary arrays
    real(KIND=dp):: DsV(size(Ulim,1),size(Ulim,1))
    real(KIND=dp):: DsU(size(Ulim,1),size(Ulim,1))
    real(KIND=dp):: hV(size(Ulim,1) ,size(Ulim,1))
    real(KIND=dp):: hU(size(Ulim,1) ,size(Ulim,1))
    integer   :: N, i,j,k, of

    H20 = 0.0_dp

    N = size(Ulim,1)
    of = N/2
    !---------------------------------------------------------------------------
    ! Construct auxiliary arrays
    !
    ! Note to the reader:
    ! The arrays are only calculated half and the other half CAN NOT BE OBTAINED
    ! BY SYMMETRY.  HOWEVER, H20 itself is antisymmetric and it is
    ! this symmetry that allows us to only calculate half of each matrix.
    hU(1:N/2,1:S)      = matmul(hlim(1:N/2,1:N/2), Ulim(1:N/2,1:S))
    hV(N/2+1:N,1:S)    = matmul(hlim(N/2+1:N,N/2+1:N), Vlim(N/2+1:N,1:S))
    DsV(1:N/2,1:S)     = matmul(Dlim(1:N/2,N/2+1:N), Vlim(N/2+1:N,1:S))
    DsU(N/2+1:N,1:S)   = matmul(Dlim(N/2+1:N,1:N/2), Ulim(1:N/2,1:S))

    !---------------------------------------------------------------------------
    ! Construct H20
    H20(S+1:N, 1:S) =                                                          &
    &                    matmul(transpose(Ulim(N/2+1:N,S+1:N)),hV(N/2+1:N, 1:S))
    H20(S+1:N, 1:S) = H20(S+1:N, 1:S) -                                        &
    &                    matmul(transpose(Vlim(1:N/2,S+1:N)),  DsV(1:N/2,1:S))
    H20(S+1:N, 1:S) = H20(S+1:N, 1:S) +                                        &
    &                   matmul(transpose(Ulim(N/2+1:N,S+1:N)),DsU(N/2+1:N, 1:S))
    H20(S+1:N, 1:S) = H20(S+1:N, 1:S) -                                        &
    &                       matmul(transpose(Vlim(1:N/2,S+1:N)), hU(1:N/2, 1:S))
    do j=1,S
      do i=S+1,N
        H20(j,i) = - H20(i,j)
      enddo
    enddo
    !---------------------------------------------------------------------------
  end function H20_sig

  function N20(Ulim,Vlim,S)
    !---------------------------------------------------------------------------
    ! Construct the gradient of the energy with respect to Lambda, N_20.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: Ulim(:,:), Vlim(:,:)
    real(KIND=dp),allocatable  :: N20(:,:)
    integer                    :: i,j,N,k
    integer, intent(in)        :: S

    if(.not.allocated(N20)) then
      N = size(Ulim,1)
      allocate(N20(N,N)) ; N20 = 0.0d0
    endif
    N = size(Ulim,1)
    N20 = 0.0d0
    if(SC) then

      N20(1:S, S+1:N) =                                                        &
      &            matmul(transpose(Ulim(1:N/2, 1:S)),Vlim(1:N/2,S+1:N))       &
      &          - matmul(transpose(Vlim(N/2+1:N,1:S)),Ulim(N/2+1:N,S+1:N))

      do j=S+1,N
        do i=1,S
          N20(j,i) = - N20(i,j)
        enddo
      enddo
    else
      do j=1,N
        do i=j+1,N
          N20(i,j) = 0.0d0
          do k=1,N
            N20(i,j) = N20(i,j) + Ulim(k,i) * Vlim(k,j) - Vlim(k,i) * Ulim(k,j)
          enddo
          N20(j,i) = - N20(i,j)
        enddo
      enddo
    endif

  end function N20

  subroutine ortho_nosig(U,V, S)
        !-----------------------------------------------------------------------
        ! Orthogonalises U and V vectors when signature is broken.
        ! Uses a simple Gramm-Schmidt routine.
        !-----------------------------------------------------------------------
        real(KIND=dp), intent(inout) :: U(:,:), V(:,:)
        real(KIND=dp)                :: overlap
        integer, intent(in)          :: S

        integer :: i,j,k, N

        N = size(U,1)
        do j=1,N
          !-------------------------------------------------------
          ! Normalise vector ( U_j )
          !                  ( V_j )
          overlap = 0.0_dp
          do i=1,N
              Overlap = Overlap + U(i,j)**2 + V(i,j)**2
          enddo
          Overlap = 1.0_dp/sqrt(overlap)
          U(1:N,j) = U(1:N,j)*overlap
          V(1:N,j) = V(1:N,j)*overlap
          !-------------------------------------------------------
          ! Orthogonalise all the rest against vector j
          do i=j+1,N
              Overlap = 0.0_dp
              do k=1,N
                  Overlap = Overlap + U(k,j) * U(k,i)
                  Overlap = Overlap + V(k,j) * V(k,i)
              enddo
              U(1:N,i) = U(1:N,i) - Overlap * U(1:N,j)
              V(1:N,i) = V(1:N,i) - Overlap * V(1:N,j)
          enddo
        enddo

  end subroutine ortho_nosig

  subroutine ortho_sig(U,V,S)
    !-----------------------------------------------------------------------
    ! Orthogonalises U and V vectors when signature is conserved.
    ! Uses a simple Gramm-Schmidt routine.
    !-----------------------------------------------------------------------
    ! S is the number of positive signature eigenvectors (U, V)
    !-----------------------------------------------------------------------
      real(KIND=dp), intent(inout) :: U(:,:), V(:,:)
      integer, intent(in)          :: S
      real(KIND=dp)                :: overlap

      integer :: i,j,k, N

      N = size(U,1)
      do j=1,S
        !-------------------------------------------------------
        ! Normalise vector ( U_j )
        !                  ( V_j )
        overlap = 0.0_dp
        do i=1,N/2
            Overlap = Overlap + U(i,j)**2
        enddo
        do i=N/2+1,N
            Overlap = Overlap + V(i,j)**2
        enddo

        Overlap      = 1.0_dp/sqrt(overlap)
        U(1:N/2,j)   = U(1:N/2,j)  *overlap
        V(N/2+1:N,j) = V(N/2+1:N,j)*overlap
        !-------------------------------------------------------
        ! Orthogonalise all the rest against vector j
        do i=j+1,N/2
            Overlap = 0.0_dp
            do k=1,N/2
                Overlap = Overlap + U(k,j) * U(k,i)
            enddo
            do k=N/2+1,N
                Overlap = Overlap + V(k,j) * V(k,i)
            enddo
            U(1:N/2,i)   = U(1:N/2  ,i) - Overlap * U(1:N/2,j)
            V(N/2+1:N,i) = V(N/2+1:N,i) - Overlap * V(N/2+1:N,j)
        enddo
      enddo

      do j=S+1,N
        !-------------------------------------------------------
        ! Normalise vector ( U_j )
        !                  ( V_j )
        overlap = 0.0_dp
        do i=1,N/2
            Overlap = Overlap + V(i,j)**2
        enddo
        do i=N/2+1,N
            Overlap = Overlap + U(i,j)**2
        enddo

        Overlap = 1.0_dp/sqrt(overlap)
        U(N/2+1:N,j)   = U(N/2+1:N,j)  *overlap
        V(1:N/2,j)     = V(1:N/2,j)*overlap
        !-------------------------------------------------------
        ! Orthogonalise all the rest against vector j
        do i=j+1,N
            Overlap = 0.0_dp
            do k=1,N/2
                Overlap = Overlap + V(k,j) * V(k,i)
            enddo
            do k=N/2+1,N
                Overlap = Overlap + U(k,j) * U(k,i)
            enddo
            U(N/2+1:N,i) = U(N/2+1:N,i) - Overlap * U(N/2+1:N,j)
            V(1:N/2,i)   = V(1:N/2,i)   - Overlap * V(1:N/2,j)
        enddo
      enddo
  end subroutine ortho_sig

  subroutine CalcQPEnergies(Ulim,Vlim,Fermi,L2,Delta)
    !---------------------------------------------------------------------------
    ! Calculate the Quasiparticle energies.
    !
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Ulim(:,:,:,:), Vlim(:,:,:,:)
    real(KIND=dp), intent(in) :: Fermi(2), L2(2)
    complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
    real(KIND=dp)             :: Vector(2*nwt,nwt), HV(2*nwt,nwt)
    integer                   :: P ,it, N,jj,i,j,k

    QuasiEnergies = 0.0_dp

    ! Construct the HFBHamiltonian
    call ConstructHFBHamiltonian(Fermi,Delta,L2,0.0_dp)

    do it=1,Iindex
        do P=1,Pindex
            N = blocksizes(P,it)
            ! Construct the eigenvector.
            Vector(  1:N  ,1:N) = Ulim(1:N,1:N,P,it)
            Vector(N+1:2*N,1:N) = Vlim(1:N,1:N,P,it)
            HV(1:2*N,1:N) = 0.0_dp
            do j=1,N
                ! Multiply every vector j with the HFB Hamiltonian
                do i=1,2*N
                    do k=1,2*N
                        HV(i,j) = HV(i,j) +                                    &
                        &                DBLE(HFBHamil(i,k,P,it)) *  Vector(k,j)
                    enddo
                enddo
            enddo
            do j=1,N
                jj = HFBColumns(j,P,it)
                QuasiEnergies(jj,P,it) = 0.0_dp
                do i=1,2*N
                    QuasiEnergies(jj,P,it) = QuasiEnergies(jj,P,it)            &
                    &                                    + Vector(i,j) * HV(i,j)
                enddo
            enddo
        enddo
    enddo
end subroutine CalcQPEnergies

subroutine checkUandV(U,V)
        !-----------------------------------------------------------------------
        ! Debugging and sanity check tool. 
        ! Check if the matrices U and V satisfy the relations to define a  
        ! correct Bogolyubov transformation, i.e. that
        !
        !      ( V^* U )
        ! W =  (       )
        !      ( U^* V )
        !
        ! is a unitary matrix, or equivalently that 
        !
        ! U^\dagger U + V^\dagger V = 1 & U U^\dagger + V^* V^T = 1
        ! U^T V       + V^T       U = 0 & U V^\dagger + V^* U^T = 0
        !
        !-----------------------------------------------------------------------
        1 format('===============================')
        2 format(' U^T U + V^T V = 1')
        3 format(' U^T V + V^T U = 0')
        4 format(' U U^T + V V^T = 1')
        5 format(' U V^T + V U^T = 0')
    
        real*8, intent(in) :: U(:,:),V(:,:)
        real*8  :: check(size(U,1),size(U,1),4)
        integer :: i,j,N
        
        check(:,:,1) = matmul(transpose(U), U) + matmul(transpose(V), V)
        check(:,:,2) = matmul(transpose(U), V) + matmul(transpose(V), U)
        check(:,:,3) = matmul( U,transpose(U)) + matmul( V,transpose(V))
        check(:,:,4) = matmul( U,transpose(V)) + matmul( V,transpose(U))
        
        print *
        print 1
        N = size(U,1)
        print 2
        !do i=1,N
        !    print ('(150f7.3)'), check(i,:,1)
        !enddo
        do j=1,N
                do i=1,N  
                        if(i.eq.j) cycle      
                        if(abs(check(i,j,1)).gt.1d-10) then
                                print *, i,j, 'check1', check(i,j,1)
                        endif
                enddo
        enddo

        print *
        print 3
        !do i=1,N
        !    print ('(150f7.3)'), check(i,:,2)
        !enddo

        do j=1,N
                do i=1,N  
                        if(abs(check(i,j,2)).gt.1d-10) then
                                print *, i,j, 'check2', check(i,j,2)
                        endif
                enddo
        enddo


        print *
        print 4
        !do i=1,N
        !    print ('(150f7.3)'), check(i,:,3)
        !enddo
        do j=1,N
                do i=1,N  
                        if(i.eq.j) cycle
                        if(abs(check(i,j,3)).gt.1d-10) then
                                print *, i,j, 'check3', check(i,j,3)
                        endif
                enddo
        enddo


        print *
        print 5
        !do i=1,N
        !    print ('(150f7.3)'), check(i,:,4)
        !enddo
        do j=1,N
                do i=1,N  
                        if(abs(check(i,j,4)).gt.1d-10) then
                                print *, i,j, 'check2', check(i,j,4)
                        endif
                enddo
        enddo


        print 1

    end subroutine checkUandV

end module GradientHFB
