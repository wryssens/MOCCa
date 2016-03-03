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
  ! Matrices that I don't want to reallocate everytime.
  real(KIND=dp),allocatable  :: Gradient(:,:,:,:) ,OldGrad(:,:,:,:)
  real(KIND=dp),allocatable  :: aN20(:,:,:,:)
  real(KIND=dp),allocatable  :: oldgradU(:,:,:,:) ,oldgradV(:,:,:,:)
  real(KIND=dp),allocatable  :: Direction(:,:,:,:), OldDir(:,:,:,:)

  procedure(H20_sig), pointer            :: H20
  procedure(constructrho_sig), pointer   :: constructrho
  procedure(constructkappa_sig), pointer :: constructkappa
  procedure(gradupdate_sig), pointer     :: gradupdate
  procedure(ortho_sig), pointer          :: ortho
  ! procedure(LN20_sig), pointer           :: LN20
contains

  subroutine GetUandV
    !----------------------------------------
    ! Get U and V values from the HFB module.
    !
    !----------------------------------------

    integer :: it,P,i,j,N, ind(2,2), Rzindex,k, S
    logical :: incolumns

    !------------------------
    ! Initialize the matrices.
    !-------------------------
    N = maxval(blocksizes)

    allocate(GradU(N,N,Pindex,Iindex)) ; GradU = 0.0_dp
    allocate(GradV(N,N,Pindex,Iindex)) ; GradV = 0.0_dp

    if(all(HFBColumns.eq.0)) then

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
                  ind(P,it) = ind(P,it)+ 1
                endif
            enddo
        enddo
    enddo

  end subroutine GetUandV

  subroutine OutHFBModule
      !-------------------------------------------------------------------------
      ! Put our solution into the HFB module correctly.
      !-------------------------------------------------------------------------
      ! TODO: put the conjugate states too
      !-------------------------------------------------------------------------
      integer :: i,j,P,it,N

      do it=1,Iindex
          do P=1,Pindex
              N = blocksizes(P,it)
              do i=1,N
                  do j=1,N
                      ! Don't forget to point HFBcolumns the correct way!
                      HFBColumns(j,P,it) = j + N
                      V(i,j+N,P,it) = GradV(i,j,P,it)
                      U(i,j+N,P,it) = GradU(i,j,P,it)
                  enddo
              enddo
              RhoHFB(1:N,1:N,P,it)   = ConstructRho(GradV(1:N,1:N,P,it))
              KappaHFB(1:N,1:N,P,it) = ConstructKappa(GradU(1:N,1:N,P,it),     &
              &                                       GradV(1:N,1:N,P,it))
          enddo
      enddo

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

  subroutine HFBFermiGradient(Fermi,L2,Delta,DeltaLN,Lipkin,Prec)
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
    3 format (' Particles: ' , 2f10.5)
    4 format (' Norm of gradient:', f10.7)

    implicit none

    !---------------------------------------------------------------------------
    ! Variables for the FindFermi API.
    real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
    complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
    complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
    logical, intent(in)                       :: Lipkin
    real(KIND=dp), intent(in)                 :: Prec

    !---------------------------------------------------------------------------
    real(KIND=dp) :: gamma(2,2), maxstep, step, L20Norm(2), LN(2)
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
    step = 0.020_dp ; OldFermi = 0.0d0
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
          H20 => H20_sig
          ConstructRho   => ConstructRho_sig
          ConstructKappa => ConstructKappa_sig
          GradUpdate     => GradUpdate_sig
          ortho          => ortho_sig
        else
          H20 => H20_nosig
          ConstructRho => ConstructRho_nosig
          ConstructKappa => ConstructKappa_nosig
          GradUpdate     => GradUpdate_nosig
          ortho          => ortho_nosig
        endif
    endif

    !---------------------------------------------------------------------------
    ! Start of the iterative solver
    do iter=1,HFBIter

      !-------------------------------------------------------------------------
      ! Construct the density and anomalous density matrix.
      ! Only necessary when LN is active, since they then contribute to the
      ! HFBHamiltonian
      do it=1,Iindex
        do P=1,Pindex
          N = blocksizes(P,it)
          RhoHFB(1:N,1:N,P,it)  = ConstructRho(  GradV(1:N,1:N,P,it))
          KappaHFB(1:N,1:N,P,it)= ConstructKappa(GradU(1:N,1:N,P,it),          &
          &                                      GradV(1:N,1:N,P,it))
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
        !l20norm(it) = 0.0d0
        do P=1,Pindex
            N = blocksizes(P,it)
            ! Save the old gradient
            Oldgrad(1:N,1:N,p,it)  = Gradient(1:N,1:N,P,it)

            Gradient(1:N,1:N,P,it) = H20(GradU(1:N,1:N,P,it),                  &
            &                            GradV(1:N,1:N,P,it),                  &
            &                            hblock(1:N,1:N,P,it),                 &
            &                            Dblock(1:N,1:N,P,it))

            aN20(1:N,1:N,P,it) = N20(GradU(1:N,1:N,P,it),GradV(1:N,1:N,P,it))

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
      if(Lipkin .and. mod(iter,10).eq.1) LN = Lncr8(Delta,DeltaLN,succes)
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
          if(abs(N20Norm(it)).gt.1d-10) then
            Fermi(it) = Fermi(it) + (Particles(it) - par(it))/N20norm(it)
          else
            Fermi(it) = Fermi(it) + 0.1*(Particles(it) - par(it))
          endif
          ! We don't take this luxury for the L2 variable, as it is kind of a
          ! nonsense to add it anyway.
          if(Lipkin) then
            L2(it) = L2(it) + 0.1*(LN(it) - L2(it))
          endif
      enddo
      !-------------------------------------------------------------------------
      ! Construct the total gradient
      do it=1,Iindex
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
          if(all((abs(GradientNorm(:,it))).lt.prec) .and.                      &
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
                gamma(P,it) = gradientnorm(p,it)/oldnorm(p,it)
                gamma(P,it) = gamma(P,it) - PR(p,it)/oldnorm(p,it)
                !---------------------------------------------------------------
                ! Update
                Direction(1:N,1:N,P,it) = Direction(1:N,1:N,P,it)  +           &
                &                             gamma(P,it) * OldDir(1:N,1:N,P,it)
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

              !-----------------------------------------------------------------
              ! Line search for a good step.
              ! In practice, this does nothing but slow the process in my
              ! experience. Fixed-step for some reason works fine.
              call Linesearch(oldgradU(1:N,1:N,P,it), oldgradV(1:N,1:N,P,it),  &
              &               Gradient(1:N,1:N,P,it), Direction(1:N,1:N,P,it), &
              &               step,   gradU(1:N,1:N,P,it), gradV(1:N,1:N,P,it),&
              &             hblock(1:N,1:N,P,it),Dblock(1:N,1:N,P,it),Fermi(it))

              OldDir(1:N,1:N,P,it)  = Direction(1:N,1:N,P,it)
              oldnorm(P,it)         = gradientnorm(P,it)
          enddo
      enddo
      !-------------------------------------------------------------------------
      ! Detect convergence or divergence?
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
 end subroutine HFBFermiGradient

  subroutine Linesearch(OldU, OldV, Grad, Direction, maxstep, Ulim, Vlim, hlim,&
    &                   Dlim, Fermi)
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

    call GradUpdate(maxstep, oldU,oldV,Direction,Ulim,Vlim)

  end subroutine Linesearch

  function ConstructRho_nosig(Vlim) result(Rho)
    !---------------------------------------------------------------------------
    ! Construct the part of the density matrix when signature is not conserved.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) ::Vlim(:,:)
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

  function ConstructRho_sig(Vlim) result(Rho)
    !---------------------------------------------------------------------------
    ! Construct the part of the density matrix when signature is conserved.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Vlim(:,:)
    real(KIND=dp)             :: Rho(size(Vlim,1),size(Vlim,1))

    integer :: N,i,j,k

    N = size(Vlim,1)
    rho = 0.0d0
    do j=1,N/2
      do i=j,N/2
        do k=1,N/2
          Rho(i,j)         = Rho(i,j)         + Vlim(i,k+N/2) * Vlim(j,k+N/2)
          Rho(i+N/2,j+N/2) = Rho(i+N/2,j+N/2) + Vlim(i+N/2,k) * Vlim(j+N/2,k)
        enddo
        Rho(j,i)         = Rho(i,j)
        Rho(j+N/2,i+N/2) = Rho(i+N/2,j+N/2)
      enddo
    enddo

  end function ConstructRho_sig

  function ConstructKappa_nosig(Ulim,Vlim) result(Kappa)
    !---------------------------------------------------------------------------
    ! Construct the part of kappa when signature is broken.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Ulim(:,:), Vlim(:,:)
    real(KIND=dp),allocatable :: Kappa(:,:)

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

  function ConstructKappa_sig(Ulim,Vlim) result(Kappa)
    !---------------------------------------------------------------------------
    ! Construct the part of kappa when signature is conserved.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Ulim(:,:), Vlim(:,:)
    real(KIND=dp),allocatable :: Kappa(:,:)

    integer :: N,i,j,k
    if(.not.allocated(Kappa)) then
      N = maxval(blocksizes)
      allocate(Kappa(N,N))
    endif

    N = size(Vlim,1)
    Kappa = 0.0d0
    do j=1,N/2
      do i=N/2+1,N
        do k=1,N/2
          Kappa(i,j)         = Kappa(i,j) + Vlim(i,k) * Ulim(j,k)
        enddo
        Kappa(j,i) = - Kappa(i,j)
      enddo
    enddo

  end function ConstructKappa_sig

  ! function Energy(hlim,Dlim,rho,kappa, Fermi)
  !   !--------------------------------------------------
  !   ! Compute the HFB energy associated with U and V.
  !   !
  !   ! E = Tr( h \rho ) - 0.5 * Tr(Delta \kappa) - \lambda <N>
  !   !
  !   !--------------------------------------------------
  !
  !   real(KIND=dp), intent(in) :: hlim(:,:), Dlim(:,:),rho(:,:), kappa(:,:)
  !   real(KIND=dp), intent(in) :: Fermi
  !   real(KIND=dp)             :: Energy
  !
  !   integer                   :: N, i,j
  !
  !   N = size(hlim,1)
  !   Energy = 0.0_dp
  !   do i=1,N
  !     do j=1,N
  !       Energy = Energy + hlim(i,j) * Rho(j,i) - 0.5*kappa(i,j)*Dlim(j,i)
  !     enddo
  !   Energy = Energy   - Fermi     * Rho(i,i)
  !   enddo
  ! end function Energy

  subroutine GradUpdate_nosig(step, U1,V1,Grad,U2,V2)
    !---------------------------------------------------------------------------
    ! Update U and V with a gradient step of size 'step'.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: U1(:,:), V1(:,:), step
    real(KIND=dp), intent(in) :: Grad(:,:)
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
    call ortho(U2,V2)

  end subroutine GradUpdate_nosig

  subroutine GradUpdate_sig(step, U1,V1,Grad,U2,V2)
    !---------------------------------------------------------------------------
    ! Update U and V with a gradient step of size 'step'.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: U1(:,:), V1(:,:), step
    real(KIND=dp), intent(in) :: Grad(:,:)
    real(KIND=dp),intent(out) :: U2(:,:), V2(:,:)
    integer                   :: i,j,P,it,N,k

    N = size(U1,1)

    U2(1:N/2  ,1:N/2)   = U1(1:N/2,1:N/2) - step *                             &
    &                     matmul(V1(1:N/2, N/2+1:N), Grad(N/2+1:N,1:N/2))
    U2(1+N/2:N,N/2+1:N) = U1(1+N/2:N,1+N/2:N) - step *                         &
    &                     matmul(V1(N/2+1:N, 1:N/2), Grad(1:N/2,N/2+1:N))
    V2(1:N/2  ,N/2+1:N)   = V1(1:N/2  ,N/2+1:N) - step *                       &
    &                      matmul(U1(1:N/2, 1:N/2), Grad(1:N/2,N/2+1:N))
    V2(N/2+1:N  ,1:N/2)   = V1(N/2+1:N  ,1:N/2) - step *                       &
    &                      matmul(U1(N/2+1:N, N/2+1:N), Grad(N/2+1:N,1:N/2))

    !Don't forget to orthonormalise
    call ortho(U2,V2)

  end subroutine GradUpdate_sig

function H20_nosig(Ulim,Vlim,hlim,Dlim) result(H20)
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

  function H20_sig(Ulim,Vlim,hlim,Dlim) result(H20)
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
    real(KIND=dp)              :: time(3)
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
    ! The commented out array calculations CAN NOT BE OBTAINED FROM THEIR
    ! symmetric counterparts. HOWEVER, H20 itself is antisymmetric and it is
    ! this symmetry that allows us to only calculate half of each matrix.
    hU(1:N/2,1:N/2)      = matmul(hlim(1:N/2,1:N/2), Ulim(1:N/2,1:N/2))
    !hU(N/2+1:N,N/2+1:N)  = matmul(hlim(N/2+1:N,N/2+1:N), Ulim(N/2+1:N,N/2+1:N))

    !hV(1:N/2, N/2+1:N)   = matmul(hlim(1:N/2,1:N/2), Vlim(1:N/2, N/2+1:N))
    hV(N/2+1:N,1:N/2)    = matmul(hlim(N/2+1:N,N/2+1:N), Vlim(N/2+1:N,1:N/2))
    !
    DsV(1:N/2,1:N/2)     = matmul(Dlim(1:N/2,N/2+1:N), Vlim(N/2+1:N,1:N/2))
    !DsV(N/2+1:N,N/2+1:N) = matmul(Dlim(N/2+1:N,1:N/2), Vlim(1:N/2,N/2+1:N))
    !
    !DsU(1:N/2,N/2+1:N)   = matmul(Dlim(1:N/2,N/2+1:N), Ulim(N/2+1:N,N/2+1:N))
    DsU(N/2+1:N,1:N/2)   = matmul(Dlim(N/2+1:N,1:N/2), Ulim(1:N/2,1:N/2))

    !---------------------------------------------------------------------------
    ! Construct H20
    H20(N/2+1:N, 1:N/2) =                                                      &
    &                matmul(transpose(Ulim(N/2+1:N,N/2+1:N)),hV(N/2+1:N, 1:N/2))
    H20(N/2+1:N, 1:N/2) = H20(N/2+1:N, 1:N/2) -                                &
    &                matmul(transpose(Vlim(1:N/2,N/2+1:N)),  DsV(1:N/2,1:N/2))
    H20(N/2+1:N, 1:N/2) = H20(N/2+1:N, 1:N/2) +                                &
    &               matmul(transpose(Ulim(N/2+1:N,N/2+1:N)),DsU(N/2+1:N, 1:N/2))
    H20(N/2+1:N, 1:N/2) = H20(N/2+1:N, 1:N/2) -                                &
    &                   matmul(transpose(Vlim(1:N/2,N/2+1:N)), hU(1:N/2, 1:N/2))
    do j=1,N/2
      do i=N/2+1,N
        H20(j,i) = - H20(i,j)
      enddo
    enddo
    !---------------------------------------------------------------------------
  end function H20_sig

  function N20(Ulim,Vlim)
    !---------------------------------------------------------------------------
    ! Construct the gradient of the energy with respect to Lambda, N_20.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: Ulim(:,:), Vlim(:,:)
    real(KIND=dp),allocatable  :: N20(:,:)
    integer :: i,j,N,k

    if(.not.allocated(N20)) then
      N = size(Ulim,1)
      allocate(N20(N,N)) ; N20 = 0.0d0
    endif
    N = size(Ulim,1)
    N20 = 0.0d0
    if(SC) then

      N20(1:N/2, N/2+1:N) =                                                    &
      &            matmul(transpose(Ulim(1:N/2, 1:N/2)),Vlim(1:N/2,N/2+1:N))   &
      &          - matmul(transpose(Vlim(N/2+1:N,1:N/2)),Ulim(N/2+1:N,N/2+1:N))

      do j=N/2+1,N
        do i=1,N/2
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

  ! function LN20_sig(Ulim,Vlim, Rho, Kappa) result(LN20)
  !   real(KIND=dp), intent(in)  :: Ulim(:,:), Vlim(:,:), Rho(:,:), Kappa(:,:)
  !   real(KIND=dp)              :: LN20(size(Ulim,1),size(Ulim,1))
  !
  !   integer :: i,j,k,N
  !
  !   !---------------------------------------------------------------------------
  !   ! Auxiliary arrays
  !   real(KIND=dp):: kappaV(size(Ulim,1),size(Ulim,1))
  !   real(KIND=dp):: kappaU(size(Ulim,1),size(Ulim,1))
  !   real(KIND=dp):: rhoV(size(Ulim,1) ,size(Ulim,1))
  !   real(KIND=dp):: rhoU(size(Ulim,1) ,size(Ulim,1))
  !
  !   N = size(Ulim,1)
  !
  !   rhoV(1:N,1:N) = matmul(rho, Vlim)
  !   rhoU(1:N,1:N) = matmul(rho, Ulim)
  !
  !   LN20 = -4.0d0 * matmul(transpose(Ulim), rhoV)                              &
  !   &    +  4.0d0 * matmul(transpose(Vlim), rhoU)
  !
  !   if(all(LN20.eq.0.0))print *, 'LN20 zero'
  ! end function LN20_sig

  subroutine ortho_nosig(U,V)
        !-----------------------------------------------------------------------
        ! Orthogonalises U and V vectors when signature is broken.
        ! Uses a simple Gramm-Schmidt routine.
        !-----------------------------------------------------------------------
        real(KIND=dp), intent(inout) :: U(:,:), V(:,:)
        real(KIND=dp)                :: overlap

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

  subroutine ortho_sig(U,V)
    !-----------------------------------------------------------------------
    ! Orthogonalises U and V vectors when signature is conserved.
    ! Uses a simple Gramm-Schmidt routine.
    !-----------------------------------------------------------------------
      real(KIND=dp), intent(inout) :: U(:,:), V(:,:)
      real(KIND=dp)                :: overlap

      integer :: i,j,k, N

      N = size(U,1)
      do j=1,N/2
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

      do j=N/2+ 1 ,N
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

end module GradientHFB