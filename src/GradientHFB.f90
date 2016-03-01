module GradientHFB

  use HFB

  implicit none

  !----------------------------------------------------------------
  ! Local storage for the U and V
  real(KIND=dp),allocatable  :: GradU(:,:,:,:), GradV(:,:,:,:)

  !----------------------------------------------------------------
  ! Store single-particle hamiltonian and delta in block-form
  real(KIND=dp), allocatable :: hblock(:,:,:,:), dblock(:,:,:,:)

  real(KIND=dp),allocatable  :: Gradient(:,:,:,:),OldGrad(:,:,:,:)
  real(KIND=dp),allocatable  :: aN20(:,:,:,:)
  real(KIND=dp),allocatable  :: oldgradU(:,:,:,:),oldgradV(:,:,:,:)
  real(KIND=dp),allocatable  :: Direction(:,:,:,:), OldDir(:,:,:,:)

  procedure(H20_sig), pointer            :: H20
  procedure(constructrho_sig), pointer   :: constructrho
  procedure(constructkappa_sig), pointer :: constructkappa
  procedure(gradupdate_sig), pointer     :: gradupdate
  procedure(ortho_sig), pointer          :: ortho
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
            if(ind(P,it).ne.N+1) then
              call stp('errorcheckind', 'N', n, 'ind', ind(P,it))
            endif
        enddo
    enddo

  end subroutine GetUandV

  subroutine OutHFBModule
      !-------------------------------------------------
      ! Put our solution into the HFB module correctly.
      !-------------------------------------------------
      ! TODO: put the conjugate states too
      !-------------------------------------------------
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

  subroutine BlockHFBHamil(Delta)
      !-------------------------------------------------
      ! Load the blockstructure of the HFBhamiltonian.
      !-------------------------------------------------

      complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
      real(KIND=dp)                             :: z(2) = 0.0_dp
      integer                                   :: N,i,j,P,it, Rzindex, Rz

      call constructHFBHamiltonian(z,Delta,z,0.0_dp)
      N = maxval(blocksizes)

      if(.not.allocated(hblock)) then
          allocate(hblock(N,N,Pindex,Iindex))
          allocate(dblock(N,N,Pindex,Iindex))
      endif

      do it=1,Iindex
          do P=1,Pindex
              N = blocksizes(P,it)
              ! positive signature
              hblock(1:N,1:N,P,it) = HFBHamil(1:N,1:N     ,P,it)
              dblock(1:N,1:N,P,it) = HFBHamil(1:N,N+1:2*N ,P,it)
          enddo
      enddo

  end subroutine BlockHFBHamil

  subroutine HFBFermiGradient(Fermi,L2,Delta,DeltaLN,Lipkin,Prec)

    implicit none

    real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
    complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
    complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
    logical, intent(in)                       :: Lipkin
    real(KIND=dp), intent(in)                 :: Prec

    real(KIND=dp)              :: gamma(2,2), maxstep, step
    real(KIND=dp)              :: N20norm(2), par(2), slope,  OldFermi(2)
    real(KIND=dp)              :: gradientnorm(2,2), oldnorm(2,2), PR(2,2)

    integer :: i,j,P,it,N, Rzindex,iter, inneriter, first=1
    logical :: converged(2)

    Particles(1) = neutrons
    Particles(2) = protons

    if(.not. allocated(GradU)) then
        N = maxval(blocksizes)
        ! Get U and V from the HFB module
        call GetUandV()
        OldgradU = GradU
        OldgradV = GradV
        allocate(Direction(N,N,pindex,Iindex)); Direction= 0.0_dp
        allocate(Gradient(N,N,Pindex,Iindex)) ; Gradient = 0.0_dp
        allocate(OldGrad(N,N,Pindex,Iindex))  ; Gradient = 0.0_dp
        allocate(aN20(N,N,Pindex,Iindex))     ; aN20     = 0.0_dp
        allocate(OldDir(N,N,Pindex,Iindex))   ; OldDir   = 0.0_dp

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

    ! Construct the matrices in block form.
    call BlockHFBHamil(Delta)

    oldnorm = 0.0d0 ; gradientnorm = 0.0d0
    step = 0.020_dp ; OldFermi = 0.0d0
    converged = .false.
    do iter=1,HFBIter
      do it=1,Iindex
          if(converged(it)) cycle
          n20norm(it) = 0.0d0
          do P=1,Pindex
              N = blocksizes(P,it)

              if(any(GradU(:,:,P,it) .eq. GradU(:,:,P,it) + 1 ))call stp('inigrad')

              Oldgrad(1:N,1:N,p,it)  = Gradient(1:N,1:N,P,it)
              Gradient(1:N,1:N,P,it) = H20(GradU(1:N,1:N,P,it),                &
              &                            GradV(1:N,1:N,P,it),                &
              &                            hblock(1:N,1:N,P,it),               &
              &                            Dblock(1:N,1:N,P,it))

              aN20(1:N,1:N,P,it) = N20(GradU(1:N,1:N,P,it),GradV(1:N,1:N,P,it))
              do j=1,N
                do i=1,N
                  N20norm(it) = N20norm(it) + aN20(i,j,P,it) * aN20(i,j,P,it)
                enddo
              enddo
          enddo

          par(it) = 0.0
          do P=1,Pindex
            N = blocksizes(P,it)
            do i=1,N
              do j=1,N
                par(it) = par(it) + GradV(i,j,P,it) * GradV(i,j,P,it)
              enddo
            enddo
          enddo

          Fermi(it) = Fermi(it) + (Particles(it) - par(it))/N20norm(it)

          do P=1,Pindex
              N = blocksizes(P,it)
              Gradient(1:N,1:N,P,it) = Gradient(1:N,1:N,P,it)                  &
              &                                 - Fermi(it) * aN20(1:N,1:N,P,it)

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
          if(all(abs(GradientNorm(:,it)).lt.prec) .and.                        &
          &     abs(Par(it) - Particles(it)).lt.Prec) then
           converged(it)=.true.
         endif
      enddo

      do it=1,Iindex
          if(converged(it)) cycle
          do P=1,Pindex
              N = blocksizes(P,it)
              Direction(1:N,1:N,P,it) = Gradient(1:N,1:N,P,it)

              if(oldnorm(P,it).ne.0.0d0) then
                gamma(P,it) = gradientnorm(p,it)/oldnorm(p,it)
                gamma(P,it) = gamma(P,it) - PR(p,it)/oldnorm(p,it)
                !gamma(P,it) = max(gamma(P,it),0.0_dp)
                Direction(1:N,1:N,P,it) = Direction(1:N,1:N,P,it)  +           &
                &                             gamma(P,it) * OldDir(1:N,1:N,P,it)
              endif

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

              oldgradU(1:N,1:N,P,it) = GradU(1:N,1:N,P,it)
              oldgradV(1:N,1:N,P,it) = GradV(1:N,1:N,P,it)

              call Linesearch(oldgradU(1:N,1:N,P,it), oldgradV(1:N,1:N,P,it),  &
              &               Gradient(1:N,1:N,P,it), Direction(1:N,1:N,P,it), &
              &               step,   gradU(1:N,1:N,P,it), gradV(1:N,1:N,P,it),&
              &             hblock(1:N,1:N,P,it),Dblock(1:N,1:N,P,it),Fermi(it))

              OldDir(1:N,1:N,P,it)  = Direction(1:N,1:N,P,it)
              oldnorm(P,it)         = gradientnorm(P,it)
          enddo
      enddo
      print *, iter, gradientnorm
      if(all(converged)) then
        print *, 'Converged', iter
        exit
      endif
   enddo
   call OutHFBModule
 end subroutine HFBFermiGradient

  subroutine Linesearch(OldU, OldV, Grad, Direction, maxstep, Ulim, Vlim, hlim, Dlim, Fermi)
    !-------------------------------------------------
    ! This does not really speed up our algorithm...
    !
    !
    !-------------------------------------------------
    real(KIND=dp), intent(in) :: OldU(:,:), OldV(:,:), hlim(:,:), Dlim(:,:)
    real(KIND=dp), intent(in) :: Grad(:,:), Fermi
    real(KIND=dp), intent(out):: Ulim(:,:), Vlim(:,:)
    real(KIND=dp), intent(inout) :: maxstep,Direction(:,:)

    real(KIND=dp), allocatable :: rho(:,:), kappa(:,:), E, step, Estart, slope
    integer :: i,j,k,N, iter, reset

    if(.not.allocated(rho)) then
      N = maxval(blocksizes)
      allocate(Rho(N,N))  ; Rho = 0.0d0
      allocate(Kappa(N,N)); Kappa=0.0d0
    endif

    N = size(OldU,1)

    call GradUpdate(maxstep, oldU,oldV,Direction,Ulim,Vlim)
    !
    ! !call stp('Search direction is not a descent direction.', 'slope', slope)
    !
    ! rho(1:N,1:N)  = ConstructRho(OldV)
    ! kappa(1:N,1:N)= ConstructKappa(OldU,OldV)
    !
    ! Estart = Energy(hlim,Dlim,rho,kappa, Fermi)
    !
    ! maxstep = maxstep * 100
    ! do iter=1,10
    !
    !   call GradUpdate(maxstep, oldU,oldV,Direction,Ulim,Vlim)
    !   rho(1:N,1:N)  = ConstructRho(Vlim)
    !   kappa(1:N,1:N)= ConstructKappa(Ulim,Vlim)
    !
    !   E = Energy(hlim,Dlim,rho,kappa, Fermi)
    !   print *, maxstep, E
    !   maxstep = maxstep * 0.8
    ! enddo
    ! print *, maxstep
    ! stop
  end subroutine Linesearch

  function ConstructRho_nosig(Vlim) result(Rho)

    real(KIND=dp), intent(in) ::Vlim(:,:)
    real(KIND=dp),allocatable :: Rho(:,:)

    integer :: N,i,j,k

    if(.not.allocated(Rho)) then
      N = maxval(blocksizes)
      allocate(Rho(N,N))
    endif

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

    real(KIND=dp), intent(in) ::Vlim(:,:)
    real(KIND=dp),allocatable :: Rho(:,:)

    integer :: N,i,j,k

    if(.not.allocated(Rho)) then
      N = maxval(blocksizes)
      allocate(Rho(N,N))
    endif

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

  function Energy(hlim,Dlim,rho,kappa, Fermi)
    !--------------------------------------------------
    ! Compute the HFB energy associated with U and V.
    !
    ! E = Tr( h \rho ) - 0.5 * Tr(Delta \kappa) - \lambda <N>
    !
    !--------------------------------------------------

    real(KIND=dp), intent(in) :: hlim(:,:), Dlim(:,:),rho(:,:), kappa(:,:)
    real(KIND=dp), intent(in) :: Fermi
    real(KIND=dp)             :: Energy

    integer                   :: N, i,j

    N = size(hlim,1)
    Energy = 0.0_dp
    do i=1,N
      do j=1,N
        Energy = Energy + hlim(i,j) * Rho(j,i) - 0.5*kappa(i,j)*Dlim(j,i)
      enddo
    Energy = Energy   - Fermi     * Rho(i,i)
    enddo
  end function Energy

  subroutine GradUpdate_nosig(step, U1,V1,Grad,U2,V2)

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
    call ortho(U2,V2)

  end subroutine GradUpdate_nosig

  subroutine GradUpdate_sig(step, U1,V1,Grad,U2,V2)

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
    ! symmetric counterparts.
    ! HOWEVER, H20 itself is antisymmetric and it is this symmetry that allows
    ! us to only calculate half of each matrix.
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
      do j=N/2+1,N
        do i=1,N/2
          N20(i,j) = 0.0d0
          do k=1,N/2
            N20(i,j) = N20(i,j) + Ulim(k,i) * Vlim(k,j)
          enddo
          do k=N/2+1,N
            N20(i,j) = N20(i,j) - Vlim(k,i) * Ulim(k,j)
          enddo

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

  subroutine ortho_nosig(U,V)
        !-------------------------
        ! Orthogonalises vectors.
        !
        !-------------------------
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
        !-------------------------
        ! Orthogonalises vectors.
        !
        !-------------------------
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

          Overlap = 1.0_dp/sqrt(overlap)
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
