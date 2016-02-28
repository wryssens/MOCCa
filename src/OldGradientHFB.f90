module GradientHFB
    !--------------------------------------------------------------------------------------
    ! Module that should be capable of solving the HFB problem by doing gradient descent.
    !
    !--------------------------------------------------------------------------------------
    !
    !
    !
    !
    !
    !

    use geninfo
    use HFB

    implicit none

    !----------------------------------------------------------------
    ! Local storage for the U and V
    real(KIND=dp),allocatable ::  GradU(:,:,:,:,:), GradV(:,:,:,:,:)
    !----------------------------------------------------------------
    ! Blocked structure for the single-particle hamiltonian and 
    ! delta
    real(KIND=dp),allocatable :: hblock(:,:,:,:,:), dblock(:,:,:,:,:)

contains

    subroutine GetUandV
        !----------------------------------------
        ! Get U and V values from the HFB module.
        !
        !----------------------------------------

        integer :: it,P,i,j,N, ind(2,2), Rzindex,k
        logical :: incolumns

        !------------------------
        ! Initialize the matrices. 
        !-------------------------
        N = maxval(blocksizes)
        Rzindex= 1 ; if(SC) Rzindex = 2

        allocate(GradU(N,N,Rzindex,Pindex,Iindex)) ; GradU = 0.0_dp
        allocate(GradV(N,N,Rzindex,Pindex,Iindex)) ; GradV = 0.0_dp

        do it=1,Iindex
            do P=1,Pindex
                N = blocksizes(P,it)

                ind(P,it) = 1
                do j=1,2*N
                    ! This loop makes sure that we go through the 
                    ! big matrix in order. Otherwise the signature
                    ! block might end up somewhere we don't expect
                    ! them.
                    incolumns = .false.
                    do k=1,N
                        if (j .eq. HFBColumns(k,P,it)) then
                            incolumns=.true.
                            exit
                        endif
                    enddo
                    if(incolumns) then
                        if(SC) then
                            if(j.le.N) then
                                GradV(1:N/2,ind(P,it),1,P,it) =  DBLE(V(  N/2+1:N  ,j,P,it))
                                GradU(1:N/2,ind(P,it),1,P,it) =  DBLE(U(      1:N/2,j,P,it))
                            else
                                GradV(1:N/2,ind(P,it)-N/2,2,P,it) =  DBLE(V(      1:N/2,j,P,it))
                                GradU(1:N/2,ind(P,it)-N/2,2,P,it) =  DBLE(U(  N/2+1:N  ,j,P,it))
                            endif
                        else
                            GradV(1:N,ind(P,it),1,P,it) = DBLE(V(1:N,j,P,it))
                            GradU(1:N,ind(P,it),1,P,it) = DBLE(U(1:N,j,P,it))
                        endif
                        ind(P,it) = ind(P,it)+ 1
                    endif
                enddo
                if(ind(P,it).ne.N+1) call stp('errorcheckind')
            enddo
        enddo

    end subroutine GetUandV

    subroutine SetUandV
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
                        !V(i,j+N,P,it) = GradV(i,j,P,it)
                        !U(i,j+N,P,it) = GradU(i,j,P,it)
                    enddo
                enddo
            enddo
        enddo

    end subroutine SetUandV

    subroutine BlockHFBHamil(Delta)
        !-------------------------------------------------
        ! Load the blockstructure of the HFBhamiltonian.
        ! but not including the Lambda and Lambda_2
        ! contribution!
        
        complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
        real(KIND=dp)                             :: z(2) = 0.0_dp
        integer                                   :: N,i,j,P,it, Rzindex, Rz

        call constructHFBHamiltonian(z,Delta,z,0.0_dp)
        N = maxval(blocksizes)

        Rzindex = 1 ; if(SC) Rzindex = 2

        if(.not.allocated(hblock)) then
            allocate(hblock(N/Rzindex,N/Rzindex,Rzindex,Pindex,Iindex))
            allocate(dblock(N/Rzindex,N/Rzindex,Rzindex,Pindex,Iindex))
        endif

        do it=1,Iindex
            do P=1,Pindex
                N = blocksizes(P,it)
                ! positive signature
                hblock(1:N/2,1:N/2,1,P,it) = HFBHamil(      1:N/2,        1:N/2,P,it)
                dblock(1:N/2,1:N/2,1,P,it) = HFBHamil(      1:N/2  ,3*N/2+1:N  ,P,it)

                ! negative signature
                hblock(1:N/2,1:N/2,2,P,it) = HFBHamil(N/2+1:N,N/2+1:N    ,P,it)
                dblock(1:N/2,1:N/2,2,P,it) = HFBHamil(N/2+1:N,  N+1:3*N/2,P,it)
            enddo
        enddo

    end subroutine BlockHFBHamil


    subroutine HFBFermiGradient(Fermi,L2,Delta,DeltaLN,Lipkin,Prec)

        real(KIND=dp), intent(inout)              :: Fermi(2), L2(2)
        complex(KIND=dp), allocatable, intent(in) :: Delta(:,:,:,:)
        complex(KIND=dp), allocatable, intent(in) :: DeltaLN(:,:,:,:)
        logical, intent(in)                       :: Lipkin
        real(KIND=dp), intent(in)                 :: Prec

        real(KIND=dp),allocatable                 :: Gradient(:,:,:,:,:)
        real(KIND=dp),allocatable                 :: oldU(:,:,:,:,:),oldV(:,:,:,:,:)
        real(KIND=dp)                             :: gradientnorm


        integer :: i,j,P,it,N, Rzindex,iter

        Rzindex = 1 ; if(SC) Rzindex = 2

        if(.not. allocated(GradU)) then
            N = maxval(blocksizes)
            ! Get U and V from the HFB module
            call GetUandV()
            OldU = GradU
            OldV = GradV
            allocate(Gradient(N/Rzindex,N/Rzindex,Rzindex,Pindex,Iindex)) ; Gradient = 0.0_dp
        endif

        ! Construct the matrices in block form.
        call BlockHFBHamil(Delta)

        do iter=1,HFBIter
            !--------------------------------------------------------------------------------------------------
            ! Construct the gradients
            do it=1,Iindex
                do P=1,Pindex
                    N = blocksizes(P,it)
                    Gradient(1:N/2,1:N/2,1,P,it) = H20(GradU(1:N/2,1:N/2,1,P,it), GradV(1:N/2,1:N/2,1,P,it),&
                        &                             hblock(1:N/2,1:N/2,1,P,it),Dblock(1:N/2,1:N/2,1,P,it))
                    do i=1,N/2
                        print *, Gradient(i,1:N/2,1,P,it)
                    enddo
                    print *
                     Gradient(1:N/2,1:N/2,2,P,it) = H20(GradU(1:N/2,1:N/2,2,P,it), GradV(1:N/2,1:N/2,2,P,it),&
                        &                            hblock(1:N/2,1:N/2,2,P,it),Dblock(1:N/2,1:N/2,2,P,it))
                    do i=1,N/2
                        print *, Gradient(i,1:N/2,2,P,it)
                    enddo
                    print *
                enddo
            enddo
            !---------------------------------------------------------------------------------------------------
            stop
            do it=1,2
                do P=1,Pindex
                    N = blocksizes(P,it)
                    oldU(1:N,1:N,1:2,P,it) = GradU(1:N,1:N,1:2,P,it)
                    oldV(1:N,1:N,1:2,P,it) = GradV(1:N,1:N,1:2,P,it)
                    call GradUpdate(0.02_dp, oldU(1:N/2,1:N/2,1,P,it),oldV(1:N/2,1:N/2,1,P,it), &
                        &                    Gradient(1:N/2,1:N/2,1,P,it),                      &
                        &                    GradU(1:N/2,1:N/2,1,P,it),GradV(1:N/2,1:N/2,1,P,it))
                    call GradUpdate(0.02_dp, oldU(1:N/2,1:N/2,2,P,it),oldV(1:N/2,1:N/2,2,P,it), &
                        &                    Gradient(1:N/2,1:N/2,2,P,it),                      &
                        &                    GradU(1:N/2,1:N/2,2,P,it),GradV(1:N/2,1:N/2,2,P,it))
                enddo
            enddo

            gradientnorm = sum(Gradient**2)
            print *, gradientnorm

        enddo
        stop

    end subroutine HFBFermiGradient

    subroutine GradUpdate(step, U1,V1,Grad,U2,V2)

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

    end subroutine GradUpdate

    subroutine ortho(U,V)    
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

    end subroutine ortho

    function H20(Ulim,Vlim,hlim, Dlim)
        !-----------------------------------------------------
        ! Computes H20 from limited U and V, using only
        ! limited hlim and Dlim.
        !-----------------------------------------------------
        real(KIND=dp),intent(in) :: Ulim(:,:), hlim(:,:)
        real(KIND=dp),intent(in) :: Vlim(:,:), Dlim(:,:)

        !-----------------------------------------------------
        ! Auxiliary arrays
        real(KIND=dp)    :: DsV(2*nwt,2*nwt), DsU(2*nwt,2*nwt)
        real(KIND=dp)    ::  hV(2*nwt,2*nwt),  hU(2*nwt,2*nwt)
        
        !-----------------------------------------------------
        ! result
        real(KIND=dp),allocatable :: H20(:,:)


        integer :: i,j,N,k

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
        allocate(H20(N,N))
        do j=1,N
            do i=1,N
                H20(i,j) = 0.0_dp
                do k=1,N
                    H20(i,j) = H20(i,j)  +  Ulim(k,i) * hV (k,j) &
                    &                    -  Vlim(k,i) * hU (k,j) &
                    &                    -  Vlim(k,i) * DsV(k,j) &
                    &                    +  Ulim(k,i) * DsU(k,j)
                enddo
                ! H20 is antisymmetric
                !H20(j,i) = -H20(i,j)
            enddo
        enddo
    end function H20


    function N20(Ulim,Vlim)
        real(KIND=dp),intent(in) :: Ulim(:,:), Vlim(:,:)

        real(KIND=dp),allocatable :: N20(:,:)
    end function N20


end module