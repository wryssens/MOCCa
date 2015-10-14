module LipkinNogami
    !----------------------------------------------------
    ! Module that treats Lipkin-Nogami calculations.
    ! Everything is based on
    !
    ! M. Bender, K. Rutz, P.-G. Reinhard and J.A. Maruhn
    ! Pairing Gaps from Nuclear Mean-Field models
    ! The European Physical Journal A
    ! July 2000, Volume 8, Issue 1, pp 59-75
    !----------------------------------------------------
    use Geninfo
    use Densities
    use Force
    use Spwfstorage

    implicit none   

    !---------------------------------------------------
    ! Response densities
    type(Densityvector) :: Response

contains

    function CalcLambda2() result(L2)
        !-----------------------------------------------
        ! Function returns the calculated Lambda2.
        !
        ! Note that this is a different calculation 
        ! from the old-school one. 
        ! In older versions of CR8, the Lambda_2 parameter
        ! was calculated with the Delta & Delta_LN and 
        ! fermi values currently in memory.
        !
        ! However, this routine assumes that a correct
        ! Fermi value has been found, the pairing problem
        ! has been solved and we can thus now calculate
        ! a Lambda_2 value for the next MEAN-FIELD
        ! iteration.
        !------------------------------------------------

        real(KIND=dp) :: L2(2), Denominator(2), Numerator(2)
        real(KIND=dp) :: u2v2(2), u4v4(2), v2
        integer       :: it, i

        !Compute the response densities
        call ComputeDensity(Response,.false.,.true.)

        ! Compute the Denominator
        Denominator = 0.0_dp
        do i=1,nwt
            it = (DensityBasis(i)%GetIsospin()+3)/2
            v2 =  DensityBasis(i)%GetOcc() 
            u2v2(it) = u2v2(it) + v2 * (1 - v2)
            u4v4(it) = u4v4(it) + (v2 * (1 - v2))**2
        enddo 
        Denominator = 8 * u2v2**2 - 16 * u4v4

        Numerator = SkyrmeResponse()! + PairingResponse() 

        L2 = Numerator/Denominator
    end function CalcLambda2

    function SkyrmeResponse() 
        !-----------------------------------------------------------------
        ! Compute the Skyrme linear response.
        !
        !
        !-----------------------------------------------------------------

        real(KIND=dp) :: RhoTot(nx,ny,nz), RhoSum(nx,ny,nz), Rho(nx,ny,nz,2)
        real(KIND=dp) :: STot(nx,ny,nz,3), SSum(nx,ny,nz), S(nx,ny,nz,3,2)
        real(KIND=dp) :: B7coef, B8coef, B12Coef, B13Coef
        real(KIND=dp) :: B8Terms(3), B12Terms(3), B13Terms(3), B16Terms(3)
        real(KIND=dp):: SkyrmeResponse(2)
        integer       :: it

        ! Total NORMAL density for the density dependent terms
        RhoTOT = Density%Rho(:,:,:,1)    + Density%Rho(:,:,:,2) 
        RhoSum = Density%Rho(:,:,:,1)**2 + Density%Rho(:,:,:,2)**2
        Rho    = Density%Rho
        if(.not.TRC) then
            S    = Density%vecs
            STot = Density%VecS(:,:,:,:,1)    + Density%VecS(:,:,:,:,2)
            SSum = sum(Density%VecS(:,:,:,:,1)**2 + Density%VecS(:,:,:,:,2)**2,4)
        endif 
        do it=1,2
            !B1 & B2 term
            SkyrmeResponse(it) = SkyrmeResponse(it) + 2*(B1 + B2) &
            &              *sum(Response%Rho(:,:,:,it)**2)
            !B3 & B4
            SkyrmeResponse(it) = SkyrmeResponse(it) + 2*(B3 + B4) &
            &             *sum(Response%Rho(:,:,:,it)*Response%tau(:,:,:,it))
            !B5 & B6
            SkyrmeResponse(it) = SkyrmeResponse(it) + 2*(B5 + B6) &
            &             *sum(Response%Rho(:,:,:,it)*Response%LapRho(:,:,:,it))
            !B7
            if(B7b .ne. 0.0_dp .or. B8b .ne. 0.0_dp) call stp('LN not ready for B7b and B8b')

            B7coef = (2 + byt3a) * (1 + byt3a)
            SkyrmeResponse(it) = SkyrmeResponse(it) + B7a &
            &             * B7Coef * sum(RhoTot**byt3a*Response%Rho(:,:,:,it)**2)
            !B8
            B8coef = (1 + byt3a) * byt3a
            B8Terms(1) = sum(RhoTot**(byt3a  )                * Response%Rho(:,:,:,it)**2)
            B8Terms(2) = sum(RhoTot**(byt3a-2) * RhoSum       * Response%Rho(:,:,:,it)**2)
            B8Terms(3) = sum(RhoTot**(byt3a-1) * Rho(:,:,:,it)* Response%Rho(:,:,:,it)**2)

            SkyrmeResponse(it) = SkyrmeResponse(it) + B8a*(2*B8Terms(1) + B8coef * B8Terms(2) + 4*byt3a * B8Terms(3))
            !B9
            SkyrmeResponse(it) = SkyrmeResponse(it) + 2*(B9 + B9q) &
            &             * sum(Response%Rho(:,:,:,it) * Response%NablaJ(:,:,:,it))
            if(.not.TRC) then
                SkyrmeResponse(it) = SkyrmeResponse(it) + 2*(B9 + B9q) &
                &         * sum(sum(Response%vecj(:,:,:,:,it)*Response%RotS(:,:,:,:,it),4))       
            endif
            !B10 & B11
            if( B10 .ne. 0.0_dp .or. B11 .ne. 0.0_dp) then 
                SkyrmeResponse(it) = SkyrmeResponse(it) + 2*(B10 + B11)*sum(Response%vecs(:,:,:,:,it)**2)
            endif

            !B12 & B13
            ! To be checked for bugs, as this is quite involved...
            ! Currently only for B12a & B13a terms
            if(B12b .ne. 0.0_dp .or. B13b .ne. 0.0_dp) call stp('LN not ready for B12b and B13b')
            if(B12a.ne.0.0_dp .or. B13a.ne. 0.0_dp) then
                B12Terms(1) = 2*      sum(RhoTot**byt3a     * sum(Response%vecs(:,:,:,:,it)**2,4))
                B12Terms(2) = 4*byt3a*sum(RhoTot**(byt3a-1) * Response%Rho(:,:,:,it)*&
                    &                 sum(STot*Response%vecs(:,:,:,:,it),4))
                B12Terms(3) = B8coef *sum(RhoTot**(byt3a-2) * sum(Stot**2,4) *Response%Rho(:,:,:,it)**2)
                SkyrmeResponse(it) = SkyrmeResponse(it) + B12a*sum(B12Terms)

                B13Terms(1) = 2*      sum(RhoTot**byt3a     * sum(Response%vecs(:,:,:,:,it)**2,4))
                B13Terms(2) = 4*byt3a*sum(RhoTot**(byt3a-1) * Response%Rho(:,:,:,it)*&
                    &                 sum(STot*Response%vecs(:,:,:,:,it),4))
                B13Terms(3) = B8coef *sum(RhoTot**(byt3a-2) * SSum *Response%Rho(:,:,:,it)**2)
                SkyrmeResponse(it) = SkyrmeResponse(it) + B13a*sum(B13Terms)
            endif

            !B14 & B15
            if(B14.ne.0.0_dp .or. B15 .ne. 0.0_dp) then
                SkyrmeResponse(it) = SkyrmeResponse(it) + 2 * (B14 + B15) * &
                &               sum(Response%JMuNu(:,:,:,:,:,it)**2)
                if(allocated(Density%VecT)) then
                    SkyrmeResponse(it) = SkyrmeResponse(it) + 2 * (B14+B15)*&
                    &             sum(Response%VecS(:,:,:,:,it) * Response%VecT(:,:,:,:,it))
                endif
            endif

            !B16 & B17
            if(B16.ne.0.0_dp .or. B17.ne.0.0_dp) then
                !                                sum J_mu,mu**2
                B16Terms(1) = 2 * (B16 + B17) * sum(Response%Jmunu(:,:,:,:,:,it)**2)
                !Extra factor of two for the interchange mu <=> nu
                !                                sum J_{mu,nu} J_{mu,nu}
                B16Terms(2) = 4 * (B16 + B17) * sum(Response%JMuNu(:,:,:,1,2,it)*Response%JMuNu(:,:,:,1,2,it) &
                    &                             + Response%JMuNu(:,:,:,1,3,it)*Response%JMuNu(:,:,:,3,1,it) &
                    &                             + Response%JMuNu(:,:,:,2,3,it)*Response%JMuNu(:,:,:,3,2,it))

                B16Terms(3) =-4 * (B16 + B17) * sum(Response%vecs(:,:,:,:,it) * Response%VecF(:,:,:,:,it))
            
                SkyrmeResponse(it) = SkyrmeResponse(it) + sum(B16Terms)
            endif

            !B18 & B19
            if(B18.ne.0.0_dp .or. B19.ne.0.0_dp) then
                SkyrmeResponse(it) = SkyrmeResponse(it) + 2*(B18 + B19) &
                &             * sum(Response%vecs(:,:,:,:,it) * Response%LapS(:,:,:,:,it))
            endif

            !B20 & B21
            if(B20.ne. 0.0_dp .or. B21 .ne. 0.0_dp) then
                SkyrmeResponse(it) = SkyrmeResponse(it) + 2*(B20+B21) & 
                &             * sum(Response%DerS(:,:,:,:,:,it)**2)
            endif
        enddo
        !Global factor dv
        SkyrmeResponse = SkyrmeResponse*dv
        !----------------------------------------------------------------------------------------------
    end function SkyrmeResponse

!     function PairingResponse()
!         !------------------------------------------------------------
!         ! Response of the pairing functional.
!         !
!         !------------------------------------------------------------

!         use Pairinginteraction

!         integer       :: it
!         real(KIND=dp) :: PairingResponse(2)
!         ! Pairing response densities
!         real(KIND=dp), allocatable, save :: Chi(:,:,:,:), ChiStar(:,:,:,:)


        

!         do it=1,2
!            ! PairingResponse(it) = 
!         enddo

!     end function PairingResponse

end module LipkinNogami