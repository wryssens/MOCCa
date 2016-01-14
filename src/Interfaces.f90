module Interfaces
!---------------------------------------------------------------------
! Module to contain all special output routines for interfacing MOCCa
! with other codes.
!
! Currently:
!  - Routine for compatibility with proj8.
!
!----------------------------------------------------------------------

use Pairing
use PairingInteraction
use BCS
use Densities

implicit none

contains 
    subroutine writeProm4(OutputFilename)
        !-------------------------------------------------------------
        ! Writes a file to OUTPUTFILENAME.proj8
        !
        !
        !
        !--------------------------------------------------------------

        character(len=*), intent(in) :: OutputFilename

        integer      :: iunit, iver,npn,npp,i,iso,nwaven,nwavep, npair, par
        integer      :: npar(2,2)
        real(KIND=dp):: Eqp(nwt)

        call get_unit(iunit)

        open(iunit, form='unformatted',file=OUTPUTFILENAME)

        !---------------------------------------------------------------
        ! Make prom4 believe it is reading EV8 wavefunctions with iver=3
        iver = 3
        write(iunit) iver
        write(iunit) 'MOCCa Prom4 Output'

        !Convert neutrons and protons to integers
        npp = floor(Protons) ; npn = floor(Neutrons)

        ! Count the number of proton & neutron wavefunctions in all parity blocks
        npar = 0
        do i=1,nwt
            iso = (HFBasis(i)%GetIsospin() + 3)/2
            par = (HFBasis(i)%GetParity()  + 3)/2
            npar(par,iso) = npar(par,iso)  + 1
        enddo
        if(sum(npar) .ne. nwt) then
            print *, npar
            call stp('Wrongly counted the wavefunctions in writeProm4.')
        endif

        write(iunit) sum(npar(:,1)), sum(npar(:,2)), nwavep, npn, npp
        write(iunit) nx,ny,nz

        ! I don't care about the iteration count and just write 0
        write(iunit) 0, dx
        ! Neither do I care about the various constraint parameters, as prom4
        ! doesnt use them anyway
        write(iunit) 0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp
        write(iunit) 0,0,0,0,0,0,0,0
        write(iunit) 0,0,0,0.0_dp,(0.0_dp, i=1,28)
        
        write(iunit) 0,0,(0.0_dp, i=1,28)
        write(iunit) 0, (0.0_dp, i=1,28)
        write(iunit) 0, (0.0_dp, i=1,28)
        ! Neither do I care about transferring the force parameters
        write(iunit) 0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp

        ! I do however care about the implementation of pairing
        if(PairingType.eq. 1) then
            if(Lipkin) then
                npair = 5
            else
                npair = 4
            endif
        elseif(PairingType.eq.2) then
            call stp("Don't try to write HFB wavefunctions to a proj8 file.")
        else
            npair = 1
        endif
        !No qp excitations for the moment so ntqp=0
        write(iunit) npair, Pairingstrength, PairingCut, PairingMu(2), alpha, Fermi, LNLambda, 0

        ! PROM4 does not care about our Deltas & Eqps
        write(iunit)  (0.0_dp , i = 1,nwt)
        write(iunit)  (0.0_dp , i = 1,nwt)

        !Write parities
        write(iunit) (HFBasis(i)%GetEnergy(),i=1,nwt)
        if(TRC) then
            write(iunit) (HFBasis(i)%GetOcc()/2.0_dp,i=1,nwt)
        else
            write(iunit) (HFBasis(i)%GetOcc(),i=1,nwt)
        endif
        write(iunit) (0.0_dp, i=1,nwt)

        ! Write the densities anyway, since it is easy. PROM4 does not care however.
        write(iunit) Density%Rho(:,:,:,1)
        write(iunit) Density%Rho(:,:,:,2)
        write(iunit) Density%Tau
        write(iunit) Density%NablaJ

        do i=1,nwt
            write(iunit) HFBasis(i)%Value%Grid(:,:,:,1,1)
            write(iunit) HFBasis(i)%Value%Grid(:,:,:,2,1)
            write(iunit) HFBasis(i)%Value%Grid(:,:,:,3,1)
            write(iunit) HFBasis(i)%Value%Grid(:,:,:,4,1)
        enddo

        close(iunit)

    end subroutine writeprom4



end module Interfaces