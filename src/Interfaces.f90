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
    subroutine writePromesse(OutputFilename)
        !--------------------------------------------------------------
        ! Writes a file to OUTPUTFILENAME that is suited for reading 
        ! with prom8 or prom4.
        !--------------------------------------------------------------
        ! Note that this routine decides on its own whether output
        ! is destined for prom4 or prom8, based on the conservation
        ! of parity.
        !--------------------------------------------------------------

        character(len=*), intent(in) :: OutputFilename

        integer      :: iunit, iver,npn,npp,i,iso,nwaven,nwavep, npair, par, kparz(nwt)
        integer      :: npar(2,2)
        real(KIND=dp):: Eqp(nwt)

        call get_unit(iunit)

        open(iunit, form='unformatted',file=OUTPUTFILENAME)

        !---------------------------------------------------------------
        ! Make prom4 believe it is reading EV8 wavefunctions with iver=3
        iver = 3
        
        write(iunit) iver
        ! Has to be 20 characters long
        write(iunit) 'MOCCa Prom Output   '

        !Convert neutrons and protons to integers
        npp = floor(Protons) ; npn = floor(Neutrons)

        ! Count the number of proton & neutron wavefunctions in all parity blocks
        npar = 0
        do i=1,nwt
            iso = (HFBasis(i)%GetIsospin() + 3)/2
            ! Prom8 takes the parities in reverse order
            par =  HFBasis(i)%GetParity()
            if(par .eq. -1) then
                par = 2
            elseif(par.eq.0) then
                par = 1
            endif
            npar(par,iso) = npar(par,iso)  + 1
        enddo

        write(iunit) sum(npar(:,1)), sum(npar(:,2)), npn, npp   
        write(iunit) nx,ny,nz

        ! I don't care about the iteration count and just write 0
        write(iunit) 0, dx
        !-------------------------------------------------------
        ! I do however care about the implementation of pairing
        if(PairingType.eq. 1) then
            if(Lipkin) then
                npair = 5
            else
                npair = 4
            endif
        elseif(PairingType.eq.2) then
            print *,"--------------------------------------------------------------" 
            print *,"Don't try to write HFB wavefunctions to a promesse file (yet)."
            print *,'--------------------------------------------------------------'
            close (iunit)
            return
        else
            npair = 1
        endif


        if(PC) then
            ! Prom8 output
            ! imtd line
            write(iunit) 0, (0.0_dp, i=1,28)
            ! qxnc line
            write(iunit) (0.0_dp, i=1,28)
            ! qxpc line
            write(iunit) (0.0_dp, i=1,28)
            ! qxtc line
            write(iunit) (0.0_dp, i=1,28)
            !Qxxn line
            write(iunit) (0.0_dp, i=1,28)
            ! Force parameters
            write(iunit) (0.0_dp, i=1,28)

            write(iunit) npair, Pairingstrength, PairingCut, PairingMu(2), 1.0_dp, alpha, Fermi, LNLambda, 0
        else
            ! Prom4 output
            ! cqcm line
            write(iunit) (0.0_dp, i=1,28)
            ! irtd line
            write(iunit) 0,0,0,0,0,0,0,0
            ! iqd line
            write(iunit) 0,0,0,0.0_dp,(0.0_dp, i=1,28)
            !icqx line
            write(iunit) 0,0,(0.0_dp, i=1,28)
            ! iq30 line
            write(iunit) 0,0,(0.0_dp, i=1,28)
            ! iq4 line
            write(iunit) 0, (0.0_dp, i=1,28)
            ! force parameters line
            write(iunit) (0.0_dp, i=1,28)

            write(iunit) npair, Pairingstrength, PairingCut, PairingMu(2), alpha,Fermi,LnLambda, 0
        endif

        ! PROM4 does not care about our Deltas & Eqps
        write(iunit)  (0.0_dp , i = 1,nwt)
        write(iunit)  (0.0_dp , i = 1,nwt)

        !Write parities
        if(PC) then
            write(iunit) (HFBasis(i)%GetParity(), i=1,nwt)
        endif

        ! Single particle energies
        write(iunit) (HFBasis(i)%GetEnergy(),i=1,nwt)
        !Occupation numbers
        if(TRC) then
            write(iunit) (HFBasis(i)%GetOcc()/2.0_dp,i=1,nwt)
        else
            write(iunit) (HFBasis(i)%GetOcc(),i=1,nwt)
        endif
        !Not yet sure what this is supposed to be.
        write(iunit) (0.0_dp, i=1,nwt)

        if(PC) then
            write(iunit) npar
        endif
 
        ! Write the densities anyway, since it is easy. PROM4 does not care however.
        write(iunit) Density%Rho(:,:,:,1)
        write(iunit) Density%Rho(:,:,:,2)
        write(iunit) Density%Tau
        write(iunit) Density%NablaJ

        ! Prom8 & Prom4 expect wavefunctions normalised to 2!
        do i=1,nwt
            write(iunit) HFBasis(i)%Value%Grid(:,:,:,1,1)*sqrt(2.0_dp)
            write(iunit) HFBasis(i)%Value%Grid(:,:,:,2,1)*sqrt(2.0_dp)
            write(iunit) HFBasis(i)%Value%Grid(:,:,:,3,1)*sqrt(2.0_dp)
            write(iunit) HFBasis(i)%Value%Grid(:,:,:,4,1)*sqrt(2.0_dp)
        enddo

        close(iunit)

    end subroutine writepromesse



end module Interfaces