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
use Derivatives

implicit none

  !-----------------------------------------------------------------------------
  ! Variables read from the files of the ancient codes, mostly used for 
  ! continuing calculations after some massage of the data.
  !-----------------------------------------------------------------------------
  real(KIND=dp), allocatable :: xkap(:,:,:,:), deltacr8(:,:,:,:)
  real(KIND=dp), allocatable :: rrt(:,:,:,:) , rrn(:,:,:,:)
  real(KIND=dp), allocatable :: dr(:,:,:,:), di(:,:,:,:)
  real(KIND=dp), allocatable :: dlnr(:,:,:,:), dlni(:,:,:,:)
  real(KIND=dp), allocatable :: cr8rho(:,:,:,:)
  real(KIND=dp), allocatable :: cr8vtau(:,:,:,:),cr8vdiv(:,:,:,:)
  integer                    :: npar(4,2)=0, nwtn = 0
  real(KIND=dp), allocatable :: Cr8A(:,:,:,:,:), cr8S(:,:,:,:,:)
  real(KIND=dp), allocatable :: Cr8U(:,:,:,:)


contains 

  subroutine ReadNil(IChan, InputFilename)
  !-----------------------------------------------------------------------------
  ! Taken from ev8.1.2.8 on 12/04/13. 
  ! Quantities that are not useful at the moment of creation are commented out
  ! and not saved.
  !-----------------------------------------------------------------------------
    use CompilationInfo
    use GenInfo
    use WaveFunctions
    use SpwfStorage
    use MeanFields

    implicit none

    integer, intent(in)                 :: IChan
    character(len=*), intent(in)        :: InputFilename
    ! Variables used for checking the input
    integer                             :: iver,nwtn,nwtp,npn,npp,mx,my,mz
    integer                             :: nwt0,itert
    real(KIND=dp)                       :: dv2,dx2
    ! Variables to put onto the Spwf
    integer                             :: kparz(nwt),Isospin
    real(KIND=dp)                       :: esp1(nwt),v2(nwt)
    type(Spinor)                        :: wf
    ! Do loop variables
    integer                             :: nw
    real(KIND=dp), dimension(nx,ny,nz)  :: WF1,WF2,WF3,WF4

    Recalc = .true.
  
    open  (IChan,form='unformatted',file=InputFileName)
    ! Read the version info and other relevant parameters.
    read  (IChan) iver
    read  (IChan) !headd (This variable is useless to remember)
    read  (IChan) nwtn,nwtp,npn,npp
    read  (IChan) mx,my,mz
    nwt0 = nwtn+nwtp

    !Checking the size of the grid and the number of Spwf.
    if   (nwt0.ne.nwt) then
      call stp ("Number of wavefunctions doesn't match file",                   &
      &         "nwt = ", nwt , ' nwt0 = ', nwt0)
    endif
    if    (nx.ne.mx)  call stp (' mx !', 'nx',nx,'mx' ,mx)
    if    (ny.ne.my)  call stp (' my !', 'ny',ny,'my', my)
    if    (nz.ne.mz)  call stp (' mz !', 'nz',nz,'mz', mz)

    !This subroutine is only compatible with versions of NIL8 & EV8 with
    ! iver = 5 or 6
    if (iver.lt.5 .or. iver.gt.6) call stp (' iver !', 'iver', iver)

    if (iver.eq.5) then
        read (IChan) itert,dx2
        dv2=2**(3)*dx2**3
        
        read (IChan) !imtd
        read (IChan) ! qxnc
        read (IChan) ! qxpc    
        read (IChan) ! qxtc
        read (IChan) ! qxxn
        read (IChan) ! t0
        read (IChan) !npair
        read (IChan) !eqp
        read (IChan) !delta
    endif
    
    if (iver.eq.6) then
        read (IChan) itert,dx2
        dv2=2**(3)*dx2**3
        
        read (IChan) !imtd,ral,epscst,cqr,cq2,rcut,acut
        !1            ,delqst,q1nst,q1pst,q1tst,q2nst,q2pst,q2tst
        read (IChan) 
        read (IChan)
        read (IChan)
        read (IChan)
        read (IChan) 
        read (IChan)
        read (IChan) 
        read (IChan) 
        read (IChan) !(eqp(i),i=1,nwt)
        read (IChan) !(delta(i),i=1,nwt)
    endif

    read (IChan) (kparz(nw),nw=1,nwt)
    read (IChan) (esp1  (nw),nw=1,nwt)
    read (IChan) v2
    read (IChan) ! v22
    read (IChan) ! npar
 
    read (IChan) Density%Rho(:,:,:,1)
    read (IChan) Density%Rho(:,:,:,2)
    read (IChan) Density%Tau
    read (IChan) Density%NablaJ

    call ChangeNumberWaveFunctions(nwt)
    
    do nw=1,nwt
      read (IChan) WF1
      read (IChan) WF2
      read (IChan) WF3
      read (IChan) WF4
      !The EV8 wavefunctions are normalised to two.
      WF1=WF1/sqrt(2.0_dp)
      WF2=WF2/sqrt(2.0_dp)
      WF3=WF3/sqrt(2.0_dp)
      WF4=WF4/sqrt(2.0_dp)

      wf = NewSpinor()
      call wf%SetComponent(1,WF1)
      call wf%SetComponent(2,WF2)
      call wf%SetComponent(3,WF3)
      call wf%SetComponent(4,WF4)
      
      if (nw.le.nwtn) then
           !Label the first nwn wavefunctions as neutron wavefunctions.
          Isospin=-1
      else
          !Label the rest as protons.
          Isospin=1
      endif
      
      HFBasis(nw) = NewWaveFunction(wf, Isospin, 1 ,kparz(nw), 1,1 )

      !Setting quantumnumbers, energies and occupation factors. 
      !(Factor two in the occupation factors)
      call HFBasis(nw)%SetOcc(2.0_dp* v2(nw))
      call HFBasis(nw)%SetEnergy(esp1(nw))
      !NIL8 and EV8 only stores the + signature WF's
      call HFBasis(nw)%SetSignature(1) 
    enddo
    close(IChan)
  end subroutine ReadNil    
  
  subroutine ReadEV4 (IChan, InputFilename)
  !-----------------------------------------------------------------------------
  ! Routine that reads EV4 .wf files. These wavefunctions do not necessarily
  ! conserve parity and thus do not carry a parity quantum number.
  !
  ! Taken from ev4.1.2.5 on 23/12/13.(Except for the option of reading from 
  ! smaller boxes) Quantities that are not useful at the moment of creation 
  ! are commented out and not saved.
  !
  ! Options:
  !   IChan is the number from which the data should be read.
  !-----------------------------------------------------------------------------
    use CompilationInfo
    use GenInfo
    use WaveFunctions
    use Moments
    use SpwfStorage

    implicit none
    
    character(len=*), intent(in) :: InputFilename
    integer, intent(in) :: IChan
    
    integer             :: nwave, nwavep, nwaven, npn, npp, Isospin,nnx,nny,nnz
    integer             :: iwave, itert,icqx,iq1ev4,iq2ev4
    real(KIND=dp)       :: v2(nwt), esp1(nwt), WF1(nx,ny,nz), WF2(nx,ny,nz)
    real(KIND=dp)       :: WF3(nx,ny,nz), WF4(nx,ny,nz),dx1,delq2,cq2=0.0_dp
    real(KIND=dp)       :: Constraints(2)
    type(Spinor)        :: wf
    type(Moment), pointer  :: Current
   
    open (IChan,form='unformatted',file=InputFileName)
    read (IChan) ! iver
    read (IChan) ! head

    read  (IChan) nwaven,nwavep,npn,npp
    read  (IChan) nnx,nny,nnz

    if((nnx.ne.nx) .or. (nny.ne.ny) .or. (nnz.ne.nz)) then
          call stp('Wrong dimensions on input. ')
    endif
    
    nwave=nwaven+nwavep        
    if(nwave.ne.nwt) call stp('Number of wavefunctions does not match')

    read (IChan) itert,dx1
    
    if (dx1.ne.dx) then
          print *, 'WARNING! DX on file is not the same as DX on input!'
    endif
    
    read (IChan) !cqcm,qcmabs,qcmcst,qcm,pentcm,ecmcst
    read (IChan) !irtd,irt2,irt3,irt4,imt2,imt3,imt4,isod
    read (IChan) !iqd,delqd,cqd,qdcst,qdd,pented,edcst
    read (IChan) icqx,iq1ev4,iq2ev4,delq2,cq2 !,rcut,
   !1            !cqx,qxcst,qxfin,qxx,pentex,excst,
   !2           cqy,qycst,qyfin,qyy,pentey,eycst,
   !3           cqz,qzcst,qzfin,qzz,pentez,ezcst,
   !4           q2cst,q2fin,g2fin
    read (IChan) !iq30,iq32,delq3,q3cut,
   !1           cq3,q30cst,q30fin,q30,pent30,e30cst,
   !2               q32cst,q32fin,q32,pent32,e32cst,q3cst,q3fin
    read (IChan) !iq4,delq4,q4cut,
   !1           cq4,q4cst,q4fin,q44,pente4,e4cst
   
    !Overriding the MOCCa constraints when the EV4 calculation was 
    ! constrained
    read (IChan) !t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,wso,wsoq
    read (IChan) !npair,gn,gp,encut,epcut,dcut,ambda,
   !1           delmax,alpha2,ntqp
    read (IChan) !(eqp  (i),i=1,nwave)
    read (IChan) !(delta(i),i=1,nwave)

    read (IChan) esp1
    read (IChan) v2
    read (IChan) !(v22 (i),i=1,nwave)

    read (IChan) Density%Rho(:,:,:,1)
    read (IChan) Density%Rho(:,:,:,2)
    read (IChan) Density%Tau
    read (IChan) Density%NablaJ

    call ChangeNumberWaveFunctions(nwt)

    do iwave=1,nwave             
      read (IChan) WF1
      read (IChan) WF2
      read (IChan) WF3
      read (IChan) WF4
      !The EV4 wavefunctions are normalised to two.
      WF1=WF1/sqrt(2.0_dp)
      WF2=WF2/sqrt(2.0_dp)
      WF3=WF3/sqrt(2.0_dp)
      WF4=WF4/sqrt(2.0_dp)

      wf= NewSpinor()
      call wf%SetComponent(1,WF1)
      call wf%SetComponent(2,WF2)
      call wf%SetComponent(3,WF3)
      call wf%SetComponent(4,WF4)
      
      if (iwave.le.nwaven) then
           !Label the first nwn wavefunctions as neutron wavefunctions.
          Isospin=-1
      else
          !Label the rest as protons.
          Isospin=1
      endif
      
      HFBasis(iwave) = NewWaveFunction(wf, Isospin, 1,0, 1,1)
      !Setting quantumnumbers, energies and occupation factors. 
      !Since these are EV4 wavefunctions, Time Simplex and 
      !Time Reversal are conserved.
      call HFBasis(iwave)%SetOcc(2*v2(iwave))
      call HFBasis(iwave)%SetEnergy(esp1(iwave))               
    enddo
  end subroutine ReadEV4
  
  subroutine ReadCr8 (IChan, InputFileName)
  !-----------------------------------------------------------------------------
  ! Reads a CR8 file and tries to continue as best as possible from it.
  ! This often needs quite some 'massaging'.
  !-----------------------------------------------------------------------------
  ! Only compatible with the latest file format, iver = 8
  !-----------------------------------------------------------------------------
    use CompilationInfo
    use GenInfo
    use WaveFunctions
    use SpwfStorage
    use Pairing,   only : Fermi, LNLambda
    use Cranking,  only : Omega, ContinueCrank, CrankValues
    use Densities, only : Recalc

    integer, intent(in) :: Ichan
    character(len=*), intent(in) :: inputFileName
    !---------------------------------------------------------------------------
    !Non-important quantities that need to get read
    integer       :: iver, npair, i, j, k, it, mulb1, mulb2, ipa, ipa2, fsize
    integer       :: Parity, Signature, Isospin, nwtp, icut
    real(KIND=dp) :: vgn, vgp, ren, rep, dcut
    real(KIND=dp) :: p(nx,ny,nz,4), trash
    !---------------------------------------------------------------------------
    ! Important quantities that are necessary for MOCCa to continue the 
    ! calculation.
    integer       :: kparz(nwt), keta(nwt)
    real(KIND=dp) :: esp1(nwt) , v2(nwt), espro(nwt)
    
    !---------------------------------------------------------------------------
    ! Signal to the density routine that it should recompute densities at the
    ! start of the iterations. This needs to be done because CR8 does not write
    ! all relevant densities to fort.12. (The time-odd ones namely...)
    ! We also print a warning for the user.
     1 format ('|---------------------------------------------------------|',/,&
    &          '| Disclaimer: Do not try to compare energies at iteration |',/,&
    &          '| 0 of MOCCa directly with the final iteration of CR8.    |',/,&
    &          '| Since CR8 does not write time-odd densities to file,    |',/,&
    &          '| MOCCa will recalculate all densities at the start.      |',/,&
    &          '| This means that the MOCCa and CR8 energy might not be   |',/,&
    &          '| exactly the same due to the density damping.            |',/,&
    &          '|---------------------------------------------------------|')
    
    Recalc = .true.
    
    print 1

    ! Do some quick checks
    if(.not. PC) call stp("MOCCa can't transform CR8 files on input.")
    if(.not. SC) call stp("MOCCa can't transform CR8 files on input.")
    if(     TRC) then
      call stp("MOCCa can't conserve Time-reversal starting from a CR8 file.")
    endif

    open  (IChan,form='unformatted',file=InputFileName)

    read(IChan) iver
    if(iver.ne.9) call stp('ReadCR8 only reads iver=9 files.', 'iver', iver)

    read(Ichan) !head line from CR8
    read(IChan) nwtn,nwtp!,npn0,npp0 in CR8

    if(nwtn + nwtp.ne.nwt) then
      call stp('Number of wavefunctions not equal to the ones on CR8 file.') 
    endif

    read(IChan) !mx,my,mz 
    read(IChan) !itert,filedx
    read(IChan) !cqx,qxcst,q0xcst,qxx,pentex,excst,
    !        1             cqy,qycst,q0ycst,qyy,pentey,eycst,
    !        2             cqz,qzcst,q0zcst,qzz,pentez,ezcst,
    !        3             q2cst,q02cst,g02cst,
    !        4             iq1,iq2,delq,cq2,q2cut
    read(IChan) !imtd,icqx,imtg
    read(Ichan) !cqr,qrtc,qrcst,qrfint,pentert,ercstt
    read(Ichan) !nforce,nfunc,ndd,ngal,njmunu,ncm2,nmass,ncoex
    read(IChan) !t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq
    
    !---------------------------------------------------------------------------
    ! Fermi energy and Lipkin-Nogami parameter for solving the pairing equations
    read(IChan) npair,vgn,vgp,ren,rep,dcut, icut,Fermi,LNLambda
    !---------------------------------------------------------------------------

    if(ContinueCrank) then
      read(IChan) trash,trash,trash,Omega(3),trash
      !Cr8's definition gets an extra minus sign!
      Omega(3) = - Omega(3) 
    else
      read(Ichan)
    endif
    read(IChan) !epair,eproj
    read(IChan) (kparz(i),i=1,nwt) !Parity of the HF-wavefunctions
    read(IChan) (keta(i),i=1,nwt)  !Signature of the HF-wavefunctions
    read(IChan) (esp1(i),i=1,nwt)  !Single-particle energies
    read(IChan) !(eqp(i),i=1,nwt)
    read(IChan) (v2(i),i=1,nwt)    ! Occupations of the canonical basis
    read(IChan) (espro(i),i=1,nwt) ! Single-particle energies in the canbasis
    read(IChan) !(icv0(i),i=1,mqp)
    ! Number of wavefunctions in parity-signature blocks.
    ! (1,it) = (P =  1, S =  1)
    ! (2,it) = (P = -1, S =  1)
    ! (3,it) = (P =  1, S = -1)
    ! (4,it) = (P = -1, S = -1) 
    read(IChan) ((npar(i,it),i=1,4),it=1,2) 

    !---------------------------------------------------------------------------
    ! Prepare all the relevant pairing matrices
    ! DeltaCR8 = Delta on file
    ! xkap     = KappaHFB on file
    ! rrn      = RhoHFB on file (positive signature part)
    ! rrt      = RhoHFB on file (negative signature part)
    ! dr       = real part of the position-space pairing field
    ! di       = imaginary part of the position-space pairing field
    ! dlnr     = real part of the position-space LN-modified pairing field
    ! dlni     = imaginary part of the position-space LN-modified pairing field
    fsize = maxval(npar)
    allocate(xkap(fsize, fsize,2,2))  ; xkap  = 0.0_dp
    allocate(Deltacr8(fsize, fsize,2,2)) ; Deltacr8 = 0.0_dp
    allocate(rrn(fsize,fsize,2,2))    ; rrn = 0.0_dp
    allocate(rrt(fsize,fsize,2,2))    ; rrt = 0.0_dp
    allocate(dr(nx,ny,nz,2), di(nx,ny,nz,2)) ; dr = 0.0_dp ; di = 0.0_dp
    allocate(dlnr(nx,ny,nz,2), dlni(nx,ny,nz,2)) ; dlnr = 0.0_dp ; dlni = 0.0_dp
    allocate(cr8rho(nx,ny,nz,2)) ; cr8rho=0.0_dp
    allocate(cr8vtau(nx,ny,nz,2)); cr8vtau=0.0_dp
    allocate(cr8vdiv(nx,ny,nz,2)); cr8vdiv=0.0_dp
    
    do it=1,2
      do ipa=1,2
        mulb1 = npar(ipa,  it)
        mulb2 = npar(ipa+2,it)
        read (IChan) !((deltacr8(i,j,ipa,it),i=1,mulb1),j=1,mulb2)
        read (IChan) ((rrn(i,j,ipa,it),i=1,mulb1),j=1,mulb1)
        read (IChan) ((rrt(i,j,ipa,it),i=1,mulb2),j=1,mulb2)
        read (IChan) ((xkap(i,j,ipa,it),i=1,mulb1),j=1,mulb2)
        read (IChan) ((deltacr8(i,j,ipa,it),i=1,mulb1),j=1,mulb2)
        read (IChan) !((uvst(i,j,ipa,it),i=1,mulb1),j=1,mulb2)
      enddo
      read (IChan) dr(:,:,:,it)
      read (IChan) di(:,:,:,it)
      read (IChan) dlnr(:,:,:,it)
      read (IChan) dlni(:,:,:,it)
      
      !-------------------------------------------------------------------------
      ! Reading the special densities
      read (IChan) !(drhor(i,it),i=1,mv)
      read (IChan) !(drhoi(i,it),i=1,mv)
      read (IChan)! (dtaur(i,it),i=1,mv)
      read (IChan)! (dtaui(i,it),i=1,mv)
      read (IChan)! (dJxxr(i,it),i=1,mv)
      read (IChan)! (dJxxi(i,it),i=1,mv)
      read (IChan)! (dJxyr(i,it),i=1,mv)
      read (IChan)! (dJxyi(i,it),i=1,mv)
      read (IChan)! (dJxzr(i,it),i=1,mv)
      read (IChan)! (dJxzi(i,it),i=1,mv)
      read (IChan)! (dJyxr(i,it),i=1,mv)
      read (IChan)! (dJyxi(i,it),i=1,mv)
      read (IChan)! (dJyyr(i,it),i=1,mv)
      read (IChan)! (dJyyi(i,it),i=1,mv)
      read (IChan)! (dJyzr(i,it),i=1,mv)
      read (IChan)! (dJyzi(i,it),i=1,mv)
      read (IChan)! (dJzxr(i,it),i=1,mv)
      read (IChan)! (dJzxi(i,it),i=1,mv)
      read (IChan)! (dJzyr(i,it),i=1,mv)
      read (IChan)! (dJzyi(i,it),i=1,mv)
      read (IChan)! (dJzzr(i,it),i=1,mv)
      read (IChan)! (dJzzi(i,it),i=1,mv)
      do i=1,4
        k=npar(i,it)
        if (k.ne.0) then
          read (IChan) 
        endif
      enddo   
    enddo
    !---------------------------------------------------------------------------
    ! Reading normal densities
    read (IChan) cr8rho(:,:,:,1)
    read (IChan) cr8rho(:,:,:,2)
    read (IChan) cr8vtau
    read (IChan) cr8vdiv
    !---------------------------------------------------------------------------
    ! Allocate space for the wavefunctions and read them
    call ChangeNumberWaveFunctions(nwt)
    do i=1,nwt
      ! Read the components
      do j=1,4
        read(IChan) p(:,:,:,j)
      enddo
      !Assign quantum numbers
      Parity=kparz(i)
      Signature=keta(i)
      if(i.le.nwtn)  Isospin=-1
      if(i.gt.nwtn)  Isospin= 1

      call HFBasis(i)%SetIsospin(Isospin)
      call HFBasis(i)%SetParity(Parity)
      call HFBasis(i)%SetSignature(Signature)
      call HFBasis(i)%SetTimeSimplex(1)   
      call HFBasis(i)%SetTimeReversal(0) 
      ! Setting the HFB occupations as the HF occupations might not always be 
      ! the best way to do things.
      call HFBasis(i)%SetOcc(v2(i))  
      call HFBasis(i)%SetEnergy(esp1(i))
      
      if(keta(i).eq.1) then                    
          call HFBasis(i)%SetGrid(1,p(:,:,:,1))
          call HFBasis(i)%SetGrid(2,p(:,:,:,2))
          call HFBasis(i)%SetGrid(3,p(:,:,:,3))
          call HFBasis(i)%SetGrid(4,p(:,:,:,4))
      else
          !The states with negative signature have a different ordering of     
          ! components in CR8.
          call HFBasis(i)%SetGrid(3,p(:,:,:,1))
          call HFBasis(i)%SetGrid(4,p(:,:,:,2))
          call HFBasis(i)%SetGrid(1,p(:,:,:,3))
          call HFBasis(i)%SetGrid(2,p(:,:,:,4))                
      endif 
    enddo
!     allocate(cr8A(nx,ny,nz,3,2), cr8S(nx,ny,nz,3,2))
!     allocate(CR8u(nx,ny,nz,2))
!     read(Ichan) Cr8u(:,:,:,1)
!     read(Ichan) Cr8u(:,:,:,2)
!     read(Ichan) cr8S(:,:,:,1,:)
!     read(Ichan) cr8S(:,:,:,2,:)
!     read(Ichan) cr8S(:,:,:,3,:)
!     read(Ichan) cr8A(:,:,:,1,:)
!     read(Ichan) cr8A(:,:,:,2,:)
!     read(Ichan) cr8A(:,:,:,3,:)

    ! Make sure the data is used in MOCCa, call the 'massaging' routine
    call prepareHFBModule
    call CR8state
    
    close(IChan)
  end subroutine ReadCr8
  
  subroutine WriteCr8 (OutputFileName)
  !-----------------------------------------------------------------------------
  ! Write a file for wave-functions that is readable by CR8.
  ! Only compatible with the latest file format, iver = 8
  !-----------------------------------------------------------------------------
    use CompilationInfo
    use GenInfo
    use WaveFunctions
    use SpwfStorage
    use Pairing,   only : Fermi, LNLambda
    use Cranking,  only : Omega, ContinueCrank, CrankValues
    use Densities, only : Recalc

    character(len=*), intent(in) :: outputFileName
    !---------------------------------------------------------------------------
    !Non-important quantities that need to get read
    real(KIND=dp) :: wf(nx,ny,nz,4), trash
    !---------------------------------------------------------------------------
    ! Important quantities that are necessary for MOCCa to continue the 
    ! calculation.
    integer       :: kparz(nwt), keta(nwt), s, p, it, i, ipa,j, k, Ichan
    real(KIND=dp) :: esp1(nwt) , v2(nwt), espro(nwt)
    
    ! Do some quick checks
    if(.not. PC) call stp("MOCCa can't transform CR8 files on input.")
    if(.not. SC) call stp("MOCCa can't transform CR8 files on input.")
    if(     TRC) then
      call stp("MOCCa can't conserve Time-reversal starting from a CR8 file.")
    endif

    Ichan = 13
    open  (IChan,form='unformatted',file=OutputFileName)

    write(Ichan) 8
    
    write(Ichan) !head line from CR8
    write(Ichan) nwn,nwp!,npn0,npp0 in CR8

    write(Ichan) !mx,my,mz 
    write(Ichan) !itert,filedx
    write(Ichan) !cqx,qxcst,q0xcst,qxx,pentex,excst,
    !        1             cqy,qycst,q0ycst,qyy,pentey,eycst,
    !        2             cqz,qzcst,q0zcst,qzz,pentez,ezcst,
    !        3             q2cst,q02cst,g02cst,
    !        4             iq1,iq2,delq,cq2,q2cut
    write(Ichan) !imtd,icqx,imtg
    write(Ichan) !cqr,qrtc,qrcst,qrfint,pentert,ercstt
    write(Ichan) !nforce,nfunc,ndd,ngal,njmunu,ncm2,nmass,ncoex
    write(Ichan) !t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq
    
    !---------------------------------------------------------------------------
    ! Fermi energy and Lipkin-Nogami parameter for solving the pairing equations
    write(Ichan) ! npair,vgn,vgp,ren,rep,dcut, icut,Fermi,LNLambda
    !---------------------------------------------------------------------------
    write(Ichan)
    
    write(Ichan) !epair,eproj
    
    npar = 0
    do i=1,nwt
        kparz(i) = HFBasis(i)%Parity
        keta(i)  = HFBasis(i)%Signature
        v2(i)    = 0.0_dp
        esp1(i)  = HFBasis(i)%Energy
        
        it = (HFBasis(i)%isospin+3)/2
        p  = kparz(i)
        s  = keta(i) 
        
        if( p.eq.1 .and. s.eq.1) then
            npar(1,it) = npar(1,it) + 1 
        elseif( p.eq.-1 .and. s.eq.1) then
            npar(2,it) = npar(1,it) + 1 
        elseif( p.eq.+1 .and. s.eq.-1) then
            npar(3,it) = npar(1,it) + 1 
        elseif( p.eq.-1 .and. s.eq.-1) then
            npar(4,it) = npar(1,it) + 1 
        endif
    enddo
    
    write(Ichan) (kparz(i),i=1,nwt) !Parity of the HF-wavefunctions
    write(Ichan) (keta(i),i=1,nwt)  !Signature of the HF-wavefunctions
    write(Ichan) (esp1(i),i=1,nwt)  !Single-particle energies
    write(Ichan) !(eqp(i),i=1,nwt)
    write(Ichan) (v2(i),i=1,nwt)    ! Occupations of the canonical basis
    write(Ichan) (esp1(i),i=1,nwt) ! Single-particle energies in the canbasis
    write(Ichan) !(icv0(i),i=1,mqp)
    ! Number of wavefunctions in parity-signature blocks.
    ! (1,it) = (P =  1, S =  1)
    ! (2,it) = (P = -1, S =  1)
    ! (3,it) = (P =  1, S = -1)
    ! (4,it) = (P = -1, S = -1) 
    write(Ichan) ((npar(i,it),i=1,4),it=1,2) 
    
    do it=1,2
      do ipa=1,2
        write(Ichan) !((deltacr8(i,j,ipa,it),i=1,mulb1),j=1,mulb2)
        write(Ichan) !((rrn(i,j,ipa,it),i=1,mulb1),j=1,mulb1)
        write(Ichan) !((rrt(i,j,ipa,it),i=1,mulb2),j=1,mulb2)
        write(Ichan) !((xkap(i,j,ipa,it),i=1,mulb1),j=1,mulb2)
        write(Ichan) !((deltacr8(i,j,ipa,it),i=1,mulb1),j=1,mulb2)
        write(Ichan) !((uvst(i,j,ipa,it),i=1,mulb1),j=1,mulb2)
      enddo
      write(Ichan) !dr(:,:,:,it)
      write(Ichan) !di(:,:,:,it)
      write(Ichan) !dlnr(:,:,:,it)
      write(Ichan) !dlni(:,:,:,it)
      
      !-------------------------------------------------------------------------
      ! Reading the special densities
      write(Ichan) !(drhor(i,it),i=1,mv)
      write(Ichan) !(drhoi(i,it),i=1,mv)
      write(Ichan)! (dtaur(i,it),i=1,mv)
      write(Ichan)! (dtaui(i,it),i=1,mv)
      write(Ichan)! (dJxxr(i,it),i=1,mv)
      write(Ichan)! (dJxxi(i,it),i=1,mv)
      write(Ichan)! (dJxyr(i,it),i=1,mv)
      write(Ichan)! (dJxyi(i,it),i=1,mv)
      write(Ichan)! (dJxzr(i,it),i=1,mv)
      write(Ichan)! (dJxzi(i,it),i=1,mv)
      write(Ichan)! (dJyxr(i,it),i=1,mv)
      write(Ichan)! (dJyxi(i,it),i=1,mv)
      write(Ichan)! (dJyyr(i,it),i=1,mv)
      write(Ichan)! (dJyyi(i,it),i=1,mv)
      write(Ichan)! (dJyzr(i,it),i=1,mv)
      write(Ichan)! (dJyzi(i,it),i=1,mv)
      write(Ichan)! (dJzxr(i,it),i=1,mv)
      write(Ichan)! (dJzxi(i,it),i=1,mv)
      write(Ichan)! (dJzyr(i,it),i=1,mv)
      write(Ichan)! (dJzyi(i,it),i=1,mv)
      write(Ichan)! (dJzzr(i,it),i=1,mv)
      write(Ichan)! (dJzzi(i,it),i=1,mv)
      do i=1,4
        k=npar(i,it)
        if (k.ne.0) then
          write(Ichan) 
        endif
      enddo   
    enddo
    !---------------------------------------------------------------------------
    ! Reading normal densities
    write(Ichan) !cr8rho(:,:,:,1)
    write(Ichan) !cr8rho(:,:,:,2)
    write(Ichan) !cr8vtau
    write(IChan) !cr8vdiv
    !---------------------------------------------------------------------------
    ! Allocate space for the wavefunctions and read them
    do i=1,nwt
      ! Read the components
      if(keta(i).eq.1) then
          wf(:,:,:,1) = HFBasis(i)%Value%Grid(:,:,:,1,1)                  
          wf(:,:,:,2) = HFBasis(i)%Value%Grid(:,:,:,2,1)                  
          wf(:,:,:,3) = HFBasis(i)%Value%Grid(:,:,:,3,1)                  
          wf(:,:,:,4) = HFBasis(i)%Value%Grid(:,:,:,4,1)                            
      else
          !The states with negative signature have a different ordering of     
          ! components in CR8.
          wf(:,:,:,3) = HFBasis(i)%Value%Grid(:,:,:,1,1)                  
          wf(:,:,:,4) = HFBasis(i)%Value%Grid(:,:,:,2,1)                  
          wf(:,:,:,1) = HFBasis(i)%Value%Grid(:,:,:,3,1)                  
          wf(:,:,:,2) = HFBasis(i)%Value%Grid(:,:,:,4,1)     
      endif 
      do j=1,4
        write(Ichan) wf(:,:,:,j)
      enddo
    enddo
    
    close(IChan)
  end subroutine writeCr8
  
  subroutine CR8State 
  !-----------------------------------------------------------------------------
  ! Subroutine that translates the data left by a CR8 run to data that is usable
  ! by MOCCa.
  ! Notably:
  ! *) dr & di   => PairingField
  ! *) rrn + rrt => RhoHFB
  ! *) xkap      => KappaHFB
  ! *) cr8rho    => Density%Rho
  ! *) cr8vtau   => Density%Tau
  ! *) cr8vdiv   => Density%NablaJ
  ! *) deltacr8  => Delta
  !-----------------------------------------------------------------------------
  ! In this routine, we will be making extensive use of the arrangement 
  ! properties of wavefunctions in CR8. Namely the HF-wavefunctions are come in 
  ! a very definite order. This order is maintained:
  ! 1) P = 1, S = 1
  ! 2) P =-1, S = 1
  ! 3) P = 1, S =-1
  ! 4) P =-1, S =-1
  ! This order is ok for every isospin, and the neutron wavefunction come before
  ! the proton wavefunctions. Note that MOCCa does not change the numbering, 
  ! and thus this order is maintained when MOCCa reads a CR8 file.
  !-----------------------------------------------------------------------------
  
  use Pairing

  integer          :: i, j, ii, jj, Count, offset, it, offseti, offsetj, P
  integer          :: ipa
  
  !-----------------------------------------------------------------------------
  !1) PairingFields constructed out of dr & di
  if(PairingType.eq.2) then
    PairingField = dcmplx(dr,di)
  endif
  !-----------------------------------------------------------------------------
  !2) RhoHFB constructed out of rrn and rrt
  if(PairingType.eq.2) then
    Count = 0
    do it=1,2
      !Positive parity, positive signature
      do i=1,blocksizes(2,it)/2
        do j=1,blocksizes(2,it)/2        
          RhoHFB(i,j,2,it) = cmplx(rrn(i,j,1,it),0.0_dp)
          RhoHFB(j,i,2,it) = RhoHFB(i,j,2,it)
        enddo
        Count = Count + 1
      enddo
      ! Negative parity, positive signature
      do i=1,blocksizes(1,it)/2
        do j=1,blocksizes(1,it)/2        
          RhoHFB(i,j,1,it) = cmplx(rrn(i,j,2,it),0.0_dp)
          RhoHFB(j,i,1,it) = RhoHFB(i,j,1,it)
        enddo
        Count = Count + 1
      enddo    
      ! Positive parity, negative signature
      do i=blocksizes(2,it)/2+1,blocksizes(2,it)
        do j=blocksizes(2,it)/2+1,blocksizes(2,it)
          RhoHFB(i,j,2,it) = cmplx(rrt(i-blocksizes(2,it)/2,j-blocksizes(2,it)/2,1,it),0.0_dp)
          RhoHFB(j,i,2,it) = RhoHFB(i,j,2,it)
        enddo
        Count = Count + 1
      enddo
      ! negative parity, negative signature
      do i=blocksizes(1,it)/2+1,blocksizes(1,it)
        do j=blocksizes(1,it)/2+1,blocksizes(1,it)
          RhoHFB(i,j,1,it) = cmplx(rrt(i-blocksizes(1,it)/2,j-blocksizes(1,it)/2,2,it),0.0_dp)
          RhoHFB(j,i,1,it) = RhoHFB(i,j,1,it)
        enddo
        Count = Count + 1
      enddo
    enddo
    
    if(Count.ne.nwt) call stp('Wrong assignment of RhoHFB from CR8 file')
    !Note, these four loops above have all been checked on correct signature
    !and parity. Only the counting check is still implemented, since that is 
    ! easy on the eyes.
  endif
  !-----------------------------------------------------------------------------
  !3) Construct KappaHFB from the xkap matrix 
  if(PairingType.eq.2) then
    do it=1,2
      ! Positive parity part (and thus opposite signatures)      
      do i=1,blocksizes(2,it)/2
        do j=1,blocksizes(2,it)/2
          KappaHFB(i,j + blocksizes(2,it)/2,2,it)   = DCMPLX(xkap(i,j,1,it),0.0_dp)
          KappaHFB(j   + blocksizes(2,it)/2,i,2,it) = - KappaHFB(i,j + blocksizes(2,it)/2,2,it)
        enddo
      enddo
      
      !Negative parity part (and thus opposite signatures)     
      do i=1,blocksizes(1,it)/2
        do j=1,blocksizes(1,it)/2
          KappaHFB(i,j + blocksizes(1,it)/2,1,it)   = DCMPLX(xkap(i,j,2,it),0.0_dp)
          KappaHFB(j   + blocksizes(1,it)/2,i,1,it) = - KappaHFB(i,j + blocksizes(1,it)/2,1,it)
        enddo
      enddo    
    enddo
  endif
  !-----------------------------------------------------------------------------
  ! 4) Read the densities
  ! *) cr8rho    => Density%Rho
  ! *) cr8vtau   => Density%Tau
  ! *) cr8vdiv   => Density%NablaJ
  Density%Rho = cr8rho
  
  !Don't forget the derivatives of the densities
  do it=1,2                                                                          
      Density%DerRho(:,:,:,1,it) = &
      & DeriveX(Density%Rho(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
      Density%DerRho(:,:,:,2,it) = &
      & DeriveY(Density%Rho(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
      Density%DerRho(:,:,:,3,it) = &
      & DeriveZ(Density%Rho(:,:,:,it), ParityInt,SignatureInt,TimeSimplexInt,1)
      Density%LapRho(:,:,:,it)   = &
      & Laplacian(Density%Rho(:,:,:,it),ParityInt,SignatureInt,TimeSimplexInt,1)
  enddo
  
  Density%Tau = cr8vtau
  Density%NablaJ = cr8vdiv
  
  !-----------------------------------------------------------------------------
  ! 5) DeltaCR8 => Delta
  if(PairingType.eq.2) then
    do it=1,2
      ! Positive parity part (and thus opposite signatures)
      do i=1,blocksizes(2,it)/2
        do j=1,blocksizes(2,it)/2
          Delta(i,j+blocksizes(2,it)/2,2,it) = Deltacr8(i,j,1,it)
          Delta(j+blocksizes(2,it)/2,i,2,it) = -  Delta(i,j+blocksizes(2,it)/2,2,it)
        enddo
      enddo
      
      !Negative parity part (and thus opposite signatures)
      do i=1,blocksizes(1,it)/2
        do j=1,blocksizes(1,it)/2
          Delta(i,j+blocksizes(1,it)/2,1,it) = DeltaCR8(i,j,2,it)
          Delta(j+blocksizes(1,it)/2,i,1,it) = -  Delta(i,j+blocksizes(1,it)/2,1,it)
        enddo
      enddo       
    enddo
  endif
 
  end subroutine CR8state
   
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
            if(Lipkin) then
                npair = 5
            else
                npair = 4
            endif
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
            write(iunit) (DensityBasis(i)%GetParity(), i=1,nwt)
        endif

        ! Single particle energies
        write(iunit) (DensityBasis(i)%GetEnergy(),i=1,nwt)
        !Occupation numbers
        if(TRC) then
            write(iunit) (DensityBasis(i)%GetOcc()/2.0_dp,i=1,nwt)
        else
            write(iunit) (DensityBasis(i)%GetOcc(),i=1,nwt)
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
            write(iunit) DensityBasis(i)%Value%Grid(:,:,:,1,1)*sqrt(2.0_dp)
            write(iunit) DensityBasis(i)%Value%Grid(:,:,:,2,1)*sqrt(2.0_dp)
            write(iunit) DensityBasis(i)%Value%Grid(:,:,:,3,1)*sqrt(2.0_dp)
            write(iunit) DensityBasis(i)%Value%Grid(:,:,:,4,1)*sqrt(2.0_dp)
        enddo

        close(iunit)

  end subroutine writepromesse
end module Interfaces
