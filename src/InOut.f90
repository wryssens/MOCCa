module InOutput

  use CompilationInfo
  use GenInfo
  use Force
  use Wavefunctions
  use Pairing
  use Transform
  use SpwfFactory
  
  implicit none
  public
  
  save
  !-----------------------------------------------------------------------------
  ! MOCCa variables
  !-----------------------------------------------------------------------------
  ! Name of the input file and output file. Default = fort.12 and fort.13
  character(len=256),public :: InputFileName="fort.12", OutputFileName="fort.13"
  !-----------------------------------------------------------------------------
  !logicals determining the symmetries of the input file
  logical ::  inTRC, inTSC,inIC,inPC,inSC
  !-----------------------------------------------------------------------------
  ! Parameters of the input file.
  integer :: filenx,fileny,filenz, filenwt
  !-----------------------------------------------------------------------------
  !Determining if MOCCa needs to write extra information to the output file, 
  !for debugging purposes. Not active at the moment.
  logical :: ExtraOutput=.false.
  !-----------------------------------------------------------------------------
  ! Whether or not to write densities to file at the end of the iterations
  logical, public    :: Pictures=.false.
  !-----------------------------------------------------------------------------
  ! Temporary storage for Kappa matrix, in order to be able to 
  ! continue HFB calculations. There needs to be some kind of temporary storage,
  ! in case the symmetries on file do not match those of the run.
  ! The array FileBlocksizes handles the effective sizes of the matrices
  complex(KIND=dp), allocatable :: InputKappa(:,:,:,:)
  integer, allocatable          :: FileBlocksizes(:,:)

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

  subroutine Input()
  !-----------------------------------------------------------------------------
  ! Master subroutine that uses the other routines to correctly read input.
  ! If allowed and necessary, MOCCa will transform the data to allow for
  ! calculations with different symmetry combinations.
  !-----------------------------------------------------------------------------
    
    1 format (' *********************************************************'  ,/,&
    &         ' * ATTENTION                                             *'  ,/,& 
    &         ' * Wavefunction file did not correspond to symmetries    *'  ,/,&
    &         ' * on input. MOCCa made the following changes :          *'  ,/,& 
    &         ' *                FILE      INPUT                        *'  ,/,& 
    &         ' * Parity          ',L1, '         ',L1,27(' '),        '*'  ,/,& 
    &         ' * Signature       ',L1, '         ',L1,27(' '),        '*'  ,/,& 
    &         ' * Time Simplex    ',L1, '         ',L1,27(' '),        '*'  ,/,&
    &         ' * Time Reversal   ',L1, '         ',L1,27(' '),        '*'  ,/,&                                                                
    &         ' * Isospin         ',L1, '         ',L1,27(' '),        '*'  ,/,&
    &         ' *                                                       *'  ,/,&
    &         ' * nx             ',i3,'       ',i3,26(' '),            '*'  ,/,&
    &         ' * ny             ',i3,'       ',i3,26(' '),            '*'  ,/,&
    &         ' * nz             ',i3,'       ',i3,26(' '),            '*'  ,/,&
    &         ' * nwt            ',i3,'       ',i3,26(' '),            '*'  ,/,&
    &         ' *********************************************************')        
    !---------------------------------------------------------------------------
    !Reading the input from data; and allocate all necessary arrays.
    call ReadRunInfo()
    !---------------------------------------------------------------------------
    !MOCCa reads the wavefunction file supplied, after determining what type 
    !of file it exactly is.
    call WFInput()
    !---------------------------------------------------------------------------
    if(    (inPC .neqv. PC)  .or. (inTSC .neqv. TSC) .or. (inSC .neqv. SC)     &
    & .or. (inTRC.neqv. TRC) .or. (inIC  .neqv. IC)) then
      !If the symmetries requested are not the ones stored on the wavefunction
      !file, transform all the wavefunctions and densities.
      call TransformInput(inTRC,inTSC,inIC,inPC,inSC,                          &
      &                   filenx,fileny,filenz,filenwt)
      print 1, inPC, PC, inSC, SC, inTSC, TSC, inTRC, TRC, inIC, IC,           &
      &       filenx,nx,fileny,ny,filenz,nz,filenwt, nwt
    endif
    ! Special mention here for the HFB anomalous density matrix.
    if(allocated(InputKappa)) then
      KappaHFB = TransformKappa(InputKappa, inPC, inIC, FileBlocksizes)
    endif
  end subroutine Input
  
  subroutine Output
  !-----------------------------------------------------------------------------
  ! Subroutine that governs the output of a wavefunction file.
  !-----------------------------------------------------------------------------
    integer           :: OutputChannel
    
    !Getting an open channel for output
    call get_unit(OutputChannel)
    
    !Writing to the output channel
    call writeMOCCa(OutputChannel)
    
  end subroutine Output
  
  subroutine WFInput()
  !----------------------------------------------------------------------------- 
  ! Subroutine that identifies the .wf file that MOCCa gets as input.
  ! Identification is based on filename, the first letters of the filename 
  ! should match the code that was used to write the file.
  !
  ! EV8 or NIL8 => ReadNil8
  ! CR8         => ReadCr8
  ! EV4         => ReadEV4
  ! MOCCa       => MOCCa
  !-----------------------------------------------------------------------------
    integer             :: InputChannel
       
    !Getting an open channel for input
    call get_unit(InputChannel)
    
    !MOCCa will try to recognise the wavefunction input it gets and call the 
    ! corresponding input subroutine
    if( &
    &     (InputFileName(1:3).eq."ev8")  &
    & .or.(InputFileName(1:4).eq."nil8") &
    & .or.(InputFileName(1:4).eq."int8") &
    & ) then
    ! MOCCa is dealing with output from either NIL8,INT8 or EV8
            call ReadNil(InputChannel)
    elseif(&
    &     (InputFileName(1:3).eq."ev4")  &
    & .or.(InputFileName(1:4).eq."int4") &
    & .or.(InputFileName(1:4).eq."d8a4") &
    & ) then
    !MOCCa is dealing with input from either EV4, D8A4 or INT4
        call ReadEv4(InputChannel)

    elseif(&
    &    ((InputFileName(1:3).eq."cr8"    ) &
    & .or.(InputFileName(1:4).eq."evcr"   ) &
    & .or.(InputFileName(1:7).eq."permcr8")) &
    & .and.(InputFileName(1:5).ne."cr8to")) then
    ! MOCCa is dealing with output from cr8, evcr8 or permcr8.
    ! (It can't be cr8to4 or cr8to1)
        call ReadCR8(InputChannel)
    elseif(&
    &      InputFilename(1:5).eq."MOCCa" ) then
        !MOCCa is dealing with its own type of input
        call ReadMOCCa(InputChannel)    
    else
      ! The input file is not clearly named; MOCCa defaults to its own 
      ! subroutine.
           call ReadMOCCa(InputChannel)
    endif
  end subroutine WFInput

  subroutine ReadRunInfo
  !-----------------------------------------------------------------------------
  ! This subroutine reads the relevant data from the user provided file.
  ! Most, if not all, of the initialisation is overseen by this routine, as 
  ! every array needs to be allocated etc...
  !-----------------------------------------------------------------------------
  ! Order of the namelists:
  ! 
  ! 1) NML = GenInfo     (module Geninfo)
  ! 2) NML = Spwfstorage (module Spwfstorage)
  ! 3) NML = Densit      (module Densities)
  ! 4) NML = Pairing     (module Pairing) 
  ! 5) NML = Derivatives (module Derivatives)
  ! 6) NML = MomentParam (module Moments)
  ! If necessary:
  !   7) NML = MomentConstraints (module Moments)
  ! 8) NML = Inandoutput (module Inoutput)
  ! 9) NML = Force       (module Force)
  ! 10)NML = Coulomb     (module Coulomb)
  ! 11)NML = Cranking    (module Cranking)
  !-----------------------------------------------------------------------------
    use Force, only       : ReadForceInfo
    use Densities, only   : ReadDensitInfo, Iniden
    use Moments, only     : ReadMomentData
    use Cranking, only    : ReadCrankingInfo
    use Coulomb, only     : ReadCoulombInfo
    use SpwfStorage, only : ReadSpwfStorageInfo
    use Mesh, only        : IniMesh
    use Cranking, only    : ReadCrankingInfo
    use Pairing, only     : ReadPairingInfo
    use Derivatives,only  : ReadDerivativesInfo
    
    implicit none
    
    NameList /InAndOutput/ InputFileName,OutputFileName,ExtraOutput,Pictures
    !--------------- Reading Input---------------------------------------   
    !Info for the GenInfo Module
    call ReadGenInfo()  
    
    !Info for the SpwfstorageModule
    call ReadSpwfStorageInfo()
    
    call ReadDensitInfo()
    
    !Allocating the densities and setting them to zero
    call Iniden
    
    !Allocating and initialising the Mesh variables
    call IniMesh
    
    !Info for the Pairing module
    call ReadPairingInfo()

    ! Info on the Derivatives Module
    call ReadDerivativesInfo()
    
    !Info on the Moments Module
    call ReadMomentData()

    !Reading the names for the in- and outputfiles.
    read (unit=*, nml=InAndOutput)
        
    ! Force namelist
    call ReadForceInfo()
    
    !Coulomb namelist
    call ReadCoulombInfo()

    !Reading the Cranking Namelist
    call ReadCrankingInfo()
    return
  end subroutine ReadRunInfo
  
  subroutine PrintInput
  !-----------------------------------------------------------------------------
  ! This subroutine prints all relevant information of the input, both user 
  ! input and wavefunction file.
  !-----------------------------------------------------------------------------
    use SpwfStorage, only: nwt, BroydenOrder
    use Derivatives, only: MaxFDOrder, MaxFDLapOrder
    use Densities, only  : DampingParam, PulayOrder, MixingScheme

    character(len=9) :: Con='Conserved', Broken='Broken   '

    1 format ( 20('-'), 'General Information ', 20('-'))
    2 format ( 'Mesh parameters' ) 
    3 format ( '   nx = ', i5 , ' ny = ' , i5 , ' nz = ' , i5)
    4 format ( '   dx = ', f20.10,' (fm  ) ')
    5 format ( '   dv = ', f5.2,' (fm^3) ')
    6 format ( 'Nucleus')
    7 format ( '    N = ', f10.5  ,'  Z = ', f10.5) 
    8 format ( 'Wavefunctions')
    9 format ( '  nwt = ', i5)
    10 format ( 'Iterative process')
    11 format ( '   dt = ', f5.2 ,' (10^{-22} s) ')
    12 format ( '  Iter= ', i5)
   121 format ( '  Imaginary Time-step method')
   122 format ( '  Nesterov optimal gradient')
    13 format ( '  Mixing Process : ' , A7, ' mixing')
   131 format ( '    Memory       : ' , i3)
   132 format ( '    Damping      : ' , f5.2) 
   133 format ( '    MixingScheme : ' , i2) 
    
    14 format ( 'Convergence Levels')
    15 format ( '  dE  =', e8.1)
    16 format ( '  dQ  =', e8.1)
    17 format ( 'Symmetries')
    18 format ( "  Parity      : ", a9) 
    19 format ( "  Signature   : ", a9) 
    20 format ( '  TimeReversal: ', a9)
    21 format ( '  TimeSimplex : ', a9)
    22 format ( '  Isospin     : ', a9)
    23 format ( 'Derivatives')
    24 format ( '  1st Order   : ', i1)
    25 format ( '  2nd Order   : ', i1)
    26 format ( '  Lagrange      ')
      
    99 format(" Force used on input file: ") 
   100 format("Name of the Force:  ", a57)
   101 format("t0 =", f12.3 , " x0  =", f12.3 , /,                             &
    &              't1 =', f12.3,  " x1  =", f12.3 , /,                        &
    &              "t2 =", f12.3 , " x2  =", f12.3 , / ,                       &
    &              "t3a=", f12.3 , " x3a =", f12.3 , " yt3a=", f12.3 ,/ ,      &
    &              "t3b=", f12.3 , " x3b =", f12.3 , ' yt3b=', f12.3 ,/,       &
    &              "te =", f12.3 , " to  =", f12.3 ,                  /,       &
    &              "wso=", f12.3 , " wsoq=", f12.3)
   102 format('Force Options' , /                                              &
    &              '  nmass   = ', i2,/                                        &
    &              '  COM1Body= ', i2,/                                        &
    &              '  COM2Body= ', i2)
    
    print 1 
    print 2
    print 3 , nx, ny, nz 
    print 4 , dx
    print 5 , dv
    print 6
    print 7 , neutrons, protons
    print 8 
    print 9 , nwt
    print 10 
    print 11, dt
    print 12, MaxIter

    call to_upper(IterType, IterType)
    select case(IterType)
    case('IMTS')
      print 121
    case('NEST')
      print 122
    end select

    if(PulayOrder.gt.1 .and. Mixingscheme.eq.1) then
      print 13, 'DIIS   '
      print 131, PulayOrder
    else
      print 13, 'Linear '  
      print 132, DampingParam
    endif
    print 133, MixingScheme
    
    print 14
    print 15, EnergyPrec
    print 16, MomentPrec
    print 17
    
    !Printing the conservation or breaking of symmetries
    if (PC)       print 18, Con
    if (.not. PC) print 18, Broken
    if (SC)       print 19, Con
    if(.not.SC)   print 19, Broken
    if (TRC)      print 20, Con
    if (.not.TRC) print 20, Broken
    if (TSC)      print 21, Con
    if(.not.TSC)  print 21, Broken
    if (IC)       print 22, Con
    if(.not.IC)   print 22, Broken
    
    print 23
    if( MaxFDOrder.ne. -1) then
      print 24, MaxFDOrder
      print 25, MaxFDLapOrder
    else
      print 26
    endif         
  end subroutine PrintInput
  
  subroutine ReadNil(IChan)
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
  
    !All symmetries are conserved in NIL/EV8
    inTRC=.true.; inIC=.true.; inPC=.true.; inSC=.true.; inTSC=.true.

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
  
  subroutine ReadEV4 (IChan)
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
    
    !---------------------------------------------------------------- -----------
    !All symmetries except parity are conserved in EV4
    inTRC=.true.; inIC=.true.; inPC=.false.; inSC=.true.; inTSC=.true.
    
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
  
  subroutine ReadCr8 (IChan)
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
    
    !---------------------------------------------------------------------------
    !All symmetries except time reversal are conserved in CR8
    inTRC=.false.; inIC=.true.; inPC=.true.; inSC=.true.; inTSC=.true.

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
   
  subroutine ReadMOCCa(Ichan)
  !-----------------------------------------------------------------------------
  ! Read info from an input file in the MOCCa format.
  !-----------------------------------------------------------------------------
    
    use geninfo
    use SpwfStorage, only : nwt, HFBasis, ChangeNumberWaveFunctions,CanBasis
    use Derivatives, only : MaxFDOrder, MaxFDLapOrder, CoulombLapOrder
    use Pairing,     only : PairingType, PairingStrength, alpha, PairingCut
    use Force
    use Cranking,    only : Omega,CrankValues,ContinueCrank
    use Moments,     only : ReadMoment, ContinueMoment
    use Densities 
    
    integer, intent(in)  :: IChan
    
    logical :: exists
    integer :: ioerror, N, Pin,Iin
    
    !File parameters to compare against
    integer       :: fileneutrons,fileprotons,i
    real(KIND=dp) :: filedx, fileintensity(4)
        
    !---------------------------------------------------------------------------
    ! Checking for a valid input file
    inquire(file=inputfilename, exist=exists)
    if(.not.exists) call stp('Input file specified does not exist!')
    open (IChan,form='unformatted',file=InputFileName)
    !---------------------------------------------------------------------------
    ! Reading the essential variables. 
    ! 1) Mesh size: nx,ny,nz,dx,dt
    ! 2) Symmetries: TRC,PC,SC,TSC,IC
    ! 3) Particles: nwt,neutrons, protons
    read(IChan, iostat=ioerror) filenx,fileny,filenz, filedx
    if(ioerror.ne.0) call stp('Input error for the mesh parameters.',          &
    &                         'ioerror = ', ioerror)
    read(IChan, iostat=ioerror) inTRC,inPC,inSC,inTSC,inIC
    if(ioerror.ne.0) call stp('Input error for the symmetries.',               &
    &                         'ioerror = ', ioerror)
    read(IChan,iostat=ioerror) filenwt,fileneutrons,fileprotons
    if(ioerror.ne.0) call stp('Input error for the wavefunction variables.',   &
    &                         'ioerror = ', ioerror)
    
    !---------------------------------------------------------------------------
    ! Checking the essential variables:
    ! a) Symmetries must not be unbroken.
    ! b) Mesh parameters & number of wavefunctions must correspond
    ! c) There should be enough space for neutrons & protons
    if(.not.inPC .and. PC) then
      call stp(' Parity broken on file, but parity conservation asked.')
    endif
    if(.not.inTSC .and. TSC) then
      call stp('Signature broken on file, but signature conservation asked.')
    endif
    if(.not.inSC .and. SC) then
      call stp('Time Simplex broken on file,'                                  &
      &      //'but Time Simplex conservation asked.')
    endif
    if(.not.inTRC .and. TRC) then
      call stp(' Time Reversal broken on file,'                                &
      &      //' but Time Reversal conservation asked.')
    endif
    
    if(filenx.ne.nx) then
      if( 2*filenx.ne.nx .or. SC ) then
        call stp(' nx is not compatible between file and data.',               &
        &        'on file: ', filenx                                           &
        &      , 'in data: ', nx     )
      endif
    endif
    if(SC.neqv.inSC .and. 2*filenx.ne.nx)then
        call stp(' nx is not compatible between file and data.',               &
        &        'on file: ', filenx                                           &
        &      , 'in data: ', nx     )
    endif
    if(fileny.ne.ny ) then
      if( 2*fileny.ne.ny .or. TSC ) then
        call stp(' ny is not compatible between file and data.',               &
        &        'on file: ', fileny                                           &
        &      , 'in data: ', ny     )
      endif
    endif
    if(TSC.neqv.inTSC .and. 2*fileny.ne.ny)then
        call stp(' ny is not compatible between file and data.',               &
        &        'on file: ', fileny                                           &
        &      , 'in data: ', ny     )
    endif
    if(filenz.ne.nz ) then
      if( 2*filenz.ne.nz .or. PC ) then
        call stp(' nz is not compatible between file and data.',               &
        &        'on file: ', filenz                                           &
        &      , 'in data: ', nz    )
      endif
    endif
    if(PC.neqv.inPC .and. 2*filenz.ne.nz)then
        call stp(' nz is not compatible between file and data.',               &
        &        'on file: ', filenz                                           &
        &      , 'in data: ', nz     )
    endif
    if( nwt.gt.filenwt .and. nwt.ne.2*filenwt) then 
      call stp('Nwt is not compatible between file and data.',                 &
      &        'In data: ', nwt, 'On file :', filenwt )
    endif
    
    N = nwt
    if(TRC) N = 2*N    
    if( protons.gt.N) call stp('Not enough spwfs to accomodate all protons!' )
    if(neutrons.gt.N) call stp('Not enough spwfs to accomodate all neutrons!')
    
    !---------------------------------------------------------------------------
    !    ! 4) WaveFunctions
    call ChangeNumberWaveFunctions(filenwt)
    do i=1,filenwt
      call HFBasis(i)%Read(IChan,filenx,fileny,filenz)
    enddo
    !---------------------------------------------------------------------------
    ! 5) Densities
    call ReadDensity(Density,Ichan,filenx,fileny,filenz, ioerror)
    if(ioerror.ne.0) then
        call stp('ReadDensity failed','Iostat', ioerror)
    endif
    !---------------------------------------------------------------------------
    ! Decide whether to cut wavefunctions or not
    call DecideToCut(nwt,filenwt)
    !---------------------------------------------------------------------------
    ! Reading the  Non-Essential Variables
    ! 6 ) Specifics of the finite difference scheme
    ! 7 )  Force  : 
    !      a) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq
    !      b) afor,hbar,hbm,xm,njmunu,nmass,COM1Body,COM2Body
    !      c) Functional Coefficients
    ! 8 )  Pairing: type, strength, alphas and PairingCut
    ! 9 )  Cranking effective, true values & intensity
    ! 10)  Constraint Values
    ! 11)  Get the initial guess for KappaHFB
    ! 12)  Get an intial guess for the Fermi levels
    ! 13)  Initial guess for Lipkin-Nogami parameter
    ! 14)  Get Momentparameters
    !---------------------------------------------------------------------------
    ! Maybe do some checks on these later?
    read(IChan, iostat=ioerror)    
    !Straight force parameters
    read(IChan,iostat=ioerror)
    !Extra parameters
    read(IChan,iostat=ioerror) 
    !Functional coefficients, not used at the moment.
    read(IChan,iostat=ioerror)
    !Pairing parameters
    read(IChan,iostat=ioerror)
    
    !Cranking Variables
    if(ContinueCrank) then
      read(IChan,iostat=ioerror) Omega, CrankValues
      if(ioerror.ne.0) then
        call stp('Did not read CrankValues correctly' , &
             &   'Iostat', ioerror)
      endif      
    else
      read(IChan,iostat=ioerror)
    endif
    
    if(PairingType.eq.2) then
      ! Read Kappa & Fermi energy only if needed
      ! But first decide on the appropriate sizes for InputKappaHFB
      PIn = 1 ; if(inPC) PIn = 2
      IIn = 1 ; if(InIc) IIn = 2
      allocate(FileBlocksizes(PIn, IIn)) ; FileBlocksizes = 0.0_dp
      if(inTRC) then
        allocate(InputKappa(2*filenwt,2*filenwt,Pin,IIn))
      else
        allocate(InputKappa(filenwt,filenwt,Pin,IIn))
      endif

      read(ICHan,iostat=ioerror) InputKappa     
      if(ioerror.ne.0) then
        call stp('Did not read InputKappa correctly' , &
             &   'Iostat', ioerror)
      endif
      read(Ichan,iostat=ioerror) FileBlockSizes
      if(ioerror.ne.0) then
        call stp('Did not read FileBlocksizes correctly' , &
             &   'Iostat', ioerror)
      endif

      read(ICHan, iostat=ioerror)  Fermi
      read(Ichan, iostat=ioerror)  LNLambda
      if(ioerror.ne.0) then
        call stp('Did not read Fermi and LNLambda correctly' , &
             &   'Iostat', ioerror)
      endif
    
    elseif(PairingType.eq.1) then
        !Read only Fermi in the BCS case
        read(IChan,iostat=ioerror)
        read(IChan,iostat=ioerror)  Fermi
        read(Ichan,iostat=ioerror) 
    else
        read(IChan,iostat=ioerror)
        read(IChan,iostat=ioerror)
        read(IChan,iostat=ioerror)
    endif    
    !Reading Moments
    if(ContinueMoment) then
      do while (.true.)
        call ReadMoment(IChan, ioerror)
        if(ioerror.eq.-1) exit !This is the End-Of-File condition
        if(ioerror.ne.0 ) then
          call stp('IO error in reading the multipole moments from file.',     &
          &        'ioerror = ', ioerror)
        endif
      enddo
    endif
    close (IChan)
  end subroutine ReadMOCCa
  
  subroutine WriteMOCCa(OChan)
  !-----------------------------------------------------------------------------
  ! Subroutine that writes the output to a MOCCa file.
  ! Note that the output of this file is not compatible with older MF codes.
  !
  !-----------------------------------------------------------------------------
  ! Essential variables are all the variables that are required to make a 
  ! calculation. The non-essential variables are all 'extra luxuries'.
  ! I hope, by making this distinction that wavefunction files will always
  ! be 'minimally functional', independent of changes to the MOCCa
  ! filestructure.
  !
  !...... Essential Variables ..................................................
  ! 1) Mesh size: nx,ny,nz,dx,dt
  ! 2) Symmetries: TRC,PC,SC,TSC,IC
  ! 3) Particles: nwt,neutrons, protons
  ! 4) WaveFunctions
  ! 5) Current densities
  !
  !..... Non-Essential Variables ...............................................
  ! 6 ) Specifics of the finite difference scheme
  ! 7 )  Force  : 
  !      a) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq
  !      b) afor,hbar,hbm,xm,njmunu,nmass,COM1Body,COM2Body
  !      c) Functional Coefficients
  ! 8 )  Pairing: type, strength, alphas and PairingCut
  ! 9 )  Omega and cranking values
  ! 10)  A guess for the kappa matrix for starting HFB calculations.
  ! 11)  The block sizes of the KappaHFB matrix
  ! 12)  Fermi energies
  ! 13)  Lipkin-Nogami parameter  
  ! 14)  Constraint Values and other parameters
  !
  !..... Testing Variables .....................................................
  !  (None at the moment) 
  !
  !..... Things to maybe include ...............................................
  !  *) Canonical basis?
  !
  !-----------------------------------------------------------------------------
    use SpwfStorage
    use CompilationInfo
    use GenInfo
    use Force
    use Cranking
    use Moments
    use Derivatives, only : MaxFDOrder, MaxFDLapOrder, CoulombLapOrder
    use Pairing,     only : PairingType, PairingStrength, alpha, PairingCut
    use Densities
    
    integer, intent(in)   :: Ochan
    integer               :: i, io
    type(Moment), pointer :: Current
    
    open (OChan,form='unformatted',file=OutputFileName)
    
    !---------------------------------------------------------------------------
    ! Essential variables
    !---------------------------------------------------------------------------     
    !Parameters of the mesh
    write(OChan,iostat=io) nx,ny,nz, dx
    !Conservation of symmetries
    write(OChan,iostat=io) TRC,PC,SC,TSC,IC
    !Number of proton, neutron and total states. 
    write(OChan,iostat=io) nwt,neutrons,protons   
    !Wavefunctions
    do i=1,nwt
      call HFBasis(i)%Write(OChan)
    enddo
    call writeDensity(Density,Ochan)
    !---------------------------------------------------------------------------
    ! Non-Essential (Continuation) variables
    !---------------------------------------------------------------------------
    !Specifics of the Finite Difference scheme
    write(OChan,iostat=io) MaxFDOrder, MaxFDLapOrder, CoulombLapOrder
    !Straight Force parameters
    write(OChan,iostat=io) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq
    !Extra parameters
    write(OChan,iostat=io) afor,hbar,hbm,COM1Body,COM2Body
    !Functional coefficients
    write(OChan,iostat=io) B1,B2,B3,B4,B5,B6,B7a,B7b,B8a,B8b,Byt3a,Byt3b,B9,   &
    &                      B9q,B10,B11,B12a,B12b,B13a,B13b,B14,B15,B16,B17,B18,&
    &                      B19,B20,B21
    !Pairing Variables
    write(OChan,iostat=io) PairingType, PairingStrength, alpha, PairingCut
    !Cranking Variables
    write(OChan,iostat=io) Omega, CrankValues
    ! HFB anomalous density matrix
    ! Make a guess only when not doing HFB
    if(PairingType.ne.2) call WriteOutKappa(PairingType)
    write(OChan,iostat=io) KappaHFB
    ! Write the dimensions of this matrix to file, and note that it has been
    ! allocated in WriteOutKappa
    write(OChan,iostat=io) blocksizes
    write(OChan,iostat=io) Fermi
    write(OChan,iostat=io) LNLambda
    !---------------------------------------------------------------------------
    ! Multipole constraint variables
    ! Special treatment: got to loop over the multipole moments and only write
    ! the constrained ones.
    Current => Root    
    do while(.true.)
      if(Current%ConstraintType.ne.0) then
        call Current%Write(Ochan)
      endif        
      if(associated(Current%Next)) then
        Current => Current%Next
      else
         exit
      endif
    enddo
    close (OChan)
    
  end subroutine WriteMOCCa

  subroutine PlotDensity()
    !---------------------------------------------------------------------
    ! Writes densities to file in formatted output, for plotting purposes.
    ! 
    !
    !---------------------------------------------------------------------
    use Densities

    integer :: iunit,i,j,k

    call get_unit(iunit)
    open(iunit, File='density.X.dat')

    do k=1,nz
      do j=1,ny
          i = nx/2+1
          if(SC) i = 1
          write(iunit, '(3f10.5)') Density%Rho(i,j,k,1), Density%Rho(i,j,k,2),&
          &                       sum(Density%Rho(i,j,k,:))

      enddo
    enddo

    close(iunit)

    call get_unit(iunit)
    open(iunit, File='density.Y.dat')

    do k=1,nz
        j = ny/2+1
        if(TSC) j = 1
        do i=1,nx
          write(iunit, '(3f10.5)') Density%Rho(i,j,k,1), Density%Rho(i,j,k,2),&
          &                       sum(Density%Rho(i,j,k,:))
        enddo
    enddo

    close(iunit)

    open(iunit, File='density.Z.dat')

    k=nz/2+1
    if(PC) k =1
       do j=1,ny
        do i=1,nx
          write(iunit, '(3f10.5)') Density%Rho(i,j,k,1), Density%Rho(i,j,k,2),&
          &                       sum(Density%Rho(i,j,k,:))
        enddo
      enddo
    close(iunit)

  end subroutine PlotDensity
  
!   subroutine PlotDensity(Direction,Coord,Isospin)
!     !---------------------------------------------------------------------------
!     ! This subroutine uses the GnuForInterface module to plot a cutthrough view 
!     ! of the density.
!     ! Direction controls the direction along which the view is "cut", Coord the 
!     ! place of cutting.
!     !    Direction=1 : Y-Z plane is plotted, with X=Coordinate
!     !    Direction=2 : X-Z plane is plotted, with Y=Coordinate
!     !    Direction=3 : Y-Z plane is plotted, with Z=Coordinate
!     ! Isospin controls which density should be plotted:
!     !    Isospin=-1 : Neutron density
!     !    Isospin= 0 : Total Density
!     !    Isospin= 1 : Proton Density
!     !---------------------------------------------------------------------------
!     ! Filenames are Picturefilename + '.X/Y/Z.' + '.eps'
!     !---------------------------------------------------------------------------   
!     use GnuFor
!     use Mesh
!     use Densities
    
!     1 format (60('-'))
!     2 format ('Creating Density plot.')
!     3 format (' Cutthrough at ', A1, ' = ' , i2 )
  
!     integer, intent(in)        :: Direction, Coord,Isospin
!     integer                    :: ierrorData,ierrorCommand,i,j,k
!     real(KIND=dp), allocatable :: ToPlot(:,:,:)
!     character(len=6)           :: xlabel,ylabel
!     character(len=256)         :: PlotFileName

!     !These ierror are useful for handling I/O errors. See the documentation of 
!     !the GnuFor module.
!     ierrorData = 0 ;   ierrorCommand=0

!     print 1
!     print 2
!     select case(Direction)
!       case (1)
!         print 3, 'X' , Coord
!       case (2)
!         print 3, 'Y' , Coord
!       case(3)
!         print 3, 'Z' , Coord
!     end select
!     print 1     
!     ! Cut along the X-direction
!     if(Direction.eq.1) then
!       xlabel='z (fm)'
!       ylabel='y (fm)'
!       PlotFileName=adjustl(adjustr(PictureFileName) // '.X')
!       !Formatting the data that should be used by gnuplot
!       allocate(ToPlot(3,ny,nz))
!       do k=1,nz
!         ToPlot(2,:,k)=MeshY
!       enddo
!       do j=1,ny
!         ToPlot(1,j,:)=MeshZ
!       enddo
!       ! Choosing the correct density to plot and putting it into the correct 
!       !formatting.
!       if(Isospin.eq.-1)then
!         do k=1,nz
!           do j=1,ny
!               ToPlot(3,j,k)=Density%Rho(Coord,j,k,1)
!           enddo
!         enddo    
!       elseif(Isospin.eq.0) then
!         do k=1,nz
!           do j=1,ny
!               ToPlot(3,j,k)=Density%Rho(Coord,j,k,1)+Density%Rho(Coord,j,k,2)
!           enddo
!         enddo    
!       else
!         do k=1,nz
!           do j=1,ny
!               ToPlot(3,j,k)=Density%Rho(Coord,j,k,2)
!           enddo
!         enddo            
!       endif
!       ! Writing the gnuplot data to a file.
!       call write_xyzgrid_data("DensityProfile", ny,nz, ToPlot, ierrorData)
!     elseif(Direction.eq.2) then
!       xlabel='z (fm)'
!       ylabel='x (fm)'
!       PlotFileName=adjustl(adjustr(PictureFileName) // '.Y')
!       !Formatting the data that should be used by gnuplot
!       allocate(ToPlot(3,nx,nz))
!       do k=1,nz
!         ToPlot(2,:,k)=MeshX
!       enddo
!       do i=1,nx
!         ToPlot(1,i,:)=MeshZ
!       enddo
!       ! Choosing the correct density to plot and putting it into the correct 
!       ! formatting.
!       if(Isospin.eq.-1)then
!         do k=1,nz
!           do i=1,nx
!               ToPlot(3,i,k)=Density%Rho(i,Coord,k,1)
!           enddo
!         enddo    
!       elseif(Isospin.eq.0) then
!         do k=1,nz
!           do i=1,nx
!               ToPlot(3,i,k)=Density%Rho(i,Coord,k,1)+Density%Rho(i,Coord,k,2)
!           enddo
!         enddo    
!       else
!         do k=1,nz
!           do i=1,nx
!               ToPlot(3,i,k)=Density%Rho(i,Coord,k,2)
!           enddo
!         enddo            
!       endif
!       ! Writing the gnuplot data to a file.
!       call write_xyzgrid_data("DensityProfile", nx,nz, ToPlot, ierrorData)    
!     elseif(Direction.eq.3) then
!       xlabel='y (fm)'
!       ylabel='x (fm)'
!       PlotFileName=adjustl(adjustr(PictureFileName) // '.Z')
!       !Formatting the data that should be used by gnuplot
!       allocate(ToPlot(3,nx,ny))
!       do j=1,ny
!         ToPlot(2,:,j)=MeshX
!       enddo
!       do i=1,nx
!         ToPlot(1,i,:)=MeshY
!       enddo
!       ! Choosing the correct density to plot and putting it into the correct 
!       ! formatting.
!       if(Isospin.eq.-1)then
!         do j=1,ny
!           do i=1,nx
!               ToPlot(3,i,j)=Density%Rho(i,j,Coord,1)
!           enddo
!         enddo    
!       elseif(Isospin.eq.0) then
!         do j=1,ny
!           do i=1,nx
!               ToPlot(3,i,j)=Density%Rho(i,j,Coord,1)+Density%Rho(i,j,Coord,2)
!           enddo
!         enddo    
!       else
!         do j=1,ny
!           do i=1,nx
!               ToPlot(3,i,j)=Density%Rho(i,j,Coord,2)
!           enddo
!         enddo            
!       endif
!       ! Writing the gnuplot data to a file.
!       call write_xyzgrid_data("DensityProfile", nx,ny, ToPlot, ierrorData)
!     else
!       call stp('Wrong Direction in subroutine PlotDensity')    
!     endif
!     if(ierrorData.eq.0) then
!       ! Writing the gnuplot command script
!       call write_xyzgrid_contour ( "Command", "DensityProfile" ,               &
!       &                      trim(plotFilename), xlabel, ylabel, ierrorCommand )
!       if(ierrorCommand.eq.0) then
!           !Run Gnuplot with the command file
!           call run_gnuplot("Command")
!       else
!           call stp('Writing gnuplot command file was unsuccesful.')
!       endif
!     else
!       call stp('Writing datafile for gnuplot was unsuccesful.')
!     endif
!     deallocate(ToPlot)
!     return
!   end subroutine PlotDensity
  
!   subroutine Visualise
!   !-----------------------------------------------------------------------------
!   ! Subroutine for drawing some density profiles.
!   !-----------------------------------------------------------------------------
!     use Mesh
!     use SpwfStorage
!     use Densities
    
!     integer :: PlotX, PlotY, PlotZ
  
!     if(.not.allocated(DensityHistory)) then
!       call iniden
!     endif
    
!     if(.not.allocated(MeshX)) then
!       call inimesh
!     endif
    
!     call UpdateDensities(1,.true.)
    
!     if(SC) then
!       PlotX = 1
!     else
!       PlotX = nx/2+1
!     endif
!     if(TSC) then
!       PlotY = 1
!     else
!       PlotY = ny/2+1
!     endif
!     if(PC) then
!       PlotZ = 1
!     else
!       PlotZ = nz/2+1
!     endif

!     call PlotDensity(1,PlotX,0)
!     call sleep(2)
!     call PlotDensity(2,PlotY,0)
!     call sleep(2)
!     call PlotDensity(3,PlotZ,0)   
!   end subroutine Visualise

end module InOutput
