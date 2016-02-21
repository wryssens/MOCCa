module InOutput

  use CompilationInfo
  use GenInfo
  use Force
  use Wavefunctions
  use Pairing
  use Transform
  use SpwfFactory
  use Interfaces
  
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
  real*8  :: filedx
  !-----------------------------------------------------------------------------
  !Determining if MOCCa needs to write extra information to the output file, 
  !for debugging purposes. Not active at the moment.
  logical :: ExtraOutput=.false.
  !-----------------------------------------------------------------------------
  ! Whether or not to write densities to file at the end of the iterations
  logical, public    :: Pictures=.false.
  !-----------------------------------------------------------------------------
  ! Temporary storage for HFB matrices, in order to be able to completely
  ! restart a calculation.
  ! The array FileBlocksizes handles the effective sizes of the matrices
  complex(KIND=dp), allocatable :: InputKappa(:,:,:,:), InputRho(:,:,:,:)
  complex(KIND=dp), allocatable :: InputU(:,:,:,:), InputV(:,:,:,:)
  complex(KIND=dp), allocatable :: InputCanTransfo(:,:,:,:)
  integer, allocatable          :: FileBlocksizes(:,:), FileHFBColumns(:,:,:)
  !-----------------------------------------------------------------------------
  ! Whether or not to write extra output files for various other codes.
  logical :: PromOutput = .false.
  !-----------------------------------------------------------------------------
  logical :: LegacyInput=.false.
  ! Convergence information
  real*8 :: fileE, filedE
  ! Force name that was used on the file
  character(len=200) :: fileforce=''
  ! Pairing information from file
  integer:: fileptype
  real*8 :: filegn,filegp
contains

  subroutine PrintFileInfo
    !---------------------------------------------------------------------------
    ! Subroutine that prints info to standardout that was found on a MOCCa
    ! file.
    !---------------------------------------------------------------------------

    1 format(20('-'), 'Information from file', 19('-'))
  100 format(60('-'))
    2 format('File:', a20)
    3 format('Mesh parameters:')
    4 format('  nx=',i5,'  ny=',i5,'  nz=',i5 )
    5 format('  dx=',f12.8, '(fm)')
    6 format(' nwt=',i5)
    7 format('Convergence:')
    8 format('   E=',f10.3)
    9 format('  dE=',e10.3)
   10 format('Skyrme parametrization:',3x, a20)
   11 format('Hartree-Fock')
   12 format('BCS pairing')
   13 format('HFB pairing')
   14 format('  Neutron pairingstrength:', f10.3)
   15 format('  Proton  pairingstrength:', f10.3)

    print 1
    print 2, trim(InputFilename)
    print 3 
    print 4, filenx,fileny,filenz
    print 5, filedx
    print 6, filenwt
    print *
    print 10, fileforce
    print *
    select case(fileptype)
    case (0)
      print 11
    case (1)
      print 12
      print 14, filegn
      print 15, filegp
    case (2)
      print 13
      print 14, filegn
      print 15, filegp
    end select
    print *
    print 7
    print 8, fileE
    print 9, filedE

  end subroutine PrintFileInfo

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
    ! Special mention here for the HFB matrices.
    if(allocated(InputKappa)) then
        if(LegacyInput) then
          KappaHFB = InputKappa
          if(InPC.neqv.PC) call stp("Don't break parity from legacy input.")
        else
          call TransformHFBMatrices(inputU, inputV, inputRho, InputKappa,        &
          &                               inPC,inIC,FileHFBColumns,FileBlocksizes)
        endif
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
    call writeMOCCa_v1(OutputChannel)

    ! Writing extra files if asked for
    if(PromOutput) call writePromesse(trim(OutputFileName)//'.prom')
    
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
            inTRC=.true.; inIC=.true.; inPC=.true.; inSC=.true.; inTSC=.true.
            call ReadNil(InputChannel, InputFilename)
    elseif(&
    &     (InputFileName(1:3).eq."ev4")  &
    & .or.(InputFileName(1:4).eq."int4") &
    & .or.(InputFileName(1:4).eq."d8a4") &
    & ) then
    !MOCCa is dealing with input from either EV4, D8A4 or INT4
      call ReadEv4(InputChannel,InputFilename)

      !All symmetries except parity are conserved in EV4
      inTRC=.true.; inIC=.true.; inPC=.false.; inSC=.true.; inTSC=.true.

    elseif(&
    &    ((InputFileName(1:3).eq."cr8"    ) &
    & .or.(InputFileName(1:4).eq."evcr"   ) &
    & .or.(InputFileName(1:7).eq."permcr8")) &
    & .and.(InputFileName(1:5).ne."cr8to")) then
      ! MOCCa is dealing with output from cr8, evcr8 or permcr8.
      ! (It can't be cr8to4 or cr8to1)
        call ReadCR8(InputChannel,InputFileName)

      !All symmetries except time reversal are conserved in CR8
      inTRC=.false.; inIC=.true.; inPC=.true.; inSC=.true.; inTSC=.true.

    elseif(&
    &      InputFilename(1:5).eq."MOCCa" ) then
        !MOCCa is dealing with its own type of input
        if(LegacyInput) then
          call ReadMOCCa_v0(InputChannel)
        else
          call ReadMOCCa_v1(InputChannel)
        endif
    else
      ! The input file is not clearly named; MOCCa defaults to its own 
      ! subroutine.
       if(LegacyInput) then
          call ReadMOCCa_v0(InputChannel)
        else
          call ReadMOCCa_v1(InputChannel)
        endif
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
    use Force, only         : ReadForceInfo
    use Densities, only     : ReadDensitInfo, Iniden
    use Moments, only       : ReadMomentData, SpecialInput
    use Cranking, only      : ReadCrankingInfo
    use Coulomb, only       : ReadCoulombInfo
    use SpwfStorage, only   : ReadSpwfStorageInfo
    use Mesh, only          : IniMesh
    use Cranking, only      : ReadCrankingInfo
    use Pairing, only       : ReadPairingInfo
    use Derivatives,only    : ReadDerivativesInfo
    use SpecialMoments,only : ReadSpecialMoments
    
    implicit none
    
    NameList /InAndOutput/ InputFileName,OutputFileName,ExtraOutput,Pictures, &
    &                      PromOutput,LegacyInput

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
    ! If signalled, read info on extra moments
    if(SpecialInput) call readSpecialMoments

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
  
  subroutine ReadMOCCa_v0(Ichan)
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
    integer       :: fileneutrons,fileprotons,i, version
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
    !---------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------
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
    ! 12)  Get an initial guess for the Fermi levels
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
      ! Only read omega from file, not the crankvalues. 
      ! Otherwise they will get overwritten
      read(IChan,iostat=ioerror) Omega !, CrankValues
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
  end subroutine ReadMOCCa_v0
  
  subroutine WriteMOCCa_v0(OChan)
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
    if(PairingType.ne.2) call GuessHFBMatrices(PairingType)
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
    
  end subroutine WriteMOCCa_v0

  subroutine WriteMOCCa_v1(OChan)
    !----------------------------------------------------------------------------
    ! Write info to an output file in the MOCCa format.
    !
    ! Current version of this routine :
    !   1 
    ! 
    ! Previous routine versions:
    !  Version            Difference
    ! ---------------------------------------------------------------------------
    !   0                 * Added version number
    !                     * Added storage of U and V HFB matrices
    !                     * Added storage of CanTransfo
    !----------------------------------------------------------------------------
    ! File format version 1
    ! ----------------------
    ! Version
    ! Convergence information: E, dE                         (*)
    ! nx,ny,nz,dx,dt
    ! TRC,PC,SC,TSC,IC
    ! nwt,neutrons,protons
    ! (nwt) Wavefunctions                                    => see Spwf module
    ! Densities                                              => See densities module
    ! MaxFDOrder, MaxFDLapOrder, CoulombLapOrder             (*)
    ! Forcename                                              (*)
    ! t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq   (*)
    ! afor,hbar,hbm,xm,njmunu,nmass,COM1Body,COM2Body        (*)
    ! Functional Coefficients                                (*)
    ! PairingType, PairingStrength, alpha, PairingCut        (*)
    ! Pairing Cutoffs                                        (*)
    ! FileBlocksizes
    ! U, V
    ! HFBColumns
    ! Rho, Kappa
    ! CanTransfo
    ! Omega, Crankvalues
    ! Moments                                                => See moments module
    ! ----------------------------------------------------------------------------
    ! * indicates that no action is done with this info at the moment
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
    use Energy

    integer, intent(in)   :: Ochan
    integer               :: i, io,P,it
    type(Moment), pointer :: Current

    open (OChan,form='unformatted',file=OutputFileName)
    ! File version
    write(Ochan, iostat=io) 1
    ! Convergence information
    write(Ochan,iostat=io) TotalEnergy, TotalEnergy-OldEnergy
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
    ! Force name
    write(Ochan,iostat=io) afor
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
    if(allocated(Pcutoffs)) then
      write(Ochan,iostat=io) Pcutoffs 
    else
      write(Ochan,iostat=io)
    endif
    write(OChan,iostat=io) Fermi, LNLambda
    ! No need to guess for the HFB matrices and variables if 
    ! PairingType = 'HFB' already.
    if(PairingType.ne.2) then
        ! Take a guess for Kappa, allocate the fileblocksizes and put
        ! other relevant matrices to zero.
        call GuessHFBMatrices(PairingType)
    endif

    write(OChan,iostat=io) blocksizes
    write(OChan,iostat=io) U,V
    write(Ochan,iostat=io) HFBColumns
    write(Ochan,iostat=io) RhoHFB,KappaHFB
    write(Ochan,iostat=io) CanTransfo

    !---------------------------------------------------------------------------
    !Cranking Variables
    write(OChan,iostat=io) Omega, CrankValues

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

end subroutine WriteMOCCa_v1

subroutine ReadMOCCa_v1(Ichan)
    !-----------------------------------------------------------------------------
    ! Read info from an input file in the MOCCa format.
    !
    ! Current version of this routine :
    !   1 
    !
    ! Previous routine versions:
    !  Version            Difference
    ! -------------------------------
    !   0                 * Added version number
    !                     * Added storage of U and V HFB matrices
    !                     * Added storage of CanTransfo
    !
    !----------------------------------------------------------------------------
    ! File format version 1
    ! ----------------------
    ! Version
    ! Convergence information: E, dE                         
    ! nx,ny,nz,dx,dt
    ! TRC,PC,SC,TSC,IC
    ! nwt,neutrons,protons
    ! (nwt) Wavefunctions                                    => see Spwf module
    ! Densities                                              => See densities module
    ! MaxFDOrder, MaxFDLapOrder, CoulombLapOrder             (*)
    ! Forcename                                              (*)
    ! t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq   (*)
    ! afor,hbar,hbm,xm,njmunu,nmass,COM1Body,COM2Body        (*)
    ! Functional Coefficients                                (*)
    ! PairingType, PairingStrength, alpha, PairingCut        (*)
    ! Pairing Cutoffs                                        (*)
    ! FileBlocksizes
    ! U, V
    ! HFBColumns
    ! Rho, Kappa
    ! CanTransfo
    ! Omega, Crankvalues
    ! Moments                                                => See moments module
    ! ----------
    ! * indicates that no action is done with this info at the moment
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
    integer :: ioerror, N, Pin,Iin, M(4),it,j,P

    !File parameters to compare against
    integer       :: fileneutrons,fileprotons,i, version
    real(KIND=dp) :: fileintensity(4)
    
    !---------------------------------------------------------------------------
    ! Checking for a valid input file
    inquire(file=inputfilename, exist=exists)
    if(.not.exists) call stp('Input file specified does not exist!')
    open (IChan,form='unformatted',file=InputFileName)

    read(IChan, iostat=ioerror) Version
    if(ioerror.ne.0) call stp('Input error for version number',                &
    &                         'ioerror', ioerror)
    if(Version.ne.1) call stp('File format has unrecognised version.',         &
    &                       'Version', Version)

    !Convergence info
    read(Ichan,iostat=ioerror) fileE, filedE

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
    !---------------------------------------------------------------------------
    ! b) Mesh parameters & number of wavefunctions must correspond
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
    if( nwt.gt.filenwt .and. TRC) then 
      call stp('Nwt is not compatible between file and data.',                 &
      &        'In data: ', nwt, 'On file :', filenwt )
    endif
    if( nwt.ne.2*filenwt .and. .not.TRC .and. inTRC) then 
      call stp('Nwt should be doubled when breaking TimeReversal.',            &
      &        'In data: ', nwt, 'On file :', filenwt )
    endif
    !---------------------------------------------------------------------------
    ! c) There should be enough space for neutrons & protons
    N = nwt
    if(TRC) N = 2*N    
    if( protons.gt.N) call stp('Not enough spwfs to accomodate all protons!' )
    if(neutrons.gt.N) call stp('Not enough spwfs to accomodate all neutrons!')
    !---------------------------------------------------------------------------
    ! 4) WaveFunctions
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
    ! 6 ) Specifics of the finite difference scheme
    !     Currently there is nothing of interest there.
    read(IChan, iostat=ioerror)
    !---------------------------------------------------------------------------
    ! 7 )  Force  : 
    !      a) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq
    !      b) afor,hbar,hbm,xm,njmunu,nmass,COM1Body,COM2Body
    !      c) Functional Coefficients
    !
    ! Force name
    read(Ichan,iostat=ioerror)  fileforce
    !Straight force parameters
    read(IChan,iostat=ioerror)
    !Extra parameters
    read(IChan,iostat=ioerror) 
    !Functional coefficients, not used at the moment.
    read(IChan,iostat=ioerror)
    !---------------------------------------------------------------------------
    ! 8) Pairing: type, strength, alphas and PairingCut
    !             Cutoffs
    !             Fermi, Lambda2
    read(IChan,iostat=ioerror) fileptype,filegn,filegp
    read(Ichan,iostat=ioerror)
    read(ICHan,iostat=ioerror)  Fermi, LNLambda
    if(ioerror.ne.0) then
        call stp('Did not read Fermi and LNLambda correctly' , &
             &   'Iostat', ioerror)
      endif
    !---------------------------------------------------------------------------
    ! 9) Reading U,V, HFBColumns, RhoHFB, KappaHFB and CanTransfo (only if needed)
    select case(PairingType)

    case(2)
        ! First decide on the appropriate sizes for U,V, Rho, Kappa and CanTransfo
        M(3) = 1 ; if(inPC) M(3) = 2
        M(4) = 1 ; if(InIc) M(4) = 2
        
        allocate(FileBlocksizes(M(3), M(4))) ; FileBlocksizes = 0

        M(1) = filenwt ; M(2) = filenwt
        if(InTRC) then
            M(1) = 2 * M(1) ; M(2) = 2 * M(2)
        endif

        allocate(InputKappa     (M(1),M(2),M(3),M(4))) 
        allocate(InputRho       (M(1),M(2),M(3),M(4)))
        allocate(InputU         (M(1),2*M(2),M(3),M(4))) 
        allocate(InputV         (M(1),2*M(2),M(3),M(4)))
        allocate(InputCanTransfo(M(1),M(2),M(3),M(4)))
        allocate(FileHFBColumns (M(1),M(3),M(4)))

        read(IChan,iostat=ioerror) Fileblocksizes
        if(ioerror.ne.0) then
            call stp('Did not read HFB blocksizes correctly','Iostat', ioerror)
        endif
        
        read(ICHan,iostat=ioerror) InputU,InputV 
        if(ioerror.ne.0) then
            call stp('Did not read U and V correctly','Iostat', ioerror)
        endif

        read(Ichan,iostat=ioerror) FileHFBColumns
        if(ioerror.ne.0) then
            call stp('Did not read HFBColumns correctly from file.')
        endif

        read(ICHan,iostat=ioerror) InputRho,InputKappa     
        if(ioerror.ne.0) then
            call stp('Did not read Rho and Kappa correctly' ,'Iostat', ioerror)
        endif
        read(ICHan,iostat=ioerror) InputCanTransfo
        if(ioerror.ne.0) then
            call stp('Did not read Rho and Kappa correctly' ,'Iostat', ioerror)
        endif
    case DEFAULT
        !--------------------------------------
        ! HF & BCS case. Don't read anything at all.
        read(Ichan,iostat=ioerror)
        read(Ichan,iostat=ioerror)
        read(Ichan,iostat=ioerror)
        read(Ichan,iostat=ioerror)
        read(Ichan,iostat=ioerror)
    end select
    !---------------------------------------------------------------------------
    ! 11) Cranking parameters: effective, true and intensity
    if(ContinueCrank) then
        ! Don't read other values, as we don't want to override them.
        read(IChan,iostat=ioerror) Omega
    else
        read(Ichan,iostat=ioerror)
    endif

    !---------------------------------------------------------------------------
    ! 10) Moment parameters
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

    close(Ichan)
end subroutine ReadMOCCa_v1

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
end module InOutput
