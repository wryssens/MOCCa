module Force
    !----------------------------------------------------------------------------------
    ! This module contains all the variables related to the Skyrme interaction used, 
    ! and that includes the physical constants, hbar, e^2 and others...
    !----------------------------------------------------------------------------------

    use Compilationinfo
    use GenInfo

    implicit none
    !-----------------------------------------------------------------------------
    ! Flag on how to calculate and print the energy of the Skyrme terms.
    ! Possible options at the moment:
    !  BTerms     => Calculate everything per coupling constant B
    !  TermByTerm => Calculate everything terms by term (DEFAULT)
    character(len=200), public   :: SkyrmeTreatment='TERMBYTERM'
    !-----------------------------------------------------------------------------
    !Name of the force
    character(len=200), public   :: afor
    !-----------------------------------------------------------------------------
    ! The diverse constants in the Skyrme force and their combinations (the Bi,Ci) 
    ! that enter the energy density.
    ! For the various definitions, see:
    !  P.Bonche, H.Flocard, P.H.Heenen; Nucl. Phys. A467 (1987); 115-135
    !-----------------------------------------------------------------------------
    !Skyrme Force parameters
    real(KIND=dp), public :: t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,te,to
    real(KIND=dp), public :: wso,wsoq
    !Functional parameters in the BFH representation
    real(KIND=dp), public :: B1,B2,B3,B4,B5,B6,B7a,B7b,B8a,B8b,Byt3a,Byt3b,B9,B9q
    real(KIND=dp), public :: B10,B11,B12a,B12b,B13a,B13b,B14
    real(KIND=dp), public :: B15,B16, B17, B18,B19,B20,B21
    !Functional parameters in the C-representation
    !Note: first component is the isoscalar constant, second the isovector part
    real(KIND=dp), public :: Crho(2),Crhosat(2), Ctau(2), Cdrho(2), CnablaJ(2)
    real(KIND=dp), public :: Cs(2), Cssat(2), Ct(2), Cf(2), Cds(2), Cnablas(2)
    !Recouplings of the tensor terms
    real(KIND=dp), public :: CJ0(2), CJ1(2), CJ2(2)

    ! J^2 terms and whether or not to take average nucleon masses
    logical               :: J2terms=.false., averagemass=.true.
    !.............................................................................
    !   Physical constants. They are set with an interaction, but have default  
    !   values.
    !    
    !   It concerns:                                                              
    !       - hbar in units of MeV * 10^{-22} s                                   
    !         builtin value:                                                      
    !             hbar   = 6.58211928                                             
    !       - m_n and m_p in units of MeV * c^{-2}                                
    !         builtin values:                                                     
    !             m_n    = 939.565379                                             
    !             m_p    = 938.272046                                             
    !       - speed of light in units of fm/(10^{-22} s)                          
    !         builtin value:                                                      
    !             c      = 29.9792458                                             
    !.............................................................................
    ! P.J. Mohr et al., Rev. Mod. Phys. 84, 2012, CODATA recommended values of     
    ! the fundamental physical constants constants:2010                           
    !                                                                             
    ! Or more practically on: http://physics.nist.gov/cuu/Constants/index.html    
    !.............................................................................                                           
    real(KIND=dp):: e2=1.43996446_dp, hbar = 6.58211928_dp, clum  =  29.9792458_dp
    real(KIND=dp):: nucleonmass(2) = (/939.565379_dp , 938.272046_dp /)
    real(KIND=dp):: hbm(2)      = 20.73551910_dp
    !-----------------------------------------------------------------------------
    ! Small remark on hbm: it is hbar^2/m and thus does not include the two that
    ! everyone is familiar with. This for historic reasons. But I hear you asking
    ! 'Wouter, then why is the default value here hbm/2?' and I'll answer that it
    ! is for input purposes that people that want to read a force can input 
    ! 20.something. In the forcereading routine, hbm is multiplied again by 2.
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! Integers governing COM-mass correction. They are intrinsically connected to 
    ! an interaction.
    !
    ! COMXBody = 
    !   0) X-body contribution is not used in the single-particle hamiltonian 
    !      and not printed in the energy summary. (And thus completely not 
    !      calculated.
    !   1) X-body contribution is not used in the single-particle hamiltonian,
    !      but calculated and added to the total energy (and printed ofcourse.)
    !   2) X-body contribution is used iteratively and ofcourse calculated, added
    !      to the total energy and printed!
    !-----------------------------------------------------------------------------
    integer, public              :: COM1Body=2, COM2Body=0

contains

    subroutine ReadForceInfo()
        !---------------------------------------------------------------------------
        ! Subroutine governing the user input of parameters regarding this module.
        !---------------------------------------------------------------------------
        NameList /Force/ afor, SkyrmeTreatment

        !Read the correct name from input
        read(unit=*, NML=Force)

        ! Make sure SkyrmeTreatment is aligned
        call to_upper(SkyrmeTreatment,Skyrmetreatment)

        ! Assigning the correct variables to the Force coefficients.    
        call ReadForce
        !Calculating all EDF coefficients
        call CalcEdfCoef
  end subroutine ReadForceInfo

  subroutine ResetConstants

    e2    = 1.43996446_dp
    hbm   = 20.73551910_dp
    nucleonmass = (/939.565379_dp , 938.272046_dp /)
    hbar = 6.58211928_dp
  end subroutine ResetConstants

  subroutine ReadForce
  !-----------------------------------------------------------------------------
  ! Subroutine that reads a force from the forces.param file that should be 
  ! provided.
  !-----------------------------------------------------------------------------

    integer            :: inunit, io
    character(len=200) :: Name, UpName, UpAfor
    logical            :: exists

    NameList /skf/ Name,t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,t3b,x3b,yt3b,te,to, &
    &                    wso,wsoq,                                         &
    &                    hbm,e2,                                           &
    &                    COM1body, COM2body,                               &
    &                    J2Terms   ,                                       &
    &                    averagemass, hbar, nucleonmass

    inquire(file='forces.param', exist=exists)
    if (.not. exists) call stp("MOCCa can't find forces.param file!")

    !Setting some default values
    t0=0.0_dp;      x0=0.0_dp
    t1=0.0_dp;      t2=0.0_dp;       x2=0.0_dp
    t3a=0.0_dp;     yt3a=0.0_dp;     t3b=0.0_dp;
    x3b=0.0_dp;     yt3b=0.0_dp
    te=0.0_dp;      to=0.0_dp
    wso=0.0_dp;     wsoq=0.0_dp 
    COM1body=0;     COM2body=0
    J2Terms=.false.
    call ResetConstants

    call get_unit(inunit)
    open(unit=inunit, file='forces.param',iostat=io)
    if (io.ne.0) call stp("Can't open forces.param!", 'iostat = ', io)
    do
        read(NML=skf, unit=inunit)
        !Creating uppercase strings
        call to_upper(afor, UpAfor)
        call to_upper(Name, UpName)
    
        if(trim(upname).eq.trim(upafor)) then
            exit
        elseif(trim(upname).eq.trim('END')) then
            call stp("Force not found on forces.param file!")
        else
             !Resetting some default values
            t0=0.0_dp;      x0=0.0_dp
            t1=0.0_dp;      t2=0.0_dp;       x2=0.0_dp
            t3a=0.0_dp;     yt3a=0.0_dp;     t3b=0.0_dp;
            x3a=0.0_dp
            x3b=0.0_dp;     yt3b=0.0_dp
            te=0.0_dp;      to=0.0_dp
            wso=0.0_dp;     wsoq=0.0_dp   
            b14=0.0_dp;     b15=0.0_dp;     b16=0.0_dp;      b17=0.0_dp
            COM1body=0;     COM2body=0
            J2Terms=.false.
            call ResetConstants
        endif
    enddo
    close(inunit)
    !Averaging masses if needed:
    if(Averagemass) nucleonmass = (nucleonmass(1) + nucleonmass(2))/2
!     !Multiply hbm by two
     hbm = hbm*2.0

  end subroutine ReadForce

  subroutine PrintForce
  !-----------------------------------------------------------------------------
  ! A subroutine that prints all the info related to the force (and the EDF) 
  ! that are used. Additionally, this prints some elementary constants that are 
  ! used.
  !-----------------------------------------------------------------------------          

    1       format(20("-"), 'Elementary Constants', 20("-"))
    3       format(22("-"), 'Skyrme Parameters', 21("-") )
    97      format("m_n            = ",f11.6, ' (MeV c^-2)',/ & 
    &              "m_p            = ",f11.6, ' (MeV c^-2)')
    98      format("e^2            = ", f11.8 ,' (sqrt(MeV fm)) ')  
    99      format("hbar^2/(2*mn)  = ", f11.8 ,' (MeV fm^2)'&
    &           ,/,"hbar^2/(2*mp)  = ", f11.8 ,' (MeV fm^2)')
    100     format("Name of the Force:  ", a57)
    101     format("t0 =", f12.3 , " x0  =", f12.3 , /, &
    &              't1 =', f12.3,  " x1  =", f12.3 , /, &
    &              "t2 =", f12.3 , " x2  =", f12.3 , / ,&
    &              "t3a=", f12.3 , " x3a =", f12.3 , " yt3a=", f12.3 ,/ ,&
    &              "t3b=", f12.3 , " x3b =", f12.3 , ' yt3b=', f12.3 ,/, &
    &              "te =", f12.3 , " to  =", f12.3 ,                  /, &
    &              "wso=", f12.3 , " wsoq=", f12.3)

    102     format('Force Options' , /  &
    &              '  COM1Body= ', i2,/ &
    &              '  COM2Body= ', i2,/ &
    &              '  J2Terms = ', l2)
    103     format(' One-body c.o.m. correction not included.')
    104     format(' One-body c.o.m. correction included perturbatively.')
    105     format(' One-body c.o.m. correction included selfconsistently.')
    106     format(' Two-body c.o.m. correction not included.')
    107     format(' Two-body c.o.m. correction included perturbatively.')
    108     format(' Two-body c.o.m. correction included selfconsistently.')

    110     format("EDF Coefficients (BFH representation)")
    111     format ( &
    &' B1 =',f13.6,'  B2 =',f13.6,'  B3 =',f13.6,/, &
    &' B4 =',f13.6,'  B5 =',f13.6,'  B6 =',f13.6,/, & 
    &' B7 =',f13.6,'  B8 =',f13.6,/,  &
    &' B9 =',f13.6,'  B9q=',f13.6,/, &
    &' B10=',f13.6,'  B11=',f13.6,'  B12=',f13.6,/, &
    &' B13=',f13.6,'  B14=',f13.6,'  B15=',f13.6,/, &
    &' B16=',f13.6,'  B17=',f13.6,'  B18=',f13.6/, &
    &' B19=',f13.6,'  B20=',f13.6,'  B21=',f13.6/,/)
    112 format ('EDF Coefficients (Isospin representation)')
    113 format (&
    &   ' C^rho_0[_0_] =',f13.6,' C^rho_1[_0_] =',f13.6,/, &
    &   ' C^rho_0[sat] =',f13.6,' C^rho_1[sat] =',f13.6,/, &
    &   ' C^s_0[_0_]   =',f13.6,' C^s_1[_0_]   =',f13.6,/, &
    &   ' C^s_0[sat]   =',f13.6,' C^s_1[sat]   =',f13.6,/, &
    &   ' C^tau_0      =',f13.6,' C^tau_1      =',f13.6,/, &
    &   ' C^Drho_0     =',f13.6,' C^Drho_1     =',f13.6,/, &
    &   ' C^divJ_0     =',f13.6,' C^divJ_1     =',f13.6,/, & 
    &   ' C^T_0        =',f13.6,' C^T_1        =',f13.6,/, &
    &   ' C^F_0        =',f13.6,' C^F_1        =',f13.6,/, &
    &   ' C^Ds_0       =',f13.6,' C^Ds_1       =',f13.6,/, &
    &   ' C^divs_0     =',f13.6,' C^divs_1     =',f13.6,/) 

    print 1
    print 98, e2
    print 97, nucleonmass
    print 99, hbm(1)/2, hbm(2)/2
    print *
    print 3
    print 100, adjustl(afor)
    print 101 , t0, x0, t1, x1, t2, x2, t3a, x3a, yt3a, t3b, x3b, yt3b, te, to,&
    &           wso, wsoq
    print *
    print 102, Com1Body, Com2body, J2Terms
    select case(COM1body)
    case(0)
        print 103
    case(1)
        print 104
    case(2)
        print 105
    end select

    select case(COM2body)
    case(0)
        print 106
    case(1)
        print 107
    case(2)
        print 108
    end select

    print *

    print 110
    print 111, B1,B2,B3,B4,B5,B6,B7a,B8a,B9,B9q,B10,B11,B12a,B13a,B14,B15,B16, &
    &          B17,B18,B19,B20,B21
    print 112
    print 113, Crho,Crhosat, Cs, Cssat, Ctau, Cdrho, CnablaJ, Ct, Cf, Cds,     &
    &          Cnablas

  end subroutine PrintForce

  subroutine CalcEDFCoef()
  !-----------------------------------------------------------------------------
  ! This subroutine calculates the Bi, the coefficients in the EDF, as a 
  ! function of the Skyrme parameters. See the formulas in the 24Mg paper and 
  ! the tensor notes.
  !-----------------------------------------------------------------------------
  ! Note for the calculation of the C's: we use the formulas in terms of B's, 
  ! because these work whether or  not we are dealing with a strict force.
  !-----------------------------------------------------------------------------
  ! Note again for the C's: they are BUGGED!
  !-----------------------------------------------------------------------------
    real(KIND=dp), parameter :: rhosat=0.16_dp

    !---------------------------------------------------------------------------
    !Calculating the B-s
    !---------------------------------------------------------------------------
    B1   = t0*(1.0_dp   + 1/2.0_dp*x0)/2.0_dp
    B2   =-t0*(1/2.0_dp + x0         )/2.0_dp

    B3   = (t1*( 1      + 1/2.0_dp*x1) + t2*(1      + 1/2.0_dp*x2))/4.0_dp
    B4   =-(t1*( 1/2.0_dp + x1       ) - t2*(1/2.0_dp +        x2))/4.0_dp

    B5   =-(3.0_dp*t1*(1  + x1/2.0_dp)-t2*(1  + x2/2.0_dp))/16.0_dp
    B6   = (3.0_dp*t1*(x1 +  1/2.0_dp)+t2*(x2 +  1/2.0_dp))/16.0_dp

    B7a  = t3a*( 1   + x3a/2.0_dp)/12.0_dp
    B8a  =-t3a*( x3a +   1/2.0_dp)/12.0_dp
    B7b  = t3b*( 1   + x3b/2.0_dp)/12.0_dp
    B8b  =-t3b*( x3b +   1/2.0_dp)/12.0_dp
    !---------------------------------------------------------------------------
    ! B9 and B9q are the constants related to the spin-orbit interaction. 
    ! Starting from a Skyrme Force, they should in principle be equal for a 
    ! true force, but apparently there are some physics cases where a difference
    ! between the isoscalar and isovector coupling is wanted.
    !---------------------------------------------------------------------------
    B9   = -wso /2.0_dp
    B9q  = -wsoq/2.0_dp

    Byt3a= yt3a
    Byt3b= yt3b

    if(.not.TRC) then
       B10  = t0*x0/4.0_dp
       B11  =-t0   /4.0_dp
       B12a = t3a*x3a/24.0_dp
       B13a =-t3a    /24.0_dp
       B12b = t3b*x3b/24.0_dp
       B13b =-t3b    /24.0_dp
    endif

    if(J2Terms) then
      B14  = -(1.0_dp/8.0_dp) * (t1*x1 + t2*x2)
      B15  =  (1.0_dp/8.0_dp) * (t1 - t2) 
      ! Tensor contribution
      B14 = B14 + (1.0_dp/4.0_dp) * (te + to)
      B15 = B15 - (1.0_dp/4.0_dp) * (te - to)
    endif
   
    B16  =-(3.0_dp/8.0_dp) * (te + to)
    B17  = (3.0_dp/8.0_dp) * (te - to)
    
!     if(J2Terms .and. (.not. TRC)) then
!       B18 =-(1.0_dp/32.0_dp)* (3 * t1 * x1 - t2 *x2)
!       B19 = (1.0_dp/32.0_dp)* (3 * t1      + t2    )
!       !Tensor contribution
!       B18 = B18 + (1.0_dp/16.0_dp) * (3*te - to)
!       B19 = B19 - (1.0_dp/16.0_dp) * (3*te + to)
!     endif
!     if(.not.TRC) then
!       B20 = (3.0_dp/16.0_dp) * (3*te - to) 
!       B21 =-(3.0_dp/16.0_dp) * (3*te + to)
!     endif

    !---------------------------------------------------------------------------
    !Calculating the C-s
    !---------------------------------------------------------------------------
    !Careful: Crho is density dependent in general. Crho contains the value at
    ! zero density, Crhosat the value at saturation density!
    !---------------------------------------------------------------------------
    Crho(1) = (3.0_dp / 8.0_dp) * t0   
    Crho(2) =-(1.0_dp / 4.0_dp) * t0 * (1.0_dp/2.0_dp + x0)

    Crhosat(1) = Crho(1)                              &
             & + 3.0_dp/48.0_dp * t3a * rhosat**byt3a &
             & + 3.0_dp/48.0_dp * t3b * rhosat**byt3b
    Crhosat(2) = Crho(2)                              &
             & - 1.0_dp/24.0_dp * t3a * (0.5_dp + x3a) * rhosat**byt3a & 
             & - 1.0_dp/24.0_dp * t3b * (0.5_dp + x3b) * rhosat**byt3b

    !Same remark for Cs: it is density dependent!
    Cs(1)   = B10 + 0.5_dp * B11
    Cs(2)=       0.5_dp * B11

    Cssat(1)= Cs(1) + (B10  + 0.5_dp*B11)                                      &
            &       + (B12a + 0.5_dp*B13a) * rhosat**byt3a &
            &       + (B12b + 0.5_dp*b13b) * rhosat**byt3b
    Cssat(2)= Cs(2) + 0.5_dp*B11                                               &
            &       + 0.5_dp*B13a  * rhosat**byt3a                             &
            &       + 0.5_dp*B13b  * rhosat**byt3b     

    Cds(1)  = B18 + 0.5_dp * B19
    Cds(2)  =       0.5_dp * B19 

    Ctau(1) = (3.0_dp/16.0_dp) * t1 + (1.0_dp/4.0_dp)*t2 * (5.0_dp/4.0_dp + x2)    
    Ctau(2) =-(1.0_dp/8.0_dp)  *(t1 * (0.5_dp + x1) - t2 * (0.5_dp + x2)) 

    Cdrho(1)=-(9.0_dp/64.0_dp) * t1 + 1.0_dp/16.0_dp * t2 * (5.0_dp/4.0_dp + x2) 
    Cdrho(2)= (1.0_dp/32.0_dp) * (3 * t1 *(0.5_dp + x1) + t2 * (0.5_dp + x2))

    CnablaJ(1) = B9 + 0.5_dp * B9q
    CNablaJ(2) =      0.5_dp * B9q

    Ct(1)   =-(B14 + 0.5_dp * B15)
    Ct(2)   =-       0.5_dp * B15

    Cf(1)   =-2.0_dp *(B16 + 0.5_dp*B17) ! additional factor -2 sign as the C
    Cf(2)   =-2.0_dp *(      0.5_dp*B17) ! refer to s*F, the b to J_ij J_ij

    Cnablas(1) = B20 + 0.5_dp * B21
    Cnablas(2) =       0.5_dp * B21

    CJ0    = -1.0_dp/3.0_dp * (Ct - 2.0_dp * Cf) 
    CJ1    = -0.5_dp        * (Ct - 0.5_dp * Cf)
    CJ2    = -                (Ct + 0.5_dp * Cf) 

    return
  end subroutine CalcEDFCoef
end module Force