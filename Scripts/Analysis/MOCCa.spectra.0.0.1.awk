#==============================================================================#
#                                                                              #
#    MOCCa.spectra.awk     release 0.0.1                                       #
#    8 November 2016                                                           #
#                                                                              #
#    Copyright W. Ryssens & M. Bender                                          #
#    Heavily inspired on cr8.spectra.1.1.2.awk by M. Bender and engineered     #                                                                              #
#    to produce the same kind of tables, only starting from MOCCa files.       #
#    Note that when I say 'kind', I do mean 'kind', and not 'identical'.       #
#                                                                              #
#    Note that Spwf.sort.awk is also required in order for this script to      #
#    function.                                                                 # 
#                                                                              #
#==============================================================================#
# TODO                                                                         #
# + implement quasiparticle basis reading & writing                            #
# + Add detection of quadrupole orientation                                    #
# + Add blocking information                                                   #
# + Add convergence info                                                       #
#==============================================================================#
#                                                                              #
#    Usage                                                                     #
#                                                                              #
#    awk -f MOCCa.spectra.awk  PREFIX.out                                      #
#                                                                              #
#    where PREFIX.out is a collection of MOCCa outputs, arranged in the desired#
#    order. (For instance by cat'ing the different files one after another).   #
#                                                                              #
#                                                                              #
#    Note that it is VERY important that the different MOCCa runs match in type#
#    DO NOT TRY TO MIX with one awk call:                                      #
#       * Calculations with different pairing options                          #
#       * Calculations with different number of SPWFs                          #
#-------------------------------------------------------------------------------
# This script generates the following files (see headers for detailed content) #
#                                                                              #
# - PREFIX.e.tab      total energies and energy of tensor terms                #
# - PREFIX.ef.tab     Fermi energies and quantities related to the LN scheme   #
# - PREFIX.[iso].qlm.tab    Multipole moments <Q_lm>  and related \beta_lm for #
#                           neutrons (n), protons(p) and total(t).             #
# - PREFIX.[iso].ql.tab     Total multipole moments <Q_l>  and related \beta_l #
#                           for neutrons (n), protons(p) and total(t).         #
#                                                                              #
# - PREFIX.[base].[iso].par=[par].sig=[sig].tab                                #   
#                                                                              #
#  where                                                                       #
#                                                                              #
#  PREFIX     name of the file read without the suffix ".out"                  #
#                                                                              #
#  [base] = hf "Hartree-Fock" basis diagonalizing the sp Hamiltonian           #
#         = ca canonical      basis diagonalizing the one-body density matrix  # 
#         = qp  quasiparticle  basis diagonalizing the q-p Hamiltonian         # 
#          (NOT IMPLEMENTED)                                                   #
#  [iso]  = n  neutrons                                                        #
#         = p  protons                                                         #
#  [par]  = +1  positive parity                                                #
#         = -1  negative parity                                                #
#         =  0  parity broken                                                  #
#  [sig]  = +1  positive signature                                             #    
#         = -1  negative signatiure                                            #
#         =  0  signature broken                                               #
#                                                                              #
#  Note that files not appropriate to the MOCCa runs will not be generated.    #
#  Meaning that                                                                #
#  - Time-reversal conserving calculations will not produce [sig]=-1 files     #
#  - Parity/signature broken calculations will only produce [par]/[sig]=0 files#
#  - Canonical and quasiparticle basis files will only be created with HFB runs#
#===============================================================================

function SortSpwfs (file, iso, basis, prefix, PC, SC, TRC, points){
    #---------------------------------------------------------------------------
    # Sort a file of SPWF info according to the symmetries of the calculation.
    # We use the extra awk script Spwf.sort.awk to sort on a specific column.
    #===========================================================================
        
    if (PC == 1) {
        # Sort the spwfs according to parity
        command = "awk -f Spwf.sort.awk 'column=4' 'points=" iqmax "' < " file;
        system(command);
        #exit
        # Moving 
        command = " mv tmp.p spwf.par=+1"
        system(command)
        
        command = " mv tmp.m spwf.par=-1"
        system(command)
    }
    else{
        command = "awk -f Spwf.sort.awk 'column=-1' 'points=" iqmax "' < " file;
        system(command)
        command = " mv tmp.zero spwf.par=0"
        system(command)
    }
    
    # Remove the original file
    system("rm " file)
    
    if (TRC == 0 && SC == 1) {
        # only distinguish between positive and negative signature when 
        # time-reversal is not conserved.
        if (PC == 1) {
            # Sort the spwfs according to signature
            # Positive spwf.parity
            file = "spwf.par=+1"
            command = "awk -f Spwf.sort.awk 'column=5' 'points=" iqmax "' < " file;
            system(command);
            
            command = " mv tmp.p spwf.par=+1.sig=+1"
            system(command)
            command = " mv tmp.m spwf.par=+1.sig=-1"
            system(command)
            
            #Negative spwf.parity
            file = "spwf.par=-1"
            command = "awk -f Spwf.sort.awk 'column=5' 'points=" iqmax "' < " file;
            system(command);
            
            command = " mv tmp.p spwf.par=-1.sig=+1"
            system(command)
            command = " mv tmp.m spwf.par=-1.sig=-1"
            system(command)
        }
        else {
            file = "spwf.neutron.par=0"
            command = "awk -f Spwf.sort.awk 'column=5' 'points=" iqmax "' < " file;
            system(command);
            
            command = " mv tmp.p spwf.par=0.sig=+1"
            system(command)
            command = " mv tmp.m spwf.par=0.sig=-1"
            system(command)
        }
    }
    else if ( TRC == 1 && SC == 1) {
        #Only signature = +1 states are printed
       if(PC == 1) {
            command =  "mv spwf.par=+1 spwf.par=+1.sig=+1"
            system(command)
            command =  "mv spwf.par=-1 spwf.par=-1.sig=+1"
            system(command)
       }
       else{
            command =  "mv spwf.par=0 spwf.par=0.sig=+1"
            system(command)
       }
    }
    else if ( SC == 0) {
        #No signature sorting
        if(PC == 1) {
            command = "awk -f Spwf.sort.awk 'column=-1' 'points=" iqmax "' < spwf.par=+1";
            system(command)
            command = " mv tmp.p spwf.par=+1.sig=0"
            system(command)
            
            command = "awk -f Spwf.sort.awk 'column=-1' 'points=" iqmax "' < spwf.par=-1";
            system(command)
            command = " mv tmp.p spwf.par=-1.sig=0"
            system(command)
    
            command="rm spwf.par=-1"
            system(command)
            command="rm spwf.par=+1"
            system(command)
          
        }
        else{
            command = "awk -f Spwf.sort.awk 'column=-1' 'points=" iqmax "' < spwf.par=0 ";
            system(command)
            command = " mv tmp.zero spwf.par=0.sig=0"
            system(command)
        }         
    }
    
    # Clearly name all of the files
    # Do this via a temporary script in order to get the correct bash shell. 
    # In a more naive way on my ubuntu installation it uses the default /bin/sh 
    # which is dash and not bash, and does not accept substring substitutions.
    command = "for f in spwf* \n do \n mv $f " prefix "." basis "." iso ".${f/spwf./}.tab \n done" 
    printf("#!/bin/bash \n") > "script.sh"
    printf(command ) >> "script.sh"
    close("script.sh")
    system("bash script.sh")
    system("rm script.sh")
}

function ReadSpwf(basis, PairingType, SC, TSC, PC, TRC, sorted)
{
    #===========================================================================
    # Depending on the pairingtype and symmetries, correctly order the Spwf info
    #===========================================================================
    
    Eind=0
    #Spenergy
    if(PairingType == "HFB" && basis =="hf")
    {
        Eind=8
    }
    if(PairingType == "HFB" && basis =="can")
    {
        Eind=5
    }
    else if (PairingType=="BCS")
    {
        Eind=6
    }
    else if (PairingType =="HF")
    {
        Eind=5
    }
                  
    #index of the state
    sorted[1] = $2
        
    # Expected value of parity
    sorted[2] = $3
    
    if( TRC != 1.0 || SC != 1.0 ){
        # Expected value of signature
        sorted[3] = $4
        #Occupation number or Rho_ii
        sorted[4] = $5
    }
    else {
        #Automatically positive signature
        sorted[3] =  1
        sorted[4] = $4
    }

    sorted[5] = $Eind
    
    # Angular momentum quantum numbers Jx, Jy, Jz and J as well as <r^2>
    k = 1
    while (k < 6){
        sorted[5+k] = $(Eind + k +1 )
        k+=1
    }

    return 
}

BEGIN{
        iq = 0
        flag[multipole]=0
        energyflag     =0
        fdenergyflag   =0
        pairingflag    =0
        angmomflag     =0
        
        HFBasisflag    =0
        CanBasisflag   =0
        
        #Only detect the first mention of HF(B)
        pairingtypeflag=1
    
}

{
        #-----------------------------------------------------------------------
        # Start to read the file with the datasets.
        # ----------------------------------------------------------------------
        # 
        # Determine filename using awk builtins 
        # (NR = number of records processed)
        # (FILENAME) = filename
        if ( NR == 1 ) {
                lengthfilename = length(FILENAME);
                prefix = substr(FILENAME,1,lengthfilename-4);
        }
        
        #  Take some parameters from the start of the calculations. 
        #  a) number of protons & neutrons
        if ( $1 == "Nucleus") {
                getline;
                neutrons=$3
                protons =$6
                mass    =neutrons+protons
                R       =1.2 * xa^(1.0/3.0)
                #print "(A,N,Z)", mass,neutrons,protons       
        }
        #  b) number of neutron & proton wave-functions
        if ( $1 == "Wavefunctions") {
                getline;
                nwt = $3
                getline;
                nwn = $3
                getline;
                nwp = $3
                #print "nwt,nwn,nwp", nwt,nwn,nwp
        }
        #   c) symmetries 
        if ( $1 == "Symmetries") {
               PC = 1
               SC = 1
               TRC= 1
               TSC= 1
               getline;
               if ($3 != "Conserved") {PC = 0}
               getline;
               if ($3 != "Conserved") {SC = 0}
               getline;
               if ($2 != "Conserved") {TRC = 0}
               getline;
               if ($3 != "Conserved") {TSC = 0}
               #print "P,Rz,TR,TS", PC, SC, TRC, TSC
        }
        #   d) number of multipole moments to consider
        if ( $1 == "Maximum" && $2 == "l" ) {
                maxmultipole=$4     
        }
        #   e) type of pairing to deal with
        if ( (pairingtypeflag == 1)  && ($1 == "BCS" || $1 == "HFB" || $1 == "Hartree-Fock;")) {
            PairingType = $1
            pairingtypeflag = 0
            if(PairingType == "Hartree-Fock;" ){
                PairingType = "HF"
            }
        }
        #   f) dealing with one-body correction or not?
        if ($1 == "COM1Body=") {
            if($2 == 0) {
                COM1flag = 0
            }
            else{
                COM1flag = 1
            }
        }
       
        #------------------------------------------------------------------------
        # This flag marks the end of a MOCCa calculation in the file
        #------------------------------------------------------------------------
        if ( $1 == "**FINAL**" ) {
            # Increment calculation counter
            iq   =  iq + 1
            iqmax=  iqmax + 1
            #Readflags
            flag[multipole]=1
            energyflag     =1
            fdenergyflag   =1
            pairingflag    =1
            if(PairingType == "HF") {
                HFBasisflag =1
                CanBasisflag=0
            }
            if(PairingType == "BCS") {
                HFBasisflag =1
                CanBasisflag=0
            }
            if(PairingType == "HFB") {
                HFBasisflag =1
                CanBasisflag=1
            }
            angmomflag=1
        }
        #-----------------------------------------------------------------------
        # Reading the SPWF info
        if( $2 == "Sp" && (CanBasisflag || HFBasisflag) ){
            getline;
            getline;
            getline;
            getline;
            getline;
            getline;
            getline;
            if(PairingType == "HFB") {
             getline;
             getline;
            }
            #----------------------------------
            # Neutron states in the HF basis
            N = 1
            if(HFBasisflag == 1){
                while( NF != 1){
                     ReadSpwf("hf", PairingType, SC, TSC, PC, TRC, line)
                     i = 1
                     while( i < NF +1 ){
                        neutronhf[iq,N,i] = line[i]
                        i+=1
                     }
                     
                    N +=1
                    getline;
                }
            }
            #----------------------------------
            # Neutron states in the canonical basis
            N = 1
            if(CanBasisflag == 1){
                getline;
                getline;
                getline;
                getline;
                getline;
                getline
                while(  NF != 1){
                    i=1 
                    ReadSpwf("can", PairingType, SC, TSC, PC, TRC, line)
                    while ( i<NF+1 ) { 
                        neutroncan[iq, N ,i]= line[i]
                        i+=1 
                    }
                    N +=1
                    getline;
                }
            }
            getline;
            getline;
            getline;
            getline;
            getline;
            getline;
            getline;
            getline;

            #----------------------------------
            # Proton states in the HF basis
            P = 1
            if(HFBasisflag == 1){
                while( NF != 1){
                    i=1
                    ReadSpwf("hf", PairingType, SC, TSC, PC, TRC, line)
                    while ( i<NF+1 ) { 
                        protonhf[iq, P ,i]= line[i]
                        i+=1 
                    }
                    P +=1
                    getline;
                }
            }
            HFBasisflag = 0
            
            #----------------------------------
            # proton states in the canonical basis
            P = 1
            if(CanBasisflag == 1){
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                while(  NF != 1){
                    i=1 
                    ReadSpwf("can", PairingType, SC, TSC, PC, TRC, line)

                    while ( i<NF+1 ) { 
                        protoncan[iq, P ,i]= line[i]
                        i+=1 
                    }
                    P +=1
                    getline;
                }
            }
            CanBasisflag = 0
        }
                   
        
        #-------------------------------------------------------------------
        # read information on multipole moments
        if ( flag[multipole] && $2 == "Multipole" && $3 == "Moments")  {
                getline; # Quantisation axis line
                getline; # Secondary axis ordering line
                getline; # Header lines
                getline; #
                getline; #  
                getline; # Number of particles
                getline; # Empty 
                getline;
                
                # Q_{lm} values
                while($1 != "RMS") {
                        if(NF != 0) {
                                ReIm = $1
                                l    = $3
                                m    = substr($4, 1, length($4)-1)
                                QN   =$5
                                QP   =$6
                                QT   =$7
                                
                                Qlm[iq,ReIm,l,m,1] = QN
                                Qlm[iq,ReIm,l,m,2] = QP
                                Qlm[iq,ReIm,l,m,3] = QT
                                
                                #printf("%03f \n",Qlm[iq,ReIm,l,m,3])
                        }
                        getline;
                }
                getline; #separator
                getline; #header
                getline;
                getline;
                # Beta_{lm} values
                while(NF != 1) {
                        if(NF != 0) {
                                #Note that MOCCa does not correctly print the
                                #beta deformation of imaginary parts of
                                #multipole moments yet
                                ReIm = "Re"
                                l  = $2
                                m  = substr($3, 1, length($3)-1)
                                BN =$4
                                BP =$5
                                BT =$6
                                
                                Beta[iq,ReIm,l,m,1] = BN
                                Beta[iq,ReIm,l,m,2] = BP
                                Beta[iq,ReIm,l,m,3] = BT
                        }
                        getline;
                }
                
                getline; #header
                getline;
                getline;

                # Ql values
                while(NF != 1) {
                        if(NF != 0) {
                                l  = substr($2, 1, length($2)-1)
                                QN =$3
                                QP =$4
                                QT =$5
                                
                                Ql[iq,l,1] = QN
                                Ql[iq,l,2] = QP
                                Ql[iq,l,3] = QT
                                
                                #printf("%03f \n",Ql[iq,l,1])
                                                               
                        }
                        getline;
                }
                getline; #header
                getline;
                getline;
                # Betal values
                while(NF != 1) {
                        if(NF != 0) {
                                l  = substr($2, 1, length($2)-1)
                                BN =$3
                                BP =$4
                                BT =$5
                                
                                Betal[iq,l,1] = BN
                                Betal[iq,l,2] = BP
                                Betal[iq,l,3] = BT
                        }
                        getline;
                }
                getline; #header
                                
                # Alternative quadrupole representations
                getline; #Cartesian header
                                                             
                QuadCart[iq,X,1] = $3
                QuadCart[iq,Y,1] = $4
                QuadCart[iq,Z,1] = $5
                getline;
                QuadCart[iq,X,2] = $3
                QuadCart[iq,Y,2] = $4
                QuadCart[iq,Z,2] = $5
                getline;
                QuadCart[iq,X,3] = $3
                QuadCart[iq,Y,3] = $4
                QuadCart[iq,Z,3] = $5
                
                #print QuadCart[iq,X,1],QuadCart[iq,X,2],QuadCart[iq,X,3]
                
                # (Q,Gamma representation)
                Q0[iq,1]     = $3
                Gamma[iq,1]  = $4
                getline;
                Q0[iq,2]     = $3
                Gamma[iq,2]  = $4
                getline;
                Q0[iq,3]     = $3
                Gamma[iq,3]  = $4

                # (iq1,iq2 representation)
                getline;
                iq1[iq,1]  = $3
                iq2[iq,1]  = $4
                getline;
                iq1[iq,2]  = $3
                iq2[iq,2]  = $4
                getline;
                iq1[iq,3]  = $3
                iq2[iq,3]  = $4
                #Reset flag
                flag[multipole]=0
        }
        
        if ( angmomflag && $2 == "Angular") {
                getline;
                getline;
                getline;
                getline;
                
                if( SC ==  0) {
                    Jx[iq]     = $3
                    OmegaX[iq] = $5
                    getline;
                }
                else{
                    Jx[iq]     = 0
                    OmegaX[iq] = 0
                }
                
                if( TSC == 0) {
                    Jy[iq]     = $3
                    OmegaY[iq] = $5
                    getline;
                }
                else{
                    Jy[iq]     = 0
                    OmegaY[iq] = 0
                }
                
                Jz[iq]     = $3
                OmegaZ[iq] = $5          
                angmomflag=0     
        }

        if ( pairingflag && $2 == "Pairing") {
                getline;
                getline;
                Fermi[iq,1] = $4
                Fermi[iq,2] = $5
                getline;
                getline;
                Dispersion[iq,1] = $2
                Dispersion[iq,2] = $3
                getline;
                #Act if next line is LN-enabled
                if($1 == "Lambda_2"){
                    Lambda2[iq,1] = $2
                    Lambda2[iq,2] = $3
                    getline;
                }
                pairingflag=0
        }
        if ( fdenergyflag && $2 == "Energies")  {
                # The finite difference energy
                getline; 
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline; #Kinetic energy is this line
                getline;
                getline;
                
                if( PairingType != "HF") {
                    getline;                    
                }   

                #Check if we are dealing with a LN calculation or not
                if($1 =="Lipkin-Nogami") {
                        getline ;
                        getline ;
                        getline ;
                }
                else{
                        getline ;
                }
                
                if(COM1flag == 1) {
                    getline;
                    getline;
                }
                
                getline;
                
                getline;
                Energy[iq,2] = $3 #Finite difference energy
                getline;
                getline;
                getline;
                getline;
                fdenergyflag = 0
        }
                
        if ( energyflag && $2 == "Lagrange")  {
                # Get all the things related to the energy
                getline; 
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                getline;
                
                Kinetic[iq,1] = $3
                Kinetic[iq,2] = $6
                Kinetic[iq,3] = $8
                
                getline;
                getline;
                
                if( PairingType != "HF") {
                    Pairing[iq,1] = $3
                    Pairing[iq,2] = $6
                    Pairing[iq,3] = $8
                    getline;                    
                }   

                #Check if we are dealing with a LN calculation or not
                if($1 =="Lipkin-Nogami") {
                        getline ;
                        ELN[iq,1] = $3
                        ELN[iq,2] = $6
                        ELN[iq,3] = $8
                        getline ;
                        getline ;
                }
                else{
                        getline ;
                        ELN[iq,1] = 0.0
                        ELN[iq,2] = 0.0
                        ELN[iq,3] = 0.0
                }
                Coulomb[iq,direct] = $2
                Coulomb[iq,exch]   = $4

                if(COM1flag == 1) {
                    getline;
                    getline;
                    
                    COM1[iq,1] = $3
                    COM1[iq,2] = $6
                    COM1[iq,3] = $8
                    
                }
                
                getline;
                getline;
                getline;
                
                Energy[iq,1] = $2 #Lagrange energy      
                getline;
                Energy[iq,2] = $2 #Finite difference energy
                getline;
                Energy[iq,3] = $3 #SPWF energy
                getline;
                getline;
                getline;
                Energy[iq,4] = $2 #Routhian
                energyflag = 0
        }        
}
#
#-------------------------------------------------------------------------------
#    All MOCCa files have been read. 
#    Now the collected data can be printed into tables
#
END{
  nn = savenn;
  np = savenp;
    
    #---------------------------------------------------------------------------
    # total energy
    # and closely related observables
    #if ( calc == "pes" ) {
        print "!       E(func)     E_FD     E(sp)  E(func)noLN      eLN(n)     eLN(p)      OmegaX      Jx        OmegaY      Jy        OmegaZ      Jz" > "tmp.e.tab";
    #}
    iq=1;
    while ( iq < iqmax + 1 ) {
        enoln = Energy[iq,1] - ELN[iq,3];
#        if ( calc == "pes" ) {
            printf("%3.0f %10.3f %10.3f %10.3f   %10.3f %10.3f   %10.6f %8.3f %12.5f %8.3f %12.5f %8.3f %12.5f \n",
               iq,Energy[iq,1],Energy[iq,2],Energy[iq,3],enoln,ELN[iq,1],ELN[iq,2],OmegaX[iq],Jx[iq],OmegaY[iq],Jy[iq],OmegaZ[iq],Jz[iq]) >> "tmp.e.tab"; 
#        }
        iq += 1;
    }
    close("tmp.e.tab");
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    # Deformations of all the Re/Im Qlm detected for neutrons, protons and total 
    iq = 1;
#    if ( calc == "pes" ) {
        printf("!    ") > "tmp.n.qlm.tab";
        printf("!    ") > "tmp.p.qlm.tab";
        printf("!    ") > "tmp.t.qlm.tab";
        #Make the header by looking which Qlm entries are initialised
#        if ( calc == "pes" ) {
            l = 1
            while ( l < maxmultipole +1) {
                m = 0
                while ( m < l+1 ) {
                    if( Qlm[iq,"Re",l,m,3] != "" ) {
                        printf("   ReQ%1.0f%1.0f       Beta%1.0f%1.0f  ",l,m,l,m)> "tmp.n.qlm.tab";
                        printf("   ReQ%1.0f%1.0f       Beta%1.0f%1.0f  ",l,m,l,m)> "tmp.p.qlm.tab";
                        printf("   ReQ%1.0f%1.0f       Beta%1.0f%1.0f  ",l,m,l,m)> "tmp.t.qlm.tab";
                    }
                    m+=1
                }
                l+=1 
            }
            printf("\n")>>"tmp.n.qlm.tab"
            printf("\n")>>"tmp.p.qlm.tab"  
            printf("\n")>>"tmp.t.qlm.tab"     
#        }
#    }
    while ( iq < iqmax + 1 ) {
#        if ( calc == "pes" ) {
            printf("%3.0f", iq) >> "tmp.n.qlm.tab"
            printf("%3.0f", iq) >> "tmp.p.qlm.tab" 
            printf("%3.0f", iq) >> "tmp.t.qlm.tab"  
            l = 1
            while ( l < maxmultipole +1) {
                m = 0
                while ( m < l+1 ) {
                    if( Qlm[iq,"Re",l,m,3] != "" ) {
                        printf("%12.3e %10.3f", Qlm[iq,"Re",l,m,1],Beta[iq,"Re",l,m,1])>> "tmp.n.qlm.tab"
                        printf("%12.3e %10.3f", Qlm[iq,"Re",l,m,2],Beta[iq,"Re",l,m,2])>> "tmp.p.qlm.tab"
                        printf("%12.3e %10.3f", Qlm[iq,"Re",l,m,3],Beta[iq,"Re",l,m,3])>> "tmp.t.qlm.tab"
                    }
                    if( Qlm[iq,"Im",l,m,3] != "" ) {
                        printf("%12.3e %10.3f", Qlm[iq,"Im",l,m,1],Beta[iq,"Re",l,m,1])>> "tmp.n.qlm.tab"
                        printf("%12.3e %10.3f", Qlm[iq,"Im",l,m,2],Beta[iq,"Re",l,m,2])>> "tmp.p.qlm.tab"
                        printf("%12.3e %10.3f", Qlm[iq,"Im",l,m,3],Beta[iq,"Re",l,m,3])>> "tmp.t.qlm.tab"
                    }
                    m +=1
                }
                l+=1 
            }
            printf("\n")>>"tmp.n.qlm.tab"
            printf("\n")>>"tmp.p.qlm.tab"  
            printf("\n")>>"tmp.t.qlm.tab"   
#        }    
        iq += 1;
    }
    close("tmp.n.qlm.tab");
    close("tmp.p.qlm.tab");
    close("tmp.t.qlm.tab");
    #---------------------------------------------------------------------------
    # Deformations of all the Ql detected for neutrons, protons and total 
    iq = 1;
#    if ( calc == "pes" ) {
        printf("!     ") > "tmp.n.ql.tab";
        printf("!     ") > "tmp.p.ql.tab";
        printf("!     ") > "tmp.t.ql.tab";
        #Make the header by looking which Qlm entries are initialised
#        if ( calc == "pes" ) {
            l = 1
            while ( l < maxmultipole +1) {
                
                if( Ql[iq,l,3] != "" ) {
                    printf("  ReQ%1.0f         Beta%1.0f   ",l,l)> "tmp.n.ql.tab";
                    printf("  ReQ%1.0f         Beta%1.0f   ",l,l)> "tmp.p.ql.tab";
                    printf("  ReQ%1.0f         Beta%1.0f   ",l,l)> "tmp.t.ql.tab";
                }
                l+=1 
            }
            printf("\n")>>"tmp.n.ql.tab"
            printf("\n")>>"tmp.p.ql.tab"  
            printf("\n")>>"tmp.t.ql.tab"     
#        }
#    }
    while ( iq < iqmax + 1 ) {
#        if ( calc == "pes" ) {
            printf("%3.0f", iq) >> "tmp.n.ql.tab"
            printf("%3.0f", iq) >> "tmp.p.ql.tab" 
            printf("%3.0f", iq) >> "tmp.t.ql.tab"  
            l = 1
            while ( l < maxmultipole +1) {
                if( Ql[iq,l,3] != "" ) {
                    printf("%12.3e %10.3f", Ql[iq,l,1],Betal[iq,l,1])>> "tmp.n.ql.tab"
                    printf("%12.3e %10.3f", Ql[iq,l,2],Betal[iq,l,2])>> "tmp.p.ql.tab"
                    printf("%12.3e %10.3f", Ql[iq,l,3],Betal[iq,l,3])>> "tmp.t.ql.tab"
                }
                l+=1 
            }
            
            printf("\n")>>"tmp.n.ql.tab"
            printf("\n")>>"tmp.p.ql.tab"  
            printf("\n")>>"tmp.t.ql.tab"   
#        }    
        iq += 1;
    }
    close("tmp.n.ql.tab");
    close("tmp.p.ql.tab");
    close("tmp.t.ql.tab");
    #---------------------------------------------------------------------------
    #------------ Fermi energies and information on Lipkin-Nogami scheme--------
#    if ( calc == "pes" ) {
        print "!      eF_n     eF_p   lambda_2n lambda_2p  eLN(n)   eLN(p) <DeltaN^2> <DeltaZ^2>" > "tmp.ef.tab";
#    }
    iq = 1;
    while ( iq < iqmax + 1 ) {
#        if ( calc == "pes" ) {
          printf("%3.0f %8.3f %8.3f  %8.3f %8.3f  %8.3f %8.3f   %8.3f  %8.3f\n",
               iq,Fermi[iq,1],Fermi[iq,2],Lambda2[iq,1],Lambda2[iq,2],ELN[iq,1], ELN[iq,2], Dispersion[iq,1], Dispersion[iq,2]) >> "tmp.ef.tab"; 
#        }
        iq += 1;
    }
    close("tmp.ef.tab")
    
    #----------------------------------------------
    # Deciding on the header for the SPWF files
    #
    header = "*             index     <P>       <Rz>       v^2      Esp       <Jx>      <Jy>      <Jz>       J        <r^2> \n"
       
    #---------------------------------------------------------------------------
    # Neutron HF basis
    printf(header) > "tmp.n.hf.tab"
    
    iq = 1;
    while ( iq < iqmax +1  ) {
        N = 1;
        while ( neutronhf[iq,N,1] != "" ){
#            if ( calc == "pes" ) {
                i = 1
                printf("%4i %4i", N, iq) >> "tmp.n.hf.tab"
                while( neutronhf[iq,N,i] != "" ) {
                    printf("%10.3f", neutronhf[iq,N,i] ) >> "tmp.n.hf.tab"
                    i +=1
                }
                
                printf("\n" ) >> "tmp.n.hf.tab"
#            }
              
            N+=1
        }
        iq+=1
    }
    close("tmp.n.hf.tab")

    #---------------------------------------------------------------------------
    # Neutron canonical basis (if HFB is active)
    if(PairingType == "HFB") {
        printf(header) > "tmp.n.can.tab"
        iq = 1;
        while ( iq < iqmax +1  ) {
            N = 1;
            while ( neutroncan[1,N,1] != "" ){
#                if ( calc == "pes" ) {
                    i = 1
                    printf("%4i %4i", N, iq) >> "tmp.n.can.tab"
                    while( neutroncan[iq,N,i] != "" ) {
                        printf("%10.3f", neutroncan[iq,N,i] ) >> "tmp.n.can.tab"
                        i +=1
#                    }
                printf("\n" ) >> "tmp.n.can.tab"
                }
                N+=1
            }
          
            iq+=1
        }
        close("tmp.n.can.tab")
    }
    #---------------------------------------------------------------------------
    # Proton HF basis
    iq = 1;
    printf(header) > "tmp.p.hf.tab"
    while ( iq < iqmax +1 ) {
        P = 1;
        while ( protonhf[iq,P,1] != ""){
#            if ( calc == "pes" ) {
                i = 1
                printf("%4i %4i", P, iq) >> "tmp.p.hf.tab"
                while( protonhf[iq,P,i] != "" ) {
                    printf("%10.3f", protonhf[iq,P,i] ) >> "tmp.p.hf.tab"
                    i +=1
                }
                
                printf("\n" ) >> "tmp.p.hf.tab"
#            }
            P+=1
        }
        iq+=1
    }
    close("tmp.p.hf.tab")
    
    #---------------------------------------------------------------------------
    # Proton canonical basis
    if(PairingType == "HFB"){
        printf(header) > "tmp.p.can.tab"
        iq = 1;
        while ( iq < iqmax +1 ) {
            P = 1;
            while ( protoncan[iq,P,1] != "" ){
#                if ( calc == "pes" ) {
                    i = 1
                    printf("%4i %4i", P, iq) >> "tmp.p.can.tab"
                    while( protoncan[iq,P,i] != "" ) {
                        printf("%10.3f", protoncan[iq,P,i] ) >> "tmp.p.can.tab"
                        i +=1
                    }
                    
                    printf("\n" ) >> "tmp.p.can.tab"
#                }
                  
                P+=1
            }
            iq+=1
        }
    close("tmp.p.can.tab")    
    }
#    exit
#    
#    #---------------------------------------------------------------------------
#    #Sort all of the SPWFs into blocks by their quantum number
    SortSpwfs("tmp.n.hf.tab", "neutron", "hf", prefix, PC, SC, TRC,iqmax) ;
    SortSpwfs("tmp.p.hf.tab", "proton" , "hf", prefix, PC, SC, TRC,iqmax) ;
    if(PairingType == "HFB") {
        SortSpwfs("tmp.n.can.tab", "neutron", "ca", prefix, PC, SC, TRC,iqmax) ;
        SortSpwfs("tmp.p.can.tab", "proton" , "ca", prefix, PC, SC, TRC,iqmax) ;
    }
    
    #---------------------------------------------------------------------------
    # Move all the temporary files to correctly named ones. 
    # Do this via a temporary script in order to get the correct bash shell. 
    # In a more naive way on my ubuntu installation it uses the default /bin/sh 
    # which is dash and not bash, and does not accept substring substitutions.
    
    printf("#!/bin/bash \n")               >  "script.sh"
    printf("for f in tmp*tab \n")          >> "script.sh"
    printf("do \n")                        >> "script.sh"
    printf("mv $f ${f/tmp/" prefix"}  \n") >> "script.sh"
    printf("done \n")                      >> "script.sh"
    
    close("script.sh")
    system("bash script.sh")
    system("rm script.sh")
    
    
}


