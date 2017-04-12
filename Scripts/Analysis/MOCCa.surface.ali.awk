#==============================================================================#
#                                                                              #
#    MOCCa.ali.surface.awk     release 0.0.1                                   #
#    29 November 2016                                                          #
#                                                                              #
#==============================================================================#
#                                                                              #
# Awk script that gets the energy out of a MOCCa file together with the values #
# of the ali-rotation angles for up to two quasiparticles.                     # 
#                                                                              #
# Usage                                                                        #
#       awk -f MOCCa.surface.ali.awk                                           #
#                                                                              #
# Where (l1,m1) indicate the first multipole moment to use, while (l2,m2)      #
# indicates the second.                                                        #
#                                                                              #
#==============================================================================#
BEGIN{
  iflag = 0;
}

{
   
  if ( $2 =="The"  && $3 == "HFBasis" && $8=="alirotated"){
    basis = 1 ;
    getline;
    getline;
    ali1 = $3
    ali3 = $4
    getline;
    ali2 = $3
    ali4 = $4
  }
  
  if ( $1 == "Alirotation" && basis == 0) {
    getline ;
    getline ;   
    alix = $2
    getline ;
    aliy = $2
  }


  if ( $1 == "**FINAL**" ) {
    iflag   = 1;
  }
  if ( iflag == 1 ) {
    if ( $1 == "Lagrange:" ) {
        ee = $2;
    }
    if ( $1 == "J_z" ) {
        jz = $3;
    }
    if ( $1 == "J_x" ) {
        jx = $3;
    }
    if ( $1 == "Re" && $2 == "Q^m_{" && $3 == "1" && $4 == "0}"){
        for (i=0 ; i < 12 ; i+=1){
            getline ;
        }
        Qm10 = $5    
    }
    if ( $1 == "Re" && $2 == "Q^m_{" && $3 == "1" && $4 == "1}"){
        for (i=0 ; i < 12 ; i+=1){
            getline ;
        }
        Qm11 = $5    
    }
    
  }
}
END{
   if( ee != 0 && basis==0) {
    # Don't include the file if it is empty and does not contain a final energy
    printf("  %8.3f %8.3f %12.3f %8.3f %8.3f %8.4f %8.4f\n",alix,aliy,ee,jx,jz, Qm10, Qm11);
   }
   if( ee != 0 && basis==1) {
    # Don't include the file if it is empty and does not contain a final energy
    printf("  %8.3f %8.3f %8.3f %8.3f %12.3f %8.3f %8.3f %8.4f %8.4f\n",ali1,ali2,ali3,ali4,ee,jx,jz, Qm10, Qm11);
   }
}
