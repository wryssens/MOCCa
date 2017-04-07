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
   if ( $1 == "Alirotation") {
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
  }
}
END{
   if( ee != 0) {
    # Don't include the file if it is empty and does not contain a final energy
    printf("  %8.3f %8.3f %12.3f %8.3f %8.3f\n",alix,aliy,ee,jx,jz);
   }
}
