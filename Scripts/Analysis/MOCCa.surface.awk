#==============================================================================#
#                                                                              #
#    MOCCa.surface.awk     release 0.0.1                                       #
#    29 November 2016                                                          #
#                                                                              #
#==============================================================================#
#                                                                              #
# Awk script that gets the energy out of a MOCCa file together with a          #
# value of a certain multipole moment (or related quantity) as a entry into    #
# a table for plotting an energy surface.                                      #
#                                                                              #
# Usage                                                                        #
#       awk -f MOCCa.surface.awk "l1=XXXX" "l2=YYYY" "m1=ZZZZ" "m2=QQQQ"       #
#              "leg=1"                                                         #
#                                                                              #
# Where (l1,m1) indicate the first multipole moment to use, while (l2,m2)      #
# indicates the second.                                                        #
#                                                                              #
#==============================================================================#
BEGIN{
  iflag = 0;
}

{
  if ( $1 == "**FINAL**" ) {
    iflag   = 1;
  }
  if ( iflag == 1 ) {
  
    if (leg == 0) {
        if ( $1 == "Beta_{" && $2 == l1 && $3 == m1 "}"  ) {
            beta1 = $6
        }
        if ( $1 == "Beta_{" && $2 == l2 && $3 == m2 "}"  ) {
            beta2 = $6
        }
        
    }
    if ( leg == 1) {
        if ($1 == "iq" && $2 == "T") {
            beta1 = $3
            beta2 = $4
        }
    }
    if ( $1 == "Lagrange:" ) {
            ee = $2;
    }
  }
}
END{
   if( ee != 0) {
    # Don't include the file if it is empty and does not contain a final energy
    printf("  %8.3f %8.3f %12.3f \n",beta1,beta2,ee);
   }
}
