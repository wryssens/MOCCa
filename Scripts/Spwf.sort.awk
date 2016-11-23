#==============================================================================#
#                                                                              #
#    Small auxiliary script for MOCCa.spectra.awk                              #
#    14 November 2016                                                          #
#                                                                              #
#                                                                              #
#    Copyright W. Ryssens & M. Bender                                          #
#    Heavily inspired on cr8.spectra.1.1.2.awk by M. Bender and engineered     #                                                   #
#    to produce the same kind of tables, only starting from MOCCa files.       #
#                                                                              #
#                                                                              #
#===============================================================================
BEGIN{

}
{
        #-----------------------------------------------------------------------
        # Start to read the file with the datasets.
        # Notice that this ignores any lines that are not spwf
        # ----------------------------------------------------------------------
        if( column == -1 && NF !=1) {
            if( $2 > qmax) {
                qmax = $2
            }
            if(Nz[$2] == "") {
                Nz[$2] = 0 
            }
            Nz[$2] += 1 
            zero[$2, Nz[$2]] = $0
        }        
        else if( NF != 1 ) {
            if( $2 > qmax) {
            qmax = $2
            }
            if( $column == 1.0) {
                
                if(NP[$2] == "") {
                    NP[$2] = 0 
                }
                NP[$2] += 1
                plus[$2, NP[$2]] = $0
            }
            else if( $column == -1.0) {
                
                if(NM[$2] == "") {
                    NM[$2] = 0 
                }
                NM[$2] += 1
                min[$2, NM[$2]] = $0
            }                      
        }                         
}

END{
    if(column != -1) {
        i = 1
        while( NP[1] + 1 > i ){
            iq = 1
            while ( iq < qmax +1){
                   printf(plus[iq,i]) >> "tmp.p"
                   printf("\n") >> "tmp.p"
                   iq +=1
            }
            printf("*\n") >> "tmp.p"
            i +=1 
        }
        
        i = 1
        while( NM[1] + 1 > i ){
            iq = 1
            while (iq < qmax +1 ){
                   printf(min[iq,i]) >> "tmp.m"
                   printf("\n")      >> "tmp.m"
                   iq +=1
            }
            printf("*\n") >> "tmp.m"
            i +=1 
        }

        close("tmp.p")
        close("tmp.m")
    }
    else {
        i = 1
        while( Nz[1] + 1 > i ){
            iq = 1
            while ( iq < qmax +1){
                   printf(zero[iq,i]) >> "tmp.zero"
                   printf("\n") >> "tmp.zero"
                   iq +=1
            }
            printf("*\n") >> "tmp.zero"
            i +=1 
        }
    
    }
       
      
}


