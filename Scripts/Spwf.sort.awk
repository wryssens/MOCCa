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
#
# These flags track whether or not the last line written to file was a 
# commentline or not, in order to avoid having multiple comment lines one after
# another. 
#
countplus=0
countmin=0
}
{
        #-----------------------------------------------------------------------
        # Start to read the file with the datasets.
        # ----------------------------------------------------------------------
                        
        if( $column == 1.0 ) {
            printf($0)   >> "tmp.p"
            printf("\n") >> "tmp.p"
            countplus += 1
            if(countplus == points){
                countplus = 0
                printf("*\n") >> "tmp.p"
            }
        }
        else if ( $column == -1.0) {
            printf($0)   >> "tmp.m"
            printf("\n") >> "tmp.m"
            countmin += 1
            if(countmin == points){
                countmin = 0
                printf("*\n") >> "tmp.m"
            }
        }
                       
        close("tmp.p")
        close("tmp.m")
}


