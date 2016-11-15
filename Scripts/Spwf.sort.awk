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
commentplus=1
commentmin =1
}

{
        #-----------------------------------------------------------------------
        # Start to read the file with the datasets.
        # ----------------------------------------------------------------------
                        
        if( $column == 1.0 ) {
            printf($0)   >> "tmp.p"
            printf("\n") >> "tmp.p"
            commentplus=0
        }
        else if ( $column == -1.0) {
            printf($0)   >> "tmp.m"
            printf("\n") >> "tmp.m"
            commentmin=0
        }
        else if ( NF == 1 && commentmin==0) {
            printf($0)   >> "tmp.m"
            printf("\n") >> "tmp.m"
            commentmin=1
        }    
        else if ( NF == 1 && commentplus==0) {    
            printf($0)   >> "tmp.p"
            printf("\n") >> "tmp.p"
            commentplus = 1
        }
                
        close("tmp.p")
        close("tmp.m")
}


