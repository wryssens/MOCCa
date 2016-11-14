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
{
        #-----------------------------------------------------------------------
        # Start to read the file with the datasets.
        # ----------------------------------------------------------------------
                        
        if( $column == 1.0 ) {
            printf($0)   >> "tmp.p"
            printf("\n") >> "tmp.p"
        }
        else if ( $column == -1.0) {
            printf($0)   >> "tmp.m"
            printf("\n") >> "tmp.m"
        }
        
        close("tmp.p")
        close("tmp.m")
}


