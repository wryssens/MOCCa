!-------------------------------------------------------------------------------
! This module contains all the routines that allow MOCCa to output graphs 
! using Gnuplot. It is entirely optional, MOCCa can run on systems that do not 
! have Gnuplot available, 
! but users will miss out on automatically generated fancy graphs.
!
! This module was entirely written by prof. J. Burkardt and taken from
!    http://people.sc.fsu.edu/~jburkardt/f_src/gnufor/gnufor.html
! It is used here under the terms of the GNU LGPL license, which can be found at
!    http://www.gnu.org/licenses/lgpl.html
!-------------------------------------------------------------------------------
! Some small additions have been made to facilitate the inclusion into the MOCCa
! program. These changes include:
!   - subroutine write_xyzgrid_contour now takes extra arguments (xlabel,ylabel)
!     which determine the labeling of the axes.
!   - changed all "set term x11" to "set term postscript". Gnuplot no longer 
!     opens a window but prints output to file.
!   - added "Outputname" argument to all command-writing subroutines. 
!     This variable determines the name of the output postscript file.
!   - changed "end" statements to "end subroutine <subroutinename>" 
!     for use with the ifort compiler.
!   - subroutine GetUnit was moved to module GenInfo
!   - In general: suppressed output to standard out when this module runs 
!     correctly. Kept the output when errors are encountered.
!-------------------------------------------------------------------------------

module GnuFor

    use CompilationInfo
    use GenInfo

    implicit none

    public
    
    contains
    
    subroutine run_gnuplot ( command_filename )
        !*****************************************************************************
        !
        !! RUN_GNUPLOT runs GNUPLOT with a given command file.
        !
        !  Discussion:
        !
        !    The GNUPLOT program must be available.  To check whether
        !    this is so, try typing
        !
        !      which gnuplot
        !
        !    If the response is
        !
        !      gnuplot: command not found
        !
        !    then you're going to have to make GNUPLOT available.
        !
        !    You may need to set the environment variable GNUTERM:
        !
        !      setenv GNUTERM x11
        !
        !    so that GNUPLOT automatically displays to your X window terminal.
        !
        !
        !    This routine expects that there is a text file containing the appropriate
        !    commands to GNUPLOT to display your picture.  There are a number of
        !    routines in this package that will do this for simple plotting tasks.
        !    Most of them require that you also set up a file of data to be plotted.
        !
        !    Once this routine invokes GNUPLOT, a graphics window should open
        !    up, and the FORTRAN program will pause.  Hitting RETURN should advance
        !    to the next picture, or terminate the window at the end, allowing the
        !    FORTRAN routine to proceed.
        !
        !
        !    You can look at the data and command files created by the routines.
        !    Moreover, you can easily modify the command file to change the options
        !    used in GNUPLOT, and then run GNUPLOT interactively, as in:
        !
        !      gnuplot commands
        !
        !    In particular, if you want a PostScript version of your graphics files,
        !    insert the command "set term postscript" at the beginning of the command
        !    file and run gnuplot as follows:
        !
        !      gnuplot commands > mypicture.ps
        !
        !    You will also have to hit RETURN once for each plot that is made.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    21 February 2001
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) COMMAND_FILENAME, the name of the
        !    command file.
        !------------------------------------------------------------------------------
        ! Change by WR: suppressed some output to standard out.
        !------------------------------------------------------------------------------
         

        character ( len = 255 ) command
        character ( len = * ) command_filename
        integer ( kind = 4 ) status
        integer ( kind = 4 ) system

        !call timestamp ( )
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'GNUFOR:'
!        write ( *, '(a)' ) '  GNUPLOT / FORTRAN90 command interface.'
        !
        !  Issue a command to the system that will start GNUPLOT, using
        !  the file we just wrote as input.
        !
        !  The "&" will run GNUPLOT in the background, so the FORTRAN program
        !  can continue execution independently, and the PERSIST switch tells
        !  GNUPLOT that if there are multiple plots, they should each go in
        !  a separate window.
        !
        !  Thanks to Morag Am-Shallem for suggesting these improvements.
        !  17 October 2007
        !
        write ( command, * ) 'gnuplot -persist ' // trim ( command_filename ) // ' &'

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'GNUFOR:'
!        write ( *, '(a)' ) '  Issuing the command:"' // trim ( command ) // '".'
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) '  Press RETURN to proceed.'

        status = system ( trim ( command ) )

        if ( status /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GNUFOR - Fatal error!'
        write ( *, '(a)' ) '  An error code was returned when the GNUPLOT command'
        write ( *, '(a)' ) '  was issued.  Perhaps GNUPLOT is not in your path.'
        write ( *, '(a)' ) '  Type "which gnuplot" to check this.'
        stop
        end if
        !
        !  Terminate.
        !
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'GNUFOR:'
!        write ( *, '(a)' ) '  Normal end of execution.'

!        write ( *, '(a)' ) ' '
!        call timestamp ( )

        return
    end subroutine run_gnuplot
    subroutine timestamp ( )

        !*****************************************************************************80
        !
        !! TIMESTAMP prints the current YMDHMS date as a time stamp.
        !
        !  Example:
        !
        !    31 May 2001   9:45:54.872 AM
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    06 August 2005
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    None
        !
         

        character ( len = 8 ) ampm
        integer ( kind = 4 ) d
        integer ( kind = 4 ) h
        integer ( kind = 4 ) m
        integer ( kind = 4 ) mm
        character ( len = 9  ), parameter, dimension(12) :: month = (/ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)
        integer ( kind = 4 ) n
        integer ( kind = 4 ) s
        integer ( kind = 4 ) values(8)
        integer ( kind = 4 ) y

        call date_and_time ( values = values )

        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

        if ( h < 12 ) then
        ampm = 'AM'
        else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
        ampm = 'Noon'
        else
        ampm = 'PM'
        end if
        else
        h = h - 12
        if ( h < 12 ) then
        ampm = 'PM'
        else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
        else
        ampm = 'AM'
        end if
        end if
        end if

        write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

        return
    end subroutine timestamp
    subroutine write_vector_data ( data_filename, n, x, y, dx, dy, ierror )

        !*****************************************************************************80
        !
        !! WRITE_VECTOR_DATA writes vector data to a file, for plotting by GNUPLOT.
        !
        !  Discussion:
        !
        !    Each vector is described by 4 values, X, Y, dX, dY, indicating that
        !    a vector is to be drawn from (X,Y) to (X+dX,Y+dY).
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    22 February 2001
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Input, integer ( kind = 4 ) N, the number of vectors.
        !
        !    Input, real ( kind = dp ) X(N), Y(N), DX(N), DY(N), the vector data.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
         

        integer ( kind = 4 ) n

        character ( len = * ) data_filename
        real ( kind = dp ) dx(n)
        real ( kind = dp ) dy(n)
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        real ( kind = dp ) x(n)
        real ( kind = dp ) y(n)

        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_VECTOR_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = data_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_VECTOR_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        do i = 1, n
        write ( file_unit, * ) x(i), y(i), dx(i), dy(i)
        end do

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_VECTOR_DATA:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT vector data file "' // &
!        trim ( data_filename ) // '"'

        return
    end subroutine write_vector_data
    subroutine write_vector_plot ( command_filename, data_filename, OutputFilename,    ierror )
        !*****************************************************************************80
        !
        !! WRITE_VECTOR_PLOT writes GNUPLOT commands to plot vectors.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    24 June 2011
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) COMMAND_FILENAME, the name of the
        !    command file.
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !----------------------------------------------------------------------------
        ! Changes by WR:
        !    Output changed from opening a window to writing to file
        !    Added argument OutputFilename, determines the name of outputfile
        !----------------------------------------------------------------------------
         

        character ( len = * ) command_filename
        character ( len = * ) data_filename
        character ( len = * ) OutputFilename !Added by WR
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        !
        !  Write the data file.
        !
        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_VECTOR_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = command_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_VECTOR_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        write ( file_unit, '(a)' ) 'set term postscript' !Added by WR
        write ( file_unit, '(a)' ) 'set output "' // trim(OutputFilename) // '.ps"' !Added by WR
    
        write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
        write ( file_unit, '(a)' ) 'set xlabel "x"'
        write ( file_unit, '(a)' ) 'set ylabel "y"'
        write ( file_unit, '(a)' ) 'set style data vector'
        write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename )
        write ( file_unit, '(a)' ) 'pause -1'
        write ( file_unit, '(a)' ) 'q'

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_VECTOR_PLOT:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT table plots command file "' // &
!        trim ( command_filename ) // '"'

        return
    end subroutine write_vector_plot
    subroutine write_xy_data ( data_filename, n, x, y, ierror )

        !*****************************************************************************80
        !
        !! WRITE_XY_DATA writes X(1:N), Y(1:N) data to a file.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    23 February 2001
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Input, integer ( kind = 4 ) N, the number of data items.
        !
        !    Input, real ( kind = dp ) X(N), Y(N), the X and Y data
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
         

        integer ( kind = 4 ) n

        character ( len = * ) data_filename
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        real ( kind = dp ) x(n)
        real ( kind = dp ) y(n)

        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = data_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        do i = 1, n
        write ( file_unit, * ) x(i), y(i)
        end do

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_XY_DATA:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT XY data file "' // &
!        trim ( data_filename ) // '"'

        return
    end subroutine write_xy_data
    subroutine write_xy_plot ( command_filename, data_filename, OutputFilename,XLabel,Ylabel, ierror )
        !*****************************************************************************80
        !
        !! WRITE_XY_PLOT writes GNUPLOT commands to plot X(1:N), Y(1:N) data.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    23 June 2011
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) COMMAND_FILENAME, the name of the
        !    command file.
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
        !----------------------------------------------------------------------------
        ! Changes by WR:
        !    Output changed from opening a window to writing to file
        !    Added argument OutputFilename, determines the name of outputfile
        !----------------------------------------------------------------------------
         

        character ( len = * ) command_filename
        character ( len = * ) data_filename
        character ( len = * ) OutputFilename    !Added by WR
        character ( len = * ) Xlabel    !Added by WR
        character ( len = * ) Ylabel   !Added by WR
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        !
        !  Write the data file.
        !
        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = command_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if
    
        write ( file_unit, '(a)' ) 'set term png' !Added by WR
        write ( file_unit, '(a)' ) 'set output "' // trim(OutputFilename) // '.png"' !Added by WR
        write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
        write ( file_unit, '(a)' ) 'set xlabel "' // trim(xlabel) // '"'
        write ( file_unit, '(a)' ) 'set ylabel "' // trim(ylabel) // '"'
        write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename ) // &
        '" using 1:2 with points'
        write ( file_unit, '(a)' ) 'pause -1'
        write ( file_unit, '(a)' ) 'q'

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_XY_PLOT:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT XY plot command file "' // &
!        trim ( command_filename ) // '"'

        return
    end subroutine write_xy_plot
    subroutine write_xyy_data ( data_filename, lda, nrow, ncol, x, ierror )

        !*****************************************************************************80
        !
        !! WRITE_XYY_DATA writes a table of data to a file, for plotting by GNUPLOT.
        !
        !  Discussion:
        !
        !    The first column of data is assumed to be the independent variable, X.
        !    Separate plots are made of X versus all the other columns of data.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    21 February 2001
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Input, integer ( kind = 4 ) LDA, the leading dimension of X.
        !
        !    Input, integer ( kind = 4 ) NROW, NCOL, the dimensions of X.
        !
        !    Input, real ( kind = dp ) X(LDA,NCOL), the NROW by NCOL data to be written.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
         

        integer ( kind = 4 ) lda
        integer ( kind = 4 ) ncol

        character ( len = * ) data_filename
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        integer ( kind = 4 ) nrow
        real ( kind = dp ) x(lda,ncol)

        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = data_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        do i = 1, nrow
        write ( file_unit, * ) x(i,1:ncol)
        end do

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_XYY_DATA:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT XYY data file "' // &
!        trim ( data_filename ) // '"'

        return
    end subroutine write_xyy_data
    subroutine write_xyy_plots ( command_filename, data_filename, &
        ncol, OutputFilename,ierror )

        !*****************************************************************************80
        !
        !! WRITE_XYY_PLOTS writes GNUPLOT commands to make multiple (X,Y) plots.
        !
        !  Discussion:
        !
        !    The first column of data is assumed to be the independent variable, X.
        !    Separate plots are made of X versus all the other columns of data.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    24 June 2011
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) COMMAND_FILENAME, the name of the
        !    command file.
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Input, integer ( kind = 4 ) NCOL, the number of columns of data.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
        !----------------------------------------------------------------------------
        ! Changes by WR:
        !    Output changed from opening a window to writing to file
        !    Added argument OutputFilename, determines the name of outputfile
        !----------------------------------------------------------------------------
         

        character ( len = * ) command_filename
        character ( len = * ) data_filename
        character ( len = * ) OutputFilename !Added by WR
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        integer ( kind = 4 ) ncol
        !
        !  Write the data file.
        !
        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = command_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        write ( file_unit, '(a)' ) 'set term postscript' ! Added by WR
        write ( file_unit, '(a)' ) 'set output "' // trim(OutputFilename) // '.ps"' !Added by WR
        write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
        write ( file_unit, '(a)' ) 'set xlabel "x"'
        write ( file_unit, '(a)' ) 'set ylabel "y"'
        do i = 2, ncol
        write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename ) // &
        '" using ', i, ' with lines'
        write ( file_unit, '(a)' ) 'pause -1'
        end do
        write ( file_unit, '(a)' ) 'q'

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_XYY_PLOTS:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT XYY plots command file "' // &
!        trim ( command_filename ) // '"'

        return
    end subroutine write_xyy_plots
    subroutine write_xyz_data ( data_filename, n, x, y, z, ierror )

        !*****************************************************************************80
        !
        !! WRITE_XYZ_DATA writes X(1:N), Y(1:N), Z(1:N) data to a file.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    23 February 2001
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Input, integer ( kind = 4 ) N, the number of data items.
        !
        !    Input, real ( kind = dp ) X(N), Y(N), Z(N), the X, Y, Z data
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
         

        integer ( kind = 4 ) n

        character ( len = * ) data_filename
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        real ( kind = dp ) x(n)
        real ( kind = dp ) y(n)
        real ( kind = dp ) z(n)

        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = data_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        do i = 1, n
        write ( file_unit, * ) x(i), y(i), z(i)
        end do

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_XYZ_DATA:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT XYZ data file "' // &
!        trim ( data_filename ) // '"'

        return
    end subroutine write_xyz_data
    subroutine write_xyz_plot ( command_filename, data_filename, OutputFilename, ierror )

        !*****************************************************************************80
        !
        !! WRITE_XYZ_PLOT writes commands to plot parametric (X,Y,Z) data.
        !
        !  Discussion:
        !
        !    This routine tries to write a command file suitable for displaying
        !    a 3D arc specified by points (X,Y,Z).  A grid data file, containing
        !    values of X, Y and Z, will also be needed.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    22 February 2001
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) COMMAND_FILENAME, the name of the
        !    command file.
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
        !----------------------------------------------------------------------------
        ! Changes by WR:
        !    Output changed from opening a window to writing to file
        !    Added argument OutputFilename, determines the name of outputfile
        !----------------------------------------------------------------------------
         

        character ( len = * ) command_filename
        character ( len = * ) data_filename
        character ( len = * ) OutputFilename !Added by WR
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        !
        !  Write the data file.
        !
        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZ_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = command_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZ_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        write ( file_unit, '(a)' ) 'set term postscript' !Added by WR
        write ( file_unit, '(a)' ) 'set output "' // trim(OutputFilename) // '.ps"' !Added by WR
        write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
        write ( file_unit, '(a)' ) 'set xlabel "x"'
        write ( file_unit, '(a)' ) 'set ylabel "y"'
        write ( file_unit, '(a)' ) 'set parametric'
        write ( file_unit, '(a)' ) 'splot "' // trim ( data_filename ) // &
        '" using 1:2:3 with lines'
        write ( file_unit, '(a)' ) 'pause -1'
        write ( file_unit, '(a)' ) 'q'

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_XYZ_PLOT:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT SPLOT command file "' // &
!        trim ( command_filename ) // '"'

        return
    end subroutine write_xyz_plot
    subroutine write_xyzgrid_contour ( command_filename, data_filename,OutputFilename, xlabel, ylabel, ierror )

        !*****************************************************************************80
        !
        !! WRITE_XYZGRID_CONTOUR writes commands to plot contours of Z(X,Y).
        !
        !  Discussion:
        !
        !    This routine tries to write a command file suitable for displaying
        !    contours of Z(X,Y) gridded data.  A grid data file, containing values
        !    of X, Y and Z, will also be needed.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    24 June 2011
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) COMMAND_FILENAME, the name of the
        !    command file.
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !----------------------------------------------------------------------------
        ! Changes by WR:
        !    Output changed from opening a window to writing to file
        !    Added argument OutputFilename, determines the name of outputfile
        !     Added xlabel and ylabel for labelling the axes
        !    Removed everything related to table.txt, I don't get it.
        !     Commented "nosurface" line. I like 3d plots :)
        !----------------------------------------------------------------------------
         

        character ( len = * ) command_filename
        character ( len = * ) data_filename
        character ( len = * ) xlabel,ylabel !Added by WR
        character ( len = * ) OutputFilename !Added by WR
        integer ( kind = 4 ) file_unit

        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        !
        !  Write the data file.
        !
        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = command_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        write ( file_unit, '(a)' ) 'set term eps' !Added by WR
        write ( file_unit, '(a)' ) 'set output "' // trim(OutputFilename) // '.eps'
        write ( file_unit, '(a)' ) 'set xlabel "' // trim( xlabel ) // '"' !Changed by WR
        write ( file_unit, '(a)' ) 'set ylabel "' // trim( ylabel ) // '"' !Changed by WR
        write ( file_unit, '(a)' ) 'set parametric'
        !write ( file_unit, '(a)' ) 'set nosurface'
        write ( file_unit, '(a)' ) 'set contour'
        write ( file_unit, '(a)' ) 'set cntrparam levels auto 15'
        write ( file_unit, '(a)' ) 'set cntrparam bspline' 
        write ( file_unit, '(a)' ) 'unset key'
        write ( file_unit, '(a)' ) 'set size ratio -1'
        write ( file_unit, '(a)' ) ' unset surface'
        write ( file_unit, '(a)' ) ' set view map' 
        write ( file_unit, '(a)' ) 'splot "' // trim ( data_filename ) // &
        '" using 1:2:3  with lines'
    
        !
        !
        !write ( file_unit, '(a)' ) ' replot '
        !write ( file_unit, '(a)' ) 'q'

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR:'
!        write ( *, '(a)' ) &
!        '  Wrote the GNUPLOT XYZGRID contour plot command file "' // &
!        trim ( command_filename ) // '"'

        return
    end subroutine write_xyzgrid_contour
    subroutine write_xyzgrid_data ( data_filename, mx, my, xyz, ierror )

        !*****************************************************************************80
        !
        !! WRITE_XYZGRID_DATA writes a file of XYZ grid data.
        !
        !  Discussion:
        !
        !    It is assumed that values of Z are available on a regular mx by my grid
        !    of (X,Y) points.
        !
        !    The form of the data file requires that all the data for a given value
        !    of Y be listed, followed by a blank line, followed by the data for
        !    another value of Y.
        !
        !  Example:
        !
        !    Here is a grid data file for a 3 by 3 grid, with Z = X + Y.
        !
        !    0.0 0.0 0.0
        !    1.0 0.0 1.0
        !    2.0 0.0 2.0
        !
        !    0.0 1.0 1.0
        !    1.0 1.0 2.0
        !    2.0 1.0 3.0
        !
        !    0.0 2.0 2.0
        !    1.0 2.0 3.0
        !    2.0 2.0 4.0
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    23 February 2001
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Input, integer ( kind = 4 ) mx, my, the dimensions of the grid.
        !
        !    Input, real ( kind = dp ) XYZ(3,mx,my), the XYZ grid data to be written.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
         

        integer ( kind = 4 ) mx
        integer ( kind = 4 ) my

        character ( len = * ) data_filename
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        integer ( kind = 4 ) j
        real ( kind = dp ) xyz(3,mx,my)

        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = data_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        do j = 1, my
        do i = 1, mx
        write ( file_unit, * ) xyz(1:3,i,j)
        end do
        write ( file_unit, '(a)' )
        end do

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_XYZGRID_DATA:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT XYZ grid data file "' // &
!        trim ( data_filename ) // '"'

        return
    end subroutine write_xyzgrid_data
    subroutine write_xyzgrid_surface ( command_filename, data_filename, OutputFilename, ierror )

        !*****************************************************************************80
        !
        !! WRITE_XYZGRID_SURFACE writes a file of GNUPLOT commands to plot a 3D surface.
        !
        !  Discussion:
        !
        !    This routine tries to write a command file suitable for displaying
        !    a surface Z(X,Y).  A grid data file, containing values of X, Y and Z,
        !    will also be needed.  The routine WRITE_XYZGRID_DATA can write this file.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    22 February 2001
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, character ( len = * ) COMMAND_FILENAME, the name of the
        !    command file.
        !
        !    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
        !
        !    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
        !
        !----------------------------------------------------------------------------
        ! Changes by WR:
        !    Output changed from opening a window to writing to file
        !    Added argument OutputFilename, determines the name of outputfile
        !----------------------------------------------------------------------------
        character ( len = * ) command_filename
        character ( len = * ) data_filename
        character ( len = * ) OutputFilename ! Added by WR
        integer ( kind = 4 ) file_unit
        integer ( kind = 4 ) ierror
        integer ( kind = 4 ) ios
        !
        !  Write the data file.
        !
        ierror = 0

        call get_unit ( file_unit )

        if ( file_unit == 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_SURFACE - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
        end if

        open ( unit = file_unit, file = command_filename, status = 'replace', &
        iostat = ios )

        if ( ios /= 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_SURFACE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
        end if

        write ( file_unit, '(a)' ) 'set term postscript' !Added by WR
        write ( file_unit, '(a)' ) 'set output "' // trim(OutputFilename) // '.ps"' !Added by WR
        write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
        write ( file_unit, '(a)' ) 'set xlabel "x"'
        write ( file_unit, '(a)' ) 'set ylabel "y"'
        write ( file_unit, '(a)' ) 'set parametric'
        write ( file_unit, '(a)' ) 'set hidden3d'
        write ( file_unit, '(a)' ) 'set contour'
        write ( file_unit, '(a)' ) 'splot "' // trim ( data_filename ) // &
        '" using 1:2:3 with lines'
        write ( file_unit, '(a)' ) 'pause -1'
        write ( file_unit, '(a)' ) 'q'

        close ( unit = file_unit )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'WRITE_SURFACE_COMMANDS:'
!        write ( *, '(a)' ) '  Wrote the GNUPLOT surface plot command file "' // &
!        trim ( command_filename ) // '"'

        return
    end subroutine write_xyzgrid_surface

end module GnuFor
