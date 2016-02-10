subroutine read_input_file( fname )
!
! Reads an _in file and parses it.

  use constants
  use geom
  use options

  implicit none

  character(len=*), intent(in)  :: fname
  character(len=132)            :: line_buffer
  character(len=132)            :: line_buffer2
  integer                       :: ii, readstat, pound_ndx, ndx

  ! Open the file for reading.
  open( 200, file = trim( fname ) )

  ! Read the file line-by-line.
  do ii = 1, 1000

     ! Read next line.
     read( 200, '(A)', iostat=readstat ) line_buffer

     ! Left adjust and trim the buffer.
     line_buffer = trim( adjustl( line_buffer ) )

     ! If we're at the end of the file, return.
     if ( readstat < 0 ) then
        return
     endif

     ! Ignore comment lines and blank lines.
     if ( index( line_buffer, "#" ) /= 1 .and. len( trim( adjustl( line_buffer ) ) ) > 0 ) then

        ! Strip out any whitespace.
        call strip_all_whitespace( line_buffer )

        ! Strip out any trailing comments.
        pound_ndx = index( line_buffer, "#" )
        if (pound_ndx /= 0 ) line_buffer = line_buffer(1:pound_ndx-1)

        ! To upper case.
        line_buffer2 = line_buffer
        call to_upper( line_buffer )

        ! Parse the string for data.

           if ( index( line_buffer, 'REPORT_EVERY_N_STEPS=') ==  1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) report_every_n_steps
              cycle
           endif

           if ( index( line_buffer, 'CHECK_CONSISTENCY_FLAG=') ==  1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) check_consistency_flag
              cycle
           endif

           if ( index( line_buffer, 'N=') ==  1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) n
              cycle
           endif

           if ( index( line_buffer, 'MX=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) mx
              cycle
           endif

           if ( index( line_buffer, 'MY=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) my
              cycle
           endif

           if ( index( line_buffer, 'MZ=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) mz
              cycle
           endif

           if ( index( line_buffer, 'TIMESTEPS_BETWEEN_WRITES=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) timesteps_between_writes
              cycle
           endif

           if ( index( line_buffer, 'DT=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) dt
              cycle
           endif

           if ( index( line_buffer, 'CFL=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) CFL
              cycle
           endif

           if ( index( line_buffer, 'T_FINAL=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer(ndx+1:), * ) t_final
              cycle
           endif

           if ( index( line_buffer, 'FNAME_RUNNAME=' ) == 1 ) then
              ndx           = index( line_buffer, '=' )
              fname_runname = trim( adjustl( line_buffer2(ndx+1:) ) )
              cycle
           endif

           if ( index( line_buffer, 'FNAME_INIT=' ) == 1 ) then
              ndx           = index( line_buffer, '=' )
              fname_init    = trim( adjustl( line_buffer2(ndx+1:) ) )
              cycle
           endif

           if ( index( line_buffer, 'LX=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) Lx
              cycle
           endif

           if ( index( line_buffer, 'LY=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) Ly
              cycle
           endif

           if ( index( line_buffer, 'LZ=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) Lz
              cycle
           endif

           if ( index( line_buffer, 'MU=' ) == 1 ) then
              ndx = index( line_buffer, '=' )
              read( line_buffer2(ndx+1:), *) mu
              cycle
           endif

        ! If we get here, throw an error.
        write(*,*) 'Error: Unclassifiable statement at line ', ii
        write(*,*) line_buffer
        stop

     else

        ! Do nothing.  This line is either blank or #-commented.
        ! XXX: This is broken.  If a user adds comments at the end of a line of meaningful input, this breaks.
        ! I should figure out a way to test for this.

     endif

  enddo

end subroutine read_input_file

subroutine strip_all_whitespace( str )

  implicit none

  character(len=*), intent(inout)  :: str
  character(len=len( str ))        :: str2
  integer                          :: ii, counter

  counter = 1
  do ii = 1, len( str )
     if ( str(ii:ii) .ne. ' ' ) then
        str2(counter:counter) = str(ii:ii)
        counter               = counter + 1
     endif
  enddo

  str2(counter:) = ' '
  str            = str2

end subroutine strip_all_whitespace

subroutine to_upper( strIn )
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)

     implicit none

     character(len=*), intent(inout) :: strIn
     character(len=len(strIn))       :: strOut
     integer                         :: i,j

     do i = 1, len( strIn )
          j = iachar( strIn(i:i) )
          if (j >= iachar( "a" ) .and. j <= iachar( "z" ) ) then
               strOut(i:i) = achar( iachar( strIn(i:i) )-32 )
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

     strIn = strOut

end subroutine to_upper
