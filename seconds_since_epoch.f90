subroutine seconds_since_epoch( timestamp )
! Computes the fractional seconds since midnight on January 1st, 1970 in UTC.

  implicit none

  ! Fractional seconds since midnight on January 1st, 1970 in UTC.
  real, intent(out)       :: timestamp

  ! Used to compute timestamp.
  integer, dimension(1:8) :: current_time_array
  integer                 :: todays_day
  integer                 :: epoch_day

  ! Get the current time, in component form.
  call date_and_time( VALUES=current_time_array )

  ! Compute the number of days from some time in the distant past to the Unix
  ! Epoch and today.
  call julian_day( 1970, 1, 1, epoch_day )
  call julian_day( current_time_array(1), current_time_array(2), &
                   current_time_array(3), todays_day )

  ! Compute the current time, in fractional epoch seconds, in UTC.  Note
  ! that the 4th component of the current time is an offset to be added
  ! to convert from UTC to the local time, so we have to subtract it to
  ! go back to UTC.
  timestamp = (todays_day - epoch_day) * 86400. + &
              current_time_array(5) * 3600. + &
              (current_time_array(6) - current_time_array(4)) * 60. + &
              current_time_array(7) + &
              current_time_array(8) / 1000.
return

end subroutine seconds_since_epoch

subroutine julian_day( year, month, day, jday )
! Computes the Julian day for the supplyed year, month, and day of month.
! The computed day is relative to some date in the distant past, such that
! January 1st, 1970 is Julian day 2,440,588.
!
! Adapted from the equation given in Communications of the ACM, Volume 11,
! Issue 10 of October 1968, page 657.  See "Letters to the editor: a machine
! algorithm for processing calendar dates" by Henry F. Fliegel and Thomas
! C. Van Flandern.

  implicit none

  integer, intent(in)  :: year
  integer, intent(in)  :: month
  integer, intent(in)  :: day

  integer, intent(out) :: jday

  ! Magic!
  jday = day - 32075 + 1461 * (year + 4800 + (month - 14) / 12) / 4 +  &
         367 * (month - 2 - ((month - 14) / 12) * 12) / 12 - &
         3 * ((year + 4900 + (month - 14) / 12) / 100) / 4

  return

end subroutine julian_day
