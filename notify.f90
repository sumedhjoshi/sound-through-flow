subroutine notify( nstring )
! This subroutine will pre-pend a datestring to the message nstring:
! 29-Mar-2013 12:41:33: nstring

  implicit none

  character(len=*), intent(in) :: nstring
  character(len=100)           :: date, time
  integer, dimension(8)        :: T
  character(len=10)            :: zone
  character(len=3)             :: month

  call date_and_time( date, time, zone, T )

  select case (T(2))
     case (1)
        month = 'Jan'
     case (2)
        month = 'Feb'
     case (3)
        month = 'Mar'
     case (4)
        month = 'Apr'
     case (5)
        month = 'May'
     case (6)
        month = 'Jun'
     case (7)
        month = 'Jul'
     case (8)
        month = 'Aug'
     case (9)
        month = 'Sep'
     case (10)
        month = 'Oct'
     case (11)
        month = 'Nov'
     case (12)
        month = 'Dec'
   end select

   write(*,*) date(7:8), '-', month, '-', date(1:4), ' ', time(1:2), ':', time(3:4), ':', time(5:6), ':  ', nstring

end subroutine notify
