module timing
!
! Contains variables related to timing.
  implicit none

  ! Timers.
  real :: timer_sync_ranks  ! Time spent doing rank synchronization.
  real :: timer_dgemm       ! Time spent computing derivatives.
  real :: timer_shuffle     ! Time spent transposing data.

end module timing
