module options
!
! Contains execution options.
! The values set here are defaults.

  implicit none

  ! Files.
  character(len=100)      :: fname_init
  character(len=100)      :: fname_runname

  ! Execution options.
  integer :: report_every_n_steps = 10 ! Number of time-steps between every report.
  logical :: check_consistency_flag = .true.
  integer :: timesteps_between_writes = 1
  real    :: CFL = 1.0

end module options
