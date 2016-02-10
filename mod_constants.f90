module constants
!
! Contains constants related to the grid.
  implicit none

  ! Time-stepping constants.
  real            :: dt
  real            :: t_final
  integer         :: nt

  ! Grid constants.
  integer         :: n                 ! Number of collocation points in each direction per subdomain
  integer         :: mx                ! Number of elements in x direction
  integer         :: my                ! Number of elements in y direction
  integer         :: mz                ! Number of elements in z direction
  integer         :: r                 ! Number of total grid points.

  ! Constants related to the local grid.
  integer         :: r_loc             ! Number of grid points per rank.
  integer         :: mx_loc            ! Number of x-elements per rank.
  integer         :: my_loc            ! Number of y-elements per rank.
  integer         :: mz_loc            ! Number of z-elements per rank.
  integer         :: nx_loc            ! Number of grid points on this rank in x.
  integer         :: ny_loc            ! Number of grid points on this rank in y.
  integer         :: nz_loc            ! Number of grid points on this rank in z.

  ! Constants related to physics.
  real            :: mu                ! Dynamic viscosity.

  ! Constants related to the grid/simulation.
  real, parameter :: big_value  = 1e25 ! If any variable exceeds this norm the code reports a warning for instability.
  real, parameter :: metric_tol = 1e-6 ! Lowest allowable value for the norm of the determinant of the metric terms.

  ! Some parallelism constants.
  !
  ! NOTE:  MPI's interface uses 4-byte integers.  explicitly define them as
  !        such rather than relying upon compiler defaults/command lines.
  integer(kind=4) :: nprocs, rank, root

  ! Adams-Bashforth coefficients.
  real            :: a1, a2, a3

end module constants
