subroutine compute_dxi( q, Dq, direction )
! Computes the derivative in the grid directions D_xi1, Dxi2, or Dxi3, which
! correspond to (logically) the x, y, and z directions.
!
  use constants
  use legendre
  use timing
  use mpi

  implicit none

  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(inout) :: q
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(out)   :: Dq
  integer, intent(in)                                            :: direction
  real                                                           :: tstart

  ! Direction:
  !
  !   3: logically z.
  !   2: logically y.
  !   1: logically x.

  ! Assuming x1, y1, and z1 are the 1D grids, the 3D grids are constructed as:
  ! X = kron( 1_{nz * ny}, x1 );
  ! Y = kron( kron( 1_nz, y1 ), 1_nx );
  ! Z = kron( z1, 1_{ny * nx} );
  !
  ! Thus we can use perfect shuffles to permute the arrays until the desired
  ! leading dimension is at the end (at least this is my hope).

  ! Note that a perfect_shuffle( p, q, x ) for a vector x splits x
  ! into p stacks of q elements each.  So, if x is built as
  ! x = xp kron xq, a pq perfect shuffle returns x = xq kron xp.

  ! Permute the input array into the appropriate coordinate-first indexing.
  tstart = MPI_Wtime()
  if ( direction == 3 ) then
     call perfect_shuffle( nz_loc, nx_loc * ny_loc, q )
  elseif ( direction == 2 ) then
     call perfect_shuffle( nz_loc * ny_loc, nx_loc, q )
     !call perfect_shuffle( nx_loc * ny_loc, nz_loc, q )
  elseif ( direction == 1 ) then
     !call perfect_shuffle( nx_loc, nz_loc * ny_loc, q )
     ! Do nothing, x already leads.
  endif
  timer_shuffle = timer_shuffle + ( MPI_Wtime() - tstart )

  tstart = MPI_Wtime()

  ! Compute the derivative.
  call DGEMM( 'n', 'n', n, r_loc / n, n, 1.0, D, n, q, n, 0.0, Dq, n )

  timer_dgemm = timer_dgemm + ( MPI_Wtime() - tstart )

  ! Permute the input array into the appropriate coordinate-first indexing.
  tstart = MPI_Wtime()
  if ( direction == 3 ) then
     call perfect_shuffle( nx_loc * ny_loc, nz_loc, q )
     call perfect_shuffle( nx_loc * ny_loc, nz_loc, Dq )
  elseif ( direction == 2 ) then
     call perfect_shuffle( nx_loc, nz_loc * ny_loc, q )
     call perfect_shuffle( nx_loc, nz_loc * ny_loc, Dq )
     !call perfect_shuffle( nz_loc, nx_loc * ny_loc, q )
     !call perfect_shuffle( nz_loc, nx_loc * ny_loc, Dq )
  elseif ( direction == 1 ) then
     !call perfect_shuffle( nz_loc * ny_loc, nx_loc, q )
     !call perfect_shuffle( nz_loc * ny_loc, nx_loc, Dq )
     ! Do nothing, x already leads.
  endif
  timer_shuffle = timer_shuffle + ( MPI_Wtime() - tstart )

end subroutine compute_dxi
