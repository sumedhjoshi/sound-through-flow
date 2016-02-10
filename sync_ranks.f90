subroutine sync_ranks( left, right, nbuffer )
!
! Takes the left/right arrays of dim nbuffer and
! sends the left to the rank on the left and the right
! to the rank on the right.  On exit, left contains information
! from the left rank, and right contains information from the
! right rank.

  use constants, only: rank, nprocs
  use mpi

  implicit none

  integer, intent(in)                       :: nbuffer
  real, dimension(1:nbuffer), intent(inout) :: left
  real, dimension(1:nbuffer), intent(inout) :: right

  real, dimension(1:nbuffer)              :: from_left, from_right
  integer(kind=4)                         :: tag2left = 1, tag2right = 2
  integer(kind=4)                         :: ierr
  integer(kind=4)                         :: req2left_sent, req2left_recv, req2right_sent, req2right_recv
  integer(kind=4), dimension(1:2)         :: allbutN, allbut0
  integer(kind=4), dimension(1:4)         :: middleranks
  integer(kind=4)                         :: nbuffer_4, rank_4

  ! Do some casting.
  nbuffer_4 = int( nbuffer, kind=4 )
  rank_4    = int( rank,    kind=4 )

  ! Post the pass to the right.
  if ( rank > 0 ) then
     call MPI_IRECV( from_left, nbuffer, MPI_DOUBLE_PRECISION, rank - 1, tag2right, MPI_COMM_WORLD, req2right_recv, ierr )
  endif
  if ( rank < nprocs - 1 ) then
     call MPI_ISEND( right, nbuffer, MPI_DOUBLE_PRECISION, rank + 1, tag2right, MPI_COMM_WORLD, req2right_sent, ierr )
  endif

  ! Post the pass to the left.
  if ( rank < nprocs - 1 ) then
     call MPI_IRECV( from_right, nbuffer, MPI_DOUBLE_PRECISION, rank + 1, tag2left, MPI_COMM_WORLD, req2left_recv, ierr )
  endif
  if ( rank > 0 ) then
     call MPI_ISEND( left, nbuffer, MPI_DOUBLE_PRECISION, rank - 1, tag2left, MPI_COMM_WORLD, req2left_sent, ierr )
  endif

  ! Wait on the left/right sends.
  if ( rank == 0 .and. nprocs > 1 ) then
     allbutN(1) = req2right_sent
     allbutN(2) = req2left_recv
     call MPI_WAITALL( int(2, kind=4), allbutN(1:2), MPI_STATUSES_IGNORE, ierr )
  elseif ( rank == nprocs - 1 .and. nprocs > 1 ) then
     allbut0(1) = req2left_sent
     allbut0(2) = req2right_recv
     call MPI_WAITALL( int(2, kind=4), allbut0(1:2), MPI_STATUSES_IGNORE, ierr )
  elseif ( nprocs > 1 ) then
     middleranks(1) = req2right_sent
     middleranks(2) = req2left_recv
     middleranks(3) = req2left_sent
     middleranks(4) = req2right_recv
     call MPI_WAITALL( int(4, kind=4) , middleranks(1:4), MPI_STATUSES_IGNORE, ierr )
  endif

  ! Update the buffers.
  left  = from_left
  right = from_right

end subroutine sync_ranks

