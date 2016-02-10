module parallel_linear_algebra
!
! Contains functions for doing parallel linear algebra on distributed arrays.

  implicit none

contains

  ! Compute a parallel distributed Euclidean norm.
  function pnorm2( x ) result( normx )

     use mpi, only: MPI_SUM, MPI_COMM_WORLD, MPI_IN_PLACE, MPI_DOUBLE_PRECISION

     implicit none

     real, dimension(:), intent(in) :: x
     real                           :: normx
     integer                        :: ierr
     real, external                 :: DDOT

     normx = DDOT( size(x,1), x, 1, x, 1 )
     call MPI_ALLREDUCE( MPI_IN_PLACE, normx, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     normx = sqrt( normx )

  end function pnorm2

  ! Compute a parallel distributed dot product.
  function pdot_product( x, y ) result( xdoty )

     use mpi, only: MPI_SUM, MPI_COMM_WORLD, MPI_IN_PLACE, MPI_DOUBLE_PRECISION

     implicit none

     real, dimension(:), intent(in) :: x, y
     real                           :: xdoty
     real                           :: iixdoty
     integer                        :: ierr
     real, external                 :: DDOT

     iixdoty = DDOT( size(x,1), x, 1, y, 1 )
     call MPI_ALLREDUCE( iixdoty, xdoty, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

  end function pdot_product

  ! Parallel distributed assignment.
  function assign_value_at( x, a, z, ndx, val )

     implicit none

     integer, intent(in)                     :: a, z      ! ownership range of this rank.
     real, dimension(1:z-a+1), intent(inout) :: x         ! local block of the distributed array.
     integer, intent(in)                   :: ndx       ! index in the full array we want to assign a value to.
     real, intent(in)                      :: val       ! value we want to assign.
     real, dimension(1:z-a+1)  :: assign_value_at
     integer                   :: local_ndx

     ! If the modification is in this rank's range, do it, otherwise return the same array.
     if ( ndx >= a .and. ndx <= z ) then

        ! Get the local index and assign.
        local_ndx = ndx - a + 1
        x(local_ndx) = val

     endif

     ! Return.
     assign_value_at = x

  end function assign_value_at

  ! Compute a parallel distributed MAX.
  function pmaxval( x ) result( xmax )

     use mpi, only: MPI_MAX, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION

     implicit none

     real, dimension(:), intent(in) :: x
     real                           :: xmax
     real                           :: iixmax
     integer                        :: ierr

     iixmax = maxval( x )
     xmax = 0.0
     call MPI_ALLREDUCE( iixmax, xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )

  end function pmaxval

end module parallel_linear_algebra
