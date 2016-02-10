! This function enforces C0 continuity between elements.
! Cq is resulting grid function that is continuous across
! the element interfaces.
subroutine enforce_continuity( q, direction )
! Computes the x derivative of a grid function.
!
  use constants
  use MPI
  use timing

  implicit none

  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(inout) :: q
  integer, intent(in)                                            :: direction
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )                :: Cq
  real, dimension( 1:nx_loc, 1:ny_loc, 1 )                       :: q_left, q_right
  integer                                                        :: ii, ndx
  real                                                           :: tstart

  ! Copy the input buffer into the output buffer.
  Cq = q
  ! Depending on the direction, enforce C0 continuity in different ways.
  select case ( direction )
     case ( 1 ) ! 1 means the logically x direction.

        do ii = 1, mx - 1
           ndx = n * ( ii - 1 ) + n
           Cq( ndx, :, : )     = ( Cq( ndx, :, : ) + Cq( ndx + 1, :, : ) ) / 2.0
           Cq( ndx + 1, :, : ) = Cq( ndx, :, : )
        enddo

     case ( 2 )

        do ii = 1, my - 1
           ndx = n * ( ii - 1 ) + n
           Cq( :, ndx, : ) = ( Cq( :, ndx, : ) + Cq( :, ndx + 1, : ) ) / 2.0
           Cq( :, ndx + 1, : ) = Cq( :, ndx, : )
        enddo

     case ( 3 ) ! 3 means the logically z direction.


        ! Send/receive function values to the left and the right.
        q_left  = q( :, :, 1:1 )
        q_right = q( :, :, nz_loc:nz_loc )
        tstart = MPI_Wtime()
        call sync_ranks( q_left, q_right, nx_loc * ny_loc )
        timer_sync_ranks = timer_sync_ranks + ( MPI_Wtime() - tstart )

        ! Average over the intra-rank interfaces.
        do ii = 1, mz_loc - 1
           ndx = n * ( ii - 1 ) + n
           Cq( :, :, ndx )     = ( Cq( :, :, ndx ) + Cq( :, :, ndx + 1 ) ) / 2.0
           Cq( :, :, ndx + 1 ) = Cq( :, :, ndx )
        enddo

        ! Average over the inter-rank interfaces.
        if ( rank > 0 ) then
           Cq( :, :, 1:1 ) = ( Cq( :, :, 1:1 ) + q_left ) / 2.0
        endif
        if ( rank < nprocs - 1 ) then
           Cq( :, :, nz_loc:nz_loc ) = ( Cq( :, :, nz_loc:nz_loc ) + q_right ) / 2.0
        endif

  end select

  ! Update output.
  q = Cq

end subroutine enforce_continuity
