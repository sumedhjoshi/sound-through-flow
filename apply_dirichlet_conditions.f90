subroutine apply_dirichlet_conditions( Df, f, direction, side )
!
! Applies to a derivative operator Df the homogenous Dirichlet boundary values
! as given by the function f.

  use constants

  implicit none

  real, dimension(1:nx_loc, 1:ny_loc, 1:nz_loc), intent(in)    :: f
  real, dimension(1:nx_loc, 1:ny_loc, 1:nz_loc), intent(inout) :: Df
  integer, intent(in)                                          :: direction, side ! [1,2,3] for direction (x,y,z) and [0,1] for side (lower, upper).

  ! Enforce homogenous Dirichlet boundary conditions in x and y ( these are easy ).
  if ( direction == 2 ) then
     if ( side == 0 ) Df( :, 1, : )      = f(:, 1, : )
     if ( side == 1 ) Df( :, ny_loc, : ) = f(:, ny_loc, : )
  endif

  if ( direction == 1 ) then
     if ( side == 0 ) Df( 1, :, : )      = f(1, :, : )
     if ( side == 1 ) Df( nx_loc, :, : ) = f(nx_loc, :, : )
  endif

  ! Enforce homogenous Dirichlet boundary conditions in z.
  if ( direction == 3 ) then

     if ( side == 0 ) then
        if ( rank == 0 )          Df( :, :, 1 )      = f(:, :, 1 )
     endif

     if ( side == 1 ) then
        if ( rank == nprocs - 1 ) Df( :, :, nz_loc ) = f( :, :, nz_loc )
     endif
  endif

end subroutine apply_dirichlet_conditions
