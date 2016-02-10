subroutine compute_laplacian( q, Lq )
!
! Computes the Laplacian Lu of the scalar field u.

  use constants
  use legendre
  use derivatives
  use mesh_deformation_maps

  implicit none

  ! Input variables.
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(out) :: Lq
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(in)  :: q
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_x, q_y, q_z

  ! Compute the divergence of the gradient.
  !call compute_gradient( q, q_x, q_y, q_z )
  !call compute_divergence( q_x, q_y, q_z, Lq )

  ! Compute the Laplacian.
  Lq = compute_dx( compute_dx( q ) ) + &
       compute_dy( compute_dy( q ) ) + &
       compute_dz( compute_dz( q ) )

end subroutine compute_laplacian
