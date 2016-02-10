subroutine compute_gradient( q, q_x, q_y, q_z )
!
! Computes the gradient of the scalar field v
! and returns components (vx,vz).
!
! This computes the gradient on the local block
! and returns the gradient on just that piece.

  use constants
  use legendre
  use mesh_deformation_maps
  use derivatives
  use mpi

  implicit none

  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(in)  :: q
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(out) :: q_x, q_y, q_z

  ! Compute the gradient.
  q_x = compute_dx( q )
  q_y = compute_dy( q )
  q_z = compute_dz( q )

end subroutine compute_gradient
