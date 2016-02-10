subroutine compute_divergence( qx, qy, qz, divq )
!
! Computes the divergence of vector field (qx, qy, qz).

  use constants
  use legendre
  use mesh_deformation_maps
  use derivatives

  implicit none

  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc), intent(in)  :: qx, qy, qz
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc), intent(out) :: divq
!  real, external :: compute_dx, compute_dy, compute_dz

  divq = compute_dx( qx ) + compute_dy( qy ) + compute_dz( qz )

end subroutine compute_divergence
