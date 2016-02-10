subroutine setup_deformed_derivatives
!
! Computes metric terms for mapping derivatives
! of grid functions from physical space into the
! master element.

  use constants
  use geom
  use legendre
  use mesh_deformation_maps

  implicit none

  ! Allocate the mesh deformation maps.
  allocate( x_xi1(1:nx_loc, 1:ny_loc, 1:nz_loc ), x_xi2( 1:nx_loc, 1:ny_loc, 1:nz_loc ), x_xi3( 1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( y_xi1(1:nx_loc, 1:ny_loc, 1:nz_loc ), y_xi2( 1:nx_loc, 1:ny_loc, 1:nz_loc ), y_xi3( 1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( z_xi1(1:nx_loc, 1:ny_loc, 1:nz_loc ), z_xi2( 1:nx_loc, 1:ny_loc, 1:nz_loc ), z_xi3( 1:nx_loc, 1:ny_loc, 1:nz_loc ) )

  ! Compute the metric derivatives.
  call compute_dxi( x, x_xi1, 1 )
  call compute_dxi( x, x_xi2, 2 )
  call compute_dxi( x, x_xi3, 3 )
  call compute_dxi( y, y_xi1, 1 )
  call compute_dxi( y, y_xi2, 2 )
  call compute_dxi( y, y_xi3, 3 )
  call compute_dxi( z, z_xi1, 1 )
  call compute_dxi( z, z_xi2, 2 )
  call compute_dxi( z, z_xi3, 3 )

  !write(*,*) norm2( x_xi1 ), norm2( x_xi2), norm2( x_xi3 )
  !write(*,*) norm2( y_xi1 ), norm2( y_xi2), norm2( y_xi3 )
  !write(*,*) norm2( z_xi1 ), norm2( z_xi2), norm2( z_xi3 )

end subroutine setup_deformed_derivatives
