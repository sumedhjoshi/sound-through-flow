module derivatives

  implicit none

contains

   ! For now, all of these functions assume a diagonal
   ! Jacobian matrix, meaning all the elements are
   ! rectangular.
   function compute_dx( q ) result( q_x )
   ! Computes the x derivative of a grid function.
   !
     use constants
     use legendre
     use mesh_deformation_maps

     implicit none

     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(in)  :: q
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_x
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi1
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi2
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi3

     ! Compute the gradient of the grid function.
     call compute_dxi( q, q_xi1, 1 )
     call compute_dxi( q, q_xi2, 2 )
     call compute_dxi( q, q_xi3, 3 )

     ! Assemble the derivative in x using metric terms.
     q_x = q_xi1 / x_xi1

     ! Enfore C0 continuity of the derivative.
     call enforce_continuity( q_x, 1 )
     !call enforce_continuity( q_x, 2 )
     !call enforce_continuity( q_x, 3 )

   end function compute_dx

   function compute_dy( q ) result( q_y )
   ! Computes the y derivative of a grid function.
   !
     use constants
     use legendre
     use mesh_deformation_maps

     implicit none

     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(in)  :: q
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_y
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi1
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi2
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi3

     ! Compute the gradient of the grid function.
     call compute_dxi( q, q_xi1, 1 )
     call compute_dxi( q, q_xi2, 2 )
     call compute_dxi( q, q_xi3, 3 )

     ! Assemble the derivative in y using metric terms.
     q_y = q_xi2 / y_xi2

     ! Enfore C0 continuity of the derivative.
     !call enforce_continuity( q_y, 1 )
     call enforce_continuity( q_y, 2 )
     !call enforce_continuity( q_y, 3 )

   end function compute_dy

   function compute_dz( q ) result( q_z )
   ! Computes the z derivative of a grid function.
   !
     use constants
     use legendre
     use mesh_deformation_maps

     implicit none

     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ), intent(in)  :: q
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_z
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi1
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi2
     real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc )              :: q_xi3

     ! Compute the gradient of the grid function.
     call compute_dxi( q, q_xi1, 1 )
     call compute_dxi( q, q_xi2, 2 )
     call compute_dxi( q, q_xi3, 3 )

     ! Assemble the derivative in z using metric terms.
     q_z = q_xi3 / z_xi3

     ! Enfore C0 continuity of the derivative.
     !call enforce_continuity( q_z, 1 )
     !call enforce_continuity( q_z, 2 )
     call enforce_continuity( q_z, 3 )

   end function compute_dz

end module derivatives
