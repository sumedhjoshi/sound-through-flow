module mesh_deformation_maps

  implicit none

  ! Mesh deformation maps; derivatives of coordinates (x,y,z) relative to (xi1,xi2,xi3).
  real, allocatable, dimension(:,:,:) :: x_xi1, x_xi2, x_xi3
  real, allocatable, dimension(:,:,:) :: y_xi1, y_xi2, y_xi3
  real, allocatable, dimension(:,:,:) :: z_xi1, z_xi2, z_xi3

end module mesh_deformation_maps
