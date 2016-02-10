module geom
!
! Contains the grid and domain geometry.

  implicit none

  real, allocatable, dimension(:,:,:) :: x, y, z
  real, dimension(1:2)                :: Lx, Ly, Lz

end module geom
