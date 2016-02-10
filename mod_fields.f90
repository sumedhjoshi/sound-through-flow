module fields
!
! Contains constants related to the grid.
  implicit none

  ! Acoustic fields at current, previous, and previous-previous times.
  real, allocatable, dimension(:,:,:) :: ux0, ux1
  real, allocatable, dimension(:,:,:) :: uy0, uy1
  real, allocatable, dimension(:,:,:) :: uz0, uz1
  real, allocatable, dimension(:,:,:) :: s0, s1

  ! Wave operator A evaluated at current previous, previous-previous, and previous-previous-previous times.
  real, allocatable, dimension(:,:,:) :: Aux1, Aux2, Aux3
  real, allocatable, dimension(:,:,:) :: Auy1, Auy2, Auy3
  real, allocatable, dimension(:,:,:) :: Auz1, Auz2, Auz3
  real, allocatable, dimension(:,:,:) :: As1, As2, As3

  ! Static hydrodynamic field.
  real, allocatable, dimension(:,:,:) :: vx, vy, vz, rho, beta

  ! Derivatives of the static hydrodynamic field.
  real, allocatable, dimension(:,:,:) :: vGvx, vGvy, vGvz
  real, allocatable, dimension(:,:,:) :: vx_x, vx_y, vx_z !
  real, allocatable, dimension(:,:,:) :: vy_x, vy_y, vy_z !
  real, allocatable, dimension(:,:,:) :: vz_x, vz_y, vz_z ! The tensor gradient of v.

  real, allocatable, dimension(:,:,:) :: rho_x, rho_y, rho_z
  real, allocatable, dimension(:,:,:) :: beta_x, beta_y, beta_z

  ! Work arrays.
  real, allocatable, dimension(:,:,:) :: work1, work2, work3

  ! Temporary array for visualization.
  real, allocatable, dimension(:,:,:) :: viz_array

end module fields
