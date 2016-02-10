module legendre
!
! This module defines the Legendre collocation points as well as the
! differentiation matrices associated to those points.

  implicit none

  real, allocatable, dimension(:)   :: points    ! Array with Gauss-Lobatto-Legendre collocation points
  real, allocatable, dimension(:)   :: wg        ! Weights for filtering.
  real, allocatable, dimension(:,:) :: D, D2, D3 ! 1st, 2nd and 3rd order differentiation matrices

end module legendre
