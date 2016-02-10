subroutine internal_mesher
!
! Builds a cartesian high-order element-based 3D grid.

  use constants
  use geom
  use legendre
  use mesh_deformation_maps

  implicit none

  real :: a, b, w, hx, hy, hz
  real, dimension( 1:nx_loc ) :: x1
  real, dimension( 1:ny_loc ) :: y1
  real, dimension( 1:nz_loc ) :: z1
  real, dimension( 1:n )      :: one_element
  integer                     :: ii, jj, iistart, iiend

  ! Get the domain bounds in z for this rank.
  w = ( Lz(2) - Lz(1) ) / nprocs
  a = Lz(1) + rank * w
  b = Lz(1) + rank * w + w

  ! Get the element sizes in each direction.
  hz = w / mz_loc
  hx = ( Lx(2) - Lx(1) ) / mx_loc
  hy = ( Ly(2) - Ly(1) ) / my_loc

  ! Build one-dimensional versions of the grid functions.
  one_element = ( points + 1 ) * hx / 2.0
  do ii = 1, mx_loc
     iistart = n * ( ii - 1 ) + 1
     iiend   = n * ( ii - 1 ) + n
     x1(iistart:iiend) = Lx(1) + hx * ( ii - 1 ) + one_element
  enddo
  one_element = ( points + 1 ) * hy / 2.0
  do ii = 1, my_loc
     iistart = n * ( ii - 1 ) + 1
     iiend   = n * ( ii - 1 ) + n
     y1(iistart:iiend) = Ly(1) + hy * ( ii - 1 ) + one_element
  enddo
  one_element = ( points + 1 ) * hz / 2.0
  do ii = 1, mz_loc
     iistart = n * ( ii - 1 ) + 1
     iiend   = n * ( ii - 1 ) + n
     !z1(iistart:iiend) = Lz(1) + hz * ( ii - 1 ) + one_element
     z1(iistart:iiend) = a + hz * ( ii - 1 ) + one_element
  enddo

  ! Extrude the 1D grids onto the 3D grid.
  do ii = 1, ny_loc
     do jj = 1, nz_loc
        x( :, ii, jj ) = x1
     enddo
  enddo
  do ii = 1, nz_loc
     do jj = 1, nx_loc
        y( jj, :, ii ) = y1
     enddo
  enddo
  do ii = 1, nx_loc
     do jj = 1, ny_loc
        z( ii, jj, : ) = z1
     enddo
  enddo

end subroutine internal_mesher
