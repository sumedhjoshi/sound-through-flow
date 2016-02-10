subroutine quad( n, x, w, ndim )

  implicit none

  ! constants.
  integer, parameter                   :: nn = 1000
  real, parameter                      :: small = 1.0e-30

  ! subroutine parameters.
  integer, intent(in)                  :: n
  integer, intent(in)                  :: ndim
  real, dimension(0:ndim), intent(in)  :: x
  real, dimension(0:ndim), intent(out) :: w

  ! local variables.
  real, dimension(0:nn)                :: al1
  integer                              :: k   ! loop index

  !  **PD-modify**: Also store values of 0 to ndim Legendre
  !  polynomials at all collocation points.
  !  determine the Gauss Quadrature weighting factors

  do k = 0, n
     call legen( al1, n, x(k), nn )

     ! Calculate weighting factor
     w(k) = 2. / (n * (n + 1) * al1(n) * al1(n) + small)
  enddo

  return
end subroutine quad


subroutine legen( al, n, xc, ndim )
  ! ---------------------------------------------------------------------
  ! Calculates values of all Legendre polynomials (and immediately highest
  ! one in hierarchy) at a given collocation point xc.

  implicit none

  ! subroutine parameters.
  integer, intent(in)                  :: ndim
  integer, intent(in)                  :: n
  real, intent(in)                     :: xc

  real, dimension(0:ndim), intent(out) :: al

  ! local variables.
  integer                              :: k      ! loop index
  integer                              :: kp, km ! k plus/minus one, respectively

  al(0)  = 1.
  al(1)  = xc

  do k = 1, n
     kp     = k + 1
     km     = k - 1
     al(kp) = (2. * k + 1.) * xc * al(k) / kp - k * al(km) / kp
  enddo

  return
end subroutine legen

