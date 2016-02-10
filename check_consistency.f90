subroutine check_consistency
!
! This routine will check the accurancy of derivatives
! computed by the differentiation matrices on this grid.
! Besides being a good way to debug code, this is also
! pretty much the only way to assess how accurate
! the results are going to be.

! XXX: only works on cartesian grids.

  use constants
  use fields
  use geom
  use parallel_linear_algebra

  implicit none

  real                                            :: pi
  character(len=128)                              :: str_l, str_x, str_y, str_z
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ) :: f, f_x_a, f_y_a, f_z_a
  real, dimension( 1:nx_loc, 1:ny_loc, 1:nz_loc ) :: Dxf, Dyf, Dzf, DivGradf
  real                                            :: len_x, len_y, len_z
  real                                            :: err_x, err_y, err_z
  integer                                         :: lambda

  ! Compute the domain length.
  len_x = Lx(2) - Lx(1)
  len_y = Ly(2) - Ly(1)
  len_z = Lz(2) - Lz(1)

  ! Get pi.
  pi = 2.0 * acos( 0.0 )

  if ( rank == root ) call notify( 'Checking error in gradient (err_x, err_y, err_z)' )

  ! Loop over a few modes.
  do lambda = 1, 10

     ! Set the function and its analytic gradient.
     f     = sin( lambda * pi * x / len_x ) * sin( lambda * pi * y / len_y ) * sin( lambda * pi * z / len_z )
     f_x_a = cos( lambda * pi * x / len_x ) * sin( lambda * pi * y / len_y ) * sin( lambda * pi * z / len_z ) * lambda * pi / len_x
     f_y_a = sin( lambda * pi * x / len_x ) * cos( lambda * pi * y / len_y ) * sin( lambda * pi * z / len_z ) * lambda * pi / len_y
     f_z_a = sin( lambda * pi * x / len_x ) * sin( lambda * pi * y / len_y ) * cos( lambda * pi * z / len_z ) * lambda * pi / len_z

     ! Compute the gradient.
     call compute_gradient( f, Dxf, Dyf, Dzf )

     ! Compute the error.
     err_x = pnorm2( reshape( f_x_a - Dxf, (/ nx_loc*ny_loc*nz_loc /) ) )
     err_y = pnorm2( reshape( f_y_a - Dyf, (/ nx_loc*ny_loc*nz_loc /) ) )
     err_z = pnorm2( reshape( f_z_a - Dzf, (/ nx_loc*ny_loc*nz_loc /) ) )

     ! Notify the user of the error.
     write( str_l, '(I3.0)' ) lambda
     write( str_x, '(D16.4)' ) err_x
     write( str_y, '(D16.4)' ) err_y
     write( str_z, '(D16.4)' ) err_z
     if ( rank == root ) then
        call notify( '   Mode '// trim( str_l ) // ' : '  &
                               // trim( adjustl( str_x ) ) &
                               // ', ' // trim( adjustl( str_y ) ) &
                               // ', ' // trim( adjustl( str_z ) ) )
     endif

     ! Compute the divergence.
     f = x**( 4 + 1 ) + y**( 4 + 1) + z**( 4 + 1)
     call compute_gradient( f, Dxf, Dyf, Dzf )
     call compute_divergence( Dxf, Dyf, Dzf, DivGradf )

     ! Set the exact solution.
     !f = -1.0 * (( lambda * pi / len_x )**2 + ( lambda * pi / len_y )**2 + ( lambda * pi / len_z )**2 ) * &
     !   ( f - x**2 - y**2 - z**2 ) + 3.0 * 2.0 * ( x + y + z )
     f = 4 * (4 + 1) * ( x**(4-1) + y**(4-1) + z**(4-1) )

     ! Compute the error.
     err_x = pnorm2( reshape( f - DivGradf, (/ nx_loc*ny_loc*nz_loc /) ) ) / pnorm2( reshape( f, (/ nx_loc*ny_loc*nz_loc /) ) )

     ! Notify the user of the error.
     write( str_l, '(I3.0)' ) lambda
     write( str_x, '(D16.4)' ) err_x
     if ( rank == root ) then
        call notify( '   Mode '// trim( str_l ) // ' : err divergence : '  &
                               // trim( adjustl( str_x ) ) )
     endif

  enddo
end subroutine check_consistency
