subroutine apply_wave_operator( s, ux, uy, uz, C, Mux, Muy, Muz )
!
! This routine takes four vectors (density, three velocity components)
! and computes the wave operator C(s,ux,uy,uz), M(s,ux,uy,uz).
  use constants
  use derivatives
  use fields

  implicit none

  real, dimension(1:nx_loc, 1:ny_loc, 1:nz_loc), intent(in)  :: s, ux, uy, uz
  real, dimension(1:nx_loc, 1:ny_loc, 1:nz_loc), intent(out) :: C, Mux, Muy, Muz ! The continuity and momentum equations.

  ! The idea is that:
  !
  !   s_new = s_old + dt * C( s_old )
  !   u_new = u_old + dt * M( s_old )

  ! Compute the divergence of the acoustic velocity field.
  call compute_divergence( ux, uy, uz, work1 )

  ! Compute the continuity update.
  C = - ( ux * rho_x + uy * rho_y + uz * rho_z ) / rho - rho * work1 / rho

  ! Compute the gradient of condensation.
  call compute_gradient( s, work1, work2, work3 )

  ! For reference:
  !      If G(vx,vy,vz) is the gradient of the vector field v, then the i-th row of G is
  !      the gradient of the scalar field vi.

  ! Add the density gradient and background nonlinear term.
  if ( rank == root ) call notify('      Computing s * v * grad( v ) ' )
  Mux = -s * vGvx
  Muy = -s * vGvy
  Muz = -s * vGvz

  ! Compute the u * grad(v) in each momentum component.
  if ( rank == root ) call notify('      Computing u * grad( v ).' )
  Mux = Mux - ( ux * vx_x + uy * vx_y + uz * vx_z )
  Muy = Muy - ( ux * vy_x + uy * vy_y + uz * vy_z )
  Muz = Muz - ( ux * vz_x + uy * vz_y + uz * vz_z )

  ! Compute the v * grad(u) term in each momentum component.
  if ( rank == root ) call notify('      Computing v * grad( u ) + laplacian( u ).' )
  call compute_gradient( ux, work1, work2, work3 )
  Mux = Mux - vx * work1 + vy * work2 + vz * work3
  Mux = Mux + ( mu / rho ) * ( compute_dx( work1 ) + compute_dy( work2 ) + compute_dz( work3 ) )
  call compute_gradient( uy, work1, work2, work3 )
  Muy = Muy - vx * work1 + vy * work2 + vz * work3
  Muy = Muy + ( mu / rho ) * ( compute_dx( work1 ) + compute_dy( work2 ) + compute_dz( work3 ) )
  call compute_gradient( uz, work1, work2, work3 )
  Muz = Muz - vx * work1 + vy * work2 + vz * work3
  Muz = Muz + ( mu / rho ) * ( compute_dx( work1 ) + compute_dy( work2 ) + compute_dz( work3 ) )

  ! Compute the ( s * grad(rho) * v ) v / rho term.
  if ( rank == root ) call notify('      Computing ( s * grad( rho ) * v ) v.' )
  Mux = Mux - s * ( rho_x * vx + rho_y * vy + rho_z * vz ) * vx / rho
  Muy = Muy - s * ( rho_x * vx + rho_y * vy + rho_z * vz ) * vy / rho
  Muz = Muz - s * ( rho_x * vx + rho_y * vy + rho_z * vz ) * vz / rho

  ! Compute the ( grad(s) * v ) v term.
  if ( rank == root ) call notify('      Computing ( grad( s ) * v ) v.' )
  call compute_gradient( s, work1, work2, work3 )
  Mux = Mux - ( work1 * vx + work2 * vy + work3 * vz ) * vx
  Muy = Muy - ( work1 * vx + work2 * vy + work3 * vz ) * vy
  Muz = Muz - ( work1 * vx + work2 * vy + work3 * vz ) * vz

  ! Compute the grad( beta * s ) term.
  if ( rank == root ) call notify('      Computing gradient( beta * s ).' )
  Mux = Mux - ( s * beta_x + beta * work1 ) / rho
  Muy = Muy - ( s * beta_y + beta * work2 ) / rho
  Muz = Muz - ( s * beta_z + beta * work3 ) / rho

  ! Apply the boundary conditions.
  if ( rank == root ) call notify('      Computing boundary conditions. ')
  call apply_dirichlet_conditions( C  , 0 * s, 1, 0 )
  call apply_dirichlet_conditions( C  , 0 * s, 2, 0 )
  call apply_dirichlet_conditions( C  , 0 * s, 3, 0 )
  call apply_dirichlet_conditions( C  , 0 * s, 1, 1 )
  call apply_dirichlet_conditions( C  , 0 * s, 2, 1 )
  call apply_dirichlet_conditions( C  , 0 * s, 3, 1 )
  !call apply_dirichlet_conditions( Mux, 0 * ux )
  !call apply_dirichlet_conditions( Muy, 0 * uy )
  !call apply_dirichlet_conditions( Muz, 0 * uz )

end subroutine apply_wave_operator
