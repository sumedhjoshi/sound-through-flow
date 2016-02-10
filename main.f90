program sem_acoustics

  ! Declare modules.
  use constants
  use derivatives
  use geom
  use HDF5
  use fields
  use legendre
  use mesh_deformation_maps
  use options
  use parallel_linear_algebra
  use timing

  !$ use omp_lib, only:        omp_get_max_threads
  use mpi

  implicit none

  ! General use variables.
  real                                              :: pi
  integer                                           :: ii

  ! Declare variables for reading command line inputs.
  integer                                           :: iargc
  integer                                           :: n_arg
  character(len=100)                                :: infile_name
  character(len=64)                                 :: caststr
  character(len=5)                                  :: procstr

  ! HDF5 and MPI status variables.
  integer(kind=4)                                   :: ierr, hdf5_status

  ! Variables for meauring time.
  real                                              :: time_start
  real                                              :: time_end
  real                                              :: solver_start
  real, allocatable, dimension(:)                   :: elapsed_time_cpu
  real                                              :: CPUtimeStart, CPUtimeEnd, CPUtimeElapsed
  integer(kind=8)                                   :: clock_start, clock_stop, clock_rate

  ! Variables for computing CFL.
  real                                              :: iidx, iidt
  real, dimension(1)                                :: max_dt

  ! Initialize MPI.
  call MPI_INIT( ierr )

  ! Get my rank and total number of processes.
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

  ! Initialize the HDF5 library.
  call h5open_f( hdf5_status )

  ! Determine system clock rate.
  call system_clock( COUNT_MAX=clock_rate )
  call system_clock( COUNT_RATE=clock_rate )

  ! Get pi.
  pi = 2.0 * acos( 0.0 )

  ! Write the process number to a string for display.
  write( procstr, '(I3)' ) rank
  procstr            = trim(procstr)

  ! Commence timing.
  if ( rank == root ) call notify( '   ' )

  ! Read command line inputs for input file.
  n_arg = iargc()
  if (n_arg < 1) then

     if ( rank == root ) then
        write(*,*) 'Error: too few input arguments:', n_arg
        write(*,'(a)') '  Usage: sem_acoustics input_file '
     endif

     stop
  else

     ! Get the input and log file names.
     call getarg( 1, infile_name )
     infile_name          = trim( infile_name )

  endif

  ! Present welcome message.
  if ( rank == root ) then
     call notify( '==============================================================================================' )
     call notify( ' ' )
     call notify( '   Spectral Finite Element Sound Propagation through an Incompressible Flow ' )
     call notify( ' ' )
     call notify( '==============================================================================================' )
     call notify( ' ' )
  endif

  ! Read the input file.
  if ( rank == root ) call notify( 'Attempting to read from file ' // infile_name )
  call read_input_file( infile_name )

  ! Configure some of the input variables.
  if ( nprocs > 1 ) then
     ! NOTE: nproc's kind is implicitly set by MPI's interface, so we have to
     !       explicitly cast nsubx so that we're comparing apples to apples.
     if ( mod( int( mz, 4 ), nprocs ) .ne. 0 ) then
        if ( rank == root ) then
           write(*,*) 'Error: incompatible execution parameters: '
           write(*,*) '      nprocs = ', nprocs , ' must divide the number of z-subdomains = ', mz
        endif
        stop
     endif
  endif

  ! Compute some secondary input variables.
  r      = n * n * n * mx * my * mz
  r_loc  = n * n * n * mx * my * mz / nprocs
  mx_loc = mx
  my_loc = my
  mz_loc = mz / nprocs
  nx_loc = n * mx_loc
  ny_loc = n * my_loc
  nz_loc = n * mz_loc

  ! Allocate the spectral differentiation matrices.
  allocate( D(1:n, 1:n), D2(1:n, 1:n), D3(1:n, 1:n) )

  ! Display input information.
  if ( rank == root ) call display_inputs()

  ! Display the estimate for memory that will be required.
  if ( rank == root ) call estimate_memory_footprint()

  ! Generate array with collocation points.
  if ( rank == root ) call notify( 'Generating collocation points.' )
  allocate( points(n), wg(n) )
  call jacobl( n - 1, 0., 0., points, n )

  ! Generate differentiation matrices.
  if ( rank == root ) call notify( 'Generating differentiation matrices.' )
  call derv( n - 1, points, D, D2, D3, n - 1 )

  ! Generate weights for Legendre polynomials and filter matrix.
  if ( rank == root ) call notify( 'Generating quadrature weights.' )
  call quad( n - 1, points, wg, n - 1 )

  ! Allocate the arrays that hold the field variables.
  allocate( ux0(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( ux1(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( uy0(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( uy1(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( uz0(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( uz1(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( s0(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( s1(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( vx(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( vy(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( vz(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( rho(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( beta(1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  allocate( vGvx(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vGvy(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vGvz(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( work1(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( work2(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( work3(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vx_x(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vx_y(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vx_z(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vy_x(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vy_y(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vy_z(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vz_x(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vz_y(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( vz_z(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( beta_x(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( beta_y(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( beta_z(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( rho_x(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( rho_y(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( rho_z(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( As1(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( As2(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( As3(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Aux1(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Aux2(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Aux3(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Auy1(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Auy2(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Auy3(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Auz1(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Auz2(1:nx_loc, 1:ny_loc, 1:nz_loc) )
  allocate( Auz3(1:nx_loc, 1:ny_loc, 1:nz_loc) )

  ! Read the mesh from the mesh file.
  allocate( x( 1:nx_loc, 1:ny_loc, 1:nz_loc ), &
            y( 1:nx_loc, 1:ny_loc, 1:nz_loc ), &
            z( 1:nx_loc, 1:ny_loc, 1:nz_loc ) )
  if ( rank == root ) call notify( 'Reading initial conditions file.' )
  call read_initfile_data( fname_init )
  !call internal_mesher()

  ! Compute the deformation maps.
  if ( rank == root ) call notify( 'Computing the deformation maps.' )
  call setup_deformed_derivatives()

  ! Set some time-stepping parameters.
  nt = floor( t_final / dt ) + 1

  ! Set the Adams-Bashforth coefficients for the first time-step.
  a1 = 1.0
  a2 = 0.0
  a3 = 0.0


!  ! Set the initial conditions.
!  s0 = exp( -( x**2 + y**2 + z**2 ) / 25.0**2 )
!  s1 = exp( -( x**2 + y**2 + z**2 ) / 25.0**2 )
!
!  ! Set the hydrodynamic flow.
!  rho  = 1.02
!  beta = rho * 343.0**2

  ! Compute the CFL constraint for this flow (this only works for cartesian grids so far).
  max_dt = 1.0e9
  do ii = 1, nx_loc - 1
     iidx = abs( x(ii + 1, 1, 1 ) - x(ii, 1, 1) )
     if ( iidx > 1.0e-8 ) then
        iidt = iidx / sqrt( beta( ii + 1, 1, 1) / rho( ii + 1, 1, 1) )
        if ( iidt < max_dt(1) ) max_dt(1) = iidt
     endif
  enddo
  do ii = 1, ny_loc - 1
     iidx = abs( y(1, ii + 1, 1 ) - y(1, ii, 1) )
     if ( iidx > 1.0e-8 ) then
        iidt = iidx / sqrt( beta( 1, ii + 1, 1) / rho( 1, ii + 1, 1) )
        if ( iidt < max_dt(1) ) max_dt(1) = iidt
     endif
  enddo
  do ii = 1, nz_loc - 1
     iidx = abs( z(1, 1, ii + 1 ) - z(1, 1, ii + 1) )
     if ( iidx > 1.0e-8 ) then
        iidt = iidx / sqrt( beta( 1, 1, ii + 1) / rho( 1, 1, ii + 1) )
        if ( iidt < max_dt(1) ) max_dt(1) = iidt
     endif
  enddo
  max_dt(1) = pmaxval( max_dt )
  if ( rank == root ) then

     write( caststr, '(D16.8)' ) max_dt
     call notify( 'Maximum timstep due to CFL constraint: ' // caststr )

  endif

  if ( max_dt(1) < dt ) then
     if ( rank == root ) then
        call notify( ' ' )
        call notify( '   WARNING: the current time-step violates CFL constraint' )
        write( caststr, '(D16.8)' ) dt
        call notify( '      Current time-step   : ' // caststr )
        write( caststr, '(D16.8)' ) CFL
        call notify( '      CFL number requested: ' // caststr )
        write( caststr, '(D16.8)' ) max_dt / CFL
        call notify( '      New time-step       : ' // caststr )
        call notify( ' ' )
     endif

     ! Fix the time-step to satisfy the CFL constraint.
     dt = max_dt(1) / CFL

     ! Set some time-stepping parameters.
     nt = floor( t_final / dt ) + 1

  endif

  ! If the number of time-steps per write exceeds the number of time-steps, set it equal.
  if ( timesteps_between_writes > nt ) timesteps_between_writes = nt

  ! Compute necessary derivatives of the hydrodynamic flow.
  if ( rank == root ) call notify( 'Compute derivatives of the background flow.' )
  call compute_gradient( vx, work1, work2, work3 )
  vGvx = vx * work1 + vy * work2 + vz * work3 ! v * grad( vx )
  call compute_gradient( vy, work1, work2, work3 )
  vGvy = vx * work1 + vy * work2 + vz * work3 ! v * grad( vy )
  call compute_gradient( vz, work1, work2, work3 )
  vGvz = vx * work1 + vy * work2 + vz * work3 ! v * grad( vz )

  ! Compute derivatives of background flow.
  call compute_gradient( rho, rho_x, rho_y, rho_z )
  call compute_gradient( beta, beta_x, beta_y, beta_z )
  call compute_gradient( vx, vx_x, vx_y, vx_z )
  call compute_gradient( vy, vy_x, vy_y, vy_z )
  call compute_gradient( vz, vz_x, vz_y, vz_z )

  ! If asked, check for accuracy of the gradient.
  if ( check_consistency_flag ) then
     call check_consistency()
  endif
  if ( rank == root ) call notify('Done checking consistency.')

  ! Start measuring time.
  time_start = MPI_Wtime()
  call seconds_since_epoch( solver_start )
  call cpu_time( CPUtimeStart )

  ! Open the output field file.
  if ( rank == root ) then
     call notify('Writing to field file: ' // trim( fname_runname ) // '.h5' )
     call open_field_file( trim( fname_runname ) // '.h5' )
  endif

  ! Write the HDF5 header.
  call write_field_header( solver_start )

  ! Start time-stepping.
  write( caststr, '(I10)') nt
  if ( rank == root ) call notify( 'Starting time loop.' )
  if ( rank == root ) call notify( 'Number of time steps: ' // caststr )
  do ii = 2, nt

     ! Start the timer.
     call system_clock( COUNT=clock_start )

     ! Notify the user.
     if ( mod( ii, report_every_n_steps ) == 0 ) then
        write( caststr, '(I10)') ii
        if ( rank == root ) call notify('   Starting time-step ' // caststr )
     endif

     ! Update AB coefficients if necessary.
     if (ii == 3) then
        a1 = 1.5
        a2 = -0.5
        a3 = 0.0
     elseif ( ii .ge. 4 ) then
        a1 = 23.0/12.0
        a2 = -4.0/3.0
        a3 = 5.0/12.0
     endif

     ! Update the field.
     s0  =  s1 + dt * ( a1 * As1  + a2 * As2  + a3 * As3 )
     ux0 = ux1 + dt * ( a1 * Aux1 + a2 * Aux2 + a3 * Aux3 )
     uy0 = uy1 + dt * ( a1 * Auy1 + a2 * Auy2 + a3 * Auy3 )
     uz0 = uz1 + dt * ( a1 * Auz1 + a2 * Auz2 + a3 * Auz3 )

     ! Check to see if the code may be going unstable.
     if ( norm2( s0 ) > big_value ) then
        write(*,*) 'Warning: rank ', rank, ' reports ||s|| of ', norm2( s0 ), '; simulation may be unstable.'
     endif

     ! Time advance.
     s1  = s0
     ux1 = ux0
     uy1 = uy0
     uz1 = uz0
     As3 = As2
     As2 = As1
     Aux3 = Aux2
     Aux2 = Aux1
     Auy3 = Auy2
     Auy2 = Auy1
     Auz3 = Auz2
     Auz2 = Auz1
     call apply_wave_operator( s0, ux0, uy0, uz0, As1, Aux1, Auy1, Auz1 )

     ! If asked write to disk, or if it is the last time-step, write to disk.
     if ( mod( ii, timesteps_between_writes ) == 0 ) then
        call write_field_variables( ii )
     endif

     ! Stop the timer.
     call system_clock( COUNT=clock_stop )

     ! How long did this timestep take?
     if ( rank == root ) then
        call write_field_timestep( ii, real( clock_stop - clock_start ) / clock_rate )
     endif

  enddo

  ! Stop timer.
  time_end = MPI_Wtime()

  ! Report total time and time per time-step.
  if ( rank == root ) then
     call notify( ' ' )
     write( caststr, '(F16.8)' ) time_end - time_start
     call notify( 'Total wall-clock seconds to solution: ' // caststr )
     write( caststr, '(F16.8)' ) ( time_end - time_start)  / nt
     call notify( 'Wall-clock Seconds per timestep     : ' // caststr )
  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  ! Report some timers.

  ! Report timing for rank synchronization.
  write( caststr, '(F16.8)' ) ( time_end - time_start ) / nt
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  if ( rank == root ) then
     call notify( ' ' )
     call notify( 'Time per step: ' )
     call notify( ' ' )
  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  call notify( '   Rank ' // procstr // ': ' // caststr )
  call sleep( 1 )
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  ! Report timing for rank synchronization.
  write( caststr, '(F16.8)' ) timer_sync_ranks / real( nt )
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  if ( rank == root ) then
     call notify( ' ' )
     call notify( 'Time per step spent in communication (rank synchronization): ' )
     call notify( ' ' )
  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  call notify( '   Rank ' // procstr // ': ' // caststr )
  call sleep( 1 )
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  ! Report timing for rank synchronization.
  write( caststr, '(F16.8)' ) timer_dgemm / real( nt )
  if ( rank == root ) then
     call notify( ' ' )
     call notify( 'Time per step spent in DGEMM: ' )
     call notify( ' ' )
  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  call notify( '   Rank ' // procstr // ': ' // caststr )
  call sleep( 1 )
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  ! Report timing for rank synchronization.
  write( caststr, '(F16.8)' ) timer_shuffle / real( nt )
  if ( rank == root ) then
     call notify( ' ' )
     call notify( 'Time per step spent in computing perfect shuffles (data transposition): ' )
     call notify( ' ' )
  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  call notify( '   Rank ' // procstr // ': ' // caststr )
  call sleep( 1 )
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  ! Get the times for all CPUs.
  call cpu_time( CPUtimeEnd )
  CPUtimeElapsed = CPUtimeEnd - CPUtimeStart
  allocate( elapsed_time_cpu( 1:nprocs ) )
  call MPI_GATHER( CPUtimeElapsed,   1, MPI_DOUBLE_PRECISION, &
                   elapsed_time_cpu, 1, MPI_DOUBLE_PRECISION, &
                   root, MPI_COMM_WORLD, ierr )

  ! Report the times into the HDF5 file.
  if ( rank == root ) call write_field_execution_stats( elapsed_time_cpu, time_end - time_start )

  ! Close the field file.
  if ( rank == root ) call close_field_file

  ! Finalize MPI.
  call MPI_FINALIZE( ierr )

end program sem_acoustics
