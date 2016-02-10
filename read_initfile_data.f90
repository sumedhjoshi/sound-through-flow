subroutine read_initfile_data( fname_init )
! Reads and parses an _init file, which sets the initial values of all the
! field variables.

  use constants, only:         n, mx, my, mz, rank, r_loc, root
  use fields
  use geom, only:              x, y, z
  use HDF5
  use mpi

  implicit none

  character(len=100), intent(in) :: fname_init

  ! Set some HDF5 data-set and file identifiers.
  integer(hid_t) :: file_id
  integer(hid_t) :: grid_x_id
  integer(hid_t) :: grid_y_id
  integer(hid_t) :: grid_z_id
  integer(hid_t) :: ic_ux_id
  integer(hid_t) :: ic_uy_id
  integer(hid_t) :: ic_uz_id
  integer(hid_t) :: ic_s_id
  integer(hid_t) :: env_rho_id
  integer(hid_t) :: env_beta_id
  integer(hid_t) :: env_vx_id
  integer(hid_t) :: env_vy_id
  integer(hid_t) :: env_vz_id

  ! Set the dataspace dimension.
  integer(hsize_t), dimension(3) :: data_dimension

  integer(kind=4) :: error_flag

  ! Set up a temporary buffer to read the data into.
  real, allocatable, dimension(:,:,:) :: read_buffer

  ! Set the data dimensions.
  data_dimension(1) = n * mx
  data_dimension(2) = n * my
  data_dimension(3) = n * mz

  ! Allocate a buffer.
  if ( rank == root ) then
     allocate( read_buffer( 1:n * mx, 1:n * my, 1:n * mz ) )
  endif

  ! Open the file.
  call h5fopen_f( fname_init, H5F_ACC_RDWR_F, file_id, error_flag )

  ! Open all the data-sets we'll need.
  call h5dopen_f( file_id, '/grid/x', grid_x_id, error_flag )
  call h5dopen_f( file_id, '/grid/y', grid_y_id, error_flag )
  call h5dopen_f( file_id, '/grid/z', grid_z_id, error_flag )
  call h5dopen_f( file_id, '/ic/ux',  ic_ux_id,  error_flag )
  call h5dopen_f( file_id, '/ic/uy',  ic_uy_id,  error_flag )
  call h5dopen_f( file_id, '/ic/uz',  ic_uz_id,  error_flag )
  call h5dopen_f( file_id, '/ic/s',   ic_s_id,   error_flag )
  call h5dopen_f( file_id, '/environment/rho',   env_rho_id,  error_flag )
  call h5dopen_f( file_id, '/environment/beta',  env_beta_id, error_flag )
  call h5dopen_f( file_id, '/environment/vx',    env_vx_id,   error_flag )
  call h5dopen_f( file_id, '/environment/vy',    env_vy_id,   error_flag )
  call h5dopen_f( file_id, '/environment/vz',    env_vz_id,   error_flag )

  ! Root reads each array and scatters to all processes..

  ! Read the grid.
  if ( rank == root ) call h5dread_f( grid_x_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    x, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )

  if ( rank == root ) call h5dread_f( grid_y_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    y, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )

  if ( rank == root ) call h5dread_f( grid_z_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    z, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )

  ! Read the initial conditions.
  if ( rank == root ) call h5dread_f( ic_ux_id,  H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    ux0, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )
  ux1 = ux0

  if ( rank == root ) call h5dread_f( ic_uy_id,  H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    uy0, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )
  uy1 = uy0

  if ( rank == root ) call h5dread_f( ic_uz_id,  H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    uz0, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )
  uz1 = uz0

  if ( rank == root ) call h5dread_f( ic_s_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    s0, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )
  s1 = s0

  ! Read the environment.
  if ( rank == root ) call h5dread_f( env_rho_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    rho, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )

  if ( rank == root ) call h5dread_f( env_beta_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    beta, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )

  if ( rank == root ) call h5dread_f( env_vx_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    vx, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )

  if ( rank == root ) call h5dread_f( env_vy_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    vy, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )

  if ( rank == root ) call h5dread_f( env_vz_id, H5T_IEEE_F64LE, read_buffer, data_dimension, error_flag )
  call MPI_SCATTER( read_buffer, r_loc, MPI_DOUBLE_PRECISION, &
                    vz, r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_flag )

  ! Close all the data-sets we just read.
  call h5dclose_f( grid_x_id, error_flag )
  call h5dclose_f( grid_y_id, error_flag )
  call h5dclose_f( grid_z_id, error_flag )
  call h5dclose_f( ic_ux_id, error_flag )
  call h5dclose_f( ic_uy_id, error_flag )
  call h5dclose_f( ic_uz_id, error_flag )
  call h5dclose_f( ic_s_id, error_flag )
  call h5dclose_f( env_rho_id, error_flag )
  call h5dclose_f( env_beta_id, error_flag )
  call h5dclose_f( env_vx_id, error_flag )
  call h5dclose_f( env_vy_id, error_flag )
  call h5dclose_f( env_vz_id, error_flag )

  ! Close the file.
  call h5fclose_f( file_id, error_flag )

end subroutine read_initfile_data
