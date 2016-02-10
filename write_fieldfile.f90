! NOTE: The calls to h5fflush_f() do not guarantee correct operation in the
!       case of a single writer, multiple reader scenario - it merely reduces
!       the window where the readers will see something different than they
!       expect.
!
!       THIS IS JUST A BANDAID ATTEMPTING TO HIDE THE FACT THAT FUNDAMENTALLY
!       HDF5 CANNOT HANDLE A SINGLE WRITER WITH MULTIPLE READERS IN ALL CASES.
!
!       It is quite possible that the reader can see an invalid view of the
!       file as it is being written, which can be seen when the time per
!       timestep is small and the number of timesteps between writes is
!       tiny.  In very limited testing, the invalid view was simply that
!       the /field/ group appeared to not contain any sub-groups, nor the
!       timestep count dataset.
!
!       The Single Writer, Multiple Reader (SWMR) interface being developed
!       by the HDF5 group (as of spring 2015) will not solve this problem.
!       According to section 3 of the most recent version of their user's
!       guide:
!
!     https://www.hdfgroup.org/HDF5/docNewFeatures/UG-HDF5-SWMR-20130629-v3.pdf
!
!       Groups and datasets cannot be added to a file accessed for SWMR.
!       Currently, the field file is written such that each timestep written
!       lives within its own group.

subroutine open_field_file( field_file_name )
! Opens the field file for writing throughout the solver's life time.  If
! the specified field file already exists, it will be overwritten.  Nothing
! is written to the field file.
!
! This routine opens both the scalar and grid data spaces.

  use constants, only: n, mx, my, mz
  use HDF5, only:      H5F_ACC_TRUNC_F, H5S_SCALAR_F, &
                       hid_t, hsize_t, &
                       h5fcreate_f, &
                       h5screate_f, h5screate_simple_f
  use io, only:        field_file_id, &
                       field_grid_dataspace_id, field_scalar_dataspace_id

  implicit none

  character(len=*), intent(in)        :: field_file_name

  integer(kind=4)                     :: data_rank
  integer(hsize_t), dimension(3)      :: data_dimensions
  integer(hid_t)                      :: hdf5_status

  ! create our output file, overwriting any existing one.
  call h5fcreate_f( trim (field_file_name ), &
                    H5F_ACC_TRUNC_F, field_file_id, hdf5_status )

  ! create a scalar dataspace that we'll use to create scalar datatypes
  ! with.  each of the non-grid variables are tied to this dataspace.
  call h5screate_f( H5S_SCALAR_F, field_scalar_dataspace_id, hdf5_status )

  ! create a 3D dataspace that we'll use for grid-interface variables
  ! (coordinates, velocities, density, etc).
  data_rank          = 3
  data_dimensions(1) = n * mx
  data_dimensions(2) = n * my
  data_dimensions(3) = n * mz
  call h5screate_simple_f( data_rank, data_dimensions, field_grid_dataspace_id, &
                           hdf5_status )

end subroutine open_field_file

subroutine write_field_header( solver_start )
! Writes the field file's header to disk.  The solver's configuration is
! collected onto the root rank and written to the output file.  This routine
! must be called before write_field_timestamp() or write_field_variables(), and
! after open_field_file() is called.
!
! This routine creates the number timesteps dataset identifier.

  use constants, only:               r, r_loc, rank, root
  use geom, only:                    x, y, z
  use HDF5, only:                    H5F_SCOPE_LOCAL_F, H5T_STD_I32LE, &
                                     hid_t, &
                                     h5dcreate_f, &
                                     h5fflush_f, &
                                     h5gclose_f, h5gcreate_f
  use io, only:                      field_file_id, &
                                     field_scalar_dataspace_id, &
                                     field_name_field, field_name_number_steps, &
                                     field_number_steps_dataset_id, field_scalar_dataspace_id
  use mpi, only:                     MPI_COMM_WORLD, MPI_DOUBLE_PRECISION

  implicit none

  real, intent(in)                :: solver_start

  integer(kind=4)                 :: ierr
  real, allocatable, dimension(:) :: x4disk, y4disk, z4disk

  ! HDF5 identifier and status code.
  integer(hid_t)                  :: group_id
  integer(hid_t)                  :: hdf5_status

  ! allocate space for a consolidated grid.
  allocate( x4disk(1:r), y4disk(1:r), z4disk(1:r) )

  ! pull the entirety of the grid onto the root rank for writing.
  call MPI_GATHER( x, r_loc, MPI_DOUBLE_PRECISION, &
                   x4disk,     r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHER( y, r_loc, MPI_DOUBLE_PRECISION, &
                   y4disk,     r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHER( z, r_loc, MPI_DOUBLE_PRECISION, &
                   z4disk,     r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )

  ! permute the grid into x/y/z-first indexing for Paraview.
  ! we need to swap the x and z coordinate to make this happen.
  ! I don't think this is a perfect shuffle, though.

  if ( rank == root ) then
     ! write out the configuration.
     call write_field_header_config

     ! write out the grid.
     call write_field_header_grid( x4disk, y4disk, z4disk )

     ! write out the execution configuration.
     call write_field_header_execution_config( solver_start )

     ! create a group to hold our time-dependent field variables.
     call h5gcreate_f( field_file_id, field_name_field, group_id, hdf5_status )
     call h5gclose_f( group_id, hdf5_status )

     ! create a scalar for the number of timesteps we have written out.  this
     ! wil be updated each time an additional step is written out to disk.
     call h5dcreate_f( field_file_id, field_name_number_steps, H5T_STD_I32LE, &
                       field_scalar_dataspace_id, field_number_steps_dataset_id, hdf5_status )

     ! flush the contents of the file as they exist so that others see the
     ! current contents.
     call h5fflush_f( field_file_id, H5F_SCOPE_LOCAL_F, hdf5_status )

  endif

  deallocate( x4disk, y4disk, z4disk )

end subroutine write_field_header

subroutine write_field_header_config
! Writes out the solver's configuration to the field file's header.

  use constants, only:         dt, mu, t_final
  use HDF5, only:              H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                               hid_t, hsize_t, &
                               h5dclose_f, h5dcreate_f, h5dwrite_f, &
                               h5gclose_f, h5gcreate_f, &
                               h5sclose_f, h5screate_f, h5screate_simple_f
  use io, only:                field_file_id, &
                               field_scalar_dataspace_id, &
                               field_name_config, &
                               field_name_mu, &
                               field_name_dt, field_name_tend, field_name_timesteps_per_write
  use options, only:           timesteps_between_writes

  implicit none

  integer(hsize_t), dimension(2) :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t) :: dataset_id
  integer(hid_t) :: group_id
  integer(hid_t) :: hdf5_status

  ! create the group to hold our configuration variables.
  call h5gcreate_f( field_file_id, field_name_config, group_id, hdf5_status )
  call h5gclose_f( group_id, hdf5_status )

  ! create a scalar for the time step parameter, dt.
  call h5dcreate_f( field_file_id, field_name_dt, H5T_IEEE_F64LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, dt, data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar for the time step threshold, tend.
  call h5dcreate_f( field_file_id, field_name_tend, H5T_IEEE_F64LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, t_final, data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar for the number of time steps per field update,
  ! timesteps_per_write.
  call h5dcreate_f( field_file_id, field_name_timesteps_per_write, H5T_STD_I32LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_INTEGER, int( timesteps_between_writes, 4 ), data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar for the viscosity, mu.
  call h5dcreate_f( field_file_id, field_name_mu, H5T_IEEE_F64LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, mu, data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! write out the GMRES configuration parameters.
  !call write_field_header_gmres

end subroutine write_field_header_config

subroutine write_field_header_grid( full_cx, full_cy, full_cz )
! Writes out the solver's grid to the field file's header.

  use constants, only:   n, r, mx, my, mz
  use HDF5, only:        H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                         hid_t, hsize_t, &
                         h5dclose_f, h5dcreate_f, h5dwrite_f, &
                         h5gclose_f, h5gcreate_f, &
                         h5screate_f, h5screate_simple_f
  use io, only:          field_file_id, &
                         field_grid_dataspace_id, field_scalar_dataspace_id, &
                         field_name_grid, field_name_n, field_name_mx, field_name_my, &
                         field_name_mz, field_name_x, field_name_y, field_name_z

  implicit none

  ! XXX: these shouldn't be passed in, but collected when we migrate to the
  !      parallel interface.
  real, dimension(1:r), intent(in)  :: full_cx
  real, dimension(1:r), intent(in)  :: full_cy
  real, dimension(1:r), intent(in)  :: full_cz

  integer(hsize_t), dimension(3)      :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                      :: dataset_id
  integer(hid_t)                      :: group_id
  integer(hid_t)                      :: hdf5_status

  ! specify the dimensions of our grid variables, full_cx and full_cz.  note
  ! that these are ignored when specified with our scalar variables.
  data_dimensions(1) = n * mx
  data_dimensions(2) = n * my
  data_dimensions(3) = n * mz

  ! create the grid group that holds all of the variables we're
  ! writing out.
  call h5gcreate_f( field_file_id, field_name_grid, group_id, hdf5_status )
  call h5gclose_f( group_id, hdf5_status )

  ! create a scalar for the collocation dimension, n.
  call h5dcreate_f( field_file_id, field_name_n, H5T_STD_I32LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_INTEGER, int( n, 4 ), data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar for the number of x sub-domains, mx.
  call h5dcreate_f( field_file_id, field_name_mx, H5T_STD_I32LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_INTEGER, int( mx, 4 ), data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar grid variable for each of the grid's x positions,
  ! cx4disk.
  call h5dcreate_f( field_file_id, field_name_x, H5T_IEEE_F64LE, &
                    field_grid_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, full_cx, data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar for the number of y sub-domains, my.
  call h5dcreate_f( field_file_id, field_name_my, H5T_STD_I32LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_INTEGER, int( my, 4 ), data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar grid variable for each of the grid's y positions,
  ! cz4disk.
  call h5dcreate_f( field_file_id, field_name_y, H5T_IEEE_F64LE, &
                    field_grid_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, full_cy, data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar for the number of z sub-domains, mz.
  call h5dcreate_f( field_file_id, field_name_mz, H5T_STD_I32LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_INTEGER, int( mz, 4 ), data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar grid variable for each of the grid's z positions,
  ! cz4disk.
  call h5dcreate_f( field_file_id, field_name_z, H5T_IEEE_F64LE, &
                    field_grid_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, full_cz, data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

end subroutine write_field_header_grid

! XXX: change the other elapsed_*_time variables
subroutine write_field_header_execution_config( solver_start )
! Writes out the solver's execution configuration to the field file's header.
!
! The wall time step dataset identifier is created in this routine.

  use constants, only:   nprocs
  use HDF5, only:        H5P_DATASET_CREATE_F, H5S_UNLIMITED_F, &
                         H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32LE, &
                         hid_t, hsize_t, &
                         h5dclose_f, h5dcreate_f, h5dwrite_f, &
                         h5gclose_f, h5gcreate_f, &
                         h5pclose_f, h5pcreate_f, h5pset_chunk_f, &
                         h5sclose_f, h5screate_f, h5screate_simple_f
  use io, only:          field_file_id, &
                         field_length1_vector_dataspace_id, field_scalar_dataspace_id, &
                         field_name_execution, field_name_number_ranks, field_name_number_threads, &
                         field_name_cpu_time, field_name_wall_time, field_name_start_time, &
                         field_name_total_wall_time, &
                         field_name_step_wall_time, &
                         field_wall_time_step_dataset_id
!$ use omp_lib, only:    omp_get_max_threads
  use options, only:     timesteps_between_writes

  implicit none

  real, intent(in)               :: solver_start

  ! number of OpenMP threads.
  integer(kind=4)                :: number_threads

  ! rank and dimension vectors for the time step vector.
  ! XXX: is this the right type for data_rank?
  integer(kind=4)                :: data_rank
  integer(hsize_t), dimension(1) :: current_dimensions
  integer(hsize_t), dimension(1) :: maximum_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                 :: dataset_id
  integer(hid_t)                 :: chunked_property_id
  integer(hid_t)                 :: extensible_dataspace_id
  integer(hid_t)                 :: group_id
  integer(hid_t)                 :: hdf5_status

  ! identify the maximum number of threads OpenMP can execute with.
  number_threads    = 1
  !$ number_threads = omp_get_max_threads()

  ! create groups to hold our execution measurements.
  call h5gcreate_f( field_file_id, field_name_execution, group_id, hdf5_status )
  call h5gclose_f( group_id, hdf5_status )
  call h5gcreate_f( field_file_id, field_name_wall_time, group_id, hdf5_status )
  call h5gclose_f( group_id, hdf5_status )

  ! create a scalar for the number of MPI ranks, nprocs.
  call h5dcreate_f( field_file_id, field_name_number_ranks, H5T_STD_I32LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_INTEGER, int( nprocs, 4 ), current_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a scalar for the number of OpenMP threads, number_threads.
  call h5dcreate_f( field_file_id, field_name_number_threads, H5T_STD_I32LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_INTEGER, number_threads, current_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! create a 1D dataspace with a fixed length of 1 that can grow to
  ! arbitrary lengths.
  data_rank             = 1
  current_dimensions(1) = 1
  maximum_dimensions(1) = H5S_UNLIMITED_F
  call h5screate_simple_f( data_rank, current_dimensions, extensible_dataspace_id, &
                           hdf5_status, maximum_dimensions )

  ! create a property that specifies our chunk size for the 1D vector
  ! dataspace.  we size each chunk to correspond to the number of timesteps
  ! between field updates (this decision is fairly arbitrary).
  data_rank             = 1
  current_dimensions(1) = timesteps_between_writes
  call h5pcreate_f( H5P_DATASET_CREATE_F, chunked_property_id, hdf5_status )
  call h5pset_chunk_f( chunked_property_id, data_rank, current_dimensions, hdf5_status )

  ! create datasets with the 1D dataspace that we will extend in the future
  ! (time, GMRES iterations, and numeric errors per step).  these datasets
  ! have to be chunked to support extension.
!  call h5dcreate_f( field_file_id, field_name_gmres_diffusion_iters, H5T_STD_I32LE, &
!                    extensible_dataspace_id, field_gmres_diffusion_iterations_dataset_id, hdf5_status, &
!                    chunked_property_id )
!  call h5dcreate_f( field_file_id, field_name_gmres_poisson_iters, H5T_STD_I32LE, &
!                    extensible_dataspace_id, field_gmres_poisson_iterations_dataset_id, hdf5_status, &
!                    chunked_property_id )
!  call h5dcreate_f( field_file_id, field_name_gmres_viscous_iters, H5T_STD_I32LE, &
!                    extensible_dataspace_id, field_gmres_viscous_iterations_dataset_id, hdf5_status, &
!                    chunked_property_id )
!  call h5dcreate_f( field_file_id, field_name_l2_poisson_error, H5T_IEEE_F64LE, &
!                    extensible_dataspace_id, field_l2_poisson_error_dataset_id, hdf5_status, &
!                    chunked_property_id )
!  call h5dcreate_f( field_file_id, field_name_l2_poisson_schur_error, H5T_IEEE_F64LE, &
!                    extensible_dataspace_id, field_l2_poisson_schur_error_dataset_id, hdf5_status, &
!                    chunked_property_id )
!  call h5dcreate_f( field_file_id, field_name_linf_diffusion_error, H5T_IEEE_F64LE, &
!                    extensible_dataspace_id, field_linf_diffusion_error_dataset_id, hdf5_status, &
!                    chunked_property_id )
!  call h5dcreate_f( field_file_id, field_name_linf_divergence_error, H5T_IEEE_F64LE, &
!                    extensible_dataspace_id, field_linf_divergence_error_dataset_id, hdf5_status, &
!                    chunked_property_id )
!  call h5dcreate_f( field_file_id, field_name_linf_viscous_x_error, H5T_IEEE_F64LE, &
!                    extensible_dataspace_id, field_linf_viscous_x_error_dataset_id, hdf5_status, &
!                    chunked_property_id )
!  call h5dcreate_f( field_file_id, field_name_linf_viscous_z_error, H5T_IEEE_F64LE, &
!                    extensible_dataspace_id, field_linf_viscous_z_error_dataset_id, hdf5_status, &
!                    chunked_property_id )
  call h5dcreate_f( field_file_id, field_name_step_wall_time, H5T_IEEE_F64LE, &
                    extensible_dataspace_id, field_wall_time_step_dataset_id, hdf5_status, &
                    chunked_property_id )

  ! create a 1x1 memory space used to extend the 1D dataspace.
  data_rank             = 1
  current_dimensions(1) = 1
  call h5screate_simple_f( data_rank, current_dimensions, field_length1_vector_dataspace_id, hdf5_status )

  ! close the extensible dataspace.
  call h5pclose_f( chunked_property_id, hdf5_status )
  call h5sclose_f( extensible_dataspace_id, hdf5_status )

  ! create a scalar for the simulation's start time, in unix seconds.
  call h5dcreate_f( field_file_id, field_name_start_time, H5T_IEEE_F64LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, solver_start, current_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

end subroutine write_field_header_execution_config

subroutine write_field_timestep( timestep_index, elapsed_time )
! Writes out timing information for a single timestep.

  use HDF5, only:    H5S_SELECT_SET_F, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5F_SCOPE_LOCAL_F, &
                     hid_t, hsize_t, &
                     h5dget_space_f, h5dset_extent_f, h5dwrite_f, &
                     h5fflush_f, &
                     h5sselect_hyperslab_f
  use io, only:      field_length1_vector_dataspace_id, &
                     field_wall_time_step_dataset_id
  use options, only: timesteps_between_writes

  implicit none

  integer, intent(in)            :: timestep_index
  real, intent(in)               :: elapsed_time

  ! HDF5 identifiers and status code.
  integer(hsize_t), dimension(1) :: data_dimensions
  integer(hsize_t), dimension(1) :: count, offset
  integer(hid_t)                 :: extensible_dataspace_id
  integer(hid_t)                 :: hdf5_status

  ! extend the datasets' extent by one element to include this timestep.
  data_dimensions(1) = timestep_index
  call h5dset_extent_f( field_wall_time_step_dataset_id, data_dimensions, hdf5_status )

  ! set the dataspace to a single scalar after the last timestep, and write
  ! out the new values.
  offset(1) = timestep_index - 1
  count(1)  = 1

  ! timestep wall time.
  call h5dget_space_f( field_wall_time_step_dataset_id, extensible_dataspace_id, hdf5_status )
  call h5sselect_hyperslab_f( extensible_dataspace_id, H5S_SELECT_SET_F, &
                              offset, count, hdf5_status )
  call h5dwrite_f( field_wall_time_step_dataset_id, H5T_NATIVE_DOUBLE, elapsed_time, &
                   count, hdf5_status, field_length1_vector_dataspace_id, extensible_dataspace_id )

  ! periodically, flush the updates to disk.  we only do this when we're going
  ! to write field variables so we don't unnecessarily introduce overhead
  ! without much gain.  while individual step timings/error are useful, they
  ! don't have to be available immediately.
  if (mod( timestep_index, timesteps_between_writes ) == 0) then
     call h5fflush_f( field_wall_time_step_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )
  end if

end subroutine write_field_timestep

subroutine write_field_variables( timestep_index )
! Writes out the field variables ux, uz, and rho, for a single time.  The
! field is collected onto the root rank and written to the output file.  This
! can only be called after write_field_header() has been called.

  use constants, only:          r, r_loc, rank, root
  use fields, only:             s0, ux0, uy0, uz0
  use HDF5, only:               H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_IEEE_F64LE, H5F_SCOPE_LOCAL_F, &
                                hid_t, hsize_t, &
                                h5dclose_f, h5dcreate_f, h5dwrite_f, &
                                h5fflush_f, &
                                h5gclose_f, h5gcreate_f
  use io, only:                 field_file_id, &
                                field_grid_dataspace_id, field_number_steps_dataset_id, &
                                field_name_step_base, field_name_s, field_name_ux, field_name_uy, field_name_uz
  use mpi, only:                MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use options, only:            timesteps_between_writes

  implicit none

  integer, intent(in)             :: timestep_index

  integer(kind=4)                 :: ierr
  real, allocatable, dimension(:) :: ux4disk, uy4disk, uz4disk, s4disk
  character(len=8)                :: step_number_string !XXX: review this name

  ! dummy dimension variable for writing datasets.  all of the configuration
  ! values written are scalars, though a vector of dimensions must be supplied
  ! to each h5dwrite_f() call.  note that this value is ignored by the HDF5
  ! library.
  integer(hsize_t), dimension(1)  :: data_dimensions = (/ 1 /)

  ! HDF5 identifiers and status code.
  integer(hid_t)                  :: dataset_id
  integer(hid_t)                  :: group_id
  integer(hid_t)                  :: hdf5_status

  ! Allocate some arrays for file I/O.
  allocate( ux4disk(1:r), uy4disk(1:r), uz4disk(1:r), s4disk(1:r) )

  call MPI_GATHER( ux0, r_loc, MPI_DOUBLE_PRECISION, &
                   ux4disk,      r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHER( uy0, r_loc, MPI_DOUBLE_PRECISION, &
                   uy4disk,      r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHER( uz0, r_loc, MPI_DOUBLE_PRECISION, &
                   uz4disk,      r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHER( s0,  r_loc, MPI_DOUBLE_PRECISION, &
                   s4disk,       r_loc, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr )

  if ( rank == root ) then

     ! create a group that holds just this timestep's field variables.
     write( step_number_string, '(I0)' ) timestep_index / timesteps_between_writes
     call h5gcreate_f( field_file_id, trim( field_name_step_base ) // step_number_string, &
                       group_id, hdf5_status )

     ! create a scalar grid variable for each of the grid's x
     ! velocities, ux4disk.
     call h5dcreate_f( group_id, field_name_ux, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call h5dwrite_f( dataset_id, &
                      H5T_NATIVE_DOUBLE, ux4disk, data_dimensions, &
                      hdf5_status )
     call h5dclose_f( dataset_id, hdf5_status )

     ! create a scalar grid variable for each of the grid's y
     ! velocities, uy4disk.
     call h5dcreate_f( group_id, field_name_uy, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call h5dwrite_f( dataset_id, &
                      H5T_NATIVE_DOUBLE, uy4disk, data_dimensions, &
                      hdf5_status )
     call h5dclose_f( dataset_id, hdf5_status )

     ! create a scalar grid variable for each of the grid's z
     ! velocities, uz4disk.
     call h5dcreate_f( group_id, field_name_uz, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call h5dwrite_f( dataset_id, &
                      H5T_NATIVE_DOUBLE, uz4disk, data_dimensions, &
                      hdf5_status )
     call h5dclose_f( dataset_id, hdf5_status )

     ! create a scalar grid variable for each of the grid's densities
     ! rho4disk.
     call h5dcreate_f( group_id, field_name_s, H5T_IEEE_F64LE, &
                       field_grid_dataspace_id, dataset_id, hdf5_status )
     call h5dwrite_f( dataset_id, &
                       H5T_NATIVE_DOUBLE, s4disk, data_dimensions, &
                       hdf5_status )
     call h5dclose_f( dataset_id, hdf5_status )

     ! flush the step we just constructed to disk.
     call h5fflush_f( group_id, H5F_SCOPE_LOCAL_F, hdf5_status )

     ! close this timestep's group and update the number of timesteps in
     ! the output file.
     call h5gclose_f( group_id, hdf5_status )
     call h5dwrite_f( field_number_steps_dataset_id, &
                      H5T_NATIVE_INTEGER, int( timestep_index / timesteps_between_writes, 4 ), data_dimensions, &
                      hdf5_status )

     ! flush the number of steps we've written to disk.
     call h5fflush_f( field_number_steps_dataset_id, H5F_SCOPE_LOCAL_F, hdf5_status )

  endif

  deallocate( ux4disk, uy4disk, uz4disk, s4disk )

end subroutine write_field_variables

subroutine write_field_execution_stats( elapsed_cpu_time, elapsed_total_time )
! Writes out the solver's execution statistics.
!
! XXX: pull the gather of CPU time into this routine when parallel writing
!      is implemented.

  use constants, only: nprocs
  use HDF5, only:      H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, &
                       hid_t, hsize_t, &
                       h5dclose_f, h5dcreate_f, h5dwrite_f, &
                       h5sclose_f, h5screate_simple_f
  use io, only:        field_file_id, &
                       field_scalar_dataspace_id, &
                       field_name_cpu_time, &
                       field_name_total_wall_time

  implicit none

  real, intent(in)               :: elapsed_cpu_time
  real, intent(in)               :: elapsed_total_time

  integer(kind=4)                :: data_rank
  integer(hsize_t), dimension(2) :: data_dimensions

  ! HDF5 identifiers and status code.
  integer(hid_t)                 :: dataset_id
  integer(hid_t)                 :: vector_dataspace_id
  integer(hid_t)                 :: hdf5_status

  ! create a 1D dataspace to hold timings from each of the MPI ranks.
  data_rank          = 1
  data_dimensions(1) = nprocs
  data_dimensions(2) = 1
  call h5screate_simple_f( data_rank, data_dimensions, vector_dataspace_id, hdf5_status )

  ! create a vector of scalars for the elapsed CPU time for each of the MPI
  ! ranks, log_cpu_time.
  call h5dcreate_f( field_file_id, field_name_cpu_time, H5T_IEEE_F64LE, &
                    vector_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_cpu_time, data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

  ! close the 1D dataspace.
  call h5sclose_f( vector_dataspace_id, hdf5_status )

  ! create a scalar for the total wall time elapsed.
  call h5dcreate_f( field_file_id, field_name_total_wall_time, H5T_IEEE_F64LE, &
                    field_scalar_dataspace_id, dataset_id, hdf5_status )
  call h5dwrite_f( dataset_id, &
                   H5T_NATIVE_DOUBLE, elapsed_total_time, data_dimensions, &
                   hdf5_status )
  call h5dclose_f( dataset_id, hdf5_status )

end subroutine write_field_execution_stats

subroutine close_field_file
! Closes the field file.  Once closed, no further I/O may be performed on the
! field file without opening it again.  The timestep dataset identifiers
! (field_number_steps_dataset_id, field_wall_time_step_dataset_id,
! and field_step_numeric_error_dataset_id), the various dataspace identifiers
! (field_grid_dataspace_id, field_length1_vector_dataspace_id, and
! field_scalar_dataspace_id), and the file identifier (field_file_id) are
! closed.

  use HDF5, only:      hid_t, &
                       h5dclose_f, &
                       h5fclose_f, &
                       h5sclose_f
  use io, only:        field_file_id, &
                       field_grid_dataspace_id, field_length1_vector_dataspace_id, &
                       field_scalar_dataspace_id, field_number_steps_dataset_id, &
                       field_wall_time_step_dataset_id

  implicit none

  integer(hid_t) :: hdf5_status

  ! close the datasets used to maintain the timestep count, timing, and errors.
  call h5dclose_f( field_number_steps_dataset_id, hdf5_status )
  call h5dclose_f( field_wall_time_step_dataset_id, hdf5_status )

  ! close the grid, scalar, and vector scalar dataspaces.
  call h5sclose_f( field_grid_dataspace_id, hdf5_status )
  call h5sclose_f( field_length1_vector_dataspace_id, hdf5_status )
  call h5sclose_f( field_scalar_dataspace_id, hdf5_status )

  ! close the file so that all of the pending writes are flushed out.
  call h5fclose_f( field_file_id, hdf5_status )

end subroutine close_field_file
