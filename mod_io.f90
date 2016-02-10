module io
! Contains global state for writing the field file as well as constant field
! file field names.

  use HDF5, only:                 hid_t

  implicit none

  ! Identifiers that are used throughout the lifetime of the solver when
  ! writing to the field file.
  integer(hid_t)               :: field_file_id

  integer(hid_t)               :: field_grid_dataspace_id
  integer(hid_t)               :: field_length1_vector_dataspace_id
  integer(hid_t)               :: field_scalar_dataspace_id

  integer(hid_t)               :: field_number_steps_dataset_id
  integer(hid_t)               :: field_wall_time_step_dataset_id
!  integer(hid_t)               :: field_gmres_diffusion_iterations_dataset_id
!  integer(hid_t)               :: field_gmres_poisson_iterations_dataset_id
!  integer(hid_t)               :: field_gmres_viscous_iterations_dataset_id
!  integer(hid_t)               :: field_l2_poisson_error_dataset_id
!  integer(hid_t)               :: field_l2_poisson_kernel_error_dataset_id
!  integer(hid_t)               :: field_l2_poisson_schur_error_dataset_id
!  integer(hid_t)               :: field_l2_schur_kernel_error_dataset_id
!  integer(hid_t)               :: field_linf_diffusion_error_dataset_id
!  integer(hid_t)               :: field_linf_divergence_error_dataset_id
!  integer(hid_t)               :: field_linf_viscous_x_error_dataset_id
!  integer(hid_t)               :: field_linf_viscous_z_error_dataset_id

  ! Paths to the grid's datasets.

  ! Configuration parameters.
  character(len=14), parameter :: field_name_config                        = "/configuration"
!  character(len=17), parameter :: field_name_uL                            = "/configuration/uL"
!  character(len=17), parameter :: field_name_uC                            = "/configuration/uC"
!  character(len=23), parameter :: field_name_facrobin                      = "/configuration/facrobin"
!  character(len=27), parameter :: field_name_facrobin_ppe                  = "/configuration/facrobin_ppe"
  character(len=24), parameter :: field_name_mu                            = "/configuration/viscosity"
!  character(len=27), parameter :: field_name_filter_order                  = "/configuration/filter_order"
  character(len=22), parameter :: field_name_dt                            = "/configuration/delta_t"
  character(len=29), parameter :: field_name_tend                          = "/configuration/time_threshold"
  character(len=30), parameter :: field_name_timesteps_per_write           = "/configuration/steps_per_write"

!  ! GMRES-specific configuration parameters.
!  character(len=20), parameter :: field_name_config_gmres                  = "/configuration/gmres"
!  character(len=38), parameter :: field_name_poisson_tolerance             = "/configuration/gmres/poisson_tolerance"
!  character(len=43), parameter :: field_name_poisson_max_iters             = "/configuration/gmres/poisson_max_iterations"
!  character(len=36), parameter :: field_name_poisson_restart               = "/configuration/gmres/poisson_restart"
!  character(len=38), parameter :: field_name_viscous_tolerance             = "/configuration/gmres/viscous_tolerance"
!  character(len=43), parameter :: field_name_viscous_max_iters             = "/configuration/gmres/viscous_max_iterations"

  ! Grid parameters.
  character(len=5),  parameter :: field_name_grid                          = "/grid"
  character(len=7),  parameter :: field_name_n                             = "/grid/n"
  character(len=8),  parameter :: field_name_mx                            = "/grid/mx"
  character(len=8),  parameter :: field_name_my                            = "/grid/my"
  character(len=8),  parameter :: field_name_mz                            = "/grid/mz"
  character(len=7),  parameter :: field_name_x                             = "/grid/x"
  character(len=7),  parameter :: field_name_y                             = "/grid/y"
  character(len=7),  parameter :: field_name_z                             = "/grid/z"

  ! Field variables.
  !
  ! NOTE: s, ux, and uz are relative names since the group that contains
  !       them is dynamically constructed at run-time from
  !       field_name_step_base.
  character(len=6),  parameter :: field_name_field                         = "/field"
  character(len=19), parameter :: field_name_number_steps                  = "/field/number_steps"
  character(len=20), parameter :: field_name_step_base                     = "/field/step" ! additional length since we append the timestep count to it
  character(len=1),  parameter :: field_name_s                             = "s"
  character(len=2),  parameter :: field_name_ux                            = "ux"
  character(len=2),  parameter :: field_name_uy                            = "uy"
  character(len=2),  parameter :: field_name_uz                            = "uz"

  ! Execution statistics.
  character(len=10), parameter :: field_name_execution                     = "/execution"
  character(len=23), parameter :: field_name_number_ranks                  = "/execution/number_ranks"
  character(len=25), parameter :: field_name_number_threads                = "/execution/number_threads"
!  character(len=16), parameter :: field_name_error                         = "/execution/error"
!  character(len=43), parameter :: field_name_gmres_diffusion_iters         = "/execution/error/gmres_diffusion_iterations"
!  character(len=41), parameter :: field_name_gmres_poisson_iters           = "/execution/error/gmres_poisson_iterations"
!  character(len=41), parameter :: field_name_gmres_viscous_iters           = "/execution/error/gmres_viscous_iterations"
!  character(len=27), parameter :: field_name_l2_poisson_error              = "/execution/error/l2_poisson"
!  character(len=34), parameter :: field_name_l2_poisson_kernel_error       = "/execution/error/l2_poisson_kernel"
!  character(len=33), parameter :: field_name_l2_poisson_schur_error        = "/execution/error/l2_poisson_schur"
!  character(len=40), parameter :: field_name_l2_poisson_schur_kernel_error = "/execution/error/l2_poisson_schur_kernel"
!  character(len=31), parameter :: field_name_linf_diffusion_error          = "/execution/error/linf_diffusion"
!  character(len=32), parameter :: field_name_linf_divergence_error         = "/execution/error/linf_divergence"
!  character(len=31), parameter :: field_name_linf_viscous_x_error          = "/execution/error/linf_viscous_x"
!  character(len=31), parameter :: field_name_linf_viscous_z_error          = "/execution/error/linf_viscous_z"
!  character(len=41), parameter :: field_name_total_numeric_error           = "/execution/error/total_norm_numeric_error"
  character(len=19), parameter :: field_name_cpu_time                      = "/execution/cpu_time"
  character(len=20), parameter :: field_name_wall_time                     = "/execution/wall_time"
  character(len=26), parameter :: field_name_start_time                    = "/execution/wall_time/start"
  character(len=26), parameter :: field_name_total_wall_time               = "/execution/wall_time/total"
!  character(len=26), parameter :: field_name_setup_wall_time               = "/execution/wall_time/setup"
!  character(len=31), parameter :: field_name_null_basis_wall_time          = "/execution/wall_time/null_basis"
!  character(len=31), parameter :: field_name_null_error_wall_time          = "/execution/wall_time/null_error"
  character(len=27), parameter :: field_name_step_wall_time                = "/execution/wall_time/steps"

end module io
