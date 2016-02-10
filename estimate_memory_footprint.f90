subroutine estimate_memory_footprint
! This subroutine outputs some information about the simulation being run.

  use constants
  !$ use omp_lib, only:        omp_get_max_threads

  implicit none

  character(len=16) :: caststr
  integer           :: count_dhydro, count_metric, count_grid, count_fields, &
                       count_hydro, count_work, count_wave, &
                       reals_per_function, bytes_per_real
  real              :: footprint_in_bytes, footprint_in_Gib, footprint_in_Gib_per_rank

  ! Set some constants.
  bytes_per_real = 8

  ! Set the number of reals to store a grid function.
  reals_per_function = r

  ! Set some counters for the number of grid functions.
  count_grid   = 3      ! 3 grid directions to store.
  count_metric = 9      ! 9 metric functions.
  count_fields = 4 * 2  ! 4 fields (ux,uy,uz,s) and 2 lagging timesteps.
  count_hydro  = 5      ! 5 hydrodynamic variables (ux,uy,uz,rho,beta).
  count_dhydro = 18     ! 18 arrays to store derivatives of hydrodynamic variables.
  count_work   = 3      ! 3 work arrays.
  count_wave   = 12     ! 12 arrays to store the wave operator.

  ! Count up the memory required in bytes.
  footprint_in_bytes        = real (bytes_per_real ) * real ( reals_per_function ) &
                                                                  * ( count_metric + &
                                                                      count_fields + &
                                                                      count_hydro + &
                                                                      count_grid + &
                                                                      count_wave + &
                                                                      count_work + &
                                                                      count_dhydro )
  footprint_in_Gib          = footprint_in_bytes / 1024.0 / 1024.0 / 1024.0
  footprint_in_Gib_per_rank = footprint_in_Gib / nprocs

  ! Display GLL grid parameters.
  call notify( 'Estimated memory footprint:' )
  call notify( '' )

  write( caststr,'(I16.0)' ) count_dhydro + count_grid + count_metric + &
                             count_fields + count_hydro + count_work + count_wave
  call notify( ' Grid-scale arrays to store     : ' // caststr )

  write( caststr,'(F16.8)' ) real( reals_per_function * bytes_per_real ) / 1024.0 / 1024.0 / 1024.0
  call notify( ' Memory per array (GB)          : ' // caststr )

  write( caststr,'(F16.8)' ) footprint_in_Gib
  call notify( ' Memory footprint (GB)          : ' // caststr )

  write( caststr,'(F16.8)' ) footprint_in_Gib_per_rank
  call notify( ' Memory footprint per rank (GB) : ' // caststr )

  call notify( ' ' )

end subroutine estimate_memory_footprint
