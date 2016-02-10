subroutine display_inputs
! This subroutine outputs some information about the simulation being run.

  use constants
  !$ use omp_lib, only:        omp_get_max_threads

  implicit none

  integer           :: number_openmp_threads
  character(len=16) :: caststr

  ! Display GLL grid parameters.
  call notify( '' )
  call notify( 'Grid information:' )
  call notify( '' )

  write( caststr,'(I12.0)' ) n
  call notify( ' n       = ' // caststr )

  write( caststr,'(I12.0)' ) mx
  call notify( ' nsubx   = ' // caststr )

  write( caststr,'(I12.0)' ) my
  call notify( ' nsuby   = ' // caststr )

  write( caststr,'(I12.0)' ) mz
  call notify( ' nsubz   = ' // caststr )

  write( caststr,'(I12.0)' ) r
  call notify( ' ntotal  = ' // caststr )

  call notify( ' ' )
  call notify( 'Per processor: ' )
  call notify( ' ' )

  write( caststr,'(I12.0)' ) r_loc
  call notify( ' ntotal  = ' // caststr )

  write( caststr,'(I12.0)' ) mx_loc
  call notify( ' nsubx   = ' // caststr )

  write( caststr,'(I12.0)' ) my_loc
  call notify( ' nsuby   = ' // caststr )

  write( caststr,'(I12.0)' ) mz_loc
  call notify( ' nsubz   = ' // caststr )

  write( caststr,'(I12.0)' ) nx_loc
  call notify( ' nx      = ' // caststr )

  write( caststr,'(I12.0)' ) ny_loc
  call notify( ' ny      = ' // caststr )

  write( caststr,'(I12.0)' ) nz_loc
  call notify( ' nz      = ' // caststr )


  call notify( '' )
  call notify( 'Parallelization parameters:' )
  call notify( '' )

  write( caststr,'(I12.0)' ) nprocs
  call notify( ' number of MPI processes  = ' // trim( adjustl( caststr ) ) )

  number_openmp_threads    = 1
  !$ number_openmp_threads = omp_get_max_threads()

  write( caststr, '(I12.0)' ) number_openmp_threads
  call notify( ' number of openMP threads = ' // trim( adjustl( caststr ) ) )
  call notify( ' ' )

end subroutine display_inputs
