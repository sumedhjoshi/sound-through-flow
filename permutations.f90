subroutine permute_z( n, mx, mz, u, trans )
! This subroutine computes the permutation to/from the standard coordinates to z coordinates.
!
! Inputs.
!
! n  [1x1]      : number of GLL points in one direction in one subdomain.
! mx [1x1]      : number of subdomains in the x direction.
! mz [1x1]      : number of subdomains in the z direction.
! trans [1x1]       : booleans that determine whether we apply Px or its transpose (and Pz).
! u [n^2*m^2x1]     : the vector that we're permuting.
!
! Outputs.
!
! permute_z [n^2*m^2x1] : the permuted vector.
!
! 22 July 2013
! Sumedh Joshi

  implicit none

  integer, intent(in)                         :: n, mx, mz
  real, dimension(1:n*mx*mz*n), intent(inout) :: u
  logical, intent(in)                         :: trans
  integer                                     :: zsub, xsub, ii, jj, counter
  integer                                     :: index_zcoords, index_standard
  real, dimension(1:n*mx*mz*n)                :: swap

  ! Compute the z-permutation.
  if (trans .eqv. .false.) then
     counter = 1

     do xsub = 1, mx
        do jj = 1, n
           do zsub = 1, mz
              do ii = 1, n

                 ! Compute the indices.
                 index_standard      = (xsub - 1) * n * n + &
                                       (zsub - 1) * n * n * mx + &
                                       (ii - 1) * n + jj
                 index_zcoords       = counter
                 counter             = counter + 1
                 swap(index_zcoords) = u(index_standard)

              enddo
           enddo
        enddo
     enddo
  endif

  if (trans .eqv. .true.) then
     counter = 1

     do xsub = 1, mx
        do jj = 1, n
           do zsub = 1, mz
              do ii = 1, n

                 ! Compute the indices.
                 index_standard       = (xsub - 1) * n * n + &
                                        (zsub - 1) * n * n * mx + &
                                        (ii - 1) * n + jj
                 index_zcoords        = counter
                 counter              = counter + 1
                 swap(index_standard) = u(index_zcoords)

              enddo
           enddo
        enddo
     enddo
  endif

  ! Update u from the swapped value.
  u = swap

end subroutine permute_z

subroutine permute_x( n, mx, mz, u, trans )
! This subroutine computes the permutation to/from the standard coordinates to x coordinates.
!
! Inputs.
!
! n  [1x1]      : number of GLL points in one direction in one subdomain.
! mx [1x1]      : number of subdomains in the x direction.
! mz [1x1]      : number of subdomains in the z direction.
! trans [1x1]       : booleans that determine whether we apply Px or its transpose (and Pz).
! u [n^2*m^2x1]     : the vector that we're permuting.
!
! Outputs.
!
! permute_x [n^2*m^2x1] : the permuted vector.
!
! 22 July 2013
! Sumedh Joshi

  implicit none

  integer, intent(in)                        :: n,mx,mz
  real, dimension(1:n*mz*mx*n),intent(inout) :: u
  logical, intent(in)                        :: trans
  integer                                    :: zsub, xsub, ii, jj, counter
  integer                                    :: index_xcoords, index_standard
  real, dimension(1:n*mx*mz*n)               :: swap

  ! Compute the x-permutation.
  counter = 1
  do zsub = 1,mz
     do jj = 1,n
        do xsub = 1,mx
           do ii = 1,n

              ! Compute the indices.
              index_standard = (zsub - 1) * n * n * mx + &
                               (xsub - 1) * n * n + &
                               (jj   - 1) * n + ii
              index_xcoords  = counter
              counter        = counter + 1

              ! Compute (standard --> x) or (x --> standard).
              if(trans .eqv. .false.) then
                 swap(index_xcoords) = u(index_standard)
              else
                 swap(index_standard) = u(index_xcoords)
              endif

           enddo
        enddo
     enddo
  enddo

  ! Update x from the swapped value.
  u = swap

end subroutine permute_x

subroutine perfect_shuffle_rows( p, q, r, x )
! x is a matrix with pq rows and r columns.  This perfect shuffles the rows of x.
! Splits x into p stacks of q elements and deals one from each stack q times.  number_of_rows(x) = p*q.

  implicit none

  integer, intent(in)                     :: p, q, r
  real, dimension(1:p*q,r), intent(inout) :: x
  integer                                 :: ii, jj, ndx, counter
  real, dimension(1:p*q,r)                :: xnew

  counter = 1
  do jj = 1, q
     do ii = 1, p
        ndx             = (ii - 1) * q + jj
        xnew(counter,:) = x(ndx,:)
        counter         = counter + 1
     enddo
  enddo

  ! Returned the shuffled vector.
  x = xnew

end subroutine perfect_shuffle_rows

subroutine permute_by( n, p, x )
! Permutes x with the permutation vector p so that
! x <-- x(p).

  implicit none

  integer, intent(in)                   :: n
  real, dimension(1:n), intent(inout)   :: x
  integer, dimension(1:n), intent(in)   :: p
  integer                               :: ii
  real, dimension(1:n)                  :: xnew

  do ii = 1,n
     xnew(ii) = x(p(ii))
  enddo

  ! Returned the shuffled vector.
  x = xnew

end subroutine permute_by

subroutine perfect_shuffle( p, q, x )
! Splits x into p stacks of q elements and deals one from each stack q times.  length(x) = p*q.

  implicit none

  integer, intent(in)                   :: p, q
  real, dimension(1:p*q), intent(inout) :: x
  integer                               :: ii, jj, ndx, counter
  real, dimension(1:p*q)                :: xnew

  !xnew    = 0.
  counter = 1
  do jj = 1, q
     do ii = 1, p
        ndx           = (ii - 1) * q + jj
        xnew(counter) = x(ndx)
        counter       = counter + 1
     enddo
  enddo

  ! Returned the shuffled vector.
  x = xnew

end subroutine perfect_shuffle

