submodule(mpi_domain_types) mpi_domain_extended_neighbors
   implicit none(type, external)
   logical :: DEBUG = .false.

contains

   module subroutine determine_extended_neighbors(self)
      !! Finds the ranks of all 26 nearest neighbors in a 3x3x3 stencil, centered on self%rank
      class(mpi_domain_t), intent(in out) :: self
      integer :: temp_coords(3)
      integer :: i, j, k, ierr, d
      logical :: oob

      do k = -1, 1
         do j = -1, 1
            do i = -1, 1

               if (i == 0 .and. j == 0 .and. k == 0) then
                  ! centre point = self
                  self%neighbors_extended(i, j, k) = self%rank
                  cycle
               end if

               ! compute candidate coords
               temp_coords = self%coords + [i, j, k]

               ! check out-of-bounds on non-periodic dims
               oob = .false.
               do d = 1, 3
                  if (.not. self%periodic(d) &
                      .and. (temp_coords(d) < 0 &
                             .or. temp_coords(d) >= self%dims(d))) then
                     oob = .true.
                     exit
                  end if
               end do

               if (oob) then
                  self%neighbors_extended(i, j, k) = MPI_PROC_NULL
               else
                  call MPI_Cart_rank(self%comm, temp_coords, self%neighbors_extended(i, j, k), ierr)
                  if (ierr /= MPI_SUCCESS) call self%abort("MPI_Cart_rank failed in finding extended neighbors.. Exiting..")
               end if

            end do
         end do
      end do

      do k = -1, 1
         do j = -1, 1
            do i = -1, 1
               if (i == 0 .and. j == 0 .and. k == 0) then
                  if (DEBUG) write (*, '(A,3(I2,1X),A,I0)') "Center (", i, j, k, "): Self = ", self%rank
               else
                  if (DEBUG) write (*, '(A,3(I2,1X),A,I0)') "Pos (", i, j, k, "): Neighbor = ", &
                     self%neighbors_extended(i, j, k)
               end if
            end do
         end do
      end do
   end subroutine determine_extended_neighbors

   !> After filling neighbours_extended, this asserts each entry is exactly what MPI says.
   module subroutine test_extended_neighbors(self)
      class(mpi_domain_t), intent(in) :: self
      integer :: i, j, k, d, ierr
      integer :: tmpc(3), expected, actual
      logical :: oob

      do k = -1, 1
         do j = -1, 1
            do i = -1, 1

               ! actual value we stored
               actual = self%neighbors_extended(i, j, k)

               ! compute what MPI_Cart_rank *should* give (or NULL if OOB on non-periodic)
               tmpc = self%coords + [i, j, k]
               oob = .false.
               do d = 1, 3
                  if (.not. self%periodic(d) .and. &
                      (tmpc(d) < 0 .or. tmpc(d) >= self%dims(d))) then
                     oob = .true.
                     exit
                  end if
                  if (self%periodic(d)) then
                     ! wrap around
                     tmpc(d) = mod(tmpc(d), self%dims(d))
                  end if
               end do

               if (oob) then
                  expected = MPI_PROC_NULL
               else
                  call MPI_Cart_rank(self%comm, tmpc, expected, ierr)
                  if (ierr /= MPI_SUCCESS) call MPI_Abort(self%comm, ierr)
               end if

               if (actual /= expected) then
                  write (*, *) "[Rank", self%rank, &
                     "] mismatch at (", i, ",", j, ",", k, &
                     "): got", actual, "expected", expected
                  call MPI_Abort(self%comm, 1)
               end if

            end do
         end do
      end do
      print *, "Successfully tested extended neighbors"
   end subroutine test_extended_neighbors
end submodule mpi_domain_extended_neighbors
