!> Extension of the simple 1D MPI implementation.
!> - It leaves the topology generation to MPI
!> - Can handle periodic boundaries with a boolean flag
!> - Heuristics may reorder the original comm (CPU topology) into a more efficient one
!> - Figures out who the neighbours are of each rank automatically
program sendrecv_3D
   use mpi
   use grid_module, only: initialize_MPI_grid, &
      left => west, right => east, down => south, up => north, back => low, front => high, &
      comm_cart, rank => my_rank
   implicit none

   integer :: ierr
   integer :: status(MPI_STATUS_SIZE)
   real(kind=4), allocatable :: array(:, :, :)
   integer :: num_cells
   integer :: tag
   integer :: i

   call initialize_MPI_grid()

   ! Create a local array with halo regions
   num_cells = 1
   allocate (array(0:num_cells + 1, 0:num_cells + 1, 0:num_cells + 1), source=real(rank, kind=4))

   tag = 1

   do i = 1, 1
      ! Perform halo exchanges in all 3 directions
      ! Exchange X direction: First to the right, followed by to the left
      call MPI_Sendrecv(array(num_cells, :, :), (num_cells + 2) * (num_cells + 2), MPI_REAL, right, tag, &
                        array(0, :, :), (num_cells + 2) * (num_cells + 2), MPI_REAL, left, tag, &
                        comm_cart, status, ierr)

      call MPI_Sendrecv(array(1, :, :), (num_cells + 2) * (num_cells + 2), MPI_REAL, left, tag, &
                        array(num_cells + 1, :, :), (num_cells + 2) * (num_cells + 2), MPI_REAL, right, tag, &
                        comm_cart, status, ierr)

      ! Exchange Y direction (up/down)
      call MPI_Sendrecv(array(:, num_cells, :), (num_cells + 2) * (num_cells + 2), MPI_REAL, up, tag, &
                        array(:, 0, :), (num_cells + 2) * (num_cells + 2), MPI_REAL, down, tag, &
                        comm_cart, status, ierr)

      call MPI_Sendrecv(array(:, 1, :), (num_cells + 2) * (num_cells + 2), MPI_REAL, down, tag, &
                        array(:, num_cells + 1, :), (num_cells + 2) * (num_cells + 2), MPI_REAL, up, tag, &
                        comm_cart, status, ierr)

      ! Exchange Z direction (front/back)
      call MPI_Sendrecv(array(:, :, num_cells), (num_cells + 2) * (num_cells + 2), MPI_REAL, front, tag, &
                        array(:, :, 0), (num_cells + 2) * (num_cells + 2), MPI_REAL, back, tag, &
                        comm_cart, status, ierr)

      call MPI_Sendrecv(array(:, :, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, back, tag, &
                        array(:, :, num_cells + 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, front, tag, &
                        comm_cart, status, ierr)
   end do

   call MPI_Finalize(ierr)

end program sendrecv_3D
