!> Simple 1D MPI implementation of a blocking halo update
program sendrecv_1D
   use mpi
   implicit none

   integer :: ierr, rank, comsize
   integer :: left, right
   integer :: tag
   integer :: status(MPI_STATUS_SIZE)
   real(kind=4), allocatable :: array_halo_exchange(:)
   integer :: num_cells

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

   left = rank - 1
   if (left < 0) left = comsize - 1
   right = rank + 1
   if (right >= comsize) right = 0

   tag = 1

   num_cells = 4
   ! Create a halo region of 1
   allocate (array_halo_exchange(0:num_cells + 1), source=real(rank, kind=4))

   ! Send right value to the right neighbor and receive left halo from the left neighbour
   call MPI_Sendrecv(array_halo_exchange(num_cells), 1, MPI_REAL, right, tag, &
                     array_halo_exchange(0), 1, MPI_REAL, left, tag, &
                     MPI_COMM_WORLD, status, ierr)

   ! Send left value to the left neighbor and receive right halo from the right neighbour
   call MPI_Sendrecv(array_halo_exchange(1), 1, MPI_REAL, left, tag, &
                     array_halo_exchange(num_cells + 1), 1, MPI_REAL, right, tag, &
                     MPI_COMM_WORLD, status, ierr)

   ! Check the values; 1:num_cells will be equal to rank, 0 will be rank-1, num_cells+1 will be rank+1
   if (any(int(array_halo_exchange(1:num_cells)) /= rank)) error stop "Internal array values updated erroneously."
   if (int(array_halo_exchange(0)) /= left) error stop "Left halo cell not updated correctly"
   if (int(array_halo_exchange(num_cells + 1)) /= right) error stop "Right halo cell not updated correctly"

   call MPI_Finalize(ierr)
end program sendrecv_1D
