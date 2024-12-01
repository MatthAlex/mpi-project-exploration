!> Extension of the simple 1D MPI implementation.
!> - It leaves the topology generation to MPI
!> - Can handle periodic boundaries with a boolean flag
!> - Heuristics may reorder the original comm (CPU topology) into a more efficient one
!> - Figures out who are the neighbours of each rank automatically
program sendrecv_3D
   use mpi
   implicit none

   integer :: ierr, rank, comsize, comm
   integer :: ndims, dim_size(3), coords(3)
   logical :: periodic(1:3), reorder
   integer :: left, right, up, down, front, back
   integer :: comm_cart
   integer :: status(MPI_STATUS_SIZE)
   real(kind=4), allocatable :: array(:, :, :)
   integer :: num_cells
   integer :: tag
   integer :: i

   call MPI_Init(ierr)
   comm = MPI_COMM_WORLD
   call MPI_Comm_size(comm, comsize, ierr)
   call MPI_Comm_rank(comm, rank, ierr)

   ! Dimensions of the Cartesian grid (3D)
   ndims = 3
   ! Let MPI_Dims_create decide the dimensions automatically
   dim_size = [0, 0, 0]

   call MPI_Dims_create(comsize, ndims, dim_size, ierr)

   ! Enable/disable periodic boundary conditions
   periodic(1) = .false. ! X
   periodic(2) = .false. ! Y
   periodic(3) = .false. ! Z
   ! If MPI can find a better topology, it can reorder the original Comm
   reorder = .true.

   ! Create the Cartesian communicator
   call MPI_Cart_create(comm, ndims, dim_size, periodic, reorder, comm_cart, ierr)

   ! Get the coordinates of this process in the grid
   call MPI_Cart_coords(comm_cart, rank, ndims, coords, ierr)

   ! Create a local array with halo regions
   num_cells = 1
   allocate (array(0:num_cells + 1, 0:num_cells + 1, 0:num_cells + 1), source=real(rank, kind=4))

   ! Compute the neighbors in each direction
   call MPI_Cart_shift(comm_cart, 0, 1, left, right, ierr) ! X direction neighbors
   call MPI_Cart_shift(comm_cart, 1, 1, down, up, ierr) ! Y direction neighbors
   call MPI_Cart_shift(comm_cart, 2, 1, back, front, ierr) ! Z direction neighbors

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
