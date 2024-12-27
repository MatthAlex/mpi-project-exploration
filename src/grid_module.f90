module grid_module
   use mpi
   implicit none
   private :: get_neighbors
   public :: initialize_MPI_grid

   !> Cartesian directions
   integer, parameter :: X_DIR = 0, Y_DIR = 1, Z_DIR = 2
   !> Number of cores alongside the cardinal Cartesian directions
   integer :: dims(3) = [0, 0, 0]
   !> Number of dimensions for the Cartesian topology decomposition
   integer, parameter :: ndims = 3
   !> Periodic Boundaries
   logical, parameter :: periods(3) = [.true., .true., .true.]
   logical, parameter :: reorder = .true.

   !> This rank
   integer, public :: my_rank
   !> Number of MPI processes
   integer, public :: comsize
   !> Ranks of Nearest Neighbours
   integer, public :: west, east, north, south, low, high
   !> Cartesian MPI communicator
   integer, public :: comm_cart
   !> Cartesian MPI coordinates
   integer, public :: my_coordinates(ndims)

contains

   !> Initializes the MPI grid with a 3D Cartesian topology.
   !> Entry point of MPI functionality.
   subroutine initialize_MPI_grid()
      integer :: original_comm, ierr, original_comsize
      original_comm = MPI_COMM_WORLD

      !> Probe the original communicator size
      call MPI_Comm_size(original_comm, original_comsize, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Comm_size"

      call MPI_Dims_create(original_comsize, ndims, dims, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error creating dimensions"

      call MPI_Cart_create(original_comm, ndims, dims, periods, reorder, comm_cart, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error creating Cartesian communicator"

      !> Probe the Cartesian communicator size
      call MPI_Comm_size(comm_cart, comsize, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Comm_size"
      if (comsize /= original_comsize) error stop "Size mismatch between old and new communicators"

      call MPI_Comm_rank(comm_cart, my_rank, ierr)
      if (my_rank == 0) print *, "Dimensions", dims
      if (ierr /= MPI_SUCCESS) error stop "Error deciding ranks in Cartesian"

      ! Populate core Cartesian coordinates
      call MPI_Cart_coords(comm_cart, my_rank, ndims, my_coordinates, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error getting cartesian coordinates"

      ! Populate Nearest Neighbours
      call get_neighbors()
   end subroutine initialize_MPI_grid

   !> Performs a displacement 1 shift across all directions. These are the nearest 6 neighbours
   subroutine get_neighbors()
      integer :: ierr
      ! Call each shift along X, Y, Z axes, and get all 6 neighbours. Shift is + and -
      call MPI_Cart_shift(comm=comm_cart, direction=X_DIR, disp=1, rank_source=west, rank_dest=east, ierror=ierr)
      call MPI_Cart_shift(comm=comm_cart, direction=Y_DIR, disp=1, rank_source=south, rank_dest=north, ierror=ierr)
      call MPI_Cart_shift(comm=comm_cart, direction=Z_DIR, disp=1, rank_source=low, rank_dest=high, ierror=ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Cartesian shift"
   end subroutine get_neighbors

end module grid_module
