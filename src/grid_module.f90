module grid_module
   use mpi, only: MPI_Cart_coords, MPI_Cart_create, MPI_Cart_shift, MPI_Comm_rank, MPI_Comm_size, MPI_Dims_create
   use mpi, only: MPI_COMM_WORLD, MPI_SUCCESS
   use enums, only: X_DIR, Y_DIR, Z_DIR, D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH, PERIODIC
   implicit none (type, external)
   private
   public :: initialize_MPI_grid

   !> Number of cores alongside the cardinal Cartesian directions
   integer :: dims(3) = [0, 0, 0]
   !> Number of dimensions for the Cartesian topology decomposition
   integer, parameter :: ndims = 3
   !> Periodic Boundaries
   logical :: is_periodic(3)
   logical, parameter :: reorder = .true.

   !> This rank
   integer, public :: my_rank
   !> Number of MPI processes
   integer, public :: comsize
   !> Ranks of Nearest Neighbours
   integer, public :: west, east, north, south, low, high
   !> Cartesian MPI communicator
   integer, public :: comm_cart
   !> Cartesian MPI coordinates. Coordinates are 0-based.
   integer, public :: my_coordinates(ndims)

contains

   !> Initializes a 3D Cartesian communicator using MPI, setting
   !> dimensions, periodicity, ranks, and nearest neighbors.
   !>
   !> 1. Extract the size of MPI_COMM_WORLD.
   !> 2. Compute how to factor that into X, Y, Z dimensions with MPI_Dims_create.
   !> 3. Set periodic boundaries with set_periodic_boundaries.
   !> 4. Create the 3D Cartesian communicator (comm_cart).
   !> 5. Determine rank, coordinates, and neighbors in the new communicator.
   subroutine initialize_MPI_grid()
      use lib_parameters, only: boundaries
      integer :: original_comm, ierr, original_comsize
      original_comm = MPI_COMM_WORLD

      ! Probe the size of the original communicator
      call MPI_Comm_size(original_comm, original_comsize, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Comm_size"

      ! Compute dims(3) that factor original_comsize into X, Y, Z processes
      call MPI_Dims_create(original_comsize, ndims, dims, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error creating dimensions"

      ! Set `is_periodic` to boundary-appropriate values
      call set_periodic_boundaries(boundaries)

      ! Create the new Cartesian communicator
      call MPI_Cart_create(original_comm, ndims, dims, is_periodic, reorder, comm_cart, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error creating Cartesian communicator"

      ! Probe the size of the Cartesian communicator
      call MPI_Comm_size(comm_cart, comsize, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Comm_size"
      if (comsize /= original_comsize) error stop "Size mismatch between old and new communicators"

      ! Rank in the new communicator
      call MPI_Comm_rank(comm_cart, my_rank, ierr)
      if (my_rank == 0) print *, "Dimensions", dims
      if (ierr /= MPI_SUCCESS) error stop "Error deciding ranks in Cartesian"

      ! Get the 3D coordinates for this rank
      call MPI_Cart_coords(comm_cart, my_rank, ndims, my_coordinates, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error getting cartesian coordinates"

      ! Determine rank's nearest neighbors in all 6 directions
      call get_neighbors()
   end subroutine initialize_MPI_grid

   !> Finds the ranks of the 6 nearest neighbors by shifting ±1 in X, Y, and Z.
   !>
   !> Sets module variables: west, east, north, south, low, high
   subroutine get_neighbors()
      integer :: ierr
      ! Call each shift along X, Y, Z axes, and get all 6 neighbours. Shift is ±1
      call MPI_Cart_shift(comm=comm_cart, direction=X_DIR, disp=1, rank_source=west, rank_dest=east, ierror=ierr)
      call MPI_Cart_shift(comm=comm_cart, direction=Y_DIR, disp=1, rank_source=south, rank_dest=north, ierror=ierr)
      call MPI_Cart_shift(comm=comm_cart, direction=Z_DIR, disp=1, rank_source=low, rank_dest=high, ierror=ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Cartesian shift"
   end subroutine get_neighbors

   !> Sets the periodic boundaries based on an array of boundary types.
   subroutine set_periodic_boundaries(bc_types)
      integer, intent(in) :: bc_types(6)
      ! Initialize to No Periodic Boundaries by default
      is_periodic = .false.
      ! x-direction is periodic if both west and east are periodic
      is_periodic(X_DIR) = (bc_types(D_WEST) == PERIODIC .and. bc_types(D_EAST) == PERIODIC)

      ! y-direction is periodic if both south and north are periodic
      is_periodic(Y_DIR) = (bc_types(D_SOUTH) == PERIODIC .and. bc_types(D_NORTH) == PERIODIC)

      ! z-direction is periodic if both low and high are periodic
      is_periodic(Z_DIR) = (bc_types(D_LOW) == PERIODIC .and. bc_types(D_HIGH) == PERIODIC)

   end subroutine set_periodic_boundaries

end module grid_module
