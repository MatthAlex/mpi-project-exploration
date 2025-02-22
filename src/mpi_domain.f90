!> Manages a 3D Cartesian grid for MPI-based simulations.
!> - Creates and configures MPI Cartesian communicator.
!> - Defines rank coordinates and nearest neighbors.
module mpi_domain
   use mpi_f08, only: MPI_Cart_coords, MPI_Cart_create, MPI_Cart_shift, MPI_Comm_rank, MPI_Comm_size, MPI_Dims_create
   use mpi_f08, only: MPI_COMM_WORLD, MPI_SUCCESS, MPI_Comm
   implicit none(type, external)
   private
   public :: create_mpi_domain

   integer :: ndims
      !! Number of dimensions for the Cartesian topology decomposition
   integer :: dims(3)
      !! Number of cores alongside the cardinal Cartesian directions
   logical :: is_periodic(3)
      !! Logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension.
   logical, parameter :: reorder = .true.
      !! Core ranking may be reordered (true) or not (false) (logical).

   integer, protected, public :: my_rank
      !! This rank
   integer, protected, public :: comsize
      !! Number of MPI processes
   integer, protected, public :: west, east, north, south, low, high
      !! Ranks of Nearest Neighbours
   integer, protected, public :: neighbor_ranks(6)
      !! Array of neighbor ranks for each of the 6 faces in 3D.

   type(MPI_Comm) :: original_comm = MPI_COMM_WORLD
      !! Original MPI communicator
   type(MPI_Comm), protected, public :: comm_cart
      !! Cartesian MPI communicator
   integer, protected, public :: my_coord(3)
      !! Cartesian MPI coordinates. Coordinates use 0-based indexing.

contains

   !> Initializes a 3D Cartesian communicator using MPI, setting
   !> dimensions, periodicity, ranks, and nearest neighbors.
   !>
   !> 1. Extract the size of MPI_COMM_WORLD.
   !> 2. Compute how to factor that into X, Y, Z dimensions with MPI_Dims_create.
   !> 3. Set periodic boundaries with set_periodic_boundaries.
   !> 4. Create the 3D Cartesian communicator (comm_cart).
   !> 5. Determine rank, coordinates, and neighbors in the new communicator.
   subroutine create_mpi_domain()
      use parameters, only: boundaries, core_decomposition, dim_decomposition
      integer :: ierr, original_comsize

      ! Probe the size of the original communicator
      call MPI_Comm_size(original_comm, original_comsize, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Comm_size"

      if (any(core_decomposition < 0)) error stop "Core decomposition can't have negative values."
      ndims = dim_decomposition
      dims = core_decomposition

      ! Creates a division of cores in a Cartesian ndims-dimensional X(, Y(, Z)) grid
      call MPI_Dims_create(nnodes=original_comsize, ndims=ndims, dims=dims, ierror=ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error creating dimensions"

      ! Set `is_periodic` to boundary-appropriate values
      call set_periodic_boundaries(boundaries)

      ! Makes a new communicator to which Cartesian topology information has been attached.
      call MPI_Cart_create(original_comm, ndims, dims, is_periodic, reorder, comm_cart, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error creating Cartesian communicator"

      ! Probe the size of the Cartesian communicator
      call MPI_Comm_size(comm_cart, comsize, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Comm_size"
      if (comsize /= original_comsize) print *, "Size mismatch between old and new communicators"

      ! Rank in the new communicator
      call MPI_Comm_rank(comm_cart, my_rank, ierr)
      if (my_rank == 0) print *, "Dimensions", dims
      if (ierr /= MPI_SUCCESS) error stop "Error deciding ranks in Cartesian"

      ! Get the 3D coordinates for this rank
      call MPI_Cart_coords(comm_cart, my_rank, ndims, my_coord, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error getting cartesian coordinates"

      ! Determine rank's nearest neighbors in all 6 directions
      call get_neighbors()
   end subroutine create_mpi_domain

   !> Finds the ranks of the 6 nearest neighbors by shifting ±1 in X, Y, and Z.
   !>
   !> Sets module variables: west, east, north, south, low, high
   subroutine get_neighbors()
      use enums, only: X_DIR, Y_DIR, Z_DIR
      integer :: ierr
      ! Call each shift along X, Y, Z axes, and get all 6 neighbours. Shift is ±1
      call MPI_Cart_shift(comm=comm_cart, direction=X_DIR, disp=1, rank_source=west, rank_dest=east, ierror=ierr)
      call MPI_Cart_shift(comm=comm_cart, direction=Y_DIR, disp=1, rank_source=south, rank_dest=north, ierror=ierr)
      call MPI_Cart_shift(comm=comm_cart, direction=Z_DIR, disp=1, rank_source=low, rank_dest=high, ierror=ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Cartesian shift"

      neighbor_ranks = [west, east, south, north, low, high]
   end subroutine get_neighbors

   !> Sets the periodic boundaries for the MPI Cartesian communicator, based on an array of boundary types.
   subroutine set_periodic_boundaries(bc_types)
      use enums, only: X_DIR, Y_DIR, Z_DIR, D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH, PERIODIC
      integer, intent(in) :: bc_types(6)
      ! Initialize to No Periodic Boundaries by default
      is_periodic = .false.
      ! `X_DIR`,..,`Z_DIR` is 0-based index where`is_periodic` is 1-based
      ! x-direction is periodic if both west and east are periodic
      is_periodic(X_DIR + 1) = (bc_types(D_WEST) == PERIODIC .and. bc_types(D_EAST) == PERIODIC)

      ! y-direction is periodic if both south and north are periodic
      is_periodic(Y_DIR + 1) = (bc_types(D_SOUTH) == PERIODIC .and. bc_types(D_NORTH) == PERIODIC)

      ! z-direction is periodic if both low and high are periodic
      is_periodic(Z_DIR + 1) = (bc_types(D_LOW) == PERIODIC .and. bc_types(D_HIGH) == PERIODIC)

   end subroutine set_periodic_boundaries

end module mpi_domain
