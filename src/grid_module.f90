module grid_module
   use mpi, only: MPI_Cart_coords, MPI_Cart_create, MPI_Cart_shift, MPI_Comm_rank, MPI_Comm_size, MPI_Dims_create
   use mpi, only: MPI_COMM_WORLD, MPI_SUCCESS
   implicit none (type, external)
   private
   public :: initialize_MPI_grid

   !> Cartesian directions
   integer, parameter :: X_DIR = 0, Y_DIR = 1, Z_DIR = 2
   !> Number of cores alongside the cardinal Cartesian directions
   integer :: dims(3) = [0, 0, 0]
   !> Number of dimensions for the Cartesian topology decomposition
   integer, parameter :: ndims = 3
   !> Periodic Boundaries
   logical :: periods(3)
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

      ! Set `periods` to boundary-appropriate values
      call set_periodic_boundaries(boundaries)

      ! Create the new Cartesian communicator
      call MPI_Cart_create(original_comm, ndims, dims, periods, reorder, comm_cart, ierr)
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

      ! Determine our nearest neighbors in all 6 directions
      call get_neighbors()
   end subroutine initialize_MPI_grid

   !> Finds the ranks of the 6 nearest neighbors by shifting ±1 in X, Y, and Z.
   !>
   !> Sets module variables: west, east, north, south, low, high
   subroutine get_neighbors()
      integer :: ierr
      ! Call each shift along X, Y, Z axes, and get all 6 neighbours. Shift is + and -
      call MPI_Cart_shift(comm=comm_cart, direction=X_DIR, disp=1, rank_source=west, rank_dest=east, ierror=ierr)
      call MPI_Cart_shift(comm=comm_cart, direction=Y_DIR, disp=1, rank_source=south, rank_dest=north, ierror=ierr)
      call MPI_Cart_shift(comm=comm_cart, direction=Z_DIR, disp=1, rank_source=low, rank_dest=high, ierror=ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Cartesian shift"
   end subroutine get_neighbors

   !> Sets the periodic boundaries based on an array of boundary codes.
   !>
   !> bc_types   6-element array of boundary codes for each face:
   !> bc_types(1)=west, bc_types(2)=east, bc_types(3)=south,
   !> bc_types(4)=north, bc_types(5)=low, bc_types(6)=high.
   !> A code of 0 indicates periodic in that face/dimension.
   subroutine set_periodic_boundaries(bc_types)
      integer, intent(in) :: bc_types(6)
      ! Initialize to No Periodic Boundaries by default
      periods = .false.
      ! x-direction is periodic if west and east are 0
      periods(1) = (bc_types(1) == 0 .and. bc_types(2) == 0)

      ! y-direction is periodic if south and north are 0
      periods(2) = (bc_types(3) == 0 .and. bc_types(4) == 0)

      ! z-direction is periodic if low and high are 0
      periods(3) = (bc_types(5) == 0 .and. bc_types(6) == 0)

   end subroutine set_periodic_boundaries

end module grid_module
