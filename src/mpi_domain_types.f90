module mpi_domain_types
   use mpi_f08, only: MPI_Cart_coords, MPI_Cart_create, MPI_Cart_shift, MPI_Comm_rank, MPI_Comm_size, MPI_Dims_create, MPI_Cart_rank
   use mpi_f08, only: MPI_COMM_WORLD, MPI_SUCCESS, MPI_Comm, MPI_PROC_NULL, MPI_Abort
   implicit none(type, external)
   private

   type, public :: mpi_domain_t
      private ! Make components private by default
      type(MPI_Comm) :: comm
         !! Cartesian communicator
      integer :: rank = -1
         !! Rank in Cartesian communicator
      integer :: size = -1
         !! Size of Cartesian communicator
      integer :: ndims = 3
         !! Number of dimensions
      integer :: dims(3) = 0
         !! Core decomposition per dimension
      integer :: coords(3) = -1
         !! Coordinates of this rank
      logical :: periodic(3) = .false.
         !! Periodicity per dimension
      integer :: neighbors(6) = MPI_PROC_NULL
         !! Ranks of [W, E, S, N, L, H] neighbors
      integer :: neighbors_extended(-1:1, -1:1, -1:1) = MPI_PROC_NULL
      logical :: reorder = .true.
         !! Core ranking may be reordered (`true`) or not (`false`)
      logical, public :: is_boundary_face(6) = .false.
         !! `True` if face is a physical boundary
      logical, public :: is_interior = .true.
         !! `True` if rank has no physical boundary faces
   contains
      procedure, public :: initialize => initialize_mpi_domain
      procedure, public :: get_communicator => get_domain_communicator
      procedure, public :: get_rank => get_domain_rank
      procedure, public :: get_neighbors => get_domain_neighbors
      procedure, public :: get_size => get_domain_size
      procedure, public :: abort => abort_mpi_processes
      procedure, public :: determine_extended_neighbors
      procedure, public :: test_extended_neighbors

      procedure, private :: determine_neighbors
      procedure, private :: set_periodicity
      procedure, private :: check_physical_boundaries
   end type mpi_domain_t

   interface
      module subroutine determine_extended_neighbors(self)
         class(mpi_domain_t), intent(inout) :: self
      end subroutine determine_extended_neighbors

      module subroutine test_extended_neighbors(self)
         class(mpi_domain_t), intent(in) :: self
      end subroutine test_extended_neighbors
   end interface
contains

   !> Subroutine to initialize the type instance
   module subroutine initialize_mpi_domain(self, requested_dims, boundary_conditions, parent_comm)
      class(mpi_domain_t), intent(inout) :: self
      integer, intent(in) :: requested_dims(3)
      integer, intent(in) :: boundary_conditions(6)
      type(MPI_Comm), intent(in), optional :: parent_comm

      type(MPI_Comm) :: comm_parent
      integer :: ierr, parent_size

      comm_parent = MPI_COMM_WORLD
      if (present(parent_comm)) comm_parent = parent_comm

      ! 1. Get original communicator size
      call MPI_Comm_size(comm_parent, parent_size, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Comm_size"

      ! 2. Determine dimensions (MPI_Dims_create)
      ! Creates a division of cores in a Cartesian ndims-dimensional X(, Y(, Z)) grid
      self%ndims = 3 ! Assuming 3D for now
      self%dims = requested_dims
      call MPI_Dims_create(nnodes=parent_size, ndims=self%ndims, dims=self%dims, ierror=ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error creating dimensions"

      ! 3. Determine periodicity from inputs
      call self%set_periodicity(boundary_conditions)

      ! 4. Create Cartesian communicator
      call MPI_Cart_create(comm_parent, self%ndims, self%dims, self%periodic, self%reorder, self%comm, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error creating Cartesian communicator"

      ! 5. Get rank, size, and coordinates in the new communicator
      call MPI_Comm_rank(self%comm, self%rank, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error deciding ranks in Cartesian"
      call MPI_Comm_size(self%comm, self%size, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Comm_size"
      if (self%size /= parent_size) print *, "Size mismatch between old and new communicators"
      call MPI_Cart_coords(self%comm, self%rank, self%ndims, self%coords, ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error getting cartesian coordinates"

      ! 6. Determine neighbors
      call self%determine_neighbors()

      ! 7. Determine which faces are physical boundaries
      call self%check_physical_boundaries()

      if (self%rank == 0) then
         print *, "MPI Domain Initialized. Dimensions:", self%dims
      end if

   end subroutine initialize_mpi_domain

   !> Sets the periodic boundaries for the MPI Cartesian communicator, based on an array of boundary types.
   module subroutine set_periodicity(self, bc_types)
      use enums, only: X_DIR, Y_DIR, Z_DIR, D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH, PERIODIC
      class(mpi_domain_t), intent(inout) :: self
      integer, intent(in) :: bc_types(6)

      ! `X_DIR`,..,`Z_DIR` is 0-based index where `is_periodic` is 1-based
      ! An axis is periodic if both directions are periodic
      self%periodic = .false.
      self%periodic(X_DIR + 1) = (bc_types(D_WEST) == PERIODIC .and. bc_types(D_EAST) == PERIODIC)
      self%periodic(Y_DIR + 1) = (bc_types(D_SOUTH) == PERIODIC .and. bc_types(D_NORTH) == PERIODIC)
      self%periodic(Z_DIR + 1) = (bc_types(D_LOW) == PERIODIC .and. bc_types(D_HIGH) == PERIODIC)
   end subroutine set_periodicity

   !> Finds the ranks of the 6 nearest neighbors by shifting ±1 in X, Y, and Z.
   module subroutine determine_neighbors(self)
      use enums, only: X_DIR, Y_DIR, Z_DIR
      class(mpi_domain_t), intent(in out) :: self
      integer :: ierr, west, east, south, north, low, high
      ! Logic from original get_neighbors
      call MPI_Cart_shift(comm=self%comm, direction=X_DIR, disp=1, rank_source=west, rank_dest=east, ierror=ierr)
      call MPI_Cart_shift(comm=self%comm, direction=Y_DIR, disp=1, rank_source=south, rank_dest=north, ierror=ierr)
      call MPI_Cart_shift(comm=self%comm, direction=Z_DIR, disp=1, rank_source=low, rank_dest=high, ierror=ierr)
      if (ierr /= MPI_SUCCESS) error stop "Error in Cartesian shift"

      self%neighbors = [west, east, south, north, low, high]
   end subroutine determine_neighbors

   module subroutine check_physical_boundaries(self)
      class(mpi_domain_t), intent(inout) :: self
      ! Logic from original determine_rank_boundaries
      self%is_boundary_face = .false.
      where (self%neighbors == MPI_PROC_NULL)
         self%is_boundary_face = .true.
      end where
      self%is_interior = .not. any(self%is_boundary_face)
   end subroutine check_physical_boundaries

   ! --- Implement simple getter functions ---
   module function get_domain_communicator(self) result(comm)
      class(mpi_domain_t), intent(in) :: self
      type(MPI_Comm) :: comm
      comm = self%comm
   end function get_domain_communicator

   module function get_domain_rank(self) result(rank)
      class(mpi_domain_t), intent(in) :: self
      integer :: rank
      rank = self%rank
   end function get_domain_rank

   module function get_domain_neighbors(self) result(neighbors_array)
      class(mpi_domain_t), intent(in) :: self
      integer :: neighbors_array(6)
      neighbors_array = self%neighbors
   end function get_domain_neighbors

   module function get_domain_size(self) result(comm_size)
      class(mpi_domain_t), intent(in) :: self
      integer :: comm_size
      comm_size = self%size
   end function get_domain_size

   !> Aborts the MPI processes cleanly
   module subroutine abort_mpi_processes(self, msg)
      class(mpi_domain_t), intent(in) :: self
      character(len=*), intent(in) :: msg
      integer :: ierr
      print *, msg
      call MPI_Abort(self%get_communicator(), ierr)
   end subroutine abort_mpi_processes
end module mpi_domain_types
