module mpi_domain_types
   use lib_mpi_enums, only: X_DIR, Y_DIR, Z_DIR, D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH, PERIODIC
   use mpi_f08, only: MPI_Cart_coords, MPI_Cart_create, MPI_Cart_shift, MPI_Comm_rank, MPI_Comm_size, MPI_Dims_create, MPI_Cart_rank
   use mpi_f08, only: MPI_COMM_WORLD, MPI_SUCCESS, MPI_Comm, MPI_PROC_NULL, MPI_Abort
   implicit none(type, external)
   private

   integer, parameter :: ABORT_ERRORCODE = 111
      !! Errorcode highlighting internal library failures.

   type, public :: mpi_domain_t
      private
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
      logical :: reorder = .true.
         !! Allow MPI to reorder ranks for topology-aware placement (NUMA, socket affinity).
      logical, public :: is_boundary_face(6) = .false.
         !! `True` if face is a physical boundary. Public to avoid getter overhead in hot loops.
      logical, public :: is_interior = .true.
         !! `True` if rank has no physical boundary faces. Public to avoid getter overhead in hot loops.
   contains
      procedure, public :: initialize => initialize_mpi_domain
      procedure, public :: get_communicator => get_domain_communicator
      procedure, public :: get_rank => get_domain_rank
      procedure, public :: get_neighbors => get_domain_neighbors
      procedure, public :: get_size => get_domain_size
      procedure, public :: get_coords => get_decomposition_coords
      procedure, public :: get_dims => get_domain_dims
      procedure, public :: get_periodicity => get_periodic_dims
      procedure, public :: coords_to_rank
      procedure, public :: abort => abort_mpi_processes

      procedure, private :: determine_neighbors
      procedure, private :: set_periodicity
      procedure, private :: check_physical_boundaries
      procedure :: log_message => domain_log_message

   end type mpi_domain_t

contains

   module subroutine initialize_mpi_domain(self, requested_dims, boundary_conditions, parent_comm)
      !! Initializes the MPI Cartesian domain instance.
      !!
      !! The initialization order matters:
      !! 1. The parent communicator is seeded into `comm` immediately so that `abort()`
      !!    works even if a later step fails.
      !! 2. `MPI_Dims_create` resolves the core decomposition from the requested dims.
      !! 3. Periodicity is derived from the boundary conditions and validated for
      !!    per-axis symmetry before the Cartesian communicator is created.
      !! 4. After communicator creation, rank metadata, neighbors, and physical
      !!    boundary faces are determined.
      class(mpi_domain_t), intent(inout) :: self
      integer, intent(in) :: requested_dims(3)
      integer, intent(in) :: boundary_conditions(6)
      type(MPI_Comm), intent(in), optional :: parent_comm

      type(MPI_Comm) :: comm_parent, comm_cart
      integer :: ierr, parent_size, rc

      comm_parent = MPI_COMM_WORLD
      if (present(parent_comm)) comm_parent = parent_comm
      self%comm = comm_parent

      call MPI_Comm_size(comm_parent, parent_size, ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Comm_size failed..")

      ! Decomposes cores in a Cartesian 3-dimensional X(, Y(, Z)) grid
      self%ndims = 3 ! Hardcoded to 3D problems - extending it to 2D/1D not in scope
      self%dims = requested_dims
      call MPI_Dims_create(nnodes=parent_size, ndims=self%ndims, dims=self%dims, ierror=ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Dims_create failed..")

      rc = validate_periodic_bcs_symmetric(boundary_conditions)
      if (rc == -1) call self%abort("ERROR: MPI: Invalid BC: X-axis periodic must be set on both West and East, or neither")
      if (rc == -2) call self%abort("ERROR: MPI: Invalid BC: Y-axis periodic must be set on both South and North, or neither")
      if (rc == -3) call self%abort("ERROR: MPI: Invalid BC: Z-axis periodic must be set on both Low and High, or neither")

      call self%set_periodicity(boundary_conditions)

      call MPI_Cart_create(comm_parent, self%ndims, self%dims, self%periodic, self%reorder, comm_cart, ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Cart_create failed..")

      self%comm = comm_cart

      call MPI_Comm_rank(self%comm, self%rank, ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Comm_rank failed..")
      call MPI_Comm_size(self%comm, self%size, ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Comm_size failed..")
      call MPI_Cart_coords(self%comm, self%rank, self%ndims, self%coords, ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Cart_coords failed..")

      call self%determine_neighbors()

      call self%check_physical_boundaries()

   end subroutine initialize_mpi_domain

   module subroutine set_periodicity(self, bc_types)
      !! Derives per-axis periodicity from the six face boundary conditions by setting `self%periodic(3)`.
      !! An axis is periodic only if both opposing faces are set to `PERIODIC`.
      !! Note: `X_DIR`/`Y_DIR`/`Z_DIR` are 0-based (MPI convention); `periodic` is 1-based.
      class(mpi_domain_t), intent(inout) :: self
      integer, intent(in) :: bc_types(6)

      ! `X_DIR`,..,`Z_DIR` is 0-based index where `is_periodic` is 1-based
      ! An axis is periodic if both directions are periodic
      self%periodic = .false.
      self%periodic(X_DIR + 1) = (bc_types(D_WEST) == PERIODIC .and. bc_types(D_EAST) == PERIODIC)
      self%periodic(Y_DIR + 1) = (bc_types(D_SOUTH) == PERIODIC .and. bc_types(D_NORTH) == PERIODIC)
      self%periodic(Z_DIR + 1) = (bc_types(D_LOW) == PERIODIC .and. bc_types(D_HIGH) == PERIODIC)
   end subroutine set_periodicity

   module subroutine determine_neighbors(self)
      !! Finds the ranks of the 6 nearest neighbors by shifting ±1 in X, Y, and Z. Populates `self%neighbors(6)`.
      class(mpi_domain_t), intent(in out) :: self
      integer :: ierr, west, east, south, north, low, high
      call MPI_Cart_shift(comm=self%comm, direction=X_DIR, disp=1, rank_source=west, rank_dest=east, ierror=ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Cart_shift failed..")
      call MPI_Cart_shift(comm=self%comm, direction=Y_DIR, disp=1, rank_source=south, rank_dest=north, ierror=ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Cart_shift failed..")
      call MPI_Cart_shift(comm=self%comm, direction=Z_DIR, disp=1, rank_source=low, rank_dest=high, ierror=ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Cart_shift failed..")

      self%neighbors = [west, east, south, north, low, high]
   end subroutine determine_neighbors

   pure module subroutine check_physical_boundaries(self)
      !! Marks which faces are physical boundaries and whether this rank is fully interior.
      !! A face is a physical boundary when its neighbor is `MPI_PROC_NULL`.
      !! A rank is interior when none of its six faces are physical boundaries —
      !! this is the fast path that skips boundary condition application entirely.
      !! Populates `self%is_boundary_face` and `self%is_interior`.
      class(mpi_domain_t), intent(inout) :: self
      self%is_boundary_face = .false.
      where (self%neighbors == MPI_PROC_NULL)
         self%is_boundary_face = .true.
      end where
      self%is_interior = .not. any(self%is_boundary_face)
   end subroutine check_physical_boundaries

   pure module function get_domain_communicator(self) result(comm)
      class(mpi_domain_t), intent(in) :: self
      type(MPI_Comm) :: comm
      comm = self%comm
   end function get_domain_communicator

   pure module function get_domain_rank(self) result(rank)
      class(mpi_domain_t), intent(in) :: self
      integer :: rank
      rank = self%rank
   end function get_domain_rank

   pure module function get_domain_neighbors(self) result(neighbors_array)
      class(mpi_domain_t), intent(in) :: self
      integer :: neighbors_array(6)
      neighbors_array = self%neighbors
   end function get_domain_neighbors

   pure module function get_domain_size(self) result(comm_size)
      class(mpi_domain_t), intent(in) :: self
      integer :: comm_size
      comm_size = self%size
   end function get_domain_size

   pure module function get_decomposition_coords(self) result(coordinates)
      class(mpi_domain_t), intent(in) :: self
      integer :: coordinates(3)
      coordinates = self%coords
   end function get_decomposition_coords

   pure module function get_domain_dims(self) result(requested_dims)
      class(mpi_domain_t), intent(in) :: self
      integer :: requested_dims(3)
      requested_dims = self%dims
   end function get_domain_dims

   pure module function get_periodic_dims(self) result(periodic_dims)
      class(mpi_domain_t), intent(in) :: self
      logical :: periodic_dims(3)
      periodic_dims = self%periodic
   end function get_periodic_dims

   module function coords_to_rank(self, coords) result(rank)
      class(mpi_domain_t), intent(in) :: self
      integer, intent(in) :: coords(3)
      integer :: rank
      integer :: ierr
      call MPI_Cart_rank(self%comm, coords, rank, ierr)
      if (ierr /= MPI_SUCCESS) call self%abort("ERROR: MPI: Cart_rank failed..")
   end function coords_to_rank

   module subroutine abort_mpi_processes(self, msg, errorcode)
      !! Logs a message and aborts all MPI processes.
      !! Uses `ABORT_ERRORCODE` (111) by default; callers may pass a specific error code
      !! via the optional `errorcode` argument for finer-grained diagnostics.
      class(mpi_domain_t), intent(in) :: self
      character(len=*), intent(in) :: msg
      integer, intent(in), optional :: errorcode
      integer :: errcode_
      call self%log_message(msg)
      errcode_ = ABORT_ERRORCODE
      if (present(errorcode)) errcode_ = errorcode
      call MPI_Abort(comm=self%comm, errorcode=errcode_)
   end subroutine abort_mpi_processes

   subroutine domain_log_message(self, msg)
      class(mpi_domain_t), intent(in) :: self
      character(len=*), intent(in) :: msg
      character(len=256) :: formatted_msg

      write (formatted_msg, "(A,I0,A,A)") "[", self%rank, "] ", trim(msg)
      print *, trim(formatted_msg)
   end subroutine domain_log_message

   pure function validate_periodic_bcs_symmetric(bc_types) result(rc)
      !! Returns 0 if periodic BCs are symmetric per axis, or a negative axis code if not.
      !! Only checks symmetry; does not validate that individual BC values are in range.
      integer, intent(in) :: bc_types(6)
      integer :: rc
      rc = 0

      ! X-axis
      if ((bc_types(D_WEST) == PERIODIC) .neqv. (bc_types(D_EAST) == PERIODIC)) then
         rc = -1
         return
      end if

      ! Y-axis
      if ((bc_types(D_SOUTH) == PERIODIC) .neqv. (bc_types(D_NORTH) == PERIODIC)) then
         rc = -2
         return
      end if

      ! Z-axis
      if ((bc_types(D_LOW) == PERIODIC) .neqv. (bc_types(D_HIGH) == PERIODIC)) then
         rc = -3
         return
      end if

   end function validate_periodic_bcs_symmetric
end module mpi_domain_types
