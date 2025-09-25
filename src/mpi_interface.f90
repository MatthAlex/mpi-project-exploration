module mpi_interface
   !! This module acts an an intermediate layer for legacy codes that can't readily integrate the OOP object.
   use mpi_f08, only: MPI_Comm, MPI_Cart_rank, MPI_SUCCESS
   use mpi_domain_types, only: mpi_domain_t
   implicit none(type, external)
   private
   public :: create_mpi_domain, get_mpi_domain_handle

   type(mpi_domain_t), target :: mpi_domain
   !! The single, authoritative instance the MPI object

   integer, allocatable, public :: PTop(:, :, :)
   !! Legacy Processor Topology - forced to support

   integer, parameter, public :: NDIM = 3
   !! Number of dimensions
   integer, parameter, public :: MASTER = 0
   !! Master Rank 0
   integer, protected, public :: myrank, num_ranks
   integer, protected, public :: mpi_coords(NDIM)
   !! Rank coordinates
   integer, protected, public :: mpi_dims(NDIM)
   !! Decomposition dimensions
   type(MPI_Comm), protected, public :: cart_comm
   !! The cartesian communicator

contains

   subroutine create_mpi_domain(user_specified_dims, boundary_conditions)
      !! Initializes the MPI Cartesian decomposition. This is the ONLY setup routine the legacy code should call.
      integer, intent(in) :: user_specified_dims(NDIM)
      ! These are the legacy [nprocx, nprocy, nprocz]
      integer, intent(in) :: boundary_conditions(NDIM * 2)
      integer :: i, j, k, ierr
      integer :: coords_buffer(0:NDIM - 1) ! MPI uses 0-based indexing for this call

      ! Step 1: Initialize our modern, reliable MPI domain object.
      ! Let MPI_Dims_create figure out the decomposition if any are 0.
      call mpi_domain%initialize(requested_dims=user_specified_dims, &
            boundary_conditions=boundary_conditions)

      ! Get MPI attributes from MPI domain.
      myrank = mpi_domain%get_rank()
      num_ranks = mpi_domain%get_size()
      cart_comm = mpi_domain%get_communicator()
      mpi_dims = mpi_domain%get_dims()
      mpi_coords = mpi_domain%get_coords()

      ! Step 2: Build the legacy PTop array as a compatibility layer.
      ! The dimensions are now known from the initialized mpi_domain.
      ! PTop is 1-based coords system.
      allocate(PTop(mpi_dims(1), mpi_dims(2), mpi_dims(3)))

      ! For every position in the process grid...
      do k = 1, mpi_dims(3)
         do j = 1, mpi_dims(2)
            do i = 1, mpi_dims(1)
               ! ...ask MPI what rank lives at these 0-based coordinates.
               coords_buffer = [i - 1, j - 1, k - 1]
               call MPI_Cart_rank(cart_comm, coords_buffer, PTop(i, j, k), ierr)
               if (ierr /= MPI_SUCCESS) call mpi_domain%abort("Failed to build PTop array")
            end do
         end do
      end do
   end subroutine create_mpi_domain

   function get_mpi_domain_handle() result(domain_ptr)
      !! Returns a pointer to the single, authoritative instance of the MPI domain.
      type(mpi_domain_t), pointer :: domain_ptr

      ! Point the result directly to our private, module-level instance.
      domain_ptr => mpi_domain

   end function get_mpi_domain_handle
end module mpi_interface
