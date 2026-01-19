!> Extension of the simple 1D MPI implementation.
!> - It leaves the topology generation to MPI
!> - Can handle periodic boundaries with a boolean flag
!> - Heuristics may reorder the original comm (CPU topology) into a more efficient one
!> - Figures out who the neighbours are of each rank automatically
program test_sendrecv_3D
   use mpi_f08, only: MPI_Init, MPI_Barrier, MPI_Finalize, MPI_WTime, MPI_Comm
   use lib_mpi_precision, only: sp, dp
   use mpi_domain_types, only: mpi_domain_t
   use lib_mpi_parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z, iterations, boundaries, core_decomposition
   use mpi_halo, only: update_mpi_halo
   use test_halo, only: check_halo_real
   use test_boundary, only: check_boundary_real
   use boundary, only: update_boundaries
   implicit none(type, external)

   type(mpi_domain_t) :: domain
   type(MPI_Comm) :: comm_cart
   integer :: ierr
   real(kind=sp), allocatable :: array(:, :, :), array_smol(:, :, :)
   integer :: i
   integer :: rank
   real(dp) :: start, finish
   real(sp) :: dir_val(6)
   !! Dirichlet values

   call MPI_Init(ierror=ierr)
   call domain%initialize(core_decomposition, boundaries)
   rank = domain%get_rank()
   comm_cart = domain%get_communicator()

   dir_val = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp, 5.0_sp, 6.0_sp]

   ! Create a local array with halo regions
   allocate (array(0:nx + 1, 0:ny + 1, 0:nz + 1), source=real(rank, kind=sp))
   ! Create a local array with halo regions of different size
   allocate (array_smol(0:nx / 4 + 1, 0:ny / 4 + 1, 0:nz / 4 + 1), source=real(rank, kind=sp))

   start = MPI_WTime()
   do i = 1, iterations
      if (rank == 0) print *, "Array"
      call update_mpi_halo(domain=domain, array=array)
      call update_boundaries(domain=domain, array=array, bc_types=boundaries, dirichlet_values=dir_val)

      call check_halo_real(domain=domain, array=array)
      call check_boundary_real(domain=domain, array=array, bc_types=boundaries, dirichlet_values=dir_val)

      call MPI_Barrier(comm=comm_cart, ierror=ierr)
      if (rank == 0) print *, "Array smol"
      call update_mpi_halo(domain=domain, array=array_smol)
      call update_boundaries(domain=domain, array=array_smol, bc_types=boundaries, dirichlet_values=dir_val)

      call check_halo_real(domain=domain, array=array_smol)
      call check_boundary_real(domain=domain, array=array_smol, bc_types=boundaries, dirichlet_values=dir_val)

      call MPI_Barrier(comm=comm_cart, ierror=ierr)
   end do
   finish = MPI_WTime()
   if (rank == 0) print '(A,1X,F8.2,A)', "MPI TIME:", finish - start, "s"

   ! Since the tests above pass
   if (rank == 0) print *, "Success!"

   call MPI_Finalize(ierror=ierr)

end program test_sendrecv_3D
