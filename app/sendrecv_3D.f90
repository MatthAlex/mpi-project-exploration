!> Extension of the simple 1D MPI implementation.
!> - It leaves the topology generation to MPI
!> - Can handle periodic boundaries with a boolean flag
!> - Heuristics may reorder the original comm (CPU topology) into a more efficient one
!> - Figures out who the neighbours are of each rank automatically
program sendrecv_3D
   use mpi
   use precision, only: sp
   use grid_module, only: initialize_MPI_grid, rank => my_rank, comm_cart
   use lib_parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z, iterations
   use lib_mpi_halo, only: update_mpi_halo
   use test_halo, only: check_halo_real
   use test_boundary, only: check_boundary_real
   use boundary, only: determine_rank_boundaries, apply_boundaries
   implicit none

   integer :: ierr
   real(kind=sp), allocatable :: array(:, :, :)
   integer :: i

   call MPI_Init(ierr)
   call initialize_MPI_grid()

   call determine_rank_boundaries()

   ! Create a local array with halo regions
   allocate (array(0:nx + 1, 0:ny + 1, 0:nz + 1), source=real(rank, kind=sp))

   do i = 1, iterations

      call update_mpi_halo(array=array)
      call apply_boundaries(array=array, bc_types=boundaries)

      call check_halo_real(array=array)
      call check_boundary_real(array=array)

      call MPI_Barrier(comm=comm_cart, ierror=ierr)
   end do

   ! Since the tests above pass
   if (rank == 0) print*, "Success!"

   call MPI_Finalize(ierr)

end program sendrecv_3D
