!> Extension of the simple 1D MPI implementation.
!> - It leaves the topology generation to MPI
!> - Can handle periodic boundaries with a boolean flag
!> - Heuristics may reorder the original comm (CPU topology) into a more efficient one
!> - Figures out who the neighbours are of each rank automatically
program sendrecv_3D
   use mpi
   use grid_module, only: initialize_MPI_grid, rank => my_rank
   use lib_parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z, iterations
   use lib_mpi_halo, only: update_mpi_halo
   implicit none

   integer :: ierr
   real(kind=4), allocatable :: array(:, :, :)
   integer :: i

   call initialize_MPI_grid()

   ! Create a local array with halo regions
   allocate (array(0:nx + 1, 0:ny + 1, 0:nz + 1), source=real(rank, kind=4))

   do i = 1, iterations

      call update_mpi_halo(array=array)

   end do

   call MPI_Finalize(ierr)

end program sendrecv_3D
