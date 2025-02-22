!> Extension of the simple 1D MPI implementation.
!> - It leaves the topology generation to MPI
!> - Can handle periodic boundaries with a boolean flag
!> - Heuristics may reorder the original comm (CPU topology) into a more efficient one
!> - Figures out who the neighbours are of each rank automatically
program sendrecv_3D
   use mpi_f08, only: MPI_Init, MPI_Barrier, MPI_Finalize, MPI_WTime
   use precision, only: sp, dp
   use mpi_domain, only: create_mpi_domain, rank => my_rank, comm_cart
   use parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z, iterations, boundaries
   use mpi_halo, only: update_mpi_halo
   use test_halo, only: check_halo_real
   use test_boundary, only: check_boundary_real
   use boundary, only: determine_rank_boundaries, update_boundary_conditions
   implicit none (type, external)

   integer :: ierr
   real(kind=sp), allocatable :: array(:, :, :), array_smol(:, :, :)
   integer :: i
   real(dp) :: start, finish

   call MPI_Init(ierr)
   call create_mpi_domain()

   call determine_rank_boundaries()

   ! Create a local array with halo regions
   allocate (array(0:nx + 1, 0:ny + 1, 0:nz + 1), source=real(rank, kind=sp))
   ! Create a local array with halo regions of different size
   allocate (array_smol(0:nx/4 + 1, 0:ny/4 + 1, 0:nz/4 + 1), source=real(rank, kind=sp))

   start = MPI_WTime()
   do i = 1, iterations

      call update_mpi_halo(array=array)
      call update_boundary_conditions(array=array, bc_types=boundaries)

      call check_halo_real(array=array)
      call check_boundary_real(array=array)

      call update_mpi_halo(array=array_smol)
      call update_boundary_conditions(array=array_smol, bc_types=boundaries)

      call check_halo_real(array=array_smol)
      call check_boundary_real(array=array_smol)

      call MPI_Barrier(comm=comm_cart, ierror=ierr)
   end do
   finish = MPI_WTime()
   if (rank == 0) print '(A,1X,F8.2,A)', "MPI TIME:", finish - start, "s"

   ! Since the tests above pass
   if (rank == 0) print *, "Success!"

   call MPI_Finalize(ierr)

end program sendrecv_3D
