program test_mpi_extended_neighbors
   use mpi_f08, only: MPI_Init, MPI_Finalize
   use mpi_domain_types, only: mpi_domain_t
   implicit none(type, external)

   type(mpi_domain_t) :: domain
   integer :: ierr, my_rank, my_size, i

   call MPI_Init(ierror=ierr)
   call domain%initialize(requested_dims=[0, 0, 0], boundary_conditions=[0, 0, 1, 1, 2, 2])
   my_rank = domain%get_rank()
   my_size = domain%get_size()
   call domain%determine_extended_neighbors()
   do i = 0, my_size-1
      if (my_rank == i) call domain%test_extended_neighbors()
   end do
   call MPI_Finalize(ierror=ierr)

end program test_mpi_extended_neighbors
