program test_invalid_bc
   use lib_mpi_enums, only: PERIODIC, DIRICHLET
   use mpi_domain_types, only: mpi_domain_t
   use mpi_f08, only: MPI_Init, MPI_Finalize
   implicit none

   type(mpi_domain_t) :: domain
   integer :: invalid_bc(6)

   call MPI_Init()

   ! This should abort
   invalid_bc = [PERIODIC, DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]
   call domain%initialize([0, 0, 0], invalid_bc)

   ! Should never reach here
   call domain%abort("TEST FAILED: Invalid BC was not caught")

end program test_invalid_bc
