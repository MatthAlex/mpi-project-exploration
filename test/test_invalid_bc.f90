program test_invalid_bc
   use mpi_f08, only: MPI_Init, MPI_Finalize
   use mpi_domain_types, only: mpi_domain_t
   use lib_mpi_enums, only: PERIODIC, DIRICHLET
   implicit none

   type(mpi_domain_t) :: domain
   integer :: ierr
   integer :: invalid_bc(6)

   call MPI_Init(ierror=ierr)

   ! This should abort
   invalid_bc = [PERIODIC, DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]
   call domain%initialize([0, 0, 0], invalid_bc)

   ! Should never reach here
   error stop "TEST FAILED: Invalid BC was not caught"

end program test_invalid_bc
