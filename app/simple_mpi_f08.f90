program simple_mpi_f08
   use mpi_f08
   implicit none

   integer :: ierr
   integer :: rank, size
   type(MPI_Status) :: status

   ! Initialize MPI environment
   call MPI_Init(ierr)
   if (ierr /= MPI_SUCCESS) error stop "MPI initialization failed"

   ! Get the number of processes
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   if (ierr /= MPI_SUCCESS) error stop "Failed to get communicator size"

   ! Get the rank of this process
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   if (ierr /= MPI_SUCCESS) error stop "Failed to get rank"

   ! Print information about this process
   print '(a,i0,a,i0)', 'Hello from process ', rank, ' of ', size

   ! Synchronize all processes
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) error stop "Barrier failed"

   ! Finalize MPI Environment
   call MPI_Finalize(ierr)
   if (ierr /= MPI_SUCCESS) error stop "MPI finalization failed"

end program simple_mpi_f08
