program simple_mpi_f08
   use mpi_f08, only: MPI_Init, MPI_SUCCESS, MPI_COMM_WORLD, MPI_Status, MPI_Sendrecv, MPI_INTEGER
   use mpi_f08, only: MPI_Comm_size, MPI_Comm_rank
   implicit none

   integer :: ierr
   integer :: rank, size
   type(MPI_Status) :: status
   integer :: left, right, sendval, recvval
   integer, parameter :: tag = 1

   ! Initialize MPI environment
   call MPI_Init(ierr)
   if (ierr /= MPI_SUCCESS) error stop "MPI initialization failed"

   ! Get the number of processes
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   if (ierr /= MPI_SUCCESS) error stop "Failed to get communicator size"

   ! Get the rank of this process
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   if (ierr /= MPI_SUCCESS) error stop "Failed to get rank"

   ! Periodic boundary
   left = rank - 1
   if (left < 0) left = size - 1
   right = rank + 1
   if (right >= size) right = 0

   ! Print information about this process
   print '(3(a,1x,i0),1x,i0)', 'Hello from process ', rank, ' of ', size, ". My neighbours are:", left, right

   ! Data to send is the rank of the current process
   sendval = rank

   ! Perform Sendrecv operation: send rank to right neighbor, receive from left neighbor
   call MPI_Sendrecv(sendval, 1, MPI_INTEGER, right, tag, & ! Send to right
                     recvval, 1, MPI_INTEGER, left, tag, &  ! Receive from left
                     MPI_COMM_WORLD, status, ierr)
   if (ierr /= MPI_SUCCESS) error stop "Sendrecv failed"
   if (recvval /= left) error stop "Wrong value received from left"

! Perform Sendrecv operation: send rank to left neighbor, receive from right neighbor
   call MPI_Sendrecv(sendval, 1, MPI_INTEGER, left, tag, & ! Send to left
                     recvval, 1, MPI_INTEGER, right, tag, &  ! Receive from right
                     MPI_COMM_WORLD, status, ierr)
   if (ierr /= MPI_SUCCESS) error stop "Sendrecv failed"
   if (recvval /= right) error stop "Wrong value received from right"

   ! Synchronize all processes
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) error stop "Barrier failed"

   ! Finalize MPI Environment
   call MPI_Finalize(ierr)
   if (ierr /= MPI_SUCCESS) error stop "MPI finalization failed"

end program simple_mpi_f08
