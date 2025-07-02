program test_simple_mpi_f08
   use mpi_f08, only: MPI_Init, MPI_SUCCESS, MPI_Status, MPI_Sendrecv, MPI_INTEGER
   use mpi_f08, only: MPI_Finalize, MPI_Barrier
   use mpi_f08, only: MPI_Comm_size, MPI_Comm_rank, MPI_Comm
   use mpi_domain_types, only: mpi_domain_t
   use lib_mpi_enums, only: D_WEST, D_EAST
   implicit none(type, external)

   integer :: ierr
   integer :: rank, size
   type(MPI_Status) :: status
   integer :: left, right, sendval, recvval
   integer, parameter :: tag = 1
   type(mpi_domain_t) :: domain
   type(MPI_Comm) :: comm
   integer :: neighbors(6)

   call MPI_Init(ierr)

   ! Initialize MPI by letting it choose the size
   ! Also set periodic boundaries in X
   call domain%initialize([4, 0, 0], [0, 0, 1, 1, 2, 2])

   ! Get only the required information
   rank = domain%get_rank()
   size = domain%get_size()
   comm = domain%get_communicator()
   neighbors = domain%get_neighbors()

   left = neighbors(D_WEST)
   right = neighbors(D_EAST)

   ! Old periodic boundary implementation
   ! ! Periodic boundary
   ! left = rank - 1
   ! if (left < 0) left = size - 1
   ! right = rank + 1
   ! if (right >= size) right = 0

   ! Print information about this process
   print '(3(a,1x,i0),1x,i0)', 'Hello from process ', rank, ' of ', size, ". My neighbours are:", left, right

   ! Data to send is the rank of the current process
   sendval = rank

   ! Perform Sendrecv operation: send rank to right neighbor, receive from left neighbor
   call MPI_Sendrecv(sendval, 1, MPI_INTEGER, right, tag, &  ! Send to right
                     recvval, 1, MPI_INTEGER, left, tag, &  ! Receive from left
                     comm, status, ierr)
   if (ierr /= MPI_SUCCESS) call domain%abort("Sendrecv failed")
   if (recvval /= left) call domain%abort("Wrong value received from left")

! Perform Sendrecv operation: send rank to left neighbor, receive from right neighbor
   call MPI_Sendrecv(sendval, 1, MPI_INTEGER, left, tag, &  ! Send to left
                     recvval, 1, MPI_INTEGER, right, tag, &  ! Receive from right
                     comm, status, ierr)
   if (ierr /= MPI_SUCCESS) call domain%abort("Sendrecv failed")
   if (recvval /= right) call domain%abort("Wrong value received from right")

   ! Synchronize all processes
   call MPI_Barrier(comm, ierr)
   if (ierr /= MPI_SUCCESS) call domain%abort("Barrier failed")

   ! Finalize MPI Environment
   call MPI_Finalize(ierr)
   if (ierr /= MPI_SUCCESS) call domain%abort("MPI finalization failed")

end program test_simple_mpi_f08
