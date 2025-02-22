!> Encapsulates MPI halo update functionality
!> - static transfer buffers - avoid allocation every update
!> - module procedure to generalise subroutine calling
!> - independent dimensionality
module lib_mpi_halo
   use mpi_f08, only: MPI_Sendrecv, MPI_STATUS, MPI_REAL, MPI_INTEGER
   use precision, only: sp
   use grid_module, only: west, east, south, north, low, high, comm_cart
   use lib_parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   implicit none(type, external)
   private
   public :: update_mpi_halo, init_arrays

   integer :: ierr
   type(MPI_Status) :: status

   !> Buffers for integer transfers - avoids allocation
   integer, dimension(ny + 2, nz + 2) :: buffer_send_x_int, buffer_rcv_x_int
   integer, dimension(nx + 2, nz + 2) :: buffer_send_y_int, buffer_rcv_y_int
   integer, dimension(nx + 2, ny + 2) :: buffer_send_z_int, buffer_rcv_z_int

   !> Buffers for real transfers - avoids reallocation
   real(kind=sp), allocatable :: buffer_send_x(:, :), buffer_rcv_x(:, :)
   real(kind=sp), allocatable :: buffer_send_y(:, :), buffer_rcv_y(:, :)

   !> Direction face-surface sizes
   integer, parameter :: X_FACE_SIZE = (ny + 2) * (nz + 2)
   integer, parameter :: Y_FACE_SIZE = (nx + 2) * (ny + 2)
   integer, parameter :: Z_FACE_SIZE = (nx + 2) * (ny + 2)

   integer, parameter :: tag = 1

   interface update_mpi_halo
      module procedure update_mpi_halo_real
      module procedure update_mpi_halo_integer
   end interface

contains

   subroutine init_arrays()
      allocate (buffer_send_x(ny + 2, nz + 2), buffer_rcv_x(ny + 2, nz + 2))
      allocate (buffer_send_y(nx + 2, nz + 2), buffer_rcv_y(nx + 2, nz + 2))
   end subroutine init_arrays

   !> Perform halo exchanges with nearest neighbors, in all 3 directions, for an array of kind sp
   subroutine update_mpi_halo_real(array)
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x(:, :) = array(ubound(array, 1) - 1, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_REAL, east, tag, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_REAL, west, tag, &
                        comm_cart, status, ierr)
      array(lbound(array, 1), :, :) = buffer_rcv_x

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y(:, :) = array(:, ubound(array, 2) - 1, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_REAL, north, tag, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_REAL, south, tag, &
                        comm_cart, status, ierr)
      array(:, lbound(array, 2), :) = buffer_rcv_y

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      call MPI_Sendrecv(array(1, 1, nz + 1), Z_FACE_SIZE, MPI_REAL, high, tag, &
                        array(1, 1, 1), Z_FACE_SIZE, MPI_REAL, low, tag, &
                        comm_cart, status, ierr)

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x(:, :) = array(2, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_REAL, west, tag, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_REAL, east, tag, &
                        comm_cart, status, ierr)
      array(nx + 2, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y(:, :) = array(:, 2, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_REAL, south, tag, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_REAL, north, tag, &
                        comm_cart, status, ierr)
      array(:, ny + 2, :) = buffer_rcv_y

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      call MPI_Sendrecv(array(1, 1, 2), Z_FACE_SIZE, MPI_REAL, low, tag, &
                        array(1, 1, nz + 2), Z_FACE_SIZE, MPI_REAL, high, tag, &
                        comm_cart, status, ierr)

   end subroutine update_mpi_halo_real

   !> Perform halo exchanges in all 3 directions for an INTEGER array
   subroutine update_mpi_halo_integer(array)
      integer, intent(in out) :: array(0:nx + 1, 0:ny + 1, 0:nz + 1)

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x_int = array(nx, :, :)
      call MPI_Sendrecv(buffer_send_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, east, tag, &
                        buffer_rcv_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, west, tag, &
                        comm_cart, status, ierr)
      array(0, :, :) = buffer_rcv_x_int

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y_int = array(:, ny, :)
      call MPI_Sendrecv(buffer_send_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, north, tag, &
                        buffer_rcv_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, south, tag, &
                        comm_cart, status, ierr)
      array(:, 0, :) = buffer_rcv_y_int

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      buffer_send_z_int = array(:, :, nz)
      call MPI_Sendrecv(buffer_send_z_int(1, 1), Z_FACE_SIZE, MPI_INTEGER, high, tag, &
                        buffer_rcv_z_int(1, 1), Z_FACE_SIZE, MPI_INTEGER, low, tag, &
                        comm_cart, status, ierr)
      array(:, :, 0) = buffer_rcv_z_int

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x_int = array(1, :, :)
      call MPI_Sendrecv(buffer_send_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, west, tag, &
                        buffer_rcv_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, east, tag, &
                        comm_cart, status, ierr)
      array(nx + 1, :, :) = buffer_rcv_x_int

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y_int = array(:, 1, :)
      call MPI_Sendrecv(buffer_send_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, south, tag, &
                        buffer_rcv_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, north, tag, &
                        comm_cart, status, ierr)
      array(:, ny + 1, :) = buffer_rcv_y_int

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      buffer_send_z_int = array(:, :, 1)
      call MPI_Sendrecv(buffer_send_z_int(1, 1), Z_FACE_SIZE, MPI_INTEGER, low, tag, &
                        buffer_rcv_z_int(1, 1), Z_FACE_SIZE, MPI_INTEGER, high, tag, &
                        comm_cart, status, ierr)
      array(:, :, nz + 1) = buffer_rcv_z_int

   end subroutine update_mpi_halo_integer

end module lib_mpi_halo
