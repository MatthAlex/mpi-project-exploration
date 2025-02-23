!> Handles halo region exchanges in a 3D MPI grid.
!> - Provides static buffers for efficient transfers.
!> - Supports real and integer array updates.
module mpi_halo
   use mpi_f08, only: MPI_Sendrecv, MPI_STATUS, MPI_REAL, MPI_INTEGER
   use precision, only: sp
   use mpi_domain, only: west, east, south, north, low, high, comm_cart
   use parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   implicit none(type, external)
   private
   public :: update_mpi_halo

   integer, parameter :: TAG_X = 1
   integer, parameter :: TAG_Y = 2
   integer, parameter :: TAG_Z = 3

   interface update_mpi_halo
      module procedure update_mpi_halo_real
      module procedure update_mpi_halo_integer
   end interface

contains
   !> Perform halo exchanges with nearest neighbors, in all 3 directions, for an array of kind sp
   subroutine update_mpi_halo_real(array)
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
      real(kind=sp), allocatable :: buffer_send_x(:, :), buffer_rcv_x(:, :)
      real(kind=sp), allocatable :: buffer_send_y(:, :), buffer_rcv_y(:, :)
      type(MPI_Status) :: status
      integer :: cells_x, cells_y, cells_z
         !! Size of the domain without halo region
      integer :: X_FACE_SIZE, Y_FACE_SIZE, Z_FACE_SIZE
      integer :: ierr

      cells_x = size(array, 1) - 2
      cells_y = size(array, 2) - 2
      cells_z = size(array, 3) - 2
      X_FACE_SIZE = (cells_y + 2) * (cells_z + 2)
      Y_FACE_SIZE = (cells_x + 2) * (cells_z + 2)
      Z_FACE_SIZE = (cells_x + 2) * (cells_y + 2)
      allocate (buffer_send_x(cells_x + 2, cells_z + 2), buffer_rcv_x(cells_y + 2, cells_z + 2))
      allocate (buffer_send_y(cells_y + 2, cells_z + 2), buffer_rcv_y(cells_x + 2, cells_z + 2))

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x(:, :) = array(ubound(array, 1) - 1, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_REAL, east, TAG_X, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_REAL, west, TAG_X, &
                        comm_cart, status, ierr)
      array(lbound(array, 1), :, :) = buffer_rcv_x

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y(:, :) = array(:, ubound(array, 2) - 1, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_REAL, north, TAG_Y, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_REAL, south, TAG_Y, &
                        comm_cart, status, ierr)
      array(:, lbound(array, 2), :) = buffer_rcv_y

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      call MPI_Sendrecv(array(1, 1, cells_z + 1), Z_FACE_SIZE, MPI_REAL, high, TAG_Z, &
                        array(1, 1, 1), Z_FACE_SIZE, MPI_REAL, low, TAG_Z, &
                        comm_cart, status, ierr)

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x(:, :) = array(2, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_REAL, west, TAG_X, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_REAL, east, TAG_X, &
                        comm_cart, status, ierr)
      array(cells_x + 2, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y(:, :) = array(:, 2, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_REAL, south, TAG_Y, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_REAL, north, TAG_Y, &
                        comm_cart, status, ierr)
      array(:, cells_y + 2, :) = buffer_rcv_y

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      call MPI_Sendrecv(array(1, 1, 2), Z_FACE_SIZE, MPI_REAL, low, TAG_Z, &
                        array(1, 1, cells_z + 2), Z_FACE_SIZE, MPI_REAL, high, TAG_Z, &
                        comm_cart, status, ierr)

   end subroutine update_mpi_halo_real

   !> Perform halo exchanges in all 3 directions for an INTEGER array
   subroutine update_mpi_halo_integer(array)
      integer, contiguous, intent(in out) :: array(:, :, :)
      integer, allocatable :: buffer_send_x_int(:, :), buffer_rcv_x_int(:, :)
      integer, allocatable :: buffer_send_y_int(:, :), buffer_rcv_y_int(:, :)
      type(MPI_Status) :: status
      integer :: cells_x, cells_y, cells_z
         !! Size of the domain without halo region
      integer :: X_FACE_SIZE, Y_FACE_SIZE, Z_FACE_SIZE
      integer :: ierr

      cells_x = size(array, 1) - 2
      cells_y = size(array, 2) - 2
      cells_z = size(array, 3) - 2
      X_FACE_SIZE = (cells_y + 2) * (cells_z + 2)
      Y_FACE_SIZE = (cells_x + 2) * (cells_z + 2)
      Z_FACE_SIZE = (cells_x + 2) * (cells_y + 2)
      allocate (buffer_send_x_int(cells_x + 2, cells_z + 2), buffer_rcv_x_int(cells_y + 2, cells_z + 2))
      allocate (buffer_send_y_int(cells_y + 2, cells_z + 2), buffer_rcv_y_int(cells_x + 2, cells_z + 2))

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x_int(:, :) = array(cells_x + 1, :, :)
      call MPI_Sendrecv(buffer_send_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, east, TAG_X, &
                        buffer_rcv_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, west, TAG_X, &
                        comm_cart, status, ierr)
      array(1, :, :) = buffer_rcv_x_int

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y_int(:, :) = array(:, cells_y + 1, :)
      call MPI_Sendrecv(buffer_send_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, north, TAG_Y, &
                        buffer_rcv_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, south, TAG_Y, &
                        comm_cart, status, ierr)
      array(:, 1, :) = buffer_rcv_y_int

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      call MPI_Sendrecv(array(1, 1, cells_z), Z_FACE_SIZE, MPI_INTEGER, high, TAG_Z, &
                        array(1, 1, 1), Z_FACE_SIZE, MPI_INTEGER, low, TAG_Z, &
                        comm_cart, status, ierr)

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x_int(:, :) = array(2, :, :)
      call MPI_Sendrecv(buffer_send_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, west, TAG_X, &
                        buffer_rcv_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, east, TAG_X, &
                        comm_cart, status, ierr)
      array(cells_x + 2, :, :) = buffer_rcv_x_int

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y_int(:, :) = array(:, 2, :)
      call MPI_Sendrecv(buffer_send_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, south, TAG_Y, &
                        buffer_rcv_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, north, TAG_Y, &
                        comm_cart, status, ierr)
      array(:, cells_y + 2, :) = buffer_rcv_y_int

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      call MPI_Sendrecv(array(1, 1, 1), Z_FACE_SIZE, MPI_INTEGER, low, TAG_Z, &
                        array(1, 1, cells_z + 1), Z_FACE_SIZE, MPI_INTEGER, high, TAG_Z, &
                        comm_cart, status, ierr)

   end subroutine update_mpi_halo_integer

end module mpi_halo
