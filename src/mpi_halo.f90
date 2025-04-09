!> Handles halo region exchanges in a 3D MPI grid.
!> - Provides static buffers for efficient transfers.
!> - Supports real and integer array updates.
module mpi_halo
   use mpi_domain_types, only: mpi_domain_t
   use mpi_f08, only: MPI_Sendrecv, MPI_STATUS, MPI_REAL, MPI_INTEGER, MPI_Comm
   use precision, only: sp
   use parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   use enums, only: D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH
   implicit none(type, external)
   private
   public :: update_mpi_halo

   integer, parameter :: TAG_X = 1, TAG_Y = 2, TAG_Z = 3
      !! MPI tags for X, Y, Z directions

   interface update_mpi_halo
      !! Generic interface that wraps the real(sp) and integer routines
      module procedure update_mpi_halo_real
      module procedure update_mpi_halo_integer
   end interface

contains
   !> Perform halo exchanges with nearest neighbors, in all 3 directions, for an array of kind sp
   subroutine update_mpi_halo_real(domain, array)
      class(mpi_domain_t), intent(in) :: domain
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
      real(kind=sp), allocatable :: buffer_send_x(:, :), buffer_rcv_x(:, :)
      real(kind=sp), allocatable :: buffer_send_y(:, :), buffer_rcv_y(:, :)
      type(MPI_Status) :: status
      type(MPI_Comm) :: comm_cart
      integer :: cells_x, cells_y, cells_z
         !! Size of the domain without halo region
      integer :: X_FACE_SIZE, Y_FACE_SIZE, Z_FACE_SIZE
      integer :: neighbors(6)
      integer :: ierr

      cells_x = size(array, 1) - 2
      cells_y = size(array, 2) - 2
      cells_z = size(array, 3) - 2
      X_FACE_SIZE = (cells_y + 2) * (cells_z + 2)
      Y_FACE_SIZE = (cells_x + 2) * (cells_z + 2)
      Z_FACE_SIZE = (cells_x + 2) * (cells_y + 2)
      allocate (buffer_send_x(cells_y + 2, cells_z + 2), buffer_rcv_x(cells_y + 2, cells_z + 2))
      allocate (buffer_send_y(cells_x + 2, cells_z + 2), buffer_rcv_y(cells_x + 2, cells_z + 2))

      neighbors = domain%get_neighbors()
      comm_cart = domain%get_communicator()

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x(:, :) = array(ubound(array, 1) - 1, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_REAL, neighbors(D_EAST), TAG_X, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_REAL, neighbors(D_WEST), TAG_X, &
                        comm_cart, status, ierr)
      array(lbound(array, 1), :, :) = buffer_rcv_x

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y(:, :) = array(:, ubound(array, 2) - 1, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_REAL, neighbors(D_NORTH), TAG_Y, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_REAL, neighbors(D_SOUTH), TAG_Y, &
                        comm_cart, status, ierr)
      array(:, lbound(array, 2), :) = buffer_rcv_y

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      call MPI_Sendrecv(array(1, 1, cells_z + 1), Z_FACE_SIZE, MPI_REAL, neighbors(D_HIGH), TAG_Z, &
                        array(1, 1, 1), Z_FACE_SIZE, MPI_REAL, neighbors(D_LOW), TAG_Z, &
                        comm_cart, status, ierr)

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x(:, :) = array(2, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_REAL, neighbors(D_WEST), TAG_X, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_REAL, neighbors(D_EAST), TAG_X, &
                        comm_cart, status, ierr)
      array(cells_x + 2, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y(:, :) = array(:, 2, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_REAL, neighbors(D_SOUTH), TAG_Y, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_REAL, neighbors(D_NORTH), TAG_Y, &
                        comm_cart, status, ierr)
      array(:, cells_y + 2, :) = buffer_rcv_y

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      call MPI_Sendrecv(array(1, 1, 2), Z_FACE_SIZE, MPI_REAL, neighbors(D_LOW), TAG_Z, &
                        array(1, 1, cells_z + 2), Z_FACE_SIZE, MPI_REAL, neighbors(D_HIGH), TAG_Z, &
                        comm_cart, status, ierr)

   end subroutine update_mpi_halo_real

   !> Perform halo exchanges in all 3 directions for an INTEGER array
   subroutine update_mpi_halo_integer(domain, array)
      class(mpi_domain_t), intent(in) :: domain
      integer, contiguous, intent(in out) :: array(:, :, :)
      integer, allocatable :: buffer_send_x_int(:, :), buffer_rcv_x_int(:, :)
      integer, allocatable :: buffer_send_y_int(:, :), buffer_rcv_y_int(:, :)
      type(MPI_Status) :: status
      type(MPI_Comm) :: comm_cart
      integer :: cells_x, cells_y, cells_z
         !! Size of the domain without halo region
      integer :: X_FACE_SIZE, Y_FACE_SIZE, Z_FACE_SIZE
      integer :: neighbors(6)
      integer :: ierr

      cells_x = size(array, 1) - 2
      cells_y = size(array, 2) - 2
      cells_z = size(array, 3) - 2
      X_FACE_SIZE = (cells_y + 2) * (cells_z + 2)
      Y_FACE_SIZE = (cells_x + 2) * (cells_z + 2)
      Z_FACE_SIZE = (cells_x + 2) * (cells_y + 2)
      allocate (buffer_send_x_int(cells_y + 2, cells_z + 2), buffer_rcv_x_int(cells_y + 2, cells_z + 2))
      allocate (buffer_send_y_int(cells_x + 2, cells_z + 2), buffer_rcv_y_int(cells_x + 2, cells_z + 2))

      neighbors = domain%get_neighbors()
      comm_cart = domain%get_communicator()

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x_int(:, :) = array(cells_x + 1, :, :)
      call MPI_Sendrecv(buffer_send_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, neighbors(D_EAST), TAG_X, &
                        buffer_rcv_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, neighbors(D_WEST), TAG_X, &
                        comm_cart, status, ierr)
      array(1, :, :) = buffer_rcv_x_int

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y_int(:, :) = array(:, cells_y + 1, :)
      call MPI_Sendrecv(buffer_send_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, neighbors(D_NORTH), TAG_Y, &
                        buffer_rcv_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, neighbors(D_SOUTH), TAG_Y, &
                        comm_cart, status, ierr)
      array(:, 1, :) = buffer_rcv_y_int

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      call MPI_Sendrecv(array(1, 1, cells_z), Z_FACE_SIZE, MPI_INTEGER, neighbors(D_HIGH), TAG_Z, &
                        array(1, 1, 1), Z_FACE_SIZE, MPI_INTEGER, neighbors(D_LOW), TAG_Z, &
                        comm_cart, status, ierr)

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x_int(:, :) = array(2, :, :)
      call MPI_Sendrecv(buffer_send_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, neighbors(D_WEST), TAG_X, &
                        buffer_rcv_x_int(1, 1), X_FACE_SIZE, MPI_INTEGER, neighbors(D_EAST), TAG_X, &
                        comm_cart, status, ierr)
      array(cells_x + 2, :, :) = buffer_rcv_x_int

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y_int(:, :) = array(:, 2, :)
      call MPI_Sendrecv(buffer_send_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, neighbors(D_SOUTH), TAG_Y, &
                        buffer_rcv_y_int(1, 1), Y_FACE_SIZE, MPI_INTEGER, neighbors(D_NORTH), TAG_Y, &
                        comm_cart, status, ierr)
      array(:, cells_y + 2, :) = buffer_rcv_y_int

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      call MPI_Sendrecv(array(1, 1, 1), Z_FACE_SIZE, MPI_INTEGER, neighbors(D_LOW), TAG_Z, &
                        array(1, 1, cells_z + 1), Z_FACE_SIZE, MPI_INTEGER, neighbors(D_HIGH), TAG_Z, &
                        comm_cart, status, ierr)

   end subroutine update_mpi_halo_integer

end module mpi_halo
