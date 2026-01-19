!> Handles halo region exchanges in a 3D MPI grid.
!> - Provides static buffers for efficient transfers.
!> - Supports real and integer array updates.
module mpi_halo
   use mpi_domain_types, only: mpi_domain_t
   use mpi_f08, only: MPI_Sendrecv, MPI_STATUS, MPI_REAL, MPI_INTEGER, MPI_Comm
   use lib_mpi_precision, only: sp
   use lib_mpi_enums, only: D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH
   implicit none(type, external)
   private
   public :: update_mpi_halo

   integer, parameter :: TAG_X = 1, TAG_Y = 2, TAG_Z = 3
      !! MPI tags for X, Y, Z directions

   interface update_mpi_halo
      !! Generic interface that wraps the real(sp) and integer routines
      module procedure update_mpi_halo_real_facade
      module procedure update_mpi_halo_integer_facade
   end interface update_mpi_halo

contains

   subroutine update_mpi_halo_real_facade(domain, array)
      class(mpi_domain_t), intent(in) :: domain
      real(sp), contiguous, intent(inout) :: array(:, :, :)

      call update_mpi_halo_real(domain, array, size(array, 1), size(array, 2), size(array, 3))
   end subroutine update_mpi_halo_real_facade

   subroutine update_mpi_halo_integer_facade(domain, array)
      class(mpi_domain_t), intent(in) :: domain
      integer, contiguous, intent(inout) :: array(:, :, :)

      call update_mpi_halo_integer(domain, array, size(array, 1), size(array, 2), size(array, 3))
   end subroutine update_mpi_halo_integer_facade

   !> Perform halo exchanges with nearest neighbors, in all 3 directions, for an array of kind sp
   subroutine update_mpi_halo_real(domain, array, nx, ny, nz)
      class(mpi_domain_t), intent(in) :: domain
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
      integer, intent(in) :: nx, ny, nz
      real(kind=sp) :: buffer_send_x(ny, nz), buffer_rcv_x(ny, nz)
      real(kind=sp) :: buffer_send_y(nx, nz), buffer_rcv_y(nx, nz)
      type(MPI_Status) :: status
      type(MPI_Comm) :: comm_cart
      integer :: X_FACE_SIZE, Y_FACE_SIZE, Z_FACE_SIZE
      integer :: neighbors(6)
      integer :: ierr

      X_FACE_SIZE = ny * nz
      Y_FACE_SIZE = nx * nz
      Z_FACE_SIZE = nx * ny

      buffer_rcv_x = 0.0_sp
      buffer_rcv_y = 0.0_sp

      neighbors = domain%get_neighbors()
      comm_cart = domain%get_communicator()

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x(:, :) = array(nx - 1, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_REAL, neighbors(D_EAST), TAG_X, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_REAL, neighbors(D_WEST), TAG_X, &
                        comm_cart, status, ierr)
      array(1, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y(:, :) = array(:, ny - 1, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_REAL, neighbors(D_NORTH), TAG_Y, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_REAL, neighbors(D_SOUTH), TAG_Y, &
                        comm_cart, status, ierr)
      array(:, 1, :) = buffer_rcv_y

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      call MPI_Sendrecv(array(1, 1, nz - 1), Z_FACE_SIZE, MPI_REAL, neighbors(D_HIGH), TAG_Z, &
                        array(1, 1, 1), Z_FACE_SIZE, MPI_REAL, neighbors(D_LOW), TAG_Z, &
                        comm_cart, status, ierr)

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x(:, :) = array(2, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_REAL, neighbors(D_WEST), TAG_X, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_REAL, neighbors(D_EAST), TAG_X, &
                        comm_cart, status, ierr)
      array(nx, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y(:, :) = array(:, 2, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_REAL, neighbors(D_SOUTH), TAG_Y, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_REAL, neighbors(D_NORTH), TAG_Y, &
                        comm_cart, status, ierr)
      array(:, ny, :) = buffer_rcv_y

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      call MPI_Sendrecv(array(1, 1, 2), Z_FACE_SIZE, MPI_REAL, neighbors(D_LOW), TAG_Z, &
                        array(1, 1, nz), Z_FACE_SIZE, MPI_REAL, neighbors(D_HIGH), TAG_Z, &
                        comm_cart, status, ierr)

   end subroutine update_mpi_halo_real

   !> Perform halo exchanges in all 3 directions for an INTEGER array
   subroutine update_mpi_halo_integer(domain, array, nx, ny, nz)
      class(mpi_domain_t), intent(in) :: domain
      integer, contiguous, intent(in out) :: array(:, :, :)
      integer, intent(in) :: nx, ny, nz
      integer :: buffer_send_x(ny, nz), buffer_rcv_x(ny, nz)
      integer :: buffer_send_y(nx, nz), buffer_rcv_y(nx, nz)
      type(MPI_Status) :: status
      type(MPI_Comm) :: comm_cart
      integer :: X_FACE_SIZE, Y_FACE_SIZE, Z_FACE_SIZE
      integer :: neighbors(6)
      integer :: ierr

      X_FACE_SIZE = ny * nz
      Y_FACE_SIZE = nx * nz
      Z_FACE_SIZE = nx * ny

      buffer_rcv_x = 0
      buffer_rcv_y = 0

      neighbors = domain%get_neighbors()
      comm_cart = domain%get_communicator()

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x(:, :) = array(nx - 1, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_INTEGER, neighbors(D_EAST), TAG_X, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_INTEGER, neighbors(D_WEST), TAG_X, &
                        comm_cart, status, ierr)
      array(1, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y(:, :) = array(:, ny - 1, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_INTEGER, neighbors(D_NORTH), TAG_Y, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_INTEGER, neighbors(D_SOUTH), TAG_Y, &
                        comm_cart, status, ierr)
      array(:, 1, :) = buffer_rcv_y

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      call MPI_Sendrecv(array(1, 1, nz - 1), Z_FACE_SIZE, MPI_INTEGER, neighbors(D_HIGH), TAG_Z, &
                        array(1, 1, 1), Z_FACE_SIZE, MPI_INTEGER, neighbors(D_LOW), TAG_Z, &
                        comm_cart, status, ierr)

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x(:, :) = array(2, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), X_FACE_SIZE, MPI_INTEGER, neighbors(D_WEST), TAG_X, &
                        buffer_rcv_x(1, 1), X_FACE_SIZE, MPI_INTEGER, neighbors(D_EAST), TAG_X, &
                        comm_cart, status, ierr)
      array(nx, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y(:, :) = array(:, 2, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), Y_FACE_SIZE, MPI_INTEGER, neighbors(D_SOUTH), TAG_Y, &
                        buffer_rcv_y(1, 1), Y_FACE_SIZE, MPI_INTEGER, neighbors(D_NORTH), TAG_Y, &
                        comm_cart, status, ierr)
      array(:, ny, :) = buffer_rcv_y

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      call MPI_Sendrecv(array(1, 1, 2), Z_FACE_SIZE, MPI_INTEGER, neighbors(D_LOW), TAG_Z, &
                        array(1, 1, nz), Z_FACE_SIZE, MPI_INTEGER, neighbors(D_HIGH), TAG_Z, &
                        comm_cart, status, ierr)

   end subroutine update_mpi_halo_integer

end module mpi_halo
