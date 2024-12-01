!> Encapsulates MPI halo update functionality
!> - static transfer buffers - avoid allocation every update
!> - module procedure to generalise subroutine calling
!> TODO: assumes cube of num_cells size
module lib_mpi_halo
   use mpi
   use grid_module, only: west, east, south, north, low, high, comm_cart
   use lib_parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   implicit none
   public :: update_mpi_halo

   integer :: ierr, rank, comm
   integer :: status(MPI_STATUS_SIZE)

   !> Buffers for integer transfers - avoids allocation
   integer, dimension(ny + 2, nz + 2) :: buffer_send_x_int
   integer, dimension(ny + 2, nz + 2) :: buffer_rcv_x_int
   integer, dimension(nx + 2, nz + 2) :: buffer_send_y_int
   integer, dimension(nx + 2, nz + 2) :: buffer_rcv_y_int
   integer, dimension(nx + 2, ny + 2) :: buffer_send_z_int
   integer, dimension(nx + 2, ny + 2) :: buffer_rcv_z_int
   !> Buffers for real transfers - avoids allocation
   real, dimension(ny + 2, nz + 2) :: buffer_send_x
   real, dimension(ny + 2, nz + 2) :: buffer_rcv_x
   real, dimension(nx + 2, nz + 2) :: buffer_send_y
   real, dimension(nx + 2, nz + 2) :: buffer_rcv_y
   real, dimension(nx + 2, ny + 2) :: buffer_send_z
   real, dimension(nx + 2, ny + 2) :: buffer_rcv_z
   integer, parameter :: tag = 1

   interface update_mpi_halo
      module procedure update_mpi_halo_real
      module procedure update_mpi_halo_integer
   end interface

contains

   !> Perform halo exchanges in all 3 directions for a REAL array
   subroutine update_mpi_halo_real(array)
      real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(in out) :: array

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x = array(nx, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), (ny + 2) * (nz + 2), MPI_REAL, east, tag, &
                        buffer_rcv_x(1, 1), (ny + 2) * (nz + 2), MPI_REAL, west, tag, &
                        comm_cart, status, ierr)
      array(0, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y = array(:, ny, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), (nx + 2) * (nz + 2), MPI_REAL, north, tag, &
                        buffer_rcv_y(1, 1), (nx + 2) * (nz + 2), MPI_REAL, south, tag, &
                        comm_cart, status, ierr)
      array(:, 0, :) = buffer_rcv_y

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      buffer_send_z = array(:, :, nz)
      call MPI_Sendrecv(buffer_send_z(1, 1), (nx + 2) * (ny + 2), MPI_REAL, high, tag, &
                        buffer_rcv_z(1, 1), (nx + 2) * (ny + 2), MPI_REAL, low, tag, &
                        comm_cart, status, ierr)
      array(:, :, 0) = buffer_rcv_z

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x = array(1, :, :)
      call MPI_Sendrecv(buffer_send_x(1, 1), (ny + 2) * (nz + 2), MPI_REAL, west, tag, &
                        buffer_rcv_x(1, 1), (ny + 2) * (nz + 2), MPI_REAL, east, tag, &
                        comm_cart, status, ierr)
      array(nx + 1, :, :) = buffer_rcv_x

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y = array(:, 1, :)
      call MPI_Sendrecv(buffer_send_y(1, 1), (nx + 2) * (nz + 2), MPI_REAL, south, tag, &
                        buffer_rcv_y(1, 1), (nx + 2) * (nz + 2), MPI_REAL, north, tag, &
                        comm_cart, status, ierr)
      array(:, ny + 1, :) = buffer_rcv_y

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      buffer_send_z = array(:, :, 1)
      call MPI_Sendrecv(buffer_send_z(1, 1), (nx + 2) * (ny + 2), MPI_REAL, low, tag, &
                        buffer_rcv_z(1, 1), (nx + 2) * (ny + 2), MPI_REAL, high, tag, &
                        comm_cart, status, ierr)
      array(:, :, nz + 1) = buffer_rcv_z

   end subroutine update_mpi_halo_real

   !> Perform halo exchanges in all 3 directions for an INTEGER array
   subroutine update_mpi_halo_integer(array)
      integer, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(in out) :: array

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_x_int = array(nx, :, :)
      call MPI_Sendrecv(buffer_send_x_int(1, 1), (ny + 2) * (nz + 2), MPI_INTEGER, east, tag, &
                        buffer_rcv_x_int(1, 1), (ny + 2) * (nz + 2), MPI_INTEGER, west, tag, &
                        comm_cart, status, ierr)
      array(0, :, :) = buffer_rcv_x_int

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_y_int = array(:, ny, :)
      call MPI_Sendrecv(buffer_send_y_int(1, 1), (nx + 2) * (nz + 2), MPI_INTEGER, north, tag, &
                        buffer_rcv_y_int(1, 1), (nx + 2) * (nz + 2), MPI_INTEGER, south, tag, &
                        comm_cart, status, ierr)
      array(:, 0, :) = buffer_rcv_y_int

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      buffer_send_z_int = array(:, :, nz)
      call MPI_Sendrecv(buffer_send_z_int(1, 1), (nx + 2) * (ny + 2), MPI_INTEGER, high, tag, &
                        buffer_rcv_z_int(1, 1), (nx + 2) * (ny + 2), MPI_INTEGER, low, tag, &
                        comm_cart, status, ierr)
      array(:, :, 0) = buffer_rcv_z_int

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_x_int = array(1, :, :)
      call MPI_Sendrecv(buffer_send_x_int(1, 1), (ny + 2) * (nz + 2), MPI_INTEGER, west, tag, &
                        buffer_rcv_x_int(1, 1), (ny + 2) * (nz + 2), MPI_INTEGER, east, tag, &
                        comm_cart, status, ierr)
      array(nx + 1, :, :) = buffer_rcv_x_int

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_y_int = array(:, 1, :)
      call MPI_Sendrecv(buffer_send_y_int(1, 1), (nx + 2) * (nz + 2), MPI_INTEGER, south, tag, &
                        buffer_rcv_y_int(1, 1), (nx + 2) * (nz + 2), MPI_INTEGER, north, tag, &
                        comm_cart, status, ierr)
      array(:, ny + 1, :) = buffer_rcv_y_int

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      buffer_send_z_int = array(:, :, 1)
      call MPI_Sendrecv(buffer_send_z_int(1, 1), (nx + 2) * (ny + 2), MPI_INTEGER, low, tag, &
                        buffer_rcv_z_int(1, 1), (nx + 2) * (ny + 2), MPI_INTEGER, high, tag, &
                        comm_cart, status, ierr)
      array(:, :, nz + 1) = buffer_rcv_z_int

   end subroutine update_mpi_halo_integer

end module lib_mpi_halo
