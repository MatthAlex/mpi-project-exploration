!> Encapsulates MPI halo update functionality
!> - static transfer buffers - avoid allocation every update
!> - module procedure to generalise subroutine calling
!> TODO: assumes cube of num_cells size
module lib_mpi_halo
   use mpi
   use grid_module, only: west, east, south, north, low, high, comm_cart
   use lib_parameters, only: num_cells
   implicit none
   public :: update_mpi_halo

   integer :: ierr, rank, comm
   integer :: status(MPI_STATUS_SIZE)

   !> Buffers for integer transfers - avoids allocation
   integer, dimension(num_cells + 2, num_cells + 2) :: buffer_send_int
   integer, dimension(num_cells + 2, num_cells + 2) :: buffer_rcv_int
   !> Buffers for real transfers - avoids allocation
   real, dimension(num_cells + 2, num_cells + 2) :: buffer_send
   real, dimension(num_cells + 2, num_cells + 2) :: buffer_rcv
   integer, parameter :: tag = 1

   interface update_mpi_halo
      module procedure update_mpi_halo_real
      module procedure update_mpi_halo_integer
   end interface

contains

   !> Perform halo exchanges in all 3 directions for a REAL array
   subroutine update_mpi_halo_real(array)
      real, dimension(0:num_cells + 1, 0:num_cells + 1, 0:num_cells + 1), intent(in out) :: array

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send = array(num_cells, :, :)
      call MPI_Sendrecv(buffer_send(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, east, tag, &
                        buffer_rcv(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, west, tag, &
                        comm_cart, status, ierr)
      array(0, :, :) = buffer_rcv

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send = array(:, num_cells, :)
      call MPI_Sendrecv(buffer_send(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, north, tag, &
                        buffer_rcv(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, south, tag, &
                        comm_cart, status, ierr)
      array(:, 0, :) = buffer_rcv

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      buffer_send = array(:, :, num_cells)
      call MPI_Sendrecv(buffer_send(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, high, tag, &
                        buffer_rcv(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, low, tag, &
                        comm_cart, status, ierr)
      array(:, :, 0) = buffer_rcv

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send = array(1, :, :)
      call MPI_Sendrecv(buffer_send(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, west, tag, &
                        buffer_rcv(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, east, tag, &
                        comm_cart, status, ierr)
      array(num_cells + 1, :, :) = buffer_rcv

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send = array(:, 1, :)
      call MPI_Sendrecv(buffer_send(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, south, tag, &
                        buffer_rcv(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, north, tag, &
                        comm_cart, status, ierr)
      array(:, num_cells + 1, :) = buffer_rcv

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      buffer_send = array(:, :, 1)
      call MPI_Sendrecv(buffer_send(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, low, tag, &
                        buffer_rcv(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, high, tag, &
                        comm_cart, status, ierr)
      array(:, :, num_cells + 1) = buffer_rcv

   end subroutine update_mpi_halo_real

   !> Perform halo exchanges in all 3 directions for an INTEGER array
   subroutine update_mpi_halo_integer(array)
      integer, dimension(0:num_cells + 1, 0:num_cells + 1, 0:num_cells + 1), intent(in out) :: array

      ! Exchange X direction; I send East halo to East neighbour, I receive West halo from West Neighbour
      buffer_send_int = array(num_cells, :, :)
      call MPI_Sendrecv(buffer_send_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, east, tag, &
                        buffer_rcv_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, west, tag, &
                        comm_cart, status, ierr)
      array(0, :, :) = buffer_rcv_int

      ! Exchange Y direction; I send North halo to North neighbour, I receive South halo from South Neighbour
      buffer_send_int = array(:, num_cells, :)
      call MPI_Sendrecv(buffer_send_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, north, tag, &
                        buffer_rcv_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, south, tag, &
                        comm_cart, status, ierr)
      array(:, 0, :) = buffer_rcv_int

      ! Exchange Z direction; I send High halo to High neighbour, I receive Low halo from Low Neighbour
      buffer_send_int = array(:, :, num_cells)
      call MPI_Sendrecv(buffer_send_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, high, tag, &
                        buffer_rcv_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, low, tag, &
                        comm_cart, status, ierr)
      array(:, :, 0) = buffer_rcv_int

      ! Exchange X direction; I send West halo to West neighbour, I receive East halo from East Neighbour
      buffer_send_int = array(1, :, :)
      call MPI_Sendrecv(buffer_send_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, west, tag, &
                        buffer_rcv_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, east, tag, &
                        comm_cart, status, ierr)
      array(num_cells + 1, :, :) = buffer_rcv_int

      ! Exchange Y direction; I send South halo to South neighbour, I receive North halo from North Neighbour
      buffer_send_int = array(:, 1, :)
      call MPI_Sendrecv(buffer_send_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, south, tag, &
                        buffer_rcv_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, north, tag, &
                        comm_cart, status, ierr)
      array(:, num_cells + 1, :) = buffer_rcv_int

      ! Exchange Z direction; I send Low halo to Low neighbour, I receive High halo from High Neighbour
      buffer_send_int = array(:, :, 1)
      call MPI_Sendrecv(buffer_send_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, low, tag, &
                        buffer_rcv_int(1, 1), (num_cells + 2) * (num_cells + 2), MPI_REAL, high, tag, &
                        comm_cart, status, ierr)
      array(:, :, num_cells + 1) = buffer_rcv_int

   end subroutine update_mpi_halo_integer

end module lib_mpi_halo
