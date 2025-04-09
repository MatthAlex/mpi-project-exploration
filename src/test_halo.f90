module test_halo
   use mpi_domain_types, only: mpi_domain_t
   use parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   use precision, only: sp
   use enums, only: D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH
   implicit none(type, external)
   private
   public :: check_halo_real
contains

   !! Check that all* values of each halo face are equal to the Neighbouring Rank for REAL array
   !! Caveat: Corners and edges are excluded, as order of SendRecv may seed values to those points,
   !! that don't belong to the actual neighbours.
   !! Reminder: we only test the non-physical halo regions
   subroutine check_halo_real(domain, array)
      class(mpi_domain_t), intent(in) :: domain
      real(kind=sp), dimension(:, :, :), intent(in) :: array
      integer :: neighbors(6)
      logical :: is_bc_face(6)
      integer :: face
      neighbors = domain%get_neighbors()
      is_bc_face = domain%is_boundary_face(:)

      do face = 1, 6
         if (is_bc_face(face)) cycle

         select case (face)
         case(D_WEST)
            if (any(int(array(1, 2:size(array, 2)-1, 2:size(array, 3)-1)) /= neighbors(D_WEST))) &
               call domain%abort("TEST HALO: Not OK in west")
         case(D_EAST)
            if (any(int(array(size(array, 1), 2:size(array, 2)-1, 2:size(array, 3)-1)) /= neighbors(D_EAST))) &
               call domain%abort("TEST HALO: Not OK in east")
         case(D_SOUTH)
            if (any(int(array(2:size(array, 1)-1, 1, 2:size(array, 3)-1)) /= neighbors(D_SOUTH))) &
               call domain%abort("TEST HALO: Not OK in south")
         case(D_NORTH)
            if (any(int(array(2:size(array, 1)-1, size(array, 2), 2:size(array, 3)-1)) /= neighbors(D_NORTH))) &
               call domain%abort("TEST HALO: Not OK in north")
         case(D_LOW)
            if (any(int(array(2:size(array, 1)-1, 2:size(array, 2)-1, 1)) /= neighbors(D_LOW))) &
               call domain%abort("TEST HALO: Not OK in low")
         case(D_HIGH)
            if (any(int(array(2:size(array, 1)-1, 2:size(array, 2)-1, size(array, 3))) /= neighbors(D_HIGH))) &
               call domain%abort("TEST HALO: Not OK in high")
         end select
      end do
   end subroutine check_halo_real
end module test_halo
