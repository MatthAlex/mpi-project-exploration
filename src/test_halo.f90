module test_halo
   use lib_parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   use grid_module, only: west, east, south, north, low, high, rank => my_rank
   implicit none
   private
   public :: check_halo_real
contains

   !! Check that all* values of each halo face are equal to the Neighbouring Rank for REAL array
   !! Caveat: Corners and edges are excluded, as order of SendRecv may seed values to those points,
   !! that don't belong to the actual neighbours.
   !! Reminder: we only test the non-physical halo regions
   subroutine check_halo_real(array)
      use boundary, only: is_bc_face
      real, dimension(0:, 0:, 0:), intent(in) :: array
      if (any(int(array(0, 1:ny, 1:nz)) /= west .and. .not. is_bc_face(1))) &
         error stop "NOT OK in west"

      if (any(int(array(nx + 1, 1:ny, 1:nz)) /= east .and. .not. is_bc_face(2))) &
         error stop "NOT OK in east"

      if (any(int(array(1:nx, 0, 1:nz)) /= south .and. .not. is_bc_face(3))) &
         error stop "NOT OK in south"

      if (any(int(array(1:nx, ny + 1, 1:nz)) /= north .and. .not. is_bc_face(4))) &
         error stop "NOT OK in north"

      if (any(int(array(1:nx, 1:ny, 0)) /= low .and. .not. is_bc_face(5))) &
         error stop "NOT OK in low"

      if (any(int(array(1:nx, 1:ny, nz + 1)) /= high .and. .not. is_bc_face(6))) &
         error stop "NOT OK in high"
   end subroutine check_halo_real
end module test_halo
