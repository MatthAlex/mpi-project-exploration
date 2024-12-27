module checks
   use lib_parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   use grid_module, only: west, east, south, north, low, high, rank => my_rank
   implicit none
   private
   public :: check_halo_real
contains

   !! Check that all* values of each halo face are equal to the Neighbouring Rank for REAL array
   !! Caveat: Corners and edges are excluded, as order of SendRecv may seed values to those points,
   !! that don't belong to the actual neighbours.
   subroutine check_halo_real(array)
      real, dimension(0:, 0:, 0:), intent(in) :: array
      if (any(int(array(0, 1:ny, 1:nz)) /= west)) &
         error stop "NOT OK in west"

      if (any(int(array(nx + 1, 1:ny, 1:nz)) /= east)) &
         error stop "NOT OK in east"

      if (any(int(array(1:nx, 0, 1:nz)) /= south)) &
         error stop "NOT OK in south"

      if (any(int(array(1:nx, ny + 1, 1:nz)) /= north)) &
         error stop "NOT OK in north"

      if (any(int(array(1:nx, 1:ny, 0)) /= low)) &
         error stop "NOT OK in low"

      if (any(int(array(1:nx, 1:ny, nz + 1)) /= high)) &
         error stop "NOT OK in high"
   end subroutine check_halo_real
end module checks
