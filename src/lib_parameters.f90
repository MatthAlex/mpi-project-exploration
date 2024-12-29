module lib_parameters
   implicit none
   public
   integer, parameter :: num_cells = 5
   integer, parameter :: num_cells_x = num_cells+1
   integer, parameter :: num_cells_y = num_cells+2
   integer, parameter :: num_cells_z = num_cells+3

   integer, parameter :: iterations = 2

   !> 0 is Periodic, 1 is Dirichlet, 2 is Neumann
   integer, parameter :: boundaries(6) = [1, 1, 2, 2, 0, 0]
end module lib_parameters
