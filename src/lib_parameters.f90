module lib_parameters
   implicit none
   public
   integer, parameter :: num_cells = 500
   integer, parameter :: num_cells_x = num_cells+1
   integer, parameter :: num_cells_y = num_cells+2
   integer, parameter :: num_cells_z = num_cells+3

   integer, parameter :: iterations = 1
end module lib_parameters
