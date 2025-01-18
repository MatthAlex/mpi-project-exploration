module lib_parameters
   use precision, only: sp
   use enums, only: PERIODIC, DIRICHLET, NEUMANN
   implicit none(type, external)
   private
   public :: num_cells_x, num_cells_y, num_cells_z
   public :: iterations, boundaries, dirichlet_value

   integer, parameter :: num_cells = 5
   integer, parameter :: num_cells_x = num_cells
   integer, parameter :: num_cells_y = num_cells
   integer, parameter :: num_cells_z = num_cells

   integer, parameter :: iterations = 2

   !> Boundaries for each face: West, East, South, North, Low, High
   integer, parameter :: boundaries(6) = [DIRICHLET, DIRICHLET, NEUMANN, NEUMANN, PERIODIC, PERIODIC]

   !> Dirichlet Boundary Condition value
   real(kind=sp), parameter :: dirichlet_value = -1.0_sp
end module lib_parameters
