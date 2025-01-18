module lib_parameters
   use precision, only: sp
   use enums, only: PERIODIC, DIRICHLET, NEUMANN
   implicit none(type, external)
   private
   public :: num_cells_x, num_cells_y, num_cells_z
   public :: iterations, boundaries, dirichlet_value, core_decomposition, dim_decomposition

   integer, parameter :: num_cells = 5
   integer, parameter :: num_cells_x = num_cells
   integer, parameter :: num_cells_y = num_cells
   integer, parameter :: num_cells_z = num_cells

   integer, parameter :: iterations = 2

   !> Boundaries for each face: West, East, South, North, Low, High
   integer, parameter :: boundaries(6) = [DIRICHLET, DIRICHLET, NEUMANN, NEUMANN, PERIODIC, PERIODIC]
   !> Number of dimensions for decomposition
   integer, parameter :: dim_decomposition = 3
   !> MPI Core decomposition. If 0, the decomposition is automatic/left to MPI
   integer, parameter :: core_decomposition(dim_decomposition) = 0

   !> Dirichlet Boundary Condition value
   real(kind=sp), parameter :: dirichlet_value = -1.0_sp
end module lib_parameters
