!> Stores simulation parameters for a 3D MPI grid.
!> - Defines grid dimensions, iterations, and boundary conditions.
module parameters
   use precision, only: sp
   use enums, only: PERIODIC, DIRICHLET, NEUMANN
   implicit none(type, external)
   private
   public :: num_cells_x, num_cells_y, num_cells_z
   public :: iterations, boundaries, dirichlet_value, core_decomposition, dim_decomposition

   integer, parameter :: num_cells = 16
      !! Rank-level number of cells in computational domain
   integer, parameter :: num_cells_x = num_cells
      !! Rank-level number of cells in computational domain in X direction
   integer, parameter :: num_cells_y = num_cells
      !! Rank-level number of cells in computational domain in Y direction
   integer, parameter :: num_cells_z = num_cells
      !! Rank-level number of cells in computational domain in Z direction

   integer, parameter :: iterations = 2
      !! Number of time steps

   integer, parameter :: boundaries(6) = [DIRICHLET, DIRICHLET, NEUMANN, NEUMANN, PERIODIC, PERIODIC]
      !! Boundaries for each face: West, East, South, North, Low, High
   integer, parameter :: dim_decomposition = 3
      !! Number of dimensions for decomposition
   integer, parameter :: core_decomposition(dim_decomposition) = 0
      !! MPI Core decomposition. If 0, the decomposition is automatic/left to MPI

   real(kind=sp), parameter :: dirichlet_value = -1.0_sp
      !! Dirichlet Boundary Condition value
end module parameters
