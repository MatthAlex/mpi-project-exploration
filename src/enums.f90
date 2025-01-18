module enums
   implicit none(type, external)
   private
   !> Cartesian directions
   integer, public, parameter :: X_DIR = 0, Y_DIR = 1, Z_DIR = 2
   !> Named constants for the 6 directions in a 3D domain.
   integer, public, parameter :: D_WEST = 1, D_EAST = 2, D_SOUTH = 3, D_NORTH = 4, D_LOW = 5, D_HIGH = 6
   !> Named constants for Boundary Conditions.
   integer, public, parameter :: PERIODIC = 0, DIRICHLET = 1, NEUMANN = 2
end module enums
