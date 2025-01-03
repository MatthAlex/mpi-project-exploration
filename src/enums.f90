module enums
   implicit none
   public
   !> Named constants for the 6 directions in a 3D domain.
   integer, parameter :: D_WEST = 1, D_EAST = 2, D_SOUTH = 3, D_NORTH = 4, D_LOW = 5, D_HIGH = 6
   !> Named constants for Boundary Conditions.
   integer, parameter :: PERIODIC = 0, DIRICHLET = 1, NEUMANN = 2

end module enums
