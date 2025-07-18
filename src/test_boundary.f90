!> Helper module that defines how to assert validity for a boundary update
module test_boundary
   use mpi_domain_types, only: mpi_domain_t
   use lib_mpi_parameters, only: nx => num_cells_x, ny => num_cells_y, nz => num_cells_z, boundaries
   use lib_mpi_precision, only: sp
   use lib_mpi_enums, only: D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH, PERIODIC, DIRICHLET, NEUMANN
   implicit none(type, external)
   private
   public :: check_boundary_real

contains

   !> Check that boundary faces have correct values based on their boundary condition type
   !> - Periodic (0): Handled by test_halo
   !> - Dirichlet (1): Should be constant value 1.0
   !> - Neumann (2): Should match adjacent interior cell
   subroutine check_boundary_real(domain, array)
      use lib_mpi_parameters, only: dirichlet_value
      class(mpi_domain_t), intent(in) :: domain
      real(kind=sp), dimension(:, :, :), intent(in) :: array
      real(kind=sp), parameter :: tolerance = 1.0e-10_sp
      logical :: is_bc_face(6)

      is_bc_face = domain%is_boundary_face(:)

      ! West face
      if (is_bc_face(D_WEST)) then
         select case (boundaries(D_WEST))
         case (DIRICHLET)
            if (any(abs(array(1, 2:size(array, 2) - 1, 2:size(array, 3) - 1) - dirichlet_value) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: West face Dirichlet")
            end if
         case (NEUMANN)
            if (any(abs(array(1, 2:size(array, 2) - 1, 2:size(array, 3) - 1) - &
                        array(2, 2:size(array, 2) - 1, 2:size(array, 3) - 1)) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: West face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary: Boundary Condition index not recognized... Exiting..")
         end select
      end if

      ! East face
      if (is_bc_face(D_EAST)) then
         select case (boundaries(D_EAST))
         case (DIRICHLET)
            if (any(abs(array(ubound(array, 1), 2:size(array, 2) - 1, 2:size(array, 3) - 1) - dirichlet_value) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: East face Dirichlet")
            end if
         case (NEUMANN)
            if (any(abs(array(ubound(array, 1), 2:size(array, 2) - 1, 2:size(array, 3) - 1) - &
                        array(ubound(array, 1) - 1, 2:size(array, 2) - 1, 2:size(array, 3) - 1)) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: East face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary: Boundary Condition index not recognized... Exiting..")
         end select
      end if

      ! South face
      if (is_bc_face(D_SOUTH)) then
         select case (boundaries(D_SOUTH))
         case (DIRICHLET)
            if (any(abs(array(2:size(array, 1) - 1, 1, 2:size(array, 3) - 1) - dirichlet_value) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: South face Dirichlet")
            end if
         case (NEUMANN)
            if (any(abs(array(2:size(array, 1) - 1, 1, 2:size(array, 3) - 1) - &
                        array(2:size(array, 1) - 1, 2, 2:size(array, 3) - 1)) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: South face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary: Boundary Condition index not recognized... Exiting..")
         end select
      end if

      ! North face
      if (is_bc_face(D_NORTH)) then
         select case (boundaries(D_NORTH))
         case (DIRICHLET)
            if (any(abs(array(2:size(array, 1) - 1, ubound(array, 2), 2:size(array, 3) - 1) - dirichlet_value) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: North face Dirichlet")
            end if
         case (NEUMANN)
            if (any(abs(array(2:size(array, 1) - 1, ubound(array, 2), 2:size(array, 3) - 1) - &
                        array(2:size(array, 1) - 1, ubound(array, 2) - 1, 2:size(array, 3) - 1)) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: North face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary: Boundary Condition index not recognized... Exiting..")
         end select
      end if

      ! Low face
      if (is_bc_face(D_LOW)) then
         select case (boundaries(D_LOW))
         case (DIRICHLET)
            if (any(abs(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, 1) - dirichlet_value) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: Low face Dirichlet")
            end if
         case (NEUMANN)
            if (any(abs(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, 1) - &
                        array(2:size(array, 1) - 1, 2:size(array, 2) - 1, 2)) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: Low face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary: Boundary Condition index not recognized... Exiting..")
         end select
      end if

      ! High face
      if (is_bc_face(D_HIGH)) then
         select case (boundaries(D_HIGH))
         case (DIRICHLET)
            if (any(abs(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, ubound(array, 3)) - dirichlet_value) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: High face Dirichlet")
            end if
         case (NEUMANN)
            if (any(abs(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, ubound(array, 3)) - &
                        array(2:size(array, 1) - 1, 2:size(array, 2) - 1, ubound(array, 3) - 1)) > tolerance)) then
               call domain%abort("TEST Boundary: Not OK: High face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary: Boundary Condition index not recognized... Exiting..")
         end select
      end if

   end subroutine check_boundary_real

end module test_boundary
