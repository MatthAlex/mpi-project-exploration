!> Helper module that defines how to assert validity for a boundary update
module test_boundary
   use mpi_domain_types, only: mpi_domain_t
   use lib_mpi_precision, only: sp
   use lib_mpi_enums, only: D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH, PERIODIC, DIRICHLET, NEUMANN
   implicit none(type, external)
   private
   public :: check_boundary_real, check_boundary_integer

contains

   !> Check that boundary faces have correct values based on their boundary condition type
   !> - Periodic (0): Handled by test_halo
   !> - Dirichlet (1): Should be constant value 1.0
   !> - Neumann (2): Should match adjacent interior cell
   subroutine check_boundary_real(domain, array, bc_types, dirichlet_values)
      class(mpi_domain_t), intent(in) :: domain
      real(kind=sp), dimension(:, :, :), intent(in) :: array
      !! Assumed-shape array extends from 1:n. `nx(|y|z)` *includes* halos
      integer, intent(in) :: bc_types(6)
      real(kind=sp), intent(in) :: dirichlet_values(6)
      real(kind=sp), parameter :: tolerance = 1.0e-6_sp
      logical :: is_bc_face(6)

      is_bc_face = domain%is_boundary_face(:)

      ! West face
      if (is_bc_face(D_WEST)) then
         select case (bc_types(D_WEST))
         case (DIRICHLET)
            if (any(abs(array(1, 2:size(array, 2) - 1, 2:size(array, 3) - 1) - dirichlet_values(D_WEST)) > tolerance)) then
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
         select case (bc_types(D_EAST))
         case (DIRICHLET)
            if (any(abs(array(ubound(array, 1), 2:size(array, 2) - 1, 2:size(array, 3) - 1) - &
               dirichlet_values(D_EAST)) > tolerance)) then
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
         select case (bc_types(D_SOUTH))
         case (DIRICHLET)
            if (any(abs(array(2:size(array, 1) - 1, 1, 2:size(array, 3) - 1) - dirichlet_values(D_SOUTH)) > tolerance)) then
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
         select case (bc_types(D_NORTH))
         case (DIRICHLET)
            if (any(abs(array(2:size(array, 1) - 1, ubound(array, 2), 2:size(array, 3) - 1) - &
               dirichlet_values(D_NORTH)) > tolerance)) then
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
         select case (bc_types(D_LOW))
         case (DIRICHLET)
            if (any(abs(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, 1) - dirichlet_values(D_LOW)) > tolerance)) then
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
         select case (bc_types(D_HIGH))
         case (DIRICHLET)
            if (any(abs(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, ubound(array, 3)) - &
               dirichlet_values(D_HIGH)) > tolerance)) then
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

   subroutine check_boundary_integer(domain, array, bc_types, dirichlet_values)
      class(mpi_domain_t), intent(in) :: domain
      integer, dimension(:, :, :), intent(in) :: array
      integer, intent(in) :: bc_types(6)
      real(kind=sp), intent(in) :: dirichlet_values(6)
      logical :: is_bc_face(6)

      is_bc_face = domain%is_boundary_face(:)

      ! West face
      if (is_bc_face(D_WEST)) then
         select case (bc_types(D_WEST))
         case (DIRICHLET)
            if (any(array(1, 2:size(array, 2) - 1, 2:size(array, 3) - 1) /= int(dirichlet_values(D_WEST)))) then
               call domain%abort("TEST Boundary INT: Not OK: West face Dirichlet")
            end if
         case (NEUMANN)
            if (any(array(1, 2:size(array, 2) - 1, 2:size(array, 3) - 1) /= &
                  array(2, 2:size(array, 2) - 1, 2:size(array, 3) - 1))) then
               call domain%abort("TEST Boundary INT: Not OK: West face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary INT: BC index not recognized")
         end select
      end if

      ! East face
      if (is_bc_face(D_EAST)) then
         select case (bc_types(D_EAST))
         case (DIRICHLET)
            if (any(array(size(array, 1), 2:size(array, 2) - 1, 2:size(array, 3) - 1) /= int(dirichlet_values(D_EAST)))) then
               call domain%abort("TEST Boundary INT: Not OK: East face Dirichlet")
            end if
         case (NEUMANN)
            if (any(array(size(array, 1), 2:size(array, 2) - 1, 2:size(array, 3) - 1) /= &
                  array(size(array, 1) - 1, 2:size(array, 2) - 1, 2:size(array, 3) - 1))) then
               call domain%abort("TEST Boundary INT: Not OK: East face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary INT: BC index not recognized")
         end select
      end if

      ! South face
      if (is_bc_face(D_SOUTH)) then
         select case (bc_types(D_SOUTH))
         case (DIRICHLET)
            if (any(array(2:size(array, 1) - 1, 1, 2:size(array, 3) - 1) /= int(dirichlet_values(D_SOUTH)))) then
               call domain%abort("TEST Boundary INT: Not OK: South face Dirichlet")
            end if
         case (NEUMANN)
            if (any(array(2:size(array, 1) - 1, 1, 2:size(array, 3) - 1) /= &
                  array(2:size(array, 1) - 1, 2, 2:size(array, 3) - 1))) then
               call domain%abort("TEST Boundary INT: Not OK: South face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary INT: BC index not recognized")
         end select
      end if

      ! North face
      if (is_bc_face(D_NORTH)) then
         select case (bc_types(D_NORTH))
         case (DIRICHLET)
            if (any(array(2:size(array, 1) - 1, size(array, 2), 2:size(array, 3) - 1) /= int(dirichlet_values(D_NORTH)))) then
               call domain%abort("TEST Boundary INT: Not OK: North face Dirichlet")
            end if
         case (NEUMANN)
            if (any(array(2:size(array, 1) - 1, size(array, 2), 2:size(array, 3) - 1) /= &
                  array(2:size(array, 1) - 1, size(array, 2) - 1, 2:size(array, 3) - 1))) then
               call domain%abort("TEST Boundary INT: Not OK: North face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary INT: BC index not recognized")
         end select
      end if

      ! Low face
      if (is_bc_face(D_LOW)) then
         select case (bc_types(D_LOW))
         case (DIRICHLET)
            if (any(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, 1) /= int(dirichlet_values(D_LOW)))) then
               call domain%abort("TEST Boundary INT: Not OK: Low face Dirichlet")
            end if
         case (NEUMANN)
            if (any(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, 1) /= &
                  array(2:size(array, 1) - 1, 2:size(array, 2) - 1, 2))) then
               call domain%abort("TEST Boundary INT: Not OK: Low face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary INT: BC index not recognized")
         end select
      end if

      ! High face
      if (is_bc_face(D_HIGH)) then
         select case (bc_types(D_HIGH))
         case (DIRICHLET)
            if (any(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, size(array, 3)) /= int(dirichlet_values(D_HIGH)))) then
               call domain%abort("TEST Boundary INT: Not OK: High face Dirichlet")
            end if
         case (NEUMANN)
            if (any(array(2:size(array, 1) - 1, 2:size(array, 2) - 1, size(array, 3)) /= &
                  array(2:size(array, 1) - 1, 2:size(array, 2) - 1, size(array, 3) - 1))) then
               call domain%abort("TEST Boundary INT: Not OK: High face Neumann")
            end if
         case default
            call domain%abort("TEST Boundary INT: BC index not recognized")
         end select
      end if
   end subroutine check_boundary_integer
end module test_boundary
