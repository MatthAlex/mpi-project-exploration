!> Module for applying boundary conditions in a 3D MPI code.
!>
!> This module defines:
!> - Named constants for each of the 6 directions in 3D
!> - Internal variables for boundary faces
!> - Routines to determine which faces belong to physical domain boundaries
!> - Routines to apply various boundary conditions (Dirichlet, Neumann, etc.)
module boundary
   use mpi_f08, only: MPI_PROC_NULL, MPI_Finalize
   use precision, only: sp
   use grid_module, only: rank => my_rank, neighbor_ranks
   use lib_parameters, only: boundaries, nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   use enums, only: D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH, PERIODIC, DIRICHLET, NEUMANN
   implicit none(type, external)
   private
   public :: determine_rank_boundaries, update_boundary_conditions

   logical, protected, public :: is_bc_face(6)
      !! Indicates if the local rank is a physical boundary, for each face.
   logical, protected, public :: is_rank_inside
      !! Indicates if this rank is completely inside (no boundary faces).
   logical, parameter :: DEBUG = .true.
      !! Print verbose debug statements

contains

   !> Determines which faces of this rank are physical boundaries.
   !>
   !> Sets:
   !> - `is_bc_face(i) = .true.` if the face i has `MPI_PROC_NULL` neighbor
   !> - `rank_is_inside = .true.` if no face is a boundary
   subroutine determine_rank_boundaries()
      is_bc_face = .false.

      ! Set boundary face booleans for this rank
      where (neighbor_ranks == MPI_PROC_NULL)
         is_bc_face = .true.
      end where

      is_rank_inside = .not. any(is_bc_face)
   end subroutine determine_rank_boundaries

   !> Applies boundary conditions to each boundary face of the local domain,
   !> only when the face belongs to the physical boundary.
   subroutine update_boundary_conditions(array, bc_types)
      use lib_parameters, only: dirichlet_value
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
         !! Source array to update boundary conditions
      integer, intent(in) :: bc_types(6)
         !! The boundary conditions for each face/direction
      integer :: face
         !! Local loop variable - West, East, South, North, Low, High

      if (is_rank_inside) return

      do face = 1, 6
         if (.not. is_bc_face(face)) cycle

         select case (bc_types(face))
         case (PERIODIC)
            ! Periodic is handled by MPI
            cycle
         case (DIRICHLET)
            call apply_dirichlet_bc(array, face, constant_value=dirichlet_value)
         case (NEUMANN)
            call apply_neumann_bc(array, face)
         case default
            call apply_custom_bc(array, face)
         end select
      end do

   end subroutine update_boundary_conditions

   !> Applies a Dirichlet boundary condition by setting boundary cells to a constant.
   subroutine apply_dirichlet_bc(array, face, constant_value)
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
      integer, intent(in) :: face
      real(kind=sp), intent(in) :: constant_value
         !! Dirichlet constant value. Static across all faces

      select case (face)
      case (D_WEST)
         if (DEBUG) write (*, '(2(A,1X,I0))') "BC DIRICHLET: Updating face", D_WEST, " for rank:", rank
         array(1, :, :) = constant_value
      case (D_EAST)
         if (DEBUG) write (*, '(2(A,1X,I0))') "BC DIRICHLET: Updating face", D_EAST, " for rank:", rank
         array(ubound(array, 1), :, :) = constant_value
      case (D_SOUTH)
         if (DEBUG) write (*, '(2(A,1X,I0))') "BC DIRICHLET: Updating face", D_SOUTH, " for rank:", rank
         array(:, 1, :) = constant_value
      case (D_NORTH)
         if (DEBUG) write (*, '(2(A,1X,I0))') "BC DIRICHLET: Updating face", D_NORTH, " for rank:", rank
         array(:, ubound(array, 2), :) = constant_value
      case (D_LOW)
         if (DEBUG) write (*, '(2(A,1X,I0))') "BC DIRICHLET: Updating face", D_LOW, " for rank:", rank
         array(:, :, 1) = constant_value
      case (D_HIGH)
         if (DEBUG) write (*, '(2(A,1X,I0))') "BC DIRICHLET: Updating face", D_HIGH, " for rank:", rank
         array(:, :, ubound(array, 3)) = constant_value
      end select
   end subroutine apply_dirichlet_bc

   !> Applies a Neumann boundary by copying data from the adjacent interior cell.
   subroutine apply_neumann_bc(array, face)
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
      integer, intent(in) :: face

      select case (face)
      case (D_WEST); array(1, :, :) = array(2, :, :)
      case (D_EAST); array(ubound(array, 1), :, :) = array(ubound(array, 1) - 1, :, :)
      case (D_SOUTH); array(:, 1, :) = array(:, 2, :)
      case (D_NORTH); array(:, ubound(array, 2), :) = array(:, ubound(array, 2) - 1, :)
      case (D_LOW); array(:, :, 1) = array(:, :, 2)
      case (D_HIGH); array(:, :, ubound(array, 3)) = array(:, :, ubound(array, 3) - 1)
      end select
   end subroutine apply_neumann_bc

   subroutine apply_custom_bc(array, face)
      real(kind=sp), intent(inout) :: array(:, :, :)
      integer, intent(in) :: face
      integer :: ierr

      print *, "Custom BC Not implemented yet. Exiting.."
      call MPI_Finalize(ierr)
   end subroutine apply_custom_bc

end module boundary
