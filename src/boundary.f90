!> Module for applying boundary conditions in a 3D MPI code.
!>
!> This module defines:
!> - Named constants for each of the 6 directions in 3D
!> - Internal variables for boundary faces
!> - Routines to determine which faces belong to physical domain boundaries
!> - Routines to apply various boundary conditions (Dirichlet, Neumann, etc.)
module boundary
   use mpi
   use grid_module, only: west, east, south, north, low, high
   use lib_parameters, only: boundaries, nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   implicit none
   private
   public :: determine_rank_boundaries, apply_boundaries
   public :: is_rank_inside, is_bc_face

   !> Named constants for the 6 directions in a 3D domain.
   integer, parameter :: D_WEST = 1, D_EAST = 2, D_SOUTH = 3, D_NORTH = 4, D_LOW = 5, D_HIGH = 6

   !> Module array of neighbor ranks for each of the 6 faces in 3D.
   integer :: neighbor_ranks(6)

   !> Whether the local rank has a physical boundary on each face.
   logical :: is_bc_face(6)

   !> Indicates if this rank is completely inside (no boundary faces).
   logical :: is_rank_inside

contains

   !> Determines which faces of this rank are physical boundaries.
   !>
   !> Sets:
   !> - `is_bc_face(i) = .true.` if the face i has `MPI_PROC_NULL` neighbor
   !> - `rank_is_inside = .true.` if no face is a boundary
   subroutine determine_rank_boundaries()
      neighbor_ranks = [west, east, south, north, low, high]
      is_bc_face = .false.

      ! Set boundary face booleans for this rank
      where (neighbor_ranks == MPI_PROC_NULL)
         is_bc_face = .true.
      end where

      is_rank_inside = .not. any(is_bc_face)
   end subroutine determine_rank_boundaries

   !> Applies boundary conditions to each boundary face of the local domain,
   !> only when the face belongs to the physical boundary.
   subroutine apply_boundaries(array, bc_types)
      real, intent(inout) :: array(:, :, :)
      integer, intent(in) :: bc_types(6)
      integer :: face

      if (is_rank_inside) return
      ! face is West, East, South, North, Low, High
      do face = 1, 6
         if (.not. is_bc_face(face)) cycle

         select case (bc_types(face))
         case (0) ! Periodic is handled by MPI
            cycle
         case (1) ! Dirichlet
            call apply_dirichlet(array, face, constant_value=1.0)
         case (2) ! Von Neumann
            call apply_neumann(array, face)
         case default ! placeholder
            call apply_custom_bc(array, face)
         end select
      end do

   end subroutine apply_boundaries

   !> Applies a Dirichlet boundary condition by setting boundary cells to a constant.
   subroutine apply_dirichlet(array, face, constant_value)
      real, intent(inout) :: array(:, :, :)
      integer, intent(in) :: face
      real, intent(in) :: constant_value ! this value could be an array of values

      select case (face)
      case (D_WEST); array(1, :, :) = constant_value
      case (D_EAST); array(ubound(array, 1), :, :) = constant_value
      case (D_SOUTH); array(:, 1, :) = constant_value
      case (D_NORTH); array(:, ubound(array, 2), :) = constant_value
      case (D_LOW); array(:, :, 1) = constant_value
      case (D_HIGH); array(:, :, ubound(array, 3)) = constant_value
      end select
   end subroutine apply_dirichlet

   !> Applies a Neumann boundary by copying data from the adjacent interior cell.
   subroutine apply_neumann(array, face)
      real, intent(inout) :: array(:, :, :)
      integer, intent(in) :: face

      select case (face)
      case (D_WEST); array(1, :, :) = array(2, :, :)
      case (D_EAST); array(ubound(array, 1), :, :) = array(ubound(array, 1) - 1, :, :)
      case (D_SOUTH); array(:, 1, :) = array(:, 2, :)
      case (D_NORTH); array(:, ubound(array, 2), :) = array(:, ubound(array, 2) - 1, :)
      case (D_LOW); array(:, :, 1) = array(:, :, 2)
      case (D_HIGH); array(:, :, ubound(array, 3)) = array(:, :, ubound(array, 3) - 1)
      end select
   end subroutine apply_neumann

   subroutine apply_custom_bc(array, face)
      real, intent(inout) :: array(:, :, :)
      integer :: face, ierr

      print *, "Custom BC Not implemented yet. Exiting.."
      call MPI_Finalize(ierr)
   end subroutine apply_custom_bc

end module boundary
