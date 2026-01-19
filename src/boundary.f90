!> Applies boundary conditions in a 3D MPI simulation.
!> - Detects physical boundary faces per rank.
!> - Supports periodic, Dirichlet, and Neumann conditions.
module boundary
   use mpi_f08, only: MPI_Abort, MPI_COMM_WORLD
   use lib_mpi_precision, only: sp
   use mpi_domain_types, only: mpi_domain_t
   use lib_mpi_parameters, only: boundaries, nx => num_cells_x, ny => num_cells_y, nz => num_cells_z
   use lib_mpi_enums, only: D_WEST, D_EAST, D_SOUTH, D_NORTH, D_LOW, D_HIGH, PERIODIC, DIRICHLET, NEUMANN
   implicit none(type, external)
   private
   public :: update_boundaries

   logical, parameter :: DEBUG = .true.
      !! Print verbose debug statements

   !! Generic interfaces that wrap the real(sp) and integer routines
   interface update_boundaries
      module procedure update_boundary_conditions_real
      module procedure update_boundary_conditions_int
   end interface update_boundaries

   interface apply_dirichlet_bc
      module procedure apply_dirichlet_bc_real
      module procedure apply_dirichlet_bc_int
   end interface apply_dirichlet_bc

   interface apply_neumann_bc
      module procedure apply_neumann_bc_real
      module procedure apply_neumann_bc_int
   end interface apply_neumann_bc

contains
   !> Applies boundary conditions to physical boundary faces only.
   !>
   !> Cell layout (1-indexed in callee due to assumed-shape remapping):
   !>   Index:    1     2    3   ...  n-1   n
   !>   Cell:   ghost | interior cells    | ghost
   !>
   !> Face ownership (mutually exclusive):
   !>   - neighbor /= `MPI_PROC_NULL`  →  halo exchange owns the ghost cell
   !>   - neighbor == `MPI_PROC_NULL`  →  boundary condition owns the ghost cell
   subroutine update_boundary_conditions_real(domain, array, bc_types, dirichlet_values)
      class(mpi_domain_t), intent(in) :: domain
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
         !! Source array to update boundary conditions
      integer, intent(in) :: bc_types(6)
         !! The boundary conditions for each face/direction
      integer :: face
         !! Local loop variable - West, East, South, North, Low, High
      real(kind=sp), intent(in) :: dirichlet_values(6)

      ! Check if rank is internal to domain - no boundaries
      if (domain%is_interior) return

      do face = 1, 6
         ! Check if face is on the boundary
         if (.not. domain%is_boundary_face(face)) cycle

         select case (bc_types(face))
         case (PERIODIC)
            ! Periodic is handled by MPI (neighbor is not MPI_PROC_NULL)
            ! This case technically shouldn't be reached if is_boundary_face is true.
            cycle
         case (DIRICHLET)
            call apply_dirichlet_bc(domain, array, face, constant_value=dirichlet_values(face))
         case (NEUMANN)
            call apply_neumann_bc(domain, array, face)
         case default
            call domain%abort("BC Not implemented yet. Exiting..")
         end select
      end do

   end subroutine update_boundary_conditions_real

   !> Applies boundary conditions to physical boundary faces only.
   !>
   !> Cell layout (1-indexed in callee due to assumed-shape remapping):
   !>   Index:    1     2    3   ...  n-1   n
   !>   Cell:   ghost | interior cells    | ghost
   !>
   !> Face ownership (mutually exclusive):
   !>   - neighbor /= `MPI_PROC_NULL`  →  halo exchange owns the ghost cell
   !>   - neighbor == `MPI_PROC_NULL`  →  boundary condition owns the ghost cell
   subroutine update_boundary_conditions_int(domain, array, bc_types, dirichlet_values)
      class(mpi_domain_t), intent(in) :: domain
      integer, contiguous, intent(in out) :: array(:, :, :)
         !! Source array to update boundary conditions
      integer, intent(in) :: bc_types(6)
         !! The boundary conditions for each face/direction
      real(kind=sp), intent(in) :: dirichlet_values(6)
      integer :: face
         !! Local loop variable - West, East, South, North, Low, High

      ! Check if rank is internal to domain - no boundaries
      if (domain%is_interior) return

      do face = 1, 6
         ! Check if face is on the boundary
         if (.not. domain%is_boundary_face(face)) cycle

         select case (bc_types(face))
         case (PERIODIC)
            ! Periodic is handled by MPI (neighbor is not MPI_PROC_NULL)
            ! This case technically shouldn't be reached if is_boundary_face is true.
            cycle
         case (DIRICHLET)
            call apply_dirichlet_bc(domain, array, face, constant_value=int(dirichlet_values(face)))
         case (NEUMANN)
            call apply_neumann_bc(domain, array, face)
         case default
            call domain%abort("BC Not implemented yet. Exiting..")
         end select
      end do

   end subroutine update_boundary_conditions_int

   !> Applies a Dirichlet boundary condition by setting boundary cells to a constant.
   subroutine apply_dirichlet_bc_real(domain, array, face, constant_value)
      class(mpi_domain_t), intent(in) :: domain
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
      integer, intent(in) :: face
      real(kind=sp), intent(in) :: constant_value
         !! Per-face Dirichlet values
      integer :: rank

      rank = domain%get_rank()
      select case (face)
      case (D_WEST)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_WEST, " for rank:", rank
         array(1, :, :) = constant_value
      case (D_EAST)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_EAST, " for rank:", rank
         array(ubound(array, 1), :, :) = constant_value
      case (D_SOUTH)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_SOUTH, " for rank:", rank
         array(:, 1, :) = constant_value
      case (D_NORTH)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_NORTH, " for rank:", rank
         array(:, ubound(array, 2), :) = constant_value
      case (D_LOW)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_LOW, " for rank:", rank
         array(:, :, 1) = constant_value
      case (D_HIGH)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_HIGH, " for rank:", rank
         array(:, :, ubound(array, 3)) = constant_value
      case default
         call domain%abort("Dirichlet boundary condition direction index not standard... Exiting..")
      end select
   end subroutine apply_dirichlet_bc_real

      !> Applies a Dirichlet boundary condition by setting boundary cells to a constant.
   subroutine apply_dirichlet_bc_int(domain, array, face, constant_value)
      class(mpi_domain_t), intent(in) :: domain
      integer, contiguous, intent(in out) :: array(:, :, :)
      integer, intent(in) :: face
      integer, intent(in) :: constant_value
         !! Per-face Dirichlet values
      integer :: rank

      rank = domain%get_rank()
      select case (face)
      case (D_WEST)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_WEST, " for rank:", rank
         array(1, :, :) = constant_value
      case (D_EAST)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_EAST, " for rank:", rank
         array(ubound(array, 1), :, :) = constant_value
      case (D_SOUTH)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_SOUTH, " for rank:", rank
         array(:, 1, :) = constant_value
      case (D_NORTH)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_NORTH, " for rank:", rank
         array(:, ubound(array, 2), :) = constant_value
      case (D_LOW)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_LOW, " for rank:", rank
         array(:, :, 1) = constant_value
      case (D_HIGH)
         if (DEBUG) write (*, "(2(A,1X,I0))") "BC DIRICHLET: Updating face", D_HIGH, " for rank:", rank
         array(:, :, ubound(array, 3)) = constant_value
      case default
         call domain%abort("Dirichlet boundary condition direction index not standard... Exiting..")
      end select
   end subroutine apply_dirichlet_bc_int

   !> Applies a Neumann boundary by copying data from the adjacent interior cell.
   subroutine apply_neumann_bc_real(domain, array, face)
      class(mpi_domain_t), intent(in) :: domain
      real(kind=sp), contiguous, intent(in out) :: array(:, :, :)
      integer, intent(in) :: face

      select case (face)
      case (D_WEST); array(1, :, :) = array(2, :, :)
      case (D_EAST); array(ubound(array, 1), :, :) = array(ubound(array, 1) - 1, :, :)
      case (D_SOUTH); array(:, 1, :) = array(:, 2, :)
      case (D_NORTH); array(:, ubound(array, 2), :) = array(:, ubound(array, 2) - 1, :)
      case (D_LOW); array(:, :, 1) = array(:, :, 2)
      case (D_HIGH); array(:, :, ubound(array, 3)) = array(:, :, ubound(array, 3) - 1)
      case default
         call domain%abort("Neumann boundary condition direction index not standard... Exiting..")
      end select
   end subroutine apply_neumann_bc_real

      !> Applies a Neumann boundary by copying data from the adjacent interior cell.
   subroutine apply_neumann_bc_int(domain, array, face)
      class(mpi_domain_t), intent(in) :: domain
      integer, contiguous, intent(in out) :: array(:, :, :)
      integer, intent(in) :: face

      select case (face)
      case (D_WEST); array(1, :, :) = array(2, :, :)
      case (D_EAST); array(ubound(array, 1), :, :) = array(ubound(array, 1) - 1, :, :)
      case (D_SOUTH); array(:, 1, :) = array(:, 2, :)
      case (D_NORTH); array(:, ubound(array, 2), :) = array(:, ubound(array, 2) - 1, :)
      case (D_LOW); array(:, :, 1) = array(:, :, 2)
      case (D_HIGH); array(:, :, ubound(array, 3)) = array(:, :, ubound(array, 3) - 1)
      case default
         call domain%abort("Neumann boundary condition direction index not standard... Exiting..")
      end select
   end subroutine apply_neumann_bc_int

end module boundary
