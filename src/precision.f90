!> Sets precision levels for numerical computations.
!> - Defines single (sp) and double (dp) precision types.
module precision
   use, intrinsic :: iso_fortran_env, only: real32, real64
   implicit none(type, external)
   private
   public :: sp, dp
   integer, parameter :: sp = real32
      !! Single-precision real kind - from ISO_FORTRAN_ENV module
   integer, parameter :: dp = real64
      !! Double-precision real kind - from ISO_FORTRAN_ENV module
end module precision
