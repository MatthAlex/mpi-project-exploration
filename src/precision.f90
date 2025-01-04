module precision
   use iso_fortran_env, only: real32, real64
   implicit none
   public
   integer, parameter :: sp = real32
   integer, parameter :: dp = real64
end module precision
