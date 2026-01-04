module mod_kinds
  ! use mod_kinds, only: dp
  use iso_fortran_env

implicit none
  ! single precision (optional)
  integer, parameter :: sp = selected_real_kind(6, 37)
  ! double precision: ~15 digits, range ~1e300+
  integer, parameter :: dp = selected_real_kind(15, 300)
  ! quadruple precision (~33 decimal digits, if supported)
  integer, parameter :: qp = selected_real_kind(33, 4931)

  integer, parameter :: i4 = int32
  integer, parameter :: i8 = int64

end module mod_kinds
