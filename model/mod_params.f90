      module mod_params
       use mod_kinds, only: dp
       implicit none
       public
       real(dp), parameter :: pi        = 3.1415926535897932384626433832795_dp
       real(dp), parameter :: tolerance = 1.0e-12_dp
      end module mod_params
