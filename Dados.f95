module Dados
  implicit none
  integer, parameter            :: d = Selected_Real_kind(15,300), s = Selected_Real_kind(7,38), id = Selected_Int_Kind(8), &
  sd = Selected_Int_Kind(15);
  integer(sd), parameter        :: nmax = 1*10**7
  real(d)                       :: tmax,  u = 0.0_d
  real(d), parameter            :: ld = 1.0_d,  pi = 3.1415926535897932_d, z0=ld/2.0_d
  real(d)                       :: vxc, vyc, xr, dr(2), r
  real(d)                       :: x, y, theta, v, vx, vy, x0, y0, theta0, t
  integer                       :: l, k, velo, raio
  logical                       :: disco(2), track
  character(len=25), parameter  :: path = "/home/henry/A/Data/Raios/"
end module Dados
