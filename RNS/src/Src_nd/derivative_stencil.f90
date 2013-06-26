module derivative_stencil_module

  implicit none

  public

  ! 4th-order first derivative on face i-1/2 is dot_product(FD4, u(i-2:i+1))/h,
  !    where u is cell-centered
  double precision,dimension(-2:1),parameter :: FD4 = &
       (/  1.d0/24.d0, -9.d0/8.d0, 9.d0/8.d0, -1.d0/24.d0  /)

end module derivative_stencil_module
