  subroutine update_inlet_xlo_2d(lo,hi,qin,t,dx)
    integer, intent(in) :: lo(2), hi(2)
    double precision, intent(in) :: t, dx(2)
    double precision, intent(inout) :: qin(nqin,lo(2):hi(2))
  end subroutine update_inlet_xlo_2d
