  subroutine update_inlet_xlo_1d(lo,hi,qin,t,dx)
    integer, intent(in) :: lo(1), hi(1)
    double precision, intent(in) :: t, dx(1)
    double precision, intent(inout) :: qin(nqin)
  end subroutine update_inlet_xlo_1d
