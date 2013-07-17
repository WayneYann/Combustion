  subroutine update_inlet_xlo_2d(lo,hi,qin,t)
    integer, intent(in) :: lo(2), hi(2)
    double precision, intent(in) :: t
    double precision, intent(inout) :: qin(nqin,lo(2):hi(2))
    return
  end subroutine update_inlet_xlo_2d
