  subroutine update_inlet_xlo_3d(lo,hi,qin,t)
    integer, intent(in) :: lo(3), hi(3)
    double precision, intent(in) :: t
    double precision, intent(inout) :: qin(nqin,lo(3):hi(3),lo(3):hi(3))
    return
  end subroutine update_inlet_xlo_3d
