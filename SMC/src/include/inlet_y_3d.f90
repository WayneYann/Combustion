  subroutine update_inlet_ylo_3d(lo,hi,qin,t,dx)
    integer, intent(in) :: lo(3), hi(3)
    double precision, intent(in) :: t, dx(3)
    double precision, intent(inout) :: qin(nqin,lo(1):hi(1),lo(3):hi(3))
    return
  end subroutine update_inlet_ylo_3d
