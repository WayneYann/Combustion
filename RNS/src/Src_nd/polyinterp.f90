module polyinterp_module

  implicit none

  double precision, dimension(-2:2), parameter :: cg1 = &
       (/ -(1.d0+4.d0*sqrt(3.d0))*11.d0/3456.d0, &
       &  (1.d0+2.d0*sqrt(3.d0))*188.d0/3456.d0, &
       &  3102.d0/3456.d0, &
       &  (1.d0-2.d0*sqrt(3.d0))*188.d0/3456.d0, &
       &  (-1.d0+4.d0*sqrt(3.d0))*11.d0/3456.d0  /)

  double precision, dimension(-2:2), parameter :: cg2 = &
       (/ (-1.d0+4.d0*sqrt(3.d0))*11.d0/3456.d0, &
       &  (1.d0-2.d0*sqrt(3.d0))*188.d0/3456.d0, &
       &  3102.d0/3456.d0, &
       &  (1.d0+2.d0*sqrt(3.d0))*188.d0/3456.d0, &
       &  -(1.d0+4.d0*sqrt(3.d0))*11.d0/3456.d0  /)

  private

  public :: cc2face_1d !, cc2gaussptface_2d

contains

  subroutine cc2face_1d(lo, hi, u, ulo, uhi, uf, flo, fhi)
    integer, intent(in) :: lo, hi, ulo, uhi, flo, fhi
    double precision, intent(in) ::  u(ulo:uhi)
    double precision             :: uf(flo:fhi)
    integer :: i
    do i=lo,hi
       uf(i) = -0.0625d0*(u(i-2)+u(i+1)) + 0.5626d0*(u(i-1)+u(i))
    end do
  end subroutine cc2face_1d


!  subroutine cc2gaussptface_2d(lo, hi, 

end module polyinterp_module
