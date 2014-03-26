module convert_module

  implicit none

  private

  public :: cellavg2cc_1d

contains

  subroutine cellavg2cc_1d(lo, hi, ua, alo, ahi, uc, clo, chi)
    integer, intent(in) :: lo, hi, alo, ahi, clo, chi
    double precision, intent(in) :: ua(alo:ahi)
    double precision, intent(out) :: uc(clo:chi)
    
    integer :: i
    double precision, parameter :: b = -1.d0/24.d0, d = 13.d0/12.d0

    do i=lo,hi
       uc(i) = b*ua(i-1) + d*ua(i) + b*ua(i+1) 
    end do

  end subroutine cellavg2cc_1d

end module convert_module
