module convert_module

  implicit none

  private

  public :: cellavg2cc_1d, cc2cellavg_1d

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

  subroutine cc2cellavg_1d(lo, hi, uc, clo, chi, ua, alo, ahi)
    integer, intent(in) :: lo, hi, alo, ahi, clo, chi
    double precision, intent(in) :: uc(clo:chi)
    double precision             :: ua(alo:ahi)
    
    integer :: i
    double precision, parameter :: b = 1.d0/24.d0, d = 11.d0/12.d0

    do i=lo,hi
       ua(i) = b*(uc(i-1)+uc(i+1)) + d*uc(i) 
    end do

  end subroutine cc2cellavg_1d

end module convert_module
