module convert_2d_module

  implicit none

  private

  public :: cellavg2cc_2d

contains

  subroutine cellavg2cc_2d(lo, hi, ua, alo, ahi, uc, clo, chi)
    integer, intent(in) :: lo(2), hi(2), alo(2), ahi(2), clo(2), chi(2)
    double precision, intent(in) :: ua(alo(1):ahi(1),alo(2):ahi(2))
    double precision             :: uc(clo(1):chi(1),clo(2):chi(2))
    
    integer :: i, j
    double precision, parameter :: b = -1.d0/24.d0, d = 13.d0/12.d0, sevenSixth=7.d0/6.d0

    if (lo(1) .eq. hi(1)) then

       i = lo(1)
       do j=lo(2),hi(2)
          uc(i,j) = b*(ua(i,j-1)+ua(i,j+1)) + d*ua(i,j) 
       end do

    else if (lo(2) .eq. hi(2)) then

       j = lo(2)
       do i=lo(1),hi(1)
          uc(i,j) = b*(ua(i-1,j)+ua(i+1,j)) + d*ua(i,j) 
       end do

    else

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uc(i,j) = b*(ua(i-1,j-1)+ua(i+1,j-1)+ua(i-1,j+1)+ua(i+1,j+1)) &
                  + sevenSixth*ua(i,j)
          end do
       end do

    end if

  end subroutine cellavg2cc_2d

end module convert_2d_module
