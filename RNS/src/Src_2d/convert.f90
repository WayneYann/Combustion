module convert_module

  implicit none

  private

  public :: cellavg2cc_2d, cc2cellavg_2d, cellavg2cc_1d, cc2cellavg_1d

contains

  subroutine cellavg2cc_2d(lo, hi, ua, alo, ahi, uc, clo, chi)
    integer, intent(in) :: lo(2), hi(2), alo(2), ahi(2), clo(2), chi(2)
    double precision, intent(in) :: ua(alo(1):ahi(1),alo(2):ahi(2))
    double precision             :: uc(clo(1):chi(1),clo(2):chi(2))

    integer :: i, j
    double precision, parameter :: b = -1.d0/24.d0, sevenSixth=7.d0/6.d0

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          uc(i,j) = b*(ua(i,j-1)+ua(i-1,j)+ua(i+1,j)+ua(i,j+1)) &
               + sevenSixth*ua(i,j)
       end do
    end do
  end subroutine cellavg2cc_2d

  subroutine cc2cellavg_2d(lo, hi, uc, clo, chi, ua, alo, ahi)
    integer, intent(in) :: lo(2), hi(2), alo(2), ahi(2), clo(2), chi(2)
    double precision, intent(in) :: uc(clo(1):chi(1),clo(2):chi(2))
    double precision             :: ua(alo(1):ahi(1),alo(2):ahi(2))
    
    integer :: i, j
    double precision, parameter :: b = 1.d0/24.d0, fiveSixth=5.d0/6.d0

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ua(i,j) = b*(uc(i,j-1)+uc(i-1,j)+uc(i+1,j)+uc(i,j+1)) &
               + fiveSixth*uc(i,j)
       end do
    end do
  end subroutine cc2cellavg_2d

  subroutine cellavg2cc_1d(lo, hi, ua, alo, ahi, uc, clo, chi, dir)
    integer, intent(in) :: lo(2), hi(2), alo(2), ahi(2), clo(2), chi(2)
    double precision, intent(in) :: ua(alo(1):ahi(1),alo(2):ahi(2))
    double precision             :: uc(clo(1):chi(1),clo(2):chi(2))
    integer, intent(in) :: dir

    integer :: i, j
    double precision, parameter :: b = -1.d0/24.d0, d = 13.d0/12.d0

    if (dir .eq. 1) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uc(i,j) = b*(ua(i-1,j)+ua(i+1,j)) + d*ua(i,j) 
          end do
       end do

    else if (dir .eq. 2) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uc(i,j) = b*(ua(i,j-1)+ua(i,j+1)) + d*ua(i,j) 
          end do
       end do

    else

       call bl_error("Error: cellavg2cc_1d: unknown dir")

    end if
  end subroutine cellavg2cc_1d

  subroutine cc2cellavg_1d(lo, hi, uc, clo, chi, ua, alo, ahi, dir)
    integer, intent(in) :: lo(2), hi(2), alo(2), ahi(2), clo(2), chi(2)
    double precision, intent(in) :: uc(clo(1):chi(1),clo(2):chi(2))
    double precision             :: ua(alo(1):ahi(1),alo(2):ahi(2))
    integer, intent(in) :: dir
    
    integer :: i, j
    double precision, parameter :: b = 1.d0/24.d0, d = 11.d0/12.d0

    if (dir .eq. 1) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ua(i,j) = b*(uc(i-1,j)+uc(i+1,j)) + d*uc(i,j) 
          end do
       end do

    else if (dir .eq. 2) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ua(i,j) = b*(uc(i,j-1)+uc(i,j+1)) + d*uc(i,j) 
          end do
       end do

    else

       call bl_error("Error: cc2cellavg_1d: unknown dir")

    end if
  end subroutine cc2cellavg_1d

end module convert_module
