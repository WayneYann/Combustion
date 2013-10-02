module convert_3d_module

  implicit none

  private

  public :: cellavg2cc_3d, cc2cellavg_3d

contains

  subroutine cellavg2cc_3d(lo, hi, ua, alo, ahi, uc, clo, chi, idir)
    integer, intent(in) :: lo(3), hi(3), alo(3), ahi(3), clo(3), chi(3)
    double precision, intent(in) :: ua(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    double precision             :: uc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    integer, intent(in), optional :: idir

    integer :: i, j, k, dm
    double precision, parameter :: b = -1.d0/24.d0, d = 13.d0/12.d0

    if (present(idir)) then
       dm = idir
    else
       dm = 0
    end if

    if (dm .eq. 1) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                uc(i,j,k) = b*(ua(i-1,j,k)+ua(i+1,j,k)) + d*ua(i,j,k)
             end do
          end do
       end do

    else if (dm .eq. 2) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                uc(i,j,k) = b*(ua(i,j-1,k)+ua(i,j+1,k)) + d*ua(i,j,k)
             end do
          end do
       end do

    else if (dm .eq. 3) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                uc(i,j,k) = b*(ua(i,j,k-1)+ua(i,j,k+1)) + d*ua(i,j,k)
             end do
          end do
       end do

    else

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                uc(i,j,k) = b*(ua(i,j,k-1)+ua(i,j-1,k)+ua(i-1,j,k) &
                     +         ua(i+1,j,k)+ua(i,j+1,k)+ua(i,j,k+1)) &
                     + 1.25d0*ua(i,j,k)
             end do
          end do
       end do

    end if

  end subroutine cellavg2cc_3d

  subroutine cc2cellavg_3d(lo, hi, uc, clo, chi, ua, alo, ahi)
    integer, intent(in) :: lo(3), hi(3), alo(3), ahi(3), clo(3), chi(3)
    double precision, intent(in) :: uc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    double precision             :: ua(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    
    integer :: i, j, k
    double precision, parameter :: b = 1.d0/24.d0, threeFourth=0.75d0

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ua(i,j,k) = b*(uc(i,j,k-1)+uc(i,j-1,k)+uc(i-1,j,k) &
                  +         uc(i+1,j,k)+uc(i,j+1,k)+uc(i,j,k+1)) &
                  + threeFourth*uc(i,j,k)
          end do
       end do
    end do

  end subroutine cc2cellavg_3d

end module convert_3d_module
