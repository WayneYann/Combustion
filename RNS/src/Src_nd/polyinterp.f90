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

  public :: cc2xface_2d, cc2yface_2d

contains

  ! subroutine cc2face_1d(lo, hi, u, ulo, uhi, uf, flo, fhi)
  !   integer, intent(in) :: lo, hi, ulo, uhi, flo, fhi
  !   double precision, intent(in) ::  u(ulo:uhi)
  !   double precision             :: uf(flo:fhi)
  !   integer :: i
  !   do i=lo,hi
  !      uf(i) = -0.0625d0*(u(i-2)+u(i+1)) + 0.5625d0*(u(i-1)+u(i))
  !   end do
  ! end subroutine cc2face_1d

  subroutine cc2xface_2d(clo,chi,u,ulo,uhi,u1,u2,flo,fhi)
    integer, intent(in) :: clo(2),chi(2),ulo(2),uhi(2),flo(2),fhi(2)
    double precision, intent(in ) :: u (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision, intent(out) :: u1(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: u2(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i, j
    double precision :: uf(clo(1):chi(1)+1,clo(2)-2:chi(2)+2)

    do  j=clo(2)-2, chi(2)+2
     do i=clo(1)   ,chi(1)+1
          uf(i,j) = -0.0625d0*(u(i-2,j)+u(i+1,j)) + 0.5625d0*(u(i-1,j)+u(i,j))
       end do
    end do

    do    j=clo(2),chi(2)
       do i=clo(1),chi(1)+1
          u1(i,j) = cg1(-2)*uf(i,j-2)+cg1(-1)*uf(i,j-1)+cg1(0)*uf(i,j) &
               +    cg1( 1)*uf(i,j+1)+cg1( 2)*uf(i,j+2)
          u2(i,j) = cg2(-2)*uf(i,j-2)+cg2(-1)*uf(i,j-1)+cg2(0)*uf(i,j) &
               +    cg2( 1)*uf(i,j+1)+cg2( 2)*uf(i,j+2)
       end do
    end do

  end subroutine cc2xface_2d

  subroutine cc2yface_2d(clo,chi,u,ulo,uhi,u1,u2,flo,fhi)
    integer, intent(in) :: clo(2),chi(2),ulo(2),uhi(2),flo(2),fhi(2)
    double precision, intent(in ) :: u (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision, intent(out) :: u1(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: u2(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i, j
    double precision :: uf(clo(1)-2:chi(1)+2,clo(2):chi(2)+1)

    do  j=clo(2)  ,chi(2)+1
     do i=clo(1)-2,chi(1)+2
          uf(i,j) = -0.0625d0*(u(i,j-2)+u(i,j+1)) + 0.5625d0*(u(i,j-1)+u(i,j))
       end do
    end do

    do    j=clo(2),chi(2)+1
       do i=clo(1),chi(1)
          u1(i,j) = cg1(-2)*uf(i-2,j)+cg1(-1)*uf(i-1,j)+cg1(0)*uf(i,j) &
               +    cg1( 1)*uf(i+1,j)+cg1( 2)*uf(i+2,j)
          u2(i,j) = cg2(-2)*uf(i-2,j)+cg2(-1)*uf(i-1,j)+cg2(0)*uf(i,j) &
               +    cg2( 1)*uf(i+1,j)+cg2( 2)*uf(i+2,j)
       end do
    end do

  end subroutine cc2yface_2d

end module polyinterp_module
