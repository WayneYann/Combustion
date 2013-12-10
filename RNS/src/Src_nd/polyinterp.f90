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

  public :: cc2xface_2d, cc2yface_2d, cc2zgauss_3d, cc2zface_3d

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


  subroutine cc2zgauss_3d(lo,hi, uc, uclo, uchi, u1, u2, uzlo,uzhi)
    integer, intent(in) :: lo(3), hi(3), uclo(3), uchi(3), uzlo(3), uzhi(3)
    double precision, intent(in ) :: uc(uclo(1):uchi(1),uclo(2):uchi(2),uclo(3):uchi(3))
    double precision, intent(out) :: u1(uzlo(1):uzhi(1),uzlo(2):uzhi(2),uzlo(3):uzhi(3))
    double precision, intent(out) :: u2(uzlo(1):uzhi(1),uzlo(2):uzhi(2),uzlo(3):uzhi(3))

    integer :: i, j, k

    do k      =lo(3),hi(3)
       do j   =lo(2),hi(2)
          do i=lo(1),hi(1)
             u1(i,j,k) = cg1(-2)*uc(i,j,k-2)+cg1(-1)*uc(i,j,k-1)+cg1(0)*uc(i,j,k) &
                  +      cg1( 1)*uc(i,j,k+1)+cg1( 2)*uc(i,j,k+2)
             u2(i,j,k) = cg2(-2)*uc(i,j,k-2)+cg2(-1)*uc(i,j,k-1)+cg2(0)*uc(i,j,k) &
                  +      cg2( 1)*uc(i,j,k+1)+cg2( 2)*uc(i,j,k+2)
          end do
       end do
    end do

  end subroutine cc2zgauss_3d


  subroutine cc2zface_3d(lo,hi, u, ulo, uhi, u1, u2, u3, u4, flo,fhi)
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),flo(3),fhi(3)
    double precision, intent(in ) :: u (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision, intent(out) :: u1(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    double precision, intent(out) :: u2(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    double precision, intent(out) :: u3(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    double precision, intent(out) :: u4(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i, j, k
    double precision, allocatable :: uz(:,:,:), uy1(:,:,:), uy2(:,:,:)
    
    allocate(uz (lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3):hi(3)+1))
    allocate(uy1(lo(1)-2:hi(1)+2,lo(2)  :hi(2)  ,lo(3):hi(3)+1))
    allocate(uy2(lo(1)-2:hi(1)+2,lo(2)  :hi(2)  ,lo(3):hi(3)+1))

    do k      =lo(3)  ,hi(3)+1
       do j   =lo(2)-2,hi(2)+2
          do i=lo(1)-2,hi(1)+2
             uz(i,j,k) = -0.0625d0*(u(i,j,k-2)+u(i,j,k+1)) + 0.5625d0*(u(i,j,k-1)+u(i,j,k))
          end do
       end do
    end do

    do k      =lo(3)  ,hi(3)+1
       do j   =lo(2)  ,hi(2)
          do i=lo(1)-2,hi(1)+2
             uy1(i,j,k) = cg1(-2)*uz(i,j-2,k)+cg1(-1)*uz(i,j-1,k)+cg1(0)*uz(i,j,k) &
                  +       cg1( 1)*uz(i,j+1,k)+cg1( 2)*uz(i,j+2,k)
             uy2(i,j,k) = cg2(-2)*uz(i,j-2,k)+cg2(-1)*uz(i,j-1,k)+cg2(0)*uz(i,j,k) &
                  +       cg2( 1)*uz(i,j+1,k)+cg2( 2)*uz(i,j+2,k)
          end do
       end do
    end do

    do k      =lo(3)  ,hi(3)+1
       do j   =lo(2)  ,hi(2)
          do i=lo(1)  ,hi(1)
             u1(i,j,k) = cg1(-2)*uy1(i-2,j,k)+cg1(-1)*uy1(i-1,j,k)+cg1(0)*uy1(i,j,k) &
                  +      cg1( 1)*uy1(i+1,j,k)+cg1( 2)*uy1(i+2,j,k)
             u2(i,j,k) = cg2(-2)*uy1(i-2,j,k)+cg2(-1)*uy1(i-1,j,k)+cg2(0)*uy1(i,j,k) &
                  +      cg2( 1)*uy1(i+1,j,k)+cg2( 2)*uy1(i+2,j,k)

             u3(i,j,k) = cg1(-2)*uy2(i-2,j,k)+cg1(-1)*uy2(i-1,j,k)+cg1(0)*uy2(i,j,k) &
                  +      cg1( 1)*uy2(i+1,j,k)+cg1( 2)*uy2(i+2,j,k)
             u4(i,j,k) = cg2(-2)*uy2(i-2,j,k)+cg2(-1)*uy2(i-1,j,k)+cg2(0)*uy2(i,j,k) &
                  +      cg2( 1)*uy2(i+1,j,k)+cg2( 2)*uy2(i+2,j,k)
          end do
       end do
    end do

    deallocate(uz,uy1,uy2)

  end subroutine cc2zface_3d

end module polyinterp_module
