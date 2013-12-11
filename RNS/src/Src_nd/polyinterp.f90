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

  double precision, dimension(-2:2), parameter :: dg1 = &
       (/ (27.d0+5.d0*sqrt(3.d0))/432.d0, &
       &  (-270.d0-92.d0*sqrt(3.d0))/432.d0, &
       &  174.d0*sqrt(3.d0)/432.d0, &
       &  (270.d0-92.d0*sqrt(3.d0))/432.d0, &
       &  (5.d0*sqrt(3.d0)-27.d0)/432.d0 /)

  double precision, dimension(-2:2), parameter :: dg2 = &
       (/ (27.d0-5.d0*sqrt(3.d0))/432.d0, &
       &  (92.d0*sqrt(3.d0)-270.d0)/432.d0, &
       &  -174.d0*sqrt(3.d0)/432.d0, &
       &  (270.d0+92.d0*sqrt(3.d0))/432.d0, &
       &  (-27.d0-5.d0*sqrt(3.d0))/432.d0 /)

  private

  public :: cc2xface_2d, cc2yface_2d, cc2DyXface_2d, cc2DxYface_2d, &
       cc2zface_3d, cc2zgauss_3d, cc2DzGauss_3d, cc2xygauss_3d, &
       cc2xyGaussDx_3d, cc2xyGaussDy_3d

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

  subroutine cc2xface_2d(lo,hi,u,ulo,uhi,uf1,uf2,flo,fhi,uc1,uc2,clo,chi)
    integer, intent(in) :: lo(2),hi(2),ulo(2),uhi(2),flo(2),fhi(2),clo(2),chi(2)
    double precision, intent(in ) :: u  (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision, intent(out) :: uf1(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: uf2(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: uc1(clo(1):chi(1),clo(2):chi(2))
    double precision, intent(out) :: uc2(clo(1):chi(1),clo(2):chi(2))

    integer :: i, j

    do j   =lo(2)  ,hi(2)
       do i=lo(1)-2,hi(1)+2
          uc1(i,j) = cg1(-2)*u(i,j-2)+cg1(-1)*u(i,j-1)+cg1(0)*u(i,j) &
               +     cg1( 1)*u(i,j+1)+cg1( 2)*u(i,j+2)
          uc2(i,j) = cg2(-2)*u(i,j-2)+cg2(-1)*u(i,j-1)+cg2(0)*u(i,j) &
               +     cg2( 1)*u(i,j+1)+cg2( 2)*u(i,j+2)          
       end do
    end do

    do    j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          uf1(i,j) = -0.0625d0*(uc1(i-2,j)+uc1(i+1,j)) + 0.5625d0*(uc1(i-1,j)+uc1(i,j))
          uf2(i,j) = -0.0625d0*(uc2(i-2,j)+uc2(i+1,j)) + 0.5625d0*(uc2(i-1,j)+uc2(i,j))
       end do
    end do

  end subroutine cc2xface_2d

  subroutine cc2yface_2d(lo,hi,u,ulo,uhi,uf1,uf2,flo,fhi,uc1,uc2,clo,chi)
    integer, intent(in) :: lo(2),hi(2),ulo(2),uhi(2),flo(2),fhi(2),clo(2),chi(2)
    double precision, intent(in ) :: u  (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision, intent(out) :: uf1(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: uf2(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: uc1(clo(1):chi(1),clo(2):chi(2))
    double precision, intent(out) :: uc2(clo(1):chi(1),clo(2):chi(2))

    integer :: i, j

    do j   =lo(2)-2,hi(2)+2
       do i=lo(1)  ,hi(1)
          uc1(i,j) = cg1(-2)*u(i-2,j)+cg1(-1)*u(i-1,j)+cg1(0)*u(i,j) &
               +     cg1( 1)*u(i+1,j)+cg1( 2)*u(i+2,j)
          uc2(i,j) = cg2(-2)*u(i-2,j)+cg2(-1)*u(i-1,j)+cg2(0)*u(i,j) &
               +     cg2( 1)*u(i+1,j)+cg2( 2)*u(i+2,j)          
       end do
    end do

    do    j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          uf1(i,j) = -0.0625d0*(uc1(i,j-2)+uc1(i,j+1)) + 0.5625d0*(uc1(i,j-1)+uc1(i,j))
       end do
    end do

    do    j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          uf2(i,j) = -0.0625d0*(uc2(i,j-2)+uc2(i,j+1)) + 0.5625d0*(uc2(i,j-1)+uc2(i,j))
       end do
    end do

  end subroutine cc2yface_2d


  subroutine cc2DyXface_2d(lo,hi,u,ulo,uhi,uf1,uf2,flo,fhi,uc1,uc2,clo,chi)
    integer, intent(in) :: lo(2),hi(2),ulo(2),uhi(2),flo(2),fhi(2),clo(2),chi(2)
    double precision, intent(in ) :: u  (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision, intent(out) :: uf1(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: uf2(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: uc1(clo(1):chi(1),clo(2):chi(2))
    double precision, intent(out) :: uc2(clo(1):chi(1),clo(2):chi(2))

    integer :: i, j

    do j   =lo(2)  ,hi(2)
       do i=lo(1)-2,hi(1)+2
          uc1(i,j) = dg1(-2)*u(i,j-2)+dg1(-1)*u(i,j-1)+dg1(0)*u(i,j) &
               +     dg1( 1)*u(i,j+1)+dg1( 2)*u(i,j+2)
          uc2(i,j) = dg2(-2)*u(i,j-2)+dg2(-1)*u(i,j-1)+dg2(0)*u(i,j) &
               +     dg2( 1)*u(i,j+1)+dg2( 2)*u(i,j+2)          
       end do
    end do

    do    j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          uf1(i,j) = -0.0625d0*(uc1(i-2,j)+uc1(i+1,j)) + 0.5625d0*(uc1(i-1,j)+uc1(i,j))
          uf2(i,j) = -0.0625d0*(uc2(i-2,j)+uc2(i+1,j)) + 0.5625d0*(uc2(i-1,j)+uc2(i,j))
       end do
    end do

  end subroutine cc2DyXface_2d


  subroutine cc2DxYface_2d(lo,hi,u,ulo,uhi,uf1,uf2,flo,fhi,uc1,uc2,clo,chi)
    integer, intent(in) :: lo(2),hi(2),ulo(2),uhi(2),flo(2),fhi(2),clo(2),chi(2)
    double precision, intent(in ) :: u  (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision, intent(out) :: uf1(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: uf2(flo(1):fhi(1),flo(2):fhi(2))
    double precision, intent(out) :: uc1(clo(1):chi(1),clo(2):chi(2))
    double precision, intent(out) :: uc2(clo(1):chi(1),clo(2):chi(2))

    integer :: i, j

    do j   =lo(2)-2,hi(2)+2
       do i=lo(1)  ,hi(1)
          uc1(i,j) = dg1(-2)*u(i-2,j)+dg1(-1)*u(i-1,j)+dg1(0)*u(i,j) &
               +     dg1( 1)*u(i+1,j)+dg1( 2)*u(i+2,j)
          uc2(i,j) = dg2(-2)*u(i-2,j)+dg2(-1)*u(i-1,j)+dg2(0)*u(i,j) &
               +     dg2( 1)*u(i+1,j)+dg2( 2)*u(i+2,j)          
       end do
    end do

    do    j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          uf1(i,j) = -0.0625d0*(uc1(i,j-2)+uc1(i,j+1)) + 0.5625d0*(uc1(i,j-1)+uc1(i,j))
       end do
    end do

    do    j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          uf2(i,j) = -0.0625d0*(uc2(i,j-2)+uc2(i,j+1)) + 0.5625d0*(uc2(i,j-1)+uc2(i,j))
       end do
    end do

  end subroutine cc2DxYface_2d


  subroutine cc2zface_3d(lo,hi, u, ulo, uhi, uf, flo,fhi)
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3), flo(3), fhi(3)
    double precision, intent(in ) :: u (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision, intent(out) :: uf(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i, j, k

    do k      =lo(3),hi(3)+1
       do j   =lo(2),hi(2)
          do i=lo(1),hi(1)
             uf(i,j,k) = -0.0625d0*(u(i,j,k-2)+u(i,j,k+1)) + 0.5625d0*(u(i,j,k-1)+u(i,j,k))
          end do
       end do
    end do
       
  end subroutine cc2zface_3d


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


  subroutine cc2DzGauss_3d(lo,hi, uc, uclo, uchi, u1, u2, uzlo,uzhi)
    integer, intent(in) :: lo(3), hi(3), uclo(3), uchi(3), uzlo(3), uzhi(3)
    double precision, intent(in ) :: uc(uclo(1):uchi(1),uclo(2):uchi(2),uclo(3):uchi(3))
    double precision, intent(out) :: u1(uzlo(1):uzhi(1),uzlo(2):uzhi(2),uzlo(3):uzhi(3))
    double precision, intent(out) :: u2(uzlo(1):uzhi(1),uzlo(2):uzhi(2),uzlo(3):uzhi(3))

    integer :: i, j, k

    do k      =lo(3),hi(3)
       do j   =lo(2),hi(2)
          do i=lo(1),hi(1)
             u1(i,j,k) = dg1(-2)*uc(i,j,k-2)+dg1(-1)*uc(i,j,k-1)+dg1(0)*uc(i,j,k) &
                  +      dg1( 1)*uc(i,j,k+1)+dg1( 2)*uc(i,j,k+2)
             u2(i,j,k) = dg2(-2)*uc(i,j,k-2)+dg2(-1)*uc(i,j,k-1)+dg2(0)*uc(i,j,k) &
                  +      dg2( 1)*uc(i,j,k+1)+dg2( 2)*uc(i,j,k+2)
          end do
       end do
    end do

  end subroutine cc2DzGauss_3d
  

  subroutine cc2xyGauss_3d(lo,hi,u,ulo,uhi,u1,u2,u3,u4,glo,ghi)
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3), glo(3), ghi(3)
    double precision, intent(in ) :: u (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision, intent(out) :: u1(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))
    double precision, intent(out) :: u2(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))
    double precision, intent(out) :: u3(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))
    double precision, intent(out) :: u4(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))

    integer :: i,j,k
    double precision :: t1(lo(1)-2:hi(1)+2), t2(lo(1)-2:hi(1)+2)

    do k      =lo(3),hi(3)
       do j   =lo(2),hi(2)

          do i=lo(1)-2,hi(1)+2
             t1(i) = cg1(-2)*u(i,j-2,k)+cg1(-1)*u(i,j-1,k)+cg1(0)*u(i,j,k) &
                  +  cg1( 1)*u(i,j+1,k)+cg1( 2)*u(i,j+2,k)
             t2(i) = cg2(-2)*u(i,j-2,k)+cg2(-1)*u(i,j-1,k)+cg2(0)*u(i,j,k) &
                  +  cg2( 1)*u(i,j+1,k)+cg2( 2)*u(i,j+2,k)
          end do

          do i=lo(1),hi(1)
             u1(i,j,k) = cg1(-2)*t1(i-2)+cg1(-1)*t1(i-1)+cg1(0)*t1(i) &
                  +      cg1( 1)*t1(i+1)+cg1( 2)*t1(i+2)
             u2(i,j,k) = cg2(-2)*t1(i-2)+cg2(-1)*t1(i-1)+cg2(0)*t1(i) &
                  +      cg2( 1)*t1(i+1)+cg2( 2)*t1(i+2)
          end do

          do i=lo(1),hi(1)
             u3(i,j,k) = cg1(-2)*t2(i-2)+cg1(-1)*t2(i-1)+cg1(0)*t2(i) &
                  +      cg1( 1)*t2(i+1)+cg1( 2)*t2(i+2)
             u4(i,j,k) = cg2(-2)*t2(i-2)+cg2(-1)*t2(i-1)+cg2(0)*t2(i) &
                  +      cg2( 1)*t2(i+1)+cg2( 2)*t2(i+2)
          end do

       end do
    end do

  end subroutine cc2xyGauss_3d


  subroutine cc2xyGaussDx_3d(lo,hi,u,ulo,uhi,ux1,ux2,ux3,ux4,xlo,xhi)
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),xlo(3),xhi(3)
    double precision, intent(in ) :: u  (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision, intent(out) :: ux1(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    double precision, intent(out) :: ux2(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    double precision, intent(out) :: ux3(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    double precision, intent(out) :: ux4(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))

    integer :: i, j, k
    double precision :: t1(lo(1)-2:hi(1)+2), t2(lo(1)-2:hi(1)+2)

    do k      =lo(3),hi(3)
       do j   =lo(2),hi(2)

          do i=lo(1)-2,hi(1)-2
             t1(i) = cg1(-2)*u(i,j-2,k)+cg1(-1)*u(i,j-1,k)+cg1(0)*u(i,j,k) &
                  +  cg1( 1)*u(i,j+1,k)+cg1( 2)*u(i,j+2,k)
             t2(i) = cg2(-2)*u(i,j-2,k)+cg2(-1)*u(i,j-1,k)+cg2(0)*u(i,j,k) &
                  +  cg2( 1)*u(i,j+1,k)+cg2( 2)*u(i,j+2,k)
          end do

          do i=lo(1),hi(1)
             ux1(i,j,k) = dg1(-2)*t1(i-2)+dg1(-1)*t1(i-1)+dg1(0)*t1(i) &
                  +       dg1( 1)*t1(i+1)+dg1( 2)*t1(i+2)
             ux2(i,j,k) = dg2(-2)*t1(i-2)+dg2(-1)*t1(i-1)+dg2(0)*t1(i) &
                  +       dg2( 1)*t1(i+1)+dg2( 2)*t1(i+2)
          end do

          do i=lo(1),hi(1)
             ux3(i,j,k) = dg1(-2)*t2(i-2)+dg1(-1)*t2(i-1)+dg1(0)*t2(i) &
                  +       dg1( 1)*t2(i+1)+dg1( 2)*t2(i+2)
             ux4(i,j,k) = dg2(-2)*t2(i-2)+dg2(-1)*t2(i-1)+dg2(0)*t2(i) &
                  +       dg2( 1)*t2(i+1)+dg2( 2)*t2(i+2)
          end do

       end do
    end do

  end subroutine cc2xyGaussDx_3d

  subroutine cc2xyGaussDy_3d(lo,hi,u,ulo,uhi,uy1,uy2,uy3,uy4,ylo,yhi)
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),ylo(3),yhi(3)
    double precision, intent(in ) :: u  (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision, intent(out) :: uy1(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    double precision, intent(out) :: uy2(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    double precision, intent(out) :: uy3(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    double precision, intent(out) :: uy4(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))

    integer :: i, j, k
    double precision :: t1(lo(1)-2:hi(1)+2), t2(lo(1)-2:hi(1)+2)

    do k      =lo(3),hi(3)
       do j   =lo(2),hi(2)

          do i=lo(1)-2,hi(1)-2
             t1(i) = dg1(-2)*u(i,j-2,k)+dg1(-1)*u(i,j-1,k)+dg1(0)*u(i,j,k) &
                  +  dg1( 1)*u(i,j+1,k)+dg1( 2)*u(i,j+2,k)
             t2(i) = dg2(-2)*u(i,j-2,k)+dg2(-1)*u(i,j-1,k)+dg2(0)*u(i,j,k) &
                  +  dg2( 1)*u(i,j+1,k)+dg2( 2)*u(i,j+2,k)
          end do

          do i=lo(1),hi(1)
             uy1(i,j,k) = cg1(-2)*t1(i-2)+cg1(-1)*t1(i-1)+cg1(0)*t1(i) &
                  +       cg1( 1)*t1(i+1)+cg1( 2)*t1(i+2)
             uy2(i,j,k) = cg2(-2)*t1(i-2)+cg2(-1)*t1(i-1)+cg2(0)*t1(i) &
                  +       cg2( 1)*t1(i+1)+cg2( 2)*t1(i+2)
          end do

          do i=lo(1),hi(1)
             uy3(i,j,k) = cg1(-2)*t2(i-2)+cg1(-1)*t2(i-1)+cg1(0)*t2(i) &
                  +       cg1( 1)*t2(i+1)+cg1( 2)*t2(i+2)
             uy4(i,j,k) = cg2(-2)*t2(i-2)+cg2(-1)*t2(i-1)+cg2(0)*t2(i) &
                  +       cg2( 1)*t2(i+1)+cg2( 2)*t2(i+2)
          end do

       end do
    end do

  end subroutine cc2xyGaussDy_3d

end module polyinterp_module
