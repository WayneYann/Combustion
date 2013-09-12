module weno_module

  implicit none

  double precision, parameter :: epsw = 1.d-6, b1=13.d0/12.d0, oneSixth=1.d0/6.d0

  ! Three types of third-order coefficients for converting cell averages to first Gauss point values
  double precision, dimension(-2:0), parameter :: L3_cg1 = &
       (/  -sqrt(3.d0)/12.d0,  1.d0/sqrt(3.d0),  1.d0-sqrt(3.d0)/4.d0  /)
  double precision, dimension(-1:1), parameter :: C3_cg1 = &
       (/   sqrt(3.d0)/12.d0,  1.d0,  -sqrt(3.d0)/12.d0  /)
  double precision, dimension( 0:2), parameter :: R3_cg1 = &
       (/  1+sqrt(3.d0)/4.d0,  -1.d0/sqrt(3.d0),  sqrt(3.d0)/12.d0  /)

  ! Three types of third-order coefficients for converting cell averages to second Gauss point values
  double precision, dimension(-2:0), parameter :: L3_cg2 = &
       (/  sqrt(3.d0)/12.d0,  -1.d0/sqrt(3.d0),  1.d0+sqrt(3.d0)/4.d0  /)
  double precision, dimension(-1:1), parameter :: C3_cg2 = &
       (/   -sqrt(3.d0)/12.d0,  1.d0,  sqrt(3.d0)/12.d0  /)
  double precision, dimension( 0:2), parameter :: R3_cg2 = &
       (/  1-sqrt(3.d0)/4.d0,  1.d0/sqrt(3.d0),  -sqrt(3.d0)/12.d0  /)

  ! Optimial weights for first Gauss point
  double precision, dimension(-2:0), parameter :: d_g1 = &
       (/ (210.d0+sqrt(3.d0))/1080.d0,  11.d0/18.d0,  (210.d0-sqrt(3.d0))/1080.d0  /)

  ! Optimial weights for second Gauss point
  double precision, dimension(-2:0), parameter :: d_g2 = &
       (/ (210.d0-sqrt(3.d0))/1080.d0,  11.d0/18.d0,  (210.d0+sqrt(3.d0))/1080.d0  /)

  ! Fifth-order coefficients for converting cell averages to two Gauss point values
  double precision, dimension(-2:2), parameter :: cg1 = &
       (/ -(1.d0+70.d0*sqrt(3.d0))/4320.d0, (4.d0+500.d0*sqrt(3.d0))/4320.d0, &
       4314.d0/4320.d0, (4.d0-500.d0*sqrt(3.d0))/4320.d0, (-1.d0+70.d0*sqrt(3.d0))/4320.d0 /)
  double precision, dimension(-2:2), parameter :: cg2 = &
       (/ (-1.d0+70.d0*sqrt(3.d0))/4320.d0, (4.d0-500.d0*sqrt(3.d0))/4320.d0, &
       4314.d0/4320.d0, (4.d0+500.d0*sqrt(3.d0))/4320.d0, -(1.d0+70.d0*sqrt(3.d0))/4320.d0 /)

  ! v_{i-1/2} = sum(cc4*v(i-2:i+1)) + O(h^4), where v(i-2:i+1) are cell averages
  double precision, dimension(-2:1), parameter :: cc4 = &
       (/  -1.d0/12.d0,  7.d0/12.d0,  7.d0/12.d0,  -1.d0/12.d0  /)

  private

  public :: weno5, cellavg2gausspt_1d, cellavg2face_1d, cellavg2gausspt_2d

contains

  subroutine weno5(v, vp, vm, vg1, vg2)
    double precision, intent(in)  :: v(-2:2)
    double precision, intent(out), optional :: vp , vm   ! v_{i+1/2} & v_{i-1/2}
    double precision, intent(out), optional :: vg1, vg2  ! at two Gauss points

    double precision :: vpr(-2:0), vmr(-2:0), v1r(-2:0), v2r(-2:0)
    double precision :: beta(-2:0), alpha(-2:0), alpha1

    beta(-2) = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta(-1) = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta( 0) = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2

    beta(-2) = 1.d0/(epsw+beta(-2))**2
    beta(-1) = 1.d0/(epsw+beta(-1))**2
    beta( 0) = 1.d0/(epsw+beta( 0))**2

    if (present(vp)) then
       alpha(-2) =      beta(-2)
       alpha(-1) = 6.d0*beta(-1)
       alpha( 0) = 3.d0*beta( 0)
       alpha1 = 1.d0/(alpha(-2) + alpha(-1) + alpha(0))

       vpr(-2) = 2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
       vpr(-1) =     -v(-1) + 5.d0*v( 0) +  2.d0*v(1)
       vpr( 0) = 2.d0*v( 0) + 5.d0*v( 1) -       v(2)
       
       vp = oneSixth*alpha1*(alpha(-2)*vpr(-2) + alpha(-1)*vpr(-1) + alpha(0)*vpr(0))
    end if

    if (present(vm)) then
       alpha(-2) = 3.d0*beta(-2)
       alpha(-1) = 6.d0*beta(-1)
       alpha( 0) =      beta( 0)
       alpha1 = 1.d0/(alpha(-2) + alpha(-1) + alpha(0))
       
       vmr(-2) =      -v(-2) + 5.d0*v(-1) + 2.d0*v(0)
       vmr(-1) =  2.d0*v(-1) + 5.d0*v(0 ) -      v(1) 
       vmr( 0) = 11.d0*v( 0) - 7.d0*v(1 ) + 2.d0*v(2)
       
       vm = oneSixth*alpha1*(alpha(-2)*vmr(-2) + alpha(-1)*vmr(-1) + alpha(0)*vmr(0))
    end if

    if (present(vg1)) then
       alpha(-2) = d_g1(-2)*beta(-2)
       alpha(-1) = d_g1(-1)*beta(-1)
       alpha( 0) = d_g1( 0)*beta( 0)
       alpha1 = 1.d0/(alpha(-2) + alpha(-1) + alpha(0))
       
       v1r(-2) = L3_cg1(-2)*v(-2) + L3_cg1(-1)*v(-1) + L3_cg1(0)*v(0)
       v1r(-1) = C3_cg1(-1)*v(-1) + C3_cg1( 0)*v( 0) + C3_cg1(1)*v(1)
       v1r( 0) = R3_cg1( 0)*v( 0) + R3_cg1( 1)*v( 1) + R3_cg1(2)*v(2)

       vg1 = alpha1*(alpha(-2)*v1r(-2) + alpha(-1)*v1r(-1) + alpha(0)*v1r(0))
    end if

    if (present(vg2)) then
       alpha(-2) = d_g2(-2)*beta(-2)
       alpha(-1) = d_g2(-1)*beta(-1)
       alpha( 0) = d_g2( 0)*beta( 0)
       alpha1 = 1.d0/(alpha(-2) + alpha(-1) + alpha(0))
       
       v2r(-2) = L3_cg2(-2)*v(-2) + L3_cg2(-1)*v(-1) + L3_cg2(0)*v(0)
       v2r(-1) = C3_cg2(-1)*v(-1) + C3_cg2( 0)*v( 0) + C3_cg2(1)*v(1)
       v2r( 0) = R3_cg2( 0)*v( 0) + R3_cg2( 1)*v( 1) + R3_cg2(2)*v(2)

       vg2 = alpha1*(alpha(-2)*v2r(-2) + alpha(-1)*v2r(-1) + alpha(0)*v2r(0))
    end if

    return
  end subroutine weno5


  subroutine cellavg2gausspt_1d(lo,hi, u, ulo,uhi, u1, u2, u12lo,u12hi)
    integer, intent(in) :: lo, hi, ulo, uhi, u12lo, u12hi
    double precision, intent(in) :: u(ulo:uhi)
    double precision :: u1(u12lo:u12hi), u2(u12lo:u12hi)

    integer :: i

    do i=lo,hi
       u1(i) = cg1(-2)*u(i-2) + cg1(-1)*u(i-1) + cg1(0)*u(i) + cg1(1)*u(i+1) + cg1(2)*u(i+2)
       u2(i) = cg2(-2)*u(i-2) + cg2(-1)*u(i-1) + cg2(0)*u(i) + cg2(1)*u(i+1) + cg2(2)*u(i+2)
    end do
  end subroutine cellavg2gausspt_1d

  subroutine cellavg2face_1d(lo, hi, u, ulo, uhi, uf, flo, fhi)
    integer, intent(in) :: lo, hi, ulo, uhi, flo, fhi
    double precision, intent(in) ::  u(ulo:uhi)
    double precision             :: uf(flo:fhi)
    integer :: i
    do i=lo,hi
       uf(i) = cc4(-2)*u(i-2) + cc4(-1)*u(i-1) + cc4(0)*u(i) + cc4(1)*u(i+1)
    end do
  end subroutine cellavg2face_1d


  subroutine cellavg2gausspt_2d(lo, hi, u, ulo, uhi, ug, glo, ghi)
    integer, intent(in) :: lo(2), hi(2), ulo(2), uhi(2), glo(2), ghi(2)
    double precision, intent(in) :: u (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision             :: ug(glo(1):ghi(1),glo(2):ghi(2),4)

    integer :: i, j, g, gg
    double precision, allocatable :: ugy(:,:,:)
    
    allocate(ugy(lo(1)-2:hi(1)+2,lo(2):hi(2),2))

    do i=lo(1)-2,hi(1)+2
       call cellavg2gausspt_1d(lo(2),hi(2), u(i,:), ulo(2),uhi(2), &
            ugy(i,:,1), ugy(i,:,2), ulo(2),uhi(2))
    end do

    do g=1,2
       gg = 2*(g-1)
       do j=lo(2),hi(2)
          call cellavg2gausspt_1d(lo(1),hi(1), ugy(:,j,g), lo(1)-2,hi(1)+2, &
               ug(:,j,gg+1), ug(:,j,gg+2), lo(1),hi(1))
       end do
    end do

    deallocate(ugy)

  end subroutine cellavg2gausspt_2d

end module weno_module
