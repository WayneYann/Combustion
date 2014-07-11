module weno_module

  implicit none

  double precision, parameter :: epsw = 1.d-6, b1=13.d0/12.d0, oneSixth=1.d0/6.d0

  ! Three types of third-order coefficients for converting cell averages to cell center point values
  double precision, dimension(-2:0), parameter :: L3_cc = &
       (/ -1.d0/24.d0,  2.d0/24.d0, 23.d0/24.d0 /)
  double precision, dimension(-1:1), parameter :: C3_cc = &  ! this is actually 4th-order
       (/ -1.d0/24.d0, 26.d0/24.d0, -1.d0/24.d0 /)
  double precision, dimension( 0:2), parameter :: R3_cc = &
       (/ 23.d0/24.d0,  2.d0/24.d0, -1.d0/24.d0 /)
  
  ! Optimal weights for cc
  double precision, parameter :: gamma_p_lr = 9.d0/214.d0
  double precision, parameter :: gamma_p_cc = 98.d0/107.d0
  double precision, parameter :: gamma_m_lr = 9.d0/67.d0
  double precision, parameter :: gamma_m_cc = 49.d0/67.d0
  double precision, parameter :: sigma_p = 107.d0/40.d0
  double precision, parameter :: sigma_m = 67.d0/40.d0

  ! Three types of third-order coefficients for converting cell averages to first Gauss point values
  double precision, dimension(-2:0), parameter :: L3_cg1 = &
       (/  -sqrt(3.d0)/12.d0,  1.d0/sqrt(3.d0),  1.d0-sqrt(3.d0)/4.d0  /)
  double precision, dimension(-1:1), parameter :: C3_cg1 = &
       (/   sqrt(3.d0)/12.d0,  1.d0,  -sqrt(3.d0)/12.d0  /)
  double precision, dimension( 0:2), parameter :: R3_cg1 = &
       (/  1.d0+sqrt(3.d0)/4.d0,  -1.d0/sqrt(3.d0),  sqrt(3.d0)/12.d0  /)

  ! Three types of third-order coefficients for converting cell averages to second Gauss point values
  double precision, dimension(-2:0), parameter :: L3_cg2 = &
       (/  sqrt(3.d0)/12.d0,  -1.d0/sqrt(3.d0),  1.d0+sqrt(3.d0)/4.d0  /)
  double precision, dimension(-1:1), parameter :: C3_cg2 = &
       (/   -sqrt(3.d0)/12.d0,  1.d0,  sqrt(3.d0)/12.d0  /)
  double precision, dimension( 0:2), parameter :: R3_cg2 = &
       (/  1.d0-sqrt(3.d0)/4.d0,  1.d0/sqrt(3.d0),  -sqrt(3.d0)/12.d0  /)

  ! Optimal weights for first Gauss point
  double precision, dimension(-2:0), parameter :: d_g1 = &
       (/ (210.d0+sqrt(3.d0))/1080.d0,  11.d0/18.d0,  (210.d0-sqrt(3.d0))/1080.d0  /)

  ! Optimal weights for second Gauss point
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

  ! given cell averages, compute derivatives at two Gauss point
  double precision, dimension(-2:2), parameter :: derg1 = &
       (/ (9.d0+2.d0*sqrt(3.d0))/108.d0,  -(36.d0+13.d0*sqrt(3.d0))/54.d0,  &
       12.d0*sqrt(3.d0)/27.d0,  (36.d0-13.d0*sqrt(3.d0))/54.d0,  (-9.d0+2.d0*sqrt(3.d0))/108.d0 /)
  double precision, dimension(-2:2), parameter :: derg2 = &
       (/ (9.d0-2.d0*sqrt(3.d0))/108.d0,  (-36.d0+13.d0*sqrt(3.d0))/54.d0,  &
       -12.d0*sqrt(3.d0)/27.d0,  (36.d0+13.d0*sqrt(3.d0))/54.d0,  -(9.d0+2.d0*sqrt(3.d0))/108.d0 /)


  interface cellavg2gausspt_2d
     module procedure cellavg2gausspt_2d_v1
     module procedure cellavg2gausspt_2d_v2
  end interface


  private

  public :: weno5, vweno5, weno5_center,  &
       cellavg2gausspt_1d, cellavg2gausspt_2d, cellavg2gausspt_2d_v1, &
       cellavg2gausspt_2d_v2, cellavg2gausspt_3d, &
       cellavg2dergausspt_1d, cellavg2dergausspt_2d, cellavg2face_1d

contains

  subroutine weno5(v, vp, vm, vg1, vg2)
    double precision, intent(in)  :: v(-2:2)
    double precision, intent(out), optional :: vp , vm   ! v_{i+1/2} & v_{i-1/2}
    double precision, intent(out), optional :: vg1, vg2  ! at two Gauss points

    double precision :: vr_2, vr_1, vr_0
    double precision :: beta_2, beta_1, beta_0
    double precision :: alpha_2, alpha_1, alpha_0
    double precision :: alpha1

    beta_2 = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta_2 = 1.d0/(epsw+beta_2)**2

    beta_1 = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta_1 = 1.d0/(epsw+beta_1)**2

    beta_0 = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
    beta_0 = 1.d0/(epsw+beta_0)**2

    if (present(vp)) then
       alpha_2 =      beta_2
       alpha_1 = 6.d0*beta_1
       alpha_0 = 3.d0*beta_0
       alpha1 = 1.d0/(alpha_2 + alpha_1 + alpha_0)

       vr_2 = 2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
       vr_1 =     -v(-1) + 5.d0*v( 0) +  2.d0*v(1)
       vr_0 = 2.d0*v( 0) + 5.d0*v( 1) -       v(2)
       
       vp = oneSixth*alpha1*(alpha_2*vr_2 + alpha_1*vr_1 + alpha_0*vr_0)
    end if

    if (present(vm)) then
       alpha_2 = 3.d0*beta_2
       alpha_1 = 6.d0*beta_1
       alpha_0 =      beta_0
       alpha1 = 1.d0/(alpha_2 + alpha_1 + alpha_0)
       
       vr_2 =      -v(-2) + 5.d0*v(-1) + 2.d0*v(0)
       vr_1 =  2.d0*v(-1) + 5.d0*v(0 ) -      v(1) 
       vr_0 = 11.d0*v( 0) - 7.d0*v(1 ) + 2.d0*v(2)
       
       vm = oneSixth*alpha1*(alpha_2*vr_2 + alpha_1*vr_1 + alpha_0*vr_0)
    end if

    if (present(vg1)) then
       alpha_2 = d_g1(-2)*beta_2
       alpha_1 = d_g1(-1)*beta_1
       alpha_0 = d_g1( 0)*beta_0
       alpha1 = 1.d0/(alpha_2 + alpha_1 + alpha_0)
       
       vr_2 = L3_cg1(-2)*v(-2) + L3_cg1(-1)*v(-1) + L3_cg1(0)*v(0)
       vr_1 = C3_cg1(-1)*v(-1) + C3_cg1( 0)*v( 0) + C3_cg1(1)*v(1)
       vr_0 = R3_cg1( 0)*v( 0) + R3_cg1( 1)*v( 1) + R3_cg1(2)*v(2)

       vg1 = alpha1*(alpha_2*vr_2 + alpha_1*vr_1 + alpha_0*vr_0)
    end if

    if (present(vg2)) then
       alpha_2 = d_g2(-2)*beta_2
       alpha_1 = d_g2(-1)*beta_1
       alpha_0 = d_g2( 0)*beta_0
       alpha1 = 1.d0/(alpha_2 + alpha_1 + alpha_0)
       
       vr_2 = L3_cg2(-2)*v(-2) + L3_cg2(-1)*v(-1) + L3_cg2(0)*v(0)
       vr_1 = C3_cg2(-1)*v(-1) + C3_cg2( 0)*v( 0) + C3_cg2(1)*v(1)
       vr_0 = R3_cg2( 0)*v( 0) + R3_cg2( 1)*v( 1) + R3_cg2(2)*v(2)

       vg2 = alpha1*(alpha_2*vr_2 + alpha_1*vr_1 + alpha_0*vr_0)
    end if

    return
  end subroutine weno5


  subroutine vweno5(lo, hi, v, vlo, vhi, glo, ghi, vp, vm, vg1, vg2)
    integer, intent(in) :: lo, hi, vlo, vhi, glo, ghi
    double precision, intent(in)  :: v(vlo:vhi)
    double precision, dimension(glo-1:ghi+1), intent(out), optional :: vp, vm  ! v_{i+1/2} & v_{i-1/2}
    double precision, dimension(glo  :ghi  ), intent(out), optional :: vg1, vg2  ! at two Gauss points

    integer :: i
    double precision, dimension(lo:hi) :: vr_2, vr_1, vr_0
    double precision, dimension(lo:hi) :: alpha_2, alpha_1, alpha_0
    double precision, dimension(lo:hi) :: beta_2, beta_1, beta_0

    !DEC$ SIMD
    do i=lo,hi
       beta_2(i) = b1*(v(i-2)-2.d0*v(i-1)+v(i))**2 + 0.25d0*(v(i-2)-4.d0*v(i-1)+3.d0*v(i))**2
       beta_2(i) = 1.d0/(epsw+beta_2(i))**2
       
       beta_1(i) = b1*(v(i-1)-2.d0*v(i  )+v(i+1))**2 + 0.25d0*(v(i-1)-v(i+1))**2
       beta_1(i) = 1.d0/(epsw+beta_1(i))**2
       
       beta_0(i) = b1*(v(i  )-2.d0*v(i+1)+v(i+2))**2 + 0.25d0*(3.d0*v(i)-4.d0*v(i+1)+v(i+2))**2
       beta_0(i) = 1.d0/(epsw+beta_0(i))**2
    end do

    if (present(vp)) then
       !DEC$ SIMD
       do i=lo,hi
          alpha_2(i) =      beta_2(i)
          alpha_1(i) = 6.d0*beta_1(i)
          alpha_0(i) = 3.d0*beta_0(i)
          
          vr_2(i) = 2.d0*v(i-2) - 7.d0*v(i-1) + 11.d0*v(i  )
          vr_1(i) =     -v(i-1) + 5.d0*v(i  ) +  2.d0*v(i+1)
          vr_0(i) = 2.d0*v(i  ) + 5.d0*v(i+1) -       v(i+2)

          vp(i) = oneSixth*(alpha_2(i)*vr_2(i) + alpha_1(i)*vr_1(i) + alpha_0(i)*vr_0(i)) &
               / (alpha_2(i) + alpha_1(i) + alpha_0(i))

          alpha_2(i) = 3.d0*beta_2(i)
          alpha_1(i) = 6.d0*beta_1(i)
          alpha_0(i) =      beta_0(i)
          
          vr_2(i) =      -v(i-2) + 5.d0*v(i-1) + 2.d0*v(i  )
          vr_1(i) =  2.d0*v(i-1) + 5.d0*v(i  ) -      v(i+1) 
          vr_0(i) = 11.d0*v(i  ) - 7.d0*v(i+1) + 2.d0*v(i+2)
          
          vm(i) = oneSixth*(alpha_2(i)*vr_2(i) + alpha_1(i)*vr_1(i) + alpha_0(i)*vr_0(i)) &
               / (alpha_2(i) + alpha_1(i) + alpha_0(i))
       end do
    endif

    if (present(vg1)) then
       !DEC$ SIMD
       do i=glo,ghi
          alpha_2(i) = d_g1(-2)*beta_2(i)
          alpha_1(i) = d_g1(-1)*beta_1(i)
          alpha_0(i) = d_g1( 0)*beta_0(i)
          
          vr_2(i) = L3_cg1(-2)*v(i-2) + L3_cg1(-1)*v(i-1) + L3_cg1(0)*v(i  )
          vr_1(i) = C3_cg1(-1)*v(i-1) + C3_cg1( 0)*v(i  ) + C3_cg1(1)*v(i+1)
          vr_0(i) = R3_cg1( 0)*v(i  ) + R3_cg1( 1)*v(i+1) + R3_cg1(2)*v(i+2)
          
          vg1(i) = (alpha_2(i)*vr_2(i) + alpha_1(i)*vr_1(i) + alpha_0(i)*vr_0(i)) &
               / (alpha_2(i) + alpha_1(i) + alpha_0(i))

          alpha_2(i) = d_g2(-2)*beta_2(i)
          alpha_1(i) = d_g2(-1)*beta_1(i)
          alpha_0(i) = d_g2( 0)*beta_0(i)
       
          vr_2(i) = L3_cg2(-2)*v(i-2) + L3_cg2(-1)*v(i-1) + L3_cg2(0)*v(i  )
          vr_1(i) = C3_cg2(-1)*v(i-1) + C3_cg2( 0)*v(i  ) + C3_cg2(1)*v(i+1)
          vr_0(i) = R3_cg2( 0)*v(i  ) + R3_cg2( 1)*v(i+1) + R3_cg2(2)*v(i+2)

          vg2(i) = (alpha_2(i)*vr_2(i) + alpha_1(i)*vr_1(i) + alpha_0(i)*vr_0(i)) &
               / (alpha_2(i) + alpha_1(i) + alpha_0(i))
       end do

    end if

  end subroutine vweno5


  subroutine weno5_center(v, vc)
    double precision, intent(in) :: v(-2:2)
    double precision, intent(out) :: vc  ! at cell center

    double precision :: alpha_2, alpha_1, alpha_0
    double precision :: beta_2, beta_1, beta_0
    double precision :: v_2, v_1, v_0

    beta_2 = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta_2 = 1.d0/(epsw+beta_2)**2

    beta_1 = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta_1 = 1.d0/(epsw+beta_1)**2

    beta_0 = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
    beta_0 = 1.d0/(epsw+beta_0)**2

    v_2 = L3_cc(-2)*v(-2) + L3_cc(-1)*v(-1) + L3_cc(0)*v(0)
    v_1 = C3_cc(-1)*v(-1) + C3_cc( 0)*v( 0) + C3_cc(1)*v(1)
    v_0 = R3_cc( 0)*v( 0) + R3_cc( 1)*v( 1) + R3_cc(2)*v(2)

    ! plus
    alpha_2 = gamma_p_lr*beta_2
    alpha_1 = gamma_p_cc*beta_1
    alpha_0 = gamma_p_lr*beta_0
    vc = sigma_p*(alpha_2*v_2 + alpha_1*v_1 + alpha_0*v_0)/(alpha_2 + alpha_1 + alpha_0)

    ! minus
    alpha_2 = gamma_m_lr*beta_2
    alpha_1 = gamma_m_cc*beta_1
    alpha_0 = gamma_m_lr*beta_0
    vc = vc-sigma_m*(alpha_2*v_2 + alpha_1*v_1 + alpha_0*v_0)/(alpha_2 + alpha_1 + alpha_0)
  end subroutine weno5_center


  subroutine cellavg2gausspt_1d(lo,hi, u, ulo,uhi, u1, u2, glo,ghi)
    integer, intent(in) :: lo, hi, ulo, uhi, glo, ghi
    double precision, intent(in)  :: u(ulo:uhi)
    double precision, intent(out) :: u1(glo:ghi), u2(glo:ghi)
    integer :: i
    do i=lo,hi
       u1(i) = cg1(-2)*u(i-2) + cg1(-1)*u(i-1) + cg1(0)*u(i) + cg1(1)*u(i+1) + cg1(2)*u(i+2)
       u2(i) = cg2(-2)*u(i-2) + cg2(-1)*u(i-1) + cg2(0)*u(i) + cg2(1)*u(i+1) + cg2(2)*u(i+2)
    end do
  end subroutine cellavg2gausspt_1d

  subroutine cellavg2gausspt_2d_v1(lo, hi, u, ulo, uhi, ug, glo, ghi)
    integer, intent(in) :: lo(2), hi(2), ulo(2), uhi(2), glo(2), ghi(2)
    double precision, intent(in)  :: u (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision, intent(out) :: ug(glo(1):ghi(1),glo(2):ghi(2),4)

    integer :: i, j, g, gg
    double precision, allocatable :: ugy(:,:,:)
    
    allocate(ugy(lo(1)-2:hi(1)+2,lo(2):hi(2),2))

    do i=lo(1)-2,hi(1)+2
       call cellavg2gausspt_1d(lo(2),hi(2), u(i,:), ulo(2),uhi(2), &
            ugy(i,:,1), ugy(i,:,2), lo(2),hi(2))
    end do

    do g=1,2
       gg = 2*(g-1)
       do j=lo(2),hi(2)
          call cellavg2gausspt_1d(lo(1),hi(1), ugy(:,j,g), lo(1)-2,hi(1)+2, &
               ug(:,j,gg+1), ug(:,j,gg+2), glo(1),ghi(1))
       end do
    end do

    deallocate(ugy)

  end subroutine cellavg2gausspt_2d_v1

  subroutine cellavg2gausspt_2d_v2(lo, hi, u, ulo, uhi, u1, u2, u3, u4, glo, ghi)
    integer, intent(in) :: lo(2), hi(2), ulo(2), uhi(2), glo(2), ghi(2)
    double precision,intent(in)  :: u (ulo(1):uhi(1),ulo(2):uhi(2))
    double precision,intent(out) :: u1(glo(1):ghi(1),glo(2):ghi(2))
    double precision,intent(out) :: u2(glo(1):ghi(1),glo(2):ghi(2))
    double precision,intent(out) :: u3(glo(1):ghi(1),glo(2):ghi(2))
    double precision,intent(out) :: u4(glo(1):ghi(1),glo(2):ghi(2))

    integer :: i, j
    double precision, allocatable :: ugy(:,:,:)
    
    allocate(ugy(lo(1)-2:hi(1)+2,lo(2):hi(2),2))

    do i=lo(1)-2,hi(1)+2
       call cellavg2gausspt_1d(lo(2),hi(2), u(i,:), ulo(2),uhi(2), &
            ugy(i,:,1), ugy(i,:,2), lo(2),hi(2))
    end do

    do j=lo(2),hi(2)
       call cellavg2gausspt_1d(lo(1),hi(1), ugy(:,j,1), lo(1)-2,hi(1)+2, &
            u1(:,j), u2(:,j), glo(1),ghi(1))
    end do

    do j=lo(2),hi(2)
       call cellavg2gausspt_1d(lo(1),hi(1), ugy(:,j,2), lo(1)-2,hi(1)+2, &
            u3(:,j), u4(:,j), glo(1),ghi(1))
    end do

    deallocate(ugy)

  end subroutine cellavg2gausspt_2d_v2


  subroutine cellavg2gausspt_3d(lo, hi, u, ulo, uhi, ug, glo, ghi)
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3), glo(3), ghi(3)
    double precision, intent(in)  :: u (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision, intent(out) :: ug(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),8)

    integer :: i, j, k, tlo(2), thi(2)
    double precision, allocatable :: ugz(:,:,:,:)
    
    allocate(ugz(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3):hi(3),2))

    do j=lo(2)-2,hi(2)+2
    do i=lo(1)-2,hi(1)+2
       call cellavg2gausspt_1d(lo(3),hi(3), u(i,j,:), ulo(3),uhi(3), &
            ugz(i,j,:,1), ugz(i,j,:,2), lo(3),hi(3))
    end do
    end do

    tlo(1) = lo(1)-2
    tlo(2) = lo(2)-2
    thi(1) = hi(1)+2
    thi(2) = hi(2)+2
    do k=lo(3),hi(3)
       call cellavg2gausspt_2d(lo(1:2), hi(1:2), ugz(:,:,k,1), tlo, thi, &
            ug(:,:,k,1:4), glo(1:2), ghi(1:2))
       call cellavg2gausspt_2d(lo(1:2), hi(1:2), ugz(:,:,k,2), tlo, thi, &
            ug(:,:,k,5:8), glo(1:2), ghi(1:2))
    end do

    deallocate(ugz)

  end subroutine cellavg2gausspt_3d


  subroutine cellavg2dergausspt_1d(lo, hi, u, ulo, uhi, du1, du2, glo, ghi)
    integer, intent(in) :: lo, hi, ulo, uhi, glo, ghi
    double precision, intent(in)  :: u(ulo:uhi)
    double precision, intent(out) :: du1(glo:ghi), du2(glo:ghi)
    integer :: i
    do i=lo,hi
       du1(i) = derg1(-2)*u(i-2) + derg1(-1)*u(i-1) + derg1(0)*u(i) + derg1(1)*u(i+1) + derg1(2)*u(i+2)
       du2(i) = derg2(-2)*u(i-2) + derg2(-1)*u(i-1) + derg2(0)*u(i) + derg2(1)*u(i+1) + derg2(2)*u(i+2)
    end do
  end subroutine cellavg2dergausspt_1d


  subroutine cellavg2dergausspt_2d(lo, hi, u, ulo, uhi, du1, du2, du3, du4, dulo, duhi, der_dir)
    integer, intent(in) :: lo(2), hi(2), ulo(2), uhi(2),dulo(2),duhi(2), der_dir
    double precision,intent(in)  ::   u( ulo(1): uhi(1), ulo(2): uhi(2))
    double precision,intent(out) :: du1(dulo(1):duhi(1),dulo(2):duhi(2))
    double precision,intent(out) :: du2(dulo(1):duhi(1),dulo(2):duhi(2))
    double precision,intent(out) :: du3(dulo(1):duhi(1),dulo(2):duhi(2))
    double precision,intent(out) :: du4(dulo(1):duhi(1),dulo(2):duhi(2))

    integer :: i, j
    double precision, allocatable :: ugy(:,:,:)

    allocate(ugy(lo(1)-2:hi(1)+2,lo(2):hi(2),2))

    if (der_dir .eq. 1) then 

       do i=lo(1)-2,hi(1)+2
          call cellavg2gausspt_1d(lo(2),hi(2), u(i,:), ulo(2),uhi(2), &
               ugy(i,:,1), ugy(i,:,2), lo(2),hi(2))
       end do

       do j=lo(2),hi(2)      ! d./dx
          call cellavg2dergausspt_1d(lo(1),hi(1), ugy(:,j,1), lo(1)-2,hi(1)+2, &
               du1(:,j), du2(:,j), dulo(1),duhi(1))
       end do
       do j=lo(2),hi(2)      ! d./dx
          call cellavg2dergausspt_1d(lo(1),hi(1), ugy(:,j,2), lo(1)-2,hi(1)+2, &
               du3(:,j), du4(:,j), dulo(1),duhi(1))
       end do
       
    else                     ! d./dy

       do i=lo(1)-2,hi(1)+2
          call cellavg2dergausspt_1d(lo(2),hi(2), u(i,:), ulo(2),uhi(2), &
               ugy(i,:,1), ugy(i,:,2), lo(2),hi(2))
       end do       

       do j=lo(2),hi(2)
          call cellavg2gausspt_1d(lo(1),hi(1), ugy(:,j,1), lo(1)-2,hi(1)+2, &
               du1(:,j), du2(:,j), dulo(1),duhi(1))
       end do
       do j=lo(2),hi(2)
          call cellavg2gausspt_1d(lo(1),hi(1), ugy(:,j,2), lo(1)-2,hi(1)+2, &
               du3(:,j), du4(:,j), dulo(1),duhi(1))
       end do

    end if

    deallocate(ugy)

  end subroutine cellavg2dergausspt_2d


  subroutine cellavg2face_1d(lo, hi, u, ulo, uhi, uf, flo, fhi)
    integer, intent(in) :: lo, hi, ulo, uhi, flo, fhi
    double precision, intent(in) ::  u(ulo:uhi)
    double precision, intent(out):: uf(flo:fhi)
    integer :: i
    do i=lo,hi
       uf(i) = cc4(-2)*u(i-2) + cc4(-1)*u(i-1) + cc4(0)*u(i) + cc4(1)*u(i+1)
    end do
  end subroutine cellavg2face_1d

end module weno_module

