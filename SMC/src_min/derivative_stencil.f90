module derivative_stencil_module

  implicit none

  public

  integer, parameter :: narrow=1, wide=2, S3D=2
  integer, save :: stencil, stencil_ng

  ! for 8th-order first derivatives
  double precision,dimension(4),parameter :: D8 = (/ 0.8d0, -0.2d0, 4.d0/105.d0, -1.d0/280.d0 /)

  ! 6th-order first derivative
  double precision,dimension(3),parameter :: D6 = (/ 0.75d0, -0.15d0, 1.d0/60.d0 /)

  ! 4th-order first derivative
  double precision,dimension(2),parameter :: D4 = (/ 2.d0/3.d0, -1.d0/12.d0 /)

  ! 3rd-order first derivative
  ! The stencil is slighly right-based with one cell on the left and two on the right
  double precision,dimension(-1:2),parameter :: DR3 = (/ -1.d0/3.d0, -0.5d0, 1.d0, -1.d0/6.d0 /)
  ! The stencil is slighly left-based with two cells on the left and on on the right
  double precision,dimension(-2:1),parameter :: DL3 = (/ 1.d0/6.d0, -1.d0, 0.5d0, 1.d0/3.d0 /)
  
  ! 3rd-order first derivative
  ! The stencil is right-based with no cell on the left and three on the right
  double precision,dimension(0:3),parameter :: DRB = (/ -11.d0/6.d0, 3.d0, -1.5d0, 1.d0/3.d0 /)
  ! The stencil is left-based with no cell on the right and three on the left
  double precision,dimension(-3:0),parameter :: DLB = (/ -1.d0/3.d0, 1.5d0, -3.d0, 11.d0/6.d0 /)

  ! coefficients for 8th-order stencil of second derivatives
  ! d(a*du/dx)/dx = H_{i+1/2} - H_{i-1/2},
  ! where H = a.M.u
  double precision, private, parameter :: M8_47 = 683.d0/10080.d0, M8_48 = -1.d0/224.d0
  double precision, save, dimension(8,8) :: M8

  ! coefficients for 6th-order stencil of second derivatives
  double precision, private, parameter :: M6_36 = 1.d0/90.d0
  double precision, save, dimension(6,6) :: M6

  ! coefficients for 4th-order stencil of second derivatives
  double precision, save, dimension(4,4) :: M4

  ! coefficients for 2nd-order stencil of second derivatives
  double precision, save, dimension(2,2) :: M2

  ! coefficients for 2nd-order biased stencil of second derivatives
  ! Given u and a on cell 1, 2, 3, 4, 
  !      the second derivate at cell 1 is give by a.B.u
  double precision, private, parameter :: BB2_33 = 2.d0
  double precision, save, dimension(4,4) :: BLB, BRB

contains
  
  subroutine stencil_init

    use bl_error_module
    use probin_module, only : stencil_type

    integer :: i, j

    if (trim(stencil_type) == "narrow") then
       stencil = narrow
    else if (trim(stencil_type) == "S3D" .or. trim(stencil_type) == "wide") then
       stencil = wide
    else
       call bl_error("unknow stencil_type")
    end if

    stencil_ng = 4

    ! 8th-order
    M8(1,1) = 5.d0/336.d0 + M8_48
    M8(2,1) = -11.d0/560.d0 - 2.d0*M8_48
    M8(3,1) = -1.d0/280.d0
    M8(4,1) = 17.d0/1680.d0 + 2.d0*M8_48
    M8(5,1) = -M8_48
    M8(6,1) = 0.d0
    M8(7,1) = 0.d0
    M8(8,1) = 0.d0

    M8(1,2) = -83.d0/3600.d0 - M8_47/5.d0 - 14.d0*M8_48/5.d0
    M8(2,2) = -31.d0/360.d0 + M8_47 + 3.d0*M8_48
    M8(3,2) = 1097.d0/5040.d0 - 2.d0*M8_47 + 6.d0*M8_48
    M8(4,2) = -319.d0/2520.d0 + 2.d0*M8_47 - 8.d0*M8_48
    M8(5,2) = -M8_47
    M8(6,2) = -139.d0/25200.d0 + M8_47/5.d0 + 9.d0*M8_48/5.d0
    M8(7,2) = 0.d0
    M8(8,2) = 0.d0

    M8(1,3) = 299.d0/50400.d0 + 2.d0*M8_47/5.d0 + 13.d0*M8_48/5.d0
    M8(2,3) = 41.d0/200.d0 - 9.d0*M8_47/5.d0 + 4.d0*M8_48/5.d0
    M8(3,3) = -1349.d0/10080.d0 + 3.d0*M8_47 - 12.d0*M8_48
    M8(4,3) = -919.d0/5040.d0 - 2.d0*M8_47 + 6.d0*M8_48
    M8(5,3) = 65.d0/224.d0 + 7.d0*M8_48
    M8(6,3) = -467.d0/25200.d0 + 3.d0*M8_47/5.d0 - 18.d0*M8_48/5.d0
    M8(7,3) = 503.d0/50400.d0 - M8_47/5.d0 - 4.d0*M8_48/5.d0
    M8(8,3) = 0.d0

    M8(1,4) = 17.d0/12600.d0 - M8_47/5.d0 - 4.d0*M8_48/5.d0
    M8(2,4) = -5927.d0/50400.d0 + 4.d0*M8_47/5.d0 - 9.d0*M8_48/5.d0
    M8(3,4) = -887.d0/5040.d0 - M8_47 + 6.d0*M8_48
    M8(4,4) = -445.d0/2016.d0
    M8(5,4) = -583.d0/720.d0 + M8_47 - 6.d0*M8_48
    M8(6,4) = -3613.d0/50400.d0 - 4.d0*M8_47/5.d0 + 9.d0*M8_48/5.d0
    M8(7,4) = -17.d0/600.d0 + M8_47/5.d0 + 4.d0*M8_48/5.d0
    M8(8,4) = -1.d0/1120.d0

    M8(1,5) = -M8(8,4)
    M8(2,5) = -M8(7,4)
    M8(3,5) = -M8(6,4)
    M8(4,5) = -M8(5,4)
    M8(5,5) = -M8(4,4)
    M8(6,5) = -M8(3,4)
    M8(7,5) = -M8(2,4)
    M8(8,5) = -M8(1,4)

    M8(1,6) = 0.d0
    M8(2,6) = -M8(7,3)
    M8(3,6) = -M8(6,3)
    M8(4,6) = -M8(5,3)
    M8(5,6) = -M8(4,3)
    M8(6,6) = -M8(3,3)
    M8(7,6) = -M8(2,3)
    M8(8,6) = -M8(1,3)

    M8(1,7) = 0.d0
    M8(2,7) = 0.d0
    M8(3,7) = -M8(6,2)
    M8(4,7) = -M8(5,2)
    M8(5,7) = -M8(4,2)
    M8(6,7) = -M8(3,2)
    M8(7,7) = -M8(2,2)
    M8(8,7) = -M8(1,2)

    M8(1,8) = 0.d0
    M8(2,8) = 0.d0
    M8(3,8) = 0.d0
    M8(4,8) = -M8(5,1)
    M8(5,8) = -M8(4,1)
    M8(6,8) = -M8(3,1)
    M8(7,8) = -M8(2,1)
    M8(8,8) = -M8(1,1)

    ! 6th-order
    M6(1,1) = -11.d0/180.d0 + M6_36
    M6(2,1) = 7.d0/60.d0 - 3.d0 * M6_36
    M6(3,1) = -1.d0/15.d0 + 3.d0 * M6_36
    M6(4,1) = -M6_36
    M6(5,1) = 0.d0
    M6(6,1) = 0.d0
    
    M6(1,2) = 1.d0/9.d0 - 2.d0 * M6_36
    M6(2,2) = -1.d0/120.d0 + 5.d0 * M6_36
    M6(3,2) = -11.d0/60.d0 - 3.d0 * M6_36
    M6(4,2) = 83.d0/360.d0 - M6_36
    M6(5,2) = -1.d0/90.d0 + M6_36
    M6(6,2) = 0.d0
    
    M6(1,3) = -1.d0/18.d0 + M6_36
    M6(2,3) = -17.d0/90.d0 - 2.d0 * M6_36
    M6(3,3) = -101.d0/360.d0 
    M6(4,3) = -137.d0/180.d0 + 2.d0 * M6_36
    M6(5,3) = -5.d0/72.d0 - M6_36
    M6(6,3) = -1.d0/180.d0
    
    M6(1,4) = -M6(6,3)
    M6(2,4) = -M6(5,3)
    M6(3,4) = -M6(4,3)
    M6(4,4) = -M6(3,3)
    M6(5,4) = -M6(2,3)
    M6(6,4) = -M6(1,3)
    
    M6(1,5) = 0.d0
    M6(2,5) = -M6(5,2)
    M6(3,5) = -M6(4,2)
    M6(4,5) = -M6(3,2)
    M6(5,5) = -M6(2,2)
    M6(6,5) = -M6(1,2)
    
    M6(1,6) = 0.d0
    M6(2,6) = 0.d0
    M6(3,6) = M6_36
    M6(4,6) = -M6(3,1)
    M6(5,6) = -M6(2,1)
    M6(6,6) = -M6(1,1)

    ! 4th-order
    M4(1,1) = 1.d0/8.d0
    M4(2,1) = -1.d0/6.d0
    M4(3,1) = 1.d0/8.d0
    M4(4,1) = 0.d0
    
    M4(1,2) = -1.d0/6.d0
    M4(2,2) = -3.d0/8.d0
    M4(3,2) = -2.d0/3.d0
    M4(4,2) = -1.d0/24.d0
    
    M4(1,3) = -M4(4,2)
    M4(2,3) = -M4(3,2)
    M4(3,3) = -M4(2,2)
    M4(4,3) = -M4(1,2)
    
    M4(1,4) = -M4(4,1)
    M4(2,4) = -M4(3,1)
    M4(3,4) = -M4(2,1)
    M4(4,4) = -M4(1,1)

    ! 2nd-order stencil
    M2(1,1) = -0.5d0
    M2(2,1) = -0.5d0
    
    M2(1,2) = 0.5d0
    M2(2,2) = 0.5d0

    ! right biased 2nd-order
    BRB(1,1) = BB2_33
    BRB(2,1) = 3.5d0 - 2.*BB2_33
    BRB(3,1) = -1.5d0 + BB2_33
    BRB(4,1) = 0.d0

    BRB(1,2) = 1.5d0 - 2.d0*BB2_33
    BRB(2,2) = -9.d0 + 4.d0*BB2_33
    BRB(3,2) = 2.5d0 - 2.d0*BB2_33
    BRB(4,2) = 0.d0

    BRB(1,3) = -1.5d0 + BB2_33
    BRB(2,3) = 5.5d0 - 2.d0*BB2_33
    BRB(3,3) = BB2_33
    BRB(4,3) = 0.d0

    BRB(1,4) = 0.d0
    BRB(2,4) = 0.d0
    BRB(3,4) = -1.d0
    BRB(4,4) = 0.d0

    do j=1,4
       do i=1,4
          BLB(i,j) = BRB(5-i,5-j)
       end do
    end do

  end subroutine stencil_init

  function first_deriv_8(u) result(du)
    double precision :: du
    double precision, intent(in) :: u(-4:4)
    du =   D8(1)*(u(1)-u(-1)) &
         + D8(2)*(u(2)-u(-2)) &
         + D8(3)*(u(3)-u(-3)) &
         + D8(4)*(u(4)-u(-4))
  end function first_deriv_8

  function first_deriv_6(u) result(du)
    double precision :: du
    double precision, intent(in) :: u(-3:3)
    du =   D6(1)*(u(1)-u(-1)) &
         + D6(2)*(u(2)-u(-2)) &
         + D6(3)*(u(3)-u(-3))
  end function first_deriv_6

  function first_deriv_4(u) result(du)
    double precision :: du
    double precision, intent(in) :: u(-2:2)
    du =   D4(1)*(u(1)-u(-1)) &
         + D4(2)*(u(2)-u(-2)) 
  end function first_deriv_4

  function first_deriv_r3(u) result(du)
    double precision :: du
    double precision, intent(in) :: u(-1:2)
    du = dot_product(u,DR3)
  end function first_deriv_r3

  function first_deriv_l3(u) result(du)
    double precision :: du
    double precision, intent(in) :: u(-2:1)
    du = dot_product(u,DL3)
  end function first_deriv_l3

  function first_deriv_rb(u) result(du)
    double precision :: du
    double precision, intent(in) :: u(0:3)
    du = dot_product(u,DRB)
  end function first_deriv_rb

  function first_deriv_lb(u) result(du)
    double precision :: du
    double precision, intent(in) :: u(-3:0)
    du = dot_product(u,DLB)
  end function first_deriv_lb

end module derivative_stencil_module
