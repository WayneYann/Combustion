module derivative_stencil_module

  implicit none

  public

  integer, parameter :: compact=1, wide=2, S3D=2
  integer, save :: stencil, stencil_ng

  ! for 8th-order first derivatives
  double precision,dimension(4),parameter :: D8 = (/ 0.8d0, -0.2d0, 4.d0/105.d0, -1.d0/280.d0 /)

  ! 6th-order first derivative
  double precision,dimension(3),parameter :: D6 = (/ 0.75d0, -0.15d0, 1.d0/60.d0 /)

  ! 4th-order first derivative
  double precision,dimension(2),parameter :: D4 = (/ 2.d0/3.d0, -1.d0/12.d0 /)

  ! 3rd-order first derivative
  ! The stencil is slighly right-based with one cell on the left and two on the right
  double precision,dimension(-1:2),parameter :: D3 = (/ -1.d0/3.d0, -0.5d0, 1.d0, -1.d0/6.d0 /)
  
  ! 3rd-order first derivative
  ! The stencil is right-based with no cell on the left and three on the right
  double precision,dimension(0:3),parameter :: DB = (/ -11.d0/6.d0, 3.d0, -1.5d0, 1.d0/3.d0 /)

  ! coefficients for 8th-order stencil of second derivatives
  double precision, private, parameter :: m47 = 683.d0/10080.d0, m48 = -1.d0/224.d0
  double precision, parameter, dimension(-4:-1,-4:3) :: DD8 = reshape( source = (/ &
       5.d0/336.d0 + m48, &
       -11.d0/560.d0 - 2.d0*m48, &
       -1.d0/280.d0, &
       17.d0/1680.d0 + 2.d0*m48, &
       !
       -83.d0/3600.d0 - m47/5.d0 - 14.d0*m48/5.d0, &
       -31.d0/360.d0 + m47 + 3.d0*m48, &
       1097.d0/5040.d0 - 2.d0*m47 + 6.d0*m48, &
       -319.d0/2520.d0 + 2.d0*m47 - 8.d0*m48, &
       !
       299.d0/50400.d0 + 2.d0*m47/5.d0 + 13.d0*m48/5.d0, &
       41.d0/200.d0 - 9.d0*m47/5.d0 + 4.d0*m48/5.d0, &
       -1349.d0/10080.d0 + 3.d0*m47 - 12.d0*m48, &
       -919.d0/5040.d0 - 2.d0*m47 + 6.d0*m48, &
       !
       17.d0/12600.d0 - m47/5.d0 - 4.d0*m48/5.d0, &
       -5927.d0/50400.d0 + 4.d0*m47/5.d0 - 9.d0*m48/5.d0, &
       -887.d0/5040.d0 - m47 + 6.d0*m48, &
       -445.d0/2016.d0, &
       !
       1.d0/1120.d0, &
       17.d0/600.d0 - m47/5.d0 - 4.d0*m48/5.d0, &
       3613.d0/50400.d0 + 4.d0*m47/5.d0 - 9.d0*m48/5.d0, &
       583.d0/720.d0 - m47 + 6.d0*m48, &
       !
       0.d0, &
       -503.d0/50400.d0 + m47/5.d0 + 4.d0*m48/5.d0, &
       467.d0/25200.d0 - 3.d0*m47/5.d0 + 18.d0*m48/5.d0, &
       -65.d0/224.d0 - 7.d0*m48, &
       !
       0.d0, &
       0.d0, &
       139.d0/25200.d0 - m47/5.d0 - 9.d0*m48/5.d0, &
       m47, &
       !
       0.d0, &
       0.d0, &
       0.d0, &
       m48   &
       /), shape = (/ 4,8 /) )

contains
  
  subroutine stencil_init

    use bl_error_module
    use probin_module, only : stencil_type

    if (trim(stencil_type) == "compact") then
       stencil = compact
    else if (trim(stencil_type) == "S3D" .or. trim(stencil_type) == "wide") then
       stencil = wide
    else
       call bl_error("unknow stencil_type")
    end if

    stencil_ng = 4

  end subroutine stencil_init

end module derivative_stencil_module
