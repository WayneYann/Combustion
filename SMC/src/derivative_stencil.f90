module derivative_stencil_module

  implicit none

  public

  integer, parameter :: compact=1, wide=2, S3D=2
  integer, save :: stencil, stencil_ng

  ! for 8th-order first derivatives
  double precision, parameter :: ALP =  0.8d0
  double precision, parameter :: BET = -0.2d0
  double precision, parameter :: GAM =  4.d0/105.d0
  double precision, parameter :: DEL = -1.d0/280.d0

  ! coefficients for 8th-order stencil of second derivatives
  double precision, parameter :: m47 = 683.d0/10080.d0, m48 = -1.d0/224.d0
  double precision, parameter :: m11 = 5.d0/336.d0 + m48, &
       &                         m12 = -83.d0/3600.d0 - m47/5.d0 - 14.d0*m48/5.d0, &
       &                         m13 = 299.d0/50400.d0 + 2.d0*m47/5.d0 + 13.d0*m48/5.d0, &
       &                         m14 = 17.d0/12600.d0 - m47/5.d0 - 4.d0*m48/5.d0, &
       &                         m15 = 1.d0/1120.d0, &
       &                         m21 = -11.d0/560.d0 - 2.d0*m48, &
       &                         m22 = -31.d0/360.d0 + m47 + 3.d0*m48, &
       &                         m23 = 41.d0/200.d0 - 9.d0*m47/5.d0 + 4.d0*m48/5.d0, &
       &                         m24 = -5927.d0/50400.d0 + 4.d0*m47/5.d0 - 9.d0*m48/5.d0, &
       &                         m25 = 17.d0/600.d0 - m47/5.d0 - 4.d0*m48/5.d0, &
       &                         m26 = -503.d0/50400.d0 + m47/5.d0 + 4.d0*m48/5.d0, &
       &                         m31 = -1.d0/280.d0, &
       &                         m32 = 1097.d0/5040.d0 - 2.d0*m47 + 6.d0*m48, &
       &                         m33 = -1349.d0/10080.d0 + 3.d0*m47 - 12.d0*m48, &
       &                         m34 = -887.d0/5040.d0 - m47 + 6.d0*m48, &
       &                         m35 = 3613.d0/50400.d0 + 4.d0*m47/5.d0 - 9.d0*m48/5.d0, &
       &                         m36 = 467.d0/25200.d0 - 3.d0*m47/5.d0 + 18.d0*m48/5.d0, &
       &                         m37 = 139.d0/25200.d0 - m47/5.d0 - 9.d0*m48/5.d0, &
       &                         m41 = 17.d0/1680.d0 + 2.d0*m48, &
       &                         m42 = -319.d0/2520.d0 + 2.d0*m47 - 8.d0*m48, &
       &                         m43 = -919.d0/5040.d0 - 2.d0*m47 + 6.d0*m48, &
       &                         m44 = -445.d0/2016.d0, &
       &                         m45 = 583.d0/720.d0 - m47 + 6.d0*m48, &
       &                         m46 = -65.d0/224.d0 - 7.d0*m48

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
