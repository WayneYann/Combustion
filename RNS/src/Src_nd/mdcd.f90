module mdcd_module

  implicit none

  integer         , save :: wenop = 1
  double precision, save :: eps = 1.d-6

  double precision, parameter, private :: dsp = 0.0463783d0
  double precision, parameter, private :: dss = 0.01d0

  double precision, parameter, private :: b1=13.d0/12.d0, oneSixth=1.d0/6.d0

  double precision, dimension(0:3), parameter, private :: mdcd_Ck = &
       (/ 1.5d0*dsp+1.5d0*dss, 0.5d0-1.5d0*dsp+4.5d0*dss, 0.5d0-1.5d0*dsp-4.5d0*dss, 1.5d0*dsp-1.5d0*dss /)

  double precision, dimension(0:3), parameter, private :: weno5_Ck = &
       (/ 1.d0, 6.d0, 3.d0, 0.d0 /)

  double precision, dimension(0:3), parameter, private :: Ck = mdcd_Ck

  private

  public :: init_mdcd, mdcd

contains

  subroutine init_mdcd(weno_p_in,weno_eps_in)
    integer, intent(in) :: weno_p_in
    double precision, intent(in) :: weno_eps_in
    wenop = weno_p_in
    eps = weno_eps_in
  end subroutine init_mdcd


  subroutine mdcd(v, vl, vr)
    double precision, intent(in) :: v(-2:3)
    double precision, intent(out) :: vl, vr

    double precision :: vk(0:3), wk(0:3), IS(0:3)

    IS(0) = eps + b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    IS(1) = eps + b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    IS(2) = eps + b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
    IS(3) = eps + b1*(v( 1)-2.d0*v( 2)+v(3))**2 + 0.25d0*(5.d0*v(1)-8.d0*v(2)+3.d0*v(3))**2
    IS(3) = maxval(IS)

    wk = Ck / IS**wenop

    vk(0) = ( 2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0))*oneSixth
    vk(1) = (     -v(-1) + 5.d0*v( 0) +  2.d0*v(1))*oneSixth
    vk(2) = ( 2.d0*v( 0) + 5.d0*v( 1) -       v(2))*oneSixth
    vk(3) = (11.d0*v( 1) - 7.d0*v( 2) +  2.d0*v(3))*oneSixth

    vl = sum(wk*vk)/sum(wk)

    !-----------------------------------------------------------

    IS(0) = eps + b1*(v(3)-2.d0*v( 2)+v( 1))**2 + 0.25d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
    IS(1) = eps + b1*(v(2)-2.d0*v( 1)+v( 0))**2 + 0.25d0*(v(2)-v(0))**2
    IS(2) = eps + b1*(v(1)-2.d0*v( 0)+v(-1))**2 + 0.25d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2
    IS(3) = eps + b1*(v(0)-2.d0*v(-1)+v(-2))**2 + 0.25d0*(5.d0*v(0)-8.d0*v(-1)+3.d0*v(-2))**2
    IS(3) = maxval(IS)

    wk = Ck / IS**wenop

    vk = vk(3:0:-1)

    vr = sum(wk*vk)/sum(wk)

  end subroutine mdcd

end module mdcd_module
