module weno_module

  implicit none

  ! coefficients for converting cell averages to two Gauss point values
  double precision, dimension(-2:2), parameter :: cg1 = &
       (/ -(1.d0+70.d0*sqrt(3.d0))/4320.d0, (4.d0+500.d0*sqrt(3.d0))/4320.d0, &
       4314.d0/4320.d0, (4.d0-500.d0*sqrt(3.d0))/4320.d0, (-1.d0+70.d0*sqrt(3.d0))/4320.d0 /)
  double precision, dimension(-2:2), parameter :: cg2 = &
       (/ (-1.d0+70.d0*sqrt(3.d0))/4320.d0, (4.d0-500.d0*sqrt(3.d0))/4320.d0, &
       4314.d0/4320.d0, (4.d0+500.d0*sqrt(3.d0))/4320.d0, -(1.d0+70.d0*sqrt(3.d0))/4320.d0 /)

  ! v_{i-1/2} = sum(cc4*v(i-2:1)) + O(h^4), where v(i-2:1) are cell averages
  double precision, dimension(-2:1), parameter :: cc4 = &
       (/  -1.d0/12.d0,  7.d0/12.d0,  7.d0/12.d0,  -1.d0/12.d0  /)

  private

  public :: weno5, cellavg2gausspt_1d

contains

  subroutine weno5(v, vp, vm)
    double precision, intent(in)  :: v(-2:2)
    double precision, intent(out) :: vp, vm ! v_{i+1/2} & v_{i-1/2}

    double precision, parameter :: epsw = 1.d-6, b1=13.d0/12.d0, oneSixth=1.d0/6.d0
    double precision :: vpr(-2:0), vmr(-2:0), beta(-2:0), alpha(-2:0), alpha1

    vpr(-2) = 2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
    vpr(-1) =     -v(-1) + 5.d0*v( 0) +  2.d0*v(1)
    vpr( 0) = 2.d0*v( 0) + 5.d0*v( 1) -       v(2)

    vmr(-2) =      -v(-2) + 5.d0*v(-1) + 2.d0*v(0)
    vmr(-1) =  2.d0*v(-1) + 5.d0*v(0 ) -      v(1) 
    vmr( 0) = 11.d0*v( 0) - 7.d0*v(1 ) + 2.d0*v(2)

    beta(-2) = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta(-1) = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta( 0) = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2

    beta(-2) = 1.d0/(epsw+beta(-2))**2
    beta(-1) = 1.d0/(epsw+beta(-1))**2
    beta( 0) = 1.d0/(epsw+beta( 0))**2

    alpha(-2) =      beta(-2)
    alpha(-1) = 6.d0*beta(-1)
    alpha( 0) = 3.d0*beta( 0)
    alpha1 = 1.d0/(alpha(-2) + alpha(-1) + alpha(0))

    vp = oneSixth*alpha1*(alpha(-2)*vpr(-2) + alpha(-1)*vpr(-1) + alpha(0)*vpr(0))

    alpha(-2) = 3.d0*beta(-2)
    alpha(-1) = 6.d0*beta(-1)
    alpha( 0) =      beta( 0)
    alpha1 = 1.d0/(alpha(-2) + alpha(-1) + alpha(0))

    vm = oneSixth*alpha1*(alpha(-2)*vmr(-2) + alpha(-1)*vmr(-1) + alpha(0)*vmr(0))

    return
  end subroutine weno5

  subroutine cellavg2gausspt_1d(u, ulo, uhi, u1, u2, lo, hi)
    integer, intent(in) :: ulo, uhi, lo, hi
    double precision, intent(in) :: u(ulo:uhi)
    double precision, intent(out) :: u1(lo:hi), u2(lo:hi)

    integer :: i

    do i=lo,hi
       u1(i) = cg1(-2)*u(i-2) + cg1(-1)*u(i-1) + cg1(0)*u(i) + cg1(1)*u(i+1) + cg1(2)*u(i+2)
       u2(i) = cg2(-2)*u(i-2) + cg2(-1)*u(i-1) + cg2(0)*u(i) + cg2(1)*u(i+1) + cg2(2)*u(i+2)
    end do
  end subroutine cellavg2gausspt_1d

end module weno_module
