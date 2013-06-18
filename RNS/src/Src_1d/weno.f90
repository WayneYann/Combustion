module weno_module

  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC

  implicit none

  private

  public reconstruct

contains

  ! L and R in UL and UR are relative to face
  subroutine reconstruct(lo, hi, U, Ulo, Uhi, UL, UR)
    use eos_module, only : eos_given_RTY
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in ) ::  U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: UL( lo(1): hi(1)+1,NVAR)
    double precision, intent(out) :: UR( lo(1): hi(1)+1,NVAR)

    integer :: i, ii, ivar
    double precision :: r1(3), r2(3), r3(3)  ! right-eigenvector
    double precision :: l1(3), l2(3), l3(3)  !  left-eigenvector
    double precision :: b1, b2
    double precision :: v, rhoInv, p, c, gamc, e, ek, H, Y(NSPEC)
    double precision :: charv(-2:2,3), vp(3), vm(3)  ! characteristic variables

    do i = lo(1)-1, hi(1)+1

       rhoInv = 1.0d0/U(i,URHO)
       v      = U(i,UMX)*rhoInv
       ek     = 0.5d0*v*v

       if (NSPEC > 0) then
          Y = U(i,UFS:UFS+NSPEC-1)*rhoInv
       end if

       call eos_given_RTY(e,p,c,gamc,U(i,URHO),U(i,UTEMP),Y)

       H = e + p*rhoInv + ek

       r1(1) = 1.0d0
       r1(2) = v - c
       r1(3) = H - v*c
       
       r2(1) = 1.0d0
       r2(2) = v + c
       r2(3) = H + v*c

       r3(1) = 1.0d0
       r3(2) = v
       r3(3) = ek

       b1 = p * rhoInv / (e*c*c)
       b2 = ek * b1

       l1(1) =  0.5d0*(b2 + v/c)
       l1(2) = -0.5d0*(b1*v + 1.0d0/c)
       l1(3) =  0.5d0*b1

       l2(1) =  0.5d0*(b2 - v/c)
       l2(2) = -0.5d0*(b1*v - 1.0d0/c)
       l2(3) =  0.5d0*b1

       l3(1) = 1.d0-b2
       l3(2) = b1*v
       l3(3) = -b1

       do ii=-2, 2
          charv(ii,1) = l1(1)*U(i+ii,URHO) + l1(2)*U(i+ii,UMX) + l1(3)*U(i+ii,UEDEN) 
          charv(ii,2) = l2(1)*U(i+ii,URHO) + l2(2)*U(i+ii,UMX) + l2(3)*U(i+ii,UEDEN) 
          charv(ii,3) = l3(1)*U(i+ii,URHO) + l3(2)*U(i+ii,UMX) + l3(3)*U(i+ii,UEDEN) 
       end do

       do ivar=1,3
          call weno5(charv(:,ivar), vp(ivar), vm(ivar))
       end do

       if (i .ne. hi(1)+1) then
          UL(i+1,URHO ) = vp(1)*r1(1) + vp(2)*r2(1) + vp(3)*r3(1)
          UL(i+1,UMX  ) = vp(1)*r1(2) + vp(2)*r2(2) + vp(3)*r3(2)
          UL(i+1,UEDEN) = vp(1)*r1(3) + vp(2)*r2(3) + vp(3)*r3(3)
          UL(i+1,UTEMP) = U(i,UTEMP)
          if (NSPEC .gt. 0) then
             print *, 'reconstruct: nspec > 0'
             stop
          end if
       end if

       if (i .ne. lo(1)-1) then
          UR(i,URHO ) = vm(1)*r1(1) + vm(2)*r2(1) + vm(3)*r3(1)
          UR(i,UMX  ) = vm(1)*r1(2) + vm(2)*r2(2) + vm(3)*r3(2)
          UR(i,UEDEN) = vm(1)*r1(3) + vm(2)*r2(3) + vm(3)*r3(3)
          UR(i,UTEMP) = U(i,UTEMP)
       end if

    end do

  end subroutine reconstruct


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

end module weno_module
