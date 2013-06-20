module weno_module

  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC, NCHARV, CFS

  implicit none

  private

  public reconstruct

contains

  ! L and R in UL and UR are relative to face
  subroutine reconstruct(lo, hi, U, Ulo, Uhi, UL, UR)
    use eos_module, only : eos_given_ReY
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in ) ::  U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: UL( lo(1): hi(1)+1,NVAR)
    double precision, intent(out) :: UR( lo(1): hi(1)+1,NVAR)

    integer,parameter :: dm = 1

    integer :: i, ii, ivar, m, n
    double precision :: egv(NCHARV,NCHARV)
    double precision :: gt, b, d(NSPEC)
    double precision :: rho, v, rhoInv, p, c, T, dpdr(NSPEC), dpde, e, ek, H, Y(NSPEC)
    double precision :: charv(-2:2,NCHARV), vp(NCHARV), vm(NCHARV) ! characteristic variables
    do n=1,NVAR
       do i=lo(1), hi(1)+1
          UL(i,n) = 0.d0
          UR(i,n) = 0.d0
       end do
    end do

    do i = lo(1)-1, hi(1)+1

       rho = 0.d0
       do n=1,nspec
          Y(n) = U(i,UFS+n-1)
          rho = rho + Y(n)
       end do

       rhoInv = 1.d0/rho

       do n=1,nspec
          Y(n) = Y(n) * rhoInv
       end do

       v  = U(i,UMX)*rhoInv
       ek = 0.5d0*v*v
       e = U(i,UEDEN)*rhoInV - ek
       T = U(i,UTEMP)

       call eos_given_ReY(p,c,T,dpdr,dpde,rho,e,Y)

       H = e + p*rhoInv + ek

       gt = dpde*rhoInv
       b = gt/(c*c)
       do n=1,nspec
          d(n) = b*(ek - e + dpdr(n)/gt)
       end do

       ! assemble left vectors
       egv(1,1) = -0.5d0*(1.d0/c + b*v)
       egv(2,1) = 0.5d0*b
       do n=1,nspec
          egv(CFS+n-1,1) = 0.5d0*(v/c + d(n))
       end do

       egv(1,2) = 0.5d0*(1.d0/c - b*v)
       egv(2,2) = 0.5d0*b
       do n=1,nspec
          egv(CFS+n-1,2) = 0.5d0*(-v/c + d(n))
       end do

       do m=1,nspec
          egv(1,CFS+m-1) = Y(m)*b*v
          egv(2,CFS+m-1) = -Y(m)*b
          do n=1,nspec
             egv(CFS+n-1,CFS+m-1) = -Y(m)*d(n)
          end do
          egv(CFS+m-1,CFS+m-1) = egv(CFS+m-1,CFS+m-1) + 1.d0
       end do

       ! convert conserved variables to characteristic variables
       do n=1,NCHARV
          do ii=-2,2
             charv(ii,n) = egv(1,n)*U(i+ii,UMX) + egv(2,n)*U(i+ii,UEDEN)
             do m=1,nspec
                charv(ii,n) = charv(ii,n) + egv(CFS+m-1,n)*U(i+ii,UFS+m-1)
             end do
          end do
       end do

       do ivar=1,NCHARV
          call weno5(charv(:,ivar), vp(ivar), vm(ivar))
       end do

       ! assemble right vectors
       egv(1,1) = v - c
       egv(2,1) = H - v*c
       do n=1,nspec
          egv(CFS+n-1,1) = Y(n)
       end do

       egv(1,2) = v + c
       egv(2,2) = H + v*c
       do n=1,nspec
          egv(CFS+n-1,2) = Y(n)
       end do

       do n=1,nspec
          egv(1,CFS+n-1) = v
          egv(2,CFS+n-1) = e + ek - dpdr(n)/gt
          do m=1,nspec
             egv(CFS+m-1,CFS+n-1) = 0.d0
          end do
          egv(CFS+n-1,CFS+n-1) = 1.d0
       end do

       if (i .ne. hi(1)+1) then
          do n=1,NCHARV
             UL(i+1,UMX  ) = UL(i+1,UMX  ) + vp(n)*egv(1,n)
             UL(i+1,UEDEN) = UL(i+1,UEDEN) + vp(n)*egv(2,n)
             do m=1,nspec
                UL(i+1,UFS+m-1) = UL(i+1,UFS+m-1) + vp(n)*egv(CFS+m-1,n)
             end do
          end do

          do m=1,nspec
             UL(i+1,URHO) = UL(i+1,URHO) + UL(i+1,UFS+m-1)
          end do

          UL(i+1,UTEMP) = U(i,UTEMP)
       end if

       if (i .ne. lo(1)-1) then
          do n=1,NCHARV
             UR(i,UMX  ) = UR(i,UMX  ) + vm(n)*egv(1,n)
             UR(i,UEDEN) = UR(i,UEDEN) + vm(n)*egv(2,n)
             do m=1,nspec
                UR(i,UFS+m-1) = UR(i,UFS+m-1) + vm(n)*egv(CFS+m-1,n)
             end do
          end do

          do m=1,nspec
             UR(i,URHO) = UR(i,URHO) + UR(i,UFS+m-1)
          end do

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
