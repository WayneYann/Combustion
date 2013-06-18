module riemann_module

  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC

  implicit none

  private

  public riemann

contains

  subroutine riemann(lo, hi, UL, UR, fx)
    integer, intent(in) :: lo(1), hi(1)
    double precision, intent(in ) :: UL(lo(1):hi(1)+1,NVAR)
    double precision, intent(in ) :: UR(lo(1):hi(1)+1,NVAR)
    double precision, intent(out) :: fx(lo(1):hi(1)+1,NVAR)

    integer :: i, n
    double precision, allocatable :: fl(:,:), fr(:,:), alpha_plus(:), alpha_mins(:)

    allocate(fl(lo(1):hi(1)+1,NVAR))
    allocate(fr(lo(1):hi(1)+1,NVAR))
    allocate(alpha_plus(lo(1):hi(1)+1))
    allocate(alpha_mins(lo(1):hi(1)+1))

    do i = lo(1), hi(1)+1
       alpha_plus(i) = 0.d0
       alpha_mins(i) = 0.d0
    end do
       
    call compute_flux_and_alpha(lo, hi, UL, fl, alpha_plus, alpha_mins)
    call compute_flux_and_alpha(lo, hi, UR, fr, alpha_plus, alpha_mins)

    do n=1,NVAR
       do i = lo(1), hi(1)+1
          fx(i,n) = (alpha_plus(i) * fl(i,n) + alpha_mins(i) * fr(i,n) &
               - alpha_plus(i) * alpha_mins(i) * (UR(i,n) - UL(i,n)))  &
               / (alpha_plus(i) + alpha_mins(i))
       end do
    end do

    deallocate(fl,fr,alpha_plus,alpha_mins)
  end subroutine riemann


  subroutine compute_flux_and_alpha(lo, hi, U, F, ap, am)
    use eos_module, only : eos_get_pcg
    integer, intent(in) :: lo(1), hi(1)
    double precision, intent(in   ) ::  U(lo(1):hi(1)+1,NVAR)
    double precision, intent(out  ) ::  F(lo(1):hi(1)+1,NVAR)
    double precision, intent(inout) :: ap(lo(1):hi(1)+1)
    double precision, intent(inout) :: am(lo(1):hi(1)+1)

    integer :: i, n
    double precision :: rho, mx, rhoE
    double precision :: v, rhoInv, p, c, gamc, T, e, ek, H, xn(NSPEC)

    do i=lo(1),hi(1)+1

       rho  = U(i,URHO)
       mx   = U(i,UMX)
       rhoE = U(i,UEDEN)
       
       rhoInv = 1.0d0/rho
       v      = mx*rhoInv
       T      = U(i,UTEMP)
       
       ek = 0.5d0*v*v
       e  = rhoE*rhoInv - ek

       if (NSPEC > 0) then
          xn = u(i,UFS:UFS+NSPEC-1)*rhoInv
       end if

       call eos_get_pcg(p,c,gamc,rho,e,T,xn)

       ap(i) = max(ap(i), v+c)
       am(i) = max(am(i), c-v)

       F(i,URHO ) = mx
       F(i,UMX  ) = mx*v + p 
       F(i,UEDEN) = (rhoE + p) * v
       F(i,UTEMP) = 0.d0

       if (NSPEC .gt. 0) then
          do n=1,NSPEC
             F(i,UFS+n-1) = U(i,UFS+n-1)*v
          end do
       end if
    end do
  end subroutine compute_flux_and_alpha

end module riemann_module
