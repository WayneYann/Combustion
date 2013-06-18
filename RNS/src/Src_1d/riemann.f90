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
    double precision, allocatable :: fl(:,:),fr(:,:),alpha_plus(:),alpha_mins(:),alpha_pm(:)

    allocate(fl(lo(1):hi(1)+1,NVAR))
    allocate(fr(lo(1):hi(1)+1,NVAR))
    allocate(alpha_plus(lo(1):hi(1)+1))
    allocate(alpha_mins(lo(1):hi(1)+1))
    allocate(alpha_pm  (lo(1):hi(1)+1))

    do i = lo(1), hi(1)+1
       alpha_plus(i) = 0.d0
       alpha_mins(i) = 0.d0
    end do

    call compute_flux_and_alpha(lo, hi, UL, fl, alpha_plus, alpha_mins)
    call compute_flux_and_alpha(lo, hi, UR, fr, alpha_plus, alpha_mins)

    do i = lo(1), hi(1)+1
       alpha_pm(i) = alpha_plus(i) * alpha_mins(i) / (alpha_plus(i) + alpha_mins(i))
       alpha_plus(i) = alpha_plus(i) / (alpha_plus(i) + alpha_mins(i))
       alpha_mins(i) = 1.d0 - alpha_plus(i)
    end do

    do n=1,NVAR
       if (n.eq.UTEMP) then
          do i = lo(1), hi(1)+1
             fx(i,n) = 0.d0
          end do
       else
          do i = lo(1), hi(1)+1
             fx(i,n) = alpha_plus(i) * fl(i,n) + alpha_mins(i) * fr(i,n) &
                  - alpha_pm(i) * (UR(i,n) - UL(i,n))
          end do
       end if
    end do

    deallocate(fl,fr,alpha_plus,alpha_mins,alpha_pm)
  end subroutine riemann


  subroutine compute_flux_and_alpha(lo, hi, U, F, ap, am)
    use eos_module, only : eos_given_ReY
    integer, intent(in) :: lo(1), hi(1)
    double precision, intent(in   ) ::  U(lo(1):hi(1)+1,NVAR)
    double precision, intent(out  ) ::  F(lo(1):hi(1)+1,NVAR)
    double precision, intent(inout) :: ap(lo(1):hi(1)+1)
    double precision, intent(inout) :: am(lo(1):hi(1)+1)

    integer :: i, n
    double precision :: rho, mx, rhoE
    double precision :: v, rhoInv, p, c, gamc, T, e, ek, H, Y(NSPEC)

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
          Y = U(i,UFS:UFS+NSPEC-1)*rhoInv
       end if

       call eos_given_ReY(gamc,p,c,T,rho,e,Y)

       ap(i) = max(ap(i), c+v)
       am(i) = max(am(i), c-v)

       F(i,URHO ) = mx
       F(i,UMX  ) = mx*v + p 
       F(i,UEDEN) = (rhoE + p) * v
       F(i,UTEMP) = 0.d0

       if (NSPEC > 0) then
          do n=1,NSPEC
             F(i,UFS+n-1) = U(i,UFS+n-1)*v
          end do
       end if
    end do
  end subroutine compute_flux_and_alpha

end module riemann_module
