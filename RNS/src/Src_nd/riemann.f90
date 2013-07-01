module riemann_module

  use meth_params_module, only : ndim, NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC

  implicit none

  private

  public riemann

contains

  subroutine riemann(lo, hi, UL, UR, fx, dir)
    integer, intent(in) :: lo, hi
    integer, intent(in), optional :: dir
    double precision, intent(in ) :: UL(lo:hi+1,NVAR)
    double precision, intent(in ) :: UR(lo:hi+1,NVAR)
    double precision, intent(out) :: fx(lo:hi+1,NVAR)

    integer :: i, n
    double precision, allocatable :: fl(:,:),fr(:,:),alpha_plus(:),alpha_mins(:),alpha_pm(:)

    allocate(fl(lo:hi+1,NVAR))
    allocate(fr(lo:hi+1,NVAR))
    allocate(alpha_plus(lo:hi+1))
    allocate(alpha_mins(lo:hi+1))
    allocate(alpha_pm  (lo:hi+1))

    do i = lo, hi+1
       alpha_plus(i) = 0.d0
       alpha_mins(i) = 0.d0
    end do

    call compute_flux_and_alpha(lo, hi, UL, fl, alpha_plus, alpha_mins, dir)
    call compute_flux_and_alpha(lo, hi, UR, fr, alpha_plus, alpha_mins, dir)

    do i = lo, hi+1
       alpha_pm(i) = alpha_plus(i) * alpha_mins(i) / (alpha_plus(i) + alpha_mins(i))
       alpha_plus(i) = alpha_plus(i) / (alpha_plus(i) + alpha_mins(i))
       alpha_mins(i) = 1.d0 - alpha_plus(i)
    end do

    do n=1,NVAR
       if (n.eq.UTEMP) then
          do i = lo, hi+1
             fx(i,n) = 0.d0
          end do
       else
          do i = lo, hi+1
             fx(i,n) = alpha_plus(i) * fl(i,n) + alpha_mins(i) * fr(i,n) &
                  - alpha_pm(i) * (UR(i,n) - UL(i,n))
          end do
       end if
    end do

    deallocate(fl,fr,alpha_plus,alpha_mins,alpha_pm)
  end subroutine riemann


  subroutine compute_flux_and_alpha(lo, hi, U, F, ap, am, dir)
    use eos_module, only : eos_given_ReY
    integer, intent(in) :: lo, hi
    integer, intent(in), optional :: dir
    double precision, intent(in   ) ::  U(lo:hi+1,NVAR)
    double precision, intent(out  ) ::  F(lo:hi+1,NVAR)
    double precision, intent(inout) :: ap(lo:hi+1)
    double precision, intent(inout) :: am(lo:hi+1)

    integer :: i, n, idir
    double precision :: rho, m(3), rhoE, v(3), vn
    double precision :: rhoInv, p, c, dpdr(NSPEC), dpde, T, e, ek, Y(NSPEC)

    if (present(dir)) then
       idir = dir
    else
       idir = 1
    end if

    do i=lo,hi+1

       rho  = U(i,URHO)
       m(1) = U(i,UMX)
       if (ndim .ge. 2) then
          m(2) = U(i,UMY)
       else
          m(2) = 0.d0
       end if
       if (ndim .eq. 3) then
          m(3) = U(i,UMZ)
       else
          m(3) = 0.d0
       end if
       rhoE = U(i,UEDEN)
       
       rhoInv = 1.0d0/rho
       v      = m*rhoInv
       T      = U(i,UTEMP)
       
       ek = 0.5d0*(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
       e  = rhoE*rhoInv - ek

       Y = U(i,UFS:UFS+NSPEC-1)*rhoInv

       call eos_given_ReY(p,c,T,dpdr,dpde,rho,e,Y)

       vn = v(idir)

       ap(i) = max(ap(i), c+vn)
       am(i) = max(am(i), c-vn)

       F(i,URHO ) = rho*vn
       F(i,UMX  ) = m(1)*vn
       if (ndim .ge. 2) then
          F(i,UMY  ) = m(2)*vn
       end if
       if (ndim .eq. 3) then
          F(i,UMZ  ) = m(3)*vn
       end if
       F(i,UMX+idir-1) = F(i,UMX+idir-1) + p
       F(i,UEDEN) = (rhoE + p) * vn
       F(i,UTEMP) = 0.d0
       do n=1,NSPEC
          F(i,UFS+n-1) = U(i,UFS+n-1)*vn
       end do
    end do
  end subroutine compute_flux_and_alpha

end module riemann_module
