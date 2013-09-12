module riemann_module

  use meth_params_module, only : ndim, NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC, &
       riemann_solver, HLL_solver, JBB_solver, HLLC_solver

  implicit none

  private

  public riemann

contains

  subroutine riemann(lo, hi, UL, UR, Ulo, Uhi, flx, flo, fhi, dir)
    integer, intent(in) :: lo, hi, Ulo, Uhi, flo, fhi
    integer, intent(in), optional :: dir
    double precision, intent(in ) ::  UL(Ulo:Uhi,NVAR)
    double precision, intent(in ) ::  UR(Ulo:Uhi,NVAR)
    double precision              :: flx(flo:fhi,NVAR)

    select case (riemann_solver)
    case (HLL_solver)
       call riemann_HLL(lo, hi, UL(lo:hi+1,:), UR(lo:hi+1,:), flx(lo:hi+1,:), dir)
    case (JBB_solver)
       call riemann_JBB(lo, hi, UL(lo:hi+1,:), UR(lo:hi+1,:), flx(lo:hi+1,:), dir)
    case (HLLC_solver)
       call riemann_HLLC(lo, hi, UL(lo:hi+1,:), UR(lo:hi+1,:), flx(lo:hi+1,:), dir)
    case default
       print *, 'unknown riemann solver'
       stop
    end select

  end subroutine riemann

  subroutine riemann_HLL(lo, hi, UL, UR, flx, dir)
    integer, intent(in) :: lo, hi
    integer, intent(in), optional :: dir
    double precision, intent(in ) ::  UL(lo:hi+1,NVAR)
    double precision, intent(in ) ::  UR(lo:hi+1,NVAR)
    double precision, intent(out) :: flx(lo:hi+1,NVAR)

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
             flx(i,n) = 0.d0
          end do
       else
          do i = lo, hi+1
             flx(i,n) = alpha_plus(i) * fl(i,n) + alpha_mins(i) * fr(i,n) &
                  - alpha_pm(i) * (UR(i,n) - UL(i,n))
          end do
       end if
    end do

    deallocate(fl,fr,alpha_plus,alpha_mins,alpha_pm)
  end subroutine riemann_HLL


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
    double precision :: rhoInv, p, c, gamc, dpdr(NSPEC), dpde, T, e, ek, Y(NSPEC)

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

       call eos_given_ReY(p,c,gamc,T,dpdr,dpde,rho,e,Y)

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


  subroutine riemann_JBB(lo, hi, UL, UR, flx, dir)
!    use prob_params_module, only : physbc_lo,Symmetry
    use eos_module, only : smalld, smallp, eos_given_ReY
    integer, intent(in) :: lo, hi
    integer, intent(in), optional :: dir
    double precision, intent(in ) ::  UL(lo:hi+1,NVAR)
    double precision, intent(in ) ::  UR(lo:hi+1,NVAR)
    double precision, intent(out) :: flx(lo:hi+1,NVAR)

    integer :: i, n, idir, ivel(3), idim
    double precision :: vflag(3), dpdr(NSPEC), dpde
    double precision :: rgdnv,regdnv, pgdnv, vgdnv(3), ekgdnv
    double precision :: rl, retl, Tl, Yl(NSPEC), vl(3), el, pl, rel, rinvl
    double precision :: rr, retr, Tr, Yr(NSPEC), vr(3), er, pr, rer, rinvr
    double precision :: wl, cl, gamcl
    double precision :: wr, cr, gamcr
    double precision :: csmall, wsmall 
    double precision :: rstar, cstar, estar, pstar, ustar
    double precision :: ro, uo, po, reo, gamco, co, entho
    double precision :: sgnm, spout, spin, ushock, scr, frac
    
    double precision, parameter :: small  = 1.d-8

    if (present(dir)) then
       idir = dir
    else
       idir = 1
    end if

    call set_vel(idir, ivel, vflag)

    do i=lo,hi+1
       
       rl    = max(UL(i,URHO), smalld)
       rinvl = 1.d0/rl
       vl(1) = UL(i,ivel(1)) * vflag(1) * rinvl
       vl(2) = UL(i,ivel(2)) * vflag(2) * rinvl
       vl(3) = UL(i,ivel(3)) * vflag(3) * rinvl
       Tl    = UL(i,UTEMP)
       retl  = UL(i,UEDEN)
       Yl    = UL(i,UFS:UFS+nspec-1) * rinvl
       rel   = retl - 0.5d0*rl*(vl(1)**2+vl(2)**2+vl(3)**2)
       el    = rel * rinvl       
       call eos_given_ReY(pl, cl, gamcl, Tl, dpdr, dpde, rl, el, Yl)

       rr    = max(UR(i,URHO), smalld)
       rinvr = 1.d0/rr
       vr(1) = UR(i,ivel(1)) * vflag(1) * rinvr
       vr(2) = UR(i,ivel(2)) * vflag(2) * rinvr
       vr(3) = UR(i,ivel(3)) * vflag(3) * rinvr
       Tr    = UR(i,UTEMP)
       retr  = UR(i,UEDEN)
       Yr    = UR(i,UFS:UFS+nspec-1) * rinvr
       rer   = retr - 0.5d0*rr*(vr(1)**2+vr(2)**2+vr(3)**2)
       er    = rer * rinvr       
       call eos_given_ReY(pr, cr, gamcr, Tr, dpdr, dpde, rr, er, Yr)

       csmall = max(small, small*cl, small*cr)
       wsmall = smalld*csmall
       wl = max(wsmall, cl*rl)
       wr = max(wsmall, cr*rr)

       pstar = ((wr*pl + wl*pr) + wl*wr*(vl(1) - vr(1)))/(wl + wr)
       ustar = ((wl*vl(1) + wr*vr(1)) + (pl - pr))/(wl + wr)
       pstar = max(pstar,smallp)

       if (ustar .gt. 0.d0) then
          ro = rl
          uo = vl(1)
          po = pl
          reo = rel
          gamco = gamcl
       else if (ustar .lt. 0.d0) then
          ro = rr
          uo = vr(1)
          po = pr
          reo = rer
          gamco = gamcr
       else
          ro = 0.5d0*(rl+rr)
          uo = 0.5d0*(vl(1)+vr(1))
          po = 0.5d0*(pl+pr)
          reo = 0.5d0*(rel+rer)
          gamco = 0.5d0*(gamcl+gamcr)
       endif
       ro = max(smalld,ro)
       
       co = sqrt(abs(gamco*po/ro))
       co = max(csmall,co)
       entho = (reo/ro + po/ro)/co**2
       rstar = ro + (pstar - po)/co**2
       rstar = max(smalld,rstar)
       estar = reo + (pstar - po)*entho
       cstar = sqrt(abs(gamco*pstar/rstar))
       cstar = max(cstar,csmall)

       sgnm = sign(1.d0,ustar)
       spout = co - sgnm*uo
       spin = cstar - sgnm*ustar
       ushock = 0.5d0*(spin + spout)
       if (pstar-po .ge. 0.d0) then
          spin = ushock
          spout = ushock
       endif
       if (spout-spin .eq. 0.d0) then
          scr = small*0.5d0*(cl+cr)
       else
          scr = spout-spin
       endif
       frac = (1.d0 + (spout + spin)/scr)*0.5d0
       frac = max(0.d0,min(1.d0,frac))

       if (ustar .gt. 0.d0) then
          vgdnv(2:3) = vl(2:3)
       else if (ustar .lt. 0.d0) then
          vgdnv(2:3) = vr(2:3)
       else
          vgdnv(2:3) = 0.5d0*(vl(2:3)+vr(2:3))
       endif
       rgdnv = frac*rstar + (1.d0 - frac)*ro

       vgdnv(1) = frac*ustar + (1.d0 - frac)*uo
       pgdnv    = frac*pstar + (1.d0 - frac)*po

       regdnv = frac*estar + (1.d0 - frac)*reo
       if (spout .lt. 0.d0) then
          rgdnv    = ro
          vgdnv(1) = uo
          pgdnv    = po
          regdnv   = reo
       endif
       if (spin .ge. 0.d0) then
          rgdnv    = rstar
          vgdnv(1) = ustar
          pgdnv    = pstar
          regdnv   = estar
       endif

       pgdnv = max(pgdnv,smallp)

       ! Enforce that fluxes through a symmetry plane are hard zero.
!       if (i .eq.0 .and. physbc_lo(idir) .eq. Symmetry) vgdnv(1) = 0.d0

       flx(i,URHO) = rgdnv*vgdnv(1)

       ekgdnv = 0.d0
       do idim=1,ndim
          flx(i,ivel(idim)) = flx(i,URHO)*vgdnv(idim)
          ekgdnv = ekgdnv + 0.5d0*vgdnv(idim)**2
       end do
       flx(i,ivel(1)) = flx(i,ivel(1)) + pgdnv

       flx(i,UEDEN) = vgdnv(1)*(regdnv + rgdnv*ekgdnv + pgdnv)
       flx(i,UTEMP) = 0.d0

       if (ustar .gt. 0.d0) then
          do n=1,NSPEC
             flx(i,UFS+n-1) = flx(i,URHO) * Yl(n)
          end do
       else if (ustar .lt. 0.d0) then
          do n=1,NSPEC
             flx(i,UFS+n-1) = flx(i,URHO) * Yr(n)
          end do
       else 
          do n=1,NSPEC
             flx(i,UFS+n-1) = flx(i,URHO) * 0.5d0*(Yl(n)+Yr(n))
          end do
       end if

    end do

  end subroutine riemann_JBB


  subroutine riemann_HLLC(lo, hi, UL, UR, flx, dir)
!    use prob_params_module, only : physbc_lo,Symmetry
    use eos_module, only : smalld, smallp, eos_given_ReY
    integer, intent(in) :: lo, hi
    integer, intent(in), optional :: dir
    double precision, intent(in ) ::  UL(lo:hi+1,NVAR)
    double precision, intent(in ) ::  UR(lo:hi+1,NVAR)
    double precision, intent(out) :: flx(lo:hi+1,NVAR)

    integer :: i, n, idir, ivel(3), idim
    double precision :: vflag(3), dpdr(NSPEC), dpde
    double precision :: rl, retl, Tl, Yl(NSPEC), vl(3), el, pl, rel, rinvl
    double precision :: rr, retr, Tr, Yr(NSPEC), vr(3), er, pr, rer, rinvr
    double precision :: cl, gamcl, Smul, Sl
    double precision :: cr, gamcr, Smur, Sr
    double precision :: Sstar, rstar, vstar(3), etstar

    if (present(dir)) then
       idir = dir
    else
       idir = 1
    end if

    call set_vel(idir, ivel, vflag)

    do i=lo,hi+1
       
       rl    = max(UL(i,URHO), smalld)
       rinvl = 1.d0/rl
       vl(1) = UL(i,ivel(1)) * vflag(1) * rinvl
       vl(2) = UL(i,ivel(2)) * vflag(2) * rinvl
       vl(3) = UL(i,ivel(3)) * vflag(3) * rinvl
       Tl    = UL(i,UTEMP)
       retl  = UL(i,UEDEN)
       Yl    = UL(i,UFS:UFS+nspec-1) * rinvl
       rel   = retl - 0.5d0*rl*(vl(1)**2+vl(2)**2+vl(3)**2)
       el    = rel * rinvl       
       call eos_given_ReY(pl, cl, gamcl, Tl, dpdr, dpde, rl, el, Yl)

       rr    = max(UR(i,URHO), smalld)
       rinvr = 1.d0/rr
       vr(1) = UR(i,ivel(1)) * vflag(1) * rinvr
       vr(2) = UR(i,ivel(2)) * vflag(2) * rinvr
       vr(3) = UR(i,ivel(3)) * vflag(3) * rinvr
       Tr    = UR(i,UTEMP)
       retr  = UR(i,UEDEN)
       Yr    = UR(i,UFS:UFS+nspec-1) * rinvr
       rer   = retr - 0.5d0*rr*(vr(1)**2+vr(2)**2+vr(3)**2)
       er    = rer * rinvr       
       call eos_given_ReY(pr, cr, gamcr, Tr, dpdr, dpde, rr, er, Yr)

       ! wave speed estimates
       Sl = min(vl(1)-cl, vr(1)-cr)
       Sr = max(vl(1)+cl, vr(1)+cr)
       Smul = Sl - vl(1)
       Smur = Sr - vr(1)
       Sstar = (pr-pl+rl*vl(1)*Smul-rr*vr(1)*Smur)/(rl*Smul-rr*Smur)

       ! HLLC flux
       if (Sl .ge. 0.d0) then

          flx(i,URHO) = rl*vl(1)

          do idim=1,ndim
             flx(i,ivel(idim)) = flx(i,URHO)*vl(idim)
          end do
          flx(i,ivel(1)) = flx(i,ivel(1)) + pl

          flx(i,UEDEN) = vl(1)*(retl + pl)
          flx(i,UTEMP) = 0.d0

          do n=1,NSPEC
             flx(i,UFS+n-1) = flx(i,URHO) * Yl(n)
          end do

       else if (Sr .le. 0.d0) then

          flx(i,URHO) = rr*vr(1)

          do idim=1,ndim
             flx(i,ivel(idim)) = flx(i,URHO)*vr(idim)
          end do
          flx(i,ivel(1)) = flx(i,ivel(1)) + pr

          flx(i,UEDEN) = vr(1)*(retr + pr)
          flx(i,UTEMP) = 0.d0

          do n=1,NSPEC
             flx(i,UFS+n-1) = flx(i,URHO) * Yr(n)
          end do

       else if (Sstar .ge. 0.d0) then

          rstar = rl*Smul/(Sl-Sstar)
          vstar(1) = Sstar
          vstar(2) = vl(2)
          vstar(3) = vl(3)
          etstar = retl*rinvl+(Sstar-vl(1))*(Sstar+pl/(rl*Smul))

          flx(i,URHO) = rl*vl(1) + Sl*(rstar-rl)

          do idim=1,ndim
             flx(i,ivel(idim)) = rl*vl(1)*vl(idim) + Sl*(rstar*vstar(idim)-rl*vl(idim))
          end do
          flx(i,ivel(1)) = flx(i,ivel(1)) + pl

          flx(i,UEDEN) = vl(1)*(retl + pl) + Sl*(rstar*etstar-retl)
          flx(i,UTEMP) = 0.d0

          do n=1,NSPEC
             flx(i,UFS+n-1) = flx(i,URHO) * Yl(n)
          end do

       else

          rstar = rr*Smur/(Sr-Sstar)
          vstar(1) = Sstar
          vstar(2) = vr(2)
          vstar(3) = vr(3)
          etstar = retr*rinvr+(Sstar-vr(1))*(Sstar+pr/(rr*Smur))

          flx(i,URHO) = rr*vr(1) + Sr*(rstar-rr)

          do idim=1,ndim
             flx(i,ivel(idim)) = rr*vr(1)*vr(idim) + Sr*(rstar*vstar(idim)-rr*vr(idim))
          end do
          flx(i,ivel(1)) = flx(i,ivel(1)) + pr

          flx(i,UEDEN) = vr(1)*(retr + pr) + Sr*(rstar*etstar-retr)
          flx(i,UTEMP) = 0.d0

          do n=1,NSPEC
             flx(i,UFS+n-1) = flx(i,URHO) * Yr(n)
          end do

       end if

    end do

  end subroutine riemann_HLLC


  subroutine set_vel(idir, ivel, vflag)
    integer, intent(in) :: idir
    integer, intent(out) :: ivel(3)
    double precision, intent(out) :: vflag(3)
    if (ndim .eq. 1) then
       vflag(1) = 1.d0
       vflag(2) = 0.d0
       vflag(3) = 0.d0
       ivel(1) = UMX
       ivel(2) = UMX
       ivel(3) = UMX
    else if (ndim .eq. 2) then
       vflag(1) = 1.d0
       vflag(2) = 1.d0
       vflag(3) = 0.d0   
       if (idir .eq. 1) then
          ivel(1) = UMX
          ivel(2) = UMY
          ivel(3) = UMX
       else
          ivel(1) = UMY
          ivel(2) = UMX
          ivel(3) = UMX
       end if
    else
       vflag(1) = 1.d0
       vflag(2) = 1.d0
       vflag(3) = 1.d0 
       if (idir .eq. 1) then
          ivel(1) = UMX
          ivel(2) = UMY
          ivel(3) = UMZ
       else if (idir .eq. 2) then
          ivel(1) = UMY
          ivel(2) = UMZ
          ivel(3) = UMX
       else
          ivel(1) = UMZ
          ivel(2) = UMX
          ivel(3) = UMY
       end if
    end if
  end subroutine set_vel


end module riemann_module
