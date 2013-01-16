module riemann_module

  implicit none

  private

  public riemannus

contains

      subroutine riemannus(ql,qr,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                           gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                           uflx,uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3, &
                           ugdnv,pgdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                           idir,ilo,ihi,jlo,jhi,kc,kflux)

      use chemistry_module, only : nspec=>nspecies
      use prob_params_module, only : physbc_lo,physbc_hi,Symmetry
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, &
                                     nadv, small_dens, small_pres

      implicit none

      integer,intent(in):: qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer,intent(in):: gd_l1,gd_l2,gd_h1,gd_h2
      integer,intent(in):: uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3
      integer,intent(in):: pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
      integer,intent(in):: idir,ilo,ihi,jlo,jhi
      integer,intent(in):: kc,kflux

      double precision,intent(in ):: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision,intent(in ):: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision,intent(in )::  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision,intent(in )::  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision,intent(in )::    cav(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision,intent(in ):: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision,intent(out):: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,uflx_l3:uflx_h3,NVAR)
      double precision,intent(out):: ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
      double precision,intent(out):: pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)

      double precision, parameter:: small = 1.d-8

      integer i, j, n, nq
      integer iadv, ispec

      double precision rgdnv,v1gdnv,v2gdnv,regdnv,ustar
      double precision rl, ul, v1l, v2l, pl, rel
      double precision rr, ur, v1r, v2r, pr, rer
      double precision wl, wr, rhoetot, scr
      double precision rstar, cstar, estar, pstar
      double precision ro, uo, po, reo, co, gamco, entho
      double precision sgnm, spin, spout, ushock, frac
      double precision wsmall, csmall,qavg

      !$OMP PARALLEL DO PRIVATE(i,j,rl,ul,v1l,v2l,pl,rel,rr,ur,v1r,v2r,pr,rer,csmall,wsmall,wl,wr,pstar,ustar,ro,uo) &
      !$OMP PRIVATE(po,reo,gamco,co,entho,rstar,estar,cstar,sgnm,spout,spin,ushock,scr,frac,v1gdnv,v2gdnv,rgdnv,regdnv) &
      !$OMP PRIVATE(rhoetot,iadv,n,nq,qavg,ispec)
      do j = jlo, jhi
         do i = ilo, ihi

            rl = ql(i,j,kc,QRHO)

            ! pick left velocities based on direction
            if(idir.eq.1) then
               ul  = ql(i,j,kc,QU)
               v1l = ql(i,j,kc,QV)
               v2l = ql(i,j,kc,QW)
            elseif(idir.eq.2) then
               ul  = ql(i,j,kc,QV)
               v1l = ql(i,j,kc,QU)
               v2l = ql(i,j,kc,QW)
            else
               ul  = ql(i,j,kc,QW)
               v1l = ql(i,j,kc,QU)
               v2l = ql(i,j,kc,QV)
            endif

            pl = ql(i,j,kc,QPRES)
            rel = ql(i,j,kc,QREINT)

            rr = qr(i,j,kc,QRHO)

            ! pick right velocities based on direction
            if(idir.eq.1) then
               ur  = qr(i,j,kc,QU)
               v1r = qr(i,j,kc,QV)
               v2r = qr(i,j,kc,QW)
            elseif(idir.eq.2) then
               ur  = qr(i,j,kc,QV)
               v1r = qr(i,j,kc,QU)
               v2r = qr(i,j,kc,QW)
            else
               ur  = qr(i,j,kc,QW)
               v1r = qr(i,j,kc,QU)
               v2r = qr(i,j,kc,QV)
            endif

            pr  = qr(i,j,kc,QPRES)
            rer = qr(i,j,kc,QREINT)

            csmall = smallc(i,j)
            wsmall = small_dens*csmall
            wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
            wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

            pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
            ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
            pstar = max(pstar,small_pres)

            if (ustar .gt. 0.d0) then
               ro = rl
               uo = ul
               po = pl
               reo = rel
               gamco = gamcl(i,j)
            else if (ustar .lt. 0.d0) then
               ro = rr
               uo = ur
               po = pr
               reo = rer
               gamco = gamcr(i,j)
            else
               ro = 0.5d0*(rl+rr)
               uo = 0.5d0*(ul+ur)
               po = 0.5d0*(pl+pr)
               reo = 0.5d0*(rel+rer)
               gamco = 0.5d0*(gamcl(i,j)+gamcr(i,j))
            endif
            ro = max(small_dens,ro)
         
            co = sqrt(abs(gamco*po/ro))
            co = max(csmall,co)
            entho = (reo/ro + po/ro)/co**2
            rstar = ro + (pstar - po)/co**2
            rstar = max(small_dens,rstar)
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
               scr = small*cav(i,j)
            else
               scr = spout-spin
            endif
            frac = (1.d0 + (spout + spin)/scr)*0.5d0
            frac = max(0.d0,min(1.d0,frac))

            if (ustar .gt. 0.d0) then
               v1gdnv = v1l
               v2gdnv = v2l
            else if (ustar .lt. 0.d0) then
               v1gdnv = v1r
               v2gdnv = v2r
            else
               v1gdnv = 0.5d0*(v1l+v1r)
               v2gdnv = 0.5d0*(v2l+v2r)
            endif
            rgdnv = frac*rstar + (1.d0 - frac)*ro

            ugdnv(i,j,kc) = frac*ustar + (1.d0 - frac)*uo
            pgdnv(i,j,kc) = frac*pstar + (1.d0 - frac)*po

            regdnv = frac*estar + (1.d0 - frac)*reo
            if (spout .lt. 0.d0) then
               rgdnv = ro
               ugdnv(i,j,kc) = uo
               pgdnv(i,j,kc) = po
               regdnv = reo
            endif
            if (spin .ge. 0.d0) then
               rgdnv = rstar
               ugdnv(i,j,kc) = ustar
               pgdnv(i,j,kc) = pstar
               regdnv = estar
            endif

            ! Enforce that fluxes through a symmetry plane are hard zero.
            if (i    .eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) &
                 ugdnv(i,j,kc) = 0.d0
            if (j    .eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) &
                 ugdnv(i,j,kc) = 0.d0
            if (kflux.eq.0 .and. physbc_lo(3) .eq. Symmetry .and. idir .eq. 3) &
                 ugdnv(i,j,kc) = 0.d0

            ! Compute fluxes, order as conserved state (not q)
            uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)

            if(idir.eq.1) then
               uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
               uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v1gdnv
               uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
            elseif(idir.eq.2) then
               uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
               uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
               uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
            else
               uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
               uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v2gdnv
               uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
            endif

            rhoetot = regdnv + 0.5d0*rgdnv*(ugdnv(i,j,kc)**2 + v1gdnv**2 + v2gdnv**2)

            uflx(i,j,kflux,UEDEN) = ugdnv(i,j,kc)*(rhoetot + pgdnv(i,j,kc))

            uflx(i,j,kflux,UEINT) = ugdnv(i,j,kc)*regdnv

            do iadv = 1, nadv
               n  = UFA + iadv - 1
               nq = QFA + iadv - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
               else
                  qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
               endif
            enddo

            do ispec = 1, nspec
               n  = UFS + ispec - 1
               nq = QFS + ispec - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
               else
                  qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
               endif
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine riemannus

end module riemann_module
