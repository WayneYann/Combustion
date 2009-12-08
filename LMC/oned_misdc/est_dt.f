      subroutine est_dt(nx,vel,scal,divu,dsdt,cfl,umax,dx,dt)
      implicit none
      integer nx
      real*8    vel(-1:nx)
      real*8   scal(-1:nx,*)
      real*8   divu(0 :nx-1)
      real*8   dsdt(0 :nx-1)
      real*8  cfl
      real*8  umax
      real*8  dx
      real*8  dt

      real*8  small
      real*8  dt_start
      real*8  dt_divu
      integer i

      small   = 1.0d-8

      umax = 0.d0
      do i = 0,nx-1
        umax = max(umax,abs(vel(i)))
      end do

      dt_start = 1.0D+20
      dt = dt_start

      if (umax .gt. small) dt = min(dt,dx/umax)
      if (dt .eq. dt_start) dt = dx

      print *,'DT = DX / V = ',dt,' USING V =',umax

      dt = dt*cfl

      print *,'DT * CFL ',dt

      call est_divu_dt(divu,dsdt,scal,dx,dt_divu)

      print *,'ESTDT: UMAX = ',umax

      print *,'DT = MIN(',dt,dt_divu,')'

      dt = min(dt,dt_divu)

      end

      subroutine est_divu_dt(divu,dsdt,scal,dx,dt)
      implicit none
      include 'spec.h'
      real*8   divu(0 :nx-1)
      real*8   dsdt(0 :nx-1)
      real*8   scal(-1:nx  ,*)
      real*8  dx
      real*8  dt

      real*8  rho
      real*8  dtcell,dtcell2
      real*8  rhoij,rhominij
      real*8  rhomin
      real*8  a,b,c
      integer i

      rhomin = rho_divu_ceiling

      dt = 1.0d20
      do i= 0,nx-1
          dtcell = dt
          if (divu_ceiling_flag.eq.1) then
            if (divu(i).gt.0.0d0) then
              rho = scal(i,Density)
              if (rho .gt. rhomin) then
                dtcell = divu_dt_factor*(1.0D0-rhomin/rho)/divu(i)
              else
                dtcell = divu_dt_factor*0.5d0/divu(i)
              endif

              if (dsdt(i) .gt. 1.0d-20) then
                if (abs(rho).gt.rhomin) then
                  rhominij = rhomin
                else
                  rhominij = .9*abs(rho) 
                endif
                rhoij = abs(rho)
                a = rhoij*dsdt(i)*0.50D0
                b = rhoij*divu(i)
                c = rhominij - rhoij
                dtcell2 = 2.0D0*c/(-b-sqrt(b**2-4.0d0*a*c))
                dtcell2 = divu_dt_factor*dtcell2
                dtcell = min(dtcell,dtcell2)
              endif
            endif
            if(dtcell .le. 0.0d0)then
              print *, 'CHECK_DIVU_DT:warning: aha'
            endif
          else 
            print *,'EST_DIVU_DT only coded for flag == 1'
            stop
          endif
          dt = min(dtcell,dt)
      enddo

      end
