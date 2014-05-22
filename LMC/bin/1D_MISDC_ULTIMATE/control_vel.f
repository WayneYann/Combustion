

      subroutine activecontrol(coft,time,dt,myproc,step,restart,usetemp)

      implicit none
      include 'spec.h'

      real*8 coft,time,dt,vslope,slocal,V_new,dVmax,dVmin
      integer myproc,step,ierr,restart,usetemp
      real*8 r1,r2,r3,r4,r5,r6,r7
      real*8 b,c,d
      real*8 sest_save
      integer i1
      logical found_it

      if (step.lt.0) return

      if (restart.ne.0) then
         open(13,file='AC_HISTORY',form='formatted',
     &        status='old',iostat=ierr)
         found_it = .false.
         if (ierr .eq. 0) then
            if (myproc.eq.0) then
               print*, 'Setting active control from history file ...'
            endif
            rewind(13)
            do
c
c                 This read(13) must correspond to the below write(13)
c
               read(13,1000,iostat=ierr) i1,r1,r2,r3,r4,r5,r6,r7
               if (ierr.ne.0) goto 100
               if (i1.eq.step) then

                  found_it = .true.
                  V_in = r2
                  tbase_control = r3
                  zbase_control = r4
                  dV_control = r5
                  sest = r6
                  coft_old = r7

               endif
            enddo

         else

            if (myproc.eq.0) then
               open(13,file='AC_HISTORY',form='formatted', status='new')
            endif

         endif

 100     if (found_it .eqv. .false.) then

            if (myproc.eq.0) then
               print*, 'Setting active control to defaults ...'
            endif

         end if
         close(13)
         return
      end if

      if (time .le. tmax_control) then

         if (usetemp.eq.0 .and. coft_old .lt. 0.d0) coft_old = coft

         zbase_control = zbase_control + V_in*dt + dV_control*dt**2
         V_in_old = V_in
         V_in = V_in + dt*dV_control

         slocal = 0.5d0*(V_in_old + V_in) - (coft - coft_old)/dt
         sest_save = sest
         sest = (1.d0 - corr)*sest + corr*slocal

         V_new = 0.d0
         if (ord_control .eq. 1) then
c           Linear
            vslope = 2*((cfix-coft)/tau_control + sest - V_in)/tau_control
            V_new = V_in + dt*vslope
            
         else if (ord_control .eq. 2) then
c           Quadratic
            vslope = 3.d0*((cfix-coft)/tau_control + sest - V_in)/tau_control
            V_new = V_in + (dt-0.5d0*dt**2/tau_control)*vslope
            
         else if (ord_control .eq. -2) then
c           Quadratic
            c = (-1.5d0/tau_control**2)*((cfix-coft)/tau_control + sest - V_in)
            b = - 2*c*tau_control
            V_new = V_in + dt*(b + dt*c) 

         else if (ord_control .eq. 3) then

c           Cubic, V'(0) = dV, V'(tau)=0
c$$$            d = (-4.d0/tau_control**3)*((cfix-coft)/tau_control
c$$$     &           + sest - V_in - dV_control*tau_control/3.d0)
c$$$            c = (-.5d0/tau_control)*(dV_control + 3*d*tau_control**2)
c$$$            b = - (2.d0*c + 3.d0*d*tau_control)*tau_control
c$$$            V_new = V_in + dt*(b + dt*(c + dt*d))
c$$$            dV_control = b + c*tau_control + d*tau_control**2

c           Cubic, V'(tau)=0, V(0) = Vin, V(tau) = sest
            d = (12/tau_control**3)*((cfix-coft)/tau_control
     &           - (V_in - sest)/3 )
            c = (V_in - sest)/tau_control**2 - 2*d*tau_control
            b = - (2*c*tau_control + 3*d*tau_control**2)
            V_new = V_in + dt*(b + dt*(c + dt*d))
            dV_control = b + c*tau_control + d*tau_control**2

         else if (ord_control .eq. -3) then
c           Cubic, V'(tau)=0, V''(tau)=0
            d = (1.5d0/tau_control**3)*((cfix-coft)/tau_control
     &           + sest - V_in)
            c = - 3*d*tau_control
            b = - (2*c + 3*d*tau_control)*tau_control
            V_new = V_in + dt*(b + dt*(c + dt*d))
            dV_control = - c*tau_control - 2*d*tau_control**2
         endif

         dVmax = changeMax_control
         dVmin = changeMax_control * max(1.d0,V_in)

         V_new = MIN(MAX(V_new,V_in-dVmin),V_in+dVmax)
         V_new = MAX(0.d0,V_new)

         tbase_control = time
         if (ord_control .ne. -3) then
            dV_control = (V_new - V_in)/dt
         endif

      endif

      if (myproc.eq.0) then
         open(13,file='AC_HISTORY',form='formatted',position='append')
         write(13,1000) step,time,V_in,tbase_control,zbase_control,
c     &        dV_control,sest,coft_old,
c     &        ABS(sest-sest_save)/(1.e-8 + 0.5*(sest+sest_save))
     &        dV_control,sest,coft
         close(13)
      endif
 1000 format(i7,8g26.18)

      coft_old = coft

      end


      subroutine control_vel(scal,dx,time,dt,step,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 dx, time, dt
      integer step,lo, hi
 
      real*8 fuelmass, flameloc
      integer i, myproc, restart, usetemp

      fuelmass = 0.d0
      do i=lo,hi
         fuelmass = fuelmass + dx*scal(i,Density)*scal(i,FirstSpec)
      enddo
      flameloc = fuelmass / (scal(lo,Density)*scal(lo,FirstSpec))

      myproc = 0
      restart = 0
      usetemp = 0
      call activecontrol(flameloc,time,dt,myproc,step,restart,usetemp)
      end
