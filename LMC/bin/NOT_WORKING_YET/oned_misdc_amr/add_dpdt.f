      subroutine add_dpdt(scal,pthermo,divu,umac,dx,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8    scal(-2:nx+1,nscal)
      real*8 pthermo(-2:nx+1)
      real*8    divu(0 :nx-1)
      real*8    umac(0 :nx  )
      real*8 Y(Nspec)
      real*8 dx
      real*8 dt
      integer lo,hi,bc(2)
      
      real*8 uadv,p_lo,p_hi
      real*8 ugradp
      real*8 denom
      real*8 gamma_inv, cp, Runiv, mwmix, RWRK, dummy
      real*8 dpdt
      real*8 dpdt_max
      integer i,n
      integer ispec,IWRK
      
      dpdt_max = 0.d0
      do i=lo,hi
         uadv = 0.5d0 * (umac(i) + umac(i+1))
         if (umac(i) .ge. 0.d0) then
            p_lo = pthermo(i-1)
         else
            p_lo = pthermo(i)
         endif
         if (umac(i+1) .ge. 0.d0) then
            p_hi = pthermo(i)
         else
            p_hi = pthermo(i+1)
         endif
         ugradp = uadv * (p_hi - p_lo) / dx 
         dpdt = (pthermo(i) - Pcgs) / dt
         dpdt = dpdt - ugradp        

c     compute gamma here
         do n = 1,Nspec
            ispec = FirstSpec + n - 1
            Y(n) = scal(i,ispec)/scal(i,Density)
         enddo
C compute (Y_i/mw_i)^-1 = mean molecular weight
         call CKMMWY(Y,IWRK,RWRK,mwmix)
         call CKRP(IWRK,RWRK,Runiv,dummy,dummy) 
         call CKCPBS(scal(i,Temp),Y,IWRK,RWRK,cp)
         gamma_inv = (cp - Runiv/mwmix)/cp

C         dpdt = dpdt / 1.4d0
         dpdt = dpdt * gamma_inv
         denom = MIN(pthermo(i),Pcgs)
C         denom = pthermo(i)
         dpdt = dpdt/denom
         dpdt = dpdt_factor * dpdt
         divu(i) = divu(i) + dpdt
         dpdt_max = MAX(dpdt_max,ABS(dpdt))
      end do
      
      write(6,1000) dpdt_max
 1000 format(' DPDT norm     = ',f14.4)
      
      end


      subroutine add_dpdt_nodal(scal,pthermo,divu,vel,dx,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8    scal(-2:nx+1,nscal)
      real*8 pthermo(-2:nx+1)
      real*8    divu(0 :nx-1)
      real*8     vel(-2:nx+1)
      real*8 Y(Nspec)
      real*8 dx
      real*8 dt
      integer lo,hi,bc(2)
      
      real*8 uadv,p_lo,p_hi
      real*8 ugradp
      real*8 denom
      real*8 gamma_inv, cp, Runiv, mwmix, RWRK, dummy
      real*8 dpdt
      real*8 dpdt_max
      integer i,n
      integer ispec,IWRK
      
      dpdt_max = 0.d0
      do i=lo,hi
         uadv = vel(i)

         if (i .eq. lo .and. bc(1) .eq. 1) then
!     inflow
            p_lo = pthermo(i)
            p_hi = pthermo(i+1)
         else if (i .eq. hi .and. bc(2) .eq. 2) then
!     outflow
            p_lo = pthermo(i-1)
            p_hi = pthermo(i)
         else
            if (uadv .ge. 0) then
               p_lo = pthermo(i-1)
               p_hi = pthermo(i)
            else                  
               p_lo = pthermo(i)
               p_hi = pthermo(i+1)
            end if
         end if

         ugradp = uadv * (p_hi - p_lo) / dx 
         dpdt = (pthermo(i) - Pcgs) / dt
         dpdt = dpdt - ugradp        

c     compute gamma here
         do n = 1,Nspec
            ispec = FirstSpec + n - 1
            Y(n) = scal(i,ispec)/scal(i,Density)
         enddo
C compute (Y_i/mw_i)^-1 = mean molecular weight
         call CKMMWY(Y,IWRK,RWRK,mwmix)
         call CKRP(IWRK,RWRK,Runiv,dummy,dummy) 
         call CKCPBS(scal(i,Temp),Y,IWRK,RWRK,cp)
         gamma_inv = (cp - Runiv/mwmix)/cp

C         dpdt = dpdt / 1.4d0
         dpdt = dpdt * gamma_inv
         denom = MIN(pthermo(i),Pcgs)
C         denom = pthermo(i)
         dpdt = dpdt/denom
         dpdt = dpdt_factor * dpdt
         divu(i) = divu(i) + dpdt
         dpdt_max = MAX(dpdt_max,ABS(dpdt))
      end do
      
      write(6,1000) dpdt_max
 1000 format(' DPDT norm     = ',f14.4)
      
      end
