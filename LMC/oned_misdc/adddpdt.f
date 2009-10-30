        subroutine add_dpdt(nx,rhort,divu,umac,dx,dt)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in
        integer nx
        real*8 rhort(-1:nx)
        real*8  divu(0 :nx-1)
        real*8  umac(0 :nx  )
        real*8 dx
        real*8 dt

c       Local variables
        real*8 uadv,p_lo,p_hi
        real*8 ugradp
        real*8 denom
        real*8 dpdt
        real*8 dpdt_max
        integer i,n
        integer ispec

        dpdt_max = 0.d0
        do i = 0,nx-1
          uadv = 0.5d0 * (umac(i) + umac(i+1))
          if (umac(i) .ge. 0.d0) then
            p_lo = rhort(i-1)
          else
            p_lo = rhort(i)
          endif
          if (umac(i+1) .ge. 0.d0) then
            p_hi = rhort(i)
          else
            p_hi = rhort(i+1)
          endif
          ugradp = uadv * (p_hi - p_lo) / dx 
          dpdt = (rhort(i) - Pamb) / dt
          dpdt = dpdt - ugradp        
          dpdt = dpdt / 1.4d0
          denom = max(epsilon(Pamb)*Pamb,min(rhort(i),Pamb))
          dpdt = dpdt/denom
          dpdt = dpdt_factor * dpdt
          divu(i) = divu(i) + dpdt
          dpdt_max = max(dpdt_max,abs(dpdt))
        end do

        write(6,1000) dpdt_max
1000    format('DPDT norm     = ',f14.4)

        end
