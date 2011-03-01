      subroutine init_soln(S,time,dx)
      implicit none
      include 'spec.h'

      real*8  S(maxscal,0:nx+1), time, dx, x, xl, xh
      integer i, n, Npmf
      real*8 pmfdata(maxspec+3), mole(maxspec), mass(maxspec)

      typVal(Temp)    = 0.d0
      typVal(RhoH)    = 0.d0
      typVal(Density) = 0.d0
      do i=1,nx
         if (probtype.eq.1) then
c            x = problo + (i+0.5d0)*dx - flame_offset
c            call pmf(x,x,pmfdata,Npmf)
            xl = problo + i*dx - flame_offset
            xh = problo + (i+1)*dx - flame_offset
            call pmf(xl,xh,pmfdata,Npmf)
            if (Npmf.ne.Nspec+3) then
               print *,'mismatched pmf'
               stop
            endif
            S(Temp,i) = pmfdata(1)
            do n=1,Nspec
               mole(n) = pmfdata(3+n)
            enddo
         else
            x = problo + (i+0.5d0)*dx
            S(Temp,i) = 298.d0
            do n=1,Nspec
               mole(n) = 0.d0
            enddo
            if (x.lt.0.5d0*(problo+probhi)) then
               mole(1) = 0.3d0
               mole(4) = 0.3d0
               mole(9) = 0.4d0
               S(Temp,i) = 298.d0
            else
               mole(4) = 0.21d0
               mole(9) = 0.79d0
               S(Temp,i) = 600.d0
               S(Temp,i) = 298.d0
            endif
         endif

         call CKXTY(mole,IWRK,RWRK,mass)         
         call CKRHOY(Pcgs,S(Temp,i),mass,IWRK,RWRK,S(Density,i))
         call CKHBMS(S(Temp,i),mass,IWRK,RWRK,S(RhoH,i))

         do n=1,Nspec
            S(FirstSpec+n-1,i) = mass(n) * S(Density,i)
         enddo
         S(RhoH,i) = S(RhoH,i) * S(Density,i)
         
         typVal(Temp)    = MAX(ABS(S(Temp,i)),   typVal(Temp))
         typVal(Density) = MAX(ABS(S(Density,i)),typVal(Density))
         typVal(RhoH   ) = MAX(ABS(S(RhoH,i)),   typVal(RhoH))
      enddo
      do n=1,Nspec
         typVal(FirstSpec+n-1) = typVal(Density)
      enddo
      time = 0.d0

      end

      subroutine apply_bcs(S,time,step)
      implicit none
      include 'spec.h'
      real*8 S(maxscal,0:nx+1)
      real*8 time
      integer step,n
c     Left boundary grow cell
      S(Density,0) = S(Density,1)
      S(Temp,0) = S(Temp,1)
      do n=1,Nspec
         S(FirstSpec+n-1,0) = S(FirstSpec+n-1,1)
      enddo
      S(RhoH,0) = S(RhoH,1)
      
c     Right boundary grow cell
      S(Density,nx+1) = S(Density,nx)
      S(Temp,nx+1) = S(Temp,nx)
      do n=1,Nspec
         S(FirstSpec+n-1,nx+1) = S(FirstSpec+n-1,nx)
      enddo
      S(RhoH,nx+1) = S(RhoH,nx)
      end

