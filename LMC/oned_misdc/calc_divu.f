      subroutine calc_divu(nx,scal,beta,Ydot,divu,dx,time)
      implicit none
      include 'spec.h'

c     Quantities passed in
      integer nx
      real*8 scal(-1:nx,nscal)
      real*8 beta(-1:nx,nscal)
      real*8 divu(0:nx-1)
      real*8 Ydot(0:nx-1,nspec)
      real*8 dx, time
      
      real*8 Y(maxspec)
      real*8 hi(maxspec,-1:nx)
      real*8 cpmix(-1:nx)
      real*8 ddivu(0:nx-1,maxspec)
      real*8 mwmix(-1:nx)
      
      real*8 RWRK, sum, vel
      integer IWRK,i,n,is,do_vel

      do i = -1,nx
         do n = 1,Nspec
            is = FirstSpec + n - 1
            Y(n) = scal(i,is) / scal(i,Density)
         enddo
         call CKMMWY(Y,IWRK,RWRK,mwmix(i))
         call CKCPBS(scal(i,Temp),Y,IWRK,RWRK,cpmix(i))
         call CKHMS(scal(i,Temp),IWRK,RWRK,hi(1,i))
      enddo
      
      do_vel = 0
      call applybc(nx,vel,scal,dx,time,do_vel)
      call get_temp_visc_terms(nx,scal,beta,divu,dx)
      do i = 0,nx-1
         divu(i) = divu(i) / (cpmix(i)*scal(i,Temp))
      enddo
      call get_spec_visc_terms(nx,scal,beta,ddivu,dx)
      do i = 0,nx-1
         do n=1,Nspec
            divu(i) = divu(i)
     &           + ddivu(i,n)*invmwt(n)*mwmix(i)
     &           - (hi(n,i)/(cpmix(i)*scal(i,Temp))
     &           -   mwmix(i)*invmwt(n))*Ydot(i,n)
         enddo
      enddo
      end
