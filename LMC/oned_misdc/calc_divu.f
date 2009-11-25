      subroutine calc_divu(nx,scal,beta,Ydot,divu,dx,time)
      implicit none
      include 'spec.h'

c     Quantities passed in
      integer nx
      real*8 scal(-1:nx,nscal)
      real*8 beta(-1:nx,nscal)
      real*8 divu(0:nx-1)
      real*8 Ydot(0:nx-1,0:maxspec)
      real*8 dx, time
      
      real*8 Y(maxspec)
      real*8 cpmix,mwmix
      real*8 ddivu(0:nx-1,maxspec)
      
      real*8 RWRK,rho,T
      integer IWRK,i,n

      call get_temp_visc_terms(nx,scal,beta,divu,dx,time)
      call get_spec_visc_terms(nx,scal,beta,ddivu,dx,time)

      do i = 0,nx-1
         rho = 0.d0
         do n = 1,Nspec
            rho = rho + scal(i,FirstSpec + n - 1)
         enddo
         do n = 1,Nspec
            Y(n) = scal(i,FirstSpec + n - 1) / rho
         enddo
         T = scal(i,Temp)
         call CKMMWY(Y,IWRK,RWRK,mwmix)
         call CKCPBS(T,Y,IWRK,RWRK,cpmix)

         divu(i) = (divu(i)/(rho*cpmix) + Ydot(i,0)) / T
         do n=1,Nspec
            divu(i) = divu(i)
     &           + (ddivu(i,n) + Ydot(i,n))*invmwt(n)*mwmix/rho
         enddo
      enddo
      end
