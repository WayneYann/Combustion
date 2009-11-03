      subroutine calc_divu(nx,scal,beta,Ydot,divu,dx)
      implicit none
      include 'spec.h'

c     Quantities passed in
      integer nx
      real*8 scal(-1:nx,nscal)
      real*8 beta(-1:nx,nscal)
      real*8 divu(0:nx-1)
      real*8 Ydot(0:nx-1,nspec)
      real*8 dx
      
      real*8 Y(maxspec)
      real*8 hi(maxspec,-1:nx)
      real*8 mwmix(-1:nx)
      real*8 cpmix(-1:nx)
      real*8 ddivu(-1:nx,maxspec)
      
      real*8 RWRK, sum
      integer IWRK,i,n

      do i = -1,nx
         do n = FirstSpec,LastSpec
            Y(n) = scal(i,n) / scal(i,Density)
         enddo
         call CKMMWY(Y,IWRK,RWRK,mwmix(i))
         call CKCPBS(scal(i,Temp),Y,IWRK,RWRK,cpmix(i))
         call CKHMS(scal(i,Temp),IWRK,RWRK,hi(1,i))
      enddo
      
      call get_temp_visc_terms(nx,scal,beta,divu,dx)
      do i = 0,nx-1
         divu(i) = divu(i) / (cpmix(i)*scal(i,Temp))
      enddo
      
      call get_spec_visc_terms(nx,scal,beta,ddivu,dx)
      do i = 0,nx-1
         do n=1,Nspec
            divu(i) = divu(i)
     &           + ddivu(i,n)*invmwt(i)*mwmix(i)
     &           + (hi(n,i)/cpmix(i) - mwmix(i)*invmwt(i))*Ydot(i,n)
         enddo
         divu(i) = divu(i) / scal(i,Density)
      enddo
      end
