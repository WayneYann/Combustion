      subroutine calc_divu(scal,beta,I_R,divu,dx,time)
      implicit none
      include 'spec.h'

c     Quantities passed in
      real*8 scal(-1:nx,nscal)
      real*8 beta(-1:nx,nscal)
      real*8 divu(0:nx-1)
      real*8 I_R(0:nx-1,0:maxspec)
      real*8 dx, time
      
      real*8 Y(maxspec)
      real*8 cpmix,mwmix
      real*8 ddivu(0:nx-1,maxspec)
      
      real*8 RWRK,rho,T
      integer IWRK,i,n
      real*8 divu_max,marc

      call get_temp_visc_terms(scal,beta,divu,dx,time)
      call get_spec_visc_terms(scal,beta,ddivu,dx,time)

      marc=0.d0
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

         divu(i) = (divu(i)/(rho*cpmix) + I_R(i,0)) / T
         marc = MAX(ABS(divu(i)),marc)
         do n=1,Nspec
            divu(i) = divu(i)
     &           + (ddivu(i,n) + I_R(i,n))*invmwt(n)*mwmix/rho
         enddo
       enddo

      divu_max = ABS(divu(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu(i)))
      enddo
      print *,'*********** DIVU norm = ',divu_max,marc
c      print *,'*********** DIVU norm = ',divu_max

      end
