      subroutine divRhoDHgradY(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx, time
      
      integer i,is,n,IWRK
      real*8 dxsqinv, rho,dv,RWRK
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 hi(maxspec,-1:nx)
      real*8 Y(maxspec,-1:nx)


c     Compute Div( rhoD h_m Grad(Yi))

      dxsqinv = 1.d0/(dx*dx)
c     Get Hi, Yi at cell centers
      do i = -1,nx
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         call CKHMS(scal(i,Temp),IWRK,RWRK,hi(1,i))
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
         enddo
      enddo

c     Compute differences
      do i = 0,nx-1
         dv = 0.d0
         do n=1,Nspec
            is = FirstSpec + n - 1
            if (coef_avg_harm.eq.1) then
               beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
               beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
            else
               beta_lo = 0.5d0*(beta(i,is) + beta(i-1,is))
               beta_hi = 0.5d0*(beta(i,is) + beta(i+1,is))
            endif

C     get h_m at edges
C CEG:: maybe it would be a better idea to just use Tmid
            flux_lo = 0.5d0*(hi(n,i)+hi(n,i-1))
            flux_hi = 0.5d0*(hi(n,i+1)+hi(n,i))   

            flux_lo = beta_lo*(Y(n,i)-Y(n,i-1))*flux_lo
            flux_hi = beta_hi*(Y(n,i+1)-Y(n,i))*flux_hi
            
            dv = dv + (flux_hi - flux_lo)
         enddo
         
         visc(i) = dv*dxsqinv
      end do

      end


      subroutine divBetaHgradY(scal,beta,visc,dx,time)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx  ,*)
      real*8 beta(-1:nx  ,*)
      real*8 visc(0 :nx-1)
      real*8 dx, time
      
      integer i,is,n,IWRK
      real*8 dxsqinv, rho,dv,RWRK
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 hi(maxspec,-1:nx)
      real*8 Y(maxspec,-1:nx)
      real*8 rdmlocp(-1:nx, Nspec)

      real*8 Tloc

c     Compute Div( h_m (rhoD - lamba/cp) Grad(Yi))

      dxsqinv = 1.d0/(dx*dx)
c     Get Hi, Yi at cell centers
      do i = -1,nx
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
         enddo
c         call CKRHOPY(rho,Pcgs,Y,IWRK,RWRK,Tloc)
         print*,"ABORT: AJN Removed 3/15/11 to get GRI30 to work"
         stop
C         call CKHMS(scal(i,Temp),IWRK,RWRK,hi(1,i))
         call CKHMS(Tloc,IWRK,RWRK,hi(1,i))
         do n=1,Nspec
C            Y(n,i) = scal(i,FirstSpec+n-1)/rho
            rdmlocp(i,n) = beta(i,FirstSpec+n-1) - beta(i,RhoH)
         enddo
      enddo

c     Compute differences
      do i = 0,nx-1
         dv = 0.d0
         do n=1,Nspec
            is = n 
            if (coef_avg_harm.eq.1) then
               beta_lo= 2.d0 / (1.d0/rdmlocp(i,is)+1.d0/rdmlocp(i-1,is))
               beta_hi= 2.d0 / (1.d0/rdmlocp(i,is)+1.d0/rdmlocp(i+1,is))
            else
               beta_lo = 0.5d0*(rdmlocp(i,is) + rdmlocp(i-1,is))
               beta_hi = 0.5d0*(rdmlocp(i,is) + rdmlocp(i+1,is))
            endif

C     get h_m at edges
C CEG:: maybe it would be a better idea to just use Tmid
            flux_lo = 0.5d0*(hi(n,i)+hi(n,i-1))
            flux_hi = 0.5d0*(hi(n,i+1)+hi(n,i))   

            flux_lo = beta_lo*(Y(n,i)-Y(n,i-1))*flux_lo
            flux_hi = beta_hi*(Y(n,i+1)-Y(n,i))*flux_hi
            
            dv = dv + (flux_hi - flux_lo)
         enddo
         
         visc(i) = dv*dxsqinv
      end do

C     debugging FIXME
C$$$ 1007 FORMAT((I5,1X),12(E22.15,1X))      
C$$$      open(UNIT=11, FILE='nuln.dat', STATUS = 'REPLACE')
C$$$      do i=0,nx-1
C$$$         write(11,1007) i, (rdmlocp(i,n),n = 1,nspec),beta(i,RhoH),
C$$$     &        beta(i,FirstSpec),beta(i,RhoH)-beta(i,FirstSpec)
C$$$      enddo
C$$$      close(11)
C         stop
CCCCCCCCCCCCC

      end
