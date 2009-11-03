        subroutine get_temp_visc_terms(nx,scal,beta,visc,dx)
        implicit none
        include 'spec.h'
        integer nx
        real*8 scal(-1:nx  ,*)
        real*8 beta(-1:nx  ,*)
        real*8 visc(0 :nx-1)
        real*8 dx

        integer i,n,is,IWRK
        real*8 beta_lo,beta_hi
        real*8 flux_lo,flux_hi
        real*8 dxsqinv,RWRK
        real*8 hi(maxspec,-1:nx)
        real*8 Y(maxspec,-1:nx)

c     Note, this returns Div(lambda.Grad(T)) + rho.D.Grad(Hi).Grad(Yi)
c     (in particular, does not scale by rho.cp)

        dxsqinv = 1.d0/(dx*dx)

        do i = 0,nx-1
           call CKHMS(scal(i,Temp),IWRK,RWRK,hi(1,i))
           do n=1,Nspec
              Y(n,i) = scal(i,FirstSpec+n-1)/scal(i,Density)
           enddo
        enddo

        do i = 1,nx

c     Compute Div( lambda Grad(T) )

c     Harmonic
c           beta_lo = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i-1,Temp))
c           beta_hi = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i+1,Temp))
c     Arithmetic
           beta_lo = 0.5*(beta(i,Temp) + beta(i-1,Temp))
           beta_hi = 0.5*(beta(i,Temp) + beta(i+1,Temp))

           flux_hi = beta_hi*(scal(i+1,Temp) - scal(i  ,Temp)) 
           flux_lo = beta_lo*(scal(i  ,Temp) - scal(i-1,Temp)) 
           visc(i) =  flux_hi - flux_lo


c     Add rhoD Grad(Yi).Grad(hi) terms
           do n=1,Nspec
              is = FirstSpec + n - 1
c     Harmonic
c              beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
c              beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
c     Arithmetic
              beta_lo = 0.5*(beta(i,is) + beta(i-1,is))
              beta_hi = 0.5*(beta(i,is) + beta(i+1,is))

              flux_lo = beta_lo*(Y(n,i)-Y(n,i-1))*(hi(n,i)-hi(n,i-1))
              flux_hi = beta_hi*(Y(n,i+1)-Y(n,i))*(hi(i+1,n)-hi(i,n))

              visc(i) = visc(i) + flux_hi - flux_lo
              
           enddo

           visc(i) = visc(i)*dxsqinv

        end do
        end

