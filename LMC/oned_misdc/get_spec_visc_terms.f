        subroutine get_spec_visc_terms(nx,scal,beta,visc,dx)
        implicit none
        include 'spec.h'
        integer nx
        real*8 scal(-1:nx  ,*)
        real*8 beta(-1:nx  ,*)
        real*8 visc(0 :nx-1,*)
        real*8 dx

        integer i,n,is,IWRK
        real*8 beta_lo,beta_hi
        real*8 flux_lo,flux_hi
        real*8 dxsqinv,RWRK
        real*8 Y(-1:nx,maxspec)

c     Note, this returns Div(rho.Di.Grad(Yi) + rho.Y.Vcor)

        dxsqinv = 1.d0/(dx*dx)

        do i = 0,nx-1
           do n=1,Nspec
              Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
           enddo
        enddo

        do i = 1,nx

c     Compute Div( rho.Di.Grad(Yi) )
           do n=1,Nspec
              is = FirstSpec + n - 1
c     Harmonic
c                beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
c                beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
c     Arithmetic
              beta_lo = 0.5*(beta(i,is) + beta(i-1,is))
              beta_hi = 0.5*(beta(i,is) + beta(i+1,is))
              
              flux_hi = beta_hi*(Y(i+1,n) - Y(i  ,n)) 
              flux_lo = beta_lo*(Y(i  ,n) - Y(i-1,n)) 
              visc(i,n) =  (flux_hi - flux_lo)*dxsqinv
           enddo

        end do
        end

