      subroutine update_spec(nx,scal_old,scal_new,
     $                       aofs,alpha,beta,Rhs,dx,dt,be_cn_theta)
      implicit none
      include 'spec.h'
      integer nx
      real*8 scal_old(-1:nx  ,nscal)
      real*8 scal_new(-1:nx  ,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8    alpha(0 :nx-1)
      real*8     beta(-1:nx  ,nscal)
      real*8      Rhs(0 :nx-1,nscal)
      real*8 dx
      real*8 dt
      real*8 be_cn_theta
      
      real*8  dth,dxsqinv
      real*8  beta_lo,beta_hi
      real*8  visc_term
      integer i,n,is
      real*8 flux_lo(maxspec), flux_hi(maxspec)
      real*8 RhoYe_lo(maxspec), RhoYe_hi(maxspec)
      real*8 Y_L, Y_C, Y_R, sum_lo, sum_hi, sumRhoYe_lo, sumRhoYe_hi

      dth = be_cn_theta * dt
      dxsqinv = 1.d0/(dx*dx)
      
c*************************************************************************
c       Create RHS = time n diffusive term.
c*************************************************************************
c     Note, this returns Div(rho.Di.Grad(Yi) + rho.Y.Vcor)


c     Note, some double calc here for cell-based indexing 
      do i = 0,nx-1
         
c     Compute Div( rho.Di.Grad(Yi) ) but ensure sum spec fluxes = 0
         sum_lo = 0.d0
         sum_hi = 0.d0
         sumRhoYe_lo = 0
         sumRhoYe_hi = 0
         do n=1,Nspec
            is = FirstSpec + n - 1
            if (coef_avg_harm.eq.1) then
               beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
               beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
            else
               beta_lo = 0.5*(beta(i,is) + beta(i-1,is))
               beta_hi = 0.5*(beta(i,is) + beta(i+1,is))
            endif

            Y_L = scal_old(i-1,is) / scal_old(i-1,Density)
            Y_C = scal_old(i  ,is) / scal_old(i  ,Density)
            Y_R = scal_old(i+1,is) / scal_old(i+1,Density)
            
            flux_hi(n) = beta_hi*(Y_R - Y_C) 
            flux_lo(n) = beta_lo*(Y_C - Y_L)
            
            if (LeEQ1 .eq. 0) then
               sum_lo = sum_lo + flux_lo(n)
               sum_hi = sum_hi + flux_hi(n)
               
               RhoYe_lo(n) = .5d0*(scal_old(i-1,is)+scal_old(i,is))
               RhoYe_hi(n) = .5d0*(scal_old(i,is)+scal_old(i+1,is))

               sumRhoYe_lo = sumRhoYe_lo + RhoYe_lo(n)
               sumRhoYe_hi = sumRhoYe_hi + RhoYe_hi(n)
            endif
         enddo
         if (LeEQ1 .eq. 0) then
            do n=1,Nspec
               is = FirstSpec + n - 1
               flux_lo(n) = flux_lo(n) - sum_lo*RhoYe_lo(n)/sumRhoYe_lo
               flux_hi(n) = flux_hi(n) - sum_hi*RhoYe_hi(n)/sumRhoYe_hi
            enddo
         endif
         do n=1,Nspec
            is = FirstSpec + n - 1
            visc_term = (flux_hi(n) - flux_lo(n))*dxsqinv
            scal_new(i,is) = scal_old(i,is) + dt*aofs(i,is)
            Rhs(i,n) = visc_term * dth + scal_new(i,is)
            alpha(i) = scal_new(i,Density)
         enddo
      enddo

      end


