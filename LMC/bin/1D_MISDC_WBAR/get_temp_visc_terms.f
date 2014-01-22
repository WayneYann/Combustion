      subroutine get_temp_visc_terms(scal,beta,beta_for_Wbar,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
      real*8 beta_for_Wbar(-1:nfine  ,nscal)
      real*8 visc(-1:nfine)
      real*8 dx
      integer lo,hi

c     Compute Div(lambda.Grad(T)) + Gamma_i.Grad(Hi)

c     compute Gamma_i.Grad(Hi)
      call gamma_dot_gradh(scal,beta,beta_for_Wbar,visc,dx,lo,hi)

c     Add Div( lambda Grad(T) )
      call addDivLambdaGradT(scal,beta,visc,dx,lo,hi)

      end

      subroutine addDivLambdaGradT(scal,beta,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
      real*8 visc(-1:nfine)
      real*8 dx
      integer lo,hi
      
      integer i
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 dxsqinv

      dxsqinv = 1.d0/(dx*dx)
      do i=lo,hi
         if (coef_avg_harm.eq.1) then
            beta_lo = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i-1,Temp))
            beta_hi = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i+1,Temp))
         else
            beta_lo = 0.5*(beta(i,Temp) + beta(i-1,Temp))
            beta_hi = 0.5*(beta(i,Temp) + beta(i+1,Temp))
         endif
         
         flux_hi = beta_hi*(scal(i+1,Temp) - scal(i  ,Temp)) 
         flux_lo = beta_lo*(scal(i  ,Temp) - scal(i-1,Temp)) 
         visc(i) = visc(i) + (flux_hi - flux_lo) * dxsqinv

      enddo

      end

      subroutine gamma_dot_gradh(scal,beta,beta_for_Wbar,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
      real*8 beta_for_Wbar(-1:nfine  ,nscal)
      real*8 visc(-1:nfine)
      real*8 dx
      integer lo,hi
      
      integer i,n,is,IWRK
      real*8 beta_lo,beta_hi
      real*8 gamma_dot_gradh_lo,gamma_dot_gradh_hi
      real*8 dxsqinv,RWRK,rho,dv
      real*8 hm(Nspec,-1:nfine)
      real*8 Y(Nspec,-1:nfine)

      real*8 gamma_lo(0:nfine-1,Nspec)
      real*8 gamma_hi(0:nfine-1,Nspec)
      real*8 gamma_Wbar_lo(0:nfine-1,Nspec)
      real*8 gamma_Wbar_hi(0:nfine-1,Nspec)

c     Compute Gamma.Grad(hi) terms
      call get_spec_visc_terms(scal,beta,gamma_lo,gamma_hi,lo,hi)
      call get_spec_visc_terms_Wbar(scal,beta_for_Wbar,
     &                              gamma_Wbar_lo,gamma_Wbar_hi,lo,hi)

c     add Wbar part
      gamma_lo = gamma_lo + gamma_Wbar_lo
      gamma_hi = gamma_hi + gamma_Wbar_hi

      call adjust_spec_diffusion_fluxes(scal,gamma_lo,gamma_hi,lo,hi)

c     Get Hi at cell centers
      do i=lo-1,hi+1
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         call CKHMS(scal(i,Temp),IWRK,RWRK,hm(1,i))
      enddo

      do i=lo,hi
         dv = 0.d0
         do n=1,Nspec
            gamma_dot_gradh_lo = gamma_lo(i,n)*(hm(n,i)-hm(n,i-1))
            gamma_dot_gradh_hi = gamma_hi(i,n)*(hm(n,i+1)-hm(n,i))

            dv = dv + (gamma_dot_gradh_hi + gamma_dot_gradh_lo)*0.5d0
         enddo
         visc(i) = dv/(dx*dx)
        end do
      end
