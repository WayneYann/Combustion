      subroutine get_temp_visc_terms(scal,beta,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
      real*8 visc(-1:nfine)
      real*8 dx
      integer lo,hi

c     Compute Div(lambda.Grad(T)) + Gamma_i.Grad(Hi)

c     compute Gamma_i.Grad(Hi)
      call gamma_dot_gradh(scal,beta,visc,dx,lo,hi)

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

      subroutine gamma_dot_gradh(scal,beta,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
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

      real*8 sum_gamma_lo,sum_gamma_hi
      real*8 sumRhoY_lo,sumRhoY_hi
      real*8 RhoYe_lo,RhoYe_hi
      real*8 X(Nspec,-1:nfine)
      real*8 scal_X(-2:nfine+1,nscal)

c     Compute Gamma.Grad(hi) terms

      dxsqinv = 1.d0/(dx*dx)

c     Get Hi, Yi, Xi, rho*Xi at cell centers
      do i=lo-1,hi+1
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         call CKHMS(scal(i,Temp),IWRK,RWRK,hm(1,i))
         do n=1,Nspec
            Y(n,i) = scal(i,FirstSpec+n-1)/rho
         enddo

c     convert Y to X
         CALL CKYTX(Y(:,i),IWRK,RWRK,X(:,i))

c     compute rho*X
         do n=1,Nspec
            scal_X(i,FirstSpec+n-1) = scal(i,Density)*X(n,i)
         end do

      enddo

c     Compute differences
      do i=lo,hi
         dv = 0.d0
         sum_gamma_lo = 0.d0
         sum_gamma_hi = 0.d0
         sumRhoY_lo = 0.d0
         sumRhoY_hi = 0.d0
         do n=1,Nspec
            is = FirstSpec + n - 1
            if (coef_avg_harm.eq.1) then
               beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
               beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
            else
               beta_lo = 0.5d0*(beta(i,is) + beta(i-1,is))
               beta_hi = 0.5d0*(beta(i,is) + beta(i+1,is))
            endif

            gamma_lo(i,n) = beta_lo*(X(n,i)-X(n,i-1))
            gamma_hi(i,n) = beta_hi*(X(n,i+1)-X(n,i))

            if (LeEQ1 .eq. 0) then

c              need to correct fluxes so they add to zero on each face
c              build up the sum of species fluxes on lo and hi faces
c              this will be "rho * V_c"
               sum_gamma_lo = sum_gamma_lo + gamma_lo(i,n)
               sum_gamma_hi = sum_gamma_hi + gamma_hi(i,n)
               
c              build up the sum of rho*Y_m
c              this will be the density
               sumRhoY_lo = sumRhoY_lo+0.5d0*(scal(i-1,is)+scal(i,is))
               sumRhoY_hi = sumRhoY_hi+0.5d0*(scal(i,is)+scal(i+1,is))
               
            end if

         end do

         if (LeEQ1 .eq. 0) then
c           correct the fluxes so they add up to zero before computing visc
            do n=1,Nspec
               is = FirstSpec + n - 1

c              compute rho*Y_m on each face
               RhoYe_lo = .5d0*(scal(i-1,is)+scal(i,is))
               RhoYe_hi = .5d0*(scal(i,is)+scal(i+1,is))

c              set flux = flux - (rho*V_c)*(rho*Y_m)/rho
               gamma_lo(i,n) = gamma_lo(i,n) 
     $              - sum_gamma_lo*RhoYe_lo/sumRhoY_lo
               gamma_hi(i,n) = gamma_hi(i,n) 
     $              - sum_gamma_hi*RhoYe_hi/sumRhoY_hi

            end do
         end if

         do n=1,Nspec
            gamma_dot_gradh_lo = gamma_lo(i,n)*(hm(n,i)-hm(n,i-1))
            gamma_dot_gradh_hi = gamma_hi(i,n)*(hm(n,i+1)-hm(n,i))

            dv = dv + (gamma_dot_gradh_hi + gamma_dot_gradh_lo)*0.5d0
         enddo
         
         visc(i) = dv*dxsqinv

        end do
      end
