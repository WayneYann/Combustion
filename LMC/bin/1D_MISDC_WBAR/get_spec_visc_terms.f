      subroutine get_spec_visc_terms(scal,beta,visc,gamma_lo,
     &                               gamma_hi,dx,lo,hi)

      implicit none
      include 'spec.h'
      real*8         scal(-2:nfine+1,nscal)
      real*8         beta(-1:nfine  ,nscal)
      real*8         visc(-1:nfine  ,Nspec)
      real*8 gamma_lo( 0:nfine-1,Nspec)
      real*8 gamma_hi( 0:nfine-1,Nspec)
      real*8 dx
      integer lo,hi
      
      integer i,n,is,IWRK
      real*8 beta_lo,beta_hi,RWRK
      real*8 dxsqinv
      real*8 Y(-1:nfine,Nspec), sum_gamma_lo, sum_gamma_hi, sumRhoX_lo, sumRhoX_hi
      real*8 RhoXe_lo, RhoXe_hi
      real*8 X(-1:nfine,Nspec)
      real*8 scal_X(-2:nfine+1,nscal)

      dxsqinv = 1.d0/(dx*dx)

      do i=lo-1,hi+1

c     compute Y = rho*Y / rho
         do n=1,Nspec
            Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo

c     convert Y to X
         CALL CKYTX(Y(i,:),IWRK,RWRK,X(i,:))

c     compute rho*X
         do n=1,Nspec
            scal_X(i,FirstSpec+n-1) = scal(i,Density)*X(i,n)
         end do

      enddo

      do i=lo,hi

         sum_gamma_lo = 0.d0
         sum_gamma_hi = 0.d0
         sumRhoX_lo = 0.d0
         sumRhoX_hi = 0.d0

         do n=1,Nspec
            is = FirstSpec + n - 1

c     compute beta on edges
            if (coef_avg_harm.eq.1) then
               beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
               beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
            else
               beta_lo = 0.5d0*(beta(i,is) + beta(i-1,is))
               beta_hi = 0.5d0*(beta(i,is) + beta(i+1,is))
            endif

c     compute gamma
            gamma_hi(i,n) = beta_hi*(X(i+1,n) - X(i  ,n)) 
            gamma_lo(i,n) = beta_lo*(X(i  ,n) - X(i-1,n)) 
 
c     compute div(gamma).  If non-unity Le we overwrite this later
            visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

            if (LeEQ1 .eq. 0) then

c              need to correct fluxes so they add to zero on each face
c              build up the sum of species fluxes on lo and hi faces
c              this will be "rho * V_c"
               sum_gamma_lo = sum_gamma_lo + gamma_lo(i,n)
               sum_gamma_hi = sum_gamma_hi + gamma_hi(i,n)
               
c              build up the sum of rho*X_m
c              this will be the density
               sumRhoX_lo = sumRhoX_lo+0.5d0*(scal_X(i-1,is)+scal_X(i,is))
               sumRhoX_hi = sumRhoX_hi+0.5d0*(scal_X(i,is)+scal_X(i+1,is))
               
            end if

         enddo

         if (LeEQ1 .eq. 0) then
c           correct the fluxes so they add up to zero before computing visc
            do n=1,Nspec
               is = FirstSpec + n - 1

c              compute rho*X_m on each face
               RhoXe_lo = .5d0*(scal_X(i-1,is)+scal_X(i,is))
               RhoXe_hi = .5d0*(scal_X(i,is)+scal_X(i+1,is))

c              set flux = flux - (rho*V_c)*(rho*X_m)/rho
               gamma_lo(i,n) = gamma_lo(i,n) 
     $              - sum_gamma_lo*RhoXe_lo/sumRhoX_lo
               gamma_hi(i,n) = gamma_hi(i,n) 
     $              - sum_gamma_hi*RhoXe_hi/sumRhoX_hi

c              compute div(gamma)
               visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

            end do
         end if

      end do
        
      end


      subroutine get_spec_visc_terms_Wbar(scal,beta_for_Wbar,visc,
     &                                    gamma_Wbar_lo,gamma_Wbar_hi,
     &                                    dx,lo,hi)

      implicit none
      include 'spec.h'
      real*8          scal(-2:nfine+1,nscal)
      real*8 beta_for_Wbar(-1:nfine  ,nscal)
      real*8          visc(-1:nfine  ,Nspec)
      real*8 gamma_Wbar_lo( 0:nfine-1,Nspec)
      real*8 gamma_Wbar_hi( 0:nfine-1,Nspec)
      real*8 dx
      integer lo,hi
      
      integer i,n,is,IWRK
      real*8 beta_for_Wbar_lo,beta_for_Wbar_hi,RWRK
      real*8 dxsqinv
      real*8 Y(-1:nfine,Nspec)
      real*8 Wbar(-2:nfine+1)

      dxsqinv = 1.d0/(dx*dx)

      do i=lo-1,hi+1

c     compute Y = rho*Y / rho
         do n=1,Nspec
            Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo
            
c     convert Y to Wbar
         CALL CKMMWY(Y(i,:),IWRK,RWRK,Wbar(i))

      enddo

      do i=lo,hi

         do n=1,Nspec
            is = FirstSpec + n - 1

c     compute beta on edges
            if (coef_avg_harm.eq.1) then
               beta_for_Wbar_lo = 2.d0 / (1.d0/beta_for_Wbar(i,is)+1.d0/beta_for_Wbar(i-1,is))
               beta_for_Wbar_hi = 2.d0 / (1.d0/beta_for_Wbar(i,is)+1.d0/beta_for_Wbar(i+1,is))
            else
               beta_for_Wbar_lo = 0.5d0*(beta_for_Wbar(i,is) + beta_for_Wbar(i-1,is))
               beta_for_Wbar_hi = 0.5d0*(beta_for_Wbar(i,is) + beta_for_Wbar(i+1,is))
            endif

c     compute gamma
            gamma_Wbar_hi(i,n) = beta_for_Wbar_hi*(Wbar(i+1) - Wbar(i  )) 
            gamma_Wbar_lo(i,n) = beta_for_Wbar_lo*(Wbar(i  ) - Wbar(i-1)) 
 
c     compute div(gamma).
c     no need to conservatively correct these
c     we will correct beta grad X after the species diffusion solve
            visc(i,n) = (gamma_Wbar_hi(i,n)-gamma_Wbar_lo(i,n))*dxsqinv

         enddo

      end do

      end


      subroutine get_spec_visc_terms_Y_and_Wbar(scal,beta_for_Y,visc,
     &                                          gamma_Wbar_lo,
     &                                          gamma_Wbar_hi,
     &     dx,lo,hi)

c     compute 
c     gamma_m = beta_for_y grad Y + gamma_Wbar
c     conservatively correct this, then set visc = (1/dxsq)*div(gamma_m)

      implicit none
      include 'spec.h'
      real*8          scal(-2:nfine+1,nscal)
      real*8    beta_for_Y(-1:nfine  ,nscal)
      real*8          visc(-1:nfine  ,Nspec)
      real*8 gamma_Wbar_lo( 0:nfine-1,Nspec)
      real*8 gamma_Wbar_hi( 0:nfine-1,Nspec)
      real*8 dx
      integer lo,hi
      
      integer i,n,is,IWRK
      real*8 beta_for_Y_lo,beta_for_Y_hi,RWRK
      real*8 dxsqinv
      real*8 Y(-1:nfine,Nspec), sum_gamma_lo, sum_gamma_hi, sumRhoX_lo, sumRhoX_hi
      real*8 RhoXe_lo, RhoXe_hi
      real*8 X(-1:nfine,Nspec)
      real*8 scal_X(-2:nfine+1,nscal)
      real*8 gamma_lo( 0:nfine-1,Nspec)
      real*8 gamma_hi( 0:nfine-1,Nspec)

      dxsqinv = 1.d0/(dx*dx)

      do i=lo-1,hi+1

c     compute Y = rho*Y / rho
         do n=1,Nspec
            Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo

c     convert Y to X
         CALL CKYTX(Y(i,:),IWRK,RWRK,X(i,:))

c     compute rho*X
         do n=1,Nspec
            scal_X(i,FirstSpec+n-1) = scal(i,Density)*X(i,n)
         end do

      enddo

      do i=lo,hi

         sum_gamma_lo = 0.d0
         sum_gamma_hi = 0.d0
         sumRhoX_lo = 0.d0
         sumRhoX_hi = 0.d0

         do n=1,Nspec
            is = FirstSpec + n - 1

c     compute beta on edges
            if (coef_avg_harm.eq.1) then
               beta_for_Y_lo = 2.d0 / (1.d0/beta_for_Y(i,is)+1.d0/beta_for_Y(i-1,is))
               beta_for_Y_hi = 2.d0 / (1.d0/beta_for_Y(i,is)+1.d0/beta_for_Y(i+1,is))
            else
               beta_for_Y_lo = 0.5d0*(beta_for_Y(i,is) + beta_for_Y(i-1,is))
               beta_for_Y_hi = 0.5d0*(beta_for_Y(i,is) + beta_for_Y(i+1,is))
            endif

c     compute gamma
            gamma_hi(i,n) = beta_for_Y_hi*(Y(i+1,n) - Y(i  ,n)) 
            gamma_lo(i,n) = beta_for_Y_lo*(Y(i  ,n) - Y(i-1,n))  

            gamma_hi(i,n) = gamma_hi(i,n) + gamma_Wbar_hi(i,n)
            gamma_lo(i,n) = gamma_lo(i,n) + gamma_Wbar_lo(i,n)

c     compute div(gamma).  If non-unity Le we overwrite this later
            visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

            if (LeEQ1 .eq. 0) then

c              need to correct fluxes so they add to zero on each face
c              build up the sum of species fluxes on lo and hi faces
c              this will be "rho * V_c"
               sum_gamma_lo = sum_gamma_lo + gamma_lo(i,n)
               sum_gamma_hi = sum_gamma_hi + gamma_hi(i,n)
               
c              build up the sum of rho*X_m
c              this will be the density
               sumRhoX_lo = sumRhoX_lo+0.5d0*(scal_X(i-1,is)+scal_X(i,is))
               sumRhoX_hi = sumRhoX_hi+0.5d0*(scal_X(i,is)+scal_X(i+1,is))
               
            end if

         enddo

         if (LeEQ1 .eq. 0) then
c           correct the fluxes so they add up to zero before computing visc
            do n=1,Nspec
               is = FirstSpec + n - 1

c              compute rho*X_m on each face
               RhoXe_lo = .5d0*(scal_X(i-1,is)+scal_X(i,is))
               RhoXe_hi = .5d0*(scal_X(i,is)+scal_X(i+1,is))

c              set flux = flux - (rho*V_c)*(rho*X_m)/rho
               gamma_lo(i,n) = gamma_lo(i,n) 
     $              - sum_gamma_lo*RhoXe_lo/sumRhoX_lo
               gamma_hi(i,n) = gamma_hi(i,n) 
     $              - sum_gamma_hi*RhoXe_hi/sumRhoX_hi

c              compute div(gamma)
               visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

            end do
         end if

      end do
        
      end

