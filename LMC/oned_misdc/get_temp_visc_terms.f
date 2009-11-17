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
        real*8 dxsqinv,RWRK,rho
        real*8 hi(maxspec,-1:nx)
        real*8 Y(maxspec,-1:nx)

        real*8 cpi(9)

c     Note, this returns Div(lambda.Grad(T)) + rho.D.Grad(Hi).Grad(Yi)

        dxsqinv = 1.d0/(dx*dx)

c     Compute Hi, Yi
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

        do i = 0,nx-1

c     Compute Div( lambda Grad(T) )
           if (coef_avg_harm.eq.1) then
              beta_lo = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i-1,Temp))
              beta_hi = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i+1,Temp))
           else
              beta_lo = 0.5*(beta(i,Temp) + beta(i-1,Temp))
              beta_hi = 0.5*(beta(i,Temp) + beta(i+1,Temp))
           endif

           flux_hi = beta_hi*(scal(i+1,Temp) - scal(i  ,Temp)) 
           flux_lo = beta_lo*(scal(i  ,Temp) - scal(i-1,Temp)) 
           visc(i) =  flux_hi - flux_lo

c     Add rhoD Grad(Yi).Grad(hi) terms
           do n=1,Nspec
              is = FirstSpec + n - 1
              if (coef_avg_harm.eq.1) then
                 beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
                 beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
              else
                 beta_lo = 0.5*(beta(i,is) + beta(i-1,is))
                 beta_hi = 0.5*(beta(i,is) + beta(i+1,is))
              endif

              flux_lo = beta_lo*(Y(n,i)-Y(n,i-1))*(hi(n,i)-hi(n,i-1))
              flux_hi = beta_hi*(Y(n,i+1)-Y(n,i))*(hi(n,i+1)-hi(n,i))

              visc(i) = visc(i) + flux_hi - flux_lo
              
           enddo

           visc(i) = visc(i) * dxsqinv

        end do
        end

