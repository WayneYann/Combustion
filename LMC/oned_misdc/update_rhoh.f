      subroutine update_rhoh(scal_old,scal_new,aofs,alpha,
     &                       beta,dRhs,Rhs,dx,dt,be_cn_theta,time)
      implicit none
      include 'spec.h'
      real*8 scal_old(-1:nx  ,nscal)
      real*8 scal_new(-1:nx  ,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8    alpha(0 :nx-1)
      real*8     beta(-1:nx  ,nscal)
      real*8     dRhs(0 :nx-1)
      real*8      Rhs(0 :nx-1)
      real*8 dx,dt,be_cn_theta,time

      real*8  visc(0:nx-1)
      real*8  h_hi,h_lo,h_mid
      real*8  flux_lo,flux_hi
      real*8  dxsqinv
      real*8  beta_lo,beta_hi
      real*8  visc_term, RWRK
      integer i,n,is, IWRK

      real*8 cp,T,H,rho
      real*8 Y(Nspec)


      if (use_temp_eqn .or. use_strang) then

         call get_rhoh_visc_terms(scal_old,beta,visc,dx,time)
c     FIXME: Add NULN terms
         if (LeEQ1 .ne. 1) then
            print *,'update_rhoh does yet support non-unity Le'
            stop
         endif

         do i = 0,nx-1
            visc_term = dt*(1.d0 - be_cn_theta)*visc(i)

            scal_new(i,RhoH) = scal_old(i,RhoH) + dt*aofs(i,RhoH)
            Rhs(i) = dRhs(i) + scal_new(i,RhoH) + visc_term
            alpha(i) = scal_new(i,Density)
         enddo

      else

         call divRhoDHgradY(scal_new,beta,visc,dx,time)

         do i = 0, nx-1
C update rhoh with advection
            scal_new(i,RhoH) = scal_old(i,RhoH) + dt*aofs(i,RhoH)            

C  and set up diffusion solve for Temp            
            T = scal_new(i,Temp)
            rho = 0.d0
            do n=1,Nspec
               rho = rho + scal_new(i,FirstSpec+n-1)
            enddo
            do n=1,Nspec
               Y(n) = scal_new(i,FirstSpec+n-1)/rho
            enddo
            call CKCPBS(T,Y,IWRK,RWRK,cp)
            CALL CKHBMS(T,Y,IWRK,RWRK,H)

            Rhs(i) = dRhs(i) + scal_new(i,RhoH) 
     $           + dt*be_cn_theta*visc(i) - scal_new(i,Density)*H
     $           + scal_new(i,Density)*cp*T
            alpha(i) = scal_new(i,Density)*cp

         enddo

      endif

      end
