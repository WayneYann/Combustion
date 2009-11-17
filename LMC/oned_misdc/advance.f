      subroutine advance(nx,vel_old,vel_new,scal_old,scal_new,
     $                   Ydot_old,Ydot_new,press_old,press_new,
     $                   divu_old,divu_new,dsdt,beta_old,beta_new,
     $                   intra,dx,dt,time)

      implicit none
      include 'spec.h'
      integer nx
      real*8   vel_new(-1:nx  )
      real*8   vel_old(-1:nx  )
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
      real*8 press_new(0 :nx  )
      real*8 press_old(0 :nx  )
      real*8  Ydot_new(0 :nx-1,nspec)
      real*8  Ydot_old(0 :nx-1,nspec)
      real*8    macvel(0 :nx  )
      real*8   veledge(0 :nx  ,nscal)
      real*8      aofs(0 :nx-1,nscal)
      real*8  divu_old(0 :nx-1)
      real*8  divu_new(0 :nx-1)
      real*8  divu_tmp(0 :nx-1)
      real*8  beta_old(-1:nx,nscal)
      real*8  beta_new(-1:nx,nscal)
      real*8    mu_old(-1:nx)
      real*8    mu_new(-1:nx)
      real*8      dsdt(0 :nx-1)
      real*8        gp(0 :nx-1)
      real*8   rhohalf(0 :nx-1)
      real*8    tforce(0 :nx-1,nscal)
      real*8      visc(0 :nx-1)
      real*8     intra(0 :nx-1,nscal)
      real*8        cp(0 :nx-1)
      real*8 x
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta
      
      real*8    diff_old(0 :nx-1,nscal)
      real*8    diff_hat(0 :nx-1,nscal)
      real*8    diff_new(0 :nx-1,nscal)
      real*8   const_src(0 :nx-1,nscal)
      real*8 lin_src_old(0 :nx-1,nscal)
      real*8 lin_src_new(0 :nx-1,nscal)
      real*8    scal_tmp(0 :nx-1,nscal)
      real*8    sumh
      integer misdc
      
      integer i,n,ispec
      integer iunit
      
      real*8 divu_max
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8 rhocp_old(0:nx-1)
      real*8   pthermo(-1:nx  )
      real*8    Ydot_max, Y(maxspec)
      real*8 RWRK, cpmix, sum
      integer IWRK, is
      
      print *,'advance: at start of time step'
c     
c*****************************************************************
c     Create MAC velocities.
c*****************************************************************
c     
      
      do i = 0,nx-1
         gp(i) = (press_old(i+1) - press_old(i)) / dx
      enddo
      
c      call minmax_scal(nx,scal_old)
      print *,'... predict edge velocities'
      call pre_mac_predict(nx,vel_old,scal_old,gp,
     $                     macvel,dx,dt)
      
      call compute_pthermo(nx,scal_old,pthermo)

      do i = 0,nx-1
         divu_tmp(i) = divu_old(i) + 0.5d0 * dt * dsdt(i)
      enddo
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU norm old = ',divu_max 
      call add_dpdt(nx,pthermo,divu_tmp,macvel,dx,dt)
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU norm new = ',divu_max 
      call macproj(nx,macvel,divu_tmp,dx)

      do n = 1,nscal
         do i = 0,nx-1
            tforce(i,n) = 0.d0
c     Initialized elsewhere, reused from previous timestep where possible
c     FIXME: starting on a bad note
cc            intra(i,n) = 0.d0
c            diff_old(i,n) = 0.d0
c            diff_hat(i,n) = 0.d0
c            diff_new(i,n) = 0.d0
c            aofs(i,n)     = 0.d0
         enddo
      enddo
c
c*****************************************************************
c     
      print *,'... creating the diffusive terms with old data'
      
      call calc_diffusivities(nx,scal_old,beta_old,mu_old,dx,time)
      call get_temp_visc_terms(nx,scal_old,beta_old,diff_old(0,Temp),dx)
      call get_spec_visc_terms(nx,scal_old,beta_old,
     &                         diff_old(0,FirstSpec),dx)
      call get_rhoh_visc_terms(nx,scal_old,beta_old,diff_old(0,RhoH),dx)
      
c*****************************************************************
      
      print *,'... computing aofs with source = diff_old + intra'
      do i = 0,nx-1
         do n = 1,nscal
            tforce(i,n) = diff_old(i,n) + intra(i,n)
         enddo

c     Temp visc terms must be scaled by 1/(rho.cp)
         do n = 1,Nspec
            Y(n) = scal_old(i,FirstSpec+n-1) / scal_old(i,Density)
         enddo
         call CKCPBS(scal_old(i,Temp),Y,IWRK,RWRK,cpmix)
         rhocp_old(i) = cpmix * scal_old(i,Density)
         tforce(i,n) = tforce(i,n)/rhocp_old(i)
      enddo
      
      call scal_aofs(nx,scal_old,macvel,aofs,tforce,dx,dt,time)

      do n = 1,nscal
         do i = 0,nx-1
            tforce(i,n) = intra(i,n)
         enddo
      enddo
      do i = 0,nx-1
         tforce(i,Temp) = tforce(i,Temp)/rhocp_old(i)
      enddo
c*****************************************************************

      print *,'... update rho'
      call update_rho(nx,scal_old,scal_new,aofs,dx,dt)

      print *,'... update to temp to define new diff coeffs'
      be_cn_theta = 1.d0
      call update_temp(nx,scal_old,scal_new,aofs,
     $                 alpha,beta_old,Rhs(0,Temp),dx,dt,be_cn_theta)
      call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,Temp),
     $              dx,dt,Temp,be_cn_theta)

      call get_hmix_given_T_RhoY(nx,scal_new,dx)
      
      call calc_diffusivities(nx,scal_new,beta_new,mu_new,dx,time+dt)

c*****************************************************************

      print *,'... update to spec with new diff. coeffs, do predictor'
      call update_spec(nx,scal_old,scal_new,aofs,alpha,
     $                 beta_old,Rhs(0,FirstSpec),dx,dt,be_cn_theta)
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta)
      enddo

c     FIXME: Adjust spec flux at np1, reset scal_new, compute Le!=1 terms
      call update_rhoh(nx,scal_old,scal_new,aofs,alpha,beta_old,
     &                 Rhs(0,RhoH),dx,dt,be_cn_theta)
      call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta)

c      do i=0,nx-1
c         print *,'enth:',i,scal_new(i,FirstSpec)/scal_new(i,Density),
c     &        scal_old(i,FirstSpec)/scal_old(i,Density)
c      enddo
c      stop
      call rhoh_to_temp(nx,scal_new)

      print *,'... recompute diffusivities, then do corrector'

      call calc_diffusivities(nx,scal_new,beta_new,mu_new,dx,time+dt)

      call update_spec(nx,scal_old,scal_new,aofs,alpha,
     $                 beta_old,Rhs(0,FirstSpec),dx,dt,be_cn_theta)
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(nx,scal_new,alpha,beta_old,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta)
      enddo

c     FIXME: Adjust spec flux at np1, reset scal_new, compute Le!=1 terms
      
      call update_rhoh(nx,scal_old,scal_new,aofs,alpha,beta_old,
     &                 Rhs(0,RhoH),dx,dt,be_cn_theta)
      call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta)
      call rhoh_to_temp(nx,scal_new)


      print *,'... create new visc terms'
      call get_temp_visc_terms(nx,scal_new,beta_new,diff_new(0,Temp),dx)
      call get_spec_visc_terms(nx,scal_new,beta_new,
     &                         diff_new(0,FirstSpec),dx)
      call get_rhoh_visc_terms(nx,scal_new,beta_new,diff_new(0,RhoH),dx)

c*****************************************************************
c       Do chemistry update with const source = adv + diff_new
c*****************************************************************

        print *,'... advancing chem with const source'
        do n = 1,nscal
          do i = 0,nx-1
              const_src(i,n) =     aofs(i,n)
            lin_src_old(i,n) = diff_new(i,n)
            lin_src_new(i,n) = diff_new(i,n)
          enddo
        enddo

        call strang_chem(nx,scal_old,scal_new,
     $                   const_src,lin_src_old,lin_src_new,
     $                   intra,dt)

        print *,'MOVING ON TO MISDC LOOP '

c*****************************************************************
c       Now begin iterations.
c*****************************************************************
   
        do misdc = 1, misdc_iterMAX

c*****************************************************************

           do n=1,nscal
              do i=0,nx-1
                 scal_tmp(i,n) = scal_new(i,n)
              enddo
           enddo

           call calc_diffusivities(nx,scal_new,beta_new,mu_new,
     &          dx,time+dt)

           print *,'... create new diff. terms : diff_hat'
           be_cn_theta = 0.5d0
           call get_rhoh_visc_terms(nx,scal_new,beta_new,
     &                              diff_hat(0,RhoH),dx)
           call get_temp_visc_terms(nx,scal_new,beta_new,
     &                              diff_hat(0,Temp),dx)
           
c*****************************************************************
           
           do n = 1,nscal
              do i = 0,nx-1
                 tforce(i,n) = intra(i,n) + diff_old(i,n)
              enddo
           enddo
           
           print *,'... computing aofs with source = diff_old + intra'
           call scal_aofs(nx,scal_old,macvel,aofs,tforce,dx,dt,time)
           
c*****************************************************************
           
           do n = 1,nscal
              do i = 0,nx-1
                 tforce(i,n) = intra(i,n)
     $                + 0.5d0 * (diff_old(i,n) + diff_hat(i,n))
     $                -          diff_hat(i,n)
              enddo
           enddo
           
           print *,'Updating species,rho with advective + intra terms'
           call update_spec(nx,scal_old,scal_new,aofs,alpha,
     $                      beta_old,Rhs(0,FirstSpec),dx,dt,be_cn_theta)
           do n=1,Nspec
              is = FirstSpec + n - 1
              call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,is),
     $                      dx,dt,is,be_cn_theta)
           enddo

           print *,'... update to rhoH with new diff. coeffs'
           call update_rhoh(nx,scal_old,scal_new,aofs,alpha,beta_old,
     &                      Rhs(0,RhoH),dx,dt,be_cn_theta)
           call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,RhoH),
     $                   dx,dt,RhoH,be_cn_theta)

           print *,'... create new diff. terms for RhoH: diff_new'
           if (be_cn_theta .eq. 1.d0) then
              do i = 0,nx-1
c                 diff_new(i,RhoH) = 
c     $                (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
c     $                -aofs(i,RhoH)-tforce(i,RhoH)
              enddo
           else if (be_cn_theta .eq. 0.5d0) then
              do i = 0,nx-1
c                 diff_new(i,RhoH) = 2.d0 * (
c     $                (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
c     $                -aofs(i,RhoH)-tforce(i,RhoH)
c     $                -0.5d0*diff_old(i,RhoH))
              enddo
           else 
              print *,'OOPS - BOGUS BE_CN_THETA ',be_cn_theta
              stop
           endif
           
           print *,'... create new temp from new rhoH'
           call rhoh_to_temp(nx,scal_new)
           
           print *,'... create new diff. terms for Temp: diff_new'
c           call compute_cp(nx,cp,scal_new)
           do i = 0,nx-1
c              diff_new(i,Temp) = diff_new(i,RhoH) / 
c     $             (cp(i)*scal_new(i,Density))
           enddo

c*****************************************************************
c       Do chemistry update with const source = adv + diff_new
c*****************************************************************

           print *,'... advancing chem with const and linear sources'
           do n = 1,nscal
              do i = 0,nx-1
                 
                 const_src(i,n) = aofs(i,n)
     $                + diff_new(i,n) - diff_hat(i,n)
                 lin_src_old(i,n) = diff_old(i,n)
                 lin_src_new(i,n) = diff_hat(i,n)
                 
              enddo
           enddo
           
           call strang_chem(nx,scal_old,scal_new,
     $                      const_src,lin_src_old,lin_src_new,
     $                      intra,dt)

           sumh = 0.d0
           do i=0,nx-1
              sumh = sumh+abs(scal_new(i,RhoH)-scal_tmp(i,RhoH))
           enddo
           print *,'SUMH ',misdc,sumh
           
        enddo

c*****************************************************************
c       Now end iterations and go on to update velocity.
c*****************************************************************
        
        do i = 0,nx-1
           do n = 1,nspec
              ispec = FirstSpec-1+n
              Ydot_new(i,n) = intra(i,ispec)
           enddo
        enddo
        
c     HACKING THIS OUT
c        call calc_diffusivities(nx,scal_new,beta_new,dx,time+dt)
        call calc_divu(nx,scal_new,beta_new,Ydot_new,divu_new,dx,time)

        do i = 0,nx-1
           rhohalf(i) = 0.5d0*(scal_old(i,Density)+scal_new(i,Density))
        enddo
        
        Ydot_max = 0.d0
        do i = 0,nx-1
           dsdt(i) = (divu_new(i) - divu_old(i)) / dt
           Ydot_max = 
     $          max(abs(Ydot_new(i,1))/scal_new(i,Density),Ydot_max)
        enddo
        print *,'YDOT MAX ',Ydot_max
        
        print *,'... update velocities'
        call vel_edge_states(nx,vel_old,scal_old,gp,
     $       macvel,veledge,dx,dt)
        call update_vel(nx,vel_old,vel_new,gp,rhohalf,
     $       macvel,veledge,dx,dt)

c     FIXME: Add viscous terms
        
        print *,'...nodal projection...'
        call project(nx,vel_old,vel_new,rhohalf,divu_new,
     $       press_old,press_new,dx,dt)
        
        sumh = 0.d0
        do i=0,nx-1
           sumh = sumh+abs(scal_new(i,RhoH)-scal_old(i,RhoH))
        enddo
        print *,'SUMON ',sumh
        end

