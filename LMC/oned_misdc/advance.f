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
      real*8   veledge(0 :nx  )
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
      real*8   vel_Rhs(0:nx-1)
      real*8 rhocp_old
      real*8   pthermo(-1:nx  )
      real*8    Ydot_max, Y(maxspec)
      real*8 RWRK, cpmix, sum
      integer IWRK, is, rho_flag
      
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

      do n = 1,Nspec
         is = FirstSpec + n - 1
         do i = 0,nx-1
            intra(i,is) = Ydot_old(i,n)
         enddo
      enddo
c
c*****************************************************************
c     
      print *,'... creating the diffusive terms with old data'
      
      call calc_diffusivities(nx,scal_old,beta_old,mu_old,dx,time)
      call get_temp_visc_terms(nx,scal_old,beta_old,diff_old(0,Temp),dx)
      call get_spec_visc_terms(nx,scal_old,beta_old,
     &                         diff_old(0,FirstSpec),dx,time)
      call get_rhoh_visc_terms(nx,scal_old,beta_old,diff_old(0,RhoH),
     &                         dx,time)
      
c*****************************************************************
      
      print *,'... computing aofs with source = diff_old + intra'
      do i = 0,nx-1
c     Temp visc terms must be scaled by 1/(rho.cp) 
         do n = 1,Nspec
            Y(n) = scal_old(i,FirstSpec+n-1) / scal_old(i,Density)
         enddo
         call CKCPBS(scal_old(i,Temp),Y,IWRK,RWRK,cpmix)
         rhocp_old = cpmix * scal_old(i,Density)
         diff_old(i,n) = diff_old(i,n)/rhocp_old
         do n = 1,nscal
            tforce(i,n) = diff_old(i,n) + intra(i,n)
         enddo
      enddo
      
      call scal_aofs(nx,scal_old,macvel,aofs,tforce,dx,dt,time)

c*****************************************************************

      print *,'... update rho'
      call update_rho(nx,scal_old,scal_new,aofs,dx,dt)

c*****************************************************************
c     Either do c-n solve for new T prior to computing new 
c     coeffs, or simply start by copying from previous time step

      if (predict_temp_for_coeffs .eq. 1) then
         print *,'... update to temp w/old coeffs'
         be_cn_theta = 0.5d0
         rho_flag = 1
         call update_temp(nx,scal_old,scal_new,aofs,
     $                    alpha,beta_old,beta_new,intra(0,Temp),
     $                    Rhs(0,Temp),dx,dt,be_cn_theta,time)
         call cn_solve(nx,scal_new,alpha,beta_old,Rhs(0,Temp),
     $                 dx,dt,Temp,be_cn_theta,rho_flag)
         call get_hmix_given_T_RhoY(nx,scal_new,dx)      
         call calc_diffusivities(nx,scal_new,beta_new,mu_new,dx,time+dt)
      else
         print *,'... set new coeffs to old values for predictor'
         do n=1,nscal
            do i=-1,nx
               scal_new(i,Temp) = scal_old(i,Temp)
               beta_new(i,n) = beta_old(i,n)
            enddo
         enddo
      endif

c*****************************************************************

      print *,'... do predictor for species'
      be_cn_theta = 0.5d0
      call update_spec(nx,scal_old,scal_new,aofs,alpha,beta_old,
     &                 Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)
      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta,rho_flag)
      enddo

c     FIXME: Adjust spec flux at np1, reset scal_new, compute Le!=1 terms
      print *,'... do predictor for rhoh'
      call update_rhoh(nx,scal_old,scal_new,aofs,alpha,beta_old,
     &                 Rhs(0,RhoH),dx,dt,be_cn_theta,time)
      rho_flag = 2
      call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta,rho_flag)
      call rhoh_to_temp(nx,scal_new)

      print *,'... recompute diffusivities'
      call calc_diffusivities(nx,scal_new,beta_new,mu_new,dx,time+dt)

      print *,'... do corrector for species'
      call update_spec(nx,scal_old,scal_new,aofs,alpha,beta_old,
     $                 Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)
      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(nx,scal_new,alpha,beta_old,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta,rho_flag)
      enddo

c     FIXME: Adjust spec flux at np1, reset scal_new, compute Le!=1 terms      
      print *,'... do corrector for species'
      call update_rhoh(nx,scal_old,scal_new,aofs,alpha,beta_old,
     &                 Rhs(0,RhoH),dx,dt,be_cn_theta,time)
      rho_flag = 2
      call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta,rho_flag)
      call rhoh_to_temp(nx,scal_new)

      print *,'... create source terms for chem'
      call get_spec_visc_terms(nx,scal_new,beta_new,
     &                         diff_new(0,FirstSpec),dx,time)
      call get_rhoh_visc_terms(nx,scal_new,beta_new,diff_new(0,RhoH),
     &                         dx,time)

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

c*****************************************************************
           
           do n = 1,nscal
              do i = 0,nx-1
                 tforce(i,n) = intra(i,n) + diff_old(i,n)
              enddo
           enddo
           
           print *,'... computing aofs with source = diff_old + intra'
           call scal_aofs(nx,scal_old,macvel,aofs,tforce,dx,dt,time)
           
           call calc_diffusivities(nx,scal_new,beta_new,mu_new,
     &          dx,time+dt)

           print *,'... create new diff. terms : diff_hat'
           be_cn_theta = 0.5d0
           call get_rhoh_visc_terms(nx,scal_new,beta_new,
     &                              diff_hat(0,RhoH),dx,time)
           call get_temp_visc_terms(nx,scal_new,beta_new,
     &                              diff_hat(0,Temp),dx)
           call get_spec_visc_terms(nx,scal_new,beta_new,
     &                              diff_hat(0,FirstSpec),dx,time)
           
c*****************************************************************
           
           do n = 1,nscal
              do i = 0,nx-1
                 tforce(i,n) = intra(i,n)
     $                + 0.5d0 * (diff_old(i,n) + diff_hat(i,n))
     $                -          diff_hat(i,n)
              enddo
           enddo
           
           print *,'Updating A+D for species with intra terms'
           call update_spec(nx,scal_old,scal_new,aofs,alpha,beta_old,
     &                      Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)
           rho_flag = 2
           do n=1,Nspec
              is = FirstSpec + n - 1
              call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,is),
     $                      dx,dt,is,be_cn_theta,rho_flag)
           enddo

           print *,'... update A+D for RhoH with intra terms'
           call update_rhoh(nx,scal_old,scal_new,aofs,alpha,beta_old,
     &                      Rhs(0,RhoH),dx,dt,be_cn_theta,time)
           rho_flag = 2
           call cn_solve(nx,scal_new,alpha,beta_new,Rhs(0,RhoH),
     $                   dx,dt,RhoH,be_cn_theta,rho_flag)
           print *,'... create new temp from new rhoH'
           call rhoh_to_temp(nx,scal_new)

           print *,'... create new diff. terms for RhoH: diff_new'
           if (be_cn_theta .eq. 1.d0) then
              do i = 0,nx-1
                 diff_new(i,RhoH) = 
     $                (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $                -aofs(i,RhoH)-tforce(i,RhoH)
              enddo
           else if (be_cn_theta .eq. 0.5d0) then
              do i = 0,nx-1
                 diff_new(i,RhoH) = 2.d0 * (
     $                (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $                -aofs(i,RhoH)-tforce(i,RhoH)
     $                -0.5d0*diff_old(i,RhoH))
              enddo
           else 
              print *,'OOPS - BOGUS BE_CN_THETA ',be_cn_theta
              stop
           endif
           
           
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

        print *,'adv: hacking state'
        do i=-1,nx
           do n=1,nscal
              scal_new(i,n) = scal_old(i,n)
           enddo
        enddo
        do i=0,nx-1
           do n=1,Nspec
              Ydot_new(i,n) = Ydot_old(i,n)
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

        call vel_edge_states(nx,vel_old,scal_old(-1,Density),gp,
     $       macvel,veledge,dx,dt,time)

        call update_vel(nx,vel_old,vel_new,gp,rhohalf,
     &                  macvel,veledge,alpha,mu_old,
     &                  vel_Rhs,dx,dt,be_cn_theta,time)
        rho_flag = 1
        call cn_solve(nx,vel_new,alpha,mu_new,vel_Rhs,
     $                dx,dt,1,be_cn_theta,rho_flag)

        print *,'...nodal projection...'
        call project(nx,vel_old,vel_new,rhohalf,divu_new,
     $       press_old,press_new,dx,dt)
        
        sumh = 0.d0
        do i=0,nx-1
           sumh = sumh+abs(scal_new(i,RhoH)-scal_old(i,RhoH))
        enddo
        print *,'SUMON ',sumh
        end

