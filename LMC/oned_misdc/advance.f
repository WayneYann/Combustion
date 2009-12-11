      subroutine advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R_old,I_R_new,press_old,press_new,
     $                   divu_old,divu_new,dsdt,beta_old,beta_new,
     $                   dx,dt,time)

      implicit none
      include 'spec.h'
      real*8   vel_new(-1:nx  )
      real*8   vel_old(-1:nx  )
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
      real*8 press_new(0 :nx  )
      real*8 press_old(0 :nx  )
      real*8   I_R_new(0:nx-1,0:maxspec)
      real*8   I_R_old(0:nx-1,0:maxspec)
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
      real*8        cp(0 :nx-1)
      real*8 x
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta
      
      real*8    diff_old(0:nx-1,nscal)
      real*8    diff_hat(0:nx-1,nscal)
      real*8    diff_new(0:nx-1,nscal)
      real*8   const_src(0:nx-1,nscal)
      real*8 lin_src_old(0:nx-1,nscal)
      real*8 lin_src_new(0:nx-1,nscal)
      real*8    scal_tmp(0:nx-1,nscal)
      integer misdc
      
      integer i,n,ispec
      integer iunit
      
      real*8 divu_max
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8      dRhs(0:nx-1,0:maxspec)
      real*8   vel_Rhs(0:nx-1)
      real*8 rhocp_old
      real*8   pthermo(-1:nx  )
      real*8    Ydot_max, Y(maxspec)
      real*8 RWRK, cpmix, sum
      integer IWRK, is, rho_flag
      
      print *,'advance: at start of time step'
      be_cn_theta = 0.5d0
c     
c*****************************************************************
c     Create MAC velocities.
c*****************************************************************
c           
C CEG:: just checking
C      write(*,*)'dt = ', dt
C      call write_plt(vel_old,scal_old,press_old,dx,0,time)
C      call write_plt(vel_new,scal_new,press_new,dx,1,time)
C      stop

      do i = 0,nx-1
         gp(i) = (press_old(i+1) - press_old(i)) / dx
      enddo
      
c      call minmax_scal(nx,scal_old)
      print *,'... predict edge velocities'
      call pre_mac_predict(vel_old,scal_old,gp,
     $                     macvel,dx,dt,time)
      
      call compute_pthermo(scal_old,pthermo)

      do i = 0,nx-1
         divu_tmp(i) = divu_old(i) + 0.5d0 * dt * dsdt(i)
      enddo
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU norm old = ',divu_max 
      call add_dpdt(pthermo,divu_tmp,macvel,dx,dt)
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU norm new = ',divu_max 
      call macproj(nx,macvel,divu_tmp,dx)

c
c*****************************************************************
c     
      print *,'... creating the diffusive terms with old data'
      
      call calc_diffusivities(scal_old,beta_old,mu_old,dx,time)
      call get_temp_visc_terms(scal_old,beta_old,
     &                         diff_old(0,Temp),dx,time)
      call get_spec_visc_terms(scal_old,beta_old,
     &                         diff_old(0,FirstSpec),dx,time)
      call get_rhoh_visc_terms(scal_old,beta_old,
     &                         diff_old(0,RhoH),dx,time)
      
c*****************************************************************
      
      print *,'... computing aofs with D(old) + R(guess)'
      do i = 0,nx-1
         do n = 0,Nspec
            I_R_new(i,n) = I_R_old(i,n)
         enddo
      enddo

      do i = 0,nx-1
         do n = 1,Nspec
            Y(n) = scal_old(i,FirstSpec+n-1) / scal_old(i,Density)
         enddo
         call CKCPBS(scal_old(i,Temp),Y,IWRK,RWRK,cpmix)
         rhocp_old = cpmix * scal_old(i,Density)
         diff_old(i,Temp) = diff_old(i,Temp)/rhocp_old

         do n = 1,Nspec
            is = FirstSpec + n - 1
            tforce(i,is) = diff_old(i,is) + I_R_new(i,n)
         enddo
         tforce(i,Temp) = diff_old(i,Temp) + I_R_new(i,0)
      enddo
       
      call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)

c*****************************************************************

      print *,'... update rho'
      call update_rho(scal_old,scal_new,aofs,dx,dt)

c*****************************************************************
c     Either do c-n solve for new T prior to computing new 
c     coeffs, or simply start by copying from previous time step

      if (predict_temp_for_coeffs .eq. 1) then
         print *,'... predict temp with old coeffs'
         rho_flag = 1
         call update_temp(scal_old,scal_new,aofs,
     $                    alpha,beta_old,beta_new,I_R_new(0,0),
     $                    Rhs(0,Temp),dx,dt,be_cn_theta,time)
C CEG:: just uses RHS and overwrites snew
         call cn_solve(scal_new,alpha,beta_old,Rhs(0,Temp),
     $                 dx,dt,Temp,be_cn_theta,rho_flag)
         call get_hmix_given_T_RhoY(scal_new,dx)      

         print *,'... compute new coeffs'
         call calc_diffusivities(scal_new,beta_new,mu_new,dx,time+dt)
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
      print *,'... do predictor for species (MISDC terms=0)'
      do i=0,nx-1
         dRhs(i,0) = 0.0d0
         do n=1,Nspec
            dRhs(i,n) = dt*I_R_new(i,n)
         enddo
      enddo
      call update_spec(scal_old,scal_new,aofs,alpha,beta_old,
     &     dRhs(0,1),Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)
      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta,rho_flag)
      enddo

      print *,'... do predictor for rhoh (MISDC terms=0)'
      call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,
     &     dRhs(0,0),Rhs(0,RhoH),dx,dt,be_cn_theta,time)
      rho_flag = 2
      call cn_solve(scal_new,alpha,beta_new,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta,rho_flag)
      call rhoh_to_temp(scal_new)

      print *,'...   extract D sources'
      do i = 0,nx-1
         diff_new(i,RhoH) = (
     $        (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $        - aofs(i,RhoH) - 
     $        (1.d0-be_cn_theta)*diff_old(i,RhoH) )/be_cn_theta
         do n=1,Nspec
            is = FirstSpec + n - 1
            diff_new(i,is) = (
     $           (scal_new(i,is)-scal_old(i,is))/dt 
     $           - aofs(i,is) - I_R_new(i,n) - 
     $           (1.d0-be_cn_theta)*diff_old(i,is) )/be_cn_theta
         enddo
      enddo

      print *,'... advance chem with A+D sources, reset I_R_new'
      do n = 1,nscal
         do i = 0,nx-1
            const_src(i,n) =     aofs(i,n)
            lin_src_old(i,n) = diff_old(i,n)
            lin_src_new(i,n) = diff_new(i,n)
         enddo
      enddo
      
C      write(20,*)'nscal = ',nscal
C      write(20,*)Density,Temp,RhoH,FirstSpec,RhoRT
C      do i = 0,nx-1
C         do n = 1,nscal
C            write(20,*)i,aofs(i,n),diff_old(i,n),diff_new(i,n)
C         enddo
C      enddo

      call strang_chem(scal_old,scal_new,
     $                 const_src,lin_src_old,lin_src_new,
     $                 I_R_new,dt)

C----------------------------------------------------------------
C----------------------------------------------------------------
      do misdc = 1, misdc_iterMAX

         be_cn_theta = 1.d0

         print *,'MISDC iteration ',misdc

         print *,'... create new diff_hat from current state'
         call calc_diffusivities(scal_new,beta_new,mu_new,
     &                           dx,time+dt)
         call get_temp_visc_terms(scal_new,beta_new,
     &                            diff_hat(0,Temp),dx,time+dt)
         call get_spec_visc_terms(scal_new,beta_new,
     &                            diff_hat(0,FirstSpec),dx,time+dt)

         do i = 0,nx-1
            do n = 1,Nspec
               ispec = FirstSpec + n - 1
               tforce(i,ispec) = I_R_new(i,n)
     &              + 0.5d0*(diff_old(i,ispec)+diff_hat(i,ispec))
            enddo
            tforce(i,Temp) = I_R_new(i,0)
     &              + 0.5d0*(diff_old(i,Temp)+diff_hat(i,Temp))
         enddo
         
         print *,'... compute A with updated D+R source'
         call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)


c*****************************************************************

         print *,'... update D for species with A + R + MISDC(D)'
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
C could take out diff_old if i left theta = 0.5
C will need to think about this for the case of Le != 1
               dRhs(i,n) = dt*(I_R_new(i,n) 
     &              + 0.5d0*(diff_old(i,is) - diff_hat(i,is)))
            enddo
            dRhs(i,0) = dt*(
     &           + 0.5d0*(diff_old(i,RhoH) - diff_new(i,RhoH)))
         enddo
         call update_spec(scal_old,scal_new,aofs,alpha,beta_old,
     &        dRhs(0,1),Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)
         rho_flag = 2
         do n=1,Nspec
            is = FirstSpec + n - 1
            call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $                    dx,dt,is,be_cn_theta,rho_flag)
         enddo

         print *,'... update D for rhoh with A + R + MISDC(D)'
         call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,
     &        dRhs(0,0),Rhs(0,RhoH),dx,dt,be_cn_theta,time)
         rho_flag = 2
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,RhoH),
     $                 dx,dt,RhoH,be_cn_theta,rho_flag)
         print *,'... create new temp from new RhoH, spec'
         call rhoh_to_temp(scal_new)

         print *,'... create diff_new from RhoH and spec solutions'
         if (be_cn_theta .ne. 0.d0) then
            do i = 0,nx-1
               diff_new(i,RhoH) = (
     $              (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $              - aofs(i,RhoH) - dRhs(i,0)/dt - 
     $              (1.d0-be_cn_theta)*diff_old(i,RhoH) )/be_cn_theta
               do n=1,Nspec
                  is = FirstSpec + n - 1
                  diff_new(i,is) = (
     $                 (scal_new(i,is)-scal_old(i,is))/dt 
     $                 - aofs(i,is) - dRhs(i,n)/dt - 
     $                 (1.d0-be_cn_theta)*diff_old(i,is) )/be_cn_theta
               enddo
            enddo
         endif
 
         print *,'... advancing chem with const and linear sources'
         do n = 1,nscal
            do i = 0,nx-1
               
               const_src(i,n) = aofs(i,n)
     $              + diff_new(i,n) - diff_hat(i,n)
               lin_src_old(i,n) = diff_old(i,n)
               lin_src_new(i,n) = diff_hat(i,n)
            enddo
         enddo
          
         call strang_chem(scal_old,scal_new,
     $                    const_src,lin_src_old,lin_src_new,
     $                    I_R_new,dt)


c*****************************************************************
c       End of MISDC iterations
c*****************************************************************
      enddo
      
      
c     call calc_diffusivities(scal_new,beta_new,dx,time+dt)
      call calc_divu(scal_new,beta_new,I_R_new,divu_new,dx,time+dt)

      do i = 0,nx-1
         rhohalf(i) = 0.5d0*(scal_old(i,Density)+scal_new(i,Density))
         dsdt(i) = (divu_new(i) - divu_old(i)) / dt
      enddo         

      print *,'... update velocities'
      
      call vel_edge_states(vel_old,scal_old(-1,Density),gp,
     $                     macvel,veledge,dx,dt,time)
      
      call update_vel(vel_old,vel_new,gp,rhohalf,
     &                macvel,veledge,alpha,mu_old,
     &                vel_Rhs,dx,dt,be_cn_theta,time)
      rho_flag = 1
      call cn_solve(vel_new,alpha,mu_new,vel_Rhs,
     $              dx,dt,1,be_cn_theta,rho_flag)

      print *,'...nodal projection...'
      call project(vel_old,vel_new,rhohalf,divu_new,
     $             press_old,press_new,dx,dt)
        
      end


