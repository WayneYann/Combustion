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
      real*8 theta
      real*8 vel_theta
      
      real*8    diff_old(0:nx-1,nscal)
      real*8    diff_hat(0:nx-1,nscal)
      real*8    diff_new(0:nx-1,nscal)
      real*8   const_src(0:nx-1,nscal)
      real*8 lin_src_old(0:nx-1,nscal)
      real*8 lin_src_new(0:nx-1,nscal)
      integer misdc
      
      integer i,n,ispec
      integer iunit
      
      real*8 divu_max
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8      dRhs(0:nx-1,0:maxspec)
      real*8   vel_Rhs(0:nx-1)
      real*8 rhocp_old,  rhocp
      real*8 Tmid
      real*8   pthermo(-1:nx  )
      real*8    Ydot_max, Y(maxspec)
      real*8 RWRK, cpmix, sum
      integer IWRK, is, rho_flag
      
C CEG debugging FIXME
      real*8 ptherm(-1:nx)
      real*8 dummy
      integer j, imax, imin
      real*8  Schange(-1:nx  ,nscal)
      real*8  change_max(nscal)
      real*8  change_min(nscal)


      print *,'advance: at start of time step'
      be_cn_theta = 1.d0
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
C CEG:: this fills ghost cells for vel_old
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
      call add_dpdt(scal_old,pthermo,divu_tmp,macvel,dx,dt)
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU norm new = ',divu_max 
      call macproj(nx,macvel,divu_tmp,dx)

C CEG FIXME
C$$$            do n=0,Nspec
C$$$               do i=0,nx-1
C$$$                  I_R_new(i,n) = 0.d0
C$$$               enddo
C$$$            enddo
CCCCCCCCC

c
c*****************************************************************
c     
C CEG:: make sure to get diffusivities at time n (old time)
      call calc_diffusivities(scal_old,beta_old,mu_old,dx,time)

      if (use_strang) then
         call strang_advance(macvel,scal_old,scal_new,
     $                   I_R_old,I_R_new,beta_old,beta_new,
     $                   dx,dt,time)
      else
         if (use_radau) then
            if (use_temp_eqn) then
               call advance_radau_temp(macvel,scal_old,scal_new,
     $              I_R_old,I_R_new,beta_old,beta_new,
     $              dx,dt,time)               
            else
               call advance_radau(macvel,scal_old,scal_new,
     $              I_R_old,I_R_new,beta_old,beta_new,
     $              dx,dt,time)
            endif
         else
         if (use_temp_eqn) then
            call advance_temp(macvel,scal_old,scal_new,
     $                   I_R_old,I_R_new,beta_old,beta_new,
     $                   dx,dt,time)
         else
            print *,'... using LOBATTO quadtrature'
            print *,'... evolving WITHOUT using temp eqn'
            print *,'... creating the diffusive terms with old data'
            print *,'... set new coeffs to old values for predictor'
            
C for provisional can use limited slopes
            unlim = 0
            lim_rxns = 0
C CEG:: this should already be true, but jus tbeing pedantic
            do i=-1,nx
               scal_new(i,Temp) = scal_old(i,Temp)
            enddo

C     CEG:: each one of these functions first calls set_bc(scal_old)
C     maybe should change this
            call get_spec_visc_terms(scal_old,beta_old,
     &           diff_old(0,FirstSpec),dx,time)
            call divRhoDHgradY(scal_old,beta_old,diff_old(0,RhoH),
     &           dx,time)
            call addDivLambdaGradT(scal_old,beta_old,
     &           diff_old(0,RhoH),dx,time)
c*****************************************************************
            
            print *,'... computing aofs with D(old) + R(guess)'

            do i = 0,nx-1
               do n = 1,Nspec
                  is = FirstSpec + n - 1
                  tforce(i,is) = diff_old(i,is) + I_R_new(i,n)
               enddo
               tforce(i,RhoH) = diff_old(i,RhoH)
            enddo
CCCCCCCCCCCCCCCCCCCCCc
 1006 FORMAT(11(E22.15,1X)) 
         call compute_pthermo(scal_old,ptherm)
         open(UNIT=11, FILE='sold.dat', STATUS = 'REPLACE')
         write(11,*)'# 256 12'
         do j=0,nx-1
            do n = 1,Nspec
               Y(n) = scal_old(j,FirstSpec+n-1)/scal_old(j,Density)
            enddo
            write(11,1006) (j+.5)*dx, 
     &                     (tforce(j,FirstSpec+n-1),n=1,Nspec),
     $                     tforce(j,RhoH)
         enddo
         close(11)
CCCCCCCCCCCCCCCCCCCCCCCC            
            call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)

c*****************************************************************

            print *,'... update rho'
            call update_rho(scal_old,scal_new,aofs,dx,dt)


c*****************************************************************
            print *,'... do predictor for species (MISDC terms=0)'
C CEG:: doesn't seem to make much difference one way or the other
            theta = 1.d0
            do i=0,nx-1
               dRhs(i,0) = dt*(1.0d0 - theta)*diff_old(i,RhoH)
               do n=1,Nspec
                  dRhs(i,n) = dt*I_R_new(i,n)
C0.d0
               enddo
            enddo

C CEG:: Trying something different FIXME??? 
            print *,'... do predictor for rhoh (MISDC terms=0)'
            call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,
     &           dRhs(0,0),Rhs(0,Temp),dx,dt,theta,time)
C Implicit solve for Temp^n+1
            rho_flag = 1
            call cn_solve(scal_new,alpha,beta_old,Rhs(0,Temp),
     $           dx,dt,Temp,theta,rho_flag)
            call calc_diffusivities(scal_new,beta_new,mu_new,dx,time+dt)
CCCCCCCCCCC

            call update_spec(scal_old,scal_new,aofs,alpha,beta_old,
     &           dRhs(0,1),Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)

            rho_flag = 2
            do n=1,Nspec
               is = FirstSpec + n - 1
               call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $              dx,dt,is,be_cn_theta,rho_flag)
            enddo

C get a better estimate for rhoH(T)
            print *,'... do predictor for rhoh (MISDC terms=0)'
            call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,
     &           dRhs(0,0),Rhs(0,Temp),dx,dt,theta,time)
C Implicit solve for Temp^n+1
            rho_flag = 1
            call cn_solve(scal_new,alpha,beta_old,Rhs(0,Temp),
     $           dx,dt,Temp,theta,rho_flag)

            print *,'...   extract D sources'
C CEG;; note that neither of these 2 fns use rhoH_new
            call divRhoDHgradY(scal_new,beta_new,diff_new(0,RhoH),
     &           dx,time+dt)
            call addDivLambdaGradT(scal_new,beta_new,
     &           diff_new(0,RhoH),dx,time+dt)
C CEG; this vs below doesn't seem to make any difference
            call get_spec_visc_terms(scal_new,beta_new,
     &           diff_new(0,FirstSpec),dx,time+dt)

C$$$            if (be_cn_theta .ne. 0.d0) then
C$$$               do i = 0,nx-1
C$$$                  do n=1,Nspec
C$$$                     is = FirstSpec + n - 1
C$$$                     diff_new(i,is) = (
C$$$     $                    (scal_new(i,is)-scal_old(i,is))/dt 
C$$$     $                    - aofs(i,is) - I_R_new(i,n) - 
C$$$     $                    (1.d0-be_cn_theta)*diff_old(i,is) 
C$$$     $                    )/be_cn_theta
C$$$                  enddo
C$$$               enddo
C$$$            else
C$$$               print *,'ERROR:: be_cn_theta=0.0d0 case not coded yet'
C$$$               stop
C$$$            endif

            if (nochem_hack) then
               print *,'WARNING! doing nochem_hack--skipping reactions'
            else
               print *,'... react with A+D sources, reset I_R_new'
               do n = 1,nscal
                  do i = 0,nx-1
C                     const_src(i,n) =     aofs(i,n)
C                     lin_src_old(i,n) = diff_old(i,n)
C                     lin_src_new(i,n) = diff_new(i,n)
                     const_src(i,n) =  aofs(i,n) + diff_new(i,n)
                     lin_src_old(i,n) = 0.d0
                     lin_src_new(i,n) = 0.d0
C treating rhoh differently than the species doesn't change anything
C                     const_src(i,RhoH) =   aofs(i,RhoH)+diff_new(i,RhoH)
C                     lin_src_old(i,RhoH) = 0.d0
C                     lin_src_new(i,RhoH) = 0.d0
                  enddo
               enddo

               call strang_chem(scal_old,scal_new,
     $              const_src,lin_src_old,lin_src_new,
     $              I_R_new,dt)
            endif
CCCCCCCCCCCCCCCCCCCCCc
         open(UNIT=11, FILE='source.dat', STATUS = 'REPLACE')
         write(11,*)'# 256 12'
         do j=0,nx-1
            do n = 1,Nspec
               Y(n) = scal_old(j,FirstSpec+n-1)/scal_old(j,Density)
            enddo
            write(11,1006) (j+.5)*dx, 
     &                     (const_src(j,FirstSpec+n-1),n=1,Nspec),
     $                     const_src(j,RhoH)
         enddo
         close(11)
         open(UNIT=11, FILE='aofs.dat', STATUS = 'REPLACE')
         write(11,*)'# 256 12'
         do j=0,nx-1
            do n = 1,Nspec
               Y(n) = scal_old(j,FirstSpec+n-1)/scal_old(j,Density)
            enddo
            write(11,1006) (j+.5)*dx, 
     &                     (aofs(j,FirstSpec+n-1),n=1,Nspec),
     $                     aofs(j,RhoH)
         enddo
         close(11)
         open(UNIT=11, FILE='diff.dat', STATUS = 'REPLACE')
         write(11,*)'# 256 12'
         do j=0,nx-1
            do n = 1,Nspec
               Y(n) = scal_old(j,FirstSpec+n-1)/scal_old(j,Density)
            enddo
            write(11,1006) (j+.5)*dx, 
     &                     (diff_new(j,FirstSpec+n-1),n=1,Nspec),
     $                     diff_new(j,RhoH)
         enddo
         close(11)
CCCCCCCCCCCCCCCCCCCCCCCC            

C     CEG debugging FIXME
C     
C     Find the estimated change in S over the timestep
C     
            do n = 1,nscal
               change_max(n) = 0.d0
               change_min(n) = 0.d0
            enddo
            do i = 0,nx-1
               do n = 1,nscal
                  Schange(i,n) = scal_new(i,n) - scal_old(i,n)
                  change_max(n) = MAX(change_max(n),Schange(i,n))
                  change_min(n) = MIN(change_min(n),Schange(i,n))
               enddo 
            enddo
            write(*,*)
            write(*,*)'Change in S over the timestep'
            write(*,*)'index      min      max'
            do n = 1,nscal
               write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
            enddo
            do i = 0,nx-1
               do n = 1,nscal
                  Schange(i,n) = scal_new(i,n)
               enddo
            enddo
 1008       FORMAT((I5,1X),(E22.15,1X))      
C----------------------------------------------------------------
C----------------------------------------------------------------
            do misdc = 1, misdc_iterMAX
C use unlimited slopes here
               unlim = 1
               lim_rxns = 1

               print *,'... doing SDC iter ',misdc

               print *,'... create new diff_hat from current state'
               call calc_diffusivities(scal_new,beta_new,mu_new,
     &              dx,time+dt)

               call get_spec_visc_terms(scal_new,beta_new,
     &              diff_hat(0,FirstSpec),dx,time+dt)
               call divRhoDHgradY(scal_new,beta_new,diff_hat(0,RhoH),
     &              dx,time+dt)
               call addDivLambdaGradT(scal_new,beta_new,
     &              diff_hat(0,RhoH),dx,time+dt)

               do i = 0,nx-1
                  do n = 1,Nspec
                     ispec = FirstSpec + n - 1
                     tforce(i,ispec) = I_R_new(i,n)
     &                    + 0.5d0*(diff_old(i,ispec)+diff_hat(i,ispec))
                  enddo
                  tforce(i,RhoH) =
     &                 + 0.5d0*(diff_old(i,RhoH)+diff_hat(i,RhoH))
               enddo
               
               print *,'... compute A with updated D+R source'
               call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)

               print *,'... update rho'
               call update_rho(scal_old,scal_new,aofs,dx,dt)

c*****************************************************************

               print *,'... update D for species with A + R + MISDC(D)'
               do i=0,nx-1
                  do n=1,Nspec
                     is = FirstSpec + n - 1
                     dRhs(i,n) = dt*(I_R_new(i,n) 
     &                    + 0.5d0*(diff_old(i,is) - diff_hat(i,is)))
                  enddo
                  dRhs(i,0) = dt * 0.5d0 
     $                 * (diff_old(i,RhoH)-diff_hat(i,RhoH))
               enddo

               call update_spec(scal_old,scal_new,aofs,alpha,beta_old,
     &              dRhs(0,1),Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)
               rho_flag = 2
               do n=1,Nspec
                  is = FirstSpec + n - 1
                  call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta,rho_flag)
               enddo

               print *,'... update D for rhoh with A + R + MISDC(D)'
               call update_rhoh(scal_old,scal_new,aofs,alpha,beta_new,
     &              dRhs(0,0),Rhs(0,Temp),dx,dt,be_cn_theta,time)
C     Implicit solve for Temp^n+1
               rho_flag = 1
               call cn_solve(scal_new,alpha,beta_new,Rhs(0,Temp),
     $              dx,dt,Temp,be_cn_theta,rho_flag)

               print *,'...   extract D sources'
C     CEG;; note that neither of these 2 fns use rhoH_new
               call divRhoDHgradY(scal_new,beta_new,diff_new(0,RhoH),
     &              dx,time+dt)
               call addDivLambdaGradT(scal_new,beta_new,
     &              diff_new(0,RhoH),dx,time+dt)

               print *,'... create diff_new from RhoH & spec solutions'
               if (be_cn_theta .ne. 0.d0) then
                  do i = 0,nx-1
                     do n=1,Nspec
                        is = FirstSpec + n - 1
                        diff_new(i,is) = (
     $                       (scal_new(i,is)-scal_old(i,is))/dt 
     $                       - aofs(i,is) - dRhs(i,n)/dt - 
     $                       (1.d0-be_cn_theta)*diff_old(i,is) 
     $                       )/be_cn_theta
                     enddo
                  enddo
               else
                  print *,'ERROR:: be_cn_theta=0.0d0 case not coded yet'
                  stop
               endif
               
               if (nochem_hack) then
                  print *,'WARNING: SDC nochem_hack--skipping reactions'
               else
                  print *,'... react with const and linear sources'
                  do n = 1,nscal
                     do i = 0,nx-1
                        const_src(i,n) = aofs(i,n)
     $                       + diff_new(i,n) - diff_hat(i,n)
                        lin_src_old(i,n) = diff_old(i,n)
                        lin_src_new(i,n) = diff_hat(i,n)
                     enddo                     
                  enddo
                  
                  call strang_chem(scal_old,scal_new,
     $                 const_src,lin_src_old,lin_src_new,
     $                 I_R_new,dt)
               endif

c*****************************************************************
c     End of MISDC iterations
c*****************************************************************
C     CEG debugging FIXME
C     
C     Find the size of the correction, ie the change in S_new
C     
               do n = 1,nscal
                  change_max(n) = 0.d0
                  change_min(n) = 0.d0
               enddo
               do i = 0,nx-1
                  do n = 1,nscal
                     Schange(i,n) = scal_new(i,n) - Schange(i,n)
                     change_max(n) = MAX(change_max(n),Schange(i,n))
                     change_min(n) = MIN(change_min(n),Schange(i,n))
                  enddo 
               enddo
               write(*,*)
               write(*,*)'Size of the correction (Change in S_new)'
               write(*,*)'index      min      max'
               do n = 1,nscal
C     write(*,*)n,change_min(n),change_max(n)
                  write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
               enddo
               do i = 0,nx-1
                  do n = 1,nscal
                     Schange(i,n) = scal_new(i,n)
                  enddo
               enddo

            enddo

         endif
C end if(use temp eqn)
         endif
         endif
C end strang vs SDC
            
      call calc_diffusivities(scal_new,beta_new,mu_new,dx,time+dt)
      call calc_divu(scal_new,beta_new,I_R_new,divu_new,dx,time+dt)

      do i = 0,nx-1
         rhohalf(i) = 0.5d0*(scal_old(i,Density)+scal_new(i,Density))
         dsdt(i) = (divu_new(i) - divu_old(i)) / dt
      enddo         
C debugging FIXME
C$$$ 1007 FORMAT((I5,1X),(E22.15,1X))      
C$$$         open(UNIT=11, FILE='dsdt.dat', STATUS = 'REPLACE')
C$$$         write(11,*)'# 256 2'
C$$$         do i=0,nx-1
C$$$            write(11,1007) i, dsdt(i)
C$$$         enddo
C$$$         close(11)
C$$$         write(*,*)'divu update'
C$$$         stop
CCCCCCCCCCCCC

      print *,'... update velocities'

      vel_theta = 0.5d0
C get velocity visc terms to use as a forcing term for advection
      call get_vel_visc_terms(vel_old,mu_old,visc,dx,time)
      do i = 0, nx-1
         visc(i) = visc(i)/scal_old(i,Density)
      enddo

      call vel_edge_states(vel_old,scal_old(-1,Density),gp,
     $                     macvel,veledge,dx,dt,time,visc)
      
      call update_vel(vel_old,vel_new,gp,rhohalf,
     &                macvel,veledge,alpha,mu_old,
     &                vel_Rhs,dx,dt,vel_theta,time)

      if (initial_iter .eq. 1) then
         call get_vel_visc_terms(vel_old,mu_old,visc,dx,time)
         do i = 0, nx-1
            vel_new(i) = vel_new(i) + visc(i)*dt/rhohalf(i)
         enddo
      else
         rho_flag = 1
         call cn_solve(vel_new,alpha,mu_new,vel_Rhs,
     $        dx,dt,1,vel_theta,rho_flag)
      endif

      print *,'...nodal projection...'
C CEG:: LMC has another var initial_step, and doesn't do the proj
C       until after all num_init_iters are done
      if (initial_iter .eq. 0) then
         call project(vel_old,vel_new,rhohalf,divu_new,
     $        press_old,press_new,dx,dt)
      endif
      end

      subroutine advance_temp (macvel,scal_old,scal_new,
     $                   I_R_old,I_R_new,
     $                   beta_old,beta_new,dx,dt,time)

      implicit none
      include 'spec.h'
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
      real*8   I_R_new(0:nx-1,0:maxspec)
      real*8   I_R_old(0:nx-1,0:maxspec)
      real*8    macvel(0 :nx  )
      real*8      aofs(0 :nx-1,nscal)
      real*8  beta_old(-1:nx,nscal)
      real*8  beta_new(-1:nx,nscal)
      real*8  mu_dummy(-1:nx)
      real*8   rhohalf(0 :nx-1)
      real*8    tforce(0 :nx-1,nscal)
      real*8      visc(0 :nx-1)
      real*8        cp(0 :nx-1)
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta
      real*8 theta
      
      real*8    diff_old(0:nx-1,nscal)
      real*8    diff_new(0:nx-1,nscal)
      real*8    diff_hat(0:nx-1,nscal)
      real*8   const_src(0:nx-1,nscal)
      real*8 lin_src_old(0:nx-1,nscal)
      real*8 lin_src_new(0:nx-1,nscal)
      
      integer i,n,ispec
      integer iunit
      
      real*8 divu_max
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8      dRhs(0:nx-1,0:maxspec)
      real*8 rhocp_old, rhocp
      real*8 Tmid
      real*8   pthermo(-1:nx  )
      real*8    Ydot_max, Y(maxspec)
      real*8 RWRK, cpmix, sum
      integer IWRK, is, rho_flag
      integer misdc

C CEG debugging FIXME
      real*8 ptherm(-1:nx)
      integer j
      real*8  Schange(-1:nx  ,nscal)
      real*8  change_max(nscal)
      real*8  change_min(nscal)


      be_cn_theta = 1.0d0

      print *,'... using LOBATTO quaadrature'
      print *,'... evolving using temperature equation'
      print *,'... creating the diffusive terms with old data'

C CEG:: each one of these functions first calls set_bc(scal_old)
C   maybe should change this
      call get_temp_visc_terms(scal_old,beta_old,
     &                         diff_old(0,Temp),dx,time)
      call get_spec_visc_terms(scal_old,beta_old,
     &                         diff_old(0,FirstSpec),dx,time)
      call get_rhoh_visc_terms(scal_old,beta_old,
     &                         diff_old(0,RhoH),dx,time)
      
c*****************************************************************
      
      print *,'... computing aofs with D(old) + R(guess)'

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
         theta = 0.5d0
C CEG:: beta_new comes into advance() with the same value as 
C       beta_old (the one we just calculated)
         call update_temp(scal_old,scal_new,aofs,
     $                    alpha,beta_old,beta_new,I_R_new(0,0),
     $                    Rhs(0,Temp),dx,dt,theta,time)
C CEG:: just uses RHS and overwrites snew
C does not fill ghost cells
         call cn_solve(scal_new,alpha,beta_old,Rhs(0,Temp),
     $                 dx,dt,Temp,theta,rho_flag)

         call get_hmix_given_T_RhoY(scal_new,dx)      

         print *,'... compute new coeffs'
         call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time+dt)
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
      if (be_cn_theta .ne. 0.d0) then
         do i = 0,nx-1
            diff_new(i,RhoH) = (
     $           (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $           - aofs(i,RhoH) -
     $           (1.d0-be_cn_theta)*diff_old(i,RhoH) )/be_cn_theta
            do n=1,Nspec
               is = FirstSpec + n - 1
               diff_new(i,is) = (
     $              (scal_new(i,is)-scal_old(i,is))/dt 
     $              - aofs(i,is) - I_R_new(i,n) - 
     $              (1.d0-be_cn_theta)*diff_old(i,is) )/be_cn_theta
            enddo
         enddo
      else
         print *,'ERROR:: not set up to work with be_cn_theta=0.0d0'
         stop
      endif

      if (nochem_hack) then
         write(*,*)'WARNING! doing nochem_hack--skipping reactions'
      else
         print *,'... react with A+D sources, reset I_R_new'
         do n = 1,nscal
            do i = 0,nx-1
               const_src(i,n) =     aofs(i,n)
C CEG:: if use linear diffusion and don't do any sdc iterations then
C       the solution develops kinks.  Using only D^n+1 leads to a solution
C       that drifts far off the equation of state (dp/dt thing not
C       implemented yet)
C               lin_src_old(i,n) = diff_old(i,n)
C               lin_src_new(i,n) = diff_new(i,n)
               lin_src_old(i,n) = diff_new(i,n)
               lin_src_new(i,n) = 0.d0

            enddo
         enddo

         call strang_chem(scal_old,scal_new,
     $                    const_src,lin_src_old,lin_src_new,
     $                    I_R_new,dt)
      endif

C CEG debugging FIXME
C
C Find the estimated change in S over the timestep
C
      do n = 1,nscal
         change_max(n) = 0.d0
         change_min(n) = 0.d0
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_new(i,n) - scal_old(i,n)
            change_max(n) = MAX(change_max(n),Schange(i,n))
            change_min(n) = MIN(change_min(n),Schange(i,n))
         enddo 
      enddo
      write(*,*)
      write(*,*)'Change in S over the timestep'
      write(*,*)'index      min      max'
      do n = 1,nscal
         write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_new(i,n)
         enddo
      enddo
 1008 FORMAT((I5,1X),(E22.15,1X))      
C----------------------------------------------------------------
C----------------------------------------------------------------
      do misdc = 1, misdc_iterMAX
         print *,'... doing SDC iter ',misdc

         print *,'... create new diff_hat from current state'
         call calc_diffusivities(scal_new,beta_new,mu_dummy,
     &                           dx,time+dt)
         call get_temp_visc_terms(scal_new,beta_new,
     &                            diff_hat(0,Temp),dx,time+dt)
         call get_spec_visc_terms(scal_new,beta_new,
     &                            diff_hat(0,FirstSpec),dx,time+dt)

         do i = 0,nx-1
            do n = 1,Nspec
               Y(n) = scal_new(i,FirstSpec+n-1) / scal_new(i,Density)
            enddo
            call CKCPBS(scal_new(i,Temp),Y,IWRK,RWRK,cpmix)
            rhocp = cpmix * scal_new(i,Density)
            diff_hat(i,Temp) = diff_hat(i,Temp)/rhocp
C     save a copy of diff_new(RhoH)
            diff_hat(i,RhoH) = diff_new(i,RhoH)
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

         print *,'... update rho'
         call update_rho(scal_old,scal_new,aofs,dx,dt)

c*****************************************************************

         print *,'... update D for species with A + R + MISDC(D)'
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
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
         else
            print *,'ERROR:: not set up to work with be_cn_theta=0.0d0'
            stop
         endif
 
         if (nochem_hack) then
            write(*,*)'WARNING:: SDC nochem_hack--skipping reactions'
         else
            print *,'... react with const and linear sources'
            do n = 1,nscal
               do i = 0,nx-1
                  
                  const_src(i,n) = aofs(i,n)
     $                 + diff_new(i,n) - diff_hat(i,n)
                  lin_src_old(i,n) = diff_old(i,n)
                  lin_src_new(i,n) = diff_hat(i,n)
               enddo
            enddo
            
            call strang_chem(scal_old,scal_new,
     $           const_src,lin_src_old,lin_src_new,
     $           I_R_new,dt)
         endif

c*****************************************************************
c       End of MISDC iterations
c*****************************************************************
C     CEG debugging FIXME
C     
C     Find the size of the correction, ie the change in S_new
C     
         do n = 1,nscal
            change_max(n) = 0.d0
            change_min(n) = 0.d0
         enddo
         do i = 0,nx-1
            do n = 1,nscal
               Schange(i,n) = scal_new(i,n) - Schange(i,n)
               change_max(n) = MAX(change_max(n),Schange(i,n))
               change_min(n) = MIN(change_min(n),Schange(i,n))
            enddo 
         enddo
         write(*,*)
         write(*,*)'Size of the correction (Change in S_new)'
         write(*,*)'index      min      max'
         do n = 1,nscal
C            write(*,*)n,change_min(n),change_max(n)
            write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
         enddo
         do i = 0,nx-1
            do n = 1,nscal
               Schange(i,n) = scal_new(i,n)
            enddo
         enddo

      enddo

      end 

      subroutine advance_radau (macvel,scal_old,scal_new,
     $                          I_R_old,I_R_new,
     $                          beta_old,beta_new,dx,dt,time)

      implicit none
      include 'spec.h'
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
      real*8  scal_tmp(-1:nx  ,nscal)
      real*8   I_R_new(0:nx-1,0:maxspec)
      real*8   I_R_old(0:nx-1,0:maxspec)
      real*8   I_R_tmp(0:nx-1,0:maxspec)
      real*8    macvel(0 :nx  )
      real*8      aofs(0 :nx-1,nscal)
      real*8  beta_old(-1:nx,nscal)
      real*8  beta_new(-1:nx,nscal)
      real*8  beta_tmp(-1:nx,nscal)
      real*8  mu_dummy(-1:nx)
      real*8   rhohalf(0 :nx-1)
      real*8    tforce(0 :nx-1,nscal)
      real*8      visc(0 :nx-1)
      real*8        cp(0 :nx-1)
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta
      real*8 theta
      
      real*8    diff_old(0:nx-1,nscal)
      real*8    diff_new(0:nx-1,nscal)
      real*8    diff_hat(0:nx-1,nscal)
      real*8    diff_tmp(0:nx-1,nscal)
      real*8   const_src(0:nx-1,nscal)
      real*8     lin_src(0:nx-1,nscal)
      
      integer i,n,ispec
      integer iunit
      
      real*8 divu_max
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8      dRhs(0:nx-1,0:maxspec)
      real*8 rhocp_old, rhocp
      real*8 Tmid
      real*8   pthermo(-1:nx  )
      real*8    Ydot_max, Y(maxspec)
      real*8 RWRK, cpmix, sum
      integer IWRK, is, rho_flag
      integer misdc

      real*8  dt_m

C CEG debugging FIXME
      real*8 ptherm(-1:nx)
      integer j
      real*8  Schange(-1:nx  ,nscal)
      real*8  change_max(nscal)
      real*8  change_min(nscal)


      be_cn_theta = 1.0d0

      print *,'... using RADAU quadrature'
      print *,'... evolving using temperature equation'
      print *,'... creating the diffusive terms with old data'

C Rhoh solve assumes reasonable values inside scal_tmp
      do n = 1, nscal
         do i=-1,nx
            scal_tmp(i,n) = scal_old(i,n)
         enddo
      enddo

C     CEG:: each one of these functions first calls set_bc(scal_old)
C     maybe should change this
      call get_spec_visc_terms(scal_old,beta_old,
     &     diff_old(0,FirstSpec),dx,time)
      call divRhoDHgradY(scal_old,beta_old,diff_old(0,RhoH),
     &     dx,time)
      call addDivLambdaGradT(scal_old,beta_old,
     &     diff_old(0,RhoH),dx,time)
c*****************************************************************
      
      print *,'... computing aofs with D(old) + R(guess)'

      do i = 0,nx-1
         do n = 1,Nspec
            is = FirstSpec + n - 1
            tforce(i,is) = diff_old(i,is) + I_R_new(i,n)
         enddo
         tforce(i,RhoH) = diff_old(i,RhoH)
      enddo
      
      call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)

c*****************************************************************
C     advance Diffusion in two steps

      dt_m = dt/3.d0
      print *,'... update rho to dt+1/3'
      call update_rho(scal_old,scal_tmp,aofs,dx,dt_m)

c*****************************************************************
      print *,'... do predictor for species (MISDC terms=0)'
      theta = 1.d0
      do i=0,nx-1
         dRhs(i,0) = 0.d0
         do n=1,Nspec
            dRhs(i,n) = dt_m*I_R_new(i,n)
         enddo
      enddo

      print *,'... do predictor for rhoh (MISDC terms=0)'
C CEG don't have new beta yet, just use old for right now
      call update_rhoh(scal_old,scal_tmp,aofs,alpha,beta_old,
     &     dRhs(0,0),Rhs(0,Temp),dx,dt_m,theta,time)
C     Implicit solve for Temp^n+1
      rho_flag = 1
      call cn_solve(scal_tmp,alpha,beta_old,Rhs(0,Temp),
     $     dx,dt_m,Temp,theta,rho_flag)
      call calc_diffusivities(scal_tmp,beta_tmp,mu_dummy,dx,time+dt_m)

      call update_spec(scal_old,scal_tmp,aofs,alpha,beta_old,
     &     dRhs(0,1),Rhs(0,FirstSpec),dx,dt_m,be_cn_theta,time)

      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_tmp,alpha,beta_tmp,Rhs(0,is),
     $        dx,dt_m,is,be_cn_theta,rho_flag)
      enddo

CCCCCCCCCCCCCCCCCCCCC
      call update_rhoh(scal_old,scal_tmp,aofs,alpha,beta_old,
     &     dRhs(0,0),Rhs(0,Temp),dx,dt_m,theta,time)
C     Implicit solve for Temp^n+1
      rho_flag = 1
      call cn_solve(scal_tmp,alpha,beta_tmp,Rhs(0,Temp),
     $     dx,dt_m,Temp,theta,rho_flag)
CCCCCCCCCCCCCCCCCcc

      print *,'...   extract D sources'
C     CEG;; note that neither of these 2 fns use rhoH_new
      call divRhoDHgradY(scal_tmp,beta_tmp,diff_new(0,RhoH),
     &     dx,time+dt_m)
      call addDivLambdaGradT(scal_tmp,beta_tmp,
     &     diff_new(0,RhoH),dx,time+dt_m)
      call get_spec_visc_terms(scal_tmp,beta_tmp,
     &     diff_new(0,FirstSpec),dx,time+dt_m)

      if (nochem_hack) then
         print *,'WARNING! doing nochem_hack--skipping reactions'
         do n = 0,Nspec
            do i = 0,nx-1
               I_R_tmp(i,n) = 0.d0
            enddo
         enddo
      else
         print *,'... react with A+D sources, reset I_R_new'
         do n = 1,nscal
            do i = 0,nx-1
               const_src(i,n) = aofs(i,n) + diff_new(i,n)
               lin_src(i,n) = 0.d0
C this doesn't help iwht the osc
C               const_src(i,n) = aofs(i,n) + diff_old(i,n)
C               lin_src(i,n) = (diff_new(i,n) - diff_old(i,n))/dt_m
            enddo
         enddo

         call chem(scal_old,scal_tmp,const_src,lin_src,
     $        I_R_tmp,dt_m)
      endif

      call calc_diffusivities(scal_tmp,beta_tmp,mu_dummy,
     &     dx,time+dt_m)

C     CEG debugging FIXME
C     
C     Find the estimated change in S over the timestep
C     
      do n = 1,nscal
         change_max(n) = 0.d0
         change_min(n) = 0.d0
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_tmp(i,n) - scal_old(i,n)
            change_max(n) = MAX(change_max(n),Schange(i,n))
            change_min(n) = MIN(change_min(n),Schange(i,n))
         enddo 
      enddo
      write(*,*)
      write(*,*)'Change in S over the timestep'
      write(*,*)'index      min      max'
      do n = 1,nscal
         write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_tmp(i,n)
         enddo
      enddo
 1008 FORMAT((I5,1X),(E22.15,1X))      

CCC debugging fixme
         call get_spec_visc_terms(scal_tmp,beta_tmp,
     &        diff_tmp(0,FirstSpec),dx,time+dt/3.d0)
         call divRhoDHgradY(scal_tmp,beta_tmp,diff_tmp(0,RhoH),
     &        dx,time+dt/3.d0)
         call addDivLambdaGradT(scal_tmp,beta_tmp,
     &        diff_tmp(0,RhoH),dx,time+dt/3.d0)

CCCCCCC
      dt_m = 2.d0*dt/3.d0
      print *,'... update rho from n+1/3 to n+1'
      call update_rho(scal_tmp,scal_new,aofs,dx,dt_m)

c*****************************************************************
      print *,'... do predictor for species (MISDC terms=0)'
      do i=0,nx-1
         do n=1,Nspec
            dRhs(i,n) = dt_m*I_R_new(i,n)
         enddo
      enddo

      print *,'... do predictor for rhoh (MISDC terms=0)'
C CEG don't have beta new yet, so use beta at n+1/3
      call update_rhoh(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &     dRhs(0,0),Rhs(0,Temp),dx,dt_m,theta,time+dt/3.d0)
C     Implicit solve for Temp^n+1
      rho_flag = 1
      call cn_solve(scal_new,alpha,beta_tmp,Rhs(0,Temp),
     $     dx,dt_m,Temp,theta,rho_flag)
      call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time+dt)

      call update_spec(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &     dRhs(0,1),Rhs(0,FirstSpec),dx,dt_m,be_cn_theta,time+dt/3.d0)

      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $        dx,dt_m,is,be_cn_theta,rho_flag)
      enddo

CCCCCCCCCCCCCCCCc
      call update_rhoh(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &     dRhs(0,0),Rhs(0,Temp),dx,dt_m,theta,time+dt/3.d0)
C     Implicit solve for Temp^n+1
      rho_flag = 1
      call cn_solve(scal_new,alpha,beta_new,Rhs(0,Temp),
     $     dx,dt_m,Temp,theta,rho_flag)
CCCCCCCCCCCCCCCCCCCCC

      print *,'...   extract D sources'
C     CEG;; note that neither of these 2 fns use rhoH_new
      call divRhoDHgradY(scal_new,beta_new,diff_new(0,RhoH),
     &     dx,time+dt)
      call addDivLambdaGradT(scal_new,beta_new,
     &     diff_new(0,RhoH),dx,time+dt)
      call get_spec_visc_terms(scal_new,beta_new,
     &     diff_new(0,FirstSpec),dx,time+dt)

      if (nochem_hack) then
         print *,'WARNING! doing nochem_hack--skipping reactions'
      else
         print *,'... react with A+D sources, reset I_R_new'
         do n = 1,nscal
            do i = 0,nx-1
               const_src(i,n) = aofs(i,n) + diff_new(i,n)
               lin_src(i,n) = 0.d0
C               const_src(i,n) = aofs(i,n) + diff_tmp(i,n)
C               lin_src(i,n) = (diff_new(i,n)-diff_tmp(i,n))/dt_m
            enddo
         enddo

         call chem(scal_tmp,scal_new,const_src,lin_src,
     $        I_R_new,dt_m)
      endif

C     CEG debugging FIXME
C     
C     Find the estimated change in S over the timestep
C     
      do n = 1,nscal
         change_max(n) = 0.d0
         change_min(n) = 0.d0
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_new(i,n) - scal_tmp(i,n)
            change_max(n) = MAX(change_max(n),Schange(i,n))
            change_min(n) = MIN(change_min(n),Schange(i,n))
         enddo 
      enddo
      write(*,*)
      write(*,*)'Change in S over the timestep'
      write(*,*)'index      min      max'
      do n = 1,nscal
         write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_new(i,n)
         enddo
      enddo

C----------------------------------------------------------------
C----------------------------------------------------------------
      do misdc = 1, misdc_iterMAX
         print *,'... doing SDC iter ',misdc

         print *,'... create new diff_hat from current state'
         call calc_diffusivities(scal_tmp,beta_tmp,mu_dummy,
     &        dx,time+dt/3.d0)
         call calc_diffusivities(scal_new,beta_new,mu_dummy,
     &        dx,time+dt)

         call get_spec_visc_terms(scal_tmp,beta_tmp,
     &        diff_tmp(0,FirstSpec),dx,time+dt/3.d0)
         call divRhoDHgradY(scal_tmp,beta_tmp,diff_tmp(0,RhoH),
     &        dx,time+dt/3.d0)
         call addDivLambdaGradT(scal_tmp,beta_tmp,
     &        diff_tmp(0,RhoH),dx,time+dt/3.d0)

         call get_spec_visc_terms(scal_new,beta_new,
     &        diff_hat(0,FirstSpec),dx,time+dt)
         call divRhoDHgradY(scal_new,beta_new,diff_hat(0,RhoH),
     &        dx,time+dt)
         call addDivLambdaGradT(scal_new,beta_new,
     &        diff_hat(0,RhoH),dx,time+dt)

         do i = 0,nx-1
            do n = 1,Nspec
               ispec = FirstSpec + n - 1
               tforce(i,ispec) = 
     $              I_R_tmp(i,n)/3.d0 + 2.d0*I_R_new(i,n)/3.d0
     &              + 0.25d0*(diff_hat(i,ispec)+3.d0*diff_tmp(i,ispec))
            enddo
            tforce(i,RhoH) =
     $           0.25d0*(diff_hat(i,RhoH)+3.d0*diff_tmp(i,RhoH))
         enddo
         
         print *,'... compute A with updated D+R source'
         call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)

C
C update diffusion with 2 steps
C  
         dt_m = dt/3.d0
         print *,'... update rho'
         call update_rho(scal_old,scal_tmp,aofs,dx,dt_m)

c*****************************************************************
         print *,'... update D for species with A + R + MISDC(D)'
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
               dRhs(i,n) = dt_m*(I_R_tmp(i,n) - diff_tmp(i,is)) 
     &              + dt*(5.d0*diff_tmp(i,is) - diff_hat(i,is))/12.d0
            enddo
            dRhs(i,0) = -dt_m * diff_tmp(i,RhoH) 
     &              + dt*(5.d0*diff_tmp(i,RhoH)-diff_hat(i,RhoH))/12.d0
         enddo

         call update_spec(scal_old,scal_tmp,aofs,alpha,beta_old,
     &        dRhs(0,1),Rhs(0,FirstSpec),dx,dt_m,be_cn_theta,time)
         rho_flag = 2
         do n=1,Nspec
            is = FirstSpec + n - 1
            call cn_solve(scal_tmp,alpha,beta_tmp,Rhs(0,is),
     $           dx,dt_m,is,be_cn_theta,rho_flag)
         enddo

         print *,'... update D for rhoh with A + R + MISDC(D)'
         call update_rhoh(scal_old,scal_tmp,aofs,alpha,beta_old,
     &        dRhs(0,0),Rhs(0,Temp),dx,dt_m,be_cn_theta,time)
C     Implicit solve for Temp^n+1
         rho_flag = 1
         call cn_solve(scal_tmp,alpha,beta_tmp,Rhs(0,Temp),
     $        dx,dt_m,Temp,be_cn_theta,rho_flag)

         print *,'...   extract D sources'
C     CEG;; note that neither of these 2 fns use rhoH_new
         call divRhoDHgradY(scal_tmp,beta_tmp,diff_new(0,RhoH),
     &        dx,time+dt_m)
         call addDivLambdaGradT(scal_tmp,beta_tmp,
     &        diff_new(0,RhoH),dx,time+dt_m)
         call get_spec_visc_terms(scal_tmp,beta_tmp,
     &        diff_new(0,FirstSpec),dx,time+dt_m)
         
         if (nochem_hack) then
            print *,'WARNING: SDC nochem_hack--skipping reactions'
         else
            print *,'... react with const and linear sources'
            do n = 1,nscal
               do i = 0,nx-1
                  const_src(i,n) = aofs(i,n)
     $                 + diff_new(i,n) - diff_tmp(i,n)
     $                 + (3.d0*diff_tmp(i,n)-diff_hat(i,n))/2.d0
                  lin_src(i,n) = (diff_hat(i,n)-diff_tmp(i,n))*1.5d0/dt
               enddo                     
            enddo
            
            call chem(scal_old,scal_tmp,const_src,lin_src,
     $           I_R_tmp,dt_m)
         endif
C
C update diffusion with 2 steps
C  
         dt_m = 2.d0*dt/3.d0
         print *,'... update rho'
         call update_rho(scal_tmp,scal_new,aofs,dx,dt_m)

c*****************************************************************
         print *,'... update D for species with A + R + MISDC(D)'
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
               dRhs(i,n) = dt_m*(I_R_new(i,n) - diff_hat(i,is)) 
     &              + dt*(diff_tmp(i,is) + diff_hat(i,is))/3.d0
            enddo
            dRhs(i,0) = -dt_m * diff_hat(i,RhoH) 
     &              + dt*(diff_tmp(i,RhoH) + diff_hat(i,RhoH))/3.d0
         enddo

         call update_spec(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &        dRhs(0,1),Rhs(0,FirstSpec),dx,dt_m,be_cn_theta,
     &        time+dt/3.d0)
         rho_flag = 2
         do n=1,Nspec
            is = FirstSpec + n - 1
            call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $           dx,dt_m,is,be_cn_theta,rho_flag)
         enddo

         print *,'... update D for rhoh with A + R + MISDC(D)'
         call update_rhoh(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &        dRhs(0,0),Rhs(0,Temp),dx,dt_m,be_cn_theta,time+dt/3.d0)
C     Implicit solve for Temp^n+1
         rho_flag = 1
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,Temp),
     $        dx,dt_m,Temp,be_cn_theta,rho_flag)

         print *,'...   extract D sources'
C     CEG;; note that neither of these 2 fns use rhoH_new
         call divRhoDHgradY(scal_new,beta_new,diff_new(0,RhoH),
     &        dx,time+dt)
         call addDivLambdaGradT(scal_new,beta_new,
     &        diff_new(0,RhoH),dx,time+dt)
         call get_spec_visc_terms(scal_new,beta_new,
     &        diff_new(0,FirstSpec),dx,time+dt)
         
         if (nochem_hack) then
            print *,'WARNING: SDC nochem_hack--skipping reactions'
         else
            print *,'... react with const and linear sources'
            do n = 1,nscal
               do i = 0,nx-1
                  const_src(i,n) = aofs(i,n)
     $                 + diff_new(i,n) - diff_hat(i,n)
     $                 + diff_tmp(i,n)
                  lin_src(i,n) = (diff_hat(i,n)-diff_tmp(i,n))*1.5d0/dt
               enddo                     
            enddo
            
            call chem(scal_tmp,scal_new,const_src,lin_src,
     $           I_R_new,dt_m)
         endif

c*****************************************************************
c     End of MISDC iterations
c*****************************************************************
C     CEG debugging FIXME
C     
C     Find the size of the correction, ie the change in S_new
C     
         do n = 1,nscal
            change_max(n) = 0.d0
            change_min(n) = 0.d0
         enddo
         do n = 1,nscal
            do i = 0,nx-1
               Schange(i,n) = scal_new(i,n) - Schange(i,n)
               change_max(n) = MAX(change_max(n),Schange(i,n))
               change_min(n) = MIN(change_min(n),Schange(i,n))
            enddo 
         enddo
         write(*,*)
         write(*,*)'Size of the correction (Change in S_new)'
         write(*,*)'index      min      max'
         do n = 1,nscal
C     write(*,*)n,change_min(n),change_max(n)
            write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
         enddo
         do i = 0,nx-1
            do n = 1,nscal
               Schange(i,n) = scal_new(i,n)
            enddo
         enddo

      enddo

      do n = 1,Nspec
         do i = 0,nx-1
            I_R_new(i,n) = I_R_tmp(i,n)/3.d0 + 2.d0*I_R_new(i,n)/3.d0 
         enddo
      enddo

      end 

      subroutine advance_radau_temp (macvel,scal_old,scal_new,
     $                          I_R_old,I_R_new,
     $                          beta_old,beta_new,dx,dt,time)

      implicit none
      include 'spec.h'
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
      real*8  scal_tmp(-1:nx  ,nscal)
      real*8   I_R_new(0:nx-1,0:maxspec)
      real*8   I_R_old(0:nx-1,0:maxspec)
      real*8   I_R_tmp(0:nx-1,0:maxspec)
      real*8    macvel(0 :nx  )
      real*8      aofs(0 :nx-1,nscal)
      real*8  beta_old(-1:nx,nscal)
      real*8  beta_new(-1:nx,nscal)
      real*8  beta_tmp(-1:nx,nscal)
      real*8  mu_dummy(-1:nx)
      real*8   rhohalf(0 :nx-1)
      real*8    tforce(0 :nx-1,nscal)
      real*8      visc(0 :nx-1)
      real*8        cp(0 :nx-1)
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta
      real*8 theta
      
      real*8    diff_old(0:nx-1,nscal)
      real*8    diff_new(0:nx-1,nscal)
      real*8    diff_hat(0:nx-1,nscal)
      real*8    diff_tmp(0:nx-1,nscal)
      real*8   const_src(0:nx-1,nscal)
      real*8     lin_src(0:nx-1,nscal)
      
      integer i,n,ispec
      integer iunit
      
      real*8 divu_max
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8      dRhs(0:nx-1,0:maxspec)
      real*8 rhocp_old, rhocp
      real*8 Tmid
      real*8   pthermo(-1:nx  )
      real*8    Ydot_max, Y(maxspec)
      real*8 RWRK, cpmix, sum
      integer IWRK, is, rho_flag
      integer misdc

      real*8  dt_m, HK(Nspec)

C CEG debugging FIXME
      real*8 ptherm(-1:nx)
      integer j
      real*8  Schange(-1:nx  ,nscal)
      real*8  change_max(nscal)
      real*8  change_min(nscal)


      be_cn_theta = 1.0d0

      print *,'... using RADAU quadrature WITH TEMP EQN'
      print *,'... evolving using temperature equation'
      print *,'... creating the diffusive terms with old data'

C Rhoh solve assumes reasonable values inside scal_tmp
      do n = 1, nscal
         do i=-1,nx
            scal_tmp(i,n) = scal_old(i,n)
         enddo
      enddo

C     CEG:: each one of these functions first calls set_bc(scal_old)
C     maybe should change this
      call get_temp_visc_terms(scal_old,beta_old,
     &                         diff_old(0,Temp),dx,time)
      call get_spec_visc_terms(scal_old,beta_old,
     &                         diff_old(0,FirstSpec),dx,time)
      call get_rhoh_visc_terms(scal_old,beta_old,
     &                         diff_old(0,RhoH),dx,time)
      
c*****************************************************************
      
      print *,'... computing aofs with D(old) + R(guess)'

      do i = 0,nx-1
         do n = 1,Nspec
            Y(n) = scal_old(i,FirstSpec+n-1) / scal_old(i,Density)
         enddo
         call CKCPBS(scal_old(i,Temp),Y,IWRK,RWRK,cpmix)
         rhocp_old = cpmix * scal_old(i,Density)
         diff_old(i,Temp) = diff_old(i,Temp)/rhocp_old

         I_R_new(i,0) = 0.d0
         call CKHMS(scal_old(i,Temp),IWRK,RWRK,HK)

         do n = 1,Nspec
            is = FirstSpec + n - 1
            tforce(i,is) = diff_old(i,is) + I_R_new(i,n)
            I_R_new(i,0) = I_R_new(i,0) +
     &           I_R_new(i,n)*HK(n)
         enddo
         tforce(i,Temp) = diff_old(i,Temp) + I_R_new(i,0)/rhocp_old
      enddo
       
      call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)

c*****************************************************************

      dt_m = dt/3.d0
      print *,'... update rho to dt+1/3'
      call update_rho(scal_old,scal_tmp,aofs,dx,dt_m)

c*****************************************************************
c     Either do c-n solve for new T prior to computing new 
c     coeffs, or simply start by copying from previous time step

      if (predict_temp_for_coeffs .eq. 1) then
         print *,'... predict temp with old coeffs'
         rho_flag = 1
         theta = 0.5d0
C CEG:: beta_new comes into advance() with the same value as 
C       beta_old (the one we just calculated)
         call update_temp(scal_old,scal_tmp,aofs,
     $                    alpha,beta_old,beta_old,I_R_new(0,0),
     $                    Rhs(0,Temp),dx,dt_m,theta,time)
C CEG:: just uses RHS and overwrites snew
C does not fill ghost cells
         call cn_solve(scal_tmp,alpha,beta_old,Rhs(0,Temp),
     $                 dx,dt_m,Temp,theta,rho_flag)

         call get_hmix_given_T_RhoY(scal_tmp,dx)      

         print *,'... compute new coeffs'
         call calc_diffusivities(scal_tmp,beta_tmp,mu_dummy,dx,
     &        time+dt_m)
      else
         print *,'... set new coeffs to old values for predictor'
         do n=1,nscal
            do i=-1,nx
               scal_tmp(i,Temp) = scal_old(i,Temp)
               beta_tmp(i,n) = beta_old(i,n)
            enddo
         enddo
      endif

c*****************************************************************
      print *,'... do predictor for species (MISDC terms=0)'
      do i=0,nx-1
         dRhs(i,0) = 0.0d0
         do n=1,Nspec
            dRhs(i,n) = dt_m*I_R_new(i,n)
         enddo
      enddo

      call update_spec(scal_old,scal_tmp,aofs,alpha,beta_old,
     &     dRhs(0,1),Rhs(0,FirstSpec),dx,dt_m,be_cn_theta,time)
      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_tmp,alpha,beta_tmp,Rhs(0,is),
     $        dx,dt_m,is,be_cn_theta,rho_flag)
      enddo

      print *,'... do predictor for rhoh (MISDC terms=0)'
      call update_rhoh(scal_old,scal_tmp,aofs,alpha,beta_old,
     &     dRhs(0,0),Rhs(0,RhoH),dx,dt_m,be_cn_theta,time)
      rho_flag = 2
      call cn_solve(scal_tmp,alpha,beta_tmp,Rhs(0,RhoH),
     $              dx,dt_m,RhoH,be_cn_theta,rho_flag)

      call rhoh_to_temp(scal_tmp)

      call get_spec_visc_terms(scal_tmp,beta_tmp,
     &                         diff_new(0,FirstSpec),dx,time+dt_m)
      call get_rhoh_visc_terms(scal_tmp,beta_tmp,
     &                         diff_new(0,RhoH),dx,time+dt_m)

      if (nochem_hack) then
         print *,'WARNING! doing nochem_hack--skipping reactions'
         do n = 0,Nspec
            do i = 0,nx-1
               I_R_tmp(i,n) = 0.d0
            enddo
         enddo
      else
CCC debugging FIXME!!!
C$$$        do i=0,nx-1
C$$$            diff_hat(i,RhoH) = diff_new(i,RhoH)
C$$$            do n=1,Nspec
C$$$               dRhs(i,n) = dt*I_R_new(i,n)
C$$$            enddo
C$$$         enddo
C$$$         call update_spec(scal_old,scal_new,aofs,alpha,beta_old,
C$$$     &        dRhs(0,1),Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)
C$$$         rho_flag = 2
C$$$         do n=1,Nspec
C$$$            is = FirstSpec + n - 1
C$$$            call cn_solve(scal_new,alpha,beta_tmp,Rhs(0,is),
C$$$     $           dx,dt,is,be_cn_theta,rho_flag)
C$$$         enddo
C$$$         call get_spec_visc_terms(scal_new,beta_tmp,
C$$$     &        diff_hat(0,FirstSpec),dx,time+dt)

C$$$CCCC   
         print *,'... react with A+D sources, reset I_R_new'
         do n = 1,nscal
            do i = 0,nx-1
               const_src(i,n) = aofs(i,n) + diff_new(i,n)
               lin_src(i,n) = 0.d0

C               const_src(i,n) = aofs(i,n) + 1.5d0*diff_tmp(i,n)
C     $              - 0.5d0*diff_hat(i,n)
C               lin_src(i,n) = (diff_hat(i,n) - diff_new(i,n))*1.5d0/dt
            enddo
         enddo

         call chem(scal_old,scal_tmp,const_src,lin_src,
     $        I_R_tmp,dt_m)
      endif

      call calc_diffusivities(scal_tmp,beta_tmp,mu_dummy,
     &     dx,time+dt_m)         

C     CEG debugging FIXME
C     
C     Find the estimated change in S over the timestep
C     
      do n = 1,nscal
         change_max(n) = 0.d0
         change_min(n) = 0.d0
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_tmp(i,n) - scal_old(i,n)
            change_max(n) = MAX(change_max(n),Schange(i,n))
            change_min(n) = MIN(change_min(n),Schange(i,n))
         enddo 
      enddo
      write(*,*)
      write(*,*)'Change in S over the timestep'
      write(*,*)'index      min      max'
      do n = 1,nscal
         write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_tmp(i,n)
         enddo
      enddo
 1008 FORMAT((I5,1X),(E22.15,1X))      

c*****************************************************************

      dt_m = 2.d0*dt/3.d0
      print *,'... update rho from n+1/3 to n+1'
      call update_rho(scal_tmp,scal_new,aofs,dx,dt_m)

      if (predict_temp_for_coeffs .eq. 1) then
         print *,'... predict temp with old coeffs'
         rho_flag = 1
         theta = 0.5d0
C CEG:: beta_new comes into advance() with the same value as 
C       beta_old (the one we just calculated)
         call update_temp(scal_tmp,scal_new,aofs,
     $                    alpha,beta_tmp,beta_tmp,I_R_new(0,0),
     $                    Rhs(0,Temp),dx,dt_m,theta,time+dt/3.d0)
C CEG:: just uses RHS and overwrites snew
C does not fill ghost cells
         call cn_solve(scal_new,alpha,beta_tmp,Rhs(0,Temp),
     $                 dx,dt_m,Temp,theta,rho_flag)

         call get_hmix_given_T_RhoY(scal_new,dx)      

         print *,'... compute new coeffs'
         call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,
     &        time+dt)
      else
         print *,'... set new coeffs to old values for predictor'
         do n=1,nscal
            do i=-1,nx
               scal_new(i,Temp) = scal_tmp(i,Temp)
               beta_new(i,n) = beta_tmp(i,n)
            enddo
         enddo
      endif

c*****************************************************************
      print *,'... do predictor for species (MISDC terms=0)'
      do i=0,nx-1
         dRhs(i,0) = 0.0d0
         do n=1,Nspec
            dRhs(i,n) = dt_m*I_R_new(i,n)
         enddo
      enddo

      call update_spec(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &     dRhs(0,1),Rhs(0,FirstSpec),dx,dt_m,be_cn_theta,time+dt/3.d0)

      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $        dx,dt_m,is,be_cn_theta,rho_flag)
      enddo

      print *,'... do predictor for rhoh (MISDC terms=0)'
      call update_rhoh(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &     dRhs(0,0),Rhs(0,RhoH),dx,dt_m,be_cn_theta,time+dt/3.d0)
      rho_flag = 2
      call cn_solve(scal_new,alpha,beta_new,Rhs(0,RhoH),
     $     dx,dt_m,RhoH,be_cn_theta,rho_flag)

      call rhoh_to_temp(scal_new)

      call get_spec_visc_terms(scal_new,beta_new,
     &     diff_new(0,FirstSpec),dx,time+dt)
      call get_rhoh_visc_terms(scal_new,beta_new,
     &                         diff_new(0,RhoH),dx,time+dt)

CCCCCCCCCCC debugging FIXME
C$$$ 1006 FORMAT((I5,1X),11(E22.15,1X))      
C$$$         call compute_pthermo(scal_new,ptherm)
C$$$         open(UNIT=11, FILE='snew.dat', STATUS = 'REPLACE')
C$$$         write(11,*)'# 256 12'
C$$$         do j=0,nx-1
C$$$            do n = 1,Nspec
C$$$               Y(n) = scal_new(j,FirstSpec+n-1)/scal_new(j,Density)
C$$$            enddo
C$$$            write(11,1006) j, macvel(j),
C$$$     &                     scal_new(j,Density),
C$$$     &                     (Y(n),n=1,Nspec),
C$$$     $                     scal_new(j,RhoH),
C$$$     $                     scal_new(j,Temp),
C$$$     $                     ptherm(j)
C$$$         enddo
C$$$         close(11)
C$$$         write(*,*)'AFTER end'
CCCCCCCCCCCCC      

      if (nochem_hack) then
         print *,'WARNING! doing nochem_hack--skipping reactions'
      else
         print *,'... react with A+D sources, reset I_R_new'
         do n = 1,nscal
            do i = 0,nx-1
C               const_src(i,n) = aofs(i,n) + diff_new(i,n)
C               lin_src(i,n) = 0.d0

               const_src(i,n) = aofs(i,n) + diff_tmp(i,n)
               lin_src(i,n) = (diff_new(i,n) - diff_tmp(i,n))/dt_m
            enddo
         enddo

         call chem(scal_tmp,scal_new,const_src,lin_src,
     $        I_R_new,dt_m)
      endif

C     CEG debugging FIXME
C     
C     Find the estimated change in S over the timestep
C     
      do n = 1,nscal
         change_max(n) = 0.d0
         change_min(n) = 0.d0
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_new(i,n) - scal_tmp(i,n)
            change_max(n) = MAX(change_max(n),Schange(i,n))
            change_min(n) = MIN(change_min(n),Schange(i,n))
         enddo 
      enddo
      write(*,*)
      write(*,*)'Change in S over the timestep'
      write(*,*)'index      min      max'
      do n = 1,nscal
         write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
      enddo
      do i = 0,nx-1
         do n = 1,nscal
            Schange(i,n) = scal_new(i,n)
         enddo
      enddo

C----------------------------------------------------------------
C----------------------------------------------------------------
      do misdc = 1, misdc_iterMAX
         print *,'... doing SDC iter ',misdc

         print *,'... create new diff_hat from current state'
         call calc_diffusivities(scal_tmp,beta_tmp,mu_dummy,
     &        dx,time+dt/3.d0)
         call calc_diffusivities(scal_new,beta_new,mu_dummy,
     &        dx,time+dt)

         call get_spec_visc_terms(scal_tmp,beta_tmp,
     &        diff_tmp(0,FirstSpec),dx,time+dt/3.d0)
         call get_rhoh_visc_terms(scal_tmp,beta_tmp,
     &        diff_tmp(0,RhoH),dx,time+dt/3.d0)
         call get_temp_visc_terms(scal_tmp,beta_tmp,
     &        diff_tmp(0,Temp),dx,time+dt/3.d0)


         call get_spec_visc_terms(scal_new,beta_new,
     &        diff_hat(0,FirstSpec),dx,time+dt)
         call get_rhoh_visc_terms(scal_new,beta_new,
     &        diff_hat(0,RhoH),dx,time+dt)
         call get_temp_visc_terms(scal_new,beta_new,
     &        diff_hat(0,Temp),dx,time+dt)

         do i = 0,nx-1
            do n = 1,Nspec
               Y(n) = scal_new(i,FirstSpec+n-1) / scal_new(i,Density)
            enddo
            call CKCPBS(scal_new(i,Temp),Y,IWRK,RWRK,cpmix)
            rhocp = cpmix * scal_new(i,Density)
            diff_hat(i,Temp) = diff_hat(i,Temp)/rhocp
            do n = 1,Nspec
               Y(n) = scal_tmp(i,FirstSpec+n-1) / scal_tmp(i,Density)
            enddo
            call CKCPBS(scal_tmp(i,Temp),Y,IWRK,RWRK,cpmix)
            rhocp = cpmix * scal_tmp(i,Density)
            diff_tmp(i,Temp) = diff_tmp(i,Temp)/rhocp

            tforce(i,Temp) = 0.d0
            call CKHMS(scal_tmp(i,Temp),IWRK,RWRK,HK)
            do n = 1,Nspec
               ispec = FirstSpec + n - 1
               tforce(i,ispec) = 
     $              I_R_tmp(i,n)/3.d0 + 2.d0*I_R_new(i,n)/3.d0
     &              + 0.25d0*(diff_hat(i,ispec)+3.d0*diff_tmp(i,ispec))
               tforce(i,Temp) = tforce(i,Temp) +
     &              I_R_tmp(i,n)*HK(n)
            enddo
            tforce(i,Temp) = tforce(i,Temp)/rhocp 
     &           + 0.25d0*(diff_hat(i,Temp)+3.d0*diff_tmp(i,Temp))

         enddo
         
         print *,'... compute A with updated D+R source'
         call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)

C
C update diffusion with 2 steps
C  
         dt_m = dt/3.d0
         print *,'... update rho'
         call update_rho(scal_old,scal_tmp,aofs,dx,dt_m)

c*****************************************************************
         print *,'... update D for species with A + R + MISDC(D)'
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
               dRhs(i,n) = dt_m*(I_R_tmp(i,n) - diff_tmp(i,is)) 
     &              + dt*(5.d0*diff_tmp(i,is) - diff_hat(i,is))/12.d0
            enddo
            dRhs(i,0) = -dt_m * diff_tmp(i,RhoH) 
     &              + dt*(5.d0*diff_tmp(i,RhoH)-diff_hat(i,RhoH))/12.d0
         enddo

         call update_spec(scal_old,scal_tmp,aofs,alpha,beta_old,
     &        dRhs(0,1),Rhs(0,FirstSpec),dx,dt_m,be_cn_theta,time)
         rho_flag = 2
         do n=1,Nspec
            is = FirstSpec + n - 1
            call cn_solve(scal_tmp,alpha,beta_tmp,Rhs(0,is),
     $           dx,dt_m,is,be_cn_theta,rho_flag)
         enddo

         print *,'... update D for rhoh with A + R + MISDC(D)'
         call update_rhoh(scal_old,scal_tmp,aofs,alpha,beta_old,
     &        dRhs(0,0),Rhs(0,RhoH),dx,dt_m,be_cn_theta,time)
C     Implicit solve for Temp^n+1
         rho_flag = 2
         call cn_solve(scal_tmp,alpha,beta_tmp,Rhs(0,RhoH),
     $        dx,dt_m,RhoH,be_cn_theta,rho_flag)

CCCCCCCCCCCC
         call rhoh_to_temp(scal_tmp)
C         call calc_diffusivities(scal_tmp,beta_tmp,mu_dummy,dx,
C     &        time+dt_m)
CCCCCCCCCCCCCc

         print *,'...   extract D sources'
C     CEG;; note that neither of these 2 fns use rhoH_new
         call get_spec_visc_terms(scal_tmp,beta_tmp,
     &        diff_new(0,FirstSpec),dx,time+dt_m)
         call get_rhoh_visc_terms(scal_tmp,beta_tmp,
     &        diff_new(0,RhoH),dx,time+dt_m)

         if (nochem_hack) then
            print *,'WARNING: SDC nochem_hack--skipping reactions'
         else
            print *,'... react with const and linear sources'
            do n = 1,nscal
               do i = 0,nx-1
                  const_src(i,n) = aofs(i,n)
     $                 + diff_new(i,n) - diff_tmp(i,n)
     $                 + (3.d0*diff_tmp(i,n)-diff_hat(i,n))/2.d0
                  lin_src(i,n) = (diff_hat(i,n)-diff_tmp(i,n))*1.5d0/dt
               enddo                     
            enddo
            
            call chem(scal_old,scal_tmp,const_src,lin_src,
     $           I_R_tmp,dt_m)
         endif
C
C update diffusion with 2 steps
C  
         dt_m = 2.d0*dt/3.d0
         print *,'... update rho'
         call update_rho(scal_tmp,scal_new,aofs,dx,dt_m)

c*****************************************************************
         print *,'... update D for species with A + R + MISDC(D)'
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
               dRhs(i,n) = dt_m*(I_R_new(i,n) - diff_hat(i,is)) 
     &              + dt*(diff_tmp(i,is) + diff_hat(i,is))/3.d0
            enddo
            dRhs(i,0) = -dt_m * diff_hat(i,RhoH) 
     &              + dt*(diff_tmp(i,RhoH) + diff_hat(i,RhoH))/3.d0
         enddo

         call update_spec(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &        dRhs(0,1),Rhs(0,FirstSpec),dx,dt_m,be_cn_theta,
     &        time+dt/3.d0)
         rho_flag = 2
         do n=1,Nspec
            is = FirstSpec + n - 1
            call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $           dx,dt_m,is,be_cn_theta,rho_flag)
         enddo

         print *,'... update D for rhoh with A + R + MISDC(D)'
         call update_rhoh(scal_tmp,scal_new,aofs,alpha,beta_tmp,
     &        dRhs(0,0),Rhs(0,RhoH),dx,dt_m,be_cn_theta,time+dt/3.d0)
C     Implicit solve for Temp^n+1
         rho_flag = 2
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,RhoH),
     $        dx,dt_m,RhoH,be_cn_theta,rho_flag)

CCCCCCCCCCCC
         call rhoh_to_temp(scal_new)
C         call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,
C     &        time+dt)
CCCCCCCCCCCCCc

         print *,'...   extract D sources'
C     CEG;; note that neither of these 2 fns use rhoH_new
         call get_spec_visc_terms(scal_new,beta_new,
     &        diff_new(0,FirstSpec),dx,time+dt)
         call get_rhoh_visc_terms(scal_new,beta_new,
     &        diff_new(0,RhoH),dx,time+dt)
         
         if (nochem_hack) then
            print *,'WARNING: SDC nochem_hack--skipping reactions'
         else
            print *,'... react with const and linear sources'
            do n = 1,nscal
               do i = 0,nx-1
                  const_src(i,n) = aofs(i,n)
     $                 + diff_new(i,n) - diff_hat(i,n)
     $                 + diff_tmp(i,n)
                  lin_src(i,n) = (diff_hat(i,n)-diff_tmp(i,n))*1.5d0/dt
               enddo                     
            enddo
            
            call chem(scal_tmp,scal_new,const_src,lin_src,
     $           I_R_new,dt_m)
         endif

c*****************************************************************
c     End of MISDC iterations
c*****************************************************************
C     CEG debugging FIXME
C     
C     Find the size of the correction, ie the change in S_new
C     
         do n = 1,nscal
            change_max(n) = 0.d0
            change_min(n) = 0.d0
         enddo
         do i = 0,nx-1
            do n = 1,nscal
               Schange(i,n) = scal_new(i,n) - Schange(i,n)
               change_max(n) = MAX(change_max(n),Schange(i,n))
               change_min(n) = MIN(change_min(n),Schange(i,n))
            enddo 
         enddo
         write(*,*)
         write(*,*)'Size of the correction (Change in S_new)'
         write(*,*)'index      min      max'
         do n = 1,nscal
C     write(*,*)n,change_min(n),change_max(n)
            write(*,*)n,MAX(ABS(change_min(n)),ABS(change_max(n)))
         enddo
         do i = 0,nx-1
            do n = 1,nscal
               Schange(i,n) = scal_new(i,n)
            enddo
         enddo

      enddo

      do n = 1,Nspec
         do i = 0,nx-1
            I_R_new(i,n) = I_R_tmp(i,n)/3.d0 + 2.d0*I_R_new(i,n)/3.d0 
         enddo
      enddo

      end 

      subroutine strang_advance (macvel,scal_old,scal_new,
     $                   I_R_old,I_R_new,
     $                   beta_old,beta_new,dx,dt,time)

      implicit none
      include 'spec.h'
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
      real*8   I_R_new(0:nx-1,0:maxspec)
      real*8   I_R_old(0:nx-1,0:maxspec)
      real*8    macvel(0 :nx  )
      real*8      aofs(0 :nx-1,nscal)
      real*8  beta_old(-1:nx,nscal)
      real*8  beta_new(-1:nx,nscal)
      real*8  mu_dummy(-1:nx)
      real*8   rhohalf(0 :nx-1)
      real*8    tforce(0 :nx-1,nscal)
      real*8      visc(0 :nx-1)
      real*8        cp(0 :nx-1)
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta
      real*8 theta
      
      real*8    diff_old(0:nx-1,nscal)
      real*8   const_src(0:nx-1,nscal)
      real*8 lin_src_old(0:nx-1,nscal)
      real*8 lin_src_new(0:nx-1,nscal)
      
      integer i,n,ispec
      integer iunit
      
      real*8 divu_max
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8      dRhs(0:nx-1,0:maxspec)
      real*8 rhocp_old
      real*8 Tmid
      real*8   pthermo(-1:nx  )
      real*8    Ydot_max, Y(maxspec)
      real*8 RWRK, cpmix, sum
      integer IWRK, is, rho_flag

C CEG debugging FIXME
      real*8 ptherm(-1:nx)
      integer j
      

      be_cn_theta = 0.5d0

      if (nochem_hack) then
         print *,'WARNING! doing nochem_hack...'
         do n = 1,nscal
            do i = 0,nx-1
               I_R_new(i,n) = 0.d0
            enddo
         enddo
      else
         print *,'... react for dt/2;  set I_R_new'
         do n = 1,nscal
            do i = 0,nx-1
               const_src(i,n) =   0.d0
               lin_src_old(i,n) = 0.d0
               lin_src_new(i,n) = 0.d0
            enddo
         enddo
         call strang_chem(scal_old,scal_new,
     $                    const_src,lin_src_old,lin_src_new,
     $                    I_R_new,dt/2.d0)
         do j=0,nx-1
            do n = FirstSpec,LastSpec
               scal_old(j,n) = scal_new(j,n)
            enddo
            scal_old(j,Temp) = scal_new(j,Temp)
         enddo
      endif
c     
c*****************************************************************
c     
      print *,'... creating the diffusive terms with old data'
C CEG:: fixme??
      call calc_diffusivities(scal_old,beta_old,mu_dummy,dx,time)
      call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time)

C CEG:: each one of these functions first calls set_bc(scal_old)
C   maybe should change this
      call get_temp_visc_terms(scal_old,beta_old,
     &                         diff_old(0,Temp),dx,time)
      call get_spec_visc_terms(scal_old,beta_old,
     &                         diff_old(0,FirstSpec),dx,time)
      call get_rhoh_visc_terms(scal_old,beta_old,
     &                         diff_old(0,RhoH),dx,time)
            
      print *,'... computing aofs with D(old)'

      do i = 0,nx-1
         do n = 1,Nspec
            Y(n) = scal_old(i,FirstSpec+n-1) / scal_old(i,Density)
         enddo
         call CKCPBS(scal_old(i,Temp),Y,IWRK,RWRK,cpmix)
         rhocp_old = cpmix * scal_old(i,Density)
         diff_old(i,Temp) = diff_old(i,Temp)/rhocp_old

         do n = 1,Nspec
            is = FirstSpec + n - 1
            tforce(i,is) = diff_old(i,is)
         enddo
         tforce(i,Temp) = diff_old(i,Temp)
      enddo
       
      call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt,time)

      print *,'... update rho'
      call update_rho(scal_old,scal_new,aofs,dx,dt)

      call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time)

c*****************************************************************
c     Either do c-n solve for new T prior to computing new 
c     coeffs, or simply start by copying from previous time step
      if (predict_temp_for_coeffs .eq. 1) then
         print *,'... predict temp with old coeffs'
         rho_flag = 1
         theta = 0.5d0
         do i=0,nx-1
            dRhs(i,0) = 0.0d0
         enddo
         call update_temp(scal_old,scal_new,aofs,
     $                    alpha,beta_old,beta_new,dRhs(0,0),
     $                    Rhs(0,Temp),dx,dt,theta,time)
C CEG:: just uses RHS and overwrites snew
C does not fill ghost cells
         call cn_solve(scal_new,alpha,beta_old,Rhs(0,Temp),
     $                 dx,dt,Temp,theta,rho_flag)

         call get_hmix_given_T_RhoY(scal_new,dx)      

         print *,'... compute new coeffs'
         call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time+dt)
      else
         print *,'... set new coeffs to old values for predictor'
         do n=1,nscal
            do i=-1,nx
               scal_new(i,Temp) = scal_old(i,Temp)
               beta_new(i,n) = beta_old(i,n)
            enddo
         enddo
      endif

      print *,'... do predictor for species'
      do i=0,nx-1
         dRhs(i,0) = 0.0d0
         do n=1,Nspec
            dRhs(i,n) = 0.d0
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

C----------------------------------------------------------------
C   Corrector

      print *,'... compute new coeffs'
      call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time+dt)

      call update_spec(scal_old,scal_new,aofs,alpha,beta_old,
     &                 dRhs(0,1),Rhs(0,FirstSpec),dx,dt,
     &                 be_cn_theta,time)
      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta,rho_flag)
      enddo
         
      print *,'... do predictor for rhoh (MISDC terms=0)'
      call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,
     &                 dRhs(0,0),Rhs(0,RhoH),dx,dt,be_cn_theta,time)
      rho_flag = 2
      call cn_solve(scal_new,alpha,beta_new,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta,rho_flag)
      call rhoh_to_temp(scal_new)

      if (nochem_hack) then
         print *,'WARNING! doing nochem_hack...'
      else
         do j=0,nx-1
            do n = FirstSpec,LastSpec
               scal_old(j,n) = scal_new(j,n)
            enddo
            scal_old(j,Temp) = scal_new(j,Temp)
            scal_old(j,Density) = scal_new(j,Density)
         enddo

         call strang_chem(scal_old,scal_new,
     $                    const_src,lin_src_old,lin_src_new,
     $                    I_R_new,dt/2.d0)
      endif

CCCCCCCCCCC debugging FIXME
C$$$ 1006 FORMAT((I5,1X),11(E22.15,1X))      
C$$$         call compute_pthermo(scal_new,ptherm)
C$$$         open(UNIT=11, FILE='corr.dat', STATUS = 'REPLACE')
C$$$         write(11,*)'# 256 12'
C$$$         do j=0,nx-1
C$$$            do n = 1,Nspec
C$$$               Y(n) = scal_new(j,FirstSpec+n-1)*1.d3
C$$$            enddo
C$$$            write(11,1006) j, macvel(j)*1.d-2, 
C$$$     &                     scal_new(j,Density)*1.d3,
C$$$     &                     (Y(n),n=1,Nspec),
C$$$     $                     scal_new(j,RhoH)*1.d-1,
C$$$     $                     scal_new(j,Temp),
C$$$     $                     ptherm(j)*1.d-1
C$$$         enddo
C$$$         close(11)
C$$$         write(*,*)'AFTER end'
CCCCCCCCCCCCC      

      end
