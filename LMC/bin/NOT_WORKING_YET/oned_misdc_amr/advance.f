      subroutine advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R,press_old,press_new,
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
      real*8   I_R(0:nx-1,0:maxspec)
      real*8    macvel(0 :nx  )
      real*8   veledge(0 :nx  )
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
      real*8      visc(0 :nx-1)
      real*8 dx
      real*8 dt
      real*8 time
      real*8 vel_theta
      real*8 divu_max
      
      integer i,j
      
      real*8     alpha(0:nx-1)
      real*8   vel_Rhs(0:nx-1)
      real*8   pthermo(-1:nx  )

      real*8   I_R_divu(0:nx-1,0:maxspec)
      real*8 WDOTK(maxspec), C(maxspec), RWRK
      integer IWRK

      integer rho_flag
      
      print *,'advance: at start of time step'

c     
c*****************************************************************
c     Create MAC velocities.
c*****************************************************************
c     
      do i = 0,nx-1
         gp(i) = (press_old(i+1) - press_old(i)) / dx
      enddo

      print *,'... predict edge velocities'
c     this fills ghost cells for vel_old
      call pre_mac_predict(vel_old,scal_old,gp,macvel,dx,dt,time)
      
      call compute_pthermo(scal_old,pthermo)

      do i = 0,nx-1
         divu_tmp(i) = divu_old(i) + 0.5d0 * dt * dsdt(i)
      enddo

c     diagnostics only
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU_TMP norm before dpdt = ',divu_max 

      call add_dpdt(scal_old,pthermo,divu_tmp,macvel,dx,dt)

c     diagnostics only
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU_TMP norm after dpdt = ',divu_max 

      call macproj(nx,macvel,divu_tmp,dx)

c     compute diffusivities at time n (old time)
c     this computes rho D_m     (for species)
c                   lambda / cp (for enthalpy)
c                   lambda      (for temperature)
      call calc_diffusivities(scal_old,beta_old,mu_old,dx,time)

      if (use_strang) then

c     Strang split advance
         call strang_advance(macvel,scal_old,scal_new,
     $                   I_R,beta_old,beta_new,
     $                   dx,dt,time)

      else

c     SDC advance
         call sdc_advance(macvel,scal_old,scal_new,
     $                    I_R,beta_old,beta_new,
     $                    dx,dt,time)

      end if

      if (use_strang) then

         ! omegadot for divu computation is average omegadot
         ! over both strang calls
         I_R_divu = I_R

      else

         ! omegadot for divu computation is instantaneous
         ! value at t^{n+1}
         do i=0,nx-1
            do j=1,Nspec
               C(j) = scal_new(i,FirstSpec+j-1)*invmwt(j)
            end do
            call CKWC(scal_new(i,Temp),C,IWRK,RWRK,WDOTK)
            do j=1,Nspec
               I_R_divu(i,j) = WDOTK(j)*mwt(j)/thickFacCH
            end do
         end do

      end if


c     this computes rho D_m     (for species)
c                   lambda / cp (for enthalpy)
c                   lambda      (for temperature)           
      call calc_diffusivities(scal_new,beta_new,mu_new,dx,time+dt)
      call calc_divu(scal_new,beta_new,I_R_divu,divu_new,dx,time+dt)

      do i = 0,nx-1
         rhohalf(i) = 0.5d0*(scal_old(i,Density)+scal_new(i,Density))
         dsdt(i) = (divu_new(i) - divu_old(i)) / dt
      enddo

      print *,'... update velocities'

      vel_theta = 0.5d0

C     get velocity visc terms to use as a forcing term for advection
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
     $                 dx,dt,1,vel_theta,rho_flag,.true.,time)
      endif

      call compute_pthermo(scal_new,pthermo)

c     diagnostics only
      divu_max = ABS(divu_new(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_new(i)))
      enddo
      print *,'DIVU_NEW norm before dpdt = ',divu_max 

      call add_dpdt_nodal(scal_new,pthermo,divu_new,vel_new,dx,dt)

c     diagnostics only
      divu_max = ABS(divu_new(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_new(i)))
      enddo
      print *,'DIVU_NEW norm after dpdt = ',divu_max 

      print *,'...nodal projection...'
      if (initial_iter .eq. 0) then
         call project(vel_new,rhohalf,divu_new,
     $                press_old,press_new,dx,dt)
      endif

      end


      subroutine sdc_advance(macvel,scal_old,scal_new,I_R,
     $                       beta_old,beta_new,dx,dt,time)

      implicit none
      include 'spec.h'
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
c     in the full LMC code, I_R only needs 0:maxspec components
c     component 0 is for rhoh
c     components 1:maxspec are for rhoX
      real*8   I_R(0:nx-1,0:maxspec)
      real*8    macvel(0 :nx  )
      real*8      aofs(0 :nx-1,nscal)
      real*8  beta_old(-1:nx,nscal)
      real*8  beta_new(-1:nx,nscal)
      real*8  mu_dummy(-1:nx)
      real*8    tforce(0 :nx-1,nscal)
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta

c     storage for conservatively corrected species flux for NULN rhoh diffusion term
c     in the full LMC code, these are called fluxNULN
      real*8 spec_flux_lo(0:nx-1,maxspec)
      real*8 spec_flux_hi(0:nx-1,maxspec)

ccccccccccccccccccccccccccccccccc
c     SDC TEMPORARIES - I_R is also a temporary, but allocated above
ccccccccccccccccccccccccccccccccc

c     in the full LMC code, these are also called
c     diff_old, diff_new, and diff_hat.
c     they only contain 0:maxspec components
c     component 0 is for rhoh
c     components 1:maxspec are for rhoX
c     differential diffusion terms for rhoh are stored elsewhere (see below)
      real*8        diff_old(0:nx-1,nscal)
      real*8        diff_new(0:nx-1,nscal)
      real*8        diff_hat(0:nx-1,nscal)

c     in the full LMC code, these are called
c     div_fluxNULN_old, div_fluxNULN_new, div_fluxNULN_hat
      real*8 diffdiff_old(0:nx-1)
      real*8 diffdiff_new(0:nx-1)

c     in the full LMC code, we only need const_src
c     it will only contain 0:maxspec components
c     component 0 is for rhoh
c     components 1:maxspec are for rhoX
      real*8   const_src(0:nx-1,nscal)
      real*8 lin_src_old(0:nx-1,nscal)
      real*8 lin_src_new(0:nx-1,nscal)
      
ccccccccccccccccccccccccccccccccc
c     END SDC TEMPORARIES
ccccccccccccccccccccccccccccccccc

      integer i,n
      
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8      dRhs(0:nx-1,0:maxspec)
      integer is, rho_flag
      integer misdc

c     temporaries for setting I_R = omegadot^n at beginning of time step
      real*8 WDOTK(maxspec), C(maxspec)

      real*8 RWRK
      integer IWRK

      diffdiff_old = 0.d0
      diffdiff_new = 0.d0

      rho_flag = 2

C----------------------------------------------------------------
c     Begin initial predictor
C----------------------------------------------------------------

c     diffusion solves in predictor are regular Crank-Nicolson
      be_cn_theta = 0.5d0

c     compute diffusion term at time n
      print *,'... computing D(U^n)'
c     compute del dot rho D grad Y and make it conservative
c     save species fluxes for differential diffusion
      call get_spec_visc_terms(scal_old,beta_old,
     &                         diff_old(0,FirstSpec),
     &                         spec_flux_lo,spec_flux_hi,dx,time)
c     compute del dot lambda/cp grad h (no differential diffusion)
      call get_rhoh_visc_terms(scal_old,beta_old,
     &                         diff_old(0,RhoH),dx,time)

c     calculate differential diffusion
      if (LeEQ1 .eq. 0) then
c        calculate sum_m del dot h_m (rho D_m - lambda/cp) grad Y_m
c        we pass in conservative rho D grad Y via spec_flux
c        we take lambda / cp from beta
c        we compute h_m from the first scal argument
c        we take the gradient of Y from the second scal argument
         call get_diffdiff_terms(scal_old,scal_old,spec_flux_lo,
     $                           spec_flux_hi,beta_old,diffdiff_old,
     $                           dx,time)
      end if

c     If .true., use I_R in predictor is instantaneous value at t^n
c     If .false., use I_R^lagged = I_R^kmax from previous time step
      if (.false.) then
         do i=0,nx-1
            do n=1,Nspec
               C(n) = scal_old(i,FirstSpec+n-1)*invmwt(n)
            end do
            call CKWC(scal_old(i,Temp),C,IWRK,RWRK,WDOTK)
            do n=1,Nspec
               I_R(i,n) = WDOTK(n)*mwt(n)/thickFacCH
            end do
         end do
      end if

c     compute advective forcing term
      print *,'... computing advective forcing term = D^n + I_R^kmax'
      do i = 0,nx-1
         do n = 1,Nspec
            is = FirstSpec + n - 1
            tforce(i,is) = diff_old(i,is) + I_R(i,n)
         enddo
         tforce(i,RhoH) = diff_old(i,RhoH) + diffdiff_old(i)
      enddo

c     compute advection term
      call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt)

c     update density
      print *,'... update rho'
      call update_rho(scal_old,scal_new,aofs,dx,dt,time)

c     compute part of the RHS for the enthalpy and species
c     diffusion solves
      do i=0,nx-1
         dRhs(i,0) = 0.0d0
         do n=1,Nspec
            dRhs(i,n) = dt*I_R(i,n)
         enddo
      enddo

c     compute RHS for species diffusion solve
      call update_spec(scal_old,scal_new,aofs,alpha,beta_old,dRhs(0,1),
     &                 Rhs(0,FirstSpec),dx,dt,be_cn_theta,time)

C     update species with diffusion solve
      print *,'... do initial diffusion solve for species'
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_new,alpha,beta_old,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta,rho_flag,.false.,time)
      enddo

      if (LeEQ1 .eq. 1) then

c        simply extract D for RhoX
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
               diff_hat(i,is) = 2.d0*((scal_new(i,is)-scal_old(i,is))/dt 
     $              - aofs(i,is) - I_R(i,n) - 0.5d0*diff_old(i,is))
            enddo
         end do

      else

c        compute del dot rho D grad Y and make it conservative
c        save species fluxes for differential diffusion
         call get_spec_visc_terms(scal_new,beta_old,
     $                            diff_hat(0,FirstSpec),
     $                            spec_flux_lo,spec_flux_hi,
     $                            dx,time)

c        update species with conservative diffusion fluxes
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
               scal_new(i,is) = scal_old(i,is) + 
     $              dt*(aofs(i,is) + I_R(i,n)
     $              + 0.5d0*diff_old(i,is) + 0.5d0*diff_hat(i,is))
            end do
         end do
         
c        calculate differential diffusion
c        calculate sum_m del dot h_m (rho D_m - lambda/cp) grad Y_m
c        we pass in conservative rho D grad Y via spec_flux
c        we take lambda / cp from beta
c        we compute h_m from the first scal argument
c        we take the gradient of Y from the second scal argument
         call get_diffdiff_terms(scal_old,scal_new,spec_flux_lo,
     $                           spec_flux_hi,beta_old,diffdiff_new,
     $                           dx,time)

c        add differential diffusion to forcing for enthalpy solve
         do i=0,nx-1
            dRhs(i,0) = dRhs(i,0) 
     $           + 0.5d0*dt*(diffdiff_old(i) + diffdiff_new(i))
         end do

      end if

c     compute RHS for enthalpy diffusion solve
      call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,dRhs(0,0),
     &                 Rhs(0,RhoH),dx,dt,be_cn_theta,time)

c     update enthalpy with diffusion solve
      call cn_solve(scal_new,alpha,beta_old,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta,rho_flag,.false.,time)

c     extract D for RhoH
      do i = 0,nx-1
         diff_hat(i,RhoH) = 2.d0*((scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $        - aofs(i,RhoH) - dRhs(i,0)/dt - 0.5d0*diff_old(i,RhoH) )
      enddo

      if (nochem_hack) then
         write(*,*)'WARNING! doing nochem_hack--skipping reactions'
      else
         print *,'... react with constant sources'
         do n = 1,nscal
            do i = 0,nx-1
               const_src(i,n) = aofs(i,n) 
     $              + 0.5d0*diff_hat(i,n) + 0.5d0*diff_old(i,n)
               lin_src_old(i,n) = 0.d0
               lin_src_new(i,n) = 0.d0
            enddo
         enddo

c        add differential diffusion
         do i=0,nx-1
            const_src(i,RhoH) = const_src(i,RhoH)
     $           + 0.5d0*(diffdiff_old(i)+diffdiff_new(i))
         end do

         call strang_chem(scal_old,scal_new,
     $                    const_src,lin_src_old,lin_src_new,
     $                    I_R,dt)

      endif

C----------------------------------------------------------------
c     End initial predictor
C----------------------------------------------------------------

C----------------------------------------------------------------
c     Begin MISDC iterations
C----------------------------------------------------------------

c     diffusion solves in SDC iterations are iterative corrections
c     that have a backward Euler character
      be_cn_theta = 1.d0

      do misdc = 1, misdc_iterMAX
         print *,'... doing SDC iter ',misdc

         print *,'... compute diff_new = D(U^{n+1,k-1})'
c        this computes rho D_m     (for species)
c                      lambda / cp (for enthalpy)
c                      lambda      (for temperature) 
         call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time+dt)
c        compute del dot rho D grad Y and make it conservative
c        save species fluxes for differential diffusion
         call get_spec_visc_terms(scal_new,beta_new,
     &                            diff_new(0,FirstSpec),
     &                            spec_flux_lo,spec_flux_hi,dx,time+dt)
c        compute del dot lambda/cp grad h (no differential diffusion)
         call get_rhoh_visc_terms(scal_new,beta_new,
     &                            diff_new(0,RhoH),dx,time+dt)

c        calculate differential diffusion
         if (LeEQ1 .eq. 0) then
c           calculate sum_m del dot h_m (rho D_m - lambda/cp) grad Y_m
c           we pass in conservative rho D grad Y via spec_flux
c           we take lambda / cp from beta
c           we compute h_m from the first scal argument
c           we take the gradient of Y from the second scal argument
            call get_diffdiff_terms(scal_new,scal_new,spec_flux_lo,
     $                              spec_flux_hi,beta_new,
     $                              diffdiff_new,dx,time)
         end if

         print *,'... computing advective forcing term = D^n + I_R^k-1'
         do i = 0,nx-1
            do n = 1,Nspec
               is = FirstSpec + n - 1
               tforce(i,is) = diff_old(i,is) + I_R(i,n)
            enddo
c           really no need to recompute this since it doesn't change
            tforce(i,RhoH) = diff_old(i,RhoH) + diffdiff_old(i)
         enddo
         
         print *,'... compute A with updated D+R source'
         call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt)

         print *,'... update rho'
         call update_rho(scal_old,scal_new,aofs,dx,dt,time)

         print *,'... update D for species with A + R + MISDC(D)'
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
c              includes deferred correction term for species
               dRhs(i,n) = dt*(I_R(i,n) 
     &              + 0.5d0*(diff_old(i,is) - diff_new(i,is)))
            enddo
c           includes deferred correction term for enthalpy
c           differential diffusion will be added later
            dRhs(i,0) = dt*(
     &           + 0.5d0*(diff_old(i,RhoH) - diff_new(i,RhoH)))
         enddo
         call update_spec(scal_old,scal_new,aofs,alpha,beta_old,
     &                    dRhs(0,1),Rhs(0,FirstSpec),dx,dt,
     &                    be_cn_theta,time)

         do n=1,Nspec
            is = FirstSpec + n - 1
            call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $                    dx,dt,is,be_cn_theta,rho_flag,.false.,time)
         enddo

         if (LeEQ1 .eq. 1) then

c           simply extract D for RhoX
            do i=0,nx-1
               do n=1,Nspec
                  is = FirstSpec + n - 1
                  diff_hat(i,is) = (scal_new(i,is)-scal_old(i,is))/dt 
     $                 - aofs(i,is) - dRhs(i,n)/dt
               enddo
            enddo

         else

c           compute del dot rho D grad Y and make it conservative
c           save species fluxes for differential diffusion
            call get_spec_visc_terms(scal_new,beta_new,
     $                               diff_hat(0,FirstSpec),
     $                               spec_flux_lo,spec_flux_hi,
     $                               dx,time)

c           add differential diffusion to forcing for enthalpy solve
            do i=0,nx-1
               dRhs(i,0) = dRhs(i,0) 
     $              + 0.5d0*dt*(diffdiff_old(i) + diffdiff_new(i))
            end do

         end if

         print *,'... update D for rhoh with A + R + MISDC(D)'
         call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,
     &                    dRhs(0,0),Rhs(0,RhoH),dx,dt,be_cn_theta,time)
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,RhoH),
     $                 dx,dt,RhoH,be_cn_theta,rho_flag,.false.,time)
         print *,'... create new temp from new RhoH, spec'

c        extract D for RhoH
         do i = 0,nx-1
            diff_hat(i,RhoH) = (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $           - aofs(i,RhoH) - dRhs(i,0)/dt
         enddo
         
         if (nochem_hack) then
            write(*,*)'WARNING:: SDC nochem_hack--skipping reactions'
         else
            print *,'... react with const sources'
            do n = 1,nscal
               do i = 0,nx-1
                  const_src(i,n) = aofs(i,n)
     $                 + 0.5d0*(diff_old(i,n)+diff_new(i,n))
     $                 + diff_hat(i,n) - diff_new(i,n)

                  lin_src_old(i,n) = 0.d0
                  lin_src_new(i,n) = 0.d0
               enddo
            enddo

c           add differential diffusion
            do i=0,nx-1
               const_src(i,RhoH) = const_src(i,RhoH)
     $              + 0.5d0*(diffdiff_old(i)+diffdiff_new(i))
            end do
            
            call strang_chem(scal_old,scal_new,
     $                       const_src,lin_src_old,lin_src_new,
     $                       I_R,dt)

         endif

C----------------------------------------------------------------
c     End MISDC iterations
C----------------------------------------------------------------

      enddo

      end

      subroutine strang_advance(macvel,scal_old,scal_new,I_R,
     $                          beta_old,beta_new,dx,dt,time)

      implicit none
      include 'spec.h'
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
      real*8   I_R(0:nx-1,0:maxspec)
      real*8    macvel(0 :nx  )
      real*8      aofs(0 :nx-1,nscal)
      real*8  beta_old(-1:nx,nscal)
      real*8  beta_new(-1:nx,nscal)
      real*8  mu_dummy(-1:nx)
      real*8    tforce(0 :nx-1,nscal)
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta
      
      real*8    diff_old(0:nx-1,nscal)
      real*8    diff_tmp(0:nx-1,nscal)
      real*8   const_src(0:nx-1,nscal)
      real*8 lin_src_old(0:nx-1,nscal)
      real*8 lin_src_new(0:nx-1,nscal)

      real*8   I_R_temp(0:nx-1,0:maxspec)
      
      integer i,n
      
      real*8     alpha(0:nx-1)
      real*8       Rhs(0:nx-1,nscal)
      real*8      dRhs(0:nx-1,0:maxspec)
      real*8 rhocp
      real*8 Y(maxspec)
      real*8 RWRK, cpmix
      integer IWRK, is, rho_flag

      real*8 spec_flux_lo(0:nx-1,maxspec)
      real*8 spec_flux_hi(0:nx-1,maxspec)

      real*8 diffdiff_old(0:nx-1)
      real*8 diffdiff_new(0:nx-1)

      diffdiff_old = 0.d0
      diffdiff_new = 0.d0

      be_cn_theta = 0.5d0

      if (nochem_hack) then
         print *,'WARNING! doing nochem_hack...'
         do n = 1,nscal
            do i = 0,nx-1
               I_R(i,n) = 0.d0
            enddo
         enddo
      else
         print *,'... react for dt/2;  set I_R'
         do n = 1,nscal
            do i = 0,nx-1
               const_src(i,n) =   0.d0
               lin_src_old(i,n) = 0.d0
               lin_src_new(i,n) = 0.d0
            enddo
         enddo
         call strang_chem(scal_old,scal_new,
     $                    const_src,lin_src_old,lin_src_new,
     $                    I_R,dt/2.d0)
         do i=0,nx-1
            do n = FirstSpec,LastSpec
               scal_old(i,n) = scal_new(i,n)
            enddo
         enddo
      endif

c     we only care about updated species out of strang_chem
c     rho and rhoh remain constant
c     call the EOS to get consistent temperature
      call rhoh_to_temp(scal_old)
c     
c*****************************************************************
c     

      print *,'... creating the diffusive terms with old data'
c    compute rho^(1) D_m^(1)     (for species)
c            lambda^(1) / cp^(1) (for enthalpy)
c            lambda^(1)          (for temperature) 
      call calc_diffusivities(scal_old,beta_old,mu_dummy,dx,time)

c     compute del dot lambda grad T + rho D grad h dot grad Y
c     the rho D grad Y term is now computed conservatively
      call get_temp_visc_terms(scal_old,beta_old,
     &                         diff_old(0,Temp),dx,time)
c     compute del dot rho D grad Y and make it conservative
c     save species fluxes for differential diffusion
      call get_spec_visc_terms(scal_old,beta_old,
     &                         diff_old(0,FirstSpec),
     &                         spec_flux_lo,spec_flux_hi,dx,time)
c     compute del dot lambda/cp grad h (no differential diffusion)
      call get_rhoh_visc_terms(scal_old,beta_old,
     &                         diff_old(0,RhoH),dx,time)

c     calculate differential diffusion
      if (LeEQ1 .eq. 0) then
c        calculate sum_m del dot h_m (rho D_m - lambda/cp) grad Y_m
c        we pass in conservative rho D grad Y via spec_flux
c        we take lambda / cp from beta
c        we compute h_m from the first scal argument
c        we take the gradient of Y from the second scal argument
         call get_diffdiff_terms(scal_old,scal_old,spec_flux_lo,
     $                           spec_flux_hi,beta_old,diffdiff_old,
     $                           dx,time)
      end if
            
      print *,'... computing aofs with D(old)'

      do i = 0,nx-1
         do n = 1,Nspec
            is = FirstSpec + n - 1
            tforce(i,is) = diff_old(i,is)
         enddo
         tforce(i,RhoH) = diff_old(i,RhoH) + diffdiff_old(i)
      enddo
       
      call scal_aofs(scal_old,macvel,aofs,tforce,dx,dt)

      print *,'... update rho'
      call update_rho(scal_old,scal_new,aofs,dx,dt,time)

      do i = 0,nx-1
         do n = 1,Nspec
            Y(n) = scal_old(i,FirstSpec+n-1) / scal_old(i,Density)
         enddo
         call CKCPBS(scal_old(i,Temp),Y,IWRK,RWRK,cpmix)
         rhocp = cpmix * 
     &        (scal_old(i,Density) + scal_new(i,Density)) / 2.d0
         tforce(i,Temp) = diff_old(i,Temp)/rhocp
      end do

c*****************************************************************
c     Either do c-n solve for new T prior to computing new 
c     coeffs, or simply start by copying from previous time step
      if (predict_temp_for_coeffs .eq. 1) then
         print *,'... predict temp with old coeffs'
         rho_flag = 1
         do i=0,nx-1
            dRhs(i,0) = 0.0d0
         enddo
         call update_temp(scal_old,scal_new,aofs,
     $                    alpha,beta_old,beta_old,dRhs(0,0),
     $                    Rhs(0,Temp),dx,dt,be_cn_theta,time)
c        just uses RHS and overwrites snew
c        does not fill ghost cells
         call cn_solve(scal_new,alpha,beta_old,Rhs(0,Temp),
     $                 dx,dt,Temp,be_cn_theta,rho_flag,.false.,time)

         print *,'... compute new coeffs'
c        compute rho^(2) D_m^(2),* (for species)
c        lambda/cp (for enthalpy) won't be used
c        lambda^(1) (for temperature) won't be used
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
     &                 dRhs(0,1),Rhs(0,FirstSpec),dx,dt,
     &                 be_cn_theta,time)

      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta,rho_flag,.false.,time)
      enddo

      if (LeEQ1 .eq. 0) then

c     compute del dot rho D grad Y and make it conservative
c     save species fluxes for differential diffusion
         call get_spec_visc_terms(scal_new,beta_new,
     &                            diff_tmp(0,FirstSpec),
     &                            spec_flux_lo,spec_flux_hi,
     &                            dx,time)

c     update species with conservative diffusion fluxes
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
               scal_new(i,is) = scal_old(i,is) + 
     $              dt*(aofs(i,is)
     $              + 0.5d0*diff_old(i,is) + 0.5d0*diff_tmp(i,is))
            end do
         end do

      end if

c     this computes rho D_m                 (for species) won't be used
c                   lambda^(2),* / cp^(2),* (for enthalpy)
c                   lambda^(2),*            (for temperature) 
      call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time+dt)

      if (LeEQ1 .eq. 0) then

c     calculate differential diffusion
c     calculate sum_m del dot h_m (rho D_m - lambda/cp) grad Y_m
c     we pass in conservative rho D grad Y via spec_flux
c     we take lambda / cp from beta
c     we compute h_m from the first scal argument
c     we take the gradient of Y from the second scal argument
         call get_diffdiff_terms(scal_new,scal_new,spec_flux_lo,
     $                           spec_flux_hi,beta_new,
     $                           diffdiff_new,dx,time)
         
         do i=0,nx-1
            dRhs(i,0) = dRhs(i,0)
     $           + 0.5d0*dt*(diffdiff_old(i) + diffdiff_new(i))
         end do
         
      end if

      print *,'... do predictor for rhoh (MISDC terms=0)'
      call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,
     &                 dRhs(0,0),Rhs(0,RhoH),dx,dt,be_cn_theta,time)
      rho_flag = 2
      call cn_solve(scal_new,alpha,beta_new,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta,rho_flag,.false.,time)

      call rhoh_to_temp(scal_new)

C----------------------------------------------------------------
C   Corrector

      print *,'... compute new coeffs'
c     this computes rho^(2) D_m^(2)     (for species)
c                   lambda^(2) / cp^(2) (for enthalpy)
c                   lambda^(2)          (for temperature) 
      call calc_diffusivities(scal_new,beta_new,mu_dummy,dx,time+dt)

      do i=0,nx-1
         dRhs(i,0) = 0.0d0
         do n=1,Nspec
            dRhs(i,n) = 0.d0
         enddo
      enddo
      call update_spec(scal_old,scal_new,aofs,alpha,beta_old,
     &                 dRhs(0,1),Rhs(0,FirstSpec),dx,dt,
     &                 be_cn_theta,time)
      rho_flag = 2
      do n=1,Nspec
         is = FirstSpec + n - 1
         call cn_solve(scal_new,alpha,beta_new,Rhs(0,is),
     $                 dx,dt,is,be_cn_theta,rho_flag,.false.,time)
      enddo

      if (LeEQ1 .eq. 0) then

c     compute del dot rho D grad Y and make it conservative
c     save species fluxes for differential diffusion
         call get_spec_visc_terms(scal_new,beta_new,
     &                            diff_tmp(0,FirstSpec),
     &                            spec_flux_lo,spec_flux_hi,
     &                            dx,time)

c     update species with conservative diffusion fluxes
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec + n - 1
               scal_new(i,is) = scal_old(i,is) + 
     $              dt*(aofs(i,is)
     $              + 0.5d0*diff_old(i,is) + 0.5d0*diff_tmp(i,is))
            end do
         end do

c     calculate differential diffusion
c     calculate sum_m del dot h_m (rho D_m - lambda/cp) grad Y_m
c     we pass in conservative rho D grad Y via spec_flux
c     we take lambda / cp from beta
c     we compute h_m from the first scal argument
c     we take the gradient of Y from the second scal argument
         call get_diffdiff_terms(scal_new,scal_new,spec_flux_lo,
     $                           spec_flux_hi,beta_new,
     $                           diffdiff_new,dx,time)

         do i=0,nx-1
            dRhs(i,0) = dRhs(i,0)
     $           + 0.5d0*dt*(diffdiff_old(i) + diffdiff_new(i))
         end do
         
      end if
         
      print *,'... do predictor for rhoh (MISDC terms=0)'
      call update_rhoh(scal_old,scal_new,aofs,alpha,beta_old,
     &                 dRhs(0,0),Rhs(0,RhoH),dx,dt,be_cn_theta,time)
      rho_flag = 2
      call cn_solve(scal_new,alpha,beta_new,Rhs(0,RhoH),
     $              dx,dt,RhoH,be_cn_theta,rho_flag,.false.,time)
      call rhoh_to_temp(scal_new)

      if (nochem_hack) then
         print *,'WARNING! doing nochem_hack...'
      else
         do i=0,nx-1
            do n = FirstSpec,LastSpec
               scal_old(i,n) = scal_new(i,n)
            enddo
            scal_old(i,Temp) = scal_new(i,Temp)
            scal_old(i,Density) = scal_new(i,Density)
         enddo

         call strang_chem(scal_old,scal_new,
     $                    const_src,lin_src_old,lin_src_new,
     $                    I_R_temp,dt/2.d0)

c        we only care about updated species out of strang_chem
c        rho and rhoh remain constant
c        call the EOS to get consistent temperature
         call rhoh_to_temp(scal_new)

         I_R = I_R + I_R_temp
         I_R = I_R / 2.d0
      endif

      end
